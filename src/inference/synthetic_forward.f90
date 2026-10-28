module tiago_synthetic_forward
    !! One forward evaluation of every configured diagnostic from a VMEC
    !! equilibrium, real diagnostic geometry and tabulated kinetic profiles,
    !! written to one netCDF-4 file with a group per diagnostic (see
    !! docs/synthetic-forward.md). The profiles enter as fixed (non-latent)
    !! spline coefficients; no new physics is added here.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_c_binding, only: c_int, c_size_t
    use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_get_var, nf90_noerr, nf90_nowrite
    use tiago_netcdf_c, only: NC_NOERR, NC_CLOBBER, NC_NETCDF4, NC_GLOBAL, NC_CHAR, NC_DOUBLE, &
        nc_create, nc_close, nc_def_grp, nc_def_dim, nc_def_var, nc_put_var_double, &
        nc_put_var_text, nc_put_att_double, put_text_att, cstr
    use tiago_real_space_state, only: real_space_state_t, located_t, locate_ok
    use tiago_vmec_state, only: vmec_state_t
    use tiago_state_sampler, only: sample_state
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t, Q_S, Q_B, Q_NE, Q_TE, Q_TI, Q_ZEFF
    use tiago_kinetic_profiles, only: kinetic_profiles_t, N_PROFILES, P_NE, P_TE, P_TI, P_ZEFF
    use tiago_sight_geometry, only: chord_t, to_rpz
    use tiago_line_integrals, only: interferometer_t, xics_t
    use tiago_thomson_scattering, only: thomson_t, thomson_volume_t
    use tiago_cold_resonance, only: cold_resonance
    use tiago_profile_geometry_io, only: has_group, read_chords, read_thomson_volumes, &
        read_ece_lines, ece_line_t, read_mse_lines, mse_line_t, profile_table_t, &
        read_profile_table
    use tiago_diagnostic_types, only: mse_channel_t
    use tiago_mse_state, only: mse_state_angles
    use tiago_sensor_set, only: sensor_set_t, read_sensor_set, sensor_kind_name
    use tiago_sensor_signals, only: sensor_set_signals
    use tiago_vacuum_forward, only: vacuum_solver_t
    use tiago_sha256, only: sha256_file
    use tiago_revision, only: tiago_source_revision
    implicit none
    private

    type, public :: synth_config_t
        character(len=:), allocatable :: wout, output
        character(len=:), allocatable :: profiles, geometry !! '' = none
        character(len=:), allocatable :: sensors, coils, extcur !! '' = none
        character(len=:), allocatable :: sheet !! contravariant | covariant | conservative
        integer :: npanel = 200 !! panels per chord (4 Gauss points each)
        integer :: thomson_points = 4 !! Gauss points per scattering volume
        integer :: ece_points = 4000 !! samples per ECE line of sight
        integer :: profile_knots = 12 !! spline intervals of the profiles in s
        integer :: plasma_nphi = 64, plasma_ntheta = 64
    end type synth_config_t

    type, public :: synth_state_t
        !! Equilibrium and profiles shared by every instrument.
        class(real_space_state_t), allocatable :: state
        type(kinetic_profiles_t) :: profiles
        real(dp), allocatable :: z(:)
        real(dp) :: aminor = 0.0_dp !! VMEC Aminor_p [m]
    end type synth_state_t

    type :: named_t
        character(len=:), allocatable :: name, units, description
        real(dp), allocatable :: values(:)
    end type named_t

    public :: synth_run, synth_load, set_profile_table, sample_kinetic, ece_cold_resonance, &
        synth_write

contains

    subroutine synth_run(cfg, ierr, message)
        type(synth_config_t), intent(in) :: cfg
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        type(synth_state_t) :: st

        call synth_load(cfg, st, ierr, message)
        if (ierr /= 0) return
        call synth_write(cfg, st, ierr, message)
    end subroutine synth_run

    subroutine synth_write(cfg, st, ierr, message)
        !! Every instrument of cfg on a loaded state, written to cfg%output.
        type(synth_config_t), intent(in) :: cfg
        type(synth_state_t), intent(in) :: st
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        integer(c_int) :: ncid

        ierr = 1
        if (nc_create(cstr(cfg%output), ior(NC_CLOBBER, NC_NETCDF4), ncid) /= NC_NOERR) then
            message = 'cannot create '//cfg%output
            return
        end if
        call provenance(ncid, cfg, st)
        ierr = 0
        if (len(cfg%geometry) > 0) then
            if (has_group(cfg%geometry, 'interferometer')) &
                call run_interferometer(ncid, cfg, st, ierr, message)
            if (ierr == 0 .and. has_group(cfg%geometry, 'thomson')) &
                call run_thomson(ncid, cfg, st, ierr, message)
            if (ierr == 0 .and. has_group(cfg%geometry, 'xics')) &
                call run_xics(ncid, cfg, st, ierr, message)
            if (ierr == 0 .and. has_group(cfg%geometry, 'ece')) &
                call run_ece(ncid, cfg, st, ierr, message)
            if (ierr == 0 .and. has_group(cfg%geometry, 'mse')) &
                call run_mse(ncid, cfg, st, ierr, message)
        end if
        if (ierr == 0 .and. len(cfg%sensors) > 0) call run_magnetics(ncid, cfg, ierr, message)
        if (nc_close(ncid) /= NC_NOERR .and. ierr == 0) then
            ierr = 1
            message = 'cannot close '//cfg%output
        end if
    end subroutine synth_write

    subroutine synth_load(cfg, st, ierr, message)
        type(synth_config_t), intent(in) :: cfg
        type(synth_state_t), intent(out) :: st
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        type(profile_table_t) :: table
        type(vmec_state_t), allocatable :: vmec
        integer :: k

        allocate(vmec)
        call vmec%load(cfg%wout, ierr, message, vector=.true.)
        if (ierr /= 0) return
        call move_alloc(vmec, st%state)
        call read_aminor(cfg%wout, st%aminor, ierr, message)
        if (ierr /= 0) return
        call st%profiles%init([(real(k, dp) / cfg%profile_knots, k = 0, cfg%profile_knots)])
        allocate(st%z(0))
        if (len(cfg%profiles) == 0) return
        call read_profile_table(cfg%profiles, table, ierr, message)
        if (ierr /= 0) return
        call set_profile_table(st, table, ierr, message)
        if (ierr /= 0) message = cfg%profiles//': '//message
    end subroutine synth_load

    subroutine set_profile_table(st, table, ierr, message)
        !! Fixed profiles from a table; Zeff stays 1 when the table has none.
        type(synth_state_t), intent(inout) :: st
        type(profile_table_t), intent(in) :: table
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        integer :: code

        ierr = 0
        message = ''
        call set_one(P_NE, table%ne, 'ne')
        if (ierr == 0) call set_one(P_TE, table%te, 'te')
        if (ierr == 0) call set_one(P_TI, table%ti, 'ti')
        if (ierr == 0 .and. allocated(table%zeff)) call set_one(P_ZEFF, table%zeff, 'zeff')
    contains
        subroutine set_one(which, values, name)
            integer, intent(in) :: which
            real(dp), intent(in) :: values(:)
            character(len=*), intent(in) :: name
            call st%profiles%set_table(which, table%s, values, code)
            if (code == 1) then
                ierr = 1
                message = name//' outside its domain (positive; 1 < zeff < z_imp)'
            else if (code /= 0) then
                ierr = 1
                message = name//' needs at least as many points as '// &
                    'spline coefficients (profile_knots + 3) and distinct s'
            end if
        end subroutine set_one
    end subroutine set_profile_table

    subroutine sample_kinetic(st, x, sample, loc)
        !! Provider rows by tiago_state_sampler, kinetic rows from the
        !! profiles at s. Points that cannot be located get zero density
        !! and temperatures and Zeff = 1.
        type(synth_state_t), intent(in) :: st
        real(dp), intent(in) :: x(:, :)
        type(plasma_sample_t), intent(out) :: sample
        type(located_t), allocatable, intent(out) :: loc(:)
        real(dp), allocatable :: f(:, :), df(:, :)
        integer :: k

        call sample_state(st%state, x, sample, loc)
        allocate(f(N_PROFILES, sample%n), df(N_PROFILES, sample%n))
        call st%profiles%evaluate(st%z, sample%q(Q_S, :), f, df)
        do k = 1, sample%n
            if (loc(k)%status /= locate_ok) then
                sample%q(Q_NE:Q_TI, k) = 0.0_dp
                sample%q(Q_ZEFF, k) = 1.0_dp
                cycle
            end if
            sample%q(Q_NE, k) = f(P_NE, k)
            sample%q(Q_TE, k) = f(P_TE, k)
            sample%q(Q_TI, k) = f(P_TI, k)
            sample%q(Q_ZEFF, k) = f(P_ZEFF, k)
            sample%dq_dx(:, Q_NE, k) = df(P_NE, k) * sample%dq_dx(:, Q_S, k)
            sample%dq_dx(:, Q_TE, k) = df(P_TE, k) * sample%dq_dx(:, Q_S, k)
            sample%dq_dx(:, Q_TI, k) = df(P_TI, k) * sample%dq_dx(:, Q_S, k)
            sample%dq_dx(:, Q_ZEFF, k) = df(P_ZEFF, k) * sample%dq_dx(:, Q_S, k)
        end do
    end subroutine sample_kinetic

    subroutine run_interferometer(ncid, cfg, st, ierr, message)
        integer(c_int), intent(in) :: ncid
        type(synth_config_t), intent(in) :: cfg
        type(synth_state_t), intent(in) :: st
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        type(chord_t), allocatable :: chords(:)
        type(interferometer_t) :: diag
        real(dp), allocatable :: y(:)

        call read_chords(cfg%geometry, 'interferometer', cfg%npanel, chords, ierr, message)
        if (ierr /= 0) return
        call diag%init(chords)
        call predict(diag, st, y)
        call write_group(ncid, 'interferometer', 'chord', chord_labels(chords), ierr, message, &
            [named('line_density', 'm-2', 'int ne dl along the chord', y)])
    end subroutine run_interferometer

    subroutine run_xics(ncid, cfg, st, ierr, message)
        integer(c_int), intent(in) :: ncid
        type(synth_config_t), intent(in) :: cfg
        type(synth_state_t), intent(in) :: st
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        type(chord_t), allocatable :: chords(:)
        type(xics_t) :: diag
        real(dp), allocatable :: y(:)
        integer :: n

        call read_chords(cfg%geometry, 'xics', cfg%npanel, chords, ierr, message)
        if (ierr /= 0) return
        call diag%init(chords)
        call predict(diag, st, y)
        n = size(chords)
        call write_group(ncid, 'xics', 'chord', chord_labels(chords), ierr, message, &
            [named('brightness', '1', 'relative brightness, unit scale per chord', y(:n)), &
            named('ti', 'eV', 'emissivity-weighted Ti along the chord', y(n + 1:))])
    end subroutine run_xics

    subroutine run_thomson(ncid, cfg, st, ierr, message)
        integer(c_int), intent(in) :: ncid
        type(synth_config_t), intent(in) :: cfg
        type(synth_state_t), intent(in) :: st
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        type(thomson_volume_t), allocatable :: volumes(:)
        type(thomson_t) :: diag
        type(plasma_sample_t) :: centre
        type(located_t), allocatable :: loc(:)
        real(dp), allocatable :: y(:), x(:, :)
        character(len=:), allocatable :: labels(:)
        integer :: n, c, width

        call read_thomson_volumes(cfg%geometry, cfg%thomson_points, volumes, ierr, message)
        if (ierr /= 0) return
        call diag%init(volumes)
        call predict(diag, st, y)
        n = size(volumes)
        allocate(x(3, n))
        do c = 1, n
            x(:, c) = to_rpz(volumes(c)%center)
        end do
        call sample_kinetic(st, x, centre, loc)
        width = maxval([(len(volumes(c)%label), c = 1, n)])
        allocate(character(len=width) :: labels(n))
        do c = 1, n
            labels(c) = volumes(c)%label
        end do
        call write_group(ncid, 'thomson', 'volume', labels, ierr, message, &
            [named('ne', 'm-3', 'volume-averaged ne', y(:n)), &
            named('te', 'eV', 'density-weighted Te', y(n + 1:)), &
            named('s', '1', 'normalized toroidal flux at the volume centre', centre%q(Q_S, :)), &
            named('r_eff', 'm', 'Aminor_p sqrt(s) at the volume centre', &
            st%aminor * sqrt(max(centre%q(Q_S, :), 0.0_dp))), &
            named('inside', '1', '1 if the centre is inside the boundary', &
            merge(1.0_dp, 0.0_dp, centre%inside))])
    end subroutine run_thomson

    subroutine run_ece(ncid, cfg, st, ierr, message)
        integer(c_int), intent(in) :: ncid
        type(synth_config_t), intent(in) :: cfg
        type(synth_state_t), intent(in) :: st
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        type(ece_line_t), allocatable :: lines(:)
        real(dp), allocatable :: te(:), s(:), rpz(:, :)
        character(len=:), allocatable :: labels(:)
        integer :: n, c, width

        call read_ece_lines(cfg%geometry, lines, ierr, message)
        if (ierr /= 0) return
        n = size(lines)
        call ece_cold_resonance(st, lines, cfg%ece_points, te, s, rpz)
        width = maxval([(len(lines(c)%label), c = 1, n)])
        allocate(character(len=width) :: labels(n))
        do c = 1, n
            labels(c) = lines(c)%label
        end do
        call write_group(ncid, 'ece', 'channel', labels, ierr, message, &
            [named('te', 'eV', 'Te at the first cold resonance inside the plasma (NaN if none)', te), &
            named('s', '1', 'normalized toroidal flux at the resonance', s), &
            named('r_eff', 'm', 'Aminor_p sqrt(s) at the resonance', st%aminor * sqrt(s)), &
            named('R', 'm', 'major radius of the resonance', rpz(1, :)), &
            named('phi', 'rad', 'toroidal angle of the resonance', rpz(2, :)), &
            named('Z', 'm', 'height of the resonance', rpz(3, :))])
    end subroutine run_ece

    subroutine ece_cold_resonance(st, lines, npoints, te, s, rpz)
        !! First point along each line of sight, inside the plasma, where
        !! harmonic f_ce(|B|) = f (cold plasma, no relativistic shift), by
        !! linear interpolation between npoints samples; Te and s there.
        type(synth_state_t), intent(in) :: st
        type(ece_line_t), intent(in) :: lines(:)
        integer, intent(in) :: npoints
        real(dp), allocatable, intent(out) :: te(:), s(:), rpz(:, :)
        type(plasma_sample_t) :: smp, res
        type(located_t), allocatable :: loc(:)
        real(dp) :: x(3, npoints), d(3), b_res, unused, t, nan
        integer :: c, j

        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(te(size(lines)), s(size(lines)), rpz(3, size(lines)))
        te = nan
        s = nan
        rpz = nan
        do c = 1, size(lines)
            d = lines(c)%direction / norm2(lines(c)%direction)
            do j = 1, npoints
                x(:, j) = to_rpz(lines(c)%origin + lines(c)%length * (j - 1) / (npoints - 1) * d)
            end do
            call sample_state(st%state, x, smp, loc)
            call cold_resonance(lines(c)%frequency, real(lines(c)%harmonic, dp), b_res, unused)
            do j = 2, npoints
                if (.not. (smp%inside(j) .and. smp%inside(j - 1))) cycle
                if ((smp%q(Q_B, j) - b_res) * (smp%q(Q_B, j - 1) - b_res) > 0.0_dp) cycle
                t = (b_res - smp%q(Q_B, j - 1)) / (smp%q(Q_B, j) - smp%q(Q_B, j - 1))
                t = lines(c)%length * (j - 2 + t) / (npoints - 1)
                call sample_kinetic(st, reshape(to_rpz(lines(c)%origin + t * d), [3, 1]), res, loc)
                if (.not. res%inside(1)) exit
                te(c) = res%q(Q_TE, 1)
                s(c) = res%q(Q_S, 1)
                rpz(:, c) = res%x(:, 1)
                exit
            end do
        end do
    end subroutine ece_cold_resonance

    subroutine run_mse(ncid, cfg, st, ierr, message)
        !! Pi-polarization angle from the interior field at each channel's
        !! sight-beam midpoint (tiago_mse_state, docs/mse.md conventions,
        !! reference +Z).
        integer(c_int), intent(in) :: ncid
        type(synth_config_t), intent(in) :: cfg
        type(synth_state_t), intent(in) :: st
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        type(mse_line_t), allocatable :: lines(:)
        type(mse_channel_t), allocatable :: channels(:)
        real(dp), allocatable :: gamma(:), field(:, :), rho(:)
        character(len=:), allocatable :: labels(:)
        integer :: c, n, width
        integer(i32) :: code

        call read_mse_lines(cfg%geometry, lines, ierr, message)
        if (ierr /= 0) return
        n = size(lines)
        allocate(channels(n))
        do c = 1, n
            channels(c)%label = lines(c)%label
            channels(c)%beam_point = lines(c)%beam_origin
            channels(c)%beam_direction = lines(c)%beam_direction
            channels(c)%sight_point = lines(c)%origin
            channels(c)%sight_direction = -lines(c)%direction
            channels(c)%max_gap = lines(c)%max_gap
        end do
        call mse_state_angles(st%state, channels, gamma, code, message, field=field, rho=rho)
        ierr = int(code)
        if (ierr /= 0) then
            message = cfg%geometry//':mse: '//message
            return
        end if
        width = maxval([(len(lines(c)%label), c = 1, n)])
        allocate(character(len=width) :: labels(n))
        do c = 1, n
            labels(c) = lines(c)%label
        end do
        call write_group(ncid, 'mse', 'channel', labels, ierr, message, &
            [named('gamma', 'rad', 'pi-polarization angle at the sight-beam midpoint', gamma), &
            named('s', '1', 'normalized toroidal flux at the midpoint', rho**2), &
            named('b_x', 'T', 'Cartesian B at the midpoint', field(1, :)), &
            named('b_y', 'T', 'Cartesian B at the midpoint', field(2, :)), &
            named('b_z', 'T', 'Cartesian B at the midpoint', field(3, :))])
    end subroutine run_mse

    subroutine run_magnetics(ncid, cfg, ierr, message)
        integer(c_int), intent(in) :: ncid
        type(synth_config_t), intent(in) :: cfg
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        type(sensor_set_t) :: set
        type(vacuum_solver_t) :: solver
        real(dp), allocatable :: y(:), kinds(:)
        character(len=:), allocatable :: labels(:)
        integer :: k, width

        call read_sensor_set(cfg%sensors, set, ierr, message)
        if (ierr /= 0) return
        if (len(cfg%coils) > 0 .and. len(cfg%extcur) > 0) then
            call solver%init(cfg%coils, cfg%extcur)
        else
            call solver%init(cfg%coils)
        end if
        call solver%enable_plasma_from_vmec(cfg%wout, cfg%plasma_nphi, cfg%plasma_ntheta, &
            covariant=cfg%sheet /= 'contravariant', conservative=cfg%sheet == 'conservative')
        call sensor_set_signals(solver, set, y)
        width = maxval([(len(set%labels(k)%s), k = 1, size(set%labels))])
        allocate(character(len=width) :: labels(size(set%labels)))
        do k = 1, size(set%labels)
            labels(k) = set%labels(k)%s
        end do
        kinds = real(set%kinds, dp)
        call write_group(ncid, 'magnetics', 'signal', labels, ierr, message, &
            [named('signal', 'Wb or T m2', &
            'coil plus plasma signal (flux loops: Wb; Rogowski and B-probes: T m2)', y), &
            named('kind', '1', 'sensor kind: '//kind_list(), kinds)])
        call solver%finalize()
    contains
        function kind_list() result(text)
            character(len=:), allocatable :: text
            integer :: j
            character(len=4) :: num
            text = ''
            do j = 1, 7
                write(num, '(I0)') j
                text = text//trim(num)//'='//trim(sensor_kind_name(j))//merge(', ', '  ', j < 7)
            end do
            text = trim(text)
        end function kind_list
    end subroutine run_magnetics

    subroutine predict(diag, st, y)
        class(profile_diagnostic_t), intent(in) :: diag
        type(synth_state_t), intent(in) :: st
        real(dp), allocatable, intent(out) :: y(:)
        type(plasma_sample_t) :: smp
        type(located_t), allocatable :: loc(:)
        real(dp), allocatable :: x(:, :), dx(:, :, :), eta(:)

        allocate(eta(diag%n_nuisance()), y(diag%n_signals()))
        eta = 0.0_dp
        call diag%sample_points(eta, x, dx)
        call sample_kinetic(st, x, smp, loc)
        call diag%predict(smp, eta, y)
    end subroutine predict

    type(named_t) function named(name, units, description, values)
        character(len=*), intent(in) :: name, units, description
        real(dp), intent(in) :: values(:)
        named%name = name
        named%units = units
        named%description = description
        named%values = values
    end function named

    function chord_labels(chords) result(labels)
        type(chord_t), intent(in) :: chords(:)
        character(len=:), allocatable :: labels(:)
        integer :: c, width
        width = maxval([(len(chords(c)%label), c = 1, size(chords))])
        allocate(character(len=width) :: labels(size(chords)))
        do c = 1, size(chords)
            labels(c) = chords(c)%label
        end do
    end function chord_labels

    subroutine write_group(ncid, group, dimname, labels, ierr, message, vars)
        integer(c_int), intent(in) :: ncid
        character(len=*), intent(in) :: group, dimname
        character(len=*), intent(in) :: labels(:)
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        type(named_t), intent(in) :: vars(:)
        integer(c_int) :: gid, dimid, strid, varid, ids(size(vars))
        integer :: k

        ierr = 1
        message = 'cannot write group '//group//' of the output'
        if (nc_def_grp(ncid, cstr(group), gid) /= NC_NOERR) return
        if (nc_def_dim(gid, cstr(dimname), int(size(labels), c_size_t), dimid) /= NC_NOERR) return
        if (nc_def_dim(gid, cstr('strlen'), int(max(len(labels), 1), c_size_t), strid) &
            /= NC_NOERR) return
        if (nc_def_var(gid, cstr('label'), NC_CHAR, 2_c_int, [dimid, strid], varid) &
            /= NC_NOERR) return
        do k = 1, size(vars)
            if (nc_def_var(gid, cstr(vars(k)%name), NC_DOUBLE, 1_c_int, [dimid], ids(k)) &
                /= NC_NOERR) return
            if (put_text_att(gid, ids(k), 'units', vars(k)%units) /= NC_NOERR) return
            if (put_text_att(gid, ids(k), 'long_name', vars(k)%description) /= NC_NOERR) return
        end do
        if (size(labels) > 0 .and. len(labels) > 0) then
            if (nc_put_var_text(gid, varid, labels) /= NC_NOERR) return
        end if
        do k = 1, size(vars)
            if (nc_put_var_double(gid, ids(k), vars(k)%values) /= NC_NOERR) return
        end do
        ierr = 0
        message = ''
    end subroutine write_group

    subroutine provenance(ncid, cfg, st)
        integer(c_int), intent(in) :: ncid
        type(synth_config_t), intent(in) :: cfg
        type(synth_state_t), intent(in) :: st
        integer(c_int) :: s
        s = put_text_att(ncid, NC_GLOBAL, 'title', 'TIAGO synthetic diagnostics')
        s = put_text_att(ncid, NC_GLOBAL, 'tiago_revision', tiago_source_revision)
        s = nc_put_att_double(ncid, NC_GLOBAL, cstr('aminor_p'), NC_DOUBLE, 1_c_size_t, [st%aminor])
        s = put_text_att(ncid, NC_GLOBAL, 'plasma_sheet', cfg%sheet)
        call input('wout', cfg%wout)
        call input('profiles', cfg%profiles)
        call input('geometry', cfg%geometry)
        call input('sensors', cfg%sensors)
        call input('coils', cfg%coils)
        call input('extcur', cfg%extcur)
    contains
        subroutine input(name, path)
            character(len=*), intent(in) :: name, path
            if (len(path) == 0) return
            s = put_text_att(ncid, NC_GLOBAL, name//'_path', path)
            s = put_text_att(ncid, NC_GLOBAL, name//'_sha256', sha256_file(path))
        end subroutine input
    end subroutine provenance

    subroutine read_aminor(path, aminor, ierr, message)
        character(len=*), intent(in) :: path
        real(dp), intent(out) :: aminor
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        integer :: ncid, varid, s
        ierr = 1
        message = path//': no Aminor_p'
        aminor = 0.0_dp
        if (nf90_open(path, nf90_nowrite, ncid) /= nf90_noerr) return
        if (nf90_inq_varid(ncid, 'Aminor_p', varid) == nf90_noerr) then
            if (nf90_get_var(ncid, varid, aminor) == nf90_noerr) ierr = 0
        end if
        s = nf90_close(ncid)
        if (ierr == 0) message = ''
    end subroutine read_aminor
end module tiago_synthetic_forward
