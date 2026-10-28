program test_synthetic_forward
    !! End-to-end synthetic forward driver: geometry and profile tables read
    !! from netCDF-4 files, every instrument evaluated, results read back
    !! from the output groups.
    !! Usage: test_synthetic_forward <workdir>             (circular state)
    !!        test_synthetic_forward <workdir> <wout.nc>   (VMEC nested state)
    !! Oracles: profile tables whose logarithm is cubic in s, so the spline
    !! fit reproduces them; on concentric circles (s = r^2, |B| = B0 R0 / R)
    !! independent Simpson quadratures of int ne dl and of the XICS
    !! emissivity moments and the closed-form cold resonance R = B0 R0 / B_res;
    !! on the VMEC state Thomson points placed at X(u) of known labels u.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit, int8
    use, intrinsic :: iso_c_binding, only: c_int, c_size_t
    use tiago_netcdf_c, only: NC_NOERR, NC_NOWRITE, NC_CLOBBER, NC_NETCDF4, NC_CHAR, NC_DOUBLE, &
        nc_open, nc_create, nc_close, nc_inq_ncid, nc_def_grp, nc_def_dim, nc_def_var, &
        nc_put_var_double, nc_put_var_text, nc_inq_varid, nc_get_var_double, cstr
    use tiago_analytic_state, only: circular_state_t
    use tiago_synthetic_forward, only: synth_config_t, synth_state_t, synth_load, synth_write, &
        set_profile_table
    use tiago_profile_geometry_io, only: profile_table_t
    use tiago_sha256, only: sha256_bytes
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: n0 = 5.0e19_dp, t0 = 3.0e3_dp, ti0 = 1.5e3_dp, zeff0 = 2.0_dp
    real(dp), parameter :: e_over_me = 1.75882001076e11_dp
    integer :: failures
    character(len=1024) :: workdir, wout

    failures = 0
    call get_command_argument(1, workdir)
    call sha256_vectors()
    if (command_argument_count() == 1) then
        call analytic_case()
    else
        call get_command_argument(2, wout)
        call vmec_case()
    end if
    if (failures > 0) then
        write(error_unit, '(I0,A)') failures, ' check(s) failed'
        error stop 1
    end if
    print '(A)', 'test_synthetic_forward passed'

contains

    pure real(dp) function ne_of(s)
        real(dp), intent(in) :: s
        ne_of = n0 * exp(-0.4_dp * s - 1.1_dp * s**2 + 0.3_dp * s**3)
    end function ne_of

    pure real(dp) function te_of(s)
        real(dp), intent(in) :: s
        te_of = t0 * exp(-1.0_dp * s - 0.9_dp * s**2 + 0.2_dp * s**3)
    end function te_of

    pure real(dp) function ti_of(s)
        real(dp), intent(in) :: s
        ti_of = ti0 * exp(-0.8_dp * s - 0.5_dp * s**2 + 0.1_dp * s**3)
    end function ti_of

    subroutine sha256_vectors()
        !! FIPS 180-4 examples.
        integer(int8) :: abc(3)
        integer(int8) :: empty(0)
        abc = int([97, 98, 99], int8)
        call same('SHA-256 of ""', sha256_bytes(empty), &
            'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855')
        call same('SHA-256 of "abc"', sha256_bytes(abc), &
            'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad')
    end subroutine sha256_vectors

    subroutine analytic_case()
        type(synth_config_t) :: cfg
        type(synth_state_t) :: st
        real(dp) :: u0, h0, chord(3, 2, 2), ts(3, 6), y(8), b_res, f, x0(3)
        real(dp), allocatable :: v(:)
        integer :: k, ierr
        character(len=:), allocatable :: message

        allocate(circular_state_t :: st%state)
        call st%profiles%init([(real(k, dp) / 12, k = 0, 12)])
        allocate(st%z(0))
        call set_profile_table(st, table(), ierr, message)
        call check(ierr == 0, 'profile table accepted')
        st%aminor = 1.0_dp

        ! Chords between their boundary crossings: vertical at R = 3.4,
        ! horizontal at Z = 0.3.
        u0 = sqrt(1.0_dp - 0.16_dp)
        h0 = sqrt(1.0_dp - 0.09_dp)
        chord(:, 1, 1) = [3.4_dp, 0.0_dp, -u0]
        chord(:, 2, 1) = [3.4_dp, 0.0_dp, u0]
        chord(:, 1, 2) = [3.0_dp - h0, 0.0_dp, 0.3_dp]
        chord(:, 2, 2) = [3.0_dp + h0, 0.0_dp, 0.3_dp]
        ts(:, 1) = [3.2_dp, 0.0_dp, 0.1_dp]
        ts(:, 2) = [0.0_dp, 2.6_dp, -0.3_dp]
        ts(:, 3) = [3.0_dp * cos(1.0_dp), 3.0_dp * sin(1.0_dp), 0.5_dp]
        ts(:, 4) = [3.7_dp, 0.0_dp, 0.0_dp]
        ts(:, 5) = [2.15_dp, 0.0_dp, 0.0_dp]
        ts(:, 6) = [3.3_dp, 0.0_dp, 0.2_dp] ! finite volume along Z
        b_res = 2.5_dp * 3.0_dp / 3.2_dp
        f = 2.0_dp * e_over_me * b_res / (2.0_dp * pi)
        x0 = [4.5_dp, 0.0_dp, 0.1_dp]
        call write_geometry(trim(workdir)//'/analytic_geometry.nc', chord, ts, &
            [f], x0, 2.5_dp)

        cfg%wout = ''
        cfg%output = trim(workdir)//'/analytic_synthetic.nc'
        cfg%profiles = ''
        cfg%geometry = trim(workdir)//'/analytic_geometry.nc'
        cfg%sensors = ''
        cfg%coils = ''
        cfg%extcur = ''
        cfg%sheet = 'contravariant'
        cfg%npanel = 200
        cfg%ece_points = 4001
        call synth_write(cfg, st, ierr, message)
        call check(ierr == 0, 'synth_write on the circular state: '//message)
        if (ierr /= 0) return

        call read_var(cfg%output, 'interferometer', 'line_density', v)
        do k = 1, 2
            y(k) = simpson(chord(:, :, k), 1)
        end do
        call close('interferometer int ne dl, vertical', v(1), y(1), 1.0e-12_dp)
        call close('interferometer int ne dl, horizontal', v(2), y(2), 1.0e-12_dp)

        call read_var(cfg%output, 'xics', 'brightness', v)
        do k = 1, 2
            y(k) = simpson(chord(:, :, k), 2)
            y(2 + k) = simpson(chord(:, :, k), 3) / y(k)
        end do
        call close('XICS brightness, vertical', v(1), y(1), 1.0e-11_dp)
        call close('XICS brightness, horizontal', v(2), y(2), 1.0e-11_dp)
        call read_var(cfg%output, 'xics', 'ti', v)
        call close('XICS weighted Ti, vertical', v(1), y(3), 1.0e-11_dp)
        call close('XICS weighted Ti, horizontal', v(2), y(4), 1.0e-11_dp)

        call read_var(cfg%output, 'thomson', 'ne', v)
        do k = 1, 5
            call close('Thomson point ne', v(k), ne_of(r2(ts(:, k))), 1.0e-12_dp)
        end do
        call close('Thomson volume ne', v(6), volume_mean(ts(:, 6), 1), 1.0e-11_dp)
        call read_var(cfg%output, 'thomson', 'te', v)
        do k = 1, 5
            call close('Thomson point Te', v(k), te_of(r2(ts(:, k))), 1.0e-12_dp)
        end do
        call close('Thomson volume Te', v(6), volume_mean(ts(:, 6), 2) / &
            volume_mean(ts(:, 6), 1), 1.0e-11_dp)
        call read_var(cfg%output, 'thomson', 's', v)
        do k = 1, 6
            call close('Thomson s', v(k), r2(ts(:, k)), 1.0e-12_dp)
        end do

        call read_var(cfg%output, 'ece', 'R', v)
        call close('ECE cold resonance R', v(1), 3.2_dp, 1.0e-7_dp)
        call read_var(cfg%output, 'ece', 'te', v)
        call close('ECE Te at the resonance', v(1), te_of(0.2_dp**2 + 0.1_dp**2), 1.0e-6_dp)
    end subroutine analytic_case

    function table() result(t)
        type(profile_table_t) :: t
        integer :: k
        t%s = [(real(k, dp) / 30, k = 0, 30)]
        allocate(t%ne(31), t%te(31), t%ti(31), t%zeff(31))
        do k = 1, 31
            t%ne(k) = ne_of(t%s(k))
            t%te(k) = te_of(t%s(k))
            t%ti(k) = ti_of(t%s(k))
        end do
        t%zeff = zeff0
    end function table

    pure real(dp) function r2(x)
        !! s = r^2 on the circular state (R0 = 3, a = 1).
        real(dp), intent(in) :: x(3)
        r2 = (hypot(x(1), x(2)) - 3.0_dp)**2 + x(3)**2
    end function r2

    pure real(dp) function integrand(s, which)
        !! 1: ne; 2: XICS emissivity (Z = 18, E = 3.1 keV); 3: emissivity Ti.
        real(dp), intent(in) :: s
        integer, intent(in) :: which
        real(dp) :: eps
        select case (which)
        case (1)
            integrand = ne_of(s)
        case default
            eps = ne_of(s)**2 * (zeff0 - 1.0_dp) * exp(-3.1e3_dp / te_of(s)) &
                / (18.0_dp * 17.0_dp * sqrt(te_of(s)))
            integrand = eps
            if (which == 3) integrand = eps * ti_of(s)
        end select
    end function integrand

    real(dp) function simpson(ends, which)
        real(dp), intent(in) :: ends(3, 2)
        integer, intent(in) :: which
        integer, parameter :: m = 20000
        real(dp) :: h, x(3)
        integer :: j
        h = norm2(ends(:, 2) - ends(:, 1)) / m
        simpson = 0.0_dp
        do j = 0, m
            x = ends(:, 1) + (ends(:, 2) - ends(:, 1)) * (real(j, dp) / m)
            simpson = simpson + merge(1.0_dp, merge(4.0_dp, 2.0_dp, mod(j, 2) == 1), &
                j == 0 .or. j == m) * integrand(min(r2(x), 1.0_dp), which)
        end do
        simpson = simpson * h / 3.0_dp
    end function simpson

    real(dp) function volume_mean(c, which)
        !! Mean over Z in [c - 0.05, c + 0.05] of ne (1) or ne Te (2).
        real(dp), intent(in) :: c(3)
        integer, intent(in) :: which
        integer, parameter :: m = 2000
        real(dp) :: s
        integer :: j
        volume_mean = 0.0_dp
        do j = 0, m
            s = r2(c + [0.0_dp, 0.0_dp, 0.1_dp * (real(j, dp) / m - 0.5_dp)])
            volume_mean = volume_mean + merge(1.0_dp, merge(4.0_dp, 2.0_dp, mod(j, 2) == 1), &
                j == 0 .or. j == m) * merge(ne_of(s), ne_of(s) * te_of(s), which == 1)
        end do
        volume_mean = volume_mean / (3.0_dp * m)
    end function volume_mean

    subroutine vmec_case()
        use tiago_vmec_state, only: vmec_state_t
        type(synth_config_t) :: cfg
        type(synth_state_t) :: st
        type(vmec_state_t) :: vm
        real(dp) :: u(3), x(3), jac(3, 3), b, db(3), ts(3, 8), chord(3, 2, 1)
        real(dp), allocatable :: v(:), s(:)
        integer :: k, ierr
        character(len=:), allocatable :: message

        call vm%load(trim(wout), ierr, message)
        call check(ierr == 0, 'load '//trim(wout))
        if (ierr /= 0) return
        do k = 1, 8
            u = [0.1_dp * k, 0.8_dp * k, 0.37_dp * k]
            call vm%eval(u, x, jac, b, db)
            ts(:, k) = x
        end do
        chord(:, 1, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
        chord(:, 2, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
        call write_geometry(trim(workdir)//'/vmec_geometry.nc', chord, ts, [real(dp) ::], &
            [0.0_dp, 0.0_dp, 0.0_dp], 1.0_dp)
        call write_profiles(trim(workdir)//'/vmec_profiles.nc', table())

        cfg%wout = trim(wout)
        cfg%output = trim(workdir)//'/vmec_synthetic.nc'
        cfg%profiles = trim(workdir)//'/vmec_profiles.nc'
        cfg%geometry = trim(workdir)//'/vmec_geometry.nc'
        cfg%sensors = ''
        cfg%coils = ''
        cfg%extcur = ''
        cfg%sheet = 'contravariant'
        call synth_load(cfg, st, ierr, message)
        call check(ierr == 0, 'synth_load: '//message)
        if (ierr /= 0) return
        call synth_write(cfg, st, ierr, message)
        call check(ierr == 0, 'synth_write on the VMEC state: '//message)
        if (ierr /= 0) return
        call read_var(cfg%output, 'thomson', 's', s)
        call read_var(cfg%output, 'thomson', 'ne', v)
        do k = 1, 8
            call close('VMEC Thomson s at X(u)', s(k), (0.1_dp * k)**2, 1.0e-9_dp)
            call close('VMEC Thomson ne at X(u)', v(k), ne_of((0.1_dp * k)**2), 1.0e-9_dp)
        end do
        call read_var(cfg%output, 'thomson', 'te', v)
        do k = 1, 8
            call close('VMEC Thomson Te at X(u)', v(k), te_of((0.1_dp * k)**2), 1.0e-9_dp)
        end do
        call read_var(cfg%output, 'thomson', 'r_eff', v)
        call close('VMEC r_eff = Aminor_p sqrt(s)', v(5), st%aminor * 0.5_dp, 1.0e-9_dp)
    end subroutine vmec_case

    subroutine write_geometry(path, chords, centres, freq, origin, length)
        !! Interferometer and XICS on the same chords, Thomson point volumes
        !! (the sixth, if present, 0.1 m along Z), ECE along -x.
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: chords(:, :, :), centres(:, :), freq(:), origin(3), length
        integer(c_int) :: ncid, s
        integer :: k
        real(dp), allocatable :: dir(:, :), len_(:)

        s = nc_create(cstr(path), ior(NC_CLOBBER, NC_NETCDF4), ncid)
        call check(s == NC_NOERR, 'create '//path)
        call put_group(ncid, 'interferometer', 'chord', &
            [character(len=4) :: ('C'//achar(48 + k), k = 1, size(chords, 3))], &
            ['start', 'end  '], [chords(:, 1, :), chords(:, 2, :)], size(chords, 3))
        call put_group(ncid, 'xics', 'chord', &
            [character(len=4) :: ('X'//achar(48 + k), k = 1, size(chords, 3))], &
            ['start', 'end  '], [chords(:, 1, :), chords(:, 2, :)], size(chords, 3))
        allocate(dir(3, size(centres, 2)), len_(size(centres, 2)))
        dir = 0.0_dp
        dir(3, :) = 1.0_dp
        len_ = 0.0_dp
        if (size(centres, 2) == 6) len_(6) = 0.1_dp
        call put_thomson(ncid, centres, dir, len_)
        if (size(freq) > 0) call put_ece(ncid, freq, origin, length)
        s = nc_close(ncid)
    end subroutine write_geometry

    subroutine put_group(ncid, group, dimname, labels, names, values, n)
        !! (n, 3) matrices names(:) stored consecutively in values.
        integer(c_int), intent(in) :: ncid
        character(len=*), intent(in) :: group, dimname, labels(:), names(:)
        real(dp), intent(in) :: values(:)
        integer, intent(in) :: n
        integer(c_int) :: gid, dn, dx, dl, var, s
        integer :: k
        s = nc_def_grp(ncid, cstr(group), gid)
        s = nc_def_dim(gid, cstr(dimname), int(n, c_size_t), dn)
        s = nc_def_dim(gid, cstr('xyz'), 3_c_size_t, dx)
        s = nc_def_dim(gid, cstr('strlen'), int(len(labels), c_size_t), dl)
        do k = 1, size(names)
            s = nc_def_var(gid, cstr(trim(names(k))), NC_DOUBLE, 2_c_int, [dn, dx], var)
            s = nc_put_var_double(gid, var, values(3 * n * (k - 1) + 1:3 * n * k))
            call check(s == NC_NOERR, 'write '//group//'/'//trim(names(k)))
        end do
        s = nc_def_var(gid, cstr('label'), NC_CHAR, 2_c_int, [dn, dl], var)
        s = nc_put_var_text(gid, var, labels)
    end subroutine put_group

    subroutine put_thomson(ncid, centres, dir, len_)
        integer(c_int), intent(in) :: ncid
        real(dp), intent(in) :: centres(:, :), dir(:, :), len_(:)
        integer(c_int) :: gid, dn, dx, var, s
        s = nc_def_grp(ncid, cstr('thomson'), gid)
        s = nc_def_dim(gid, cstr('volume'), int(size(len_), c_size_t), dn)
        s = nc_def_dim(gid, cstr('xyz'), 3_c_size_t, dx)
        s = nc_def_var(gid, cstr('center'), NC_DOUBLE, 2_c_int, [dn, dx], var)
        s = nc_put_var_double(gid, var, centres)
        s = nc_def_var(gid, cstr('direction'), NC_DOUBLE, 2_c_int, [dn, dx], var)
        s = nc_put_var_double(gid, var, dir)
        s = nc_def_var(gid, cstr('length'), NC_DOUBLE, 1_c_int, [dn], var)
        s = nc_put_var_double(gid, var, len_)
        call check(s == NC_NOERR, 'write thomson')
    end subroutine put_thomson

    subroutine put_ece(ncid, freq, origin, length)
        integer(c_int), intent(in) :: ncid
        real(dp), intent(in) :: freq(:), origin(3), length
        integer(c_int) :: gid, dn, dx, var, s
        integer :: n
        n = size(freq)
        s = nc_def_grp(ncid, cstr('ece'), gid)
        s = nc_def_dim(gid, cstr('channel'), int(n, c_size_t), dn)
        s = nc_def_dim(gid, cstr('xyz'), 3_c_size_t, dx)
        s = nc_def_var(gid, cstr('frequency'), NC_DOUBLE, 1_c_int, [dn], var)
        s = nc_put_var_double(gid, var, freq)
        s = nc_def_var(gid, cstr('harmonic'), NC_DOUBLE, 1_c_int, [dn], var)
        s = nc_put_var_double(gid, var, spread(2.0_dp, 1, n))
        s = nc_def_var(gid, cstr('origin'), NC_DOUBLE, 2_c_int, [dn, dx], var)
        s = nc_put_var_double(gid, var, spread(origin, 2, n))
        s = nc_def_var(gid, cstr('direction'), NC_DOUBLE, 2_c_int, [dn, dx], var)
        s = nc_put_var_double(gid, var, spread([-1.0_dp, 0.0_dp, 0.0_dp], 2, n))
        s = nc_def_var(gid, cstr('length'), NC_DOUBLE, 1_c_int, [dn], var)
        s = nc_put_var_double(gid, var, spread(length, 1, n))
        call check(s == NC_NOERR, 'write ece')
    end subroutine put_ece

    subroutine write_profiles(path, t)
        character(len=*), intent(in) :: path
        type(profile_table_t), intent(in) :: t
        integer(c_int) :: ncid, dn, var, s
        s = nc_create(cstr(path), ior(NC_CLOBBER, NC_NETCDF4), ncid)
        s = nc_def_dim(ncid, cstr('s'), int(size(t%s), c_size_t), dn)
        s = nc_def_var(ncid, cstr('s'), NC_DOUBLE, 1_c_int, [dn], var)
        s = nc_put_var_double(ncid, var, t%s)
        s = nc_def_var(ncid, cstr('ne'), NC_DOUBLE, 1_c_int, [dn], var)
        s = nc_put_var_double(ncid, var, t%ne)
        s = nc_def_var(ncid, cstr('te'), NC_DOUBLE, 1_c_int, [dn], var)
        s = nc_put_var_double(ncid, var, t%te)
        s = nc_def_var(ncid, cstr('ti'), NC_DOUBLE, 1_c_int, [dn], var)
        s = nc_put_var_double(ncid, var, t%ti)
        call check(s == NC_NOERR, 'write '//path)
        s = nc_close(ncid)
    end subroutine write_profiles

    subroutine read_var(path, group, name, v)
        character(len=*), intent(in) :: path, group, name
        real(dp), allocatable, intent(out) :: v(:)
        integer(c_int) :: ncid, gid, var, s
        allocate(v(64))
        v = huge(1.0_dp)
        s = nc_open(cstr(path), NC_NOWRITE, ncid)
        if (s == NC_NOERR) s = nc_inq_ncid(ncid, cstr(group), gid)
        if (s == NC_NOERR) s = nc_inq_varid(gid, cstr(name), var)
        if (s == NC_NOERR) s = nc_get_var_double(gid, var, v)
        call check(s == NC_NOERR, 'read '//group//'/'//name)
        s = nc_close(ncid)
    end subroutine read_var

    subroutine check(ok, name)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: name
        if (.not. ok) then
            write(error_unit, '(A)') 'FAIL '//name
            failures = failures + 1
        end if
    end subroutine check

    subroutine same(name, a, b)
        character(len=*), intent(in) :: name, a, b
        call check(a == b, name//': '//a)
        if (a == b) print '(A)', 'ok   '//name
    end subroutine same

    subroutine close(name, a, b, rtol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: a, b, rtol
        if (abs(a - b) > rtol * abs(b) .or. a /= a) then
            write(error_unit, '(A,2ES24.16)') 'FAIL '//name//': ', a, b
            failures = failures + 1
        else
            print '(A,2ES24.16)', 'ok   '//name//': ', a, b
        end if
    end subroutine close
end program test_synthetic_forward
