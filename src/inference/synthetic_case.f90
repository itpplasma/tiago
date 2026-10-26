module tiago_synthetic_case
    !! Synthetic fixed-boundary joint-inference case for tests and the
    !! coverage study: magnetic diagnostics from DIAGNO files (as written by
    !! benchmarks/reconstruction/sensors.awk), equilibrium/coil parameters with
    !! Gaussian priors, kinetic profiles with a Gaussian-process prior, and
    !! ideal local profile probes (Thomson-like ne/Te points on one toroidal
    !! cut, CXRS-like Ti points on another) placed at given (rho, theta) of a
    !! reference equilibrium. No coils: the plasma signals are the model.
    !! place_instruments sets up the calibration study instead: a Thomson
    !! system with several lasers and a log-scale density calibration, an
    !! interferometer with a shared systematic scale, and optionally the Ti
    !! probes; setup_magnetics can keep only a subset of the magnetic signals.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    use tiago_diagnostic_types, only: bprobe_t
    use tiago_flux_loops, only: read_flux_loop_file
    use tiago_segmented_rogowski, only: read_segmented_rogowski_file
    use tiago_bprobes, only: read_bprobe_file
    use tiago_reconstruction, only: rec, KIND_FLUX, KIND_SEGROG, KIND_BPROBE, model_signals, &
        signal_value
    use tiago_joint_model, only: joint_model_t
    use tiago_point_probe, only: point_probe_t
    use tiago_thomson_scattering, only: thomson_t, thomson_volume_t
    use tiago_line_integrals, only: interferometer_t
    use tiago_sight_geometry, only: chord_t
    use tiago_profile_diagnostic, only: Q_NE, Q_TE, Q_TI, plasma_sample_t
    use tiago_real_space_state, only: located_t
    use fortnum_rng, only: rng_t, rng_normal
    implicit none
    private

    public :: setup_magnetics, place_probes, place_instruments, synthesize_magnetics, &
        synthesize_profiles

contains

    subroutine setup_magnetics(vmec_input, workdir, sensor_dir, names, prior, prior_sigma, ftol, &
            keep)
        !! rec with the given parameters (no coils), the DIAGNO signals in
        !! sensor_dir as measurements (values set by synthesize_magnetics):
        !! keep = 'all' (default), 'diamagnetic' (flux loops labelled DIA_*,
        !! the diamagnetic loops of sensors.awk) or 'none'.
        character(len=*), intent(in) :: vmec_input, workdir, sensor_dir
        character(len=*), intent(in) :: names(:)
        real(dp), intent(in) :: prior(:), prior_sigma(:), ftol
        character(len=*), intent(in), optional :: keep
        logical, allocatable :: mask(:)
        integer(i32) :: ierr
        character(len=:), allocatable :: msg
        type(bprobe_t), allocatable :: measured(:)
        integer :: k, m, n

        n = size(names)
        rec%n = n
        allocate(rec%names(n), rec%is_coil(n), rec%group(n), rec%eq_index(n), rec%x0(n), &
            rec%scale(n), rec%prior(n), rec%prior_sigma(n), rec%history(0))
        do k = 1, n
            rec%names(k) = names(k)
            rec%is_coil(k) = .false.
            rec%group(k) = 0
            rec%eq_index(k) = k
        end do
        call execute_command_line('mkdir -p "'//workdir//'"')
        call rec%eq%init(vmec_input, workdir, names)
        call rec%eq%vmec%set_input('ftol', ftol)
        call rec%eq%vmec%set_input('niter', 20000.0_dp)
        rec%rule%samples_per_segment = 4
        rec%rule%gauss = .true.
        rec%nphi = 32
        rec%ntheta = 32
        call read_flux_loop_file(sensor_dir//'/flux.diagno', rec%loops, ierr, msg)
        if (ierr /= 0) error stop 'flux: '//msg
        call read_segmented_rogowski_file(sensor_dir//'/seg.diagno', rec%segs, ierr, msg)
        if (ierr /= 0) error stop 'segrog: '//msg
        call read_bprobe_file(sensor_dir//'/probes.diagno', measured, ierr, msg, .false.)
        if (ierr /= 0) error stop 'bprobes: '//msg
        rec%n_diag_probes = size(measured)
        rec%n_fb = 0
        rec%probes = measured
        call rec%coils%init('')
        allocate(rec%extcur(0), rec%flux_m(size(rec%loops), 0), rec%seg_m(size(rec%segs), 0), &
            rec%probe_m(size(rec%probes), 0))
        call rec%plasma%init('')
        m = size(rec%loops) + size(rec%segs) + rec%n_diag_probes
        allocate(rec%meas_kind(m), rec%meas_index(m), rec%meas_value(m), rec%meas_sigma(m), &
            rec%meas_label(m))
        m = 0
        do k = 1, size(rec%loops)
            m = m + 1
            rec%meas_kind(m) = KIND_FLUX
            rec%meas_index(m) = k
            rec%meas_label(m) = rec%loops(k)%label
        end do
        do k = 1, size(rec%segs)
            m = m + 1
            rec%meas_kind(m) = KIND_SEGROG
            rec%meas_index(m) = k
            rec%meas_label(m) = rec%segs(k)%label
        end do
        do k = 1, rec%n_diag_probes
            m = m + 1
            rec%meas_kind(m) = KIND_BPROBE
            rec%meas_index(m) = k
            rec%meas_label(m) = rec%probes(k)%label
        end do
        rec%n_meas = m
        allocate(mask(m))
        mask = .true.
        if (present(keep)) then
            select case (keep)
            case ('all')
            case ('diamagnetic')
                mask = rec%meas_kind == KIND_FLUX .and. index(rec%meas_label, 'DIA_') == 1
            case ('none')
                mask = .false.
            case default
                error stop 'setup_magnetics: keep must be all, diamagnetic or none'
            end select
        end if
        rec%meas_kind = pack(rec%meas_kind, mask)
        rec%meas_index = pack(rec%meas_index, mask)
        rec%meas_label = pack(rec%meas_label, mask)
        rec%n_meas = count(mask)
        deallocate(rec%meas_value, rec%meas_sigma)
        allocate(rec%meas_value(rec%n_meas), rec%meas_sigma(rec%n_meas))
        rec%meas_value = 0.0_dp
        rec%meas_sigma = 1.0_dp
        rec%prior = prior
        rec%prior_sigma = prior_sigma
        rec%scale = max(abs(prior), 1.0e-3_dp)
        rec%x0 = prior
        rec%prior_param = pack([(k, k = 1, n)], prior_sigma > 0.0_dp)
        rec%n_prior = size(rec%prior_param)
    end subroutine setup_magnetics

    pure function identity(n) result(a)
        integer, intent(in) :: n
        real(dp) :: a(n, n)
        integer :: i
        a = 0.0_dp
        do i = 1, n
            a(i, i) = 1.0_dp
        end do
    end function identity

    subroutine place_probes(model, rho, theta_cut, phi_ne_te, phi_ti, align_sigma)
        !! Two probe sets at the current equilibrium's (rho, theta) positions:
        !! ne and Te at every rho on cut phi_ne_te, Ti on cut phi_ti, each
        !! position at both theta_cut and theta_cut + pi; synthesize_profiles
        !! sets data and covariances. align_sigma > 0 adds (dR, dZ) nuisances
        !! of the ne/Te set. Call joint setup again afterwards.
        type(joint_model_t), intent(inout) :: model
        real(dp), intent(in) :: rho(:), theta_cut, phi_ne_te, phi_ti, align_sigma
        type(point_probe_t) :: a, b
        real(dp) :: r, z, d(2, 3)
        integer :: i, k, np

        np = 2 * size(rho)
        allocate(a%x0(3, 2 * np), a%quantity(2 * np), b%x0(3, np), b%quantity(np))
        k = 0
        do i = 1, size(rho)
            call model%state%point(rho(i), theta_cut, phi_ne_te, r, z, d)
            a%x0(:, 2 * k + 1) = [r, phi_ne_te, z]
            a%x0(:, 2 * k + 2) = [r, phi_ne_te, z]
            call model%state%point(rho(i), theta_cut, phi_ti, r, z, d)
            b%x0(:, k + 1) = [r, phi_ti, z]
            k = k + 1
            call model%state%point(rho(i), theta_cut + acos(-1.0_dp), phi_ne_te, r, z, d)
            a%x0(:, 2 * k + 1) = [r, phi_ne_te, z]
            a%x0(:, 2 * k + 2) = [r, phi_ne_te, z]
            call model%state%point(rho(i), theta_cut + acos(-1.0_dp), phi_ti, r, z, d)
            b%x0(:, k + 1) = [r, phi_ti, z]
            k = k + 1
        end do
        a%quantity = [([Q_NE, Q_TE], i = 1, np)]
        b%quantity = Q_TI
        a%aligned = align_sigma <= 0.0_dp
        b%aligned = .true.
        if (allocated(model%diags)) deallocate(model%diags)
        allocate(model%diags(2))
        allocate(model%diags(1)%d, source=a)
        allocate(model%diags(2)%d, source=b)
        if (align_sigma > 0.0_dp) then
            model%diags(1)%eta_sigma = [align_sigma, align_sigma]
        else
            allocate(model%diags(1)%eta_sigma(0))
        end if
        allocate(model%diags(2)%eta_sigma(0))
        model%diags(1)%gain_sigma = 0.0_dp
        model%diags(2)%gain_sigma = 0.0_dp
        model%diags(1)%data = [(0.0_dp, i = 1, 2 * np)]
        model%diags(2)%data = [(0.0_dp, i = 1, np)]
        ! unit covariances until synthesize_profiles sets the data
        model%diags(1)%covariance = identity(2 * np)
        model%diags(2)%covariance = identity(np)
        model%diags(1)%model_sigma = [(0.0_dp, i = 1, 2 * np)]
        model%diags(2)%model_sigma = [(0.0_dp, i = 1, np)]
    end subroutine place_probes

    subroutine place_instruments(model, rho, theta_cut, phi_ts, n_laser, cal_mean, cal_sigma, &
            cal_nu, laser_sigma, n_chords, interf_sigma, with_ti, phi_ti)
        !! Calibration study at the current equilibrium. Slot 1: Thomson
        !! volumes (1 cm, vertical laser) at every rho on cut phi_ts, at
        !! theta_cut and theta_cut + pi, probed by each of n_laser lasers,
        !! with the log-scale density calibration (cal_sigma > 0; per-laser
        !! offsets for laser_sigma > 0). Slot 2 (n_chords > 0): horizontal
        !! interferometer chords on the cut pi / (2 nfp) at heights spread over
        !! the plasma, with one shared log-scale systematic of prior sd
        !! interf_sigma (> 0). Last slot (with_ti): the CXRS-like Ti probes of
        !! place_probes on cut phi_ti. synthesize_profiles sets data and
        !! covariances; call joint setup again afterwards.
        type(joint_model_t), intent(inout) :: model
        real(dp), intent(in) :: rho(:), theta_cut, phi_ts, cal_mean, cal_sigma, cal_nu, &
            laser_sigma, interf_sigma, phi_ti
        integer, intent(in) :: n_laser, n_chords
        logical, intent(in) :: with_ti
        type(thomson_t) :: ts
        type(interferometer_t) :: interf
        type(point_probe_t) :: ti
        type(thomson_volume_t), allocatable :: vol(:)
        type(chord_t), allocatable :: chords(:)
        real(dp) :: r, z, d(2, 3), phi, rmin, rmax, zmin, zmax, t
        integer :: i, j, l, k, nd, nv
        character(len=32) :: tag

        nv = 2 * size(rho)
        allocate(vol(nv * n_laser))
        k = 0
        do l = 1, n_laser
            do i = 1, size(rho)
                do j = 0, 1
                    call model%state%point(rho(i), theta_cut + j * acos(-1.0_dp), phi_ts, r, z, d)
                    k = k + 1
                    write(tag, '("L",I0,"_",I0)') l, (i - 1) * 2 + j + 1
                    vol(k)%label = trim(tag)
                    vol(k)%center = [r * cos(phi_ts), r * sin(phi_ts), z]
                    vol(k)%direction = [0.0_dp, 0.0_dp, 1.0_dp]
                    vol(k)%length = 0.01_dp
                    vol(k)%npoints = 3
                    vol(k)%laser = l
                end do
            end do
        end do
        call ts%init(vol, aligned=.true.)
        nd = 1 + merge(1, 0, n_chords > 0) + merge(1, 0, with_ti)
        if (allocated(model%diags)) deallocate(model%diags)
        allocate(model%diags(nd))
        allocate(model%diags(1)%d, source=ts)
        model%diags(1)%cal_sigma = cal_sigma
        model%diags(1)%cal_mean = cal_mean
        model%diags(1)%cal_nu = cal_nu
        model%diags(1)%laser_sigma = laser_sigma
        model%diags(1)%cal_group = ts%calibration_groups()
        k = 1
        if (n_chords > 0) then
            ! bounding box of the boundary on the chord cut
            phi = acos(-1.0_dp) / (2 * rec%eq%nfp)
            rmin = huge(1.0_dp)
            rmax = -huge(1.0_dp)
            zmin = huge(1.0_dp)
            zmax = -huge(1.0_dp)
            do i = 0, 127
                call model%state%point(1.0_dp, i * acos(-1.0_dp) / 64, phi, r, z, d)
                rmin = min(rmin, r)
                rmax = max(rmax, r)
                zmin = min(zmin, z)
                zmax = max(zmax, z)
            end do
            allocate(chords(n_chords))
            do i = 1, n_chords
                t = (i - 0.5_dp) / n_chords
                z = zmin + (0.15_dp + 0.7_dp * t) * (zmax - zmin)
                write(tag, '("chord_",I0)') i
                chords(i)%label = trim(tag)
                r = rmin - 0.1_dp * (rmax - rmin)
                chords(i)%origin = [r * cos(phi), r * sin(phi), z]
                chords(i)%direction = [cos(phi), sin(phi), 0.0_dp]
                chords(i)%length = 1.2_dp * (rmax - rmin)
                chords(i)%npanel = 40
                chords(i)%ngauss = 4
            end do
            call interf%init(chords)
            k = k + 1
            allocate(model%diags(k)%d, source=interf)
            model%diags(k)%cal_sigma = interf_sigma
        end if
        if (with_ti) then
            allocate(ti%x0(3, nv), ti%quantity(nv))
            j = 0
            do i = 1, size(rho)
                call model%state%point(rho(i), theta_cut, phi_ti, r, z, d)
                ti%x0(:, j + 1) = [r, phi_ti, z]
                call model%state%point(rho(i), theta_cut + acos(-1.0_dp), phi_ti, r, z, d)
                ti%x0(:, j + 2) = [r, phi_ti, z]
                j = j + 2
            end do
            ti%quantity = Q_TI
            ti%aligned = .true.
            k = k + 1
            allocate(model%diags(k)%d, source=ti)
        end if
        do k = 1, nd
            associate (slot => model%diags(k))
                allocate(slot%eta_sigma(0))
                slot%gain_sigma = 0.0_dp
                slot%data = [(0.0_dp, i = 1, slot%d%n_signals())]
                slot%covariance = identity(slot%d%n_signals())
                slot%model_sigma = [(0.0_dp, i = 1, slot%d%n_signals())]
            end associate
        end do
    end subroutine place_instruments

    subroutine synthesize_magnetics(model, x_true, sigma_rel, sigma_flux, sigma_segrog, sigma_bprobe, &
            rng, noise)
        !! rec measurements and sigmas from the current equilibrium (solved at
        !! the truth by the caller) at parameters x_true, with Gaussian noise;
        !! resets the model's diagonal magnetic covariance to these sigmas
        !! (joint setup factors whatever rec%meas_sigma held before).
        type(joint_model_t), intent(inout) :: model
        real(dp), intent(in) :: x_true(:), sigma_rel, sigma_flux, sigma_segrog, sigma_bprobe
        type(rng_t), intent(inout) :: rng
        logical, intent(in) :: noise
        real(dp), allocatable :: flux(:), seg(:), probe(:)
        real(dp) :: exact, floor, e
        integer :: i

        call model_signals(x_true, flux, seg, probe)
        do i = 1, rec%n_meas
            exact = signal_value(rec%meas_kind(i), rec%meas_index(i), flux, seg, probe)
            select case (rec%meas_kind(i))
            case (KIND_FLUX)
                floor = sigma_flux
            case (KIND_SEGROG)
                floor = sigma_segrog
            case default
                floor = sigma_bprobe
            end select
            rec%meas_sigma(i) = sigma_rel * abs(exact) + floor
            call rng_normal(rng, e)
            rec%meas_value(i) = exact
            if (noise) rec%meas_value(i) = exact + rec%meas_sigma(i) * e
        end do
        if (allocated(model%mag_covariance)) deallocate(model%mag_covariance)
        allocate(model%mag_covariance(rec%n_meas, rec%n_meas))
        model%mag_covariance = 0.0_dp
        do i = 1, rec%n_meas
            model%mag_covariance(i, i) = rec%meas_sigma(i)**2
        end do
        if (allocated(model%mag_model_sigma)) call model%update_covariances()
    end subroutine synthesize_magnetics

    subroutine synthesize_profiles(model, theta_true, sigma_rel, sigma_floor, rng, noise, rel)
        !! Profile-probe data at theta_true (including the gains and
        !! calibration scales of theta_true): diagonal covariance with sigma =
        !! sigma_rel |y| + sigma_floor(d) (rel(d) instead of sigma_rel if
        !! given), Gaussian noise; resets model_sigma.
        type(joint_model_t), intent(inout) :: model
        real(dp), intent(in) :: theta_true(:), sigma_rel, sigma_floor(:)
        type(rng_t), intent(inout) :: rng
        logical, intent(in) :: noise
        real(dp), intent(in), optional :: rel(:)
        real(dp), allocatable :: y(:), eta(:), dx(:, :, :), r(:)
        type(plasma_sample_t) :: sample
        type(located_t), allocatable :: loc(:)
        real(dp) :: e, gain, sr
        integer :: d, i
        logical :: ok

        allocate(r(model%n_residuals()))
        call model%residual(theta_true, r, ok)
        if (.not. ok) error stop 'synthesize_profiles: no equilibrium at the truth'
        do d = 1, size(model%diags)
            allocate(y(model%diags(d)%d%n_signals()))
            call model%predict_diagnostic(theta_true, d, y, sample, loc, eta, dx)
            call model%unpack(theta_true, gain=gain, d=d)
            y = gain * model%signal_scale(theta_true, d) * y
            sr = sigma_rel
            if (present(rel)) sr = rel(d)
            model%diags(d)%covariance = 0.0_dp * spread(y, 1, size(y)) * spread(y, 2, size(y))
            do i = 1, size(y)
                model%diags(d)%covariance(i, i) = (sr * abs(y(i)) + sigma_floor(d))**2
                call rng_normal(rng, e)
                model%diags(d)%data(i) = y(i)
                if (noise) model%diags(d)%data(i) = y(i) + sqrt(model%diags(d)%covariance(i, i)) * e
            end do
            model%diags(d)%model_sigma = 0.0_dp
            deallocate(y)
        end do
        call model%update_covariances()
    end subroutine synthesize_profiles

end module tiago_synthetic_case
