module joint_recorded_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use tiago_joint_model, only: joint_model_t
    use tiago_reconstruction, only: rec, invalidate_caches
    use tiago_synthetic_case, only: setup_magnetics, synthesize_magnetics
    use tiago_kinetic_profiles, only: P_NE, P_TE, P_TI, P_ZEFF
    use tiago_calibration, only: calibration_prior_t
    use tiago_thomson_scattering, only: thomson_volume_t
    use tiago_thomson_polychromator, only: recorded_thomson_t, polychromator_channel_t
    use tiago_dispersion_interferometer, only: dispersion_interferometer_t
    use tiago_sight_geometry, only: chord_t
    use tiago_ece_radiometer, only: ece_t, ece_channel_t
    use tiago_ece_receiver, only: ece_response_channel_t
    use tiago_profile_diagnostic, only: plasma_sample_t
    use tiago_real_space_state, only: located_t
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_rng, only: rng_t, rng_seed
    use fortnum_status, only: fortnum_status_t
    implicit none
    integer, parameter::nt = 9, ne = 4
contains
    subroutine configure_joint(model, input, sensors, work, truth)
        type(joint_model_t), intent(out)::model
        character(*), intent(in)::input, sensors, work
        real(dp), allocatable, intent(out)::truth(:)
        type(recorded_thomson_t)::ts
        type(dispersion_interferometer_t)::di
        type(ece_t)::ece
        type(calibration_prior_t)::prior, walls
        type(thomson_volume_t)::volumes(3)
        type(polychromator_channel_t)::bands(nt)
        type(chord_t)::chords(3)
        type(ece_channel_t)::rays(ne)
        type(ece_response_channel_t)::receiver(ne)
        type(rng_t)::rng
        type(fortnum_status_t)::status
        type(plasma_sample_t)::sample
        type(located_t), allocatable::locations(:)
        real(dp), allocatable::mean(:), sigma(:), r(:), eta(:), dx(:, :, :), y(:)
        real(dp)::radius, z, d(2, 3), x(3), xu(3, 3), b, db(3), wave(32), weights(32), &
                   lo, hi, c(3), noise_floor
        real(dp), parameter::pi = acos(-1._dp), laser = 532.e-9_dp
        integer::i, j, k, code, row, ix
        logical::ok
        call setup_magnetics(input, work, sensors, [character(len=64)::'PHIEDGE'], &
                             [.514386_dp], [.01_dp], 1.e-18_dp)
        call rec%eq%vmec%set_input('hot_restart', 0._dp)
        ! The shared integrator records flux in Wb on all channels.
        rec%n_meas = 3*size(rec%loops); rec%n_diag_probes = 0
        rec%meas_kind = [(1, i=1, rec%n_meas)]
        rec%meas_index = [(mod(i - 1, size(rec%loops)) + 1, i=1, rec%n_meas)]
        rec%meas_label = [character(len=64):: &
                          (rec%loops(mod(i - 1, size(rec%loops)) + 1)%label, i=1, rec%n_meas)]
        rec%meas_value = [(0._dp, i=1, rec%n_meas)]
        rec%meas_sigma = [(1._dp, i=1, rec%n_meas)]
        deallocate (rec%probes, rec%segs, rec%probe_m, rec%seg_m)
        allocate (rec%probes(0), rec%segs(0), rec%probe_m(0, 0), rec%seg_m(0, 0))
        call model%prof%init([0._dp, .4_dp, 1._dp])
        associate (g => model%prof%greville)
            call model%prof%set_profile(P_NE, 2.e18_dp*exp(-.5_dp*g), .true., &
                                        .15_dp, .5_dp, .3_dp, .6_dp, .2_dp)
            call model%prof%set_profile(P_TE, 1200*exp(-.6_dp*g), .true., &
                                        .15_dp, .5_dp, .3_dp, .6_dp, .2_dp)
            call model%prof%set_profile(P_TI, 900*exp(-.6_dp*g), .false.)
            call model%prof%set_profile(P_ZEFF, 1.5_dp + 0*g, .false.)
        end associate
        allocate (model%diags(0))
        call model%setup()
        allocate (truth(model%n)); truth = 0; truth(1) = rec%prior(1)/rec%scale(1)
        allocate (r(model%n_residuals()))
        call model%residual(truth, r, ok)
        call require(ok, 'joint initial equilibrium')
        do i = 1, 3
            call model%state%point(.2_dp + .2_dp*i, 0._dp, 0._dp, radius, z, d)
            volumes(i)%label = 'synthetic'
            volumes(i)%center = [radius, 0._dp, z]
            volumes(i)%length = .003_dp
            volumes(i)%npoints = 3
            chords(i) = chord_t('synthetic', [radius, 0._dp, -.18_dp], &
                                [0._dp, 0._dp, 1._dp], .36_dp, 8, 4)
        end do
        allocate (mean(3 + 2*nt), sigma(3 + 2*nt))
        mean = 1; mean(3) = .8_dp; sigma = .002_dp
        sigma(1) = .2_dp
        do i = 1, nt
            j = mod(i - 1, 3) + 1
            select case (j)
            case (1); lo = .82_dp; hi = .92_dp
            case (2); lo = .92_dp; hi = .97_dp
            case default; lo = .97_dp; hi = 1.08_dp
            end select
            call gauss_legendre_ab(32, lo*laser, hi*laser, wave, weights)
            bands(i)%name = 'synthetic'
            bands(i)%volume = (i - 1)/3 + 1
            bands(i)%lambda0 = laser; bands(i)%theta = pi/2
            bands(i)%energy = .5_dp; bands(i)%solid_angle = .01_dp
            bands(i)%wavelength = wave; bands(i)%weights = weights
            bands(i)%transmission = spread(1._dp, 1, 32)
            bands(i)%qe = spread(.6_dp, 1, 32)
            mean(5 + 2*(i - 1)) = 1.e-10_dp
            sigma(5 + 2*(i - 1)) = .002e-10_dp
        end do
        call make_prior(prior, mean, sigma)
        call ts%init(volumes, bands, prior, code, aligned=.true.)
        call require(code == 0, 'joint Thomson setup')
        deallocate (mean, sigma)
        mean = [1.e-5_dp, 1._dp, 1._dp, 0._dp, 0._dp, 0._dp, 1._dp, 1._dp]
        sigma = [1.e-8_dp, .003_dp, .003_dp, .003_dp, .003_dp, .2_dp, .002_dp, .002_dp]
        call make_prior(prior, mean, sigma)
        call di%init(chords, prior, [0._dp, 0._dp, 0._dp], [0._dp, 0._dp, 0._dp], code)
        call require(code == 0, 'joint interferometer setup')
        call model%state%point(.85_dp, 0._dp, 0._dp, radius, z, d)
        call model%state%eval([.60_dp, 0._dp, 0._dp], x, xu, b, db)
        deallocate (mean, sigma)
        allocate (mean(3*ne), sigma(3*ne))
        do i = 1, ne
            rays(i)%label = 'synthetic'
            rays(i)%origin = [radius, 0._dp, z]
            rays(i)%direction = [-1._dp, 0._dp, 0._dp]
            rays(i)%length = .08_dp
            rays(i)%nscan = 16
            rays(i)%frequency = 2*1.7588200107721631e11_dp*b/(2*pi)* &
                                (.97_dp + .02_dp*(i - 1))
            receiver(i)%rays = reshape([i], [1, 1])
            receiver(i)%frequency = [rays(i)%frequency]
            receiver(i)%weights = [1.e8_dp]
            receiver(i)%response = [1._dp]; receiver(i)%antenna = [1._dp]
            receiver(i)%duration = .01_dp
            mean(3*i - 2:3*i) = [.01_dp, .02_dp, 20._dp]
            sigma(3*i - 2:3*i) = [.001_dp, .0001_dp, .01_dp]
        end do
        call make_prior(prior, mean, sigma)
        deallocate (mean, sigma)
        allocate (mean(2*ne), sigma(2*ne))
        do i = 1, ne
            mean(2*i - 1:2*i) = [.2_dp, 10._dp]
            sigma(2*i - 1:2*i) = [.001_dp, .01_dp]
        end do
        call make_prior(walls, mean, sigma)
        call ece%init(rays, receiver, prior, walls, code)
        call require(code == 0, 'joint ECE setup')
        deallocate (model%diags, truth, r)
        allocate (model%diags(3))
        allocate (model%diags(1)%d, source=ts)
        allocate (model%diags(2)%d, source=di)
        allocate (model%diags(3)%d, source=ece)
        do i = 1, 3
            j = model%diags(i)%d%n_signals()
            allocate (model%diags(i)%data(j), model%diags(i)%covariance(j, j))
            model%diags(i)%data = 0
            model%diags(i)%covariance = 0
            do k = 1, j
                model%diags(i)%covariance(k, k) = 1
            end do
        end do
        call make_prior(model%mag_calibration, [1._dp, 0._dp, 0._dp], &
                        [.01_dp, 1.e-5_dp, 1.e-5_dp])
        model%mag_time = [(real((i - 1)/size(rec%loops), dp), i=1, rec%n_meas)]
        call model%setup()
        ! A common ECE gain error, with its own standard normal prior.
        row = model%diags(3)%eta_offset - model%n_rec - model%n_prof
        do i = 2, ne
            model%eta_factor(row + 3*i - 2, :) = model%eta_factor(row + 1, :)
        end do
        allocate (truth(model%n)); truth = 0; truth(1) = rec%prior(1)/rec%scale(1) + .002_dp
        truth(model%n_rec + 1:model%n_rec + model%n_prof) = .15_dp
        truth(model%diags(1)%eta_offset + 1) = 1._dp
        truth(model%diags(2)%eta_offset + 6) = .7_dp
        truth(model%diags(3)%eta_offset + 1) = 1.2_dp
        truth(model%mag_calibration_index + 3) = 1.5_dp
        allocate (r(model%n_residuals()))
        call model%residual(truth, r, ok)
        call require(ok, 'joint recorded truth')
        call rng_seed(rng, 20261028_int64, status)
        call synthesize_magnetics(model, truth(:rec%n)*rec%scale, &
                                  .002_dp, 1.e-7_dp, 1.e-7_dp, 1.e-7_dp, rng, .false.)
        ix = model%mag_calibration_index
        c = model%mag_calibration%physical(truth(ix + 1:ix + 3))
        rec%meas_value = c(1)*rec%meas_value + c(2) + c(3)*model%mag_time
        do i = 1, 3
            j = model%diags(i)%d%n_signals()
            allocate (y(j))
            call model%predict_diagnostic(truth, i, y, sample, locations, eta, dx)
            model%diags(i)%data = y
            model%diags(i)%covariance = 0
            noise_floor = 1.e-5_dp
            if (i == 1) noise_floor = 1.e-16_dp
            do k = 1, j
                model%diags(i)%covariance(k, k) = (.005_dp*abs(y(k)) + noise_floor)**2
            end do
            deallocate (y)
        end do
        call invalidate_caches()
        call model%update_covariances()
    end subroutine
    subroutine make_prior(prior, mean, sigma)
        type(calibration_prior_t), intent(out)::prior
        real(dp), intent(in)::mean(:), sigma(:)
        real(dp)::cov(size(mean), size(mean))
        integer::i, code
        cov = 0
        do i = 1, size(mean)
            cov(i, i) = sigma(i)**2
        end do
        call prior%init(mean, cov, code)
        call require(code == 0, 'synthetic calibration prior')
    end subroutine
    subroutine require(ok, message)
        logical, intent(in)::ok
        character(*), intent(in)::message
        if (.not. ok) then
            print *, message
            error stop 1
        end if
    end subroutine
end module
