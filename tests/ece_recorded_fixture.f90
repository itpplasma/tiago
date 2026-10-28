module ece_recorded_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use ece_geometry_fixture, only: setup
    use tiago_ece_state_medium, only: state_medium_t
    use tiago_ece_radiometer, only: ece_t, ece_channel_t
    use tiago_ece_receiver, only: ece_response_channel_t
    use tiago_vmec_state, only: vmec_state_t
    use tiago_kinetic_profiles, only: P_NE, P_TE
    use tiago_calibration, only: calibration_prior_t
    implicit none
contains
    subroutine configure(medium, diagnostic, nray, ncell, scale, varying)
        type(state_medium_t), intent(out) :: medium
        type(ece_t), intent(out) :: diagnostic
        integer, intent(in) :: nray, ncell
        real(dp), intent(in) :: scale
        logical, intent(in) :: varying
        type(ece_channel_t) :: channels(nray)
        type(ece_response_channel_t) :: response(nray)
        type(calibration_prior_t) :: prior, walls
        real(dp) :: mean(3*nray), covariance(3*nray, 3*nray), &
            wm(2*nray), wc(2*nray, 2*nray)
        integer :: c, code
        call setup(medium)
        select type (s => medium%state)
        type is (vmec_state_t)
            s%rmnc = scale*s%rmnc
            s%zmns = scale*s%zmns
            s%bsupvmnc = s%bsupvmnc/scale
        end select
        call medium%profiles%set_profile(P_NE, &
            spread(3.0e17_dp, 1, medium%profiles%ncoef), &
            .true., &
            sigma_f=0.1_dp, l_core=0.4_dp, l_edge=0.2_dp, s0=0.5_dp, &
            sw=0.2_dp)
        if (varying) then
            call medium%profiles%set_profile(P_TE, &
                1000*exp(-.5_dp*medium%profiles%greville), &
                .true., &
                sigma_f=0.1_dp, l_core=0.4_dp, l_edge=0.2_dp, s0=0.5_dp, &
                sw=0.2_dp)
        else
            call medium%profiles%set_profile(P_TE, &
                spread(1000.0_dp, 1, medium%profiles%ncoef), &
                .true., &
                sigma_f=0.1_dp, l_core=0.4_dp, l_edge=0.2_dp, s0=0.5_dp, &
                sw=0.2_dp)
        end if
        deallocate (medium%z)
        allocate (medium%z(medium%profiles%n_latent))
        medium%z = 0
        covariance = 0
        wc = 0
        do c = 1, nray
            channels(c)%label = 'synthetic'
            channels(c)%frequency = 2*1.7588200107721631e11_dp* &
                (.3_dp*(3.58_dp + .02_dp*real(c - 1, dp)/max(1, &
                nray - 1)))/(2*acos(-1.0_dp))
            channels(c)%origin = scale*[3.65_dp, 0.0_dp, 0.1_dp]
            channels(c)%direction = [-1.0_dp, 0.0_dp, 0.0_dp]
            channels(c)%length = scale*.15_dp
            channels(c)%nscan = ncell
            response(c)%rays = reshape([c], [1, 1])
            response(c)%frequency = [channels(c)%frequency]
            response(c)%weights = [1.0e8_dp]
            response(c)%response = [1.0_dp]
            response(c)%antenna = [1.0_dp]
            response(c)%duration = .01_dp
            mean(3*c - 2:3*c) = [0.01_dp, 0.02_dp, 20.0_dp]
            covariance(3*c - 2, 3*c - 2) = 1.0e-6_dp
            covariance(3*c - 1, 3*c - 1) = 1.0e-4_dp
            covariance(3*c, 3*c) = 1.0_dp
            wm(2*c - 1:2*c) = [.2_dp, 10.0_dp]
            wc(2*c - 1, 2*c - 1) = .01_dp**2
            wc(2*c, 2*c) = 1
        end do
        call prior%init(mean, covariance, code)
        if (code /= 0) error stop 'synthetic receiver prior'
        call walls%init(wm, wc, code)
        if (code /= 0) error stop 'synthetic wall prior'
        call diagnostic%init(channels, response, prior, walls, code)
        if (code /= 0) error stop 'synthetic ECE configuration'
    end subroutine
end module
