module ece_geometry_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_vmec_state, only: vmec_state_t
    use tiago_ece_state_medium, only: state_medium_t
    use tiago_kinetic_profiles, only: P_NE
    use tiago_vmecpp_state, only: vmecpp_state_t
    implicit none
contains
    subroutine setup(medium)
        type(state_medium_t), intent(out) :: medium
        integer :: j
        allocate (vmec_state_t :: medium%state)
        select type (s => medium%state)
        type is (vmec_state_t)
            s%ns = 5
            s%mnmax = 2
            s%mnmax_nyq = 2
            s%xm = [0.0_dp, 1.0_dp]
            s%xn = [0.0_dp, 0.0_dp]
            s%xm_nyq = s%xm
            s%xn_nyq = s%xn
            s%rho_full = [(real(j - 1, dp)/4, j=1, 5)]
            s%rho_half = [(real(j, dp)/5, j=1, 4)]
            allocate (s%rmnc(2, 5), s%zmns(2, 5), s%bmnc(2, 4), &
                      s%bsupumnc(2, 4), s%bsupvmnc(2, 4))
            s%rmnc(1, :) = 3
            s%rmnc(2, :) = s%rho_full
            s%zmns = 0
            s%zmns(2, :) = s%rho_full
            s%bmnc(1, :) = 0.9_dp
            s%bmnc(2, :) = 0.3_dp*s%rho_half
            s%bsupumnc = 0
            s%bsupvmnc = 0
            s%bsupvmnc(1, :) = 0.3_dp
        end select
        call medium%profiles%init([0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp])
        call medium%profiles%set_profile(P_NE, &
                                         1.0e19_dp*exp(-0.4_dp*medium%profiles%greville), .true., &
                                         sigma_f=0.1_dp, l_core=0.4_dp, l_edge=0.2_dp, &
                                         s0=0.5_dp, sw=0.2_dp)
        allocate (medium%z(medium%profiles%n_latent))
        medium%z = 0
    end subroutine
    subroutine setup_vmecpp(medium)
        type(state_medium_t), intent(inout)::medium
        type(vmecpp_state_t)::state
        real(dp)::coef(0:1, 0:1, 5, 12), bmnc(2, 5), bu(2, 5), bv(2, 5)
        integer::j
        coef = 0
        bmnc = 0
        bu = 0
        bv = 0
        do j = 1, 5
            coef(0, 0, j, 1) = 3
            coef(0, 1, j, 1) = sqrt(real(j - 1, dp)/4)
            coef(0, 1, j, 5) = sqrt(real(j - 1, dp)/4)
            bmnc(1, j) = 0.9_dp
            if (j > 1) bmnc(2, j) = 0.3_dp*sqrt((j - 1.5_dp)/4)
            bv(1, j) = 0.3_dp
        end do
        call state%init(coef, 1, [0.0_dp, 1.0_dp], [0.0_dp, 0.0_dp], bmnc, bu, bv)
        deallocate (medium%state)
        allocate (medium%state, source=state)
    end subroutine

end module
