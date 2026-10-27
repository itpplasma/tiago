program test_ece_geometry_limits
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use ece_geometry_fixture, only: setup
    use tiago_ece_state_medium, only: state_medium_t
    use tiago_kinetic_profiles, only: kinetic_profiles_t, P_NE
    use tiago_ece_ray, only: cold_ray_t, trace_cold_ray, ECE_O_MODE
    use fortnum_status, only: fortnum_status_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), frequency = 9.0e10_dp
    ! CODATA 2018, as in test_ece_ray; Stix, Waves in Plasmas (1992), ch. 1:
    ! perpendicular ordinary mode n^2=1-ne/ncritical. This gives the two
    ! independent analytic Hamiltonian trajectories used in PR #47's tests.
    real(dp), parameter :: critical = (2*pi*frequency)**2*8.8541878128e-12_dp* &
                           9.1093837015e-31_dp/1.602176634e-19_dp**2
    type(state_medium_t) :: medium
    type(cold_ray_t) :: ray
    type(fortnum_status_t) :: status
    real(dp) :: initial(6), expected(6), t, error, coarse, q(4), dq(3, 4), d2q(3, 3, 4)
    real(dp), allocatable :: p(:, :), dp_(:, :, :)
    integer :: n, j, kind
    do kind = 1, 2
        coarse = huge(1.0_dp)
        do n = 65, 129, 64
            call setup(medium)
            call density_fit(medium, n, kind)
            allocate (p(4, medium%n_params()), dp_(3, 4, medium%n_params()))
            if (kind == 1) then
                initial = [3.35_dp, 0.0_dp, 0.0_dp, sqrt(1 - 1.2_dp*0.35_dp), 0.0_dp, &
                           0.0_dp]
            else
                initial = [3.5_dp, 0.0_dp, 0.1_dp, 0.3_dp, 0.0_dp, sqrt(0.55_dp)]
            end if
            call medium%sample(initial(:3), q, dq, d2q, p, dp_)
            ! Use the exact shell of the interpolated density at launch.
            initial(4:) = initial(4:)*sqrt((1 - q(4)/critical)/sum(initial(4:)**2))
            if (kind == 1) then
                call trace_cold_ray(medium, frequency, ECE_O_MODE, initial, 1.8_dp, &
                                    ray, status, rtol=1.0e-11_dp, atol=1.0e-13_dp, &
                                    hmax=0.01_dp)
            else
                call trace_cold_ray(medium, frequency, ECE_O_MODE, initial, 2*pi, &
                                    ray, status, rtol=1.0e-11_dp, atol=1.0e-13_dp, &
                                    hmax=0.02_dp)
            end if
            if (status%code /= 0) then
                print *, status%msg
                error stop 'VMEC analytic-limit trace failed'
            end if
            error = 0
            do j = 1, size(ray%trace%t)
                t = ray%trace%t(j)
                expected = initial
                if (kind == 1) then
                    expected(1) = initial(1) + initial(4)*t - 1.2_dp*t*t/4
                    expected(4) = initial(4) - 1.2_dp*t/2
                else
                    expected(1) = 3 + (initial(1) - 3)*cos(t) + initial(4)*sin(t)
                    expected(3) = initial(3)*cos(t) + initial(6)*sin(t)
                    expected(4) = initial(4)*cos(t) - (initial(1) - 3)*sin(t)
                    expected(6) = initial(6)*cos(t) - initial(3)*sin(t)
                end if
                error = max(error, maxval(abs(ray%trace%y(:, j) - expected)))
            end do
            print '(A,I0,A,I0,A,ES14.6)', 'limit=', kind, ', knots=', n, &
                ', max error [m/1]=', error
            if (n == 129) then
                if (error > 5.0e-7_dp) error stop 'VMEC slab/cylinder analytic error'
                if (error > coarse/3) error stop 'VMEC profile refinement failed'
            end if
            coarse = error
            deallocate (p, dp_)
        end do
    end do
contains
    subroutine density_fit(m, n, kind)
        type(state_medium_t), intent(inout)::m
        integer, intent(in)::n, kind
        type(kinetic_profiles_t)::prof
        real(dp), allocatable::values(:), z(:)
        integer::i
        call prof%init([(real(i - 1, dp)/(n - 1), i=1, n)])
        allocate (values(prof%ncoef))
        values = critical*(0.1_dp + prof%greville)
        if (kind == 1) values = critical*1.2_dp*sqrt(max(0.01_dp, prof%greville))
        call prof%set_profile(P_NE, values, .true., sigma_f=0.1_dp, l_core=0.4_dp, &
                              l_edge=0.2_dp, s0=0.5_dp, sw=0.2_dp)
        allocate (z(prof%n_latent))
        z = 0
        call prof%fit_latent(P_NE, prof%greville, values, z)
        prof%p(P_NE)%mean = prof%coefficients(z, P_NE)
        prof%p(P_NE)%latent = .false.
        prof%n_latent = 0
        m%profiles = prof
        m%z = [real(dp)::]
    end subroutine
end program
