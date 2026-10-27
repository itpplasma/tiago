program benchmark_ece_geometry
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use ece_geometry_fixture, only: setup
    use tiago_ece_state_medium, only: state_medium_t
    use tiago_ece_ray, only: cold_ray_t, trace_cold_rays, cold_hamiltonian, ECE_O_MODE
    use fortnum_status, only: fortnum_status_t
    implicit none
    integer, parameter :: n = 256, repetitions = 4
    type(state_medium_t) :: medium
    type(cold_ray_t), allocatable :: rays(:)
    type(fortnum_status_t), allocatable :: statuses(:)
    type(fortnum_status_t) :: status
    real(dp) :: initial(6, n), frequencies(n), durations(n), q(4), dq(3, 4), d2q(3, 3, &
                                                                                 4)
    real(dp) :: h, dh(8), d2h(8, 8), expected(6), error
    real(dp), allocatable :: qp(:, :), dqp(:, :, :)
    integer :: modes(n), i
    integer(int64)::first, last, rate
    call setup(medium)
    allocate (qp(4, medium%n_params()), dqp(3, 4, medium%n_params()))
    call medium%sample([3.4_dp, 0.0_dp, 0.2_dp], q, dq, d2q, qp, dqp)
    call cold_hamiltonian([0.8_dp, 0.0_dp, 0.6_dp], q(:3), q(4), 9.0e10_dp, &
                          ECE_O_MODE, h, dh, d2h)
    initial(:3, :) = spread([3.4_dp, 0.0_dp, 0.2_dp], 2, n)
    initial(4:, :) = spread(sqrt(1 - 2*h)*[0.8_dp, 0.0_dp, 0.6_dp], 2, n)
    frequencies = 9.0e10_dp
    durations = 0.15_dp
    modes = ECE_O_MODE
    call trace_cold_rays(medium, frequencies, modes, initial, durations, rays, &
                         statuses, status)
    if (status%code /= 0) error stop 'geometry warm-up failed'
    expected = rays(1)%trace%y(:, size(rays(1)%trace%t))
    call system_clock(first, rate)
    do i = 1, repetitions
        call trace_cold_rays(medium, frequencies, modes, initial, durations, rays, &
                             statuses, status)
        if (status%code /= 0) error stop 'geometry benchmark trace failed'
    end do
    call system_clock(last)
    error = 0
    do i = 1, n
        error = max(error, maxval(abs(rays(i)%trace%y(:, &
                                                      size(rays(i)%trace%t)) - expected)))
    end do
    print '(A)', &
        'rays,repetitions,seconds,rays_per_second,max_endpoint_difference_m_or_1'
    print '(I0,A,I0,A,ES16.8,A,ES16.8,A,ES16.8)', n, ',', repetitions, ',', &
        real(last - first, dp)/real(rate, dp), ',', real(n*repetitions, dp)*real(rate, &
                                                                                 dp)/ &
        real(last - first, dp), ',', error
end program
