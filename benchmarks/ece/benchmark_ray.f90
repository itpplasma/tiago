program benchmark_ece_ray
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use ece_ray_fixture, only: polynomial_medium_t
    use tiago_ece_ray, only: cold_ray_t, trace_cold_rays, ECE_O_MODE
    use fortnum_status, only: fortnum_status_t
    implicit none
    integer, parameter :: n = 512, repetitions = 20
    real(dp), parameter :: pi = acos(-1.0_dp), frequency = 9.0e10_dp
    real(dp), &
        parameter :: critical = (2*pi*frequency)**2*8.8541878128e-12_dp &
                     *9.1093837015e-31_dp/1.602176634e-19_dp**2
    type(polynomial_medium_t) :: medium
    type(cold_ray_t), allocatable :: rays(:)
    type(fortnum_status_t), allocatable :: statuses(:)
    type(fortnum_status_t) :: status
    integer(int64) :: first, last, rate
    real(dp) :: initial(6, n), frequencies(n), durations(n), error, &
                expected(6)
    integer :: i, j, modes(n), points
    medium%p(3) = 1.0_dp
    medium%p(8) = critical
    do i = 1, n
        initial(:, i) = [-0.7_dp, 0.1_dp, 0.0_dp, 0.7_dp, 0.1_dp, 0.0_dp]
    end do
    frequencies = frequency
    durations = 2*pi
    modes = ECE_O_MODE
    call trace_cold_rays(medium, frequencies, modes, initial, durations, &
                         rays, statuses, status)
    if (status%code /= 0) error stop 'benchmark warm-up failed'
    call system_clock(first, rate)
    do j = 1, repetitions
        call trace_cold_rays(medium, frequencies, modes, initial, durations, &
                             rays, statuses, status)
        if (status%code /= 0) error stop 'benchmark batch failed'
    end do
    call system_clock(last)
    error = 0.0_dp
    points = 0
    do i = 1, n
        expected = initial(:, i)
        error = max(error, maxval(abs(rays(i)%trace%y(:, &
                                                    size(rays(i)%trace%t)) - expected)))
        points = points + size(rays(i)%trace%t)
    end do
    print '(A)', &
        'rays,repetitions,seconds,rays_per_second,max_orbit_error,accepted_points'
    print '(I0,A,I0,A,ES16.8,A,ES16.8,A,ES16.8,A,I0)', n, ',', repetitions, &
        ',', &
        real(last - first, dp)/real(rate, dp), ',', &
        real(n*repetitions, dp)*real(rate, dp)/real(last - first, dp), ',', error, &
        ',', points
end program
