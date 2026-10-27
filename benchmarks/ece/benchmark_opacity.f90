program benchmark_opacity
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use tiago_ece_opacity
    implicit none
    integer, parameter :: n = 131072, reps = 16
    real(dp), allocatable :: te(:), ne(:), b(:), f(:), a(:), j(:), n2(:), rows(:, :, :)
    integer, allocatable :: mode(:), harmonic(:), status(:)
    integer :: k, r
    integer(int64) :: first, last, rate
    real(dp) :: elapsed
    allocate (te(n), ne(n), b(n), f(n), a(n), j(n), n2(n), rows(3, 4, n), &
              mode(n), harmonic(n), status(n))
    te = 510.99895069_dp
    ne = .05_dp*8.8541878188e-12_dp*4/9.1093837139e-31_dp
    b = 2.0_dp
    harmonic = 2
    mode = ECE_X_MODE
    do k = 1, n
        f(k) = 2*1.602176634e-19_dp*2/(2*acos(-1.0_dp)*9.1093837139e-31_dp* &
                                       (1 + .001_dp*(1.0_dp + 4.0_dp*k/n)))
    end do
    call ece_opacity_batch(te, ne, b, f, mode, harmonic, a, j, n2, rows, status)
    if (any(status /= 0)) error stop 'benchmark opacity status'
    call system_clock(first, rate)
    do r = 1, reps
        call ece_opacity_batch(te, ne, b, f, mode, harmonic, a, j, n2, rows, status)
    end do
    call system_clock(last)
    elapsed = real(last - first, dp)/real(rate, dp)
    if (any(status /= 0)) error stop 'benchmark opacity status'
    print '(A)', 'samples,repetitions,seconds,samples_per_second,alpha_sum_per_m'
    print '(I0,A,I0,A,ES16.8,A,ES16.8,A,ES24.16)', n, ',', reps, ',', elapsed, ',', &
        real(n, dp)*reps/elapsed, ',', sum(a)
end program
