program benchmark_ece_transfer
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use tiago_ece_transfer, only: ece_transfer_t, thermal_transfer_batch
    implicit none
    integer, parameter :: n = 256, nr = 4096, repeats = 20
    real(dp), allocatable :: te(:, :), alpha(:, :), length(:, :), rho(:), wall(:)
    type(ece_transfer_t), allocatable :: trace(:)
    integer, allocatable :: status(:)
    integer(int64) :: start, finish, rate
    real(dp) :: elapsed, error, reference
    integer :: k
    allocate (te(n, nr), alpha(n, nr), length(n, nr), rho(nr), wall(nr), trace(nr), &
              status(nr))
    te = 1000.0_dp
    alpha = .7_dp
    length = 1.0_dp/n
    rho = 0.0_dp
    wall = 0.0_dp
    call thermal_transfer_batch(te, alpha, length, rho, wall, trace, status)
    call system_clock(start, rate)
    do k = 1, repeats
        call thermal_transfer_batch(te, alpha, length, rho, wall, trace, status)
        if (any(status /= 0)) error stop 'benchmark transfer failed'
    end do
    call system_clock(finish)
    elapsed = real(finish - start, dp)/rate
    reference = 1000*(1 - exp(-.7_dp))
    error = maxval(abs(trace%temperature - reference))
    if (error > 1.0e-9_dp) error stop 'benchmark accuracy'
    print '(A)', 'rays,cells,repeats,seconds,rays_per_second,max_error_eV'
    print '(I0,A,I0,A,I0,A,ES16.8,A,ES16.8,A,ES16.8)', &
        nr, ',', n, ',', repeats, ',', elapsed, ',', nr*real(repeats, dp)/elapsed, &
            ',', error
end program
