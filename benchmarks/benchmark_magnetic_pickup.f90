program benchmark_magnetic_pickup
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use omp_lib, only: omp_get_wtime
    use tiago_magnetic_pickup, only: magnetic_pickup_t
    use tiago_calibration, only: calibration_prior_t
    implicit none
    type(magnetic_pickup_t) :: model
    type(calibration_prior_t) :: prior
    integer(int64), allocatable :: clock(:)
    complex(dp), allocatable :: h(:)
    real(dp), allocatable :: flux(:), y(:)
    real(dp) :: cov(2, 2), start, finish, checksum
    integer :: n, j, repeat, status, m
    integer, parameter :: sizes(3) = [729, 2187, 59049]
    cov = reshape([1.e-4_dp, 0.0_dp, 0.0_dp, 1.e-12_dp], [2, 2])
    call prior%init([1.0_dp, 0.0_dp], cov, status)
    prior%record%source_id = 'synthetic'
    prior%record%record_id = 'identity electronics'
    prior%record%metadata = 'H: ideal voltage ratio, all stages identity'
    prior%valid_from = 0
    prior%valid_to = 1
    do m = 1, size(sizes)
        n = sizes(m)
        allocate (clock(n), h(n), flux(n), y(n))
        do j = 1, n
            clock(j) = int(j - 1, int64)*1000_int64
            flux(j) = 1.e-6_dp*cos(2*acos(-1.0_dp)*7*(j - 1)/n)
        end do
        h = 1
        call model%init(clock, 0_int64, prior, h, 1.e6_dp, 'Wb', 'V', 'V/V', status)
        if (status /= 0) error stop 'benchmark init'
        repeat = max(10, 2000000/n)
        call model%predict(flux, [0.0_dp, 0.0_dp], y, status)
        start = omp_get_wtime()
        do j = 1, repeat
            call model%predict(flux, [0.0_dp, 0.0_dp], y, status)
            if (status /= 0) error stop 'benchmark predict'
        end do
        finish = omp_get_wtime()
        checksum = sum(abs(y))
        write (*, '(i8,1x,i8,1x,es16.8,1x,es16.8)') &
            n, repeat, (finish - start)/repeat, checksum
        deallocate (clock, h, flux, y)
    end do
end program
