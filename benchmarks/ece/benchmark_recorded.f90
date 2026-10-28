program benchmark_ece_recorded
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use ece_recorded_fixture, only: configure
    use tiago_ece_state_medium, only: state_medium_t
    use tiago_ece_radiometer, only: ece_t
    use tiago_ece_recorded, only: recorded_ece_t, record_ece
    use fortnum_status, only: fortnum_status_t
    implicit none
    integer, parameter :: n = 64, repetitions = 4
    type(state_medium_t) :: medium
    type(ece_t) :: diagnostic
    type(recorded_ece_t) :: tape
    type(fortnum_status_t) :: status
    real(dp), allocatable :: pd(:), pb(:), eta(:), ed(:), eb(:)
    real(dp) :: yd(n), vd(n), yb(n), vb(n), elapsed(3), checksum(3)
    integer(int64) :: first, last, rate
    integer :: i
    call configure(medium, diagnostic, n, 48, 1.0_dp, .true.)
    allocate (pd(medium%n_params()), pb(medium%n_params()), &
        eta(diagnostic%n_nuisance()), ed(diagnostic%n_nuisance()), &
        eb(diagnostic%n_nuisance()))
    pd = 0
    pd(size(pd) - 5) = 1
    eta = 0
    ed = .1_dp
    yb = 1
    vb = 1000
    call record_ece(diagnostic, medium, eta, tape, status)
    if (status%code /= 0) error stop 'recorded warm-up'
    call system_clock(first, rate)
    do i = 1, repetitions
        call record_ece(diagnostic, medium, eta, tape, status)
        if (status%code /= 0) error stop 'recorded benchmark'
    end do
    call system_clock(last)
    elapsed(1) = real(last - first, dp)/rate
    checksum(1) = sum(tape%voltage)
    call system_clock(first)
    do i = 1, repetitions
        call tape%jvp(pd, ed, yd, vd, status)
        if (status%code /= 0) error stop 'recorded tangent benchmark'
    end do
    call system_clock(last)
    elapsed(2) = real(last - first, dp)/rate
    checksum(2) = sum(yd) + sum(vd)
    call system_clock(first)
    do i = 1, repetitions
        call tape%vjp(yb, vb, pb, eb, status)
        if (status%code /= 0) error stop 'recorded adjoint benchmark'
    end do
    call system_clock(last)
    elapsed(3) = real(last - first, dp)/rate
    checksum(3) = sum(pb) + sum(eb)
    print '(A)', 'operation,rays,cells,repetitions,seconds,rays_per_second,checksum'
    do i = 1, 3
        print '(I0,3(",",I0),3(",",ES18.10))', i, n, 48, repetitions, elapsed(i), &
            n*repetitions/elapsed(i), checksum(i)
    end do
end program
