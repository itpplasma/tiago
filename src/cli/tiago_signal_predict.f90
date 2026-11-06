program tiago_signal_predict
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_signal_correlation, only: signal_correlation_t
    use fortnum_status, only: fortnum_status_t
    implicit none
    type :: query_t
        integer, allocatable :: observed(:), missing(:)
        real(dp), allocatable :: values(:), location(:), scale(:, :)
        real(dp) :: df, log_observed
        type(fortnum_status_t) :: status
    end type
    type(signal_correlation_t) :: model
    type(fortnum_status_t) :: status
    type(query_t), allocatable :: queries(:)
    real(dp), allocatable :: samples(:, :), mean0(:), psi0(:, :)
    real(dp) :: kappa0, nu0
    character(len=4096) :: input, output
    integer :: u, v, d, n, cases, k, i, no, nm
    call get_command_argument(1, input)
    call get_command_argument(2, output)
    open(newunit=u, file=trim(input), status='old', action='read')
    read(u, *) d, n, cases, kappa0, nu0
    if (d < 1 .or. n < 0 .or. cases < 1) error stop 'invalid prediction dimensions'
    allocate(samples(d, n), mean0(d), psi0(d, d), queries(cases))
    read(u, *) mean0
    do i = 1, d
        read(u, *) psi0(i, :)
    end do
    do i = 1, n
        read(u, *) samples(:, i)
    end do
    call model%fit(samples, mean0, psi0, kappa0, nu0, status)
    if (status%code /= 0) error stop 'invalid signal training model'
    do k = 1, cases
        read(u, *) no
        if (no < 0 .or. no > d) error stop 'invalid observed signal count'
        allocate(queries(k)%observed(no), queries(k)%values(no))
        if (no > 0) then
            read(u, *) queries(k)%observed
            read(u, *) queries(k)%values
        end if
    end do
    close(u)
    !$omp parallel do default(shared) private(k) schedule(static)
    do k = 1, cases
        call model%predict(queries(k)%observed, queries(k)%values, &
            queries(k)%missing, queries(k)%location, queries(k)%scale, &
            queries(k)%df, queries(k)%log_observed, queries(k)%status)
    end do
    !$omp end parallel do
    open(newunit=v, file=trim(output), status='replace', action='write')
    do k = 1, cases
        if (queries(k)%status%code /= 0) error stop 'invalid signal prediction'
        nm = size(queries(k)%missing)
        write(v, *) k, nm, queries(k)%df, queries(k)%log_observed
        if (nm > 0) then
            write(v, *) queries(k)%missing
            write(v, *) queries(k)%location
            do i = 1, nm
                write(v, *) queries(k)%scale(i, :)
            end do
        end if
    end do
    close(v)
end program
