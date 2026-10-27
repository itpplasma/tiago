program benchmark_archive
    use, intrinsic :: iso_fortran_env, only: dp => real64, int8
    use tiago_rest_adapter, only: rest_adapter_t
    use tiago_signal_source, only: signal_request_t, signal_window_t, &
        SIGNAL_OK, SIGNAL_MORE
    use tiago_http_client, only: http_clock
    implicit none
    type(rest_adapter_t) :: adapter
    type(signal_request_t) :: request
    type(signal_window_t) :: window
    integer(int8), allocatable :: bytes(:)
    character(2048) :: config, payload
    integer :: unit, n, status, i
    real(dp) :: start, elapsed, error
    call get_command_argument(1, config)
    call get_command_argument(2, payload)
    call adapter%configure(trim(config), status)
    if (status /= SIGNAL_OK) stop 1
    open (newunit=unit, file=trim(payload), access='stream', form='unformatted')
    inquire (unit=unit, size=n)
    allocate (bytes(n))
    read (unit) bytes
    close (unit)
    request%channel = 'large'
    request%version = 'v1'
    request%first = 0.0_dp
    request%last = 1.0_dp
    print '(A)', 'repeat,samples,elapsed_s,max_absolute_error_eV'
    do i = 1, 6
        start = http_clock()
        call adapter%decode(request, bytes, .false., window, status)
        if (status /= SIGNAL_MORE) stop 2
        call adapter%decode(request, [integer(int8) ::], .true., window, status)
        if (status /= SIGNAL_OK) stop 3
        elapsed = http_clock() - start
        error = maxval(abs(window%samples - (4.0_dp*window%time + 3.0_dp)))
        if (error > 1e-13_dp) stop 4
        ! The first pass warms allocations and the OpenMP runtime.
        if (i > 1) print '(I0,A,I0,A,F12.8,A,ES12.4)', &
            i - 1, ',', size(window%samples), ',', elapsed, ',', error
    end do
end program
