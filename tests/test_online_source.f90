module binary_test_decoder
    use, intrinsic :: iso_fortran_env, only: int8
    use tiago_rest_adapter, only: payload_decoder_t
    use tiago_signal_source, only: signal_request_t, signal_window_t, &
        SIGNAL_ERROR, SIGNAL_OK, SIGNAL_MORE
    implicit none
    type, extends(payload_decoder_t) :: binary_decoder_t
        integer(int8) :: buffer(6)
        integer :: used = 0
    contains
        procedure :: decode
    end type
contains
    subroutine decode(self, request, bytes, final, window, status)
        class(binary_decoder_t), intent(inout) :: self
        type(signal_request_t), intent(in) :: request
        integer(int8), intent(in) :: bytes(:)
        logical, intent(in) :: final
        type(signal_window_t), intent(out) :: window
        integer, intent(out) :: status
        status = SIGNAL_ERROR
        if (size(bytes) > 6 - self%used) then
            self%used = 0
            return
        end if
        self%buffer(self%used + 1:self%used + size(bytes)) = bytes
        self%used = self%used + size(bytes)
        if (.not. final) then
            status = SIGNAL_MORE
            return
        end if
        if (self%used /= 6) then
            self%used = 0
            return
        end if
        self%used = 0
        if (any(self%buffer /= &
            [0_int8, 1_int8, 127_int8, int(-128, int8), -1_int8, 0_int8])) return
        window%time = [request%first]
        window%samples = [42.0d0]
        window%variance = [0.0d0]
        window%units = 'test'
        window%provenance%source_id = 'synthetic-binary'
        window%provenance%version = 'v1'
        window%provenance%request_key = request%key()
        window%provenance%record_id = 'binary'
        window%provenance%metadata = 'synthetic bytes including NUL and high bits'
        call window%seal()
        status = SIGNAL_OK
    end subroutine
end module

program test_online_source
    use, intrinsic :: iso_fortran_env, only: dp => real64, int8
    use tiago_signal_source, only: signal_request_t, signal_window_t, &
        cached_source_t, backend_request_t, SIGNAL_OK
    use tiago_rest_adapter, only: rest_adapter_t
    use tiago_online_source, only: online_source_t
    use tiago_signal_reduction, only: linear_reduction_t
    use tiago_http_client, only: http_get, http_policy_t, http_clock
    use binary_test_decoder, only: binary_decoder_t
    implicit none
    type(rest_adapter_t), target :: adapter
    type(binary_decoder_t), target :: binary
    type(online_source_t), target :: source
    type(cached_source_t) :: cache
    type(linear_reduction_t), target :: step
    type(signal_request_t) :: request
    type(signal_window_t) :: window, second
    type(backend_request_t) :: backend
    type(http_policy_t) :: policy
    integer(int8), allocatable :: bytes(:)
    character(2048) :: directory, base
    character(16), parameter :: failures(*) = [character(16) :: &
        'wrongversion', 'malformed', 'outside', 'null', 'badjson', 'negative', 'truncated']
    real(dp) :: start, elapsed
    integer :: status, i, unit
    call get_command_argument(1, directory)
    call adapter%configure(trim(directory)//'/adapter.json', status)
    call require(status == SIGNAL_OK, 'configure')
    source%adapter => adapter
    source%policy = adapter%policy
    request%channel = 'ok'
    request%version = 'v1'
    request%first = 0.0_dp
    request%last = 1.0_dp
    call source%read(request, window, status)
    call require(status == SIGNAL_OK, 'JSON source')
    call require(size(window%samples) == 3, 'JSON shape')
    call require(maxval(abs(window%samples - [3.0_dp, 5.0_dp, 7.0_dp])) &
        < 1e-14_dp, 'synthetic ramp and unit conversion')
    call require(maxval(abs(window%variance - 1.0_dp)) < 1e-14_dp, &
        'variance conversion')
    print '(A,ES12.4)', 'ramp_max_absolute_error_eV=', &
        maxval(abs(window%samples - [3.0_dp, 5.0_dp, 7.0_dp]))
    call cache%init(source, 3)
    request%channel = 'split'
    request%last = 2.5_dp
    call cache%read(request, window, status)
    call require(status == SIGNAL_OK, 'chunked source')
    call require(size(window%samples) == 9, 'three windows')
    call require(all(window%time(2:) > window%time(:8)), 'no duplicate boundary')
    call require(maxval(abs(window%samples - (4*window%time + 3))) &
        < 1e-13_dp, 'chunked values')
    call cache%read(request, second, status)
    call require(status == SIGNAL_OK, 'cache hit')
    call require(window%provenance%checksum == second%provenance%checksum, &
        'cache content')
    call source%read(request, second, status)
    call require(status == SIGNAL_OK, 'repeat uncached chunked read')
    request%last = 1.0_dp
    do i = 1, size(failures)
        request%channel = trim(failures(i))
        call source%read(request, window, status)
        call require(status /= SIGNAL_OK, 'reject '//request%channel)
    end do
    request%channel = 'chunked'
    call source%read(request, window, status)
    call require(status == SIGNAL_OK, 'HTTP chunk framing and decoder recovery')
    call require(size(window%samples) == 3, 'chunk-framed sample count')
    source%policy%max_response_bytes = 128
    request%channel = 'oversize'
    call source%read(request, window, status)
    call require(status /= SIGNAL_OK, 'streamed response bound')
    source%policy = adapter%policy
    request%channel = 'empty'
    call source%read(request, window, status)
    call require(status == SIGNAL_OK, 'empty array')
    call require(size(window%samples) == 0, 'empty result')
    adapter%decoder => binary
    request%channel = 'binary'
    call source%read(request, window, status)
    call require(status == SIGNAL_OK, 'binary decoder')
    call require(abs(window%samples(1) - 42.0_dp) < 1e-14_dp, 'binary value')
    nullify(adapter%decoder)
    request%channel = 'ok'
    request%version = 'v1/a b'
    call adapter%map_request(request, backend, status)
    call require(status == SIGNAL_OK, 'URL mapping')
    call require(index(backend%endpoint, 'v1%2Fa%20b') > 0, 'URL escape')
    request%version = 'v1'
    ! Process one bounded window at a time before retaining reduced observations.
    step%id = 'three-sample-mean-v1'
    step%units = 'eV'
    allocate (step%weights(1, 3), step%noise_covariance(1, 1))
    step%weights = 1.0_dp/3
    step%noise_covariance = 0.0_dp
    source%reduction => step
    request%last = 2.5_dp
    call source%read(request, window, status)
    call require(status == SIGNAL_OK, 'online reduction')
    call require(size(window%samples) == 3, 'reduced windows')
    call require(maxval(abs(window%samples - [5.0_dp, 9.0_dp, 12.0_dp])) &
        < 1e-13_dp, 'reduced analytic means')
    nullify(source%reduction)
    source%policy%max_output_samples = 1
    call source%read(request, window, status)
    call require(status /= SIGNAL_OK, 'output bound')
    source%policy = adapter%policy
    source%policy%max_windows = 1
    call source%read(request, window, status)
    call require(status /= SIGNAL_OK, 'window count bound')
    source%policy = adapter%policy
    source%policy%max_retries = 1
    request%channel = 'windowretry'
    call source%read(request, window, status)
    call require(status /= SIGNAL_OK, 'retry budget shared across windows')
    request%channel = 'always500'
    call source%read(request, window, status)
    call require(status /= SIGNAL_OK, 'retry cap across windows')
    source%policy = adapter%policy
    request%channel = 'versioned'
    deallocate (request%version)
    call source%read(request, window, status)
    call require(status == SIGNAL_OK, 'resolve then pin revision')
    call require(window%provenance%version == 'v1', 'resolved revision')
    request%version = 'v1'
    ! Decode enough samples to exercise OpenMP unit conversion, and measure it.
    source%policy = adapter%policy
    request%channel = 'large'
    request%last = 1.0_dp
    start = http_clock()
    call source%read(request, window, status)
    elapsed = http_clock() - start
    call require(status == SIGNAL_OK, 'large JSON')
    call require(size(window%samples) == 20000, 'large size')
    call require(maxval(abs(window%samples - (4*window%time+3))) &
        < 1e-13_dp, 'large values')
    print '(A,F10.6)', 'large_20000_samples_elapsed_s=', elapsed
    open (newunit=unit, file=trim(directory)//'/base.txt')
    read (unit, '(A)') base
    close (unit)
    policy%min_interval = 0.0_dp
    call http_get(trim(base)//'/stats', policy, bytes, status)
    open (newunit=unit, file=trim(directory)//'/stats-after.json', &
        access='stream', form='unformatted', status='replace')
    write (unit) bytes
    close (unit)
    print *, 'Online adapter, cache, decoder and reduction controls passed'
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print *, 'FAIL: ', message
            stop 1
        end if
    end subroutine
end program
