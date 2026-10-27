module tiago_rest_adapter
    use, intrinsic :: iso_fortran_env, only: dp => real64, int8
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use json_module, only: json_file
    use tiago_signal_source, only: signal_adapter_t, signal_request_t, &
        signal_window_t, backend_request_t, SIGNAL_OK, SIGNAL_ERROR, SIGNAL_MORE, &
        version_matches
    use tiago_http_client, only: http_policy_t
    use tiago_recorded_affine, only: recorded_affine
    implicit none
    private

    type, abstract, public :: payload_decoder_t
    contains
        procedure(payload_i), deferred :: decode
    end type
    abstract interface
        subroutine payload_i(self, request, bytes, final, window, status)
            import payload_decoder_t, signal_request_t, signal_window_t, int8
            class(payload_decoder_t), intent(inout) :: self
            type(signal_request_t), intent(in) :: request
            integer(int8), intent(in) :: bytes(:)
            logical, intent(in) :: final
            type(signal_window_t), intent(out) :: window
            integer, intent(out) :: status
        end subroutine
    end interface
    type :: channel_t
        character(:), allocatable :: logical_name, url, units
        character(:), allocatable :: time_path, value_path, variance_path
        character(:), allocatable :: version_path, version
        real(dp) :: scale = 1.0_dp, offset = 0.0_dp, time_scale = 1.0_dp
    end type
    type, extends(signal_adapter_t), public :: rest_adapter_t
        private
        type(http_policy_t), public :: policy
        class(payload_decoder_t), pointer, public :: decoder => null()
        type(channel_t), allocatable :: channels(:)
        character(:), allocatable :: source_id, buffer
        integer :: used = 0
    contains
        procedure :: configure => rest_configure
        procedure :: map_request => rest_map
        procedure :: decode => rest_decode
    end type
contains
    subroutine rest_configure(self, path, status)
        class(rest_adapter_t), intent(inout) :: self
        character(*), intent(in) :: path
        integer, intent(out) :: status
        type(json_file) :: config
        integer :: n, i, j
        logical :: found, ok
        character(40) :: prefix
        character(:), allocatable :: base
        status = SIGNAL_ERROR
        self%used = 0
        if (allocated(self%channels)) deallocate (self%channels)
        self%policy = http_policy_t()
        call config%initialize(allow_duplicate_keys=.false., &
            null_to_real_mode=1)
        call config%load(path)
        if (config%failed()) return
        ok = .true.
        call string('source_id', self%source_id)
        call config%info('channels', found=found, n_children=n)
        if (.not. found) return
        if (n < 1 .or. n > 10000) return
        allocate (self%channels(n))
        do i = 1, n
            write (prefix, '("channels(",I0,").")') i
            base = trim(prefix)
            call string(base//'logical', self%channels(i)%logical_name)
            call string(base//'url', self%channels(i)%url)
            call string(base//'units', self%channels(i)%units)
            call string(base//'time_path', self%channels(i)%time_path)
            call string(base//'value_path', self%channels(i)%value_path)
            call string(base//'version_path', self%channels(i)%version_path)
            call string(base//'version', self%channels(i)%version)
            self%channels(i)%variance_path = ''
            call config%get(base//'variance_path', &
                self%channels(i)%variance_path, found)
            if (.not. found) self%channels(i)%variance_path = ''
            call real_option(base//'scale', self%channels(i)%scale)
            call real_option(base//'offset', self%channels(i)%offset)
            call real_option(base//'time_scale', self%channels(i)%time_scale)
            if (.not. ok) exit
            if (self%channels(i)%time_scale <= 0.0_dp) ok = .false.
            if (index(self%channels(i)%url, '{first}') == 0) ok = .false.
            if (index(self%channels(i)%url, '{last}') == 0) ok = .false.
            do j = 1, i - 1
                if (self%channels(i)%logical_name == &
                    self%channels(j)%logical_name) ok = .false.
            end do
        end do
        call real_option('policy.min_interval', self%policy%min_interval)
        call real_option('policy.timeout', self%policy%timeout)
        call real_option('policy.connect_timeout', self%policy%connect_timeout)
        call real_option('policy.backoff', self%policy%backoff)
        call real_option('policy.max_backoff', self%policy%max_backoff)
        call real_option('policy.jitter', self%policy%jitter)
        call real_option('policy.max_elapsed', self%policy%max_elapsed)
        call real_option('policy.max_window', self%policy%max_window)
        call int_option('policy.max_retries', self%policy%max_retries)
        call int_option('policy.max_response_bytes', self%policy%max_response_bytes)
        call int_option('policy.max_windows', self%policy%max_windows)
        call int_option('policy.max_output_samples', self%policy%max_output_samples)
        call config%info('policy.user_agent', found=found)
        if (found) call string('policy.user_agent', self%policy%user_agent)
        if (.not. self%policy%valid()) ok = .false.
        if (config%failed()) ok = .false.
        call config%destroy()
        if (ok) then
            status = SIGNAL_OK
        else
            deallocate (self%channels)
        end if
    contains
        subroutine string(key, value)
            character(*), intent(in) :: key
            character(:), allocatable, intent(out) :: value
            logical :: exists
            call config%get(key, value, exists)
            if (.not. exists) then
                ok = .false.
                return
            end if
            if (len(value) == 0) ok = .false.
        end subroutine
        subroutine real_option(key, value)
            character(*), intent(in) :: key
            real(dp), intent(inout) :: value
            logical :: exists
            call config%info(key, found=exists)
            if (.not. exists) return
            call config%get(key, value)
            if (.not. ieee_is_finite(value)) ok = .false.
        end subroutine
        subroutine int_option(key, value)
            character(*), intent(in) :: key
            integer, intent(inout) :: value
            logical :: exists
            call config%info(key, found=exists)
            if (exists) call config%get(key, value)
        end subroutine
    end subroutine

    integer function channel_index(self, request) result(idx)
        class(rest_adapter_t), intent(in) :: self
        type(signal_request_t), intent(in) :: request
        integer :: i
        idx = 0
        if (.not. request%valid()) return
        if (.not. allocated(self%channels)) return
        do i = 1, size(self%channels)
            if (self%channels(i)%logical_name /= request%channel) cycle
            idx = i
            return
        end do
    end function

    subroutine rest_map(self, request, backend, status)
        class(rest_adapter_t), intent(in) :: self
        type(signal_request_t), intent(in) :: request
        type(backend_request_t), intent(out) :: backend
        integer, intent(out) :: status
        character(:), allocatable :: url, version
        character(40) :: number
        integer :: i
        status = SIGNAL_ERROR
        i = channel_index(self, request)
        if (i == 0) return
        version = self%channels(i)%version
        if (allocated(request%version)) then
            if (len(request%version) > 0) version = request%version
        end if
        url = self%channels(i)%url
        write (number, '(ES25.17E3)') request%first
        call replace(url, '{first}', trim(adjustl(number)))
        write (number, '(ES25.17E3)') request%last
        call replace(url, '{last}', trim(adjustl(number)))
        call replace(url, '{version}', escape(version))
        call replace(url, '{channel}', escape(request%channel))
        if (index(url, '{') > 0 .or. index(url, '}') > 0) return
        ! Escape plus signs in scientific exponents in query parameters.
        call replace(url, '+', '%2B')
        backend%endpoint = url
        backend%method = 'GET'
        backend%body = ''
        backend%request_key = request%key()
        status = SIGNAL_OK
    end subroutine

    subroutine replace(text, pattern, value)
        character(:), allocatable, intent(inout) :: text
        character(*), intent(in) :: pattern, value
        character(:), allocatable :: tail, result
        integer :: k
        tail = text
        result = ''
        do
            k = index(tail, pattern)
            if (k == 0) exit
            result = result//tail(:k - 1)//value
            tail = tail(k + len(pattern):)
        end do
        text = result//tail
    end subroutine

    function escape(text) result(encoded)
        character(*), intent(in) :: text
        character(:), allocatable :: encoded
        character(3) :: hex
        integer :: i, c
        encoded = ''
        do i = 1, len(text)
            c = iachar(text(i:i))
            if (index('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'// &
                '0123456789-._~', text(i:i)) > 0) then
                encoded = encoded//text(i:i)
            else
                write (hex, '("%",Z2.2)') c
                encoded = encoded//hex
            end if
        end do
    end function

    subroutine rest_decode(self, request, bytes, final, window, status)
        class(rest_adapter_t), intent(inout) :: self
        type(signal_request_t), intent(in) :: request
        integer(int8), intent(in) :: bytes(:)
        logical, intent(in) :: final
        type(signal_window_t), intent(out) :: window
        integer, intent(out) :: status
        type(json_file) :: response
        type(channel_t) :: ch
        real(dp), allocatable :: times(:), values(:), variance(:)
        real(dp) :: derivative, d_gain, d_offset
        character(:), allocatable :: version
        logical :: found
        integer :: idx, i, n
        status = SIGNAL_ERROR
        if (associated(self%decoder)) then
            call self%decoder%decode(request, bytes, final, window, status)
            return
        end if
        idx = channel_index(self, request)
        if (idx == 0) return
        if (.not. self%policy%valid()) return
        if (.not. allocated(self%buffer)) then
            allocate (character(self%policy%max_response_bytes) :: self%buffer)
        end if
        if (len(self%buffer) /= self%policy%max_response_bytes) then
            deallocate (self%buffer)
            allocate (character(self%policy%max_response_bytes) :: self%buffer)
            self%used = 0
        end if
        if (size(bytes) > len(self%buffer) - self%used) then
            self%used = 0
            return
        end if
        do i = 1, size(bytes)
            self%buffer(self%used + i:self%used + i) = achar(iand(int(bytes(i)), 255))
        end do
        self%used = self%used + size(bytes)
        if (.not. final) then
            status = SIGNAL_MORE
            return
        end if
        call response%initialize(allow_duplicate_keys=.false., null_to_real_mode=1)
        call response%deserialize(self%buffer(:self%used))
        self%used = 0
        if (response%failed()) return
        ch = self%channels(idx)
        call response%get(ch%time_path, times, found)
        if (.not. found) return
        call response%get(ch%value_path, values, found)
        if (.not. found) return
        call response%get(ch%version_path, version, found)
        if (.not. found) return
        if (len(version) == 0) return
        if (.not. version_matches(request, version)) return
        if (size(times) /= size(values)) return
        n = size(times)
        allocate (variance(n))
        variance = 0.0_dp
        if (len(ch%variance_path) > 0) then
            deallocate (variance)
            call response%get(ch%variance_path, variance, found)
            if (.not. found) return
            if (size(variance) /= n) return
        end if
        call response%destroy()
        if (any(variance < 0.0_dp)) return
        allocate (window%time(n), window%samples(n), window%variance(n))
        !$omp parallel do private(derivative,d_gain,d_offset) if(n>1024)
        do i = 1, n
            call recorded_affine(values(i), ch%scale, ch%offset, &
                window%samples(i), derivative, d_gain, d_offset)
            window%variance(i) = derivative*derivative*variance(i)
            window%time(i) = times(i)*ch%time_scale
        end do
        !$omp end parallel do
        window%units = ch%units
        if (.not. window%valid()) return
        if (any(window%time < request%first)) return
        if (any(window%time > request%last)) return
        window%provenance%source_id = self%source_id
        window%provenance%version = version
        window%provenance%request_key = request%key()
        window%provenance%record_id = request%channel
        window%provenance%metadata = 'REST JSON; time origin: archive-defined; '// &
            'scaling is a fixed unit conversion; zero variance if not supplied'
        call window%seal()
        status = SIGNAL_OK
    end subroutine
end module tiago_rest_adapter
