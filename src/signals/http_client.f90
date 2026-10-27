module tiago_http_client
    use, intrinsic :: iso_c_binding, only: c_ptr, c_funptr, c_int, c_long, &
        c_double, c_size_t, c_int8_t, c_char, c_null_char, c_associated, c_f_pointer, &
        c_loc, c_funloc
    use, intrinsic :: iso_fortran_env, only: dp => real64, int8
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    integer, parameter, public :: HTTP_OK = 0, HTTP_ERROR = 1, HTTP_OVERSIZE = 3
    type, public :: http_policy_t
        real(dp) :: min_interval = 0.5_dp, timeout = 30.0_dp
        real(dp) :: connect_timeout = 10.0_dp, backoff = 0.5_dp
        real(dp) :: max_backoff = 30.0_dp, jitter = 0.25_dp
        real(dp) :: max_elapsed = 120.0_dp, max_window = 60.0_dp
        integer :: max_retries = 3, max_response_bytes = 8388608
        integer :: max_windows = 10000, max_output_samples = 1000000
        character(:), allocatable :: user_agent
    contains
        procedure :: valid => policy_valid
    end type
    abstract interface
        subroutine http_consumer_i(bytes, final, context, status)
            import int8
            integer(int8), intent(in) :: bytes(:)
            logical, intent(in) :: final
            class(*), intent(inout) :: context
            integer, intent(out) :: status
        end subroutine
    end interface
    type :: receive_t
        type(c_ptr) :: bytes
        integer(c_size_t) :: used, capacity
        integer(c_int) :: overflow
        type(c_ptr) :: handle
        procedure(http_consumer_i), pointer, nopass :: consumer => null()
        class(*), pointer :: context => null()
        logical :: consumer_error = .false., delivered = .false.
    end type
    interface
        function curl_global() bind(C, name='tiago_curl_global') result(rc)
            import c_int
            integer(c_int) :: rc
        end function
        function curl_init() bind(C, name='curl_easy_init') result(handle)
            import c_ptr
            type(c_ptr) :: handle
        end function
        function curl_perform(handle) bind(C, name='curl_easy_perform') result(rc)
            import c_ptr, c_int
            type(c_ptr), value :: handle
            integer(c_int) :: rc
        end function
        subroutine curl_cleanup(handle) bind(C, name='curl_easy_cleanup')
            import c_ptr
            type(c_ptr), value :: handle
        end subroutine
        function curl_setup(handle, url, agent, timeout, connect_timeout, &
                writer, data) bind(C, name='tiago_curl_setup') result(rc)
            import c_ptr, c_funptr, c_long, c_int, c_char
            type(c_ptr), value :: handle, data
            character(c_char), intent(in) :: url(*), agent(*)
            integer(c_long), value :: timeout, connect_timeout
            type(c_funptr), value :: writer
            integer(c_int) :: rc
        end function
        function curl_info(handle, code, retry) &
                bind(C, name='tiago_curl_info') result(rc)
            import c_ptr, c_long, c_double, c_int
            type(c_ptr), value :: handle
            integer(c_long), intent(out) :: code
            real(c_double), intent(out) :: retry
            integer(c_int) :: rc
        end function
        function host_acquire(url, interval, deadline) &
                bind(C, name='tiago_http_acquire') result(gate)
            import c_ptr, c_double, c_char
            character(c_char), intent(in) :: url(*)
            real(c_double), value :: interval, deadline
            type(c_ptr) :: gate
        end function
        subroutine host_release(gate, delay) bind(C, name='tiago_http_release')
            import c_ptr, c_double
            type(c_ptr), value :: gate
            real(c_double), value :: delay
        end subroutine
        function http_clock() bind(C, name='tiago_http_clock') result(seconds)
            import c_double
            real(c_double) :: seconds
        end function
        subroutine http_sleep(seconds) bind(C, name='tiago_http_sleep')
            import c_double
            real(c_double), value :: seconds
        end subroutine
    end interface
    public :: http_get, http_clock, http_sleep
contains
    logical function policy_valid(self)
        class(http_policy_t), intent(in) :: self
        real(dp) :: settings(9)
        settings = [self%min_interval, self%timeout, self%connect_timeout, &
            self%backoff, self%max_backoff, self%jitter, self%max_elapsed, &
            self%max_window, real(self%max_response_bytes, dp)]
        policy_valid = .false.
        if (.not. all(ieee_is_finite(settings))) return
        if (any(settings < 0.0_dp)) return
        if (min(self%timeout, self%connect_timeout, self%max_elapsed, &
            self%max_window) <= 0.0_dp) return
        if (self%timeout > real(huge(0_c_long), dp)/1000.0_dp) return
        if (self%connect_timeout > real(huge(0_c_long), dp)/1000.0_dp) return
        if (self%max_retries < 0 .or. self%max_retries > 100) return
        if (self%max_response_bytes < 1 .or. self%max_windows < 1) return
        if (self%max_output_samples < 1 .or. self%jitter > 1.0_dp) return
        if (allocated(self%user_agent)) then
            if (len_trim(self%user_agent) == 0) return
            if (index(self%user_agent, achar(0)) > 0) return
            if (index(self%user_agent, achar(10)) > 0) return
            if (index(self%user_agent, achar(13)) > 0) return
        end if
        policy_valid = .true.
    end function

    function receive(data, width, count, userdata) bind(C) result(written)
        type(c_ptr), value :: data, userdata
        integer(c_size_t), value :: width, count
        integer(c_size_t) :: written, n
        integer(c_long) :: code
        real(c_double) :: retry
        integer(c_int) :: rc
        integer :: status
        type(receive_t), pointer :: state
        integer(c_int8_t), pointer :: incoming(:), buffer(:)
        written = 0_c_size_t
        call c_f_pointer(userdata, state)
        if (width /= 0_c_size_t) then
            if (count > huge(n)/width) then
                state%overflow = 1
                return
            end if
        end if
        n = width*count
        if (n > state%capacity - state%used) then
            state%overflow = 1
            return
        end if
        call c_f_pointer(data, incoming, [int(n)])
        if (associated(state%consumer)) then
            rc = curl_info(state%handle, code, retry)
            if (rc /= 0) return
            if (code >= 200 .and. code < 300) then
                state%delivered = .true.
                call state%consumer(incoming, .false., state%context, status)
                if (status /= 0) then
                    state%consumer_error = .true.
                    return
                end if
            end if
        else
            call c_f_pointer(state%bytes, buffer, [int(state%capacity)])
            buffer(int(state%used) + 1:int(state%used + n)) = incoming
        end if
        state%used = state%used + n
        written = n
    end function

    subroutine http_get(url, policy, bytes, status, response_code, retries_left, &
            deadline, consumer, context)
        character(*), intent(in) :: url
        type(http_policy_t), intent(in) :: policy
        integer(int8), allocatable, intent(out) :: bytes(:)
        integer, intent(out) :: status
        integer, optional, intent(out) :: response_code
        integer, optional, intent(inout) :: retries_left
        real(dp), optional, intent(in) :: deadline
        procedure(http_consumer_i), optional :: consumer
        class(*), optional, target, intent(inout) :: context
        integer(c_int8_t), allocatable, target :: storage(:)
        type(receive_t), target :: state
        type(c_ptr) :: handle, gate
        integer(c_int) :: rc, info_rc
        integer(c_long) :: code
        integer :: budget, attempt, consume_status
        real(dp) :: retry, delay, random, end_time, remaining, exponential
        character(:), allocatable :: agent
        status = HTTP_ERROR
        if (present(response_code)) response_code = 0
        if (.not. policy%valid()) return
        if (index(url, achar(0)) > 0) return
        if (curl_global() /= 0) return
        end_time = http_clock() + policy%max_elapsed
        if (present(deadline)) end_time = min(end_time, deadline)
        budget = policy%max_retries
        if (present(retries_left)) budget = min(budget, retries_left)
        agent = 'TIAGO/1 (+https://github.com/itpplasma/tiago)'
        if (allocated(policy%user_agent)) agent = policy%user_agent
        if (present(consumer) .neqv. present(context)) return
        if (present(consumer)) then
            allocate (storage(1))
        else
            allocate (storage(policy%max_response_bytes))
        end if
        attempt = 0
        exponential = min(policy%backoff, policy%max_backoff)
        do
            gate = host_acquire(url//c_null_char, policy%min_interval, end_time)
            if (.not. c_associated(gate)) return
            handle = curl_init()
            if (.not. c_associated(handle)) then
                call host_release(gate, 0.0_dp)
                return
            end if
            remaining = min(policy%timeout, end_time - http_clock())
            state%bytes = c_loc(storage)
            state%used = 0_c_size_t
            state%capacity = int(policy%max_response_bytes, c_size_t)
            state%overflow = 0
            state%consumer_error = .false.
            state%delivered = .false.
            state%handle = handle
            if (present(consumer)) then
                state%consumer => consumer
                state%context => context
            end if
            rc = curl_setup(handle, url//c_null_char, agent//c_null_char, &
                max(1_c_long, int(remaining*1000.0_dp, c_long)), &
                max(1_c_long, int(policy%connect_timeout*1000.0_dp, c_long)), &
                c_funloc(receive), c_loc(state))
            if (rc == 0) rc = curl_perform(handle)
            info_rc = curl_info(handle, code, retry)
            call curl_cleanup(handle)
            if (present(response_code)) response_code = int(code)
            delay = 0.0_dp
            if (state%overflow /= 0) then
                call host_release(gate, 0.0_dp)
                status = HTTP_OVERSIZE
                return
            end if
            if (rc == 0 .and. info_rc == 0) then
                if (code >= 200 .and. code < 300) then
                    call host_release(gate, 0.0_dp)
                    if (present(consumer)) then
                        call consumer([integer(int8) ::], .true., context, &
                            consume_status)
                        if (consume_status /= 0) return
                        allocate (bytes(0))
                    else
                        bytes = storage(:int(state%used))
                    end if
                    status = HTTP_OK
                    return
                end if
            end if
            if (state%consumer_error) then
                call host_release(gate, 0.0_dp)
                return
            end if
            if (present(consumer)) then
                if (rc /= 0 .and. state%delivered) then
                    call host_release(gate, 0.0_dp)
                    return
                end if
            end if
            ! Retry only transient transport failures, throttling and server faults.
            if (rc /= 0) then
                select case (rc)
                case (5, 6, 7, 18, 28, 52, 55, 56)
                case default
                    call host_release(gate, 0.0_dp)
                    return
                end select
            else if (code /= 429 .and. (code < 500 .or. code > 599)) then
                call host_release(gate, 0.0_dp)
                return
            end if
            call random_number(random)
            delay = max(retry, exponential*(1.0_dp + policy%jitter*random))
            call host_release(gate, delay)
            if (attempt >= budget) return
            if (http_clock() + delay >= end_time) return
            attempt = attempt + 1
            if (present(retries_left)) retries_left = retries_left - 1
            exponential = min(policy%max_backoff, exponential*2.0_dp)
        end do
    end subroutine
end module tiago_http_client
