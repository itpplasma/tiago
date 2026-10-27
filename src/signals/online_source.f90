module tiago_online_source
    use, intrinsic :: iso_fortran_env, only: dp => real64, int8
    use tiago_signal_source, only: signal_source_t, signal_adapter_t, &
        signal_request_t, signal_window_t, backend_request_t, &
        SIGNAL_OK, SIGNAL_ERROR, SIGNAL_MORE, version_matches, content_checksum
    use tiago_signal_reduction, only: reduction_step_t
    use tiago_http_client, only: http_policy_t, http_get, http_clock, HTTP_OK
    implicit none
    private
    type :: stream_context_t
        class(signal_adapter_t), pointer :: adapter
        class(reduction_step_t), pointer :: reduction
        type(signal_request_t) :: request
        type(signal_window_t) :: joined
        logical :: first = .true.
        integer :: total = 0, max_output_samples
    end type
    type, extends(signal_source_t), public :: online_source_t
        class(signal_adapter_t), pointer :: adapter => null()
        class(reduction_step_t), pointer :: reduction => null()
        type(http_policy_t) :: policy
    contains
        procedure :: read => online_read
    end type
contains
    subroutine online_read(self, request, window, status)
        class(online_source_t), intent(inout) :: self
        type(signal_request_t), intent(in) :: request
        type(signal_window_t), intent(out) :: window
        integer, intent(out) :: status
        type(signal_request_t) :: part
        type(backend_request_t) :: backend
        type(signal_window_t) :: discarded
        type(stream_context_t), target :: context
        integer(int8), allocatable :: bytes(:)
        integer :: retries, nwin, k, ignored
        real(dp) :: deadline, right
        status = SIGNAL_ERROR
        if (.not. request%valid()) return
        if (.not. associated(self%adapter)) return
        if (.not. self%policy%valid()) return
        deadline = http_clock() + self%policy%max_elapsed
        retries = self%policy%max_retries
        nwin = 0
        context%adapter => self%adapter
        context%reduction => self%reduction
        context%max_output_samples = self%policy%max_output_samples
        part = request
        do
            nwin = nwin + 1
            if (nwin > self%policy%max_windows) return
            if (request%last - part%first > self%policy%max_window) then
                right = part%first + self%policy%max_window
            else
                right = request%last
            end if
            if (right < part%first) return
            part%last = right
            call self%adapter%map_request(part, backend, k)
            if (k /= SIGNAL_OK) return
            if (.not. allocated(backend%method)) return
            if (backend%method /= 'GET') return
            if (.not. allocated(backend%endpoint)) return
            context%request = part
            call http_get(backend%endpoint, self%policy, bytes, k, &
                retries_left=retries, deadline=deadline, &
                consumer=consume, context=context)
            if (k /= HTTP_OK) then
                ! Finalize an interrupted decoder; its output is discarded.
                call self%adapter%decode(part, [integer(int8) ::], .true., &
                    discarded, ignored)
                return
            end if
            part%version = context%request%version
            deallocate (bytes)
            if (right == request%last) exit
            ! Closed intervals split without requesting a boundary twice.
            part%first = nearest(right, 1.0_dp)
            if (part%first <= right) return
        end do
        if (context%first) return
        if (.not. context%joined%valid()) return
        context%joined%provenance%request_key = request%key()
        call context%joined%seal()
        window = context%joined
        status = SIGNAL_OK
    end subroutine

    subroutine consume(bytes, final, state, status)
        integer(int8), intent(in) :: bytes(:)
        logical, intent(in) :: final
        class(*), intent(inout) :: state
        integer, intent(out) :: status
        type(signal_window_t) :: raw, reduced
        integer :: k
        status = SIGNAL_ERROR
        select type (state)
            type is (stream_context_t)
            call state%adapter%decode(state%request, bytes, final, raw, k)
            if (k == SIGNAL_MORE) then
                if (.not. final) status = SIGNAL_OK
                return
            end if
            if (k /= SIGNAL_OK) return
            if (.not. raw%valid()) return
            if (.not. raw%has_provenance()) return
            if (.not. version_matches(state%request, raw%provenance%version)) return
            if (content_checksum(raw) /= raw%provenance%checksum) return
            if (any(raw%time < state%request%first)) return
            if (any(raw%time > state%request%last)) return
            state%request%version = raw%provenance%version
            if (associated(state%reduction)) then
                call state%reduction%apply(raw, reduced, k)
                if (k /= SIGNAL_OK) return
            else
                reduced = raw
            end if
            if (.not. reduced%valid()) return
            if (.not. reduced%has_provenance()) return
            if (size(reduced%samples) > &
                state%max_output_samples - state%total) return
            state%total = state%total + size(reduced%samples)
            call append(state%joined, reduced, state%first, k)
            if (k /= SIGNAL_OK) return
            state%first = .false.
            status = SIGNAL_OK
        end select
    end subroutine

    subroutine append(joined, next, first, status)
        type(signal_window_t), intent(inout) :: joined
        type(signal_window_t), intent(in) :: next
        logical, intent(in) :: first
        integer, intent(out) :: status
        real(dp), allocatable :: covariance(:, :)
        integer :: n, m
        status = SIGNAL_ERROR
        if (first) then
            joined = next
            status = SIGNAL_OK
            return
        end if
        if (joined%units /= next%units) return
        if (joined%provenance%version /= next%provenance%version) return
        if (joined%provenance%source_id /= next%provenance%source_id) return
        n = size(joined%time)
        m = size(next%time)
        if (n > 0 .and. m > 0) then
            if (next%time(1) < joined%time(n)) return
        end if
        if (allocated(joined%covariance) .or. allocated(next%covariance)) then
            allocate (covariance(n + m, n + m))
            covariance = 0.0_dp
            if (allocated(joined%covariance)) covariance(:n, :n) = joined%covariance
            if (allocated(next%covariance)) covariance(n + 1:, n + 1:) = next%covariance
            call move_alloc(covariance, joined%covariance)
        end if
        joined%time = [joined%time, next%time]
        joined%samples = [joined%samples, next%samples]
        joined%variance = [joined%variance, next%variance]
        joined%valid_from = max(joined%valid_from, next%valid_from)
        joined%valid_to = min(joined%valid_to, next%valid_to)
        joined%provenance%metadata = joined%provenance%metadata// &
            new_line('a')//'window_checksum='//next%provenance%checksum
        status = SIGNAL_OK
    end subroutine
end module tiago_online_source
