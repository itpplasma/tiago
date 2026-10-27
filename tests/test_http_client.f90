program test_http_client
    use, intrinsic :: iso_fortran_env, only: dp => real64, int8
    use tiago_http_client, only: http_policy_t, http_get, http_clock, &
        HTTP_OK, HTTP_OVERSIZE
    implicit none
    type(http_policy_t) :: policy
    integer(int8), allocatable :: bytes(:)
    character(2048) :: base_file, base
    integer :: unit, status, code, budget
    real(dp) :: start, elapsed
    call get_command_argument(1, base_file)
    open (newunit=unit, file=trim(base_file)//'/base.txt')
    read (unit, '(A)') base
    close (unit)
    policy%min_interval = 0.05_dp
    policy%backoff = 0.03_dp
    policy%jitter = 0.2_dp
    policy%max_elapsed = 10.0_dp
    policy%timeout = 2.0_dp
    policy%max_response_bytes = 2048
    call http_get(trim(base)//'/retry429', policy, bytes, status, code)
    call require(status == HTTP_OK .and. code == 200, '429 retry')
    ! The server event timestamps independently check the Retry-After floor.
    call http_get(trim(base)//'/date429', policy, bytes, status)
    call require(status == HTTP_OK, 'HTTP-date retry')
    start = http_clock()
    call http_get(trim(base)//'/retry500', policy, bytes, status)
    elapsed = http_clock() - start
    call require(status == HTTP_OK, '503 retry')
    call require(elapsed >= 0.09_dp, 'exponential backoff timing')
    print '(A,F10.6)', 'backoff_elapsed_s=', elapsed
    call http_get(trim(base)//'/drop', policy, bytes, status)
    call require(status == HTTP_OK, 'connection retry')
    budget = 1
    call http_get(trim(base)//'/always500', policy, bytes, status, &
        retries_left=budget)
    call require(status /= HTTP_OK .and. budget == 0, 'retry budget')
    call http_get(trim(base)//'/missing', policy, bytes, status, code)
    call require(status /= HTTP_OK .and. code == 404, 'no retry on 404')
    policy%min_interval = 0.0_dp
    call http_get(trim(base)//'/oversize', policy, bytes, status)
    call require(status == HTTP_OVERSIZE, 'hard response cap')
    policy%timeout = 0.05_dp
    policy%max_retries = 0
    start = http_clock()
    call http_get(trim(base)//'/slow', policy, bytes, status)
    elapsed = http_clock() - start
    call require(status /= HTTP_OK, 'timeout')
    print '(A,F10.6)', 'timeout_elapsed_s=', elapsed
    call require(elapsed < 0.25_dp, 'bounded timeout')
    policy%timeout = 2.0_dp
    policy%min_interval = 0.5_dp
    call http_get(trim(base)//'/interval', policy, bytes, status)
    call require(status == HTTP_OK, 'first paced request')
    start = http_clock()
    call http_get(trim(base)//'/interval', policy, bytes, status)
    elapsed = http_clock() - start
    call require(status == HTTP_OK, 'second paced request')
    call require(elapsed >= 0.49_dp, 'default minimum interval')
    print '(A,F10.6)', 'interval_elapsed_s=', elapsed
    policy%min_interval = 0.01_dp
    policy%max_elapsed = 0.02_dp
    call http_get(trim(base)//'/deadline', policy, bytes, status)
    call require(status /= HTTP_OK, 'pacing respects overall deadline')
    call http_get('file:///etc/hosts', policy, bytes, status)
    call require(status /= HTTP_OK, 'HTTP protocols only')
    policy%min_interval = -1.0_dp
    call http_get(trim(base), policy, bytes, status)
    call require(status /= HTTP_OK, 'invalid policy')
    print *, 'HTTP client controls passed'
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
