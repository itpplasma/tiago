program test_http_parallel
    use, intrinsic :: iso_fortran_env, only: int8
    use tiago_http_client, only: http_policy_t, http_get, HTTP_OK
    implicit none
    type(http_policy_t) :: policy
    integer(int8), allocatable :: bytes(:)
    integer :: i, unit, status(4)
    character(2048) :: directory, base
    call get_command_argument(1, directory)
    open (newunit=unit, file=trim(directory)//'/base.txt')
    read (unit, '(A)') base
    close (unit)
    policy%min_interval = 0.01d0
    !$omp parallel do private(bytes) num_threads(4)
    do i = 1, 4
        call http_get(trim(base)//'/parallel', policy, bytes, status(i))
    end do
    !$omp end parallel do
    if (any(status /= HTTP_OK)) stop 1
    print *, 'Four concurrent callers completed under the shared host gate'
end program
