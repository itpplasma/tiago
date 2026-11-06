program tiago_linear_infer
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_linear_posterior, only: linear_posterior_t
    use tiago_inference_engines, only: laplace, laplace_t
    use fortnum_status, only: fortnum_status_t
    implicit none
    type :: case_t
        type(linear_posterior_t) :: model
        type(laplace_t) :: result
        logical :: converged
        integer :: parameters
    end type
    type(case_t), allocatable :: batch(:)
    type(fortnum_status_t) :: status
    real(dp), allocatable :: a(:, :), y(:), sd(:)
    real(dp) :: normalization
    character(len=4096) :: input, output
    integer :: u, v, cases, k, m, n, i
    call get_command_argument(1, input)
    call get_command_argument(2, output)
    open (newunit=u, file=trim(input), status='old', action='read')
    read (u, *) cases
    if (cases < 1) error stop 'empty linear batch'
    allocate (batch(cases))
    do k = 1, cases
        read (u, *) m, n, normalization
        if (m < 1 .or. n < 1) error stop 'invalid linear dimensions'
        batch(k)%parameters = n
        allocate (a(m, n), y(m), sd(n))
        read (u, *) sd
        do i = 1, m
            read (u, *) y(i), a(i, :)
        end do
        call batch(k)%model%configure(a, y, normalization, sd, status)
        if (status%code /= 0) error stop 'invalid linear model'
        deallocate (a, y, sd)
    end do
    close (u)
    !$omp parallel do default(shared) private(k) schedule(dynamic)
    do k = 1, cases
        call fit_case(batch(k))
    end do
    !$omp end parallel do
    open (newunit=v, file=trim(output), status='replace', action='write')
    do k = 1, cases
        if (.not. batch(k)%converged) error stop 'linear Gaussian solve failed'
        if (batch(k)%result%info /= 0 .or. .not. batch(k)%result%proper) &
            error stop 'linear Laplace failed'
        n = batch(k)%parameters
        write (v, *) k, n, batch(k)%result%log_evidence, batch(k)%result%chi2
        write (v, *) batch(k)%result%mean
        do i = 1, n
            write (v, *) batch(k)%result%cov(i, :)
        end do
    end do
    close (v)
contains
    recursive subroutine fit_case(item)
        type(case_t), intent(inout) :: item
        real(dp), allocatable :: theta(:), residual(:), jac(:, :), precision(:, :)
        logical :: ok
        integer :: info
        external :: dposv
        allocate(theta(item%parameters), residual(item%model%n_residuals()))
        theta = 0
        call item%model%residual(theta, residual, ok)
        call item%model%jacobian(theta, jac)
        precision = matmul(transpose(jac), jac)
        theta = -matmul(transpose(jac), residual)
        call dposv('U', size(theta), 1, precision, size(theta), theta, size(theta), info)
        item%converged = ok .and. info == 0
        if (item%converged) call laplace(item%model, theta, item%result)
    end subroutine
end program
