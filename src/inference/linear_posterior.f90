module tiago_linear_posterior
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_local_posterior, only: local_posterior_t
    use fortnum_rng, only: rng_t, rng_normal
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    implicit none
    private
    type, extends(local_posterior_t), public :: linear_posterior_t
        private
        real(dp), allocatable :: rows(:, :), rhs(:), prior_sd(:)
        real(dp) :: normalization = 0
    contains
        procedure :: configure
        procedure :: residual
        procedure :: jacobian
        procedure :: n_residuals
        procedure :: log_normalization
        procedure :: proper_prior
        procedure :: draw_prior
    end type
contains
    subroutine configure(self, design, data, log_likelihood_normalization, sd, status)
        class(linear_posterior_t), intent(inout) :: self
        real(dp), intent(in) :: design(:, :), data(:), sd(:)
        real(dp), intent(in) :: log_likelihood_normalization
        type(fortnum_status_t), intent(out) :: status
        integer :: m, n, i
        if (allocated(self%rows)) deallocate (self%rows, self%rhs, self%prior_sd)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid linear posterior')
        m = size(design, 1)
        n = size(design, 2)
        if (m < 1 .or. n < 1) return
        if (size(data) /= m .or. size(sd) /= n) return
        if (.not. all(ieee_is_finite(design))) return
        if (.not. all(ieee_is_finite(data))) return
        if (.not. all(ieee_is_finite(sd))) return
        if (.not. ieee_is_finite(log_likelihood_normalization)) return
        if (any(sd <= 0)) return
        allocate (self%rows(m + n, n), self%rhs(m + n))
        self%prior_sd = sd
        self%rows = 0
        self%rhs = 0
        self%rows(:m, :) = design
        self%rhs(:m) = data
        do i = 1, n
            self%rows(m + i, i) = 1/sd(i)
        end do
        self%normalization = log_likelihood_normalization &
            - 0.5_dp*n*log(2*acos(-1.0_dp)) - sum(log(sd))
        call status_set(status, FORTNUM_OK, '')
    end subroutine
    subroutine residual(self, theta, r, ok)
        class(linear_posterior_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), intent(out) :: r(:)
        logical, intent(out) :: ok
        ok = allocated(self%rows)
        if (.not. ok) return
        ok = size(theta) == size(self%rows, 2) .and. size(r) == size(self%rhs)
        if (.not. ok) return
        r = matmul(self%rows, theta) - self%rhs
    end subroutine
    subroutine jacobian(self, theta, jac)
        class(linear_posterior_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), allocatable, intent(out) :: jac(:, :)
        if (.not. allocated(self%rows)) error stop 'unconfigured linear posterior'
        if (size(theta) /= size(self%rows, 2)) error stop 'linear parameter count'
        jac = self%rows
    end subroutine
    integer function n_residuals(self) result(n)
        class(linear_posterior_t), intent(in) :: self
        n = 0
        if (allocated(self%rhs)) n = size(self%rhs)
    end function
    real(dp) function log_normalization(self) result(c)
        class(linear_posterior_t), intent(in) :: self
        c = self%normalization
    end function
    logical function proper_prior(self) result(proper)
        class(linear_posterior_t), intent(in) :: self
        proper = allocated(self%prior_sd)
        if (proper) proper = all(self%prior_sd > 0)
    end function
    subroutine draw_prior(self, rng, theta, ok)
        class(linear_posterior_t), intent(in) :: self
        type(rng_t), intent(inout) :: rng
        real(dp), intent(out) :: theta(:)
        logical, intent(out) :: ok
        integer :: i
        ok = self%proper_prior()
        if (.not. ok) return
        ok = size(theta) == size(self%prior_sd)
        if (.not. ok) return
        do i = 1, size(theta)
            call rng_normal(rng, theta(i))
        end do
        theta = self%prior_sd*theta
    end subroutine
end module
