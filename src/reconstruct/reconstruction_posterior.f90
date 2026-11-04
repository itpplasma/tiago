module tiago_reconstruction_posterior
    !! Magnetic reconstruction in scaled coordinates, using the shared engines.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_local_posterior, only: local_posterior_t
    use tiago_reconstruction, only: rec, residual_at, jacobian_at, n_rows, &
        invalidate_caches
    use fortnum_rng, only: rng_t, rng_normal
    implicit none
    private
    type, extends(local_posterior_t), public :: reconstruction_posterior_t
    contains
        procedure :: residual => magnetic_residual
        procedure :: jacobian => magnetic_jacobian
        procedure :: n_residuals => magnetic_rows
        procedure :: log_normalization => magnetic_normalization
        procedure :: proper_prior => magnetic_proper
        procedure :: draw_prior => magnetic_draw
        procedure :: exact => magnetic_exact
        procedure :: save_equilibrium => magnetic_save
    end type
contains
    subroutine magnetic_residual(self, theta, r, ok)
        class(reconstruction_posterior_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), intent(out) :: r(:)
        logical, intent(out) :: ok
        call residual_at(theta*rec%scale, r)
        ok = all(ieee_is_finite(r))
    end subroutine

    subroutine magnetic_jacobian(self, theta, jac)
        class(reconstruction_posterior_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), allocatable, intent(out) :: jac(:, :)
        integer :: i
        call jacobian_at(theta*rec%scale, jac)
        do i = 1, rec%n
            jac(:, i) = jac(:, i)*rec%scale(i)
        end do
    end subroutine

    integer function magnetic_rows(self) result(n)
        class(reconstruction_posterior_t), intent(in) :: self
        n = n_rows()
    end function

    logical function magnetic_proper(self) result(proper)
        class(reconstruction_posterior_t), intent(in) :: self
        proper = rec%n_prior == rec%n
    end function

    real(dp) function magnetic_normalization(self) result(c)
        class(reconstruction_posterior_t), intent(in) :: self
        integer :: i, k
        c = -0.5_dp*n_rows()*log(2.0_dp*acos(-1.0_dp)) - &
            sum(log(rec%meas_sigma)) - rec%n_fb*log(rec%fb_sigma)
        do i = 1, rec%n_prior
            k = rec%prior_param(i)
            c = c - log(rec%prior_sigma(k)/rec%scale(k))
        end do
    end function

    subroutine magnetic_draw(self, rng, theta, ok)
        class(reconstruction_posterior_t), intent(in) :: self
        type(rng_t), intent(inout) :: rng
        real(dp), intent(out) :: theta(:)
        logical, intent(out) :: ok
        integer :: i
        ok = self%proper_prior() .and. size(theta) == rec%n
        if (.not. ok) return
        do i = 1, rec%n
            call rng_normal(rng, theta(i))
        end do
        theta = (rec%prior + rec%prior_sigma*theta)/rec%scale
    end subroutine

    subroutine magnetic_exact(self)
        class(reconstruction_posterior_t), intent(inout) :: self
        call invalidate_caches()
    end subroutine

    subroutine magnetic_save(self, prefix, path)
        class(reconstruction_posterior_t), intent(inout) :: self
        character(len=*), intent(in) :: prefix
        character(len=:), allocatable, intent(out) :: path
        path = prefix//'.nc'
        call rec%eq%write_wout(path)
    end subroutine
end module
