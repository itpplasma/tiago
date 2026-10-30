module tiago_local_posterior
    !! Posterior contract in the coordinates used by MAP and Laplace.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_rng, only: rng_t
    implicit none
    private

    type, abstract, public :: local_posterior_t
    contains
        procedure(residual_interface), deferred :: residual
        procedure(jacobian_interface), deferred :: jacobian
        procedure(count_interface), deferred :: n_residuals
        procedure(normalization_interface), deferred :: log_normalization
        procedure(proper_interface), deferred :: proper_prior
        procedure(draw_interface), deferred :: draw_prior
        procedure :: exact => posterior_exact
        procedure :: save_equilibrium => posterior_save
    end type

    abstract interface
        subroutine residual_interface(self, theta, r, ok)
            import local_posterior_t, dp
            class(local_posterior_t), intent(inout) :: self
            real(dp), intent(in) :: theta(:)
            real(dp), intent(out) :: r(:)
            logical, intent(out) :: ok
        end subroutine
        subroutine jacobian_interface(self, theta, jac)
            import local_posterior_t, dp
            class(local_posterior_t), intent(inout) :: self
            real(dp), intent(in) :: theta(:)
            real(dp), allocatable, intent(out) :: jac(:, :)
        end subroutine
        integer function count_interface(self)
            import local_posterior_t
            class(local_posterior_t), intent(in) :: self
        end function
        real(dp) function normalization_interface(self)
            import local_posterior_t, dp
            class(local_posterior_t), intent(in) :: self
        end function
        logical function proper_interface(self)
            import local_posterior_t
            class(local_posterior_t), intent(in) :: self
        end function
        subroutine draw_interface(self, rng, theta, ok)
            import local_posterior_t, dp, rng_t
            class(local_posterior_t), intent(in) :: self
            type(rng_t), intent(inout) :: rng
            real(dp), intent(out) :: theta(:)
            logical, intent(out) :: ok
        end subroutine
    end interface
contains
    subroutine posterior_exact(self)
        class(local_posterior_t), intent(inout) :: self
    end subroutine

    subroutine posterior_save(self, prefix, path)
        class(local_posterior_t), intent(inout) :: self
        character(len=*), intent(in) :: prefix
        character(len=:), allocatable, intent(out) :: path
        path = ''
    end subroutine
end module
