module tiago_label_profile
    !! Kinetic profiles f(rho) as seen by the profile diagnostics: a value, its
    !! rho-derivative, and linear derivatives with respect to the profile's own
    !! parameters. Latent parameterizations and priors (S0 of the inference
    !! plan) extend radial_profile_t; node_profile_t is the reference
    !! implementation (C1 cubic Hermite on node values, constant beyond the
    !! end nodes).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_radial_interp, only: hermite_weights, max_weights
    implicit none
    private

    type, abstract, public :: radial_profile_t
    contains
        procedure(n_params_i), deferred :: n_params
        procedure(get_params_i), deferred :: get_params
        procedure(set_params_i), deferred :: set_params
        procedure(eval_i), deferred :: eval
        procedure(param_jvp_i), deferred :: param_jvp
        procedure(param_vjp_i), deferred :: param_vjp
    end type radial_profile_t

    abstract interface
        pure integer function n_params_i(self)
            import :: radial_profile_t
            class(radial_profile_t), intent(in) :: self
        end function n_params_i

        pure subroutine get_params_i(self, p)
            import :: radial_profile_t, dp
            class(radial_profile_t), intent(in) :: self
            real(dp), intent(out) :: p(:)
        end subroutine get_params_i

        subroutine set_params_i(self, p)
            import :: radial_profile_t, dp
            class(radial_profile_t), intent(inout) :: self
            real(dp), intent(in) :: p(:)
        end subroutine set_params_i

        pure subroutine eval_i(self, rho, f, dfdrho)
            import :: radial_profile_t, dp
            class(radial_profile_t), intent(in) :: self
            real(dp), intent(in) :: rho
            real(dp), intent(out) :: f, dfdrho
        end subroutine eval_i

        pure subroutine param_jvp_i(self, rho, dp_, df)
            import :: radial_profile_t, dp
            class(radial_profile_t), intent(in) :: self
            real(dp), intent(in) :: rho, dp_(:)
            real(dp), intent(out) :: df
        end subroutine param_jvp_i

        pure subroutine param_vjp_i(self, rho, fbar, pbar)
            !! pbar += (df/dp)^T fbar
            import :: radial_profile_t, dp
            class(radial_profile_t), intent(in) :: self
            real(dp), intent(in) :: rho, fbar
            real(dp), intent(inout) :: pbar(:)
        end subroutine param_vjp_i
    end interface

    type, extends(radial_profile_t), public :: node_profile_t
        real(dp), allocatable :: nodes(:) !! increasing rho nodes
        real(dp), allocatable :: values(:)
    contains
        procedure :: n_params => node_n_params
        procedure :: get_params => node_get_params
        procedure :: set_params => node_set_params
        procedure :: eval => node_eval
        procedure :: param_jvp => node_param_jvp
        procedure :: param_vjp => node_param_vjp
    end type node_profile_t

contains

    pure integer function node_n_params(self)
        class(node_profile_t), intent(in) :: self
        node_n_params = size(self%values)
    end function node_n_params

    pure subroutine node_get_params(self, p)
        class(node_profile_t), intent(in) :: self
        real(dp), intent(out) :: p(:)
        p = self%values
    end subroutine node_get_params

    subroutine node_set_params(self, p)
        class(node_profile_t), intent(inout) :: self
        real(dp), intent(in) :: p(:)
        self%values = p
    end subroutine node_set_params

    pure subroutine node_eval(self, rho, f, dfdrho)
        class(node_profile_t), intent(in) :: self
        real(dp), intent(in) :: rho
        real(dp), intent(out) :: f, dfdrho
        integer :: first, nw
        real(dp) :: w(max_weights), dw(max_weights)

        call hermite_weights(self%nodes, rho, .false., first, nw, w, dw)
        f = dot_product(w(:nw), self%values(first:first + nw - 1))
        dfdrho = dot_product(dw(:nw), self%values(first:first + nw - 1))
    end subroutine node_eval

    pure subroutine node_param_jvp(self, rho, dp_, df)
        class(node_profile_t), intent(in) :: self
        real(dp), intent(in) :: rho, dp_(:)
        real(dp), intent(out) :: df
        integer :: first, nw
        real(dp) :: w(max_weights), dw(max_weights)

        call hermite_weights(self%nodes, rho, .false., first, nw, w, dw)
        df = dot_product(w(:nw), dp_(first:first + nw - 1))
    end subroutine node_param_jvp

    pure subroutine node_param_vjp(self, rho, fbar, pbar)
        class(node_profile_t), intent(in) :: self
        real(dp), intent(in) :: rho, fbar
        real(dp), intent(inout) :: pbar(:)
        integer :: first, nw
        real(dp) :: w(max_weights), dw(max_weights)

        call hermite_weights(self%nodes, rho, .false., first, nw, w, dw)
        pbar(first:first + nw - 1) = pbar(first:first + nw - 1) + fbar * w(:nw)
    end subroutine node_param_vjp
end module tiago_label_profile
