module tiago_profile_tape
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t
    implicit none
    private

    type, abstract, public :: profile_tape_t
        logical :: valid = .false.
    contains
        procedure(tangent_i), deferred :: signal_jvp
        procedure(adjoint_i), deferred :: signal_vjp
    end type

    abstract interface
        subroutine tangent_i(self, parameter_dot, eta_dot, signal_dot, status)
            import profile_tape_t, dp, fortnum_status_t
            class(profile_tape_t), intent(in) :: self
            real(dp), intent(in) :: parameter_dot(:), eta_dot(:)
            real(dp), intent(out) :: signal_dot(:)
            type(fortnum_status_t), intent(out) :: status
        end subroutine

        subroutine adjoint_i(self, signal_bar, parameter_bar, eta_bar, status)
            import profile_tape_t, dp, fortnum_status_t
            class(profile_tape_t), intent(in) :: self
            real(dp), intent(in) :: signal_bar(:)
            real(dp), intent(out) :: parameter_bar(:), eta_bar(:)
            type(fortnum_status_t), intent(out) :: status
        end subroutine
    end interface
end module
