module tiago_calibrated_tape
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_DOMAIN_ERROR
    use tiago_profile_tape, only: profile_tape_t
    use tiago_recorded_diagnostic, only: calibrated_diagnostic_t
    implicit none
    private

    type, extends(profile_tape_t), public :: calibrated_tape_t
        private
        class(profile_tape_t), allocatable :: inner
        type(calibrated_diagnostic_t) :: response
        real(dp), allocatable :: raw(:), eta(:)
    contains
        procedure :: init => calibrated_record
        procedure :: signal_jvp => calibrated_tangent
        procedure :: signal_vjp => calibrated_adjoint
    end type
contains
    subroutine calibrated_record(self, diagnostic, eta, raw, inner)
        class(calibrated_tape_t), intent(out) :: self
        type(calibrated_diagnostic_t), intent(in) :: diagnostic
        real(dp), intent(in) :: eta(:), raw(:)
        class(profile_tape_t), allocatable, intent(inout) :: inner
        if (.not. allocated(inner)) return
        if (.not. inner%valid) return
        if (size(eta) /= diagnostic%n_nuisance()) return
        if (size(raw) /= diagnostic%n_signals()) return
        self%response = diagnostic
        self%raw = raw
        self%eta = eta
        call move_alloc(inner, self%inner)
        self%valid = .true.
    end subroutine

    subroutine calibrated_tangent(self, parameter_dot, eta_dot, signal_dot, status)
        class(calibrated_tape_t), intent(in) :: self
        real(dp), intent(in) :: parameter_dot(:), eta_dot(:)
        real(dp), intent(out) :: signal_dot(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: raw_dot(:)
        integer :: k
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid calibrated tape tangent')
        if (.not. self%valid) return
        if (size(eta_dot) /= size(self%eta)) return
        if (size(signal_dot) /= size(self%raw)) return
        k = size(self%eta) - self%response%n_calibration()
        allocate (raw_dot(size(self%raw)))
        call self%inner%signal_jvp(parameter_dot, eta_dot(:k), raw_dot, status)
        if (status%code /= 0) return
        call self%response%response_jvp(self%raw, self%eta(k + 1:), raw_dot, &
            eta_dot(k + 1:), signal_dot)
    end subroutine

    subroutine calibrated_adjoint(self, signal_bar, parameter_bar, eta_bar, status)
        class(calibrated_tape_t), intent(in) :: self
        real(dp), intent(in) :: signal_bar(:)
        real(dp), intent(out) :: parameter_bar(:), eta_bar(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: raw_bar(:)
        integer :: k
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid calibrated tape adjoint')
        if (.not. self%valid) return
        if (size(eta_bar) /= size(self%eta)) return
        if (size(signal_bar) /= size(self%raw)) return
        k = size(self%eta) - self%response%n_calibration()
        allocate (raw_bar(size(self%raw)))
        call self%response%response_vjp(self%raw, self%eta(k + 1:), signal_bar, &
            raw_bar, eta_bar(k + 1:))
        call self%inner%signal_vjp(raw_bar, parameter_bar, eta_bar(:k), status)
    end subroutine
end module
