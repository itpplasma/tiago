module tiago_recorded_diagnostic
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t
    use tiago_calibration, only: calibration_prior_t
    use tiago_recorded_affine, only: recorded_affine
    implicit none
    private

    type, abstract, extends(profile_diagnostic_t), public :: recorded_diagnostic_t
        type(calibration_prior_t) :: calibration
    contains
        procedure :: n_calibration => calibration_count
    end type

    type, extends(recorded_diagnostic_t), public :: calibrated_diagnostic_t
        class(profile_diagnostic_t), allocatable :: model
    contains
        procedure :: init => calibrated_init
        procedure :: n_signals => calibrated_count
        procedure :: n_nuisance => calibrated_nuisance
        procedure :: label => calibrated_label
        procedure :: sample_points => calibrated_points
        procedure :: predict => calibrated_predict
        procedure :: jvp => calibrated_jvp
        procedure :: vjp => calibrated_vjp
    end type
contains
    integer function calibration_count(self)
        class(recorded_diagnostic_t), intent(in) :: self
        calibration_count = self%calibration%n_params()
    end function

    subroutine calibrated_init(self, model, prior, status)
        class(calibrated_diagnostic_t), intent(inout) :: self
        class(profile_diagnostic_t), intent(in) :: model
        type(calibration_prior_t), intent(in) :: prior
        integer, intent(out) :: status
        status = 1
        if (prior%n_params() /= 2) return
        if (.not. allocated(prior%factor)) return
        if (allocated(self%model)) deallocate (self%model)
        allocate (self%model, source=model)
        self%calibration = prior
        status = 0
    end subroutine

    integer function calibrated_count(self)
        class(calibrated_diagnostic_t), intent(in) :: self
        calibrated_count = self%model%n_signals()
    end function

    integer function calibrated_nuisance(self)
        class(calibrated_diagnostic_t), intent(in) :: self
        calibrated_nuisance = self%model%n_nuisance() + self%n_calibration()
    end function

    function calibrated_label(self, i) result(text)
        class(calibrated_diagnostic_t), intent(in) :: self
        integer, intent(in) :: i
        character(:), allocatable :: text
        text = 'recorded:'//self%model%label(i)
    end function

    subroutine calibrated_points(self, eta, x, dx_deta)
        class(calibrated_diagnostic_t), intent(in) :: self
        real(dp), intent(in) :: eta(:)
        real(dp), allocatable, intent(out) :: x(:, :), dx_deta(:, :, :)
        real(dp), allocatable :: inner_dx(:, :, :)
        integer :: k
        k = self%model%n_nuisance()
        call self%model%sample_points(eta(:k), x, inner_dx)
        allocate (dx_deta(3, size(x, 2), self%n_nuisance()))
        dx_deta = 0.0_dp
        dx_deta(:, :, :k) = inner_dx
    end subroutine

    subroutine calibrated_predict(self, sample, eta, y)
        class(calibrated_diagnostic_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: c(2), raw(size(y)), dy, dg, db
        integer :: i, k
        k = self%model%n_nuisance()
        c = self%calibration%physical(eta(k + 1:))
        call self%model%predict(sample, eta(:k), raw)
        !$omp parallel do default(shared) private(i, dy, dg, db) if(size(y)>1024)
        do i = 1, size(y)
            call recorded_affine(raw(i), c(1), c(2), y(i), dy, dg, db)
        end do
        !$omp end parallel do
    end subroutine

    subroutine calibrated_jvp(self, sample, eta, q_dot, eta_dot, y_dot)
        class(calibrated_diagnostic_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        real(dp) :: c(2), dc(2), raw(size(y_dot)), dr(size(y_dot))
        real(dp) :: value, dy, dg, db
        integer :: i, k
        k = self%model%n_nuisance()
        c = self%calibration%physical(eta(k + 1:))
        dc = self%calibration%jvp(eta_dot(k + 1:))
        call self%model%predict(sample, eta(:k), raw)
        call self%model%jvp(sample, eta(:k), q_dot, eta_dot(:k), dr)
        !$omp parallel do default(shared) private(i, value, dy, dg, db) &
        !$omp if(size(y_dot)>1024)
        do i = 1, size(y_dot)
            call recorded_affine(raw(i), c(1), c(2), value, dy, dg, db)
            y_dot(i) = dy*dr(i) + dg*dc(1) + db*dc(2)
        end do
        !$omp end parallel do
    end subroutine

    subroutine calibrated_vjp(self, sample, eta, y_bar, q_bar, eta_bar)
        class(calibrated_diagnostic_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        real(dp) :: c(2), cbar(2), raw(size(y_bar)), rb(size(y_bar))
        real(dp) :: value, dy, dg, db
        integer :: i, k
        k = self%model%n_nuisance()
        c = self%calibration%physical(eta(k + 1:))
        call self%model%predict(sample, eta(:k), raw)
        cbar = 0.0_dp
        !$omp parallel do default(shared) private(i, value, dy, dg, db) &
        !$omp reduction(+:cbar) if(size(y_bar)>1024)
        do i = 1, size(y_bar)
            call recorded_affine(raw(i), c(1), c(2), value, dy, dg, db)
            rb(i) = dy*y_bar(i)
            cbar(1) = cbar(1) + dg*y_bar(i)
            cbar(2) = cbar(2) + db*y_bar(i)
        end do
        !$omp end parallel do
        call self%model%vjp(sample, eta(:k), rb, q_bar, eta_bar(:k))
        eta_bar(k + 1:) = self%calibration%vjp(cbar)
    end subroutine
end module tiago_recorded_diagnostic
