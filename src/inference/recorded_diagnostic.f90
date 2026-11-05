module tiago_recorded_diagnostic
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t
    use tiago_calibration, only: calibration_prior_t
    use tiago_recorded_affine, only: recorded_affine
    use tiago_recorded_log_gain, only: recorded_log_gain
    implicit none
    private

    type, abstract, extends(profile_diagnostic_t), public :: recorded_diagnostic_t
        type(calibration_prior_t) :: calibration
    contains
        procedure :: n_calibration => calibration_count
    end type

    type, extends(recorded_diagnostic_t), public :: calibrated_diagnostic_t
        class(profile_diagnostic_t), allocatable :: model
        integer, allocatable, private :: groups(:)
        integer, private :: n_groups = 1
        logical, private :: log_gains = .false., gain_only = .false.
    contains
        procedure :: init => calibrated_init
        procedure :: n_signals => calibrated_count
        procedure :: n_nuisance => calibrated_nuisance
        procedure :: label => calibrated_label
        procedure :: sample_points => calibrated_points
        procedure :: predict => calibrated_predict
        procedure :: jvp => calibrated_jvp
        procedure :: vjp => calibrated_vjp
        procedure :: response_predict
        procedure :: response_jvp
        procedure :: response_vjp
        procedure :: uses_low_fidelity_te => calibrated_low_fidelity_te
    end type
contains
    integer function calibration_count(self)
        class(recorded_diagnostic_t), intent(in) :: self
        calibration_count = self%calibration%n_params()
    end function

    subroutine calibrated_init(self, model, prior, status, groups, log_gains)
        class(calibrated_diagnostic_t), intent(inout) :: self
        class(profile_diagnostic_t), intent(in) :: model
        type(calibration_prior_t), intent(in) :: prior
        integer, intent(out) :: status
        integer, intent(in), optional :: groups(:)
        logical, intent(in), optional :: log_gains
        integer :: ng, g
        status = 1
        ng = 1
        if (present(groups)) then
            if (size(groups) /= model%n_signals()) return
            if (size(groups) == 0) return
            if (any(groups < 1)) return
            ng = maxval(groups)
            if (ng > model%n_signals()) return
            do g = 1, ng
                if (.not. any(groups == g)) return
            end do
        end if
        if (prior%n_params() /= ng .and. prior%n_params() /= 2*ng) return
        if (.not. allocated(prior%factor)) return
        if (size(prior%factor, 1) /= prior%n_params()) return
        if (size(prior%factor, 2) /= prior%n_params()) return
        if (allocated(self%model)) deallocate (self%model)
        allocate (self%model, source=model)
        self%calibration = prior
        self%n_groups = ng
        self%gain_only = prior%n_params() == ng
        if (allocated(self%groups)) deallocate (self%groups)
        allocate (self%groups(model%n_signals()))
        self%groups = 1
        if (present(groups)) self%groups = groups
        self%log_gains = .false.
        if (present(log_gains)) self%log_gains = log_gains
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

    logical function calibrated_low_fidelity_te(self)
        class(calibrated_diagnostic_t), intent(in) :: self
        calibrated_low_fidelity_te = self%model%uses_low_fidelity_te()
    end function

    subroutine calibrated_predict(self, sample, eta, y)
        class(calibrated_diagnostic_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: raw(size(y))
        integer :: k
        k = self%model%n_nuisance()
        call self%model%predict(sample, eta(:k), raw)
        call self%response_predict(raw, eta(k + 1:), y)
    end subroutine

    subroutine calibrated_jvp(self, sample, eta, q_dot, eta_dot, y_dot)
        class(calibrated_diagnostic_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        real(dp) :: raw(size(y_dot)), dr(size(y_dot))
        integer :: k
        k = self%model%n_nuisance()
        call self%model%predict(sample, eta(:k), raw)
        call self%model%jvp(sample, eta(:k), q_dot, eta_dot(:k), dr)
        call self%response_jvp(raw, eta(k + 1:), dr, eta_dot(k + 1:), y_dot)
    end subroutine

    subroutine calibrated_vjp(self, sample, eta, y_bar, q_bar, eta_bar)
        class(calibrated_diagnostic_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        real(dp) :: raw(size(y_bar)), rb(size(y_bar))
        integer :: k
        k = self%model%n_nuisance()
        call self%model%predict(sample, eta(:k), raw)
        call self%response_vjp(raw, eta(k + 1:), y_bar, rb, eta_bar(k + 1:))
        call self%model%vjp(sample, eta(:k), rb, q_bar, eta_bar(:k))
    end subroutine

    subroutine response_predict(self, raw, z, y)
        class(calibrated_diagnostic_t), intent(in) :: self
        real(dp), intent(in) :: raw(:), z(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: c(2*self%n_groups), slope(self%n_groups), dy, dg, db
        integer :: i, g
        call calibration_values(self, z, c, slope)
        !$omp parallel do default(shared) private(i, g, dy, dg, db) if(size(y)>1024)
        do i = 1, size(y)
            g = 2*self%groups(i) - 1
            call recorded_affine(raw(i), c(g), c(g + 1), y(i), dy, dg, db)
        end do
        !$omp end parallel do
    end subroutine

    subroutine response_jvp(self, raw, z, raw_dot, z_dot, y_dot)
        class(calibrated_diagnostic_t), intent(in) :: self
        real(dp), intent(in) :: raw(:), z(:), raw_dot(:), z_dot(:)
        real(dp), intent(out) :: y_dot(:)
        real(dp) :: c(2*self%n_groups), dc(2*self%n_groups), slope(self%n_groups)
        real(dp) :: value, dy, dg, db, dc_prior(size(self%calibration%mean))
        integer :: i, g
        call calibration_values(self, z, c, slope)
        dc_prior = self%calibration%jvp(z_dot)
        dc = 0
        if (self%gain_only) then
            dc(1::2) = dc_prior
        else
            dc = dc_prior
        end if
        dc(1::2) = slope*dc(1::2)
        !$omp parallel do default(shared) private(i, g, value, dy, dg, db) &
        !$omp if(size(y_dot)>1024)
        do i = 1, size(y_dot)
            g = 2*self%groups(i) - 1
            call recorded_affine(raw(i), c(g), c(g + 1), value, dy, dg, db)
            y_dot(i) = dy*raw_dot(i) + dg*dc(g) + db*dc(g + 1)
        end do
        !$omp end parallel do
    end subroutine

    subroutine response_vjp(self, raw, z, y_bar, raw_bar, z_bar)
        class(calibrated_diagnostic_t), intent(in) :: self
        real(dp), intent(in) :: raw(:), z(:), y_bar(:)
        real(dp), intent(out) :: raw_bar(:), z_bar(:)
        real(dp) :: c(2*self%n_groups), cbar(2*self%n_groups), slope(self%n_groups)
        real(dp) :: value, dy, dg, db
        integer :: i, g
        call calibration_values(self, z, c, slope)
        cbar = 0
        !$omp parallel do default(shared) private(i, g, value, dy, dg, db) &
        !$omp reduction(+:cbar) if(size(y_bar)>1024)
        do i = 1, size(y_bar)
            g = 2*self%groups(i) - 1
            call recorded_affine(raw(i), c(g), c(g + 1), value, dy, dg, db)
            raw_bar(i) = dy*y_bar(i)
            cbar(g) = cbar(g) + dg*y_bar(i)
            cbar(g + 1) = cbar(g + 1) + db*y_bar(i)
        end do
        !$omp end parallel do
        cbar(1::2) = slope*cbar(1::2)
        if (self%gain_only) then
            z_bar = self%calibration%vjp(cbar(1::2))
        else
            z_bar = self%calibration%vjp(cbar)
        end if
    end subroutine

    subroutine calibration_values(self, z, c, slope)
        class(calibrated_diagnostic_t), intent(in) :: self
        real(dp), intent(in) :: z(:)
        real(dp), intent(out) :: c(:), slope(:)
        real(dp) :: gain, c_prior(size(z))
        integer :: g
        c_prior = self%calibration%physical(z)
        c = 0
        if (self%gain_only) then
            c(1::2) = c_prior
        else
            c = c_prior
        end if
        slope = 1.0_dp
        if (.not. self%log_gains) return
        do g = 1, size(slope)
            call recorded_log_gain(c(2*g - 1), gain, slope(g))
            c(2*g - 1) = gain
        end do
    end subroutine
end module tiago_recorded_diagnostic
