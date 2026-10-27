module tiago_integrator_diagnostic
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t
    use tiago_recorded_diagnostic, only: recorded_diagnostic_t
    use tiago_calibration, only: calibration_prior_t
    use tiago_integrator_signal, only: integrator_signal
    implicit none
    private
    public :: integrator_diagnostic_t

    type, extends(recorded_diagnostic_t), public :: integrator_diagnostic_t
        class(profile_diagnostic_t), allocatable :: model
        real(dp), allocatable :: time(:)
        real(dp) :: origin = 0.0_dp
    contains
        procedure :: init => integrator_init
        procedure :: n_signals => integrator_count
        procedure :: n_nuisance => integrator_nuisance
        procedure :: label => integrator_label
        procedure :: sample_points => integrator_points
        procedure :: predict => integrator_predict
        procedure :: jvp => integrator_jvp
        procedure :: vjp => integrator_vjp
    end type
contains
    subroutine integrator_init(self, model, prior, time, origin, status)
        class(integrator_diagnostic_t), intent(inout) :: self
        class(profile_diagnostic_t), intent(in) :: model
        type(calibration_prior_t), intent(in) :: prior
        real(dp), intent(in) :: time(:), origin
        integer, intent(out) :: status
        status = 1
        if (prior%n_params() /= 3) return
        if (size(time) /= model%n_signals()) return
        if (.not. all(ieee_is_finite(time))) return
        if (.not. ieee_is_finite(origin)) return
        if (any(time < prior%valid_from)) return
        if (any(time > prior%valid_to)) return
        if (.not. allocated(prior%factor)) return
        if (allocated(self%model)) deallocate (self%model)
        allocate (self%model, source=model)
        self%calibration = prior
        self%time = time
        self%origin = origin
        status = 0
    end subroutine

    integer function integrator_count(self)
        class(integrator_diagnostic_t), intent(in) :: self
        integrator_count = self%model%n_signals()
    end function

    integer function integrator_nuisance(self)
        class(integrator_diagnostic_t), intent(in) :: self
        integrator_nuisance = self%model%n_nuisance() + self%n_calibration()
    end function

    function integrator_label(self, i) result(text)
        class(integrator_diagnostic_t), intent(in) :: self
        integer, intent(in) :: i
        character(:), allocatable :: text
        text = 'recorded:'//self%model%label(i)
    end function

    subroutine integrator_points(self, eta, x, dx_deta)
        class(integrator_diagnostic_t), intent(in) :: self
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

    subroutine integrator_predict(self, sample, eta, y)
        class(integrator_diagnostic_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: c(3), dy, dg, db, dd
        real(dp), allocatable :: raw(:)
        integer :: i, k
        allocate (raw(size(y)))
        k = self%model%n_nuisance()
        c = self%calibration%physical(eta(k + 1:))
        call self%model%predict(sample, eta(:k), raw)
        !$omp parallel do default(shared) private(i, dy, dg, db, dd) if(size(y)>1024)
        do i = 1, size(y)
            call integrator_signal(raw(i), c(1), c(2), c(3), self%time(i), &
                                   self%origin, y(i), dy, dg, db, dd)
        end do
        !$omp end parallel do
    end subroutine

    subroutine integrator_jvp(self, sample, eta, q_dot, eta_dot, y_dot)
        class(integrator_diagnostic_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        real(dp) :: c(3), dc(3)
        real(dp), allocatable :: raw(:), dr(:)
        real(dp) :: value, dy, dg, db, dd
        integer :: i, k
        allocate (raw(size(y_dot)), dr(size(y_dot)))
        k = self%model%n_nuisance()
        c = self%calibration%physical(eta(k + 1:))
        dc = self%calibration%jvp(eta_dot(k + 1:))
        call self%model%predict(sample, eta(:k), raw)
        call self%model%jvp(sample, eta(:k), q_dot, eta_dot(:k), dr)
        !$omp parallel do default(shared) private(i, value, dy, dg, db, dd) &
        !$omp if(size(y_dot)>1024)
        do i = 1, size(y_dot)
            call integrator_signal(raw(i), c(1), c(2), c(3), self%time(i), &
                                   self%origin, value, dy, dg, db, dd)
            y_dot(i) = dy*dr(i) + dg*dc(1) + db*dc(2) + dd*dc(3)
        end do
        !$omp end parallel do
    end subroutine

    subroutine integrator_vjp(self, sample, eta, y_bar, q_bar, eta_bar)
        class(integrator_diagnostic_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        real(dp) :: c(3), cbar(3)
        real(dp), allocatable :: raw(:), rb(:)
        real(dp) :: value, dy, dg, db, dd
        integer :: i, k
        allocate (raw(size(y_bar)), rb(size(y_bar)))
        k = self%model%n_nuisance()
        c = self%calibration%physical(eta(k + 1:))
        call self%model%predict(sample, eta(:k), raw)
        cbar = 0.0_dp
        !$omp parallel do default(shared) private(i, value, dy, dg, db, dd) &
        !$omp reduction(+:cbar) if(size(y_bar)>1024)
        do i = 1, size(y_bar)
            call integrator_signal(raw(i), c(1), c(2), c(3), self%time(i), &
                                   self%origin, value, dy, dg, db, dd)
            rb(i) = dy*y_bar(i)
            cbar(1) = cbar(1) + dg*y_bar(i)
            cbar(2) = cbar(2) + db*y_bar(i)
            cbar(3) = cbar(3) + dd*y_bar(i)
        end do
        !$omp end parallel do
        call self%model%vjp(sample, eta(:k), rb, q_bar, eta_bar(:k))
        eta_bar(k + 1:) = self%calibration%vjp(cbar)
    end subroutine
end module tiago_integrator_diagnostic
