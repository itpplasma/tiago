module tiago_digitizer_diagnostic
    !! Digitiser stage: maps a recorded volt-level diagnostic to ADC counts.
    !! Per-channel volts-per-count and offset are calibration latents with
    !! priors, ordered [vres_1..vres_n, voff_1..voff_n]. Predictions are
    !! expected counts as reals (no rounding) so derivatives exist; use
    !! tiago_adc for integer quantisation and saturation flags.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t
    use tiago_recorded_diagnostic, only: recorded_diagnostic_t
    use tiago_calibration, only: calibration_prior_t
    use tiago_digitizer_signal, only: digitizer_signal
    implicit none
    private
    public :: digitizer_diagnostic_t

    type, extends(recorded_diagnostic_t), public :: digitizer_diagnostic_t
        class(profile_diagnostic_t), allocatable :: model
        integer :: bits = 16
    contains
        procedure :: init => digitizer_init
        procedure :: n_signals => digitizer_count
        procedure :: n_nuisance => digitizer_nuisance
        procedure :: label => digitizer_label
        procedure :: sample_points => digitizer_points
        procedure :: predict => digitizer_predict
        procedure :: jvp => digitizer_jvp
        procedure :: vjp => digitizer_vjp
    end type
contains
    subroutine digitizer_init(self, model, prior, bits, status)
        class(digitizer_diagnostic_t), intent(inout) :: self
        class(profile_diagnostic_t), intent(in) :: model
        type(calibration_prior_t), intent(in) :: prior
        integer, intent(in) :: bits
        integer, intent(out) :: status
        integer :: n
        status = 1
        n = model%n_signals()
        if (prior%n_params() /= 2*n) return
        if (bits < 2 .or. bits > 31) return
        ! vres means are the first n prior means and must be positive.
        if (.not. all(ieee_is_finite(prior%mean))) return
        if (any(prior%mean(:n) <= 0.0_dp)) return
        if (.not. allocated(prior%factor)) return
        if (allocated(self%model)) deallocate (self%model)
        allocate (self%model, source=model)
        self%calibration = prior
        self%bits = bits
        status = 0
    end subroutine

    integer function digitizer_count(self)
        class(digitizer_diagnostic_t), intent(in) :: self
        digitizer_count = self%model%n_signals()
    end function

    integer function digitizer_nuisance(self)
        class(digitizer_diagnostic_t), intent(in) :: self
        digitizer_nuisance = self%model%n_nuisance() + self%n_calibration()
    end function

    function digitizer_label(self, i) result(text)
        class(digitizer_diagnostic_t), intent(in) :: self
        integer, intent(in) :: i
        character(:), allocatable :: text
        text = 'counts:'//self%model%label(i)
    end function

    subroutine digitizer_points(self, eta, x, dx_deta)
        class(digitizer_diagnostic_t), intent(in) :: self
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

    subroutine digitizer_predict(self, sample, eta, y)
        class(digitizer_diagnostic_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        real(dp), allocatable :: volts(:), c(:)
        real(dp) :: dvolt, dvres, dvoff
        integer :: i, k, n
        allocate (volts(size(y)), c(2*size(y)))
        k = self%model%n_nuisance()
        n = self%n_signals()
        c = self%calibration%physical(eta(k + 1:))
        call self%model%predict(sample, eta(:k), volts)
        !$omp parallel do default(shared) private(i, dvolt, dvres, dvoff) &
        !$omp if(size(y)>1024)
        do i = 1, size(y)
            call digitizer_signal(volts(i), c(i), c(n + i), y(i), dvolt, dvres, dvoff)
        end do
        !$omp end parallel do
    end subroutine

    subroutine digitizer_jvp(self, sample, eta, q_dot, eta_dot, y_dot)
        class(digitizer_diagnostic_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        real(dp), allocatable :: volts(:), dv(:), c(:), dc(:)
        real(dp) :: value, dvolt, dvres, dvoff
        integer :: i, k, n
        allocate (volts(size(y_dot)), dv(size(y_dot)), &
            c(2*size(y_dot)), dc(2*size(y_dot)))
        k = self%model%n_nuisance()
        n = self%n_signals()
        c = self%calibration%physical(eta(k + 1:))
        dc = self%calibration%jvp(eta_dot(k + 1:))
        call self%model%predict(sample, eta(:k), volts)
        call self%model%jvp(sample, eta(:k), q_dot, eta_dot(:k), dv)
        !$omp parallel do default(shared) private(i, value, dvolt, dvres, dvoff) &
        !$omp if(size(y_dot)>1024)
        do i = 1, size(y_dot)
            call digitizer_signal(volts(i), c(i), c(n + i), value, dvolt, dvres, dvoff)
            y_dot(i) = dvolt*dv(i) + dvres*dc(i) + dvoff*dc(n + i)
        end do
        !$omp end parallel do
    end subroutine

    subroutine digitizer_vjp(self, sample, eta, y_bar, q_bar, eta_bar)
        class(digitizer_diagnostic_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        real(dp), allocatable :: volts(:), rb(:), c(:), cbar(:)
        real(dp) :: value, dvolt, dvres, dvoff
        integer :: i, k, n
        allocate (volts(size(y_bar)), rb(size(y_bar)), c(2*size(y_bar)), &
            cbar(2*size(y_bar)))
        k = self%model%n_nuisance()
        n = self%n_signals()
        c = self%calibration%physical(eta(k + 1:))
        call self%model%predict(sample, eta(:k), volts)
        cbar = 0.0_dp
        !$omp parallel do default(shared) private(i, value, dvolt, dvres, dvoff) &
        !$omp reduction(+:cbar) if(size(y_bar)>1024)
        do i = 1, size(y_bar)
            call digitizer_signal(volts(i), c(i), c(n + i), value, dvolt, dvres, dvoff)
            rb(i) = dvolt*y_bar(i)
            cbar(i) = cbar(i) + dvres*y_bar(i)
            cbar(n + i) = cbar(n + i) + dvoff*y_bar(i)
        end do
        !$omp end parallel do
        call self%model%vjp(sample, eta(:k), rb, q_bar, eta_bar(:k))
        eta_bar(k + 1:) = self%calibration%vjp(cbar)
    end subroutine
end module tiago_digitizer_diagnostic
