module tiago_ece_radiometer
    !! Recorded ECE on prepared rays; see docs/ece-recorded.md.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_profile_diagnostic, only: plasma_sample_t, Q_B, Q_NE, Q_TE
    use tiago_recorded_diagnostic, only: recorded_diagnostic_t
    use tiago_calibration, only: calibration_prior_t
    use tiago_ece_ray, only: ECE_X_MODE
    use tiago_ece_opacity, only: ECE_X1_TENUOUS
    use tiago_ece_opacity_transfer, only: thermal_ece_transfer, &
        ece_opacity_transfer_t
    use tiago_ece_receiver, only: ece_receiver_t, ece_response_channel_t
    implicit none
    private

    type, public :: ece_channel_t
        character(len=:), allocatable :: label
        real(dp) :: frequency = 0.0_dp
        integer :: harmonic = 2
        real(dp) :: origin(3) = 0.0_dp
        real(dp) :: direction(3) = [-1.0_dp, 0.0_dp, 0.0_dp]
        real(dp) :: length = 1.0_dp !! Hamiltonian duration, metres
        integer :: nscan = 1000
        integer :: mode = ECE_X_MODE
    end type

    type, public, extends(recorded_diagnostic_t) :: ece_t
        type(ece_channel_t), allocatable :: channels(:)
        integer, allocatable :: first(:)
        type(ece_receiver_t) :: receiver
        type(calibration_prior_t) :: wall_prior
        real(dp), allocatable :: points(:, :), cell_length(:), parallel_cos2(:)
    contains
        procedure :: init => ece_init
        procedure :: n_signals => ece_n_signals
        procedure :: n_nuisance => ece_n_nuisance
        procedure :: n_calibration => ece_n_nuisance
        procedure :: label => ece_label
        procedure :: sample_points => ece_points
        procedure :: predict => ece_predict
        procedure :: jvp => ece_jvp
        procedure :: vjp => ece_vjp
        procedure :: record => ece_record
        procedure :: record_jvp
        procedure :: record_vjp
        procedure :: brightness
    end type
contains
    subroutine ece_init(self, channels, response, prior, wall_prior, status)
        class(ece_t), intent(inout) :: self
        type(ece_channel_t), intent(in) :: channels(:)
        type(ece_response_channel_t), intent(in) :: response(:)
        type(calibration_prior_t), intent(in) :: prior, wall_prior
        integer, intent(out) :: status
        integer :: c, i, f
        status = 1
        if (size(channels) == 0 .or. size(response) == 0) return
        if (wall_prior%n_params() /= 2*size(channels)) return
        if (.not. allocated(wall_prior%factor)) return
        do c = 1, size(channels)
            if (.not. allocated(channels(c)%label)) return
            if (channels(c)%nscan < 1) return
            if (.not. all(ieee_is_finite([channels(c)%frequency, &
                channels(c)%origin, channels(c)%direction, &
                channels(c)%length]))) return
            if (channels(c)%frequency <= 0 .or. channels(c)%length <= 0) return
            if (norm2(channels(c)%direction) <= 0) return
        end do
        call self%receiver%init(response, size(channels), prior, status)
        if (status /= 0) return
        status = 1
        do c = 1, size(response)
            do f = 1, size(response(c)%frequency)
                do i = 1, size(response(c)%rays, 1)
                    if (abs(channels(response(c)%rays(i, f))%frequency - &
                        response(c)%frequency(f)) > 0) return
                end do
            end do
        end do
        self%channels = channels
        self%calibration = prior
        self%wall_prior = wall_prior
        if (allocated(self%first)) deallocate (self%first)
        allocate (self%first(size(channels) + 1))
        self%first(1) = 1
        do c = 1, size(channels)
            self%first(c + 1) = self%first(c) + channels(c)%nscan
        end do
        if (allocated(self%points)) deallocate (self%points, self%cell_length, &
            self%parallel_cos2)
        status = 0
    end subroutine

    integer function ece_n_signals(self)
        class(ece_t), intent(in) :: self
        ece_n_signals = self%calibration%n_params()/3
    end function

    integer function ece_n_nuisance(self)
        class(ece_t), intent(in) :: self
        ece_n_nuisance = self%calibration%n_params() + self%wall_prior%n_params()
    end function

    function ece_label(self, i) result(text)
        class(ece_t), intent(in) :: self
        integer, intent(in) :: i
        character(:), allocatable :: text
        character(24) :: number
        if (i < 1 .or. i > self%n_signals()) error stop 'ECE: signal index'
        write (number, '(I0)') i
        text = 'ECE voltage '//trim(number)
    end function

    subroutine ece_points(self, eta, x, dx_deta)
        class(ece_t), intent(in) :: self
        real(dp), intent(in) :: eta(:)
        real(dp), allocatable, intent(out) :: x(:, :), dx_deta(:, :, :)
        if (.not. allocated(self%points)) error stop &
            'ECE: prepare refracted rays first'
        if (size(eta) /= self%n_nuisance()) error stop 'ECE: nuisance shape'
        x = self%points
        allocate (dx_deta(3, size(x, 2), size(eta)))
        dx_deta = 0
    end subroutine

    subroutine brightness(self, sample, eta, traces, temperature, status)
        class(ece_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        type(ece_opacity_transfer_t), intent(out) :: traces(:)
        real(dp), intent(out) :: temperature(:)
        integer, intent(out) :: status
        real(dp) :: walls(size(self%wall_prior%mean))
        integer :: c, lo, hi, code(size(self%channels)), nc
        status = 1
        if (.not. allocated(self%points)) return
        if (size(eta) /= self%n_nuisance()) return
        if (size(sample%q, 2) /= size(self%cell_length)) return
        nc = self%calibration%n_params()
        walls = self%wall_prior%physical(eta(nc + 1:))
        ! The bound makes parallel Doppler broadening small relative to the
        ! relativistic width. It is checked wherever opacity is evaluated.
        if (any(self%parallel_cos2 > 1.0e-8_dp*sample%q(Q_TE, :)/ &
            510998.95069_dp)) then
            status = 3
            return
        end if
        !$omp parallel do default(shared) private(c, lo, hi) schedule(dynamic, 1)
        do c = 1, size(self%channels)
            lo = self%first(c)
            hi = self%first(c + 1) - 1
            call thermal_ece_transfer(sample%q(Q_TE, lo:hi), &
                sample%q(Q_NE, lo:hi), sample%q(Q_B, lo:hi), &
                self%channels(c)%frequency, &
                self%channels(c)%mode, &
                self%channels(c)%harmonic, &
                self%cell_length(lo:hi), &
                walls(2*c - 1), walls(2*c), traces(c), code(c), &
                ECE_X1_TENUOUS)
            temperature(c) = traces(c)%temperature
        end do
        !$omp end parallel do
        status = maxval(code)
    end subroutine

    subroutine ece_record(self, sample, eta, y, variance, status)
        class(ece_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:), variance(:)
        integer, intent(out) :: status
        type(ece_opacity_transfer_t) :: traces(size(self%channels))
        real(dp) :: temperature(size(self%channels))
        call self%brightness(sample, eta, traces, temperature, status)
        if (status /= 0) return
        call self%receiver%evaluate(temperature, &
            eta(:self%calibration%n_params()), &
            y, variance, status)
    end subroutine

    subroutine ece_predict(self, sample, eta, y)
        class(ece_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: variance(size(y))
        integer :: status
        call self%record(sample, eta, y, variance, status)
        if (status == 3) error stop 'ECE: ray exceeds perpendicular-opacity angle bound'
        if (status /= 0) error stop 'ECE: invalid opacity regime or receiver inputs'
    end subroutine

    subroutine record_jvp(self, sample, eta, q_dot, eta_dot, length_dot, &
        y_dot, variance_dot, status)
        class(ece_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:), length_dot(:)
        real(dp), intent(out) :: y_dot(:), variance_dot(:)
        integer, intent(out) :: status
        type(ece_opacity_transfer_t) :: traces(size(self%channels))
        real(dp) :: temperature(size(self%channels)), dt(size(self%channels))
        real(dp) :: dw(size(self%wall_prior%mean))
        integer :: c, lo, hi, nc
        call self%brightness(sample, eta, traces, temperature, status)
        if (status /= 0) return
        nc = self%calibration%n_params()
        dw = self%wall_prior%jvp(eta_dot(nc + 1:))
        !$omp parallel do default(shared) private(c, lo, hi) schedule(static)
        do c = 1, size(self%channels)
            lo = self%first(c)
            hi = self%first(c + 1) - 1
            call traces(c)%jvp(q_dot(Q_TE, lo:hi), q_dot(Q_NE, lo:hi), &
                q_dot(Q_B, lo:hi), 0.0_dp, length_dot(lo:hi), &
                dw(2*c - 1), dw(2*c), dt(c))
        end do
        !$omp end parallel do
        call self%receiver%jvp(temperature, eta(:nc), dt, eta_dot(:nc), &
            y_dot, variance_dot)
    end subroutine

    subroutine ece_jvp(self, sample, eta, q_dot, eta_dot, y_dot)
        class(ece_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        real(dp) :: dv(size(y_dot))
        integer :: status
        call self%record_jvp(sample, eta, q_dot, eta_dot, &
            0*self%cell_length, y_dot, dv, status)
        if (status /= 0) error stop 'ECE: invalid tangent evaluation'
    end subroutine

    subroutine record_vjp(self, sample, eta, y_bar, variance_bar, q_bar, &
        eta_bar, length_bar, status)
        class(ece_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:), variance_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:), length_bar(:)
        integer, intent(out) :: status
        type(ece_opacity_transfer_t) :: traces(size(self%channels))
        real(dp) :: temperature(size(self%channels)), tb(size(self%channels))
        real(dp) :: wb(size(self%wall_prior%mean)), fb
        integer :: c, lo, hi, nc
        q_bar = 0
        eta_bar = 0
        length_bar = 0
        call self%brightness(sample, eta, traces, temperature, status)
        if (status /= 0) return
        nc = self%calibration%n_params()
        call self%receiver%vjp(temperature, eta(:nc), y_bar, variance_bar, &
            tb, eta_bar(:nc))
        !$omp parallel do default(shared) private(c, lo, hi, fb) schedule(static)
        do c = 1, size(self%channels)
            lo = self%first(c)
            hi = self%first(c + 1) - 1
            call traces(c)%vjp(tb(c), q_bar(Q_TE, lo:hi), q_bar(Q_NE, lo:hi), &
                q_bar(Q_B, lo:hi), fb, length_bar(lo:hi), wb(2*c - 1), &
                wb(2*c))
        end do
        !$omp end parallel do
        eta_bar(nc + 1:) = self%wall_prior%vjp(wb)
    end subroutine

    subroutine ece_vjp(self, sample, eta, y_bar, q_bar, eta_bar)
        class(ece_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        real(dp) :: lb(size(self%cell_length))
        integer :: status
        call self%record_vjp(sample, eta, y_bar, 0*y_bar, q_bar, eta_bar, lb, &
            status)
        if (status /= 0) error stop 'ECE: invalid adjoint evaluation'
    end subroutine
end module
