module tiago_harmonic_lockin
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_harmonic_mixer, only: harmonic_mixer
    implicit none
    private
    real(dp), parameter :: twopi = 2.0_dp*acos(-1.0_dp)

    type, public :: harmonic_lockin_config_t
        integer(int64) :: cadence_ns = 0, epoch_ns = 0, delay_ns = 0
        integer(int64) :: valid_from_ns = 0, valid_to_ns = -1
        integer :: samples_per_period = 0, reference_polarity = 0
        real(dp) :: frequency_hz = 0, modulation_depth = 0
        real(dp) :: detector_gain = 0, detector_offset = 0
        real(dp) :: reference_gain = 0, reference_offset = 0
        real(dp) :: reference_phase = 0, reference_tolerance = 0
        real(dp) :: iq_gain(2) = 0, iq_offset(2) = 0
        real(dp) :: minimum_normalizer = 0, alias_tolerance = 0
        real(dp) :: minimum_radius = 0
    end type

    type, public :: harmonic_lockin_output_t
        integer(int64), allocatable :: time_ns(:)
        real(dp), allocatable :: iq(:, :)
        logical, allocatable :: valid(:)
    end type

    type, public :: harmonic_lockin_t
        private
        type(harmonic_lockin_config_t) :: config
        real(dp), allocatable :: weights(:, :), reference(:)
        real(dp) :: normalizer(2) = 0, accumulator(2) = 0
        real(dp) :: weight_sum(2) = 0, detector_sum = 0
        integer :: position = 0
        integer(int64) :: last_ns = 0
        logical :: initialized = .false., started = .false.
    contains
        procedure :: initialize
        procedure :: consume
        procedure :: buffered_samples
    end type
contains
    subroutine initialize(self, config, status)
        class(harmonic_lockin_t), intent(inout) :: self
        type(harmonic_lockin_config_t), intent(in) :: config
        integer, intent(out) :: status
        real(dp) :: theta, phase, value(3), jac(3, 4), discrete(2), dc(2)
        real(dp) :: leakage(2)
        real(dp), allocatable :: weights(:, :), reference(:)
        integer :: k, n
        status = 1
        n = config%samples_per_period
        if (n < 8 .or. config%cadence_ns <= 0) return
        if (config%cadence_ns > huge(0_int64)/int(n, int64)) return
        if (config%reference_polarity /= 1 .and. config%reference_polarity /= -1) return
        if (config%valid_to_ns < config%valid_from_ns) return
        if (config%epoch_ns < config%valid_from_ns) return
        if (config%epoch_ns > config%valid_to_ns) return
        if (.not. all(ieee_is_finite([config%frequency_hz, &
            config%modulation_depth, config%detector_gain, config%detector_offset, &
            config%reference_gain, config%reference_offset, config%reference_phase, &
            config%reference_tolerance, config%iq_gain, config%iq_offset, &
            config%minimum_normalizer, config%alias_tolerance, &
            config%minimum_radius]))) return
        if (config%frequency_hz <= 0 .or. config%modulation_depth <= 0) return
        if (config%detector_gain == 0 .or. config%reference_gain <= 0) return
        if (any(config%iq_gain == 0)) return
        if (config%reference_tolerance <= 0 .or. config%alias_tolerance <= 0) return
        if (config%minimum_normalizer <= 0 .or. config%minimum_radius < 0) return
        if (abs(config%frequency_hz*real(config%cadence_ns, dp)*1e-9_dp &
            *real(n, dp) - 1.0_dp) > 64*epsilon(1.0_dp)) return
        dc = [bessel_jn(2, config%modulation_depth), &
            -bessel_jn(1, config%modulation_depth)]
        if (any(abs(dc) <= config%minimum_normalizer)) return
        allocate (weights(2, n), reference(n))
        phase = config%reference_phase - twopi*real(modulo(config%delay_ns, &
            config%cadence_ns*int(n, int64)), dp) &
            /real(config%cadence_ns*int(n, int64), dp)
        discrete = 0
        leakage = 0
        do k = 1, n
            theta = twopi*real(k - 1, dp)/real(n, dp) + phase
            call harmonic_mixer(0.0_dp, 0.0_dp, 1.0_dp, theta, value, jac)
            weights(:, k) = value(2:3)
            reference(k) = real(config%reference_polarity, dp)*sin(&
                twopi*real(k - 1, dp)/real(n, dp) + config%reference_phase)
            discrete(1) = discrete(1) + cos(config%modulation_depth*sin(theta)) &
                *weights(1, k)/real(n, dp)
            discrete(2) = discrete(2) + sin(config%modulation_depth*sin(theta)) &
                *sin(theta)/real(n, dp)
            leakage(1) = leakage(1) &
                + sin(config%modulation_depth*sin(theta))*weights(1, k)/real(n, dp)
            leakage(2) = leakage(2) &
                + cos(config%modulation_depth*sin(theta))*sin(theta)/real(n, dp)
        end do
        ! Reject grids whose harmonic aliasing exceeds the declared tolerance.
        if (any(abs(discrete - [dc(1), -dc(2)]) > config%alias_tolerance)) return
        if (any(abs(leakage) > config%alias_tolerance)) return
        call move_alloc(weights, self%weights)
        call move_alloc(reference, self%reference)
        self%normalizer = real(n, dp)*dc
        self%config = config
        self%position = 0
        self%accumulator = 0
        self%weight_sum = 0
        self%detector_sum = 0
        self%started = .false.
        self%initialized = .true.
        status = 0
    end subroutine

    subroutine consume(self, time_ns, detector, reference, output, status)
        class(harmonic_lockin_t), intent(inout) :: self
        integer(int64), intent(in) :: time_ns(:)
        real(dp), intent(in) :: detector(:), reference(:)
        type(harmonic_lockin_output_t), intent(out) :: output
        integer, intent(out) :: status
        integer :: k, p, n, m, count
        integer(int64) :: previous, dt
        real(dp) :: sample, measured, value(3), jac(3, 4), theta
        status = 1
        if (.not. self%initialized) return
        m = size(time_ns)
        if (size(detector) /= m .or. size(reference) /= m) return
        if (.not. all(ieee_is_finite(detector))) return
        if (.not. all(ieee_is_finite(reference))) return
        n = self%config%samples_per_period
        if (m > huge(m) - n) return
        dt = self%config%cadence_ns
        previous = self%last_ns
        ! Preflight the complete chunk so errors never partially advance the state.
        do k = 1, m
            if (time_ns(k) < self%config%valid_from_ns) return
            if (time_ns(k) > self%config%valid_to_ns) return
            if (k == 1 .and. .not. self%started) then
                if (time_ns(k) /= self%config%epoch_ns) return
            else
                if (previous > huge(previous) - dt) return
                if (time_ns(k) /= previous + dt) return
            end if
            previous = time_ns(k)
            p = modulo(self%position + k - 1, n) + 1
            measured = (reference(k) - self%config%reference_offset) &
                /self%config%reference_gain
            if (.not. ieee_is_finite(measured)) return
            if (abs(measured - self%reference(p)) > &
                self%config%reference_tolerance) return
            sample = (detector(k) - self%config%detector_offset) &
                /self%config%detector_gain
            if (.not. ieee_is_finite(sample)) return
        end do
        count = (self%position + m)/n
        allocate (output%iq(2, count), output%time_ns(count), output%valid(count))
        count = 0
        do k = 1, m
            p = self%position + 1
            theta = 0 ! Only the affine row is used; oscillator weights are cached.
            call harmonic_mixer(detector(k), self%config%detector_offset, &
                self%config%detector_gain, theta, value, jac)
            sample = value(1)
            self%accumulator = self%accumulator + sample*self%weights(:, p)
            self%weight_sum = self%weight_sum + self%weights(:, p)
            self%detector_sum = self%detector_sum + sample
            self%position = p
            self%last_ns = time_ns(k)
            self%started = .true.
            if (p /= n) cycle
            count = count + 1
            output%iq(:, count) = ((self%accumulator &
                - self%detector_sum*self%weight_sum/real(n, dp)) &
                /self%normalizer - self%config%iq_offset)/self%config%iq_gain
            output%time_ns(count) = time_ns(k)
            output%valid(count) = all(ieee_is_finite(output%iq(:, count)))
            output%valid(count) = output%valid(count) .and. &
                norm2(output%iq(:, count)) > self%config%minimum_radius
            self%position = 0
            self%accumulator = 0
            self%weight_sum = 0
            self%detector_sum = 0
        end do
        status = 0
    end subroutine

    integer function buffered_samples(self) result(n)
        class(harmonic_lockin_t), intent(in) :: self
        n = self%position
    end function
end module
