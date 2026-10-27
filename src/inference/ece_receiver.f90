module tiago_ece_receiver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use tiago_calibration, only: calibration_prior_t
    use tiago_recorded_affine, only: recorded_affine
    use tiago_ece_radiometer_noise, only: ece_radiometer_noise
    implicit none
    private

    type, public :: ece_response_channel_t
        ! Each column of rays refers to one frequency; rows are antenna nodes.
        integer, allocatable :: rays(:, :)
        real(dp), allocatable :: frequency(:), weights(:), response(:), antenna(:)
        real(dp) :: duration = 0.0_dp
    end type

    type, public :: ece_receiver_t
        private
        integer :: n_rays = 0
        type(ece_response_channel_t), allocatable :: channels(:)
        real(dp), allocatable :: bandwidth(:), norm(:)
        type(calibration_prior_t) :: calibration
    contains
        procedure :: init => receiver_init
        procedure :: n_nuisance => receiver_nuisance
        procedure :: effective_bandwidth
        procedure :: evaluate => receiver_evaluate
        procedure :: jvp => receiver_jvp
        procedure :: vjp => receiver_vjp
    end type
contains
    subroutine receiver_init(self, channels, n_rays, prior, status)
        class(ece_receiver_t), intent(out) :: self
        type(ece_response_channel_t), intent(in) :: channels(:)
        integer, intent(in) :: n_rays
        type(calibration_prior_t), intent(in) :: prior
        integer, intent(out) :: status
        real(dp) :: area, square_area
        integer :: c, n, nf, na
        status = 1
        n = size(channels)
        if (n == 0 .or. n_rays < 1) return
        if (prior%n_params() /= 3*n) return
        if (.not. allocated(prior%factor)) return
        do c = 1, n
            if (prior%mean(3*c - 2) <= 0.0_dp .or. prior%mean(3*c) < 0.0_dp) return
            if (.not. allocated(channels(c)%rays)) return
            if (.not. allocated(channels(c)%frequency)) return
            if (.not. allocated(channels(c)%weights)) return
            if (.not. allocated(channels(c)%response)) return
            if (.not. allocated(channels(c)%antenna)) return
            nf = size(channels(c)%frequency)
            na = size(channels(c)%antenna)
            if (nf == 0 .or. na == 0) return
            if (any(shape(channels(c)%rays) /= [na, nf])) return
            if (any(channels(c)%rays < 1) .or. any(channels(c)%rays > n_rays)) return
            if (size(channels(c)%weights) /= nf .or. &
                size(channels(c)%response) /= nf) return
            if (.not. all(ieee_is_finite(channels(c)%frequency))) return
            if (.not. all(ieee_is_finite(channels(c)%weights))) return
            if (.not. all(ieee_is_finite(channels(c)%response))) return
            if (.not. all(ieee_is_finite(channels(c)%antenna))) return
            if (.not. ieee_is_finite(channels(c)%duration)) return
            if (any(channels(c)%frequency <= 0.0_dp)) return
            if (any(channels(c)%weights < 0.0_dp) .or. &
                any(channels(c)%response < 0.0_dp)) return
            if (any(channels(c)%antenna < 0.0_dp) .or. &
                channels(c)%duration <= 0.0_dp) return
            area = dot_product(channels(c)%weights, channels(c)%response)
            square_area = dot_product(channels(c)%weights, channels(c)%response**2)
            if (.not. all(ieee_is_finite([area, square_area, &
                sum(channels(c)%antenna)]))) return
            if (area <= 0.0_dp .or. square_area <= 0.0_dp) return
            if (sum(channels(c)%antenna) <= 0.0_dp) return
            if (.not. ieee_is_finite(area*(area/square_area))) return
            if (area*(area/square_area) <= 0.0_dp) return
        end do
        self%channels = channels
        self%calibration = prior
        self%n_rays = n_rays
        allocate (self%bandwidth(n), self%norm(n))
        do c = 1, n
            self%channels(c)%antenna = channels(c)%antenna/sum(channels(c)%antenna)
            area = dot_product(channels(c)%weights, channels(c)%response)
            square_area = dot_product(channels(c)%weights, channels(c)%response**2)
            self%norm(c) = area
            self%bandwidth(c) = area*(area/square_area)
        end do
        status = 0
    end subroutine

    integer function receiver_nuisance(self)
        class(ece_receiver_t), intent(in) :: self
        receiver_nuisance = self%calibration%n_params()
    end function

    function effective_bandwidth(self) result(bandwidth)
        class(ece_receiver_t), intent(in) :: self
        real(dp), allocatable :: bandwidth(:)
        if (.not. allocated(self%bandwidth)) error stop 'ECE receiver: uninitialized'
        bandwidth = self%bandwidth
    end function

    subroutine receiver_evaluate(self, temperature, eta, voltage, variance, status)
        class(ece_receiver_t), intent(in) :: self
        real(dp), intent(in) :: temperature(:), eta(:)
        real(dp), intent(out) :: voltage(:), variance(:)
        integer, intent(out) :: status
        real(dp) :: physical(size(eta)), mean_te, dr, dg, db, row(5), noise, tf
        integer :: c, f
        status = 1
        voltage = ieee_value(0.0_dp, ieee_quiet_nan)
        variance = ieee_value(0.0_dp, ieee_quiet_nan)
        if (.not. allocated(self%channels)) return
        if (size(voltage) /= size(self%channels) .or. &
            size(variance) /= size(voltage)) return
        if (.not. valid_inputs(self, temperature, eta, physical)) return
        !$omp parallel do schedule(static) &
        !$omp private(mean_te, dr, dg, db, row, noise, tf, f)
        do c = 1, size(self%channels)
            mean_te = channel_temperature(self, c, temperature)
            call recorded_affine(mean_te, physical(3*c - 2), physical(3*c - 1), &
                                 voltage(c), dr, dg, db)
            variance(c) = 0.0_dp
            do f = 1, size(self%channels(c)%frequency)
                tf = frequency_temperature(self, c, f, temperature)
                call ece_radiometer_noise(tf, physical(3*c), physical(3*c - 2), &
                                          1.0_dp, self%channels(c)%duration, noise, row)
                variance(c) = variance(c) + noise*noise_weight(self, c, f)
            end do
        end do
        !$omp end parallel do
        if (.not. all(ieee_is_finite(voltage))) return
        if (.not. all(ieee_is_finite(variance))) return
        status = 0
    end subroutine

    subroutine receiver_jvp(self, temperature, eta, te_dot, eta_dot, voltage_dot, &
                            variance_dot)
        class(ece_receiver_t), intent(in) :: self
        real(dp), intent(in) :: temperature(:), eta(:), te_dot(:), eta_dot(:)
        real(dp), intent(out) :: voltage_dot(:), variance_dot(:)
        real(dp) :: physical(size(eta)), dc(size(eta)), mean_te, dt, value, dr, dg, &
                    db, row(5), tf, df
        integer :: c, f
        if (.not. valid_inputs(self, temperature, eta, physical)) &
            error stop 'ECE receiver: invalid inputs'
        if (size(te_dot) /= self%n_rays .or. size(eta_dot) /= size(eta)) &
            error stop 'ECE receiver: tangent input shape'
        if (size(voltage_dot) /= size(self%channels) .or. &
            size(variance_dot) /= size(voltage_dot)) &
            error stop 'ECE receiver: tangent output shape'
        dc = self%calibration%jvp(eta_dot)
        !$omp parallel do schedule(static) &
        !$omp private(mean_te, dt, value, dr, dg, db, row, tf, df, f)
        do c = 1, size(self%channels)
            mean_te = channel_temperature(self, c, temperature)
            dt = channel_temperature(self, c, te_dot)
            call recorded_affine(mean_te, physical(3*c - 2), physical(3*c - 1), &
                value, dr, dg, db)
            voltage_dot(c) = dr*dt + dg*dc(3*c - 2) + db*dc(3*c - 1)
            variance_dot(c) = 0.0_dp
            do f = 1, size(self%channels(c)%frequency)
                tf = frequency_temperature(self, c, f, temperature)
                df = frequency_temperature(self, c, f, te_dot)
                call ece_radiometer_noise(tf, physical(3*c), physical(3*c - 2), &
                                          1.0_dp, self%channels(c)%duration, value, row)
                variance_dot(c) = variance_dot(c) + noise_weight(self, c, f)* &
                                  (row(1)*df + row(2)*dc(3*c) + row(3)*dc(3*c - 2))
            end do
        end do
        !$omp end parallel do
    end subroutine

    subroutine receiver_vjp(self, temperature, eta, voltage_bar, variance_bar, &
                            te_bar, eta_bar)
        class(ece_receiver_t), intent(in) :: self
        real(dp), intent(in) :: temperature(:), eta(:), voltage_bar(:), variance_bar(:)
        real(dp), intent(out) :: te_bar(:), eta_bar(:)
        real(dp) :: physical(size(eta)), cb(size(eta)), mean_te, value, dr, dg, db, &
                    row(5), seed, tf, nb
        integer :: c, a, f, k
        if (.not. valid_inputs(self, temperature, eta, physical)) &
            error stop 'ECE receiver: invalid inputs'
        if (size(te_bar) /= self%n_rays .or. size(eta_bar) /= size(eta)) &
            error stop 'ECE receiver: adjoint output shape'
        if (size(voltage_bar) /= size(self%channels) .or. &
            size(variance_bar) /= size(voltage_bar)) &
            error stop 'ECE receiver: adjoint input shape'
        te_bar = 0.0_dp
        cb = 0.0_dp
        !$omp parallel do schedule(static) &
        !$omp private(mean_te, value, dr, dg, db, row, seed, tf, nb, a, f, k) &
        !$omp reduction(+:te_bar)
        do c = 1, size(self%channels)
            mean_te = channel_temperature(self, c, temperature)
            call recorded_affine(mean_te, physical(3*c - 2), physical(3*c - 1), &
                value, dr, dg, db)
            seed = dr*voltage_bar(c)
            cb(3*c - 2) = dg*voltage_bar(c)
            cb(3*c - 1) = db*voltage_bar(c)
            do f = 1, size(self%channels(c)%frequency)
                tf = frequency_temperature(self, c, f, temperature)
                call ece_radiometer_noise(tf, physical(3*c), physical(3*c - 2), &
                                          1.0_dp, self%channels(c)%duration, value, row)
                nb = noise_weight(self, c, f)*variance_bar(c)
                cb(3*c - 2) = cb(3*c - 2) + row(3)*nb
                cb(3*c) = cb(3*c) + row(2)*nb
                do a = 1, size(self%channels(c)%antenna)
                    k = self%channels(c)%rays(a, f)
                    te_bar(k) = te_bar(k) + seed*node_weight(self, c, a, f) + &
                                nb*row(1)*self%channels(c)%antenna(a)
                end do
            end do
        end do
        !$omp end parallel do
        eta_bar = self%calibration%vjp(cb)
    end subroutine

    logical function valid_inputs(self, temperature, eta, physical)
        class(ece_receiver_t), intent(in) :: self
        real(dp), intent(in) :: temperature(:), eta(:)
        real(dp), intent(out) :: physical(:)
        integer :: j
        valid_inputs = .false.
        if (.not. allocated(self%channels)) return
        if (size(temperature) /= self%n_rays .or. size(eta) /= self%n_nuisance()) return
        if (.not. all(ieee_is_finite(temperature))) return
        if (.not. all(ieee_is_finite(eta))) return
        if (any(temperature < 0.0_dp)) return
        physical = self%calibration%physical(eta)
        if (.not. all(ieee_is_finite(physical))) return
        do j = 1, size(self%channels)
            if (physical(3*j - 2) <= 0.0_dp .or. physical(3*j) < 0.0_dp) return
        end do
        valid_inputs = .true.
    end function

    real(dp) function node_weight(self, c, a, f) result(w)
        class(ece_receiver_t), intent(in) :: self
        integer, intent(in) :: c, a, f
        w = self%channels(c)%antenna(a)*self%channels(c)%weights(f)* &
            self%channels(c)%response(f)/self%norm(c)
    end function

    real(dp) function noise_weight(self, c, f) result(w)
        class(ece_receiver_t), intent(in) :: self
        integer, intent(in) :: c, f
        w = self%channels(c)%weights(f)*(self%channels(c)%response(f)/self%norm(c))**2
    end function

    real(dp) function frequency_temperature(self, c, f, temperature) result(t)
        class(ece_receiver_t), intent(in) :: self
        integer, intent(in) :: c, f
        real(dp), intent(in) :: temperature(:)
        integer :: a
        t = 0.0_dp
        do a = 1, size(self%channels(c)%antenna)
            t = t + self%channels(c)%antenna(a)*temperature(self%channels(c)%rays(a, f))
        end do
    end function

    real(dp) function channel_temperature(self, c, temperature) result(t)
        class(ece_receiver_t), intent(in) :: self
        integer, intent(in) :: c
        real(dp), intent(in) :: temperature(:)
        integer :: a, f
        t = 0.0_dp
        do f = 1, size(self%channels(c)%frequency)
            do a = 1, size(self%channels(c)%antenna)
                t = t + node_weight(self, c, a, &
                    f)*temperature(self%channels(c)%rays(a, f))
            end do
        end do
    end function
end module
