module tiago_magnetic_pickup
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_fft, only: fft_c2c
    use tiago_calibration, only: calibration_prior_t
    use tiago_magnetic_pickup_leaf, only: magnetic_pickup_leaf
    implicit none
    private
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer(int64), parameter :: exact_ns = 9007199254740992_int64
    type, public :: magnetic_pickup_t
        private
        type(calibration_prior_t) :: calibration
        complex(dp), allocatable :: response(:)
        real(dp) :: dt = 0
        logical :: valid = .false.
    contains
        procedure :: init => pickup_init
        procedure :: predict => pickup_predict
        procedure :: jvp => pickup_jvp
        procedure :: vjp => pickup_vjp
    end type
contains
    subroutine pickup_init(self, utc_ns, reference_ns, prior, response, &
            max_frequency, flux_units, voltage_units, &
            response_units, status)
        class(magnetic_pickup_t), intent(inout) :: self
        integer(int64), intent(in) :: utc_ns(:), reference_ns
        type(calibration_prior_t), intent(in) :: prior
        complex(dp), intent(in) :: response(:)
        real(dp), intent(in) :: max_frequency
        character(*), intent(in) :: flux_units, voltage_units, response_units
        integer, intent(out) :: status
        integer(int64) :: step, delta
        integer :: n, j
        real(dp) :: time
        self%valid = .false.
        status = 1
        n = size(utc_ns)
        if (n < 3 .or. mod(n, 2) == 0) return
        if (size(response) /= n) return
        if (flux_units /= 'Wb' .or. voltage_units /= 'V') return
        if (response_units /= 'V/V') return
        ! Nonnegative UTC nanoseconds make subtraction overflow impossible.
        if (reference_ns < 0 .or. any(utc_ns < 0)) return
        if (utc_ns(2) <= utc_ns(1)) return
        step = utc_ns(2) - utc_ns(1)
        if (step > exact_ns) return
        do j = 2, n
            if (utc_ns(j) <= utc_ns(j - 1)) return
            if (utc_ns(j) - utc_ns(j - 1) /= step) return
        end do
        if (prior%n_params() /= 2) return
        if (.not. allocated(prior%factor)) return
        if (any(shape(prior%factor) /= [2, 2])) return
        if (.not. all(ieee_is_finite(prior%factor))) return
        if (.not. all(ieee_is_finite(prior%mean))) return
        if (.not. ieee_is_finite(prior%valid_from)) return
        if (.not. ieee_is_finite(prior%valid_to)) return
        if (prior%valid_to < prior%valid_from) return
        if (.not. allocated(prior%record%record_id)) return
        if (.not. allocated(prior%record%source_id)) return
        if (len_trim(prior%record%record_id) == 0) return
        if (len_trim(prior%record%source_id) == 0) return
        if (.not. allocated(prior%record%metadata)) return
        if (len_trim(prior%record%metadata) == 0) return
        do j = 1, n
            delta = utc_ns(j) - reference_ns
            if (abs(delta) > exact_ns) return
            time = real(delta, dp)*1.e-9_dp
            if (time < prior%valid_from .or. time > prior%valid_to) return
        end do
        self%dt = real(step, dp)*1.e-9_dp
        if (.not. ieee_is_finite(max_frequency)) return
        if (max_frequency < real(n/2, dp)/(n*self%dt)) return
        if (.not. all(ieee_is_finite(real(response, dp)))) return
        if (.not. all(ieee_is_finite(aimag(response)))) return
        ! Exact conjugate pairs ensure a real operator without hidden projection.
        if (aimag(response(1)) /= 0) return
        do j = 2, n/2 + 1
            if (response(n - j + 2) /= conjg(response(j))) return
        end do
        self%calibration = prior
        self%response = response
        self%valid = .true.
        status = 0
    end subroutine

    subroutine transfer(self, z, t, dg, dd, status)
        class(magnetic_pickup_t), intent(in) :: self
        real(dp), intent(in) :: z(:)
        complex(dp), allocatable, intent(out) :: t(:), dg(:), dd(:)
        integer, intent(out) :: status
        real(dp) :: c(2), omega, tr, ti, gr, gi, dr, di
        integer :: n, j, k
        status = 1
        if (.not. self%valid) return
        if (size(z) /= 2) return
        if (.not. all(ieee_is_finite(z))) return
        c = self%calibration%physical(z)
        if (.not. all(ieee_is_finite(c))) return
        n = size(self%response)
        ! Bound the delay to half a period, avoiding ambiguous wrapped clocks.
        if (abs(c(2)) > n*self%dt/2) return
        allocate (t(n), dg(n), dd(n))
        !$omp parallel do default(shared) private(j,k,omega,tr,ti,gr,gi,dr,di) &
        !$omp if(n>1024)
        do j = 1, n
            k = j - 1
            if (k > n/2) k = k - n
            omega = 2*pi*real(k, dp)/(n*self%dt)
            call magnetic_pickup_leaf(omega, c(1), c(2), &
                real(self%response(j), dp), &
                aimag(self%response(j)), &
                tr, ti, gr, gi, dr, di)
            t(j) = cmplx(tr, ti, dp)
            dg(j) = cmplx(gr, gi, dp)
            dd(j) = cmplx(dr, di, dp)
        end do
        !$omp end parallel do
        if (.not. all(ieee_is_finite(real(t, dp)))) return
        if (.not. all(ieee_is_finite(aimag(t)))) return
        if (.not. all(ieee_is_finite(real(dg, dp)))) return
        if (.not. all(ieee_is_finite(aimag(dg)))) return
        if (.not. all(ieee_is_finite(real(dd, dp)))) return
        if (.not. all(ieee_is_finite(aimag(dd)))) return
        status = 0
    end subroutine

    logical function input_valid(self, flux, y)
        class(magnetic_pickup_t), intent(in) :: self
        real(dp), intent(in) :: flux(:), y(:)
        input_valid = .false.
        if (.not. self%valid) return
        if (size(flux) /= size(self%response)) return
        if (size(y) /= size(flux)) return
        if (.not. all(ieee_is_finite(flux))) return
        input_valid = .true.
    end function

    subroutine pickup_predict(self, flux, z, voltage, status)
        class(magnetic_pickup_t), intent(in) :: self
        real(dp), intent(in) :: flux(:), z(:)
        real(dp), intent(out) :: voltage(:)
        integer, intent(out) :: status
        complex(dp), allocatable :: x(:), t(:), dg(:), dd(:)
        voltage = 0
        status = 1
        if (.not. input_valid(self, flux, voltage)) return
        call transfer(self, z, t, dg, dd, status)
        if (status /= 0) return
        x = cmplx(flux, 0.0_dp, dp)
        call fft_c2c(x, -1)
        x = t*x
        call fft_c2c(x, 1)
        call real_output(x, voltage, status)
    end subroutine

    subroutine pickup_jvp(self, flux, z, flux_dot, z_dot, voltage_dot, status)
        class(magnetic_pickup_t), intent(in) :: self
        real(dp), intent(in) :: flux(:), z(:), flux_dot(:), z_dot(:)
        real(dp), intent(out) :: voltage_dot(:)
        integer, intent(out) :: status
        complex(dp), allocatable :: x(:), dx(:), t(:), dg(:), dd(:)
        real(dp) :: dc(2)
        voltage_dot = 0
        status = 1
        if (.not. input_valid(self, flux, voltage_dot)) return
        if (size(flux_dot) /= size(flux) .or. size(z_dot) /= 2) return
        if (.not. all(ieee_is_finite(flux_dot))) return
        if (.not. all(ieee_is_finite(z_dot))) return
        call transfer(self, z, t, dg, dd, status)
        if (status /= 0) return
        dc = self%calibration%jvp(z_dot)
        x = cmplx(flux, 0.0_dp, dp)
        dx = cmplx(flux_dot, 0.0_dp, dp)
        call fft_c2c(x, -1)
        call fft_c2c(dx, -1)
        dx = t*dx + (dg*dc(1) + dd*dc(2))*x
        call fft_c2c(dx, 1)
        call real_output(dx, voltage_dot, status)
    end subroutine

    subroutine pickup_vjp(self, flux, z, voltage_bar, flux_bar, z_bar, status)
        class(magnetic_pickup_t), intent(in) :: self
        real(dp), intent(in) :: flux(:), z(:), voltage_bar(:)
        real(dp), intent(out) :: flux_bar(:), z_bar(:)
        integer, intent(out) :: status
        complex(dp), allocatable :: x(:), b(:), t(:), dg(:), dd(:)
        real(dp) :: cb(2)
        flux_bar = 0
        z_bar = 0
        status = 1
        if (.not. input_valid(self, flux, flux_bar)) return
        if (size(voltage_bar) /= size(flux) .or. size(z_bar) /= 2) return
        if (.not. all(ieee_is_finite(voltage_bar))) return
        call transfer(self, z, t, dg, dd, status)
        if (status /= 0) return
        x = cmplx(flux, 0.0_dp, dp)
        b = cmplx(voltage_bar, 0.0_dp, dp)
        call fft_c2c(x, -1)
        call fft_c2c(b, -1)
        cb = [real(sum(conjg(b)*dg*x), dp), &
            real(sum(conjg(b)*dd*x), dp)]/size(x)
        b = conjg(t)*b
        call fft_c2c(b, 1)
        call real_output(b, flux_bar, status)
        if (status /= 0) return
        z_bar = self%calibration%vjp(cb)
        if (.not. all(ieee_is_finite(z_bar))) status = 1
    end subroutine

    subroutine real_output(x, y, status)
        complex(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer, intent(out) :: status
        status = 1
        y = 0
        if (.not. all(ieee_is_finite(real(x, dp)))) return
        if (.not. all(ieee_is_finite(aimag(x)))) return
        if (maxval(abs(aimag(x))) > 1.e-10_dp*max(1.0_dp, maxval(abs(x)))) return
        y = real(x, dp)/size(x)
        status = 0
    end subroutine
end module
