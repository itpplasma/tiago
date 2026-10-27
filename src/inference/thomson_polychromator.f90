module tiago_thomson_polychromator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_profile_diagnostic, only: plasma_sample_t, Q_NE, Q_TE
    use tiago_recorded_diagnostic, only: recorded_diagnostic_t
    use tiago_calibration, only: calibration_prior_t
    use tiago_thomson_scattering, only: thomson_t, thomson_volume_t
    use tiago_thomson_spectrum, only: thomson_spectrum_t
    use tiago_ts_photoelectrons, only: ts_photoelectrons
    use tiago_ts_recorded, only: ts_recorded
    use tiago_ts_rayleigh, only: ts_rayleigh
    implicit none
    private
    type, public :: polychromator_channel_t
        character(:), allocatable :: name
        integer :: volume = 1
        real(dp) :: lambda0 = 0.0_dp, theta = 0.0_dp
        real(dp) :: energy = 0.0_dp, solid_angle = 0.0_dp
        real(dp), allocatable :: wavelength(:), weights(:), transmission(:), qe(:)
        ! Response at the laser line for an elastic Rayleigh calibration.
        real(dp) :: elastic_response = 0.0_dp
    end type
    type, extends(recorded_diagnostic_t), public :: recorded_thomson_t
        type(thomson_t) :: sampling
        type(thomson_spectrum_t) :: spectrum
        type(polychromator_channel_t), allocatable :: channels(:)
        type(calibration_prior_t) :: gas
        logical :: rayleigh_mode = .false.
    contains
        procedure :: init => poly_init
        procedure :: n_signals => poly_count
        procedure :: n_nuisance => poly_nuisance
        procedure :: label => poly_label
        procedure :: sample_points => poly_points
        procedure :: predict => poly_predict
        procedure :: photoelectrons => poly_photons
        procedure :: photon_jvp
        procedure :: photon_vjp
        procedure :: jvp => poly_jvp
        procedure :: vjp => poly_vjp
    end type
contains
    subroutine poly_init(self, volumes, channels, prior, status, aligned, gas)
        class(recorded_thomson_t), intent(inout) :: self
        type(thomson_volume_t), intent(in) :: volumes(:)
        type(polychromator_channel_t), intent(in) :: channels(:)
        type(calibration_prior_t), intent(in) :: prior
        integer, intent(out) :: status
        logical, optional, intent(in) :: aligned
        type(calibration_prior_t), optional, intent(in) :: gas
        integer :: i, n
        logical :: alignment
        status = 1
        n = size(channels)
        if (n == 0 .or. size(volumes) == 0) return
        if (prior%n_params() /= 3 + 2*n) return
        if (.not. allocated(prior%factor)) return
        if (any(prior%mean <= 0.0_dp)) return
        do i = 1, size(volumes)
            if (.not. all(ieee_is_finite(volumes(i)%center))) return
            if (.not. all(ieee_is_finite(volumes(i)%direction))) return
            if (norm2(volumes(i)%direction) <= 0.0_dp) return
            if (.not. ieee_is_finite(volumes(i)%length)) return
            if (volumes(i)%length <= 0.0_dp) return
        end do
        do i = 1, n
            if (channels(i)%volume < 1 .or. channels(i)%volume > size(volumes)) return
            if (.not. allocated(channels(i)%name)) return
            if (.not. all(ieee_is_finite([channels(i)%lambda0, channels(i)%theta, &
                                          channels(i)%energy, channels(i)%solid_angle, channels(i)%elastic_response]))) return
            if (channels(i)%lambda0 <= 0.0_dp .or. channels(i)%energy < 0.0_dp .or. &
                channels(i)%solid_angle < 0.0_dp) return
            if (channels(i)%theta <= 0.0_dp .or. channels(i)%theta > acos(-1.0_dp)) return
            if (channels(i)%elastic_response < 0.0_dp) return
            if (.not. allocated(channels(i)%wavelength)) return
            if (.not. allocated(channels(i)%weights)) return
            if (.not. allocated(channels(i)%transmission)) return
            if (.not. allocated(channels(i)%qe)) return
            if (size(channels(i)%wavelength) == 0) return
            if (size(channels(i)%weights) /= size(channels(i)%wavelength)) return
            if (size(channels(i)%qe) /= size(channels(i)%wavelength)) return
            if (size(channels(i)%transmission) /= size(channels(i)%wavelength)) return
            if (.not. all(ieee_is_finite(channels(i)%wavelength))) return
            if (.not. all(ieee_is_finite(channels(i)%weights))) return
            if (.not. all(ieee_is_finite(channels(i)%qe))) return
            if (.not. all(ieee_is_finite(channels(i)%transmission))) return
            if (any(channels(i)%wavelength <= 0.0_dp)) return
            if (any(channels(i)%weights < 0.0_dp)) return
            if (any(channels(i)%qe < 0.0_dp) .or. any(channels(i)%qe > 1.0_dp)) return
            if (any(channels(i)%transmission < 0.0_dp) .or. &
                any(channels(i)%transmission > 1.0_dp)) return
        end do
        self%rayleigh_mode = present(gas)
        if (present(gas)) then
            if (gas%n_params() /= 3) return
            if (.not. allocated(gas%factor)) return
            if (any(gas%mean <= 0.0_dp)) return
            self%gas = gas
        end if
        alignment = .true.
        if (present(aligned)) alignment = aligned
        call self%sampling%init(volumes, alignment)
        call self%spectrum%init()
        self%channels = channels
        self%calibration = prior
        status = 0
    end subroutine
    integer function poly_count(self)
        class(recorded_thomson_t), intent(in) :: self
        poly_count = size(self%channels)
    end function
    integer function poly_nuisance(self)
        class(recorded_thomson_t), intent(in) :: self
        poly_nuisance = self%sampling%n_nuisance() + self%n_calibration()
        if (self%rayleigh_mode) poly_nuisance = poly_nuisance + 3
    end function
    function poly_label(self, i) result(text)
        class(recorded_thomson_t), intent(in) :: self
        integer, intent(in) :: i
        character(:), allocatable :: text
        text = 'Thomson V s: '//self%channels(i)%name
        if (self%rayleigh_mode) text = 'Rayleigh V s: '//self%channels(i)%name
    end function
    subroutine poly_points(self, eta, x, dx_deta)
        class(recorded_thomson_t), intent(in) :: self
        real(dp), intent(in) :: eta(:)
        real(dp), allocatable, intent(out) :: x(:, :), dx_deta(:, :, :)
        real(dp), allocatable :: inner(:, :, :)
        integer :: n
        n = self%sampling%n_nuisance()
        call self%sampling%sample_points(eta(:n), x, inner)
        allocate (dx_deta(3, size(x, 2), self%n_nuisance()))
        dx_deta = 0.0_dp
        dx_deta(:, :, :n) = inner
    end subroutine
    subroutine poly_predict(self, sample, eta, y)
        class(recorded_thomson_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        call evaluate(self, sample, eta, y)
    end subroutine
    subroutine poly_photons(self, sample, eta, y)
        class(recorded_thomson_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        call evaluate(self, sample, eta, y, photons=.true.)
    end subroutine
    subroutine poly_jvp(self, sample, eta, q_dot, eta_dot, y_dot)
        class(recorded_thomson_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        call tangent(self, sample, eta, q_dot, eta_dot, y_dot, .false.)
    end subroutine
    subroutine photon_jvp(self, sample, eta, q_dot, eta_dot, y_dot)
        class(recorded_thomson_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        call tangent(self, sample, eta, q_dot, eta_dot, y_dot, .true.)
    end subroutine
    subroutine tangent(self, sample, eta, q_dot, eta_dot, y_dot, photons)
        class(recorded_thomson_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        real(dp), allocatable :: rows(:, :, :), ce(:, :), ge(:, :), dc(:)
        real(dp) :: y(size(y_dot)), dg(3)
        logical, intent(in) :: photons
        integer :: i, k, n, m
        n = self%sampling%n_nuisance()
        m = self%n_calibration()
        allocate (rows(2, size(self%sampling%w), size(y_dot)), ce(m, size(y_dot)), ge(3, size(y_dot)))
        call evaluate(self, sample, eta, y, rows, ce, ge, photons)
        dc = self%calibration%jvp(eta_dot(n + 1:n + m))
        dg = 0.0_dp
        if (self%rayleigh_mode) dg = self%gas%jvp(eta_dot(n + m + 1:))
        do i = 1, size(y_dot)
            y_dot(i) = dot_product(ce(:, i), dc) + dot_product(ge(:, i), dg)
            do k = 1, size(self%sampling%w)
                y_dot(i) = y_dot(i) + rows(1, k, i)*q_dot(Q_NE, k) + rows(2, k, i)*q_dot(Q_TE, k)
            end do
        end do
    end subroutine
    subroutine poly_vjp(self, sample, eta, y_bar, q_bar, eta_bar)
        class(recorded_thomson_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        call reverse(self, sample, eta, y_bar, q_bar, eta_bar, .false.)
    end subroutine
    subroutine photon_vjp(self, sample, eta, y_bar, q_bar, eta_bar)
        class(recorded_thomson_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        call reverse(self, sample, eta, y_bar, q_bar, eta_bar, .true.)
    end subroutine
    subroutine reverse(self, sample, eta, y_bar, q_bar, eta_bar, photons)
        class(recorded_thomson_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        real(dp), allocatable :: rows(:, :, :), ce(:, :), ge(:, :)
        real(dp) :: y(size(y_bar))
        logical, intent(in) :: photons
        integer :: i, k, n, m
        n = self%sampling%n_nuisance()
        m = self%n_calibration()
        allocate (rows(2, size(self%sampling%w), size(y)), ce(m, size(y)), ge(3, size(y)))
        call evaluate(self, sample, eta, y, rows, ce, ge, photons)
        q_bar = 0.0_dp
        eta_bar = 0.0_dp
        do i = 1, size(y)
            do k = 1, size(self%sampling%w)
                q_bar(Q_NE, k) = q_bar(Q_NE, k) + rows(1, k, i)*y_bar(i)
                q_bar(Q_TE, k) = q_bar(Q_TE, k) + rows(2, k, i)*y_bar(i)
            end do
        end do
        eta_bar(n + 1:n + m) = self%calibration%vjp(matmul(ce, y_bar))
        if (self%rayleigh_mode) eta_bar(n + m + 1:) = self%gas%vjp(matmul(ge, y_bar))
    end subroutine
    subroutine evaluate(self, sample, eta, y, rows, ce, ge, photons)
        class(recorded_thomson_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        real(dp), optional, intent(out) :: rows(:, :, :), ce(:, :), ge(:, :)
        logical, optional, intent(in) :: photons
        real(dp), allocatable :: c(:)
        real(dp) :: gas(3), d(11), f, band, dband, weight, response, gain
        integer :: i, k, v, j, n, m
        logical :: pe
        n = self%sampling%n_nuisance()
        m = self%n_calibration()
        c = self%calibration%physical(eta(n + 1:n + m))
        if (any(c <= 0.0_dp)) error stop 'Nonpositive Thomson calibration parameter'
        if (self%rayleigh_mode) then
            gas = self%gas%physical(eta(n + m + 1:))
            if (any(gas <= 0.0_dp)) error stop 'Nonpositive Rayleigh gas parameter'
        end if
        if (.not. self%rayleigh_mode) then
            if (.not. all(ieee_is_finite(sample%q(Q_NE, :)))) error stop 'Nonfinite Thomson ne'
            if (any(sample%q(Q_NE, :) < 0.0_dp)) error stop 'Negative Thomson ne'
            if (any(sample%q(Q_TE, :) < 0.0_dp)) error stop 'Negative Thomson Te'
        end if
        pe = .false.
        if (present(photons)) pe = photons
        if (present(rows)) rows = 0.0_dp
        if (present(ce)) ce = 0.0_dp
        if (present(ge)) ge = 0.0_dp
        y = 0.0_dp
        !$omp parallel do default(shared) private(i,k,v,j,d,f,band,dband,weight,response,gain) &
        !$omp schedule(static) if(size(y)>=8)
        do i = 1, size(y)
            v = self%channels(i)%volume
            j = 4 + 2*(i - 1)
            if (self%rayleigh_mode) then
                response = self%channels(i)%elastic_response*c(3)*c(j)
                gain = c(j + 1)
                if (pe) gain = 1.0_dp
                call ts_rayleigh(gas(1), gas(2), self%channels(i)%energy, &
                                 self%channels(i)%lambda0, self%channels(i)%solid_angle, &
                                 self%sampling%volumes(v)%length, response, c(1), c(2), gain, &
                                 gas(3), f, d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), &
                                 d(9), d(10), d(11))
                if (pe) d(10) = 0.0_dp
                y(i) = f
                if (present(ge)) ge(:, i) = d([1, 2, 11])
                if (present(ce)) then
                    ce(1, i) = d(8)
                    ce(2, i) = d(9)
                    ce(3, i) = d(7)*self%channels(i)%elastic_response*c(j)
                    ce(j, i) = d(7)*self%channels(i)%elastic_response*c(3)
                    ce(j + 1, i) = d(10)
                end if
                cycle
            end if
            do k = 1, size(self%sampling%w)
                if (self%sampling%owner(k) /= v) cycle
                if (sample%q(Q_NE, k) <= 0.0_dp .and. sample%q(Q_TE, k) <= 0.0_dp) cycle
                call channel_band(self, i, sample%q(Q_TE, k), band, dband)
                if (pe) then
                    call ts_photoelectrons(sample%q(Q_NE, k), band, self%channels(i)%energy, &
                                           self%channels(i)%lambda0, self%channels(i)%solid_angle, &
                                           self%sampling%volumes(v)%length, c(3), c(1), c(2), c(j), &
                                           f, d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), d(9), d(10))
                    d(11) = 0.0_dp
                else
                    call ts_recorded(sample%q(Q_NE, k), band, self%channels(i)%energy, &
                                     self%channels(i)%lambda0, self%channels(i)%solid_angle, &
                                     self%sampling%volumes(v)%length, c(3), c(1), c(2), c(j), c(j + 1), &
                                     f, d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), d(9), d(10), d(11))
                end if
                weight = self%sampling%w(k)
                y(i) = y(i) + weight*f
                if (present(rows)) rows(:, k, i) = weight*[d(1), d(2)*dband]
                if (present(ce)) then
                    ce(1, i) = ce(1, i) + weight*d(8)
                    ce(2, i) = ce(2, i) + weight*d(9)
                    ce(3, i) = ce(3, i) + weight*d(7)
                    ce(j, i) = ce(j, i) + weight*d(10)
                    ce(j + 1, i) = ce(j + 1, i) + weight*d(11)
                end if
            end do
        end do
        !$omp end parallel do
    end subroutine
    subroutine channel_band(self, i, te, band, dband)
        class(recorded_thomson_t), intent(in) :: self
        integer, intent(in) :: i
        real(dp), intent(in) :: te
        real(dp), intent(out) :: band, dband
        real(dp) :: s(size(self%channels(i)%wavelength)), ds(size(self%channels(i)%wavelength))
        real(dp) :: w(size(s))
        call self%spectrum%cross_section(te, self%channels(i)%lambda0, &
                                         self%channels(i)%theta, self%channels(i)%wavelength, s, ds)
        w = self%channels(i)%weights*self%channels(i)%transmission*self%channels(i)%qe
        band = dot_product(w, s)
        dband = dot_product(w, ds)
    end subroutine
end module
