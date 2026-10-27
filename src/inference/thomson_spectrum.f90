module tiago_thomson_spectrum
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_quadrature, only: gauss_gen_laguerre, gauss_legendre_ab
    use fortnum_special_complex_bessel, only: bessel_k_complex
    use fortnum_status, only: fortnum_status_t
    use tiago_ts_depolarization, only: ts_depolarization
    use tiago_ts_spectral, only: ts_spectral
    use tiago_ts_bessel_chain, only: ts_bessel_chain
    implicit none
    private
    type, public :: thomson_spectrum_t
        private
        real(dp), allocatable :: nodes(:), weights(:)
    contains
        procedure :: init => spectrum_init
        procedure :: cross_section => spectrum_cross_section
        procedure :: normalized => spectrum_normalized
    end type
contains
    subroutine spectrum_init(self, order)
        class(thomson_spectrum_t), intent(inout) :: self
        integer, optional, intent(in) :: order
        integer :: n
        n = 24
        if (present(order)) n = order
        if (n < 4) error stop 'Thomson: depolarization quadrature needs order >= 4'
        if (allocated(self%nodes)) deallocate (self%nodes, self%weights)
        allocate (self%nodes(n), self%weights(n))
        call gauss_gen_laguerre(n, 0.0_dp, self%nodes, self%weights)
    end subroutine

    subroutine spectrum_cross_section(self, te, lambda0, theta, wavelength, s, dsdt)
        !! Photon differential cross-section divided by re**2 [1/m].
        !! Incident and detected polarization perpendicular to scattering plane.
        class(thomson_spectrum_t), intent(in) :: self
        real(dp), intent(in) :: te, lambda0, theta, wavelength(:)
        real(dp), intent(out) :: s(:), dsdt(:)
        complex(dp) :: bk(3)
        type(fortnum_status_t) :: status
        real(dp) :: kdot, a, depol, ddt, value, deriv, dt, dk, dq, r, mu
        integer :: j, k
        if (.not. allocated(self%nodes)) error stop 'Thomson spectrum not initialized'
        if (.not. all(ieee_is_finite([te, lambda0, theta]))) &
            error stop 'Nonfinite Thomson input'
        if (te < 10.0_dp .or. te > 20000.0_dp) error stop 'Thomson Te outside 10--20000 eV'
        if (lambda0 <= 0.0_dp .or. theta <= 0.0_dp .or. theta > acos(-1.0_dp)) &
            error stop 'Invalid Thomson wavelength or scattering angle'
        if (.not. all(ieee_is_finite(wavelength))) error stop 'Nonfinite wavelength'
        if (any(wavelength <= 0.0_dp)) error stop 'Nonpositive wavelength'
        a = 510998.95069_dp/te
        do k = 1, 3
            call bessel_k_complex(k, cmplx(a, 0.0_dp, dp), .true., bk(k), status)
            if (status%code /= 0) error stop 'Thomson Bessel evaluation failed'
        end do
        call ts_bessel_chain(te, real(bk(1), dp), real(bk(2), dp), real(bk(3), dp), kdot)
        mu = cos(theta)
        do j = 1, size(wavelength)
            r = wavelength(j)/lambda0
            depol = 0.0_dp
            ddt = 0.0_dp
            do k = 1, size(self%nodes)
                call ts_depolarization(te, r, mu, self%nodes(k), value, deriv)
                depol = depol + self%weights(k)*value
                ddt = ddt + self%weights(k)*deriv
            end do
            call ts_spectral(te, r, mu, lambda0, real(bk(2), dp), depol, value, dt, dk, dq)
            s(j) = value
            dsdt(j) = dt + dk*kdot + dq*ddt
        end do
    end subroutine

    subroutine spectrum_normalized(self, te, lambda0, theta, wavelength, s, dsdt, area)
        class(thomson_spectrum_t), intent(in) :: self
        real(dp), intent(in) :: te, lambda0, theta, wavelength(:)
        real(dp), intent(out) :: s(:), dsdt(:)
        real(dp), optional, intent(out) :: area
        real(dp), parameter :: edge(21) = [-4.0_dp, -2.0_dp, -1.0_dp, -.5_dp, &
                                           -.25_dp, -.125_dp, -.0625_dp, -.03125_dp, -.015625_dp, -.0078125_dp, &
                                           0.0_dp, .0078125_dp, .015625_dp, .03125_dp, .0625_dp, .125_dp, &
                                           .25_dp, .5_dp, 1.0_dp, 2.0_dp, 4.0_dp]
        real(dp) :: z(24), w(24), wave(480), weights(480), raw(480), dr(480)
        real(dp) :: norm, dnorm
        integer :: k, lo, hi
        do k = 1, 20
            lo = (k - 1)*24 + 1
            hi = k*24
            call gauss_legendre_ab(24, edge(k), edge(k + 1), z, w)
            wave(lo:hi) = lambda0*exp(z*sqrt(1 - cos(theta)))
            weights(lo:hi) = w*sqrt(1 - cos(theta))*wave(lo:hi)
        end do
        call self%cross_section(te, lambda0, theta, wave, raw, dr)
        norm = dot_product(weights, raw)
        dnorm = dot_product(weights, dr)
        call self%cross_section(te, lambda0, theta, wavelength, s, dsdt)
        dsdt = (dsdt - s*dnorm/norm)/norm
        s = s/norm
        if (present(area)) area = norm
    end subroutine
end module
