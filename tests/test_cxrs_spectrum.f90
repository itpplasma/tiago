module test_cxrs_spectrum_plasma
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_cxrs_spectrum, only: cxrs_plasma_t
    implicit none

    type, extends(cxrs_plasma_t) :: uniform_t
        !! Uniform plasma for x(1) > 1 m. ne is tiny so the beam is not
        !! attenuated, and Ti = 0 removes the halo unless ti is set.
        real(dp) :: ne = 1.0e6_dp, te = 1.0e3_dp, ti = 0.0_dp, zeff = 1.5_dp
        real(dp) :: n_z = 1.0e18_dp, t_z = 2.0e3_dp
        real(dp) :: flow(3) = [3.0e4_dp, -1.0e5_dp, 2.0e4_dp], field(3) = [0.3_dp, 2.0_dp, -0.4_dp]
    contains
        procedure :: eval => uniform_eval
        procedure :: impurity => uniform_impurity
    end type uniform_t
contains
    subroutine uniform_eval(self, x, ne, te, ti, zeff)
        class(uniform_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: ne, te, ti, zeff
        ne = 0.0_dp
        te = 0.0_dp
        ti = 0.0_dp
        zeff = 1.0_dp
        if (x(1) <= 1.0_dp) return
        ne = self%ne
        te = self%te
        ti = self%ti
        zeff = self%zeff
    end subroutine uniform_eval

    subroutine uniform_impurity(self, x, n_z, t_z, flow, field)
        class(uniform_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: n_z, t_z, flow(3), field(3)
        n_z = 0.0_dp
        t_z = 0.0_dp
        flow = self%flow
        field = self%field
        if (x(1) <= 1.0_dp) return
        n_z = self%n_z
        t_z = self%t_z
    end subroutine uniform_impurity
end module test_cxrs_spectrum_plasma

program test_cxrs_spectrum
    !! CXRS spectra and line moments (tiago_cxrs_spectrum).
    !! Oracles: Gauss-Hermite and Gauss-Laguerre rules on monomials; the
    !! hydrogen 2p-1s radial integral 128 sqrt(6) / 243 a0; Wigner 6j values
    !! 1/6 and -3/70; the hydrogen 2p3/2 - 2p1/2 Dirac splitting against its
    !! (Z alpha)^4 term and the 2:1 Lyman-alpha intensities; the cross-section
    !! moments against closed-form Gaussian moments for q linear in the
    !! collision energy; the brightness of a Gaussian beam seen across its
    !! axis in closed form; apparent velocity and temperature equal to the
    !! flow and T_Z for constant q; the moments of the synthetic spectrum
    !! against the channel moments with fine structure, Zeeman effect and
    !! instrument function.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_cxrs_spectrum
    use tiago_beam_neutrals, only: neutral_beam_t, beam_component_flux, beam_component_speed, &
        BEAM_COMPONENTS
    use test_cxrs_spectrum_plasma, only: uniform_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), c_light = 299792458.0_dp
    real(dp), parameter :: e_charge = 1.602176634e-19_dp, amu = 1.66053906660e-27_dp
    integer :: failures

    failures = 0
    call check_quadrature()
    call check_atomic()
    call check_cross_section_moments()
    call check_channel()
    if (failures > 0) then
        write (error_unit, '(I0,A)') failures, ' CXRS spectrum check(s) failed'
        error stop 1
    end if
    print '(A)', 'CXRS spectrum: all checks passed'

contains

    subroutine expect(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            failures = failures + 1
            write (error_unit, '(A)') 'FAIL: '//what
        end if
    end subroutine expect

    subroutine check_quadrature()
        real(dp) :: xh(16), wh(16), xl(10), wl(10)
        call gauss_hermite(16, xh, wh)
        call gauss_laguerre(10, xl, wl)
        call expect(abs(sum(wh * xh**2) - sqrt(pi) / 2) < 1.0e-13_dp, 'Gauss-Hermite x^2')
        call expect(abs(sum(wh * xh**8) - 105 * sqrt(pi) / 16) < 1.0e-12_dp, 'Gauss-Hermite x^8')
        call expect(abs(sum(wl * xl**5) / 120 - 1) < 1.0e-13_dp, 'Gauss-Laguerre x^5')
        call expect(abs(sum(wl * xl**19) / gamma(20.0_dp) - 1) < 1.0e-11_dp, 'Gauss-Laguerre x^19')
    end subroutine check_quadrature

    subroutine check_atomic()
        type(cxrs_line_t) :: line
        real(dp), allocatable :: shift(:), weight(:)
        real(dp) :: split
        call expect(abs(abs(hydrogenic_radial_integral(2, 1, 1, 0)) - 128 * sqrt(6.0_dp) / 243) &
            < 1.0e-12_dp, 'hydrogen 2p-1s radial integral')
        call expect(abs(wigner_6j(2, 2, 2, 2, 2, 2) - 1.0_dp / 6) < 1.0e-14_dp, '6j {1 1 1; 1 1 1}')
        call expect(abs(wigner_6j(4, 4, 4, 4, 4, 4) + 3.0_dp / 70) < 1.0e-14_dp, '6j {2 2 2; 2 2 2}')
        ! Leading order mu c^2 (Z alpha)^4 / 32 for n = 2, j = 3/2 - 1/2.
        split = dirac_energy(1, 2, 3, 1.00727646688_dp) - dirac_energy(1, 2, 1, 1.00727646688_dp)
        call expect(abs(split / (510998.95_dp / (1 + 5.48579909065e-4_dp / 1.00727646688_dp) &
            * 7.2973525693e-3_dp**4 / 32) - 1) < 1.0e-4_dp, 'hydrogen 2p fine structure')
        line%z = 1
        line%n_upper = 2
        line%n_lower = 1
        line%mass_amu = 1.00727646688_dp
        line%wavelength = 121.567e-9_dp
        call fine_structure_pattern(line, shift, weight)
        call expect(size(weight) == 2, 'Lyman-alpha has two components')
        call expect(abs(maxval(weight) - 2.0_dp / 3) < 1.0e-12_dp, 'Lyman-alpha 2:1 intensities')
    end subroutine check_atomic

    subroutine check_cross_section_moments()
        type(cx_rate_table_t) :: table
        real(dp) :: xh(16), wh(16), xl(10), wl(10), sigma, drift, a, m0, m1, m2, kappa, q0, q1
        real(dp) :: r0, r1, r2
        call gauss_hermite(16, xh, wh)
        call gauss_laguerre(10, xl, wl)
        q0 = 2.0e-15_dp
        q1 = 3.0e-20_dp
        table%energy = [0.0_dp, 1.0e7_dp]
        table%rate = [q0, q0 + q1 * 1.0e7_dp]
        kappa = 0.5_dp * amu / e_charge
        sigma = sqrt(2.0e3_dp * e_charge / (14.003074_dp * amu))
        drift = 2.2e6_dp
        a = 0.37_dp
        call cx_velocity_moments(table, xh, wh, xl, wl, sigma, drift, a, m0, m1, m2)
        r0 = q0 + q1 * kappa * (drift**2 + 3 * sigma**2)
        r1 = -a * 2 * q1 * kappa * drift * sigma**2
        r2 = q0 * sigma**2 + q1 * kappa * (5 * sigma**4 + drift**2 * sigma**2)
        call expect(abs(m0 / r0 - 1) < 1.0e-12_dp, 'cross-section moment 0')
        call expect(abs(m1 / r1 - 1) < 1.0e-10_dp, 'cross-section moment 1')
        call expect(abs(m2 / r2 - 1) < 1.0e-10_dp, 'cross-section moment 2')
    end subroutine check_cross_section_moments

    function test_beam() result(beam)
        type(neutral_beam_t) :: beam
        beam%origin = [0.0_dp, 0.0_dp, 0.0_dp]
        beam%direction = [1.0_dp, 0.0_dp, 0.0_dp]
        beam%energy = 60.0e3_dp
        beam%power = 2.0e6_dp
        beam%power_fraction = [0.6_dp, 0.3_dp, 0.1_dp]
        beam%width = [0.12_dp, 0.18_dp]
        beam%length = 4.0_dp
        beam%nodes = 401
    end function test_beam

    subroutine check_channel()
        type(cxrs_spectrometer_t) :: spec
        type(cxrs_line_t) :: line
        type(cxrs_channel_t) :: channel
        type(uniform_t) :: plasma
        type(cxrs_moments_t) :: mom
        type(neutral_beam_t) :: beam(1)
        real(dp) :: lam(6001), s(6001), dl, b, l0, ref, total, mean, var, e(3)
        integer :: status, j, i
        beam(1) = test_beam()
        line%rate%energy = [1.0_dp]
        line%rate%rate = [5.0e-15_dp]
        line%fine_structure = .false.
        line%zeeman = .false.
        call spec%init(line, beam, plasma, status)
        call expect(status == 0, 'spectrometer init')
        ! Vertical chord crossing the beam axis at x = 2.5 m (ev = z).
        channel%chord%origin = [2.5_dp, 0.0_dp, 1.5_dp]
        channel%chord%direction = [0.0_dp, 0.0_dp, -1.0_dp]
        channel%chord%length = 3.0_dp
        channel%chord%npanel = 300
        mom = spec%moments(plasma, channel)
        ref = 0.0_dp
        do j = 1, BEAM_COMPONENTS
            ref = ref + beam_component_flux(beam(1), j) / (beam_component_speed(beam(1), j) &
                * pi * 0.12_dp * 0.18_dp) * sqrt(pi) * 0.18_dp
        end do
        ref = ref * plasma%n_z * 5.0e-15_dp / (4 * pi)
        call expect(abs(mom%brightness / ref - 1) < 1.0e-10_dp, 'brightness of a Gaussian beam')
        call expect(abs(mom%velocity - dot_product(plasma%flow, [0.0_dp, 0.0_dp, -1.0_dp])) &
            < 1.0e-6_dp, 'velocity for constant q')
        call expect(abs(mom%temperature / plasma%t_z - 1) < 1.0e-10_dp, 'temperature for constant q')

        ! Energy-dependent q, halo, fine structure, Zeeman and instrument.
        line%rate%energy = [0.0_dp, 5.0e3_dp, 2.0e4_dp, 4.0e4_dp]
        line%rate%rate = [1.0e-16_dp, 2.0e-15_dp, 6.0e-15_dp, 4.0e-15_dp]
        line%fine_structure = .true.
        line%zeeman = .true.
        plasma%ti = 2.5e3_dp
        plasma%ne = 4.0e19_dp
        call spec%init(line, beam, plasma, status)
        channel%chord%npanel = 60
        channel%chord%direction = [0.4_dp, 0.3_dp, -1.0_dp]
        channel%chord%origin = [2.5_dp - 0.6_dp, -0.45_dp, 1.5_dp]
        channel%instrument_fwhm = 0.05e-9_dp
        mom = spec%moments(plasma, channel)
        call expect(mom%halo_fraction > 0.0_dp .and. mom%halo_fraction < 1.0_dp, 'halo fraction')
        l0 = line%wavelength
        dl = 6.0e-9_dp / (size(lam) - 1)
        do i = 1, size(lam)
            lam(i) = l0 - 3.0e-9_dp + (i - 1) * dl
        end do
        s = spec%spectrum(plasma, channel, lam)
        total = sum(s) * dl
        mean = sum(s * lam) * dl / total
        var = sum(s * (lam - mean)**2) * dl / total
        call expect(abs(total / mom%brightness - 1) < 1.0e-8_dp, 'spectrum integral')
        call expect(abs((mean / l0 - 1) * c_light - mom%velocity) < 1.0e-3_dp, 'spectrum centroid')
        b = 14.003074_dp * amu * c_light**2 * var / l0**2 / e_charge
        call expect(abs(b / mom%temperature - 1) < 1.0e-6_dp, 'spectrum variance')
        call expect(mom%temperature > mom%temperature_doppler, 'fine structure and Zeeman broaden')
        e = channel%chord%direction / norm2(channel%chord%direction)
        call expect(abs(mom%velocity - mom%velocity_local) > 1.0_dp, 'cross-section shift is present')
    end subroutine check_channel
end program test_cxrs_spectrum
