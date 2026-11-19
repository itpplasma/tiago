module tiago_cxrs_spectrum
    !! Charge-exchange recombination spectra of a hydrogen-like impurity line
    !! along straight sight lines crossing neutral beams (tiago_beam_neutrals).
    !!
    !! A fully stripped impurity of charge Z, density n_Z, temperature T_Z and
    !! flow V captures from beam component j (velocity u_j, density n_j) and
    !! from the halo (density n_h). With the effective emission rate q(E) of
    !! the observed n -> n' line, tabulated in the collision energy per amu E,
    !! the local emission of donor d is
    !!
    !!   eps_d f_d(v_l) = n_Z n_d / (4 pi) int f_Z(v) q(E(|v - u_d|)) delta(v_l - v . e) d^3v,
    !!
    !! e the unit vector of the sight line pointing away from the observer.
    !! The weight q(|v - u_d|) of the Maxwellian f_Z shifts and narrows the
    !! velocity distribution of emitting ions (cross-section effect, von
    !! Hellermann et al., Plasma Phys. Control. Fusion 37 (1995) 71). Its
    !! zeroth, first and second moments in v_l are integrated by Gauss-Hermite
    !! quadrature along u_d - V and Gauss-Laguerre quadrature in the
    !! perpendicular speed squared. Halo donors are thermal at the ion
    !! temperature and flow, so their weight is q at the mean collision energy
    !! E = 1.5 T_i (1/A_D + 1/A_Z) and carries no velocity dependence.
    !!
    !! The rest line is the zero-field fine-structure pattern of the n -> n'
    !! transition (Dirac energies, statistically populated upper j levels,
    !! hydrogenic line strengths) convolved with the normal Zeeman triplet at
    !! |B|; the pattern is centred on its intensity-weighted mean. Each sight
    !! point contributes Gaussians with its local mean and variance, and the
    !! instrument function is a Gaussian. Line moments of a channel are the
    !! emission-weighted sums over the sight line; the apparent temperature is
    !! m_Z c^2 var(lambda) / lambda_0^2.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_sight_geometry, only: chord_t, gauss_legendre01
    use tiago_beam_neutrals, only: neutral_beam_t, beam_plasma_t, beam_attenuation_t, &
        BEAM_COMPONENTS, beam_axis_nodes, beam_attenuation, beam_density, beam_frame, &
        beam_component_speed, beam_component_flux, beam_optical_depth, halo_density, &
        hydrogen_charge_exchange, beam_check
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: c_light = 299792458.0_dp
    real(dp), parameter :: e_charge = 1.602176634e-19_dp
    real(dp), parameter :: amu = 1.66053906660e-27_dp
    real(dp), parameter :: m_e = 9.1093837015e-31_dp
    real(dp), parameter :: alpha_fs = 7.2973525693e-3_dp
    real(dp), parameter :: h_planck = 6.62607015e-34_dp
    real(dp), parameter :: mu_bohr = 9.2740100783e-24_dp
    integer, parameter :: n_hermite = 16, n_laguerre = 10

    type, public :: cx_rate_table_t
        !! Effective emission rate coefficient q [m^3/s] of the observed line
        !! against collision energy [eV/amu], piecewise linear in energy and
        !! constant beyond the ends.
        real(dp), allocatable :: energy(:)
        real(dp), allocatable :: rate(:)
    end type cx_rate_table_t

    type, public :: cxrs_line_t
        integer :: z = 7 !! nuclear charge of the hydrogen-like emitter
        integer :: n_upper = 9
        integer :: n_lower = 8
        real(dp) :: mass_amu = 14.003074_dp !! emitter mass [u]
        real(dp) :: wavelength = 566.937e-9_dp !! rest wavelength of the analysis [m]
        type(cx_rate_table_t) :: rate
        logical :: fine_structure = .true.
        logical :: zeeman = .true.
    end type cxrs_line_t

    type, public :: cxrs_channel_t
        type(chord_t) :: chord !! origin at the observer, direction into the plasma
        real(dp) :: instrument_fwhm = 0.0_dp !! Gaussian instrument function [m]
    end type cxrs_channel_t

    type, public :: cxrs_moments_t
        real(dp) :: brightness = 0.0_dp !! photons m^-2 s^-1 sr^-1
        real(dp) :: halo_fraction = 0.0_dp !! halo share of the brightness
        real(dp) :: velocity = 0.0_dp !! apparent line-of-sight velocity [m/s], away from observer
        real(dp) :: temperature = 0.0_dp !! apparent T from the full variance [eV]
        real(dp) :: temperature_doppler = 0.0_dp !! apparent T, Doppler variance only [eV]
        real(dp) :: velocity_local = 0.0_dp !! emission-weighted V . e [m/s], no cross-section effect
        real(dp) :: temperature_local = 0.0_dp !! emission-weighted T_Z [eV]
        real(dp) :: location(3) = 0.0_dp !! emission-weighted position, Cartesian [m]
    end type cxrs_moments_t

    type, abstract, public, extends(beam_plasma_t) :: cxrs_plasma_t
        !! Adds the emitting impurity: density [m^-3], temperature [eV], flow
        !! velocity and magnetic field, both Cartesian [m/s], [T].
    contains
        procedure(impurity_eval), deferred :: impurity
    end type cxrs_plasma_t

    abstract interface
        subroutine impurity_eval(self, x, n_z, t_z, flow, field)
            import :: cxrs_plasma_t, dp
            class(cxrs_plasma_t), intent(in) :: self
            real(dp), intent(in) :: x(3)
            real(dp), intent(out) :: n_z, t_z, flow(3), field(3)
        end subroutine impurity_eval
    end interface

    type, public :: cxrs_spectrometer_t
        type(cxrs_line_t) :: line
        type(neutral_beam_t), allocatable :: beams(:)
        type(beam_attenuation_t), allocatable :: att(:)
        real(dp), allocatable :: pattern_shift(:) !! fine-structure offsets [m]
        real(dp), allocatable :: pattern_weight(:) !! normalized intensities
        real(dp) :: xh(n_hermite) = 0.0_dp, wh(n_hermite) = 0.0_dp
        real(dp) :: xl(n_laguerre) = 0.0_dp, wl(n_laguerre) = 0.0_dp
    contains
        procedure :: init => spectrometer_init
        procedure :: moments => spectrometer_moments
        procedure :: spectrum => spectrometer_spectrum
    end type cxrs_spectrometer_t

    public :: cx_rate, cx_velocity_moments, fine_structure_pattern, zeeman_variance, &
        hydrogenic_radial_integral, wigner_6j, dirac_energy, gauss_hermite, gauss_laguerre

contains

    pure real(dp) function cx_rate(table, e)
        !! q(E) [m^3/s] at the collision energy E [eV/amu].
        type(cx_rate_table_t), intent(in) :: table
        real(dp), intent(in) :: e
        integer :: k, n
        n = size(table%energy)
        if (e <= table%energy(1)) then
            cx_rate = table%rate(1)
        else if (e >= table%energy(n)) then
            cx_rate = table%rate(n)
        else
            k = 1
            do while (table%energy(k + 1) < e)
                k = k + 1
            end do
            cx_rate = table%rate(k) + (table%rate(k + 1) - table%rate(k)) * (e - table%energy(k)) &
                / (table%energy(k + 1) - table%energy(k))
        end if
    end function cx_rate

    pure subroutine gauss_hermite(n, x, w)
        !! Nodes and weights of int exp(-x^2) f(x) dx.
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)
        real(dp) :: z, z1, p1, p2, p3, pp
        integer :: i, j, it
        z = 0.0_dp
        do i = 1, (n + 1) / 2
            if (i == 1) then
                z = sqrt(real(2 * n + 1, dp)) - 1.85575_dp * real(2 * n + 1, dp)**(-1.0_dp / 6.0_dp)
            else if (i == 2) then
                z = z - 1.14_dp * real(n, dp)**0.426_dp / z
            else if (i == 3) then
                z = 1.86_dp * z - 0.86_dp * x(1)
            else if (i == 4) then
                z = 1.91_dp * z - 0.91_dp * x(2)
            else
                z = 2.0_dp * z - x(i - 2)
            end if
            do it = 1, 100
                p1 = pi**(-0.25_dp)
                p2 = 0.0_dp
                do j = 1, n
                    p3 = p2
                    p2 = p1
                    p1 = z * sqrt(2.0_dp / j) * p2 - sqrt(real(j - 1, dp) / j) * p3
                end do
                pp = sqrt(2.0_dp * n) * p2
                z1 = z
                z = z1 - p1 / pp
                if (abs(z - z1) <= 1.0e-14_dp) exit
            end do
            x(i) = z
            x(n + 1 - i) = -z
            w(i) = 2.0_dp / (pp * pp)
            w(n + 1 - i) = w(i)
        end do
    end subroutine gauss_hermite

    pure subroutine gauss_laguerre(n, x, w)
        !! Nodes and weights of int_0^inf exp(-x) f(x) dx.
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)
        real(dp) :: z, z1, p1, p2, p3, pp
        integer :: i, j, it
        z = 0.0_dp
        do i = 1, n
            if (i == 1) then
                z = 3.0_dp / (1.0_dp + 2.4_dp * n)
            else if (i == 2) then
                z = z + 15.0_dp / (1.0_dp + 2.5_dp * n)
            else
                z = z + (1.0_dp + 2.55_dp * (i - 2)) / (1.9_dp * (i - 2)) * (z - x(i - 2))
            end if
            do it = 1, 100
                p1 = 1.0_dp
                p2 = 0.0_dp
                do j = 1, n
                    p3 = p2
                    p2 = p1
                    p1 = ((2 * j - 1 - z) * p2 - (j - 1) * p3) / j
                end do
                pp = n * (p1 - p2) / z
                z1 = z
                z = z1 - p1 / pp
                if (abs(z - z1) <= 1.0e-14_dp * max(1.0_dp, z)) exit
            end do
            x(i) = z
            w(i) = -1.0_dp / (pp * n * p2)
        end do
    end subroutine gauss_laguerre

    pure subroutine cx_velocity_moments(table, xh, wh, xl, wl, sigma, drift, a, m0, m1, m2)
        !! Moments int f q, int f q w_l, int f q w_l^2 of the Maxwellian
        !! deviation w (per-dimension variance sigma^2) weighted with
        !! q(E(|w - drift e1|)), w_l = a w_1 + b w_2, a^2 + b^2 = 1, e1 along
        !! the donor velocity relative to the flow, drift >= 0 [m/s].
        type(cx_rate_table_t), intent(in) :: table
        real(dp), intent(in) :: xh(:), wh(:), xl(:), wl(:), sigma, drift, a
        real(dp), intent(out) :: m0, m1, m2
        real(dp), parameter :: kappa = 0.5_dp * amu / e_charge
        real(dp) :: w1, wp2, q, s1, s11, spp
        integer :: i, k
        m0 = 0.0_dp
        s1 = 0.0_dp
        s11 = 0.0_dp
        spp = 0.0_dp
        do i = 1, size(xh)
            w1 = sqrt(2.0_dp) * sigma * xh(i)
            do k = 1, size(xl)
                wp2 = 2.0_dp * sigma * sigma * xl(k)
                q = wh(i) / sqrt(pi) * wl(k) * cx_rate(table, kappa * ((w1 - drift)**2 + wp2))
                m0 = m0 + q
                s1 = s1 + q * w1
                s11 = s11 + q * w1 * w1
                spp = spp + q * wp2
            end do
        end do
        m1 = a * s1
        m2 = a * a * s11 + 0.5_dp * (1.0_dp - a * a) * spp
    end subroutine cx_velocity_moments

    pure real(dp) function dirac_energy(z, n, j2, mass_amu)
        !! Binding energy [eV] (negative) of a hydrogen-like level n, j = j2/2,
        !! Dirac formula with the reduced electron mass.
        integer, intent(in) :: z, n, j2
        real(dp), intent(in) :: mass_amu
        real(dp) :: za, k, mu
        za = z * alpha_fs
        k = 0.5_dp * (j2 + 1)
        mu = m_e * mass_amu * amu / (m_e + mass_amu * amu)
        dirac_energy = mu * c_light**2 / e_charge &
            * (1.0_dp / sqrt(1.0_dp + (za / (n - k + sqrt(k * k - za * za)))**2) - 1.0_dp)
    end function dirac_energy

    pure real(dp) function log_factorial(n)
        integer, intent(in) :: n
        log_factorial = log_gamma(real(n + 1, dp))
    end function log_factorial

    pure real(dp) function triangle(a2, b2, c2)
        !! Racah Delta(abc) for doubled arguments.
        integer, intent(in) :: a2, b2, c2
        triangle = exp(0.5_dp * (log_factorial((a2 + b2 - c2) / 2) + log_factorial((a2 - b2 + c2) / 2) &
            + log_factorial((-a2 + b2 + c2) / 2) - log_factorial((a2 + b2 + c2) / 2 + 1)))
    end function triangle

    pure logical function triad(a2, b2, c2)
        integer, intent(in) :: a2, b2, c2
        triad = c2 >= abs(a2 - b2) .and. c2 <= a2 + b2 .and. mod(a2 + b2 + c2, 2) == 0
    end function triad

    pure real(dp) function wigner_6j(j1, j2, j3, j4, j5, j6)
        !! {j1 j2 j3; j4 j5 j6} for doubled angular momenta (Racah formula).
        integer, intent(in) :: j1, j2, j3, j4, j5, j6
        integer :: t, tmin, tmax
        real(dp) :: s
        wigner_6j = 0.0_dp
        if (.not. (triad(j1, j2, j3) .and. triad(j1, j5, j6) .and. triad(j4, j2, j6) &
            .and. triad(j4, j5, j3))) return
        tmin = max(j1 + j2 + j3, j1 + j5 + j6, j4 + j2 + j6, j4 + j5 + j3) / 2
        tmax = min(j1 + j2 + j4 + j5, j2 + j3 + j5 + j6, j3 + j1 + j6 + j4) / 2
        s = 0.0_dp
        do t = tmin, tmax
            s = s + (-1)**t * exp(log_factorial(t + 1) - log_factorial(t - (j1 + j2 + j3) / 2) &
                - log_factorial(t - (j1 + j5 + j6) / 2) - log_factorial(t - (j4 + j2 + j6) / 2) &
                - log_factorial(t - (j4 + j5 + j3) / 2) - log_factorial((j1 + j2 + j4 + j5) / 2 - t) &
                - log_factorial((j2 + j3 + j5 + j6) / 2 - t) - log_factorial((j3 + j1 + j6 + j4) / 2 - t))
        end do
        wigner_6j = s * triangle(j1, j2, j3) * triangle(j1, j5, j6) * triangle(j4, j2, j6) &
            * triangle(j4, j5, j3)
    end function wigner_6j

    pure real(dp) function laguerre(k, alpha, x)
        !! Generalized Laguerre polynomial L_k^(alpha)(x).
        integer, intent(in) :: k, alpha
        real(dp), intent(in) :: x
        real(dp) :: p0, p1, p2
        integer :: i
        p0 = 1.0_dp
        if (k == 0) then
            laguerre = p0
            return
        end if
        p1 = 1.0_dp + alpha - x
        do i = 1, k - 1
            p2 = ((2 * i + 1 + alpha - x) * p1 - (i + alpha) * p0) / (i + 1)
            p0 = p1
            p1 = p2
        end do
        laguerre = p1
    end function laguerre

    pure real(dp) function hydrogenic_radial_integral(n1, l1, n2, l2)
        !! int R_n1l1 R_n2l2 r^3 dr in units of a0 / Z (Gauss-Laguerre, exact
        !! for these polynomial-times-exponential integrands).
        integer, intent(in) :: n1, l1, n2, l2
        integer, parameter :: nq = 40
        real(dp) :: x(nq), w(nq), r, c, s, s1, s2
        integer :: i
        call gauss_laguerre(nq, x, w)
        s1 = 0.0_dp
        s2 = 0.0_dp
        s = 0.0_dp
        do i = 1, nq
            r = x(i) * n1 / 2.0_dp
            s1 = s1 + w(i) * n1 / 2.0_dp * radial(n1, l1, r)**2 * r * r * exp(x(i))
            r = x(i) * n2 / 2.0_dp
            s2 = s2 + w(i) * n2 / 2.0_dp * radial(n2, l2, r)**2 * r * r * exp(x(i))
            c = 1.0_dp / n1 + 1.0_dp / n2
            r = x(i) / c
            s = s + w(i) / c * radial(n1, l1, r) * radial(n2, l2, r) * r**3 * exp(x(i))
        end do
        hydrogenic_radial_integral = s / sqrt(s1 * s2)
    contains
        pure real(dp) function radial(n, l, rr)
            integer, intent(in) :: n, l
            real(dp), intent(in) :: rr
            radial = rr**l * exp(-rr / n) * laguerre(n - l - 1, 2 * l + 1, 2.0_dp * rr / n)
        end function radial
    end function hydrogenic_radial_integral

    subroutine fine_structure_pattern(line, shift, weight)
        !! Zero-field components n l j -> n' l' j': wavelength offsets [m] from
        !! the intensity-weighted mean and normalized intensities, for
        !! statistically populated upper levels (intensity proportional to the
        !! line strength (2j+1)(2j'+1) {l j 1/2; j' l' 1}^2 max(l,l') R^2).
        type(cxrs_line_t), intent(in) :: line
        real(dp), allocatable, intent(out) :: shift(:), weight(:)
        real(dp) :: lam(4 * line%n_upper * 4), s(4 * line%n_upper * 4), de, r
        integer :: l1, l2, j1, j2, m
        m = 0
        do l1 = 0, line%n_upper - 1
            do l2 = max(l1 - 1, 0), min(l1 + 1, line%n_lower - 1)
                if (abs(l1 - l2) /= 1) cycle
                r = hydrogenic_radial_integral(line%n_upper, l1, line%n_lower, l2)
                do j1 = max(2 * l1 - 1, 1), 2 * l1 + 1, 2
                    do j2 = max(2 * l2 - 1, 1), 2 * l2 + 1, 2
                        if (abs(j1 - j2) > 2) cycle
                        m = m + 1
                        s(m) = (j1 + 1) * (j2 + 1) * wigner_6j(2 * l1, j1, 1, j2, 2 * l2, 2)**2 &
                            * max(l1, l2) * r * r
                        de = dirac_energy(line%z, line%n_upper, j1, line%mass_amu) &
                            - dirac_energy(line%z, line%n_lower, j2, line%mass_amu)
                        lam(m) = h_planck * c_light / (de * e_charge)
                    end do
                end do
            end do
        end do
        weight = s(:m) / sum(s(:m))
        shift = lam(:m) - sum(weight * lam(:m))
    end subroutine fine_structure_pattern

    pure real(dp) function zeeman_variance(line, field, e)
        !! Wavelength variance [m^2] of the normal Zeeman triplet at |B| for a
        !! sight line e: sigma components at +-lambda^2 mu_B B / (h c) carry
        !! the fraction (1 + cos^2 theta) / 2 of the intensity.
        type(cxrs_line_t), intent(in) :: line
        real(dp), intent(in) :: field(3), e(3)
        real(dp) :: b, c2
        b = norm2(field)
        zeeman_variance = 0.0_dp
        if (b <= 0.0_dp) return
        c2 = (dot_product(field, e) / (b * norm2(e)))**2
        zeeman_variance = 0.5_dp * (1.0_dp + c2) * (line%wavelength**2 * mu_bohr * b / (h_planck * c_light))**2
    end function zeeman_variance

    subroutine spectrometer_init(self, line, beams, plasma, status)
        !! status: beam_check of the first failing beam, 5 if a beam origin is
        !! inside the plasma, 6 if the rate table is empty or not increasing.
        class(cxrs_spectrometer_t), intent(inout) :: self
        type(cxrs_line_t), intent(in) :: line
        type(neutral_beam_t), intent(in) :: beams(:)
        class(cxrs_plasma_t), intent(in) :: plasma
        integer, intent(out) :: status
        real(dp), allocatable :: axis(:, :), ne(:), te(:), ti(:), zeff(:)
        integer :: b, k
        status = 6
        if (.not. allocated(line%rate%energy)) return
        if (size(line%rate%energy) < 1 .or. size(line%rate%energy) /= size(line%rate%rate)) return
        if (size(line%rate%energy) > 1) then
            if (any(line%rate%energy(2:) <= line%rate%energy(:size(line%rate%energy) - 1))) return
        end if
        self%line = line
        self%beams = beams
        if (allocated(self%att)) deallocate(self%att)
        allocate(self%att(size(beams)))
        do b = 1, size(beams)
            status = beam_check(beams(b))
            if (status /= 0) return
            axis = beam_axis_nodes(beams(b))
            allocate(ne(beams(b)%nodes), te(beams(b)%nodes), ti(beams(b)%nodes), zeff(beams(b)%nodes))
            do k = 1, beams(b)%nodes
                call plasma%eval(axis(:, k), ne(k), te(k), ti(k), zeff(k))
            end do
            status = 5
            if (ne(1) > 0.0_dp) return
            call beam_attenuation(beams(b), ne, te, zeff, self%att(b))
            deallocate(ne, te, ti, zeff)
        end do
        status = 0
        if (line%fine_structure) then
            call fine_structure_pattern(line, self%pattern_shift, self%pattern_weight)
        else
            self%pattern_shift = [0.0_dp]
            self%pattern_weight = [1.0_dp]
        end if
        call gauss_hermite(n_hermite, self%xh, self%wh)
        call gauss_laguerre(n_laguerre, self%xl, self%wl)
    end subroutine spectrometer_init

    subroutine local_emission(self, plasma, x, e, eps, mean, var, halo, tz, flow_l, field)
        !! Emission [photons m^-3 s^-1 sr^-1], mean and variance of v . e of the
        !! emitting ions at x, summed over beams, components and the halo.
        class(cxrs_spectrometer_t), intent(in) :: self
        class(cxrs_plasma_t), intent(in) :: plasma
        real(dp), intent(in) :: x(3), e(3)
        real(dp), intent(out) :: eps, mean, var, halo, tz, flow_l, field(3)
        real(dp) :: nz, flow(3), ne, te, ti, zeff, nb(BEAM_COMPONENTS), u(3), d(3), eh(3), ev(3)
        real(dp) :: rel(3), drift, a, m0, m1, m2, sigma, s0, s1, s2, src, l, nh, q, mz, ebar
        integer :: b, j
        eps = 0.0_dp
        mean = 0.0_dp
        var = 0.0_dp
        halo = 0.0_dp
        call plasma%impurity(x, nz, tz, flow, field)
        flow_l = dot_product(flow, e)
        if (.not. (nz > 0.0_dp .and. tz > 0.0_dp)) return
        call plasma%eval(x, ne, te, ti, zeff)
        mz = self%line%mass_amu * amu
        sigma = sqrt(tz * e_charge / mz)
        s0 = 0.0_dp
        s1 = 0.0_dp
        s2 = 0.0_dp
        do b = 1, size(self%beams)
            call beam_density(self%beams(b), self%att(b), x, nb)
            call beam_frame(self%beams(b), d, eh, ev)
            do j = 1, BEAM_COMPONENTS
                if (.not. nb(j) > 0.0_dp) cycle
                u = beam_component_speed(self%beams(b), j) * d
                rel = u - flow
                drift = norm2(rel)
                a = dot_product(rel, e) / max(drift, tiny(1.0_dp))
                call cx_velocity_moments(self%line%rate, self%xh, self%wh, self%xl, self%wl, sigma, &
                    drift, a, m0, m1, m2)
                s0 = s0 + nb(j) * m0
                s1 = s1 + nb(j) * m1
                s2 = s2 + nb(j) * m2
            end do
            l = dot_product(x - self%beams(b)%origin, d)
            if (l < 0.0_dp) cycle
            src = 0.0_dp
            do j = 1, BEAM_COMPONENTS
                if (self%beams(b)%power_fraction(j) <= 0.0_dp) cycle
                src = src + beam_component_flux(self%beams(b), j) &
                    * exp(-beam_optical_depth(self%beams(b), self%att(b), l, j)) &
                    * hydrogen_charge_exchange(self%beams(b)%energy / j / self%beams(b)%isotope * 1.0e-3_dp)
            end do
            call halo_density(self%beams(b), x, src, ne, te, ti, zeff, nh)
            ebar = 1.5_dp * ti * (1.0_dp / self%beams(b)%mass_amu + 1.0_dp / self%line%mass_amu)
            q = cx_rate(self%line%rate, ebar)
            halo = halo + nz * nh * q / (4.0_dp * pi)
            s0 = s0 + nh * q
            s2 = s2 + nh * q * sigma * sigma
        end do
        if (.not. s0 > 0.0_dp) return
        eps = nz * s0 / (4.0_dp * pi)
        mean = flow_l + s1 / s0
        var = s2 / s0 - (s1 / s0)**2
    end subroutine local_emission

    subroutine sight_points(chord, x, w)
        type(chord_t), intent(in) :: chord
        real(dp), allocatable, intent(out) :: x(:, :), w(:)
        real(dp), allocatable :: xg(:), wg(:)
        real(dp) :: d(3), h
        integer :: i, k, m
        call gauss_legendre01(chord%ngauss, xg, wg)
        d = chord%direction / norm2(chord%direction)
        h = chord%length / chord%npanel
        allocate(x(3, chord%npanel * chord%ngauss), w(chord%npanel * chord%ngauss))
        m = 0
        do i = 1, chord%npanel
            do k = 1, chord%ngauss
                m = m + 1
                x(:, m) = chord%origin + (i - 1 + xg(k)) * h * d
                w(m) = h * wg(k)
            end do
        end do
    end subroutine sight_points

    function spectrometer_moments(self, plasma, channel) result(out)
        !! Line moments of one channel.
        class(cxrs_spectrometer_t), intent(in) :: self
        class(cxrs_plasma_t), intent(in) :: plasma
        type(cxrs_channel_t), intent(in) :: channel
        type(cxrs_moments_t) :: out
        real(dp), allocatable :: x(:, :), w(:)
        real(dp) :: e(3), eps, mean, var, halo, tz, fl, field(3), s0, s1, s2, sz, sh, st, sv, sx(3)
        real(dp) :: mz, lam_var, doppler
        integer :: i
        call sight_points(channel%chord, x, w)
        e = channel%chord%direction / norm2(channel%chord%direction)
        s0 = 0.0_dp
        s1 = 0.0_dp
        s2 = 0.0_dp
        sz = 0.0_dp
        sh = 0.0_dp
        st = 0.0_dp
        sv = 0.0_dp
        sx = 0.0_dp
        do i = 1, size(w)
            call local_emission(self, plasma, x(:, i), e, eps, mean, var, halo, tz, fl, field)
            if (.not. eps > 0.0_dp) cycle
            s0 = s0 + w(i) * eps
            s1 = s1 + w(i) * eps * mean
            s2 = s2 + w(i) * eps * (var + mean * mean)
            sh = sh + w(i) * halo
            st = st + w(i) * eps * tz
            sv = sv + w(i) * eps * fl
            sx = sx + w(i) * eps * x(:, i)
            if (self%line%zeeman) sz = sz + w(i) * eps * zeeman_variance(self%line, field, e)
        end do
        if (.not. s0 > 0.0_dp) return
        mz = self%line%mass_amu * amu
        out%brightness = s0
        out%halo_fraction = sh / s0
        out%velocity = s1 / s0
        doppler = s2 / s0 - out%velocity**2
        out%temperature_doppler = mz * doppler / e_charge
        lam_var = self%line%wavelength**2 * doppler / c_light**2 + sum(self%pattern_weight &
            * self%pattern_shift**2) + sz / s0 + (channel%instrument_fwhm / sqrt(8.0_dp * log(2.0_dp)))**2
        out%temperature = mz * c_light**2 * lam_var / self%line%wavelength**2 / e_charge
        out%velocity_local = sv / s0
        out%temperature_local = st / s0
        out%location = sx / s0
    end function spectrometer_moments

    function spectrometer_spectrum(self, plasma, channel, lambda) result(spec)
        !! Spectral radiance [photons m^-2 s^-1 sr^-1 m^-1] at wavelengths lambda [m].
        class(cxrs_spectrometer_t), intent(in) :: self
        class(cxrs_plasma_t), intent(in) :: plasma
        type(cxrs_channel_t), intent(in) :: channel
        real(dp), intent(in) :: lambda(:)
        real(dp) :: spec(size(lambda))
        real(dp), allocatable :: x(:, :), w(:)
        real(dp) :: e(3), eps, mean, var, halo, tz, fl, field(3), lc, s2, dz, inst, l0
        real(dp) :: zs(3), zw(3), c2, b
        integer :: i, k, m
        call sight_points(channel%chord, x, w)
        e = channel%chord%direction / norm2(channel%chord%direction)
        l0 = self%line%wavelength
        inst = (channel%instrument_fwhm / sqrt(8.0_dp * log(2.0_dp)))**2
        spec = 0.0_dp
        do i = 1, size(w)
            call local_emission(self, plasma, x(:, i), e, eps, mean, var, halo, tz, fl, field)
            if (.not. eps > 0.0_dp) cycle
            lc = l0 * (1.0_dp + mean / c_light)
            s2 = l0**2 * var / c_light**2 + inst
            zs = 0.0_dp
            zw = [1.0_dp, 0.0_dp, 0.0_dp]
            b = norm2(field)
            if (self%line%zeeman .and. b > 0.0_dp) then
                c2 = (dot_product(field, e) / b)**2
                dz = l0**2 * mu_bohr * b / (h_planck * c_light)
                zs = [0.0_dp, dz, -dz]
                zw = [0.5_dp * (1.0_dp - c2), 0.25_dp * (1.0_dp + c2), 0.25_dp * (1.0_dp + c2)]
            end if
            do k = 1, size(self%pattern_shift)
                do m = 1, 3
                    if (zw(m) <= 0.0_dp) cycle
                    spec = spec + w(i) * eps * self%pattern_weight(k) * zw(m) &
                        * exp(-0.5_dp * (lambda - lc - self%pattern_shift(k) - zs(m))**2 / s2) &
                        / sqrt(2.0_dp * pi * s2)
                end do
            end do
        end do
    end function spectrometer_spectrum
end module tiago_cxrs_spectrum
