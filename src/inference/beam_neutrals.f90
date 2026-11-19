module tiago_beam_neutrals
    !! Neutral heating-beam density with energy components, Gaussian
    !! transverse profile, divergence, electron-density-dependent stopping and
    !! a thermal halo. Shared by the charge-exchange, FIDA and neutron models.
    !!
    !! Component j = 1, 2, 3 carries the energy E / j per particle. Its density
    !! at a point x with axial distance l >= 0 from the beam origin and
    !! transverse offsets (a, b) is
    !!
    !!   n_j = Gamma_j / (v_j pi w_h w_v) exp(-a^2 / w_h^2 - b^2 / w_v^2) exp(-tau_j(l)),
    !!   Gamma_j = P f_j j / E,   v_j = sqrt(2 E / (j m)),
    !!   w_k(l)^2 = w_k^2 (1 - l / F_k)^2 + (l alpha_k)^2,
    !!   tau_j(l) = int_0^l ne sigma_s(E / (j A), ne, Te, Zeff) ds   along the axis,
    !!
    !! with power P, power fractions f_j, 1/e half-widths w_k, focal lengths F_k
    !! and 1/e divergences alpha_k. sigma_s is the fit of Suzuki et al., Plasma
    !! Phys. Control. Fusion 40 (1998) 2097, eqs. (24), (26), (28), tables 2a,
    !! 2c, 3a, 3c (beam isotope column at B = 5 T; table b for He and B).
    !! The attenuation of the axis applies to the whole cross-section.
    !!
    !! Halo: beam neutrals charge-exchange with thermal ions at the rate
    !! S = sum_j n_j n_i sigma_cx(E_j / A) v_j, Janev, Reiter and Samm (2003)
    !! eq. (44). Halo atoms are ionized by electrons at nu = ne <sigma v>_e(Te)
    !! (Maxwellian average of eq. (14)) and diffuse between charge exchanges
    !! with ions at nu_cx = n_i sigma_cx(3 Ti / A) sqrt(6 Ti / m). Over their
    !! lifetime they spread by L^2 = (Ti / m) / ((nu_cx + nu) nu) per direction,
    !! so the halo density is the beam profile widened to W_k^2 = w_k^2 + 4 L^2
    !! with the same number per unit length, S / nu, evaluated with the plasma
    !! at the point. Proton-impact ionization and charge exchange of halo atoms
    !! with impurities are neglected.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_sight_geometry, only: gauss_legendre01
    use tiago_beam_stopping_suzuki, only: beam_stopping_suzuki
    use tiago_hydrogen_1s_cross_sections, only: hydrogen_1s_cross_sections
    implicit none
    private

    integer, parameter, public :: BEAM_COMPONENTS = 3
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: e_charge = 1.602176634e-19_dp
    real(dp), parameter :: amu = 1.66053906660e-27_dp
    real(dp), parameter :: m_e = 9.1093837015e-31_dp

    type, public :: neutral_beam_t
        real(dp) :: origin(3) = 0.0_dp !! Cartesian [m], outside the plasma
        real(dp) :: direction(3) = [1.0_dp, 0.0_dp, 0.0_dp] !! normalized on use
        real(dp) :: energy = 60.0e3_dp !! full energy per particle [eV]
        integer :: isotope = 2 !! beam mass number A: 1, 2 or 3
        real(dp) :: mass_amu = 2.01410177812_dp !! beam atom mass [u]
        real(dp) :: power = 2.5e6_dp !! neutral power through the origin [W]
        real(dp) :: power_fraction(BEAM_COMPONENTS) = [1.0_dp, 0.0_dp, 0.0_dp]
        real(dp) :: width(2) = [0.1_dp, 0.1_dp] !! 1/e half-widths at the origin [m]
        real(dp) :: divergence(2) = 0.0_dp !! 1/e half-angles [rad]
        real(dp) :: focus(2) = huge(1.0_dp) !! focal distances from the origin [m]
        real(dp) :: length = 4.0_dp !! axis length of the attenuation grid [m]
        integer :: nodes = 401 !! axis nodes of the attenuation grid
        integer :: impurity_z = 7 !! impurity of the stopping fit: 2, 5, 6, 7 or 8
    end type neutral_beam_t

    type, public :: beam_attenuation_t
        !! Axis nodes s_k = (k - 1) h with plasma values and the cumulative
        !! optical depth tau(k, j) (trapezoid rule).
        real(dp) :: h = 0.0_dp
        real(dp), allocatable :: tau(:, :) !! (nodes, component)
        real(dp), allocatable :: rate(:, :) !! ne sigma_s at the nodes [m^-1]
        real(dp), allocatable :: dtau_dne(:, :), dtau_dte(:, :), dtau_dzeff(:, :)
        !! partials of ne sigma_s at node k, (nodes, component) [m^-1 per unit]
    end type beam_attenuation_t

    type, abstract, public :: beam_plasma_t
        !! Plasma seen by the beam: ne [m^-3], Te, Ti [eV] and Zeff at a
        !! Cartesian point; all zero outside the plasma.
    contains
        procedure(plasma_eval), deferred :: eval
    end type beam_plasma_t

    abstract interface
        subroutine plasma_eval(self, x, ne, te, ti, zeff)
            import :: beam_plasma_t, dp
            class(beam_plasma_t), intent(in) :: self
            real(dp), intent(in) :: x(3)
            real(dp), intent(out) :: ne, te, ti, zeff
        end subroutine plasma_eval
    end interface

    public :: beam_check, beam_axis_nodes, beam_frame, beam_attenuation, beam_density, &
        beam_density_jvp, beam_density_vjp, beam_component_speed, beam_component_flux, &
        beam_densities, halo_density, stopping_cross_section, electron_ionization_rate, &
        hydrogen_charge_exchange, beam_optical_depth

    ! Suzuki et al. (1998) tables 2a, 3a: columns H, D, T at B = 5 T.
    real(dp), parameter :: a_high(10, 3) = reshape([ &
        1.27e1_dp, 1.25_dp, 4.52e-1_dp, 1.05e-2_dp, 5.47e-1_dp, -1.02e-1_dp, 3.60e-1_dp, &
        -2.98e-2_dp, -9.59e-2_dp, 4.21e-3_dp, &
        1.41e1_dp, 1.11_dp, 4.08e-1_dp, 1.05e-2_dp, 5.47e-1_dp, -4.03e-2_dp, 3.45e-1_dp, &
        -2.88e-2_dp, -9.71e-2_dp, 4.74e-3_dp, &
        1.27e1_dp, 1.26_dp, 4.49e-1_dp, 1.05e-2_dp, 5.47e-1_dp, -5.77e-3_dp, 3.36e-1_dp, &
        -2.82e-2_dp, -9.74e-2_dp, 4.87e-3_dp], [10, 3])
    real(dp), parameter :: a_low(10, 3) = reshape([ &
        -5.29e1_dp, -1.36_dp, 7.19e-2_dp, 1.37e-2_dp, 4.54e-1_dp, 4.03e-1_dp, -2.20e-1_dp, &
        6.66e-2_dp, -6.77e-2_dp, -1.48e-3_dp, &
        -6.79e1_dp, -1.22_dp, 8.14e-2_dp, 1.39e-2_dp, 4.54e-1_dp, 4.65e-1_dp, -2.73e-1_dp, &
        7.51e-2_dp, -6.30e-2_dp, -5.08e-4_dp, &
        -7.42e1_dp, -1.18_dp, 8.43e-2_dp, 1.39e-2_dp, 4.53e-1_dp, 4.91e-1_dp, -2.94e-1_dp, &
        7.88e-2_dp, -6.12e-2_dp, -1.85e-4_dp], [10, 3])
    ! B_ijk in the order 111, 112, 121, 122, 211, ..., 322 for Z = 2, 5, 6, 7, 8:
    ! tables 2b/3b for He and B, 2c/3c (1 < Zeff < 5) for C, N and O.
    integer, parameter :: impurity_charge(5) = [2, 5, 6, 7, 8]
    real(dp), parameter :: b_high(12, 5) = reshape([ &
        2.31e-1_dp, 3.43e-1_dp, -1.85e-1_dp, -1.62e-2_dp, 1.05e-1_dp, -7.03e-2_dp, &
        5.31e-2_dp, 3.42e-3_dp, -8.38e-3_dp, 4.15e-3_dp, -3.35e-3_dp, -2.21e-4_dp, &
        -7.32e-1_dp, 1.83e-2_dp, -1.55e-1_dp, -1.72e-2_dp, 3.21e-1_dp, 9.46e-3_dp, &
        3.97e-2_dp, 4.20e-3_dp, -2.04e-2_dp, -6.19e-4_dp, -2.24e-3_dp, -2.54e-4_dp, &
        -1.01_dp, -8.65e-3_dp, -1.24e-1_dp, -1.45e-2_dp, 3.91e-1_dp, 1.61e-2_dp, &
        2.98e-2_dp, 3.32e-3_dp, -2.48e-2_dp, -1.04e-3_dp, -1.52e-3_dp, -1.89e-4_dp, &
        -1.02_dp, -1.39e-2_dp, -9.79e-2_dp, -1.17e-2_dp, 3.75e-1_dp, 1.56e-2_dp, &
        2.24e-2_dp, 2.54e-3_dp, -2.26e-2_dp, -8.89e-4_dp, -1.04e-3_dp, -1.39e-4_dp, &
        -1.02_dp, -1.48e-2_dp, -6.74e-2_dp, -9.17e-3_dp, 3.59e-1_dp, 1.43e-2_dp, &
        1.39e-2_dp, 1.84e-3_dp, -2.09e-2_dp, -7.32e-4_dp, -5.02e-4_dp, -9.49e-5_dp], [12, 5])
    real(dp), parameter :: b_low(12, 5) = reshape([ &
        -7.92e-1_dp, 4.20e-2_dp, 5.30e-2_dp, -1.39e-2_dp, 3.01e-1_dp, -2.64e-2_dp, &
        -2.99e-2_dp, 6.07e-3_dp, 2.72e-4_dp, 6.11e-3_dp, 3.47e-3_dp, -9.19e-4_dp, &
        1.22e-1_dp, 5.27e-2_dp, -4.30e-4_dp, -3.18e-3_dp, -1.51e-1_dp, -3.64e-2_dp, &
        3.43e-3_dp, 1.51e-3_dp, 4.20e-2_dp, 6.92e-3_dp, -1.41e-3_dp, -2.90e-4_dp, &
        1.61e-1_dp, 5.98e-2_dp, -3.36e-3_dp, -4.26e-3_dp, -1.57e-1_dp, -3.96e-2_dp, &
        4.60e-3_dp, 2.19e-3_dp, 3.91e-2_dp, 7.11e-3_dp, -1.44e-3_dp, -3.85e-4_dp, &
        1.39e-1_dp, 6.06e-2_dp, -3.06e-3_dp, -4.55e-3_dp, -1.33e-1_dp, -3.94e-2_dp, &
        3.99e-3_dp, 2.36e-3_dp, 3.35e-2_dp, 6.90e-3_dp, -1.24e-3_dp, -4.05e-4_dp, &
        1.11e-1_dp, 5.41e-2_dp, -3.46e-4_dp, -3.68e-3_dp, -1.08e-1_dp, -3.47e-2_dp, &
        1.93e-3_dp, 1.81e-3_dp, 2.80e-2_dp, 6.04e-3_dp, -8.41e-4_dp, -3.17e-4_dp], [12, 5])

contains

    pure integer function beam_check(beam)
        !! 0 if valid; 1 bad direction, 2 bad isotope or impurity, 3 non-positive
        !! energy, power, width or grid, 4 a component energy below the 10 keV/amu
        !! limit of the stopping fit (9 keV/amu tolerated) or above 10 MeV/amu.
        type(neutral_beam_t), intent(in) :: beam
        real(dp) :: e_amu
        integer :: j
        beam_check = 1
        if (.not. norm2(beam%direction) > 0.0_dp) return
        beam_check = 2
        if (beam%isotope < 1 .or. beam%isotope > 3) return
        if (.not. any(impurity_charge == beam%impurity_z)) return
        beam_check = 3
        if (.not. (beam%energy > 0.0_dp .and. beam%power >= 0.0_dp .and. beam%mass_amu > 0.0_dp)) return
        if (.not. all(beam%width > 0.0_dp) .or. any(beam%divergence < 0.0_dp)) return
        if (.not. beam%length > 0.0_dp .or. beam%nodes < 2) return
        if (any(beam%power_fraction < 0.0_dp)) return
        beam_check = 4
        do j = 1, BEAM_COMPONENTS
            if (beam%power_fraction(j) <= 0.0_dp) cycle
            e_amu = beam%energy / j / beam%isotope * 1.0e-3_dp
            if (e_amu < 9.0_dp .or. e_amu > 1.0e4_dp) return
        end do
        beam_check = 0
    end function beam_check

    pure subroutine beam_frame(beam, d, eh, ev)
        !! Axis d, horizontal eh = (z x d) / |z x d| and vertical ev = d x eh.
        type(neutral_beam_t), intent(in) :: beam
        real(dp), intent(out) :: d(3), eh(3), ev(3)
        d = beam%direction / norm2(beam%direction)
        eh = [-d(2), d(1), 0.0_dp]
        if (norm2(eh) < 1.0e-12_dp) eh = [0.0_dp, -d(3), d(2)]
        eh = eh / norm2(eh)
        ev = [d(2) * eh(3) - d(3) * eh(2), d(3) * eh(1) - d(1) * eh(3), d(1) * eh(2) - d(2) * eh(1)]
    end subroutine beam_frame

    pure function beam_axis_nodes(beam) result(x)
        !! Cartesian nodes of the attenuation grid [m], (3, nodes).
        type(neutral_beam_t), intent(in) :: beam
        real(dp) :: x(3, beam%nodes), d(3), eh(3), ev(3), h
        integer :: k
        call beam_frame(beam, d, eh, ev)
        h = beam%length / (beam%nodes - 1)
        do k = 1, beam%nodes
            x(:, k) = beam%origin + (k - 1) * h * d
        end do
    end function beam_axis_nodes

    pure real(dp) function beam_component_speed(beam, j)
        !! Speed of component j [m/s].
        type(neutral_beam_t), intent(in) :: beam
        integer, intent(in) :: j
        beam_component_speed = sqrt(2.0_dp * beam%energy * e_charge / (j * beam%mass_amu * amu))
    end function beam_component_speed

    pure real(dp) function beam_component_flux(beam, j)
        !! Atoms per second of component j through the origin.
        type(neutral_beam_t), intent(in) :: beam
        integer, intent(in) :: j
        beam_component_flux = beam%power * beam%power_fraction(j) * j / (beam%energy * e_charge)
    end function beam_component_flux

    pure subroutine stopping_cross_section(beam, j, ne, te, zeff, sigma, dsigma)
        !! Suzuki stopping cross-section [m^2] of component j and its partials
        !! in (ne, Te, Zeff). Zero for ne <= 0 or Te <= 0 (outside the plasma).
        type(neutral_beam_t), intent(in) :: beam
        integer, intent(in) :: j
        real(dp), intent(in) :: ne, te, zeff
        real(dp), intent(out) :: sigma, dsigma(3)
        real(dp) :: e_amu, a(10), b(12)
        integer :: iz
        sigma = 0.0_dp
        dsigma = 0.0_dp
        if (.not. (ne > 0.0_dp .and. te > 0.0_dp)) return
        e_amu = beam%energy / j / beam%isotope * 1.0e-3_dp
        iz = findloc(impurity_charge, beam%impurity_z, 1)
        if (e_amu < 100.0_dp) then
            a = a_low(:, beam%isotope)
            b = b_low(:, iz)
        else
            a = a_high(:, beam%isotope)
            b = b_high(:, iz)
        end if
        call beam_stopping_suzuki(e_amu, ne, te, zeff, a(1), a(2), a(3), a(4), a(5), a(6), a(7), &
            a(8), a(9), a(10), b(1), b(2), b(3), b(4), b(5), b(6), b(7), b(8), b(9), b(10), &
            b(11), b(12), sigma, dsigma(1), dsigma(2), dsigma(3))
    end subroutine stopping_cross_section

    pure subroutine beam_attenuation(beam, ne, te, zeff, att)
        !! Cumulative optical depth of every component from plasma values at
        !! the axis nodes (beam_axis_nodes).
        type(neutral_beam_t), intent(in) :: beam
        real(dp), intent(in) :: ne(:), te(:), zeff(:)
        type(beam_attenuation_t), intent(out) :: att
        real(dp) :: sigma, ds(3)
        integer :: j, k
        att%h = beam%length / (beam%nodes - 1)
        allocate(att%tau(beam%nodes, BEAM_COMPONENTS), att%rate(beam%nodes, BEAM_COMPONENTS), &
            att%dtau_dne(beam%nodes, BEAM_COMPONENTS), &
            att%dtau_dte(beam%nodes, BEAM_COMPONENTS), att%dtau_dzeff(beam%nodes, BEAM_COMPONENTS))
        do j = 1, BEAM_COMPONENTS
            do k = 1, beam%nodes
                call stopping_cross_section(beam, j, ne(k), te(k), zeff(k), sigma, ds)
                att%rate(k, j) = ne(k) * sigma
                att%dtau_dne(k, j) = sigma + ne(k) * ds(1)
                att%dtau_dte(k, j) = ne(k) * ds(2)
                att%dtau_dzeff(k, j) = ne(k) * ds(3)
            end do
            att%tau(1, j) = 0.0_dp
            do k = 2, beam%nodes
                att%tau(k, j) = att%tau(k - 1, j) + 0.5_dp * att%h * (att%rate(k - 1, j) + att%rate(k, j))
            end do
        end do
    end subroutine beam_attenuation

    pure subroutine axial_weights(beam, att, l, c, k0)
        !! tau(l) = sum_k c(k) rate(k), k = 1 .. k0 + 1 (trapezoid up to l,
        !! linear between nodes), constant beyond the grid.
        type(neutral_beam_t), intent(in) :: beam
        type(beam_attenuation_t), intent(in) :: att
        real(dp), intent(in) :: l
        real(dp), intent(out) :: c(beam%nodes)
        integer, intent(out) :: k0
        real(dp) :: t
        c = 0.0_dp
        k0 = min(int(max(l, 0.0_dp) / att%h) + 1, beam%nodes - 1)
        t = min(max(l, 0.0_dp) / att%h - (k0 - 1), 1.0_dp)
        c(1:k0) = att%h
        c(1) = 0.5_dp * att%h
        c(k0) = c(k0) - 0.5_dp * att%h
        ! Partial panel [s_k0, s_k0 + t h]: trapezoid of the linear interpolant.
        c(k0) = c(k0) + 0.5_dp * att%h * t * (2.0_dp - t)
        c(k0 + 1) = 0.5_dp * att%h * t * t
        if (k0 == 1) c(1) = 0.5_dp * att%h * t * (2.0_dp - t)
    end subroutine axial_weights

    pure subroutine beam_geometry(beam, x, l, profile, area)
        !! Axial distance l, Gaussian factor and pi w_h w_v at x.
        type(neutral_beam_t), intent(in) :: beam
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: l, profile, area
        real(dp) :: d(3), eh(3), ev(3), rel(3), w(2), off(2)
        call beam_frame(beam, d, eh, ev)
        rel = x - beam%origin
        l = dot_product(rel, d)
        off = [dot_product(rel, eh), dot_product(rel, ev)]
        w = sqrt((beam%width * (1.0_dp - l / beam%focus))**2 + (l * beam%divergence)**2)
        profile = exp(-sum((off / w)**2))
        area = pi * w(1) * w(2)
    end subroutine beam_geometry

    pure subroutine beam_density(beam, att, x, nb)
        !! Component densities [m^-3] at a Cartesian point; zero behind the origin.
        type(neutral_beam_t), intent(in) :: beam
        type(beam_attenuation_t), intent(in) :: att
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: nb(BEAM_COMPONENTS)
        real(dp) :: l, profile, area
        integer :: j
        nb = 0.0_dp
        call beam_geometry(beam, x, l, profile, area)
        if (l < 0.0_dp) return
        do j = 1, BEAM_COMPONENTS
            if (beam%power_fraction(j) <= 0.0_dp) cycle
            nb(j) = beam_component_flux(beam, j) / (beam_component_speed(beam, j) * area) &
                * profile * exp(-beam_optical_depth(beam, att, l, j))
        end do
    end subroutine beam_density

    pure real(dp) function beam_optical_depth(beam, att, l, j)
        !! Optical depth tau_j at axial distance l [m].
        type(neutral_beam_t), intent(in) :: beam
        type(beam_attenuation_t), intent(in) :: att
        real(dp), intent(in) :: l
        integer, intent(in) :: j
        real(dp) :: c(beam%nodes)
        integer :: k0
        call axial_weights(beam, att, l, c, k0)
        beam_optical_depth = sum(c(:k0 + 1) * att%rate(:k0 + 1, j))
    end function beam_optical_depth

    pure subroutine beam_density_jvp(beam, att, x, dne, dte, dzeff, dnb)
        !! Tangent of the component densities at x for node perturbations of
        !! (ne, Te, Zeff) along the axis.
        type(neutral_beam_t), intent(in) :: beam
        type(beam_attenuation_t), intent(in) :: att
        real(dp), intent(in) :: x(3), dne(:), dte(:), dzeff(:)
        real(dp), intent(out) :: dnb(BEAM_COMPONENTS)
        real(dp) :: nb(BEAM_COMPONENTS), l, profile, area, c(beam%nodes)
        integer :: j, k0
        call beam_density(beam, att, x, nb)
        dnb = 0.0_dp
        call beam_geometry(beam, x, l, profile, area)
        if (l < 0.0_dp) return
        call axial_weights(beam, att, l, c, k0)
        do j = 1, BEAM_COMPONENTS
            dnb(j) = -nb(j) * sum(c(:k0 + 1) * (att%dtau_dne(:k0 + 1, j) * dne(:k0 + 1) &
                + att%dtau_dte(:k0 + 1, j) * dte(:k0 + 1) + att%dtau_dzeff(:k0 + 1, j) * dzeff(:k0 + 1)))
        end do
    end subroutine beam_density_jvp

    pure subroutine beam_density_vjp(beam, att, x, nb_bar, ne_bar, te_bar, zeff_bar)
        !! Adds the cotangent of the component densities at x to the node
        !! cotangents of (ne, Te, Zeff).
        type(neutral_beam_t), intent(in) :: beam
        type(beam_attenuation_t), intent(in) :: att
        real(dp), intent(in) :: x(3), nb_bar(BEAM_COMPONENTS)
        real(dp), intent(inout) :: ne_bar(:), te_bar(:), zeff_bar(:)
        real(dp) :: nb(BEAM_COMPONENTS), l, profile, area, c(beam%nodes), g
        integer :: j, k0
        call beam_density(beam, att, x, nb)
        call beam_geometry(beam, x, l, profile, area)
        if (l < 0.0_dp) return
        call axial_weights(beam, att, l, c, k0)
        do j = 1, BEAM_COMPONENTS
            g = -nb(j) * nb_bar(j)
            ne_bar(:k0 + 1) = ne_bar(:k0 + 1) + g * c(:k0 + 1) * att%dtau_dne(:k0 + 1, j)
            te_bar(:k0 + 1) = te_bar(:k0 + 1) + g * c(:k0 + 1) * att%dtau_dte(:k0 + 1, j)
            zeff_bar(:k0 + 1) = zeff_bar(:k0 + 1) + g * c(:k0 + 1) * att%dtau_dzeff(:k0 + 1, j)
        end do
    end subroutine beam_density_vjp

    pure real(dp) function hydrogen_charge_exchange(e_amu)
        !! H+ + H(1s) total charge-exchange cross-section [m^2] at E [keV/amu].
        real(dp), intent(in) :: e_amu
        real(dp) :: sp, se
        call hydrogen_1s_cross_sections(e_amu, 1.0e3_dp, hydrogen_charge_exchange, sp, se)
    end function hydrogen_charge_exchange

    pure real(dp) function electron_ionization_rate(te)
        !! Maxwellian <sigma v> [m^3/s] of e + H(1s) ionization at Te [eV].
        real(dp), intent(in) :: te
        real(dp), parameter :: ion = 13.6_dp
        integer, parameter :: npanel = 200
        real(dp), allocatable :: xg(:), wg(:)
        real(dp) :: y, e, scx, sp, se, s, h
        integer :: i, k
        electron_ionization_rate = 0.0_dp
        if (.not. te > 0.0_dp) return
        call gauss_legendre01(4, xg, wg)
        ! <sigma v> = sqrt(8 / (pi m)) T^(-3/2) int_I^inf sigma E exp(-E/T) dE,
        ! integrated in y = ln(E / I) up to E = I + 60 T.
        h = log(1.0_dp + 60.0_dp * te / ion) / npanel
        s = 0.0_dp
        do i = 1, npanel
            do k = 1, 4
                y = (i - 1 + xg(k)) * h
                e = ion * exp(y)
                call hydrogen_1s_cross_sections(1.0_dp, e, scx, sp, se)
                s = s + h * wg(k) * max(se, 0.0_dp) * e * e * exp(-e / te)
            end do
        end do
        electron_ionization_rate = sqrt(8.0_dp / (pi * m_e)) * (te * e_charge)**(-1.5_dp) &
            * s * e_charge**2
    end function electron_ionization_rate

    pure subroutine halo_density(beam, x, nb_axis_rate, ne, te, ti, zeff, nh)
        !! Halo density [m^-3] at x from the charge-exchange source per unit
        !! length nb_axis_rate = sum_j Gamma_j exp(-tau_j) sigma_cx(E_j / A) [m^2/s]
        !! (see beam_densities) and the local plasma.
        type(neutral_beam_t), intent(in) :: beam
        real(dp), intent(in) :: x(3), nb_axis_rate, ne, te, ti, zeff
        real(dp), intent(out) :: nh
        real(dp) :: d(3), eh(3), ev(3), rel(3), l, off(2), w(2), wh(2), ni, nu, nu_cx, m, len2
        nh = 0.0_dp
        if (.not. (ne > 0.0_dp .and. te > 0.0_dp .and. ti > 0.0_dp)) return
        call beam_frame(beam, d, eh, ev)
        rel = x - beam%origin
        l = dot_product(rel, d)
        if (l < 0.0_dp) return
        off = [dot_product(rel, eh), dot_product(rel, ev)]
        w = sqrt((beam%width * (1.0_dp - l / beam%focus))**2 + (l * beam%divergence)**2)
        m = beam%mass_amu * amu
        ni = ne * (beam%impurity_z - zeff) / (beam%impurity_z - 1)
        nu = ne * electron_ionization_rate(te)
        nu_cx = ni * hydrogen_charge_exchange(3.0e-3_dp * ti / beam%isotope) &
            * sqrt(6.0_dp * ti * e_charge / m)
        len2 = ti * e_charge / m / ((nu_cx + nu) * nu)
        wh = sqrt(w**2 + 4.0_dp * len2)
        nh = ni * nb_axis_rate / nu * exp(-sum((off / wh)**2)) / (pi * wh(1) * wh(2))
    end subroutine halo_density

    subroutine beam_densities(beam, plasma, x, nb, nh, status)
        !! Component and halo densities [m^-3] at Cartesian points x(3, n) for a
        !! plasma given as a callback. status: beam_check, or 5 if the plasma
        !! is not empty at the beam origin.
        type(neutral_beam_t), intent(in) :: beam
        class(beam_plasma_t), intent(in) :: plasma
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: nb(BEAM_COMPONENTS, size(x, 2)), nh(size(x, 2))
        integer, intent(out) :: status
        type(beam_attenuation_t) :: att
        real(dp) :: axis(3, beam%nodes), ne(beam%nodes), te(beam%nodes), ti(beam%nodes)
        real(dp) :: zeff(beam%nodes), d(3), eh(3), ev(3), l, src, pne, pte, pti, pzeff
        integer :: i, j, k
        nb = 0.0_dp
        nh = 0.0_dp
        status = beam_check(beam)
        if (status /= 0) return
        axis = beam_axis_nodes(beam)
        do k = 1, beam%nodes
            call plasma%eval(axis(:, k), ne(k), te(k), ti(k), zeff(k))
        end do
        status = 5
        if (ne(1) > 0.0_dp) return
        status = 0
        call beam_attenuation(beam, ne, te, zeff, att)
        call beam_frame(beam, d, eh, ev)
        do i = 1, size(x, 2)
            call beam_density(beam, att, x(:, i), nb(:, i))
            l = dot_product(x(:, i) - beam%origin, d)
            if (l < 0.0_dp) cycle
            src = 0.0_dp
            do j = 1, BEAM_COMPONENTS
                if (beam%power_fraction(j) <= 0.0_dp) cycle
                src = src + beam_component_flux(beam, j) &
                    * exp(-beam_optical_depth(beam, att, l, j)) &
                    * hydrogen_charge_exchange(beam%energy / j / beam%isotope * 1.0e-3_dp)
            end do
            call plasma%eval(x(:, i), pne, pte, pti, pzeff)
            call halo_density(beam, x(:, i), src, pne, pte, pti, pzeff, nh(i))
        end do
    end subroutine beam_densities
end module tiago_beam_neutrals
