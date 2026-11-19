module test_beam_neutrals_plasma
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_beam_neutrals, only: beam_plasma_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)

    type, extends(beam_plasma_t) :: slab_t
        !! ne = n0 sin^2(pi (x - x1) / (x2 - x1)) for x1 < x < x2, Te, Ti, Zeff constant.
        real(dp) :: x1 = 1.0_dp, x2 = 3.0_dp, n0 = 6.0e19_dp, te = 2.0e3_dp
        real(dp) :: ti = 2.5e3_dp, zeff = 1.6_dp
    contains
        procedure :: eval => slab_eval
    end type slab_t
contains
    subroutine slab_eval(self, x, ne, te, ti, zeff)
        class(slab_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: ne, te, ti, zeff
        ne = 0.0_dp
        te = 0.0_dp
        ti = 0.0_dp
        zeff = 1.0_dp
        if (x(1) <= self%x1 .or. x(1) >= self%x2) return
        ne = self%n0 * sin(pi * (x(1) - self%x1) / (self%x2 - self%x1))**2
        te = self%te
        ti = self%ti
        zeff = self%zeff
    end subroutine slab_eval
end module test_beam_neutrals_plasma

program test_beam_neutrals
    !! Neutral beam densities (tiago_beam_neutrals).
    !! Oracles: the Suzuki stopping fit against ASCOT5's independent
    !! implementation (itpplasma/ascot5 29e093f8, src/suzuki.c, D beam into
    !! D + N7+); its partial derivatives against central differences; the
    !! optical depth of a slab against composite Gauss sums of the closed-form
    !! integrand; the unattenuated density by hand; the transverse integral of
    !! n_j v_j against the attenuated particle flux with focus and divergence;
    !! the JVP against central differences of the node values and the
    !! transposition identity of the VJP; the Maxwellian ionization rate
    !! against a trapezoid sum over electron speed; the halo per unit length
    !! against its charge-exchange source over the ionization rate.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_beam_neutrals
    use tiago_sight_geometry, only: gauss_legendre01
    use test_beam_neutrals_plasma, only: slab_t, pi
    implicit none
    real(dp), parameter :: e_charge = 1.602176634e-19_dp
    integer :: failures

    failures = 0
    call check_suzuki()
    call check_attenuation()
    call check_profile()
    call check_derivatives()
    call check_ionization_rate()
    call check_halo()
    if (failures > 0) then
        write (error_unit, '(I0,A)') failures, ' beam neutral check(s) failed'
        error stop 1
    end if
    print '(A)', 'beam neutrals: all checks passed'

contains

    subroutine expect(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            failures = failures + 1
            write (error_unit, '(A)') 'FAIL: '//what
        end if
    end subroutine expect

    function d_beam(e_amu, z) result(beam)
        real(dp), intent(in) :: e_amu
        integer, intent(in) :: z
        type(neutral_beam_t) :: beam
        beam%isotope = 2
        beam%energy = 2.0e3_dp * e_amu
        beam%impurity_z = z
    end function d_beam

    subroutine check_suzuki()
        ! (E [keV/amu], ne [m^-3], Te [eV], Zeff, sigma [m^2]) from ASCOT5.
        real(dp), parameter :: ref(5, 7) = reshape([ &
            9.58_dp, 4.0e19_dp, 1500.0_dp, 1.0_dp, 9.8326043781e-20_dp, &
            30.0_dp, 4.0e19_dp, 1500.0_dp, 1.6_dp, 5.6227716794e-20_dp, &
            15.0_dp, 6.0e19_dp, 3000.0_dp, 1.6_dp, 7.9345841973e-20_dp, &
            10.0_dp, 2.0e19_dp, 800.0_dp, 2.2_dp, 1.0024612504e-19_dp, &
            20.0_dp, 1.0e20_dp, 5000.0_dp, 1.3_dp, 6.6844819620e-20_dp, &
            150.0_dp, 5.0e19_dp, 10000.0_dp, 1.8_dp, 1.9405816174e-20_dp, &
            500.0_dp, 1.0e20_dp, 20000.0_dp, 2.5_dp, 1.1984533447e-20_dp], [5, 7])
        real(dp) :: sigma, ds(3), sp, sm, dm(3), h
        integer :: i, k
        do i = 1, size(ref, 2)
            call stopping_cross_section(d_beam(ref(1, i), 7), 1, ref(2, i), ref(3, i), ref(4, i), &
                sigma, ds)
            call expect(abs(sigma / ref(5, i) - 1.0_dp) < 1.0e-9_dp, 'Suzuki against ASCOT5')
            do k = 1, 3
                block
                    real(dp) :: p(3)
                    p = ref(2:4, i)
                    h = 1.0e-5_dp * p(k)
                    p(k) = ref(1 + k, i) + h
                    call stopping_cross_section(d_beam(ref(1, i), 7), 1, p(1), p(2), p(3), sp, dm)
                    p(k) = ref(1 + k, i) - h
                    call stopping_cross_section(d_beam(ref(1, i), 7), 1, p(1), p(2), p(3), sm, dm)
                    call expect(abs((sp - sm) / (2 * h) - ds(k)) <= 1.0e-6_dp * abs(ds(k)) &
                        + 1.0e-12_dp * sigma / ref(1 + k, i), 'Suzuki partial derivative')
                end block
            end do
        end do
    end subroutine check_suzuki

    function test_beam() result(beam)
        type(neutral_beam_t) :: beam
        beam%origin = [0.0_dp, 0.0_dp, 0.0_dp]
        beam%direction = [1.0_dp, 0.0_dp, 0.0_dp]
        beam%energy = 60.0e3_dp
        beam%power = 2.0e6_dp
        beam%power_fraction = [0.6_dp, 0.3_dp, 0.1_dp]
        beam%width = [0.12_dp, 0.18_dp]
        beam%divergence = [0.01_dp, 0.02_dp]
        beam%focus = [6.0_dp, 8.0_dp]
        beam%length = 4.0_dp
        beam%nodes = 2001
    end function test_beam

    subroutine check_attenuation()
        type(neutral_beam_t) :: beam
        type(slab_t) :: slab
        real(dp) :: nb(BEAM_COMPONENTS, 1), nh(1), x(3, 1), tau, sigma, ds(3), s, ne
        real(dp), allocatable :: xg(:), wg(:)
        integer :: status, i, k, j
        beam = test_beam()
        x(:, 1) = [2.4_dp, 0.0_dp, 0.0_dp]
        call beam_densities(beam, slab, x, nb, nh, status)
        call expect(status == 0, 'beam_densities status')
        call gauss_legendre01(8, xg, wg)
        do j = 1, BEAM_COMPONENTS
            tau = 0.0_dp
            do i = 1, 200
                do k = 1, 8
                    s = slab%x1 + (i - 1 + xg(k)) * (x(1, 1) - slab%x1) / 200
                    ne = slab%n0 * sin(pi * (s - slab%x1) / (slab%x2 - slab%x1))**2
                    call stopping_cross_section(beam, j, ne, slab%te, slab%zeff, sigma, ds)
                    tau = tau + wg(k) * (x(1, 1) - slab%x1) / 200 * ne * sigma
                end do
            end do
            ! The axis passes through x, so the profile factor is one there.
            call expect(abs(nb(j, 1) / (beam_component_flux(beam, j) &
                / (beam_component_speed(beam, j) * area_at(beam, x(1, 1)))) - exp(-tau)) &
                < 1.0e-4_dp * exp(-tau), 'slab optical depth')
        end do
        call expect(tau > 0.1_dp, 'slab optical depth is not trivial')
        ! Unattenuated density before the plasma: Gamma_1 / (v_1 pi w_h w_v).
        x(:, 1) = [0.5_dp, 0.0_dp, 0.0_dp]
        call beam_densities(beam, slab, x, nb, nh, status)
        call expect(abs(nb(1, 1) / (0.6_dp * 2.0e6_dp / (60.0e3_dp * e_charge) &
            / (sqrt(2 * 60.0e3_dp * e_charge / (2.01410177812_dp * 1.66053906660e-27_dp)) &
            * area_at(beam, 0.5_dp))) - 1.0_dp) < 1.0e-12_dp, 'unattenuated density')
    end subroutine check_attenuation

    real(dp) function area_at(beam, l)
        type(neutral_beam_t), intent(in) :: beam
        real(dp), intent(in) :: l
        real(dp) :: w(2)
        w = sqrt((beam%width * (1 - l / beam%focus))**2 + (l * beam%divergence)**2)
        area_at = pi * w(1) * w(2)
    end function area_at

    subroutine check_profile()
        type(neutral_beam_t) :: beam
        type(slab_t) :: slab
        integer, parameter :: m = 81
        real(dp) :: x(3, m * m), nb(BEAM_COMPONENTS, m * m), nh(m * m), flux, tau
        real(dp) :: ne(4001), te(4001), ti(4001), zeff(4001), axis(3, 4001), h
        type(beam_attenuation_t) :: att
        integer :: i, k, status, j
        beam = test_beam()
        beam%direction = [1.0_dp, 0.3_dp, -0.2_dp]
        beam%nodes = 4001
        h = 2.4_dp / (m - 1)
        do i = 1, m
            do k = 1, m
                ! Plane x = 2.6 m with the transverse offsets on a 2.4 m square.
                x(:, (i - 1) * m + k) = [2.6_dp, 0.78_dp - 1.2_dp + (i - 1) * h, &
                    -0.52_dp - 1.2_dp + (k - 1) * h]
            end do
        end do
        call beam_densities(beam, slab, x, nb, nh, status)
        axis = beam_axis_nodes(beam)
        do k = 1, beam%nodes
            call slab%eval(axis(:, k), ne(k), te(k), ti(k), zeff(k))
        end do
        call beam_attenuation(beam, ne, te, zeff, att)
        do j = 1, BEAM_COMPONENTS
            ! Flux through the plane x = 2.6: n v d_x dA. Exp(-tau) varies across
            ! the plane because points have different axial distances, so
            ! compare with the attenuation at each point.
            flux = 0.0_dp
            do i = 1, m * m
                tau = beam_optical_depth(beam, att, dot_product(x(:, i), beam%direction) &
                    / norm2(beam%direction), j)
                flux = flux + nb(j, i) * exp(tau) * beam_component_speed(beam, j) * h * h &
                    * beam%direction(1) / norm2(beam%direction)
            end do
            ! Off-axis points have l different from the axis point, so the
            ! Gaussian widths vary over the plane; the tolerance covers that.
            call expect(abs(flux / beam_component_flux(beam, j) - 1.0_dp) < 2.0e-3_dp, &
                'transverse integral of n v')
        end do
    end subroutine check_profile

    subroutine check_derivatives()
        type(neutral_beam_t) :: beam
        type(slab_t) :: slab
        type(beam_attenuation_t) :: att, ap, am
        integer, parameter :: n = 301
        real(dp) :: axis(3, n), ne(n), te(n), ti(n), zeff(n), dne(n), dte(n), dz(n)
        real(dp) :: x(3), nbp(3), nbm(3), dnb(3), nb_bar(3), ne_bar(n), te_bar(n), z_bar(n), eps
        integer :: k
        beam = test_beam()
        beam%nodes = n
        axis = beam_axis_nodes(beam)
        do k = 1, n
            call slab%eval(axis(:, k), ne(k), te(k), ti(k), zeff(k))
            te(k) = te(k) * (1 + 0.1_dp * sin(1.0_dp * k))
            dne(k) = 1.0e19_dp * cos(0.3_dp * k)
            dte(k) = 300.0_dp * sin(0.7_dp * k)
            dz(k) = 0.2_dp * cos(1.1_dp * k)
            if (.not. ne(k) > 0.0_dp) then
                dne(k) = 0.0_dp
                dte(k) = 0.0_dp
                dz(k) = 0.0_dp
            end if
        end do
        x = [2.37_dp, 0.05_dp, -0.08_dp]
        eps = 1.0e-4_dp
        call beam_attenuation(beam, ne, te, zeff, att)
        call beam_attenuation(beam, ne + eps * dne, te + eps * dte, zeff + eps * dz, ap)
        call beam_attenuation(beam, ne - eps * dne, te - eps * dte, zeff - eps * dz, am)
        call beam_density(beam, ap, x, nbp)
        call beam_density(beam, am, x, nbm)
        call beam_density_jvp(beam, att, x, dne, dte, dz, dnb)
        call expect(all(abs((nbp - nbm) / (2 * eps) - dnb) <= 1.0e-6_dp * abs(dnb) + 1.0e-6_dp), &
            'beam density JVP against central differences')
        call expect(any(abs(dnb) > 0.0_dp), 'beam density JVP is not zero')
        nb_bar = [0.3_dp, -1.2_dp, 2.0_dp]
        ne_bar = 0.0_dp
        te_bar = 0.0_dp
        z_bar = 0.0_dp
        call beam_density_vjp(beam, att, x, nb_bar, ne_bar, te_bar, z_bar)
        call expect(abs(dot_product(nb_bar, dnb) - (dot_product(ne_bar, dne) &
            + dot_product(te_bar, dte) + dot_product(z_bar, dz))) <= 1.0e-12_dp &
            * sum(abs(nb_bar * dnb)), 'beam density VJP transposition')
    end subroutine check_derivatives

    subroutine check_ionization_rate()
        real(dp), parameter :: m_e = 9.1093837015e-31_dp
        real(dp) :: te, v, dv, e, s, sp, scx, sig, ref, f
        integer :: i, it
        real(dp), parameter :: temps(3) = [20.0_dp, 300.0_dp, 4000.0_dp]
        do it = 1, 3
            te = temps(it)
            ! <sigma v> = int sigma(v) v f(v) 4 pi v^2 dv over a Maxwellian.
            dv = sqrt(2 * 80 * te * e_charge / m_e) / 400000
            ref = 0.0_dp
            do i = 1, 400000
                v = (i - 0.5_dp) * dv
                e = 0.5_dp * m_e * v * v / e_charge
                if (e <= 13.6_dp) cycle
                block
                    use tiago_hydrogen_1s_cross_sections, only: hydrogen_1s_cross_sections
                    call hydrogen_1s_cross_sections(1.0_dp, e, scx, sp, sig)
                end block
                f = (m_e / (2 * pi * te * e_charge))**1.5_dp * exp(-e / te)
                ref = ref + max(sig, 0.0_dp) * v * f * 4 * pi * v * v * dv
            end do
            s = electron_ionization_rate(te)
            call expect(abs(s / ref - 1.0_dp) < 1.0e-6_dp, 'Maxwellian ionization rate')
        end do
    end subroutine check_ionization_rate

    subroutine check_halo()
        type(neutral_beam_t) :: beam
        integer, parameter :: m = 121
        real(dp) :: h, s, nh, ne, te, ti, zeff, src, total
        integer :: i, k
        beam = test_beam()
        ne = 5.0e19_dp
        te = 2.0e3_dp
        ti = 2.5e3_dp
        zeff = 1.6_dp
        src = 3.0e-3_dp
        h = 3.0_dp / (m - 1)
        total = 0.0_dp
        do i = 1, m
            do k = 1, m
                call halo_density(beam, [2.0_dp, -1.5_dp + (i - 1) * h, -1.5_dp + (k - 1) * h], &
                    src, ne, te, ti, zeff, nh)
                total = total + nh * h * h
            end do
        end do
        s = ne * (7 - zeff) / 6 * src / (ne * electron_ionization_rate(te))
        call expect(abs(total / s - 1.0_dp) < 1.0e-6_dp, 'halo atoms per unit length')
    end subroutine check_halo
end program test_beam_neutrals
