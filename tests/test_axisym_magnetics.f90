program test_axisym_magnetics
    !! Axisymmetric sensors (docs/axisymmetric-free-boundary.md), checked by
    !! 1. straight-segment Biot-Savart on polygons of 8192 and 16384 sides of
    !!    the same circle, Richardson-extrapolated in the number of sides;
    !! 2. a Solov'ev plasma: the volume integral of j_phi against the virtual
    !!    casing sheet K_phi = -dpsi/dn / (mu0 R) on its boundary psi = 0,
    !!    and the total current against the Rogowski signal;
    !! 3. coil-group responses, rectangle turns and the pickup length average.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_axisym_magnetics, only: axisym_sensors_t, axisym_flux_field, axisym_rectangle
    use tiago_coil_kernels, only: coil_set_t, build_coil_set
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 4.0e-7_dp * pi
    ! Solov'ev psi = c1 (R^2 - r0^2)^2 / 8 + c2 Z^2 / 2 - c3, so that
    ! Delta* psi = c1 R^2 + c2 = -mu0 R j_phi.
    real(dp), parameter :: c1 = 1.0_dp, c2 = 0.8_dp, c3 = 0.0595_dp, r0 = 1.0_dp
    integer :: failures

    failures = 0
    call check_polygon()
    call check_solovev()
    call check_coils()
    if (failures > 0) then
        write(error_unit, '(I0,A)') failures, ' axisym magnetics checks failed'
        error stop 1
    end if
    print '(A)', 'axisym magnetics: all checks passed'

contains

    subroutine check_polygon()
        real(dp), parameter :: a = 1.2_dp, h = 0.3_dp, cur = 2.5e4_dp
        real(dp), parameter :: pts(2, 4) = reshape([1.8_dp, 0.9_dp, 0.4_dp, -0.6_dp, &
            1.25_dp, 0.31_dp, 2.6_dp, 1.7_dp], [2, 4])
        real(dp) :: psi, br, bz, e(3, 2), ref(3)
        integer :: i, level
        do i = 1, size(pts, 2)
            call axisym_flux_field(pts(1, i), pts(2, i), [a], [h], [cur], psi, br, bz)
            do level = 1, 2
                e(:, level) = polygon(a, h, cur, pts(1, i), pts(2, i), 8192 * level)
            end do
            ! The inscribed polygon errs by O(n^-2).
            ref = (4 * e(:, 2) - e(:, 1)) / 3
            call close_to(psi, ref(1), 1.0e-10_dp, 'polygon psi')
            call close_to(br, ref(2), 1.0e-10_dp, 'polygon B_R')
            call close_to(bz, ref(3), 1.0e-10_dp, 'polygon B_Z')
        end do
    end subroutine check_polygon

    function polygon(a, h, cur, r, z, n) result(v)
        real(dp), intent(in) :: a, h, cur, r, z
        integer, intent(in) :: n
        real(dp) :: v(3), t, x(n + 1), y(n + 1), zz(n + 1), c(n + 1), p(3, 1), f(3, 1)
        type(coil_set_t) :: cs
        integer :: k
        do k = 1, n + 1
            t = 2 * pi * (k - 1) / n
            x(k) = a * cos(t)
            y(k) = a * sin(t)
            zz(k) = h
        end do
        c = cur
        cs = build_coil_set(x, y, zz, c)
        p(:, 1) = [r, 0.0_dp, z]
        call cs%vector_potential(p, f)
        v(1) = r * f(2, 1)
        call cs%field(p, f)
        v(2:3) = [f(1, 1), f(3, 1)]
    end function polygon

    subroutine check_solovev()
        integer, parameter :: nrho = 40, ntheta = 256
        real(dp), allocatable :: rv(:), zv(:), wv(:), rb(:), zb(:), wb(:), sv(:), sb(:)
        real(dp) :: xg(nrho), wg(nrho), theta, rho, dr, dz
        type(axisym_sensors_t) :: sensors
        integer :: i, j, k
        allocate(rv(nrho * ntheta), zv(nrho * ntheta), wv(nrho * ntheta))
        allocate(rb(ntheta), zb(ntheta), wb(ntheta))
        call gauss_legendre_unit(nrho, xg, wg)
        k = 0
        do j = 1, ntheta
            theta = 2 * pi * (j - 1) / ntheta
            rho = boundary_radius(theta)
            do i = 1, nrho
                k = k + 1
                rv(k) = r0 + rho * xg(i) * cos(theta)
                zv(k) = rho * xg(i) * sin(theta)
                wv(k) = -(c1 * rv(k)**2 + c2) / (mu0 * rv(k)) &
                    * rho**2 * xg(i) * wg(i) * 2 * pi / ntheta
            end do
            ! Boundary point and its tangent from the implicit dpsi/dtheta = 0.
            rb(j) = r0 + rho * cos(theta)
            zb(j) = rho * sin(theta)
            call boundary_tangent(theta, rho, dr, dz)
            wb(j) = -(dpsi_dr(rb(j)) * dz - dpsi_dz(zb(j)) * dr) / (mu0 * rb(j)) &
                * 2 * pi / ntheta
        end do
        call setup(sensors)
        sv = sensors%signals(rv, zv, wv)
        sb = sensors%signals(rb, zb, wb)
        do i = 1, size(sv)
            call close_to(sv(i), sb(i), 1.0e-11_dp, 'Solovev volume against sheet')
        end do
        call close_to(sum(wv), sum(wb), 1.0e-11_dp, 'Solovev current against circulation')
        call close_to(sv(size(sv)), sum(wv), 1.0e-14_dp, 'Rogowski encloses the plasma')
        call close_to(sum(wv), solovev_current(), 1.0e-11_dp, 'Solovev current closed form')
    end subroutine check_solovev

    subroutine setup(sensors)
        type(axisym_sensors_t), intent(out) :: sensors
        integer :: i
        allocate(sensors%flux(3), sensors%pickup(4), sensors%rogowski(1))
        sensors%flux(1)%r = 1.45_dp
        sensors%flux(1)%z = 0.0_dp
        sensors%flux(2)%r = 0.45_dp
        sensors%flux(2)%z = 0.2_dp
        sensors%flux(3)%r = 1.0_dp
        sensors%flux(3)%z = 0.6_dp
        do i = 1, 4
            sensors%pickup(i)%r = 1.0_dp + 0.6_dp * cos(pi * (i - 0.5_dp) / 2)
            sensors%pickup(i)%z = 0.6_dp * sin(pi * (i - 0.5_dp) / 2)
            sensors%pickup(i)%tr = -sin(pi * (i - 0.5_dp) / 2)
            sensors%pickup(i)%tz = cos(pi * (i - 0.5_dp) / 2)
            sensors%pickup(i)%nr = sensors%pickup(i)%tr
            sensors%pickup(i)%nz = sensors%pickup(i)%tz
            sensors%pickup(i)%length = 0.05_dp
            sensors%pickup(i)%points = 4
            sensors%pickup(i)%area = 0.12_dp
        end do
        sensors%rogowski(1)%r = [0.4_dp, 1.6_dp, 1.6_dp, 0.4_dp]
        sensors%rogowski(1)%z = [-0.7_dp, -0.7_dp, 0.7_dp, 0.7_dp]
    end subroutine setup

    subroutine check_coils()
        real(dp), allocatable :: rf(:), zf(:), wf(:), a(:, :), s(:), one(:)
        real(dp), allocatable :: r2(:), z2(:), w2(:)
        type(axisym_sensors_t) :: sensors
        real(dp) :: psi, br, bz, avg, x
        integer :: k, group(2 * 25 + 1)
        call axisym_rectangle(1.75_dp, 0.6_dp, 0.1_dp, 0.2_dp, 12.0_dp, 5, rf, zf, wf)
        call close_to(sum(wf), 12.0_dp, 1.0e-13_dp, 'rectangle turns')
        call axisym_rectangle(1.75_dp, -0.6_dp, 0.1_dp, 0.2_dp, -7.0_dp, 5, r2, z2, w2)
        rf = [rf, r2, 0.9_dp]
        zf = [zf, z2, 0.65_dp]
        wf = [wf, w2, 3.0_dp]
        group = [(1, k = 1, 25), (2, k = 1, 25), 2]
        call setup(sensors)
        a = sensors%response(rf, zf, wf, group, 2)
        s = sensors%signals(rf, zf, wf * merge(4.0e3_dp, -1.5e3_dp, group == 1))
        one = matmul(a, [4.0e3_dp, -1.5e3_dp])
        do k = 1, size(s)
            call close_to(one(k), s(k), 1.0e-13_dp, 'group response linearity')
        end do
        ! The Rogowski contour holds only the 3-turn filament of group 2.
        call close_to(a(size(s), 2), 3.0_dp, 1.0e-15_dp, 'Rogowski turns')
        ! Pickup length average against a 2001-point trapezoid rule.
        associate (p => sensors%pickup(2))
            avg = 0
            do k = 0, 2000
                x = -p%length / 2 + p%length * k / 2000
                call axisym_flux_field(p%r + x * p%tr, p%z + x * p%tz, rf, zf, wf, psi, br, bz)
                avg = avg + merge(0.5_dp, 1.0_dp, k == 0 .or. k == 2000) &
                    * (br * p%nr + bz * p%nz) / 2000
            end do
            s = sensors%signals(rf, zf, wf)
            call close_to(s(size(sensors%flux) + 2), p%area * avg, 1.0e-7_dp, 'pickup length average')
        end associate
    end subroutine check_coils

    real(dp) function psi_s(r, z)
        real(dp), intent(in) :: r, z
        psi_s = c1 * (r**2 - r0**2)**2 / 8 + c2 * z**2 / 2 - c3
    end function psi_s

    real(dp) function dpsi_dr(r)
        real(dp), intent(in) :: r
        dpsi_dr = c1 * r * (r**2 - r0**2) / 2
    end function dpsi_dr

    real(dp) function dpsi_dz(z)
        real(dp), intent(in) :: z
        dpsi_dz = c2 * z
    end function dpsi_dz

    real(dp) function boundary_radius(theta) result(rho)
        !! Root of psi_s along the ray from (r0, 0), by bisection then Newton.
        real(dp), intent(in) :: theta
        real(dp) :: lo, hi, f
        integer :: it
        lo = 0
        hi = 0.9_dp
        do it = 1, 60
            rho = (lo + hi) / 2
            if (psi_s(r0 + rho * cos(theta), rho * sin(theta)) < 0) then
                lo = rho
            else
                hi = rho
            end if
        end do
        do it = 1, 3
            f = psi_s(r0 + rho * cos(theta), rho * sin(theta))
            rho = rho - f / (dpsi_dr(r0 + rho * cos(theta)) * cos(theta) &
                + dpsi_dz(rho * sin(theta)) * sin(theta))
        end do
    end function boundary_radius

    subroutine boundary_tangent(theta, rho, dr, dz)
        real(dp), intent(in) :: theta, rho
        real(dp), intent(out) :: dr, dz
        real(dp) :: pr, pz, drho
        pr = dpsi_dr(r0 + rho * cos(theta))
        pz = dpsi_dz(rho * sin(theta))
        drho = -(-pr * rho * sin(theta) + pz * rho * cos(theta)) &
            / (pr * cos(theta) + pz * sin(theta))
        dr = drho * cos(theta) - rho * sin(theta)
        dz = drho * sin(theta) + rho * cos(theta)
    end subroutine boundary_tangent

    real(dp) function solovev_current() result(total)
        !! I = -(1/mu0) int (c1 R + c2/R) dR dZ with the half height
        !! zh(R) = sqrt(2 (c3 - c1 (R^2 - r0^2)^2 / 8) / c2), by 400-point
        !! Gauss-Legendre in u with R^2 - r0^2 = sqrt(8 c3/c1) sin(u).
        integer, parameter :: n = 400
        real(dp) :: x(n), w(n), u, s, r, zh, drdu
        integer :: i
        call gauss_legendre_unit(n, x, w)
        s = sqrt(8 * c3 / c1)
        total = 0
        do i = 1, n
            u = -pi / 2 + pi * x(i)
            r = sqrt(r0**2 + s * sin(u))
            drdu = s * cos(u) / (2 * r)
            zh = sqrt(2 * c3 / c2) * cos(u)
            total = total - (c1 * r + c2 / r) * 2 * zh * drdu * pi * w(i) / mu0
        end do
    end function solovev_current

    subroutine gauss_legendre_unit(n, x, w)
        !! n-point Gauss-Legendre rule on [0, 1].
        use fortnum_quadrature, only: gauss_legendre_ab
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)
        call gauss_legendre_ab(n, 0.0_dp, 1.0_dp, x, w)
    end subroutine gauss_legendre_unit

    subroutine close_to(value, reference, tolerance, label)
        real(dp), intent(in) :: value, reference, tolerance
        character(*), intent(in) :: label
        if (abs(value - reference) > tolerance * max(abs(reference), 1.0e-300_dp)) then
            write(error_unit, '(A,": ",ES23.15," vs ",ES23.15)') label, value, reference
            failures = failures + 1
        end if
    end subroutine close_to

end program test_axisym_magnetics
