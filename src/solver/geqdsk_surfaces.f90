module tiago_geqdsk_surfaces
    !! Flux-surface averages on a G-EQDSK equilibrium (docs/bootstrap.md).
    !! The surface psi_N = x is traced along rays from the magnetic axis,
    !! R = R_axis + rho(theta) cos(theta), Z = Z_axis + rho(theta) sin(theta),
    !! and averaged with the weight dl / B_p,
    !!     <A> = int A dl/B_p / int dl/B_p,
    !! by the periodic trapezoidal rule in theta, with
    !!     d rho / d theta = rho (psi_R sin - psi_Z cos) / (psi_R cos + psi_Z sin).
    !! The effective trapped fraction is the bounce integral of Lin-Liu and
    !! Miller, Phys. Plasmas 2, 1666 (1995), Eq. 1,
    !!     f_t = 1 - 3/4 <B^2> int_0^{1/B_max} lambda dlambda / <sqrt(1 - lambda B)>,
    !! after lambda = (1 - y^2) / B_max with Gauss-Legendre points in y.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use batch_interpolate, only: evaluate_batch_splines_1d
    use tiago_geqdsk_plasma, only: geqdsk_plasma_t
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    type, public :: flux_surface_t
        real(dp) :: psin = 0
        real(dp) :: r_min = 0, r_max = 0 !! extent of the surface [m]
        real(dp) :: r_geo = 0, epsilon = 0 !! (r_max + r_min) / 2 [m], (r_max - r_min) / (r_max + r_min)
        real(dp) :: f = 0 !! F = R B_phi [T m], sign of the helicity
        real(dp) :: q = 0 !! |q|
        real(dp) :: b2 = 0, b_max = 0 !! <B^2> [T^2], max |B| [T]
        real(dp) :: f_trap = 0 !! effective trapped fraction
    end type flux_surface_t

    public :: flux_surface, surface_jdotb

contains

    subroutine flux_surface(eq, psin, surf, ntheta)
        !! Averages of the surface psi_N = psin in (0, 1), ntheta rays (default 512).
        type(geqdsk_plasma_t), intent(in) :: eq
        real(dp), intent(in) :: psin
        type(flux_surface_t), intent(out) :: surf
        integer, intent(in), optional :: ntheta
        integer, parameter :: ny = 48
        real(dp), allocatable :: r(:), b(:), w(:)
        real(dp) :: theta, rho, drho, c, s, x, a, gz, y(5), dpsi, bp2, rr, zz, dr, dz
        real(dp) :: yq(ny), wq(ny), h, total
        integer :: n, j, k

        n = 512
        if (present(ntheta)) n = ntheta
        allocate(r(n), b(n), w(n))
        call evaluate_batch_splines_1d(eq%profile_spline, psin, y)
        surf%psin = psin
        surf%q = y(3)
        surf%f = eq%helicity * y(4)
        dpsi = eq%psi_b - eq%psi_axis
        do j = 1, n
            theta = 2 * pi * (j - 1) / n
            c = cos(theta)
            s = sin(theta)
            rho = ray_root()
            rr = eq%r_axis + rho * c
            zz = eq%z_axis + rho * s
            call eq%label(rr, zz, x, a, gz)
            drho = rho * (a * s - gz * c) / (a * c + gz * s)
            dr = drho * c - rho * s
            dz = drho * s + rho * c
            bp2 = (a**2 + gz**2) * dpsi**2 / rr**2
            r(j) = rr
            b(j) = sqrt(bp2 + (y(4) / rr)**2)
            w(j) = sqrt(dr**2 + dz**2) / sqrt(bp2)
        end do
        total = sum(w)
        surf%r_min = minval(r)
        surf%r_max = maxval(r)
        surf%r_geo = (surf%r_max + surf%r_min) / 2
        surf%epsilon = (surf%r_max - surf%r_min) / (surf%r_max + surf%r_min)
        surf%b2 = sum(b**2 * w) / total
        surf%b_max = maxval(b)
        ! int_0^{1/B_max} lambda dlambda / <sqrt(1 - lambda B)>
        !   = 2 / B_max^2 int_0^1 y (1 - y^2) / <sqrt(1 - (1 - y^2) B / B_max)> dy
        call gauss_legendre(ny, yq, wq)
        yq = (yq + 1) / 2
        wq = wq / 2
        h = 0
        do k = 1, ny
            h = h + wq(k) * yq(k) * (1 - yq(k)**2) &
                / (sum(sqrt(max(1 - (1 - yq(k)**2) * b / surf%b_max, 0.0_dp)) * w) / total)
        end do
        surf%f_trap = 1 - 0.75_dp * surf%b2 * 2 * h / surf%b_max**2

    contains

        real(dp) function ray_root() result(t)
            !! First rho where psi_N = psin along the ray: marching out from the
            !! axis in steps of 2 mm, so that a ray past the X-point does not
            !! bracket a crossing in the private flux region, then bisection
            !! refined with Newton steps.
            real(dp), parameter :: step = 2.0e-3_dp
            real(dp) :: lo, hi, f, gr, gzz
            integer :: it
            lo = 0
            hi = step
            do
                call eq%label(eq%r_axis + hi * c, eq%z_axis + hi * s, f, gr, gzz)
                if (f >= psin) exit
                lo = hi
                hi = hi + step
            end do
            t = (lo + hi) / 2
            do it = 1, 100
                call eq%label(eq%r_axis + t * c, eq%z_axis + t * s, f, gr, gzz)
                if (f > psin) then
                    hi = t
                else
                    lo = t
                end if
                t = t - (f - psin) / (gr * c + gzz * s)
                if (t <= lo .or. t >= hi) t = (lo + hi) / 2
                if (hi - lo < 1.0e-13_dp .or. abs(f - psin) < 1.0e-14_dp) exit
            end do
        end function ray_root
    end subroutine flux_surface

    real(dp) function surface_jdotb(eq, surf) result(jb)
        !! <j.B> = F p' + F' <B^2> / mu0 [A T / m^2] of the equilibrium on the
        !! surface, from j = F' B / mu0 + p' R e_phi with B_p = grad psi x grad phi.
        type(geqdsk_plasma_t), intent(in) :: eq
        type(flux_surface_t), intent(in) :: surf
        real(dp), parameter :: mu0 = 4.0e-7_dp * pi
        real(dp) :: y(5)
        call evaluate_batch_splines_1d(eq%profile_spline, surf%psin, y)
        jb = surf%f * y(1) + y(2) / surf%f * surf%b2 / mu0
    end function surface_jdotb

end module tiago_geqdsk_surfaces
