module tiago_flux_map
    !! Flux label s and |B| at real-space points from a stellarator-symmetric
    !! VMEC geometry, and the reverse-mode derivative of s with respect to the
    !! geometry coefficients.
    !!
    !! Geometry: VMEC++ physical product-basis coefficients coef(0:ntor,
    !! 0:mpol-1, ns, 12) (tiago_vmecpp%geometry), R = sum r_cc cos(m t) cos(n z)
    !! + r_ss sin(m t) sin(n z), Z = sum z_sc sin(m t) cos(n z) + z_cs cos(m t)
    !! sin(n z), z = nfp phi, on the full grid s_j = (j - 1) / (ns - 1).
    !! Radially, even-m coefficients are interpolated linearly in s and odd-m
    !! coefficients as rho times a linear interpolant of c / rho (rho = sqrt(s),
    !! c / rho at the axis taken from the first surface off axis), continued
    !! linearly beyond s = 1. For a point (R, phi, Z) Newton's method solves
    !! R(rho, t, phi) = R, Z(rho, t, phi) = Z for (rho, t). Implicit
    !! differentiation of these two equations gives ds/dx and, for a cotangent
    !! s_bar, coef_bar = -lambda_R dR/dcoef - lambda_Z dZ/dcoef with
    !! lambda = J^-T (2 rho s_bar, 0), J = d(R, Z)/d(rho, t).
    !!
    !! |B| comes from the wout Nyquist spectrum bmnc (half grid, linear in s,
    !! constant beyond the outermost half surface) at (s, t, phi); its
    !! dependence on the equilibrium is not differentiated here (see
    !! docs/inference.md).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter :: R_CC = 1, R_SS = 2, Z_SC = 5, Z_CS = 6

    type, public :: flux_map_t
        integer :: ns = 0, mpol = 0, ntor = 0, nfp = 1
        real(dp), allocatable :: coef(:, :, :, :)          !! (0:ntor, 0:mpol-1, ns, 12)
        real(dp), allocatable :: rho_full(:)
        ! |B| spectrum (wout layout): bmnc(mn, j) on the half grid j = 2 ... ns
        real(dp), allocatable :: xm_nyq(:), xn_nyq(:), bmnc(:, :)
    contains
        procedure :: init => map_init
        procedure :: locate => map_locate
        procedure :: s_vjp => map_s_vjp
        procedure :: point => map_point
    end type flux_map_t

    type, public :: located_t
        !! Flux coordinates of points and the data the reverse mode needs.
        integer :: n = 0
        real(dp), allocatable :: rho(:), theta(:), phi(:)
        real(dp), allocatable :: jinv(:, :, :)     !! (2, 2, n) inverse of d(R, Z)/d(rho, theta)
        real(dp), allocatable :: ds_dx(:, :)       !! (3, n) d s / d (R, phi, Z)
        real(dp), allocatable :: b(:), db_dx(:, :) !! |B| and d|B| / d (R, phi, Z)
        logical, allocatable :: converged(:)
    end type located_t

contains

    subroutine map_init(self, coef, nfp, xm_nyq, xn_nyq, bmnc)
        class(flux_map_t), intent(inout) :: self
        real(dp), intent(in) :: coef(0:, 0:, :, :)
        integer, intent(in) :: nfp
        real(dp), intent(in), optional :: xm_nyq(:), xn_nyq(:), bmnc(:, :)
        integer :: j

        self%ntor = size(coef, 1) - 1
        self%mpol = size(coef, 2)
        self%ns = size(coef, 3)
        self%nfp = nfp
        self%coef = coef
        self%rho_full = [(sqrt(real(j - 1, dp) / (self%ns - 1)), j = 1, self%ns)]
        if (present(bmnc)) then
            self%xm_nyq = xm_nyq
            self%xn_nyq = xn_nyq
            self%bmnc = bmnc
        end if
    end subroutine map_init

    pure subroutine radial_weights(self, rho, odd, j, w, dw)
        !! c(rho) = w(1) c_j + w(2) c_(j+1) (odd m: c_1 replaced by c_2) and dw/drho.
        type(flux_map_t), intent(in) :: self
        real(dp), intent(in) :: rho
        logical, intent(in) :: odd
        integer, intent(out) :: j
        real(dp), intent(out) :: w(2), dw(2)
        real(dp) :: ds, s, t, ra, rb

        ds = 1.0_dp / (self%ns - 1)
        s = rho**2
        j = min(max(int(s / ds) + 1, 1), self%ns - 1)
        t = (s - (j - 1) * ds) / ds
        if (.not. odd) then
            w = [1.0_dp - t, t]
            dw = [-2.0_dp * rho / ds, 2.0_dp * rho / ds]
        else
            ra = self%rho_full(j)
            rb = self%rho_full(j + 1)
            if (j == 1) ra = rb       ! c / rho at the axis from the first surface
            w = [rho * (1.0_dp - t) / ra, rho * t / rb]
            dw = [((1.0_dp - t) - 2.0_dp * rho**2 / ds) / ra, (t + 2.0_dp * rho**2 / ds) / rb]
        end if
    end subroutine radial_weights

    pure subroutine map_point(self, rho, theta, phi, r, z, d)
        !! R, Z at (rho, theta, phi) and d(R, Z)/d(rho, theta, phi) in d(2, 3).
        class(flux_map_t), intent(in) :: self
        real(dp), intent(in) :: rho, theta, phi
        real(dp), intent(out) :: r, z, d(2, 3)
        real(dp) :: w(2), dw(2), cm, sm, cn, sn, zeta, nn, rcc, rss, zsc, zcs
        real(dp) :: drcc, drss, dzsc, dzcs
        integer :: m, n, j, jj

        r = 0.0_dp
        z = 0.0_dp
        d = 0.0_dp
        zeta = self%nfp * phi
        do m = 0, self%mpol - 1
            call radial_weights(self, rho, mod(m, 2) == 1, j, w, dw)
            jj = j
            if (mod(m, 2) == 1 .and. j == 1) jj = 2
            cm = cos(m * theta)
            sm = sin(m * theta)
            do n = 0, self%ntor
                nn = real(self%nfp * n, dp)
                cn = cos(n * zeta)
                sn = sin(n * zeta)
                rcc = w(1) * self%coef(n, m, jj, R_CC) + w(2) * self%coef(n, m, j + 1, R_CC)
                rss = w(1) * self%coef(n, m, jj, R_SS) + w(2) * self%coef(n, m, j + 1, R_SS)
                zsc = w(1) * self%coef(n, m, jj, Z_SC) + w(2) * self%coef(n, m, j + 1, Z_SC)
                zcs = w(1) * self%coef(n, m, jj, Z_CS) + w(2) * self%coef(n, m, j + 1, Z_CS)
                drcc = dw(1) * self%coef(n, m, jj, R_CC) + dw(2) * self%coef(n, m, j + 1, R_CC)
                drss = dw(1) * self%coef(n, m, jj, R_SS) + dw(2) * self%coef(n, m, j + 1, R_SS)
                dzsc = dw(1) * self%coef(n, m, jj, Z_SC) + dw(2) * self%coef(n, m, j + 1, Z_SC)
                dzcs = dw(1) * self%coef(n, m, jj, Z_CS) + dw(2) * self%coef(n, m, j + 1, Z_CS)
                r = r + rcc * cm * cn + rss * sm * sn
                z = z + zsc * sm * cn + zcs * cm * sn
                d(1, 1) = d(1, 1) + drcc * cm * cn + drss * sm * sn
                d(2, 1) = d(2, 1) + dzsc * sm * cn + dzcs * cm * sn
                d(1, 2) = d(1, 2) - m * rcc * sm * cn + m * rss * cm * sn
                d(2, 2) = d(2, 2) + m * zsc * cm * cn - m * zcs * sm * sn
                d(1, 3) = d(1, 3) - nn * rcc * cm * sn + nn * rss * sm * cn
                d(2, 3) = d(2, 3) - nn * zsc * sm * sn + nn * zcs * cm * cn
            end do
        end do
    end subroutine map_point

    subroutine map_locate(self, x, loc)
        !! Flux coordinates of points x(3, n) = (R, phi, Z).
        class(flux_map_t), intent(in) :: self
        real(dp), intent(in) :: x(:, :)
        type(located_t), intent(out) :: loc
        integer :: k, n

        n = size(x, 2)
        loc%n = n
        allocate(loc%rho(n), loc%theta(n), loc%phi(n), loc%jinv(2, 2, n), loc%ds_dx(3, n), &
            loc%b(n), loc%db_dx(3, n), loc%converged(n))
        !$omp parallel do schedule(dynamic)
        do k = 1, n
            call locate_one(self, x(:, k), loc%rho(k), loc%theta(k), loc%jinv(:, :, k), &
                loc%ds_dx(:, k), loc%converged(k))
            loc%phi(k) = x(2, k)
            call field_strength(self, loc%rho(k), loc%theta(k), x(2, k), loc%jinv(:, :, k), &
                loc%ds_dx(:, k), loc%b(k), loc%db_dx(:, k))
        end do
        !$omp end parallel do
    end subroutine map_locate

    subroutine locate_one(self, x, rho, theta, jinv, ds_dx, converged)
        type(flux_map_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: rho, theta, jinv(2, 2), ds_dx(3)
        logical, intent(out) :: converged
        real(dp) :: r, z, d(2, 3), f(2), step(2), det, best, dist, r0, z0, rt, tt, lam, fnew(2)
        real(dp) :: rho_try, theta_try, scale
        integer :: it, i, k, back

        ! start: closest point of a coarse (rho, theta) grid
        best = huge(1.0_dp)
        rho = 0.5_dp
        theta = 0.0_dp
        do i = 1, 12
            rt = i / 10.0_dp
            do k = 0, 31
                tt = k * 2.0_dp * acos(-1.0_dp) / 32
                call self%point(rt, tt, x(2), r, z, d)
                dist = (r - x(1))**2 + (z - x(3))**2
                if (dist < best) then
                    best = dist
                    rho = rt
                    theta = tt
                end if
            end do
        end do
        call self%point(0.0_dp, 0.0_dp, x(2), r0, z0, d)
        scale = max(abs(r0), 1.0_dp)
        converged = .false.
        call self%point(rho, theta, x(2), r, z, d)
        f = [r - x(1), z - x(3)]
        do it = 1, 60
            det = d(1, 1) * d(2, 2) - d(1, 2) * d(2, 1)
            if (det == 0.0_dp) exit
            step = -[d(2, 2) * f(1) - d(1, 2) * f(2), -d(2, 1) * f(1) + d(1, 1) * f(2)] / det
            lam = 1.0_dp
            do back = 1, 30
                rho_try = rho + lam * step(1)
                theta_try = theta + lam * step(2)
                if (rho_try < 0.0_dp) then
                    rho_try = -rho_try
                    theta_try = theta_try + acos(-1.0_dp)
                end if
                call self%point(rho_try, theta_try, x(2), r, z, d)
                fnew = [r - x(1), z - x(3)]
                if (norm2(fnew) < norm2(f) .or. norm2(fnew) < 1.0e-13_dp * scale) exit
                lam = 0.5_dp * lam
            end do
            rho = rho_try
            theta = theta_try
            f = fnew
            if (norm2(f) < 1.0e-13_dp * scale) then
                converged = .true.
                exit
            end if
        end do
        call self%point(rho, theta, x(2), r, z, d)
        det = d(1, 1) * d(2, 2) - d(1, 2) * d(2, 1)
        jinv = reshape([d(2, 2), -d(2, 1), -d(1, 2), d(1, 1)], [2, 2]) / det
        ! s = rho^2: ds/dR = 2 rho drho/dR, ds/dZ = 2 rho drho/dZ,
        ! ds/dphi = -2 rho (J^-1 d(R, Z)/dphi)_rho
        ds_dx(1) = 2.0_dp * rho * jinv(1, 1)
        ds_dx(3) = 2.0_dp * rho * jinv(1, 2)
        ds_dx(2) = -2.0_dp * rho * (jinv(1, 1) * d(1, 3) + jinv(1, 2) * d(2, 3))
        theta = modulo(theta, 2.0_dp * acos(-1.0_dp))
    end subroutine locate_one

    subroutine field_strength(self, rho, theta, phi, jinv, ds_dx, b, db_dx)
        !! |B| and its spatial gradient from the Nyquist spectrum (0 without one).
        type(flux_map_t), intent(in) :: self
        real(dp), intent(in) :: rho, theta, phi, jinv(2, 2), ds_dx(3)
        real(dp), intent(out) :: b, db_dx(3)
        real(dp) :: s, ds, t, sh1, bm, dbm_ds, db_ds, db_dt, db_dp, arg, dtheta_dx(3)
        real(dp) :: r, z, d(2, 3)
        integer :: mn, j

        b = 0.0_dp
        db_dx = 0.0_dp
        if (.not. allocated(self%bmnc)) return
        s = rho**2
        ds = 1.0_dp / (self%ns - 1)
        ! half-grid surfaces j = 2 ... ns at s = (j - 1.5) ds
        sh1 = 0.5_dp * ds
        if (s <= sh1) then
            j = 2
            t = 0.0_dp
        else
            j = min(int((s - sh1) / ds) + 2, self%ns - 1)
            t = min((s - sh1 - (j - 2) * ds) / ds, 1.0_dp)
        end if
        db_ds = 0.0_dp
        db_dt = 0.0_dp
        db_dp = 0.0_dp
        do mn = 1, size(self%xm_nyq)
            bm = (1.0_dp - t) * self%bmnc(mn, j) + t * self%bmnc(mn, j + 1)
            dbm_ds = 0.0_dp
            if (s > sh1 .and. t < 1.0_dp) dbm_ds = (self%bmnc(mn, j + 1) - self%bmnc(mn, j)) / ds
            arg = self%xm_nyq(mn) * theta - self%xn_nyq(mn) * phi
            b = b + bm * cos(arg)
            db_ds = db_ds + dbm_ds * cos(arg)
            db_dt = db_dt - bm * self%xm_nyq(mn) * sin(arg)
            db_dp = db_dp + bm * self%xn_nyq(mn) * sin(arg)
        end do
        call self%point(rho, theta, phi, r, z, d)
        ! theta(R, phi, Z): dtheta/d(R, Z) = J^-1 row 2, dtheta/dphi = -(J^-1 dX/dphi)_theta
        dtheta_dx = [jinv(2, 1), -(jinv(2, 1) * d(1, 3) + jinv(2, 2) * d(2, 3)), jinv(2, 2)]
        db_dx = db_ds * ds_dx + db_dt * dtheta_dx
        db_dx(2) = db_dx(2) + db_dp
    end subroutine field_strength

    subroutine map_s_vjp(self, loc, s_bar, coef_bar)
        !! coef_bar += (ds/dcoef)^T s_bar for the located points.
        class(flux_map_t), intent(in) :: self
        type(located_t), intent(in) :: loc
        real(dp), intent(in) :: s_bar(:)
        real(dp), intent(inout) :: coef_bar(0:, 0:, :, :)
        real(dp) :: lam(2), w(2), dw(2), cm, sm, cn, sn, zeta
        integer :: k, m, n, j, jj

        !$omp parallel do private(lam, w, dw, cm, sm, cn, sn, zeta, m, n, j, jj) &
        !$omp reduction(+:coef_bar)
        do k = 1, loc%n
            if (s_bar(k) == 0.0_dp) cycle
            ! lambda = J^-T (2 rho s_bar, 0); coef_bar -= lambda . dF/dcoef
            lam = [loc%jinv(1, 1, k), loc%jinv(1, 2, k)] * 2.0_dp * loc%rho(k) * s_bar(k)
            zeta = self%nfp * loc%phi(k)
            do m = 0, self%mpol - 1
                call radial_weights(self, loc%rho(k), mod(m, 2) == 1, j, w, dw)
                jj = j
                if (mod(m, 2) == 1 .and. j == 1) jj = 2
                cm = cos(m * loc%theta(k))
                sm = sin(m * loc%theta(k))
                do n = 0, self%ntor
                    cn = cos(n * zeta)
                    sn = sin(n * zeta)
                    coef_bar(n, m, jj, R_CC) = coef_bar(n, m, jj, R_CC) - lam(1) * w(1) * cm * cn
                    coef_bar(n, m, j + 1, R_CC) = coef_bar(n, m, j + 1, R_CC) - lam(1) * w(2) * cm * cn
                    coef_bar(n, m, jj, R_SS) = coef_bar(n, m, jj, R_SS) - lam(1) * w(1) * sm * sn
                    coef_bar(n, m, j + 1, R_SS) = coef_bar(n, m, j + 1, R_SS) - lam(1) * w(2) * sm * sn
                    coef_bar(n, m, jj, Z_SC) = coef_bar(n, m, jj, Z_SC) - lam(2) * w(1) * sm * cn
                    coef_bar(n, m, j + 1, Z_SC) = coef_bar(n, m, j + 1, Z_SC) - lam(2) * w(2) * sm * cn
                    coef_bar(n, m, jj, Z_CS) = coef_bar(n, m, jj, Z_CS) - lam(2) * w(1) * cm * sn
                    coef_bar(n, m, j + 1, Z_CS) = coef_bar(n, m, j + 1, Z_CS) - lam(2) * w(2) * cm * sn
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine map_s_vjp

end module tiago_flux_map
