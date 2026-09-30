module tiago_vmec_edge
    !! The edge vector Tiago's plasma model reads, computed from a VMEC++
    !! equilibrium's geometry, and its reverse-mode derivative:
    !!
    !!     y = [ B^u_mn, B^v_mn (s = 1),  rmnc_mn, zmns_mn (s = 1) ].
    !!
    !! B^u, B^v follow VMEC++'s output stage (output_quantities.cc, or its JAX
    !! port vmecpp/autodiff_wout.py) on the two outermost half-grid surfaces:
    !! the half-grid Jacobian tau * R and the metric from the odd/even-m split
    !! of R, Z (odd m scaled by 1/sqrt(s)), B^v = phi' (1 + d lambda/d theta) / sqrt(g),
    !! B^u = (chi' - phi' d lambda/d zeta) / sqrt(g), where chi' is iota phi' or,
    !! for a prescribed current (ncurr = 1), the value that makes the surface
    !! average of B_u equal the enclosed current. Both are Fourier-analysed to
    !! VMEC's Nyquist modes and extrapolated to s = 1 as 1.5 b(ns) - 0.5 b(ns-1),
    !! as in the wout file. The boundary R, Z come from the outermost surface.
    !! Stellarator-symmetric equilibria only (VMEC++'s adjoint has the same limit).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: d_s_half_d_s_interp = 0.25_dp
    ! geometry blocks of tiago_vmecpp%geometry
    integer, parameter, public :: R_CC = 1, R_SS = 2, Z_SC = 5, Z_CS = 6, L_SC = 9, L_CS = 10

    type, public :: vmec_edge_t
        integer :: ns, mpol, ntor, nfp, ncurr
        integer :: ntheta_even, ntheta_r, nzeta, mnyq, nnyq
        integer :: mnmax, mnmax_nyq
        real(dp), allocatable :: xm(:), xn(:), xm_nyq(:), xn_nyq(:)   !! wout mode tables
        real(dp), allocatable :: theta(:), zeta(:), w_int(:)          !! reduced grid
        real(dp), allocatable :: sqrt_s_full(:), sqrt_s_half(:), s_full(:), odd_scale(:)
        real(dp), allocatable :: kernel(:, :, :)                      !! (mn_nyq, nzeta, ntheta_r)
    contains
        procedure :: init => edge_init
        procedure :: n_y => edge_n_y
        procedure :: forward => edge_forward
        procedure :: vjp => edge_vjp
    end type vmec_edge_t

    ! real-space values of R, Z, lambda on one full surface, parity-split
    type :: surface_t
        real(dp), allocatable :: r(:, :, :), z(:, :, :), l(:, :, :)   !! (nzeta, ntheta_r, 6)
    end type surface_t
    ! slots of surface_t arrays
    integer, parameter :: VE = 1, VO = 2, TE = 3, TO = 4, ZE = 5, ZO = 6

contains

    subroutine edge_init(self, ns, mpol, ntor, nfp, ncurr, ntheta_input, nzeta_input)
        !! Grid and mode sizes as VMEC++'s Sizes::computeDerivedSizes.
        class(vmec_edge_t), intent(inout) :: self
        integer, intent(in) :: ns, mpol, ntor, nfp, ncurr, ntheta_input, nzeta_input
        real(dp), allocatable :: cosmui(:, :), sinmui(:, :), cosnv(:, :), sinnv(:, :)
        real(dp) :: mscale, nscale, dmult, int_norm, sgn
        integer :: ntheta, k, l, m, n, mn, an

        if (ns < 4) error stop 'VMEC edge field needs ns >= 4'
        self%ns = ns
        self%mpol = mpol
        self%ntor = ntor
        self%nfp = nfp
        self%ncurr = ncurr
        ntheta = max(ntheta_input, 2 * mpol + 6)
        self%ntheta_even = 2 * (ntheta / 2)
        self%ntheta_r = self%ntheta_even / 2 + 1
        if (ntor == 0) then
            self%nzeta = max(nzeta_input, 1)
        else
            self%nzeta = max(nzeta_input, 2 * ntor + 4)
        end if
        self%mnyq = max(0, self%ntheta_even / 2, mpol - 1)
        self%nnyq = max(0, self%nzeta / 2, ntor)
        call mode_table(mpol, ntor, nfp, self%xm, self%xn)
        call mode_table(self%mnyq + 1, self%nnyq, nfp, self%xm_nyq, self%xn_nyq)
        self%mnmax = size(self%xm)
        self%mnmax_nyq = size(self%xm_nyq)

        self%theta = [(2.0_dp * pi * k / self%ntheta_even, k = 0, self%ntheta_r - 1)]
        self%zeta = [(2.0_dp * pi * l / self%nzeta, l = 0, self%nzeta - 1)]
        allocate(self%w_int(self%ntheta_r))
        self%w_int = 1.0_dp / (self%nzeta * (self%ntheta_r - 1))
        self%w_int(1) = self%w_int(1) / 2.0_dp
        self%w_int(self%ntheta_r) = self%w_int(self%ntheta_r) / 2.0_dp
        self%s_full = [(real(k, dp) / (ns - 1), k = 0, ns - 1)]
        self%sqrt_s_full = sqrt(self%s_full)
        self%sqrt_s_full(ns) = 1.0_dp
        self%sqrt_s_half = [(sqrt((k + 0.5_dp) / (ns - 1)), k = 0, ns - 2)]
        self%odd_scale = 1.0_dp / max(self%sqrt_s_full, sqrt(1.0_dp / (ns - 1)))

        ! Nyquist forward-DFT kernel of cos(m theta - n zeta) (output_quantities)
        int_norm = 1.0_dp / (self%nzeta * (self%ntheta_r - 1))
        allocate(cosmui(0:self%mnyq, self%ntheta_r), sinmui(0:self%mnyq, self%ntheta_r))
        do m = 0, self%mnyq
            mscale = merge(1.0_dp, sqrt(2.0_dp), m == 0)
            cosmui(m, :) = cos(m * self%theta) * mscale * int_norm
            sinmui(m, :) = sin(m * self%theta) * mscale * int_norm
            cosmui(m, 1) = cosmui(m, 1) / 2.0_dp
            cosmui(m, self%ntheta_r) = cosmui(m, self%ntheta_r) / 2.0_dp
        end do
        if (self%mnyq /= 0) cosmui(self%mnyq, :) = cosmui(self%mnyq, :) / 2.0_dp
        allocate(cosnv(self%nzeta, 0:self%nnyq), sinnv(self%nzeta, 0:self%nnyq))
        do n = 0, self%nnyq
            nscale = merge(1.0_dp, sqrt(2.0_dp), n == 0)
            cosnv(:, n) = cos(n * self%zeta) * nscale
            sinnv(:, n) = sin(n * self%zeta) * nscale
        end do
        if (self%nnyq /= 0) cosnv(:, self%nnyq) = cosnv(:, self%nnyq) / 2.0_dp
        allocate(self%kernel(self%mnmax_nyq, self%nzeta, self%ntheta_r))
        do mn = 1, self%mnmax_nyq
            m = nint(self%xm_nyq(mn))
            n = nint(self%xn_nyq(mn)) / nfp
            an = abs(n)
            sgn = sign(1.0_dp, real(n, dp))
            if (n == 0) sgn = 0.0_dp
            dmult = merge(1.0_dp, sqrt(2.0_dp), m == 0) * merge(1.0_dp, sqrt(2.0_dp), an == 0) * 0.5_dp
            if (m == 0 .or. n == 0) dmult = 2.0_dp * dmult
            do k = 1, self%ntheta_r
                do l = 1, self%nzeta
                    self%kernel(mn, l, k) = dmult * (cosnv(l, an) * cosmui(m, k) &
                        + sgn * sinnv(l, an) * sinmui(m, k))
                end do
            end do
        end do
    end subroutine edge_init

    integer function edge_n_y(self)
        class(vmec_edge_t), intent(in) :: self
        edge_n_y = 2 * self%mnmax_nyq + 2 * self%mnmax
    end function edge_n_y

    subroutine mode_table(m_size, n_size, nfp, xm, xn)
        integer, intent(in) :: m_size, n_size, nfp
        real(dp), allocatable, intent(out) :: xm(:), xn(:)
        integer :: m, n, k

        allocate(xm(n_size + 1 + (m_size - 1) * (2 * n_size + 1)))
        allocate(xn(size(xm)))
        k = 0
        do n = 0, n_size
            k = k + 1
            xm(k) = 0.0_dp
            xn(k) = real(nfp * n, dp)
        end do
        do m = 1, m_size - 1
            do n = -n_size, n_size
                k = k + 1
                xm(k) = real(m, dp)
                xn(k) = real(nfp * n, dp)
            end do
        end do
    end subroutine mode_table

    ! ---- forward ---------------------------------------------------------------

    subroutine synthesize(self, coef, j, phip, s)
        !! Real space of R, Z, lambda on full surface j (1-based), parity-split,
        !! with the odd-m 1/sqrt(s) scaling; lambda times phi'(j).
        class(vmec_edge_t), intent(in) :: self
        real(dp), intent(in) :: coef(0:, 0:, :, :), phip
        integer, intent(in) :: j
        type(surface_t), intent(out) :: s
        integer :: m, n, p
        real(dp) :: scale

        allocate(s%r(self%nzeta, self%ntheta_r, 6), s%z(self%nzeta, self%ntheta_r, 6), &
            s%l(self%nzeta, self%ntheta_r, 6))
        s%r = 0.0_dp
        s%z = 0.0_dp
        s%l = 0.0_dp
        do m = 0, self%mpol - 1
            p = merge(0, 1, mod(m, 2) == 0)   ! 0 even, 1 odd
            scale = merge(1.0_dp, self%odd_scale(j), p == 0)
            do n = 0, self%ntor
                call add(s%r, p, m, n, scale * coef(n, m, j, R_CC), 'c', 'c')
                call add(s%r, p, m, n, scale * coef(n, m, j, R_SS), 's', 's')
                call add(s%z, p, m, n, scale * coef(n, m, j, Z_SC), 's', 'c')
                call add(s%z, p, m, n, scale * coef(n, m, j, Z_CS), 'c', 's')
                call add(s%l, p, m, n, scale * phip * coef(n, m, j, L_SC), 's', 'c')
                call add(s%l, p, m, n, scale * phip * coef(n, m, j, L_CS), 'c', 's')
            end do
        end do
    contains
        subroutine add(x, parity, m, n, c, pol, tor)
            real(dp), intent(inout) :: x(:, :, :)
            integer, intent(in) :: parity, m, n
            real(dp), intent(in) :: c
            character, intent(in) :: pol, tor
            real(dp) :: pv(self%ntheta_r), pd(self%ntheta_r), tv(self%nzeta), td(self%nzeta)
            integer :: k
            if (c == 0.0_dp) return
            call basis(self, m, n, pol, tor, pv, pd, tv, td)
            do k = 1, self%ntheta_r
                x(:, k, VE + parity) = x(:, k, VE + parity) + c * pv(k) * tv
                x(:, k, TE + parity) = x(:, k, TE + parity) + c * pd(k) * tv
                x(:, k, ZE + parity) = x(:, k, ZE + parity) + c * pv(k) * td
            end do
        end subroutine add
    end subroutine synthesize

    subroutine basis(self, m, n, pol, tor, pv, pd, tv, td)
        !! Poloidal and toroidal factors and their theta and (geometric) phi derivatives.
        type(vmec_edge_t), intent(in) :: self
        integer, intent(in) :: m, n
        character, intent(in) :: pol, tor
        real(dp), intent(out) :: pv(:), pd(:), tv(:), td(:)
        real(dp) :: nn

        nn = real(self%nfp * n, dp)
        if (pol == 'c') then
            pv = cos(m * self%theta)
            pd = -m * sin(m * self%theta)
        else
            pv = sin(m * self%theta)
            pd = m * cos(m * self%theta)
        end if
        if (tor == 'c') then
            tv = cos(n * self%zeta)
            td = -nn * sin(n * self%zeta)
        else
            tv = sin(n * self%zeta)
            td = nn * cos(n * self%zeta)
        end if
    end subroutine basis

    subroutine half_field(self, a, b, h, phip_a, phip_b, phip_h, iota_h, current_h, bu, bv, &
            chip, parts)
        !! B^u, B^v on half surface h (between full surfaces a = h and b = h + 1).
        class(vmec_edge_t), intent(in) :: self
        type(surface_t), intent(in) :: a, b
        integer, intent(in) :: h
        real(dp), intent(in) :: phip_a, phip_b, phip_h, iota_h, current_h
        real(dp), intent(out) :: bu(:, :), bv(:, :), chip
        real(dp), allocatable, intent(out), optional :: parts(:, :, :)
        real(dp), allocatable :: r12(:, :), ru12(:, :), zu12(:, :), rs(:, :), zs(:, :)
        real(dp), allocatable :: tau(:, :), g(:, :), guu(:, :), guv(:, :), bul(:, :)
        real(dp) :: sh, si, so, ds, plasma_current, average
        integer :: k

        sh = self%sqrt_s_half(h)
        si = self%s_full(h)
        so = self%s_full(h + 1)
        ds = 1.0_dp / (self%ns - 1)
        r12 = hv(a%r(:, :, VE), b%r(:, :, VE), a%r(:, :, VO), b%r(:, :, VO))
        ru12 = hv(a%r(:, :, TE), b%r(:, :, TE), a%r(:, :, TO), b%r(:, :, TO))
        zu12 = hv(a%z(:, :, TE), b%z(:, :, TE), a%z(:, :, TO), b%z(:, :, TO))
        rs = ((b%r(:, :, VE) - a%r(:, :, VE)) + sh * (b%r(:, :, VO) - a%r(:, :, VO))) / ds
        zs = ((b%z(:, :, VE) - a%z(:, :, VE)) + sh * (b%z(:, :, VO) - a%z(:, :, VO))) / ds
        tau = ru12 * zs - rs * zu12 + d_s_half_d_s_interp * ( &
            b%r(:, :, TO) * b%z(:, :, VO) + a%r(:, :, TO) * a%z(:, :, VO) &
            - b%z(:, :, TO) * b%r(:, :, VO) - a%z(:, :, TO) * a%r(:, :, VO) &
            + (b%r(:, :, TE) * b%z(:, :, VO) + a%r(:, :, TE) * a%z(:, :, VO) &
            - b%z(:, :, TE) * b%r(:, :, VO) - a%z(:, :, TE) * a%r(:, :, VO)) / sh)
        g = tau * r12
        bv = hv(a%l(:, :, TE) + phip_a, b%l(:, :, TE) + phip_b, a%l(:, :, TO), b%l(:, :, TO)) / g
        bul = hv(-a%l(:, :, ZE), -b%l(:, :, ZE), -a%l(:, :, ZO), -b%l(:, :, ZO)) / g
        if (self%ncurr == 1) then
            guu = metric(a%r(:, :, TE), a%r(:, :, TO), a%r(:, :, TE), a%r(:, :, TO), &
                b%r(:, :, TE), b%r(:, :, TO), b%r(:, :, TE), b%r(:, :, TO)) &
                + metric(a%z(:, :, TE), a%z(:, :, TO), a%z(:, :, TE), a%z(:, :, TO), &
                b%z(:, :, TE), b%z(:, :, TO), b%z(:, :, TE), b%z(:, :, TO))
            guv = metric(a%r(:, :, TE), a%r(:, :, TO), a%r(:, :, ZE), a%r(:, :, ZO), &
                b%r(:, :, TE), b%r(:, :, TO), b%r(:, :, ZE), b%r(:, :, ZO)) &
                + metric(a%z(:, :, TE), a%z(:, :, TO), a%z(:, :, ZE), a%z(:, :, ZO), &
                b%z(:, :, TE), b%z(:, :, TO), b%z(:, :, ZE), b%z(:, :, ZO))
            plasma_current = 0.0_dp
            average = 0.0_dp
            do k = 1, self%ntheta_r
                plasma_current = plasma_current + self%w_int(k) * &
                    sum(guu(:, k) * bul(:, k) + guv(:, k) * bv(:, k))
                average = average + self%w_int(k) * sum(guu(:, k) / g(:, k))
            end do
            chip = (current_h - plasma_current) / average
        else
            chip = iota_h * phip_h
        end if
        bu = bul + chip / g
        if (present(parts)) then
            allocate(parts(self%nzeta, self%ntheta_r, 10))
            parts(:, :, 1) = r12
            parts(:, :, 2) = ru12
            parts(:, :, 3) = zu12
            parts(:, :, 4) = rs
            parts(:, :, 5) = zs
            parts(:, :, 6) = tau
            parts(:, :, 7) = g
            parts(:, :, 8) = bul
            if (self%ncurr == 1) then
                parts(:, :, 9) = guu
                parts(:, :, 10) = guv
            else
                parts(:, :, 9:10) = 0.0_dp
            end if
        end if
    contains
        pure function hv(ea, eb, oa, ob) result(v)
            real(dp), intent(in) :: ea(:, :), eb(:, :), oa(:, :), ob(:, :)
            real(dp) :: v(size(ea, 1), size(ea, 2))
            v = 0.5_dp * ((ea + eb) + sh * (oa + ob))
        end function hv
        pure function metric(ae_a, ao_a, be_a, bo_a, ae_b, ao_b, be_b, bo_b) result(v)
            real(dp), intent(in) :: ae_a(:, :), ao_a(:, :), be_a(:, :), bo_a(:, :)
            real(dp), intent(in) :: ae_b(:, :), ao_b(:, :), be_b(:, :), bo_b(:, :)
            real(dp) :: v(size(ae_a, 1), size(ae_a, 2))
            v = 0.5_dp * (ae_a * be_a + ae_b * be_b + si * ao_a * bo_a + so * ao_b * bo_b) &
                + 0.5_dp * sh * (ae_a * bo_a + ae_b * bo_b + be_a * ao_a + be_b * ao_b)
        end function metric
    end subroutine half_field

    subroutine edge_forward(self, coef, phip_full, phip_half, iota_half, current_half, y)
        !! y = [B^u_mn, B^v_mn, rmnc_mn, zmns_mn] at s = 1.
        class(vmec_edge_t), intent(in) :: self
        real(dp), intent(in) :: coef(0:, 0:, :, :), phip_full(:), phip_half(:)
        real(dp), intent(in) :: iota_half(:), current_half(:)
        real(dp), intent(out) :: y(:)
        type(surface_t) :: s(3)
        real(dp), allocatable :: bu(:, :), bv(:, :)
        real(dp) :: chip, weight
        integer :: q, h, nb

        nb = self%mnmax_nyq
        allocate(bu(self%nzeta, self%ntheta_r), bv(self%nzeta, self%ntheta_r))
        do q = 1, 3
            call synthesize(self, coef, self%ns - 3 + q, phip_full(self%ns - 3 + q), s(q))
        end do
        y = 0.0_dp
        do q = 1, 2
            h = self%ns - 3 + q          ! half surfaces ns-2, ns-1 (1-based)
            weight = merge(-0.5_dp, 1.5_dp, q == 1)
            call half_field(self, s(q), s(q + 1), h, phip_full(h), phip_full(h + 1), &
                phip_half(h), iota_half(h), current_half(h), bu, bv, chip)
            y(:nb) = y(:nb) + weight * analyse(bu)
            y(nb + 1:2 * nb) = y(nb + 1:2 * nb) + weight * analyse(bv)
        end do
        call edge_geometry(self, coef, y(2 * nb + 1:))
    contains
        function analyse(f) result(c)
            real(dp), intent(in) :: f(:, :)
            real(dp) :: c(self%mnmax_nyq)
            integer :: mn
            do mn = 1, self%mnmax_nyq
                c(mn) = sum(self%kernel(mn, :, :) * f)
            end do
        end function analyse
    end subroutine edge_forward

    subroutine edge_geometry(self, coef, yg)
        !! rmnc, zmns (wout modes) of the outermost surface from the product basis:
        !! cos m t cos n z = (cos(m t - n z) + cos(m t + n z)) / 2, etc.
        class(vmec_edge_t), intent(in) :: self
        real(dp), intent(in) :: coef(0:, 0:, :, :)
        real(dp), intent(out) :: yg(:)
        integer :: mn, m, n, an, j

        j = self%ns
        do mn = 1, self%mnmax
            m = nint(self%xm(mn))
            n = nint(self%xn(mn)) / self%nfp
            an = abs(n)
            if (m == 0) then
                yg(mn) = coef(an, 0, j, R_CC)
                yg(self%mnmax + mn) = -coef(an, 0, j, Z_CS)
            else if (n == 0) then
                yg(mn) = coef(0, m, j, R_CC)
                yg(self%mnmax + mn) = coef(0, m, j, Z_SC)
            else
                yg(mn) = 0.5_dp * (coef(an, m, j, R_CC) + sign(1, n) * coef(an, m, j, R_SS))
                yg(self%mnmax + mn) = 0.5_dp * (coef(an, m, j, Z_SC) - sign(1, n) * coef(an, m, j, Z_CS))
            end if
        end do
    end subroutine edge_geometry

    ! ---- reverse mode -------------------------------------------------------------

    subroutine edge_vjp(self, coef, phip_full, phip_half, iota_half, current_half, ybar, &
            coef_bar, iota_bar, current_bar)
        !! Cotangents of the geometry blocks and of the half-grid iota and current
        !! (the direct dependence; the solver's dependence is the adjoint's) for a
        !! cotangent ybar of y.
        class(vmec_edge_t), intent(in) :: self
        real(dp), intent(in) :: coef(0:, 0:, :, :), phip_full(:), phip_half(:)
        real(dp), intent(in) :: iota_half(:), current_half(:), ybar(:)
        real(dp), intent(out) :: coef_bar(0:, 0:, :, :), iota_bar(:), current_bar(:)
        type(surface_t) :: s(3), sb(3)
        real(dp), allocatable :: bu(:, :), bv(:, :), bub(:, :), bvb(:, :), parts(:, :, :)
        real(dp) :: chip, weight
        integer :: q, h, nb, mn, j

        nb = self%mnmax_nyq
        coef_bar = 0.0_dp
        iota_bar = 0.0_dp
        current_bar = 0.0_dp
        allocate(bu(self%nzeta, self%ntheta_r), bv(self%nzeta, self%ntheta_r))
        allocate(bub(self%nzeta, self%ntheta_r), bvb(self%nzeta, self%ntheta_r))
        do q = 1, 3
            call synthesize(self, coef, self%ns - 3 + q, phip_full(self%ns - 3 + q), s(q))
            allocate(sb(q)%r, sb(q)%z, sb(q)%l, mold=s(q)%r)
            sb(q)%r = 0.0_dp
            sb(q)%z = 0.0_dp
            sb(q)%l = 0.0_dp
        end do
        do q = 1, 2
            h = self%ns - 3 + q
            weight = merge(-0.5_dp, 1.5_dp, q == 1)
            call half_field(self, s(q), s(q + 1), h, phip_full(h), phip_full(h + 1), &
                phip_half(h), iota_half(h), current_half(h), bu, bv, chip, parts)
            bub = 0.0_dp
            bvb = 0.0_dp
            do mn = 1, nb
                bub = bub + weight * ybar(mn) * self%kernel(mn, :, :)
                bvb = bvb + weight * ybar(nb + mn) * self%kernel(mn, :, :)
            end do
            call half_field_vjp(self, s(q), s(q + 1), sb(q), sb(q + 1), h, bv, chip, parts, &
                bub, bvb, iota_bar(h), current_bar(h), phip_half(h))
        end do
        do q = 1, 3
            j = self%ns - 3 + q
            call synthesize_vjp(self, coef_bar, j, phip_full(j), sb(q))
        end do
        call edge_geometry_vjp(self, ybar(2 * nb + 1:), coef_bar)
    end subroutine edge_vjp

    subroutine half_field_vjp(self, a, b, ab, bb, h, bv, chip, parts, bub, bvb, iota_bar, &
            current_bar, phip_h)
        !! Reverse of half_field: accumulates into the surface cotangents ab, bb.
        class(vmec_edge_t), intent(in) :: self
        type(surface_t), intent(in) :: a, b
        type(surface_t), intent(inout) :: ab, bb
        integer, intent(in) :: h
        real(dp), intent(in) :: bv(:, :), chip, parts(:, :, :), phip_h
        real(dp), intent(inout) :: bub(:, :), bvb(:, :)
        real(dp), intent(inout) :: iota_bar, current_bar
        real(dp), allocatable :: r12(:, :), ru12(:, :), zu12(:, :), rs(:, :), zs(:, :)
        real(dp), allocatable :: tau(:, :), g(:, :), bul(:, :), guu(:, :), guv(:, :)
        real(dp), allocatable :: gb(:, :), bulb(:, :), guub(:, :), guvb(:, :), hub(:, :), hvb(:, :)
        real(dp), allocatable :: taub(:, :), r12b(:, :), ru12b(:, :), zu12b(:, :), rsb(:, :), zsb(:, :)
        real(dp), allocatable :: t2b(:, :)
        real(dp) :: sh, si, so, ds, chipb, ipb, averageb, average
        integer :: k

        sh = self%sqrt_s_half(h)
        si = self%s_full(h)
        so = self%s_full(h + 1)
        ds = 1.0_dp / (self%ns - 1)
        r12 = parts(:, :, 1)
        ru12 = parts(:, :, 2)
        zu12 = parts(:, :, 3)
        rs = parts(:, :, 4)
        zs = parts(:, :, 5)
        tau = parts(:, :, 6)
        g = parts(:, :, 7)
        bul = parts(:, :, 8)
        guu = parts(:, :, 9)
        guv = parts(:, :, 10)
        allocate(gb, bulb, guub, guvb, mold=g)
        guub = 0.0_dp
        guvb = 0.0_dp
        ! bu = bul + chip / g
        bulb = bub
        chipb = 0.0_dp
        do k = 1, self%ntheta_r
            chipb = chipb + sum(bub(:, k) / g(:, k))
        end do
        gb = -bub * chip / g**2
        if (self%ncurr == 1) then
            average = 0.0_dp
            do k = 1, self%ntheta_r
                average = average + self%w_int(k) * sum(guu(:, k) / g(:, k))
            end do
            ! chip = (current - Ip) / average
            current_bar = current_bar + chipb / average
            ipb = -chipb / average
            averageb = -chipb * chip / average
            do k = 1, self%ntheta_r
                guub(:, k) = guub(:, k) + ipb * self%w_int(k) * bul(:, k) &
                    + averageb * self%w_int(k) / g(:, k)
                bulb(:, k) = bulb(:, k) + ipb * self%w_int(k) * guu(:, k)
                guvb(:, k) = guvb(:, k) + ipb * self%w_int(k) * bv(:, k)
                bvb(:, k) = bvb(:, k) + ipb * self%w_int(k) * guv(:, k)
                gb(:, k) = gb(:, k) - averageb * self%w_int(k) * guu(:, k) / g(:, k)**2
            end do
        else
            iota_bar = iota_bar + chipb * phip_h
        end if
        ! bv = hv(lu) / g, bul = hv(lv) / g
        hvb = bvb / g
        gb = gb - bvb * bv / g
        hub = bulb / g
        gb = gb - bulb * bul / g
        call hv_vjp(hvb, ab%l(:, :, TE), bb%l(:, :, TE), ab%l(:, :, TO), bb%l(:, :, TO))
        call hv_vjp(-hub, ab%l(:, :, ZE), bb%l(:, :, ZE), ab%l(:, :, ZO), bb%l(:, :, ZO))
        ! g = tau r12; tau = ru12 zs - rs zu12 + c tau2
        taub = gb * r12
        r12b = gb * tau
        ru12b = taub * zs
        zsb = taub * ru12
        rsb = -taub * zu12
        zu12b = -taub * rs
        t2b = d_s_half_d_s_interp * taub
        ab%r(:, :, TO) = ab%r(:, :, TO) + t2b * a%z(:, :, VO)
        ab%z(:, :, VO) = ab%z(:, :, VO) + t2b * a%r(:, :, TO)
        bb%r(:, :, TO) = bb%r(:, :, TO) + t2b * b%z(:, :, VO)
        bb%z(:, :, VO) = bb%z(:, :, VO) + t2b * b%r(:, :, TO)
        ab%z(:, :, TO) = ab%z(:, :, TO) - t2b * a%r(:, :, VO)
        ab%r(:, :, VO) = ab%r(:, :, VO) - t2b * a%z(:, :, TO)
        bb%z(:, :, TO) = bb%z(:, :, TO) - t2b * b%r(:, :, VO)
        bb%r(:, :, VO) = bb%r(:, :, VO) - t2b * b%z(:, :, TO)
        ab%r(:, :, TE) = ab%r(:, :, TE) + t2b * a%z(:, :, VO) / sh
        ab%z(:, :, VO) = ab%z(:, :, VO) + t2b * a%r(:, :, TE) / sh
        bb%r(:, :, TE) = bb%r(:, :, TE) + t2b * b%z(:, :, VO) / sh
        bb%z(:, :, VO) = bb%z(:, :, VO) + t2b * b%r(:, :, TE) / sh
        ab%z(:, :, TE) = ab%z(:, :, TE) - t2b * a%r(:, :, VO) / sh
        ab%r(:, :, VO) = ab%r(:, :, VO) - t2b * a%z(:, :, TE) / sh
        bb%z(:, :, TE) = bb%z(:, :, TE) - t2b * b%r(:, :, VO) / sh
        bb%r(:, :, VO) = bb%r(:, :, VO) - t2b * b%z(:, :, TE) / sh
        ! rs, zs radial differences
        ab%r(:, :, VE) = ab%r(:, :, VE) - rsb / ds
        bb%r(:, :, VE) = bb%r(:, :, VE) + rsb / ds
        ab%r(:, :, VO) = ab%r(:, :, VO) - sh * rsb / ds
        bb%r(:, :, VO) = bb%r(:, :, VO) + sh * rsb / ds
        ab%z(:, :, VE) = ab%z(:, :, VE) - zsb / ds
        bb%z(:, :, VE) = bb%z(:, :, VE) + zsb / ds
        ab%z(:, :, VO) = ab%z(:, :, VO) - sh * zsb / ds
        bb%z(:, :, VO) = bb%z(:, :, VO) + sh * zsb / ds
        ! half averages
        call hv_vjp(r12b, ab%r(:, :, VE), bb%r(:, :, VE), ab%r(:, :, VO), bb%r(:, :, VO))
        call hv_vjp(ru12b, ab%r(:, :, TE), bb%r(:, :, TE), ab%r(:, :, TO), bb%r(:, :, TO))
        call hv_vjp(zu12b, ab%z(:, :, TE), bb%z(:, :, TE), ab%z(:, :, TO), bb%z(:, :, TO))
        if (self%ncurr == 1) then
            call metric_vjp(guub, TE, TE, 'r')
            call metric_vjp(guub, TE, TE, 'z')
            call metric_vjp(guvb, TE, ZE, 'r')
            call metric_vjp(guvb, TE, ZE, 'z')
        end if
    contains
        subroutine hv_vjp(vb, ea, eb, oa, ob)
            real(dp), intent(in) :: vb(:, :)
            real(dp), intent(inout) :: ea(:, :), eb(:, :), oa(:, :), ob(:, :)
            ea = ea + 0.5_dp * vb
            eb = eb + 0.5_dp * vb
            oa = oa + 0.5_dp * sh * vb
            ob = ob + 0.5_dp * sh * vb
        end subroutine hv_vjp
        subroutine metric_vjp(vb, first, second, which)
            !! metric(x_first, x_second) of field which ('r' or 'z'); first and
            !! second name the even slot (the odd slot follows it).
            real(dp), intent(in) :: vb(:, :)
            integer, intent(in) :: first, second
            character, intent(in) :: which
            if (which == 'r') then
                call pair(vb, a%r, b%r, ab%r, bb%r, first, second)
            else
                call pair(vb, a%z, b%z, ab%z, bb%z, first, second)
            end if
        end subroutine metric_vjp
        subroutine pair(vb, xa, xb, xab, xbb, p, q)
            real(dp), intent(in) :: vb(:, :), xa(:, :, :), xb(:, :, :)
            real(dp), intent(inout) :: xab(:, :, :), xbb(:, :, :)
            integer, intent(in) :: p, q
            ! v = 0.5 (Pe_a Qe_a + Pe_b Qe_b + si Po_a Qo_a + so Po_b Qo_b)
            !   + 0.5 sh (Pe_a Qo_a + Pe_b Qo_b + Qe_a Po_a + Qe_b Po_b)
            xab(:, :, p) = xab(:, :, p) + vb * (0.5_dp * xa(:, :, q) + 0.5_dp * sh * xa(:, :, q + 1))
            xab(:, :, q) = xab(:, :, q) + vb * (0.5_dp * xa(:, :, p) + 0.5_dp * sh * xa(:, :, p + 1))
            xab(:, :, p + 1) = xab(:, :, p + 1) + vb * (0.5_dp * si * xa(:, :, q + 1) &
                + 0.5_dp * sh * xa(:, :, q))
            xab(:, :, q + 1) = xab(:, :, q + 1) + vb * (0.5_dp * si * xa(:, :, p + 1) &
                + 0.5_dp * sh * xa(:, :, p))
            xbb(:, :, p) = xbb(:, :, p) + vb * (0.5_dp * xb(:, :, q) + 0.5_dp * sh * xb(:, :, q + 1))
            xbb(:, :, q) = xbb(:, :, q) + vb * (0.5_dp * xb(:, :, p) + 0.5_dp * sh * xb(:, :, p + 1))
            xbb(:, :, p + 1) = xbb(:, :, p + 1) + vb * (0.5_dp * so * xb(:, :, q + 1) &
                + 0.5_dp * sh * xb(:, :, q))
            xbb(:, :, q + 1) = xbb(:, :, q + 1) + vb * (0.5_dp * so * xb(:, :, p + 1) &
                + 0.5_dp * sh * xb(:, :, p))
        end subroutine pair
    end subroutine half_field_vjp

    subroutine synthesize_vjp(self, coef_bar, j, phip, sb)
        !! Transpose of synthesize for surface j.
        class(vmec_edge_t), intent(in) :: self
        real(dp), intent(inout) :: coef_bar(0:, 0:, :, :)
        integer, intent(in) :: j
        real(dp), intent(in) :: phip
        type(surface_t), intent(in) :: sb
        integer :: m, n, p
        real(dp) :: scale

        do m = 0, self%mpol - 1
            p = merge(0, 1, mod(m, 2) == 0)
            scale = merge(1.0_dp, self%odd_scale(j), p == 0)
            do n = 0, self%ntor
                coef_bar(n, m, j, R_CC) = coef_bar(n, m, j, R_CC) + scale * dot(sb%r, p, m, n, 'c', 'c')
                coef_bar(n, m, j, R_SS) = coef_bar(n, m, j, R_SS) + scale * dot(sb%r, p, m, n, 's', 's')
                coef_bar(n, m, j, Z_SC) = coef_bar(n, m, j, Z_SC) + scale * dot(sb%z, p, m, n, 's', 'c')
                coef_bar(n, m, j, Z_CS) = coef_bar(n, m, j, Z_CS) + scale * dot(sb%z, p, m, n, 'c', 's')
                coef_bar(n, m, j, L_SC) = coef_bar(n, m, j, L_SC) &
                    + scale * phip * dot(sb%l, p, m, n, 's', 'c')
                coef_bar(n, m, j, L_CS) = coef_bar(n, m, j, L_CS) &
                    + scale * phip * dot(sb%l, p, m, n, 'c', 's')
            end do
        end do
    contains
        real(dp) function dot(x, parity, m, n, pol, tor)
            real(dp), intent(in) :: x(:, :, :)
            integer, intent(in) :: parity, m, n
            character, intent(in) :: pol, tor
            real(dp) :: pv(self%ntheta_r), pd(self%ntheta_r), tv(self%nzeta), td(self%nzeta)
            integer :: k
            call basis(self, m, n, pol, tor, pv, pd, tv, td)
            dot = 0.0_dp
            do k = 1, self%ntheta_r
                dot = dot + pv(k) * sum(x(:, k, VE + parity) * tv) &
                    + pd(k) * sum(x(:, k, TE + parity) * tv) + pv(k) * sum(x(:, k, ZE + parity) * td)
            end do
        end function dot
    end subroutine synthesize_vjp

    subroutine edge_geometry_vjp(self, ygb, coef_bar)
        class(vmec_edge_t), intent(in) :: self
        real(dp), intent(in) :: ygb(:)
        real(dp), intent(inout) :: coef_bar(0:, 0:, :, :)
        integer :: mn, m, n, an, j
        real(dp) :: rb, zb

        j = self%ns
        do mn = 1, self%mnmax
            m = nint(self%xm(mn))
            n = nint(self%xn(mn)) / self%nfp
            an = abs(n)
            rb = ygb(mn)
            zb = ygb(self%mnmax + mn)
            if (m == 0) then
                coef_bar(an, 0, j, R_CC) = coef_bar(an, 0, j, R_CC) + rb
                coef_bar(an, 0, j, Z_CS) = coef_bar(an, 0, j, Z_CS) - zb
            else if (n == 0) then
                coef_bar(0, m, j, R_CC) = coef_bar(0, m, j, R_CC) + rb
                coef_bar(0, m, j, Z_SC) = coef_bar(0, m, j, Z_SC) + zb
            else
                coef_bar(an, m, j, R_CC) = coef_bar(an, m, j, R_CC) + 0.5_dp * rb
                coef_bar(an, m, j, R_SS) = coef_bar(an, m, j, R_SS) + 0.5_dp * sign(1, n) * rb
                coef_bar(an, m, j, Z_SC) = coef_bar(an, m, j, Z_SC) + 0.5_dp * zb
                coef_bar(an, m, j, Z_CS) = coef_bar(an, m, j, Z_CS) - 0.5_dp * sign(1, n) * zb
            end if
        end do
    end subroutine edge_geometry_vjp
end module tiago_vmec_edge
