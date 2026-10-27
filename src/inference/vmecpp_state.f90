module tiago_vmecpp_state
    !! VMEC++ equilibrium as a real_space_state_t (tiago_real_space_state):
    !! labels u = (rho = sqrt(s), VMEC theta, cylindrical phi), x Cartesian.
    !!
    !! The parameters are the solver's physical geometry coefficients
    !! coef(0:ntor, 0:mpol-1, ns, 12) (tiago_vmecpp%geometry, product basis
    !! R = sum r_cc cos(m t) cos(n z) + r_ss sin(m t) sin(n z),
    !! Z = sum z_sc sin(m t) cos(n z) + z_cs cos(m t) sin(n z), z = nfp phi)
    !! followed by the |B| Nyquist spectrum bmnc on the half grid (wout layout,
    !! surfaces 2 ... ns). Every coefficient is a C1 cubic Hermite interpolant
    !! in rho of its node values (tiago_radial_interp, as tiago_vmec_state for
    !! wout files), extended beyond the last node. The geometry cotangent of
    !! param_vjp is the cotangent VMEC++'s implicit adjoint accepts
    !! (tiago_equilibrium%jacobian, coef_cot); the bmnc cotangent is its
    !! bmnc_cot (see docs/inference.md).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_real_space_state, only: real_space_state_t, vector_jets_t, locate_unsupported, locate_ok
    use tiago_radial_interp, only: hermite_weights, max_weights
    use tiago_vmec_state, only: vmec_state_t
    use tiago_vmec_product_pair, only: vmec_product_pair
    implicit none
    private

    integer, parameter :: R_CC = 1, R_SS = 2, Z_SC = 5, Z_CS = 6

    type, extends(real_space_state_t), public :: vmecpp_state_t
        integer :: ns = 0, mpol = 0, ntor = 0, nfp = 1, mnmax_nyq = 0
        real(dp), allocatable :: coef(:, :, :, :)          !! (0:ntor, 0:mpol-1, ns, 12)
        real(dp), allocatable :: xm_nyq(:), xn_nyq(:)
        real(dp), allocatable :: bmnc(:, :)                !! (mnmax_nyq, ns - 1): surfaces 2 ... ns
        real(dp), allocatable :: rho_full(:), rho_half(:)
        real(dp), allocatable :: bsupumnc(:,:), bsupvmnc(:,:)
    contains
        procedure :: vector_jets => state_vector_jets
        procedure :: init => state_init
        procedure :: n_geometry => state_n_geometry
        procedure :: n_params => state_n_params
        procedure :: get_params => state_get_params
        procedure :: set_params => state_set_params
        procedure :: eval => state_eval
        procedure :: param_jvp => state_param_jvp
        procedure :: param_vjp => state_param_vjp
        procedure :: point => state_point
    end type vmecpp_state_t

contains

    subroutine state_init(self, coef, nfp, xm_nyq, xn_nyq, bmnc, bsupumnc, bsupvmnc)
        !! bmnc: wout layout (mnmax_nyq, ns) with the axis column ignored.
        class(vmecpp_state_t), intent(inout) :: self
        real(dp), intent(in) :: coef(0:, 0:, :, :), xm_nyq(:), xn_nyq(:), bmnc(:, :)
        real(dp), intent(in), optional :: bsupumnc(:,:), bsupvmnc(:,:)
        integer, intent(in) :: nfp
        integer :: j

        self%ntor = size(coef, 1) - 1
        self%mpol = size(coef, 2)
        self%ns = size(coef, 3)
        self%nfp = nfp
        self%coef = coef
        self%xm_nyq = xm_nyq
        self%xn_nyq = xn_nyq
        self%mnmax_nyq = size(xm_nyq)
        self%bmnc = bmnc(:, 2:self%ns)
        if (allocated(self%bsupumnc)) deallocate(self%bsupumnc,self%bsupvmnc)
        if (present(bsupumnc) .neqv. present(bsupvmnc)) &
            error stop 'vmecpp_state: both contravariant spectra are required'
        if (present(bsupumnc)) then
            self%bsupumnc = bsupumnc(:,2:self%ns)
            self%bsupvmnc = bsupvmnc(:,2:self%ns)
        end if
        self%rho_full = [(sqrt(real(j - 1, dp) / (self%ns - 1)), j = 1, self%ns)]
        self%rho_half = [(sqrt((j - 1.5_dp) / (self%ns - 1)), j = 2, self%ns)]
    end subroutine state_init

    pure integer function state_n_geometry(self)
        class(vmecpp_state_t), intent(in) :: self
        state_n_geometry = size(self%coef)
    end function state_n_geometry

    pure integer function state_n_params(self)
        class(vmecpp_state_t), intent(in) :: self
        state_n_params = size(self%coef) + size(self%bmnc)
        if (allocated(self%bsupumnc)) state_n_params=state_n_params+2*size(self%bmnc)
    end function state_n_params

    pure subroutine state_get_params(self, p)
        class(vmecpp_state_t), intent(in) :: self
        real(dp), intent(out) :: p(:)
        integer :: ng
        ng = size(self%coef)
        p(1:ng) = reshape(self%coef, [ng])
        p(ng + 1:ng+size(self%bmnc)) = reshape(self%bmnc, [size(self%bmnc)])
        if (allocated(self%bsupumnc)) then
            p(ng+size(self%bmnc)+1:ng+2*size(self%bmnc)) = &
                reshape(self%bsupumnc,[size(self%bmnc)])
            p(ng+2*size(self%bmnc)+1:) = reshape(self%bsupvmnc,[size(self%bmnc)])
        end if
    end subroutine state_get_params

    subroutine state_set_params(self, p)
        class(vmecpp_state_t), intent(inout) :: self
        real(dp), intent(in) :: p(:)
        integer :: ng
        ng = size(self%coef)
        self%coef = reshape(p(1:ng), shape(self%coef))
        self%bmnc = reshape(p(ng + 1:ng+size(self%bmnc)), shape(self%bmnc))
        if (allocated(self%bsupumnc)) then
            self%bsupumnc=reshape(p(ng+size(self%bmnc)+1:ng+2*size(self%bmnc)),shape(self%bmnc))
            self%bsupvmnc=reshape(p(ng+2*size(self%bmnc)+1:),shape(self%bmnc))
        end if
    end subroutine state_set_params

    pure subroutine state_point(self, rho, theta, phi, r, z, d)
        !! R, Z at labels and d(R, Z)/d(rho, theta, phi) (cylindrical).
        class(vmecpp_state_t), intent(in) :: self
        real(dp), intent(in) :: rho, theta, phi
        real(dp), intent(out) :: r, z, d(2, 3)
        real(dp) :: w(max_weights), dw(max_weights), cm, sm, cn, sn, nn, zeta
        real(dp) :: rcc, rss, zsc, zcs, drcc, drss, dzsc, dzcs
        integer :: first, nw, m, n, k1, k2

        call hermite_weights(self%rho_full, rho, .true., first, nw, w, dw)
        k1 = first
        k2 = first + nw - 1
        r = 0.0_dp
        z = 0.0_dp
        d = 0.0_dp
        zeta = self%nfp * phi
        do m = 0, self%mpol - 1
            cm = cos(m * theta)
            sm = sin(m * theta)
            do n = 0, self%ntor
                nn = real(self%nfp * n, dp)
                cn = cos(n * zeta)
                sn = sin(n * zeta)
                rcc = dot_product(w(:nw), self%coef(n, m, k1:k2, R_CC))
                rss = dot_product(w(:nw), self%coef(n, m, k1:k2, R_SS))
                zsc = dot_product(w(:nw), self%coef(n, m, k1:k2, Z_SC))
                zcs = dot_product(w(:nw), self%coef(n, m, k1:k2, Z_CS))
                drcc = dot_product(dw(:nw), self%coef(n, m, k1:k2, R_CC))
                drss = dot_product(dw(:nw), self%coef(n, m, k1:k2, R_SS))
                dzsc = dot_product(dw(:nw), self%coef(n, m, k1:k2, Z_SC))
                dzcs = dot_product(dw(:nw), self%coef(n, m, k1:k2, Z_CS))
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
    end subroutine state_point

    pure subroutine state_eval(self, u, x, dxdu, b, dbdu)
        class(vmecpp_state_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3), dxdu(3, 3), b, dbdu(3)
        real(dp) :: r, z, d(2, 3), cp, sp, wh(max_weights), dwh(max_weights), c, s, cb, dcb, ang
        integer :: fh, nh, mn

        call self%point(u(1), u(2), u(3), r, z, d)
        cp = cos(u(3))
        sp = sin(u(3))
        x = [r * cp, r * sp, z]
        dxdu(1, :) = d(1, :) * cp
        dxdu(2, :) = d(1, :) * sp
        dxdu(3, :) = d(2, :)
        dxdu(1, 3) = dxdu(1, 3) - r * sp
        dxdu(2, 3) = dxdu(2, 3) + r * cp
        call hermite_weights(self%rho_half, u(1), .true., fh, nh, wh, dwh)
        b = 0.0_dp
        dbdu = 0.0_dp
        do mn = 1, self%mnmax_nyq
            ang = self%xm_nyq(mn) * u(2) - self%xn_nyq(mn) * u(3)
            c = cos(ang)
            s = sin(ang)
            cb = dot_product(wh(:nh), self%bmnc(mn, fh:fh + nh - 1))
            dcb = dot_product(dwh(:nh), self%bmnc(mn, fh:fh + nh - 1))
            b = b + cb * c
            dbdu = dbdu + [dcb * c, -self%xm_nyq(mn) * cb * s, self%xn_nyq(mn) * cb * s]
        end do
    end subroutine state_eval

    pure subroutine state_param_jvp(self, u, dp_, dx, db)
        class(vmecpp_state_t), intent(in) :: self
        real(dp), intent(in) :: u(3), dp_(:)
        real(dp), intent(out) :: dx(3), db
        real(dp) :: w(max_weights), dw(max_weights), wh(max_weights), dwh(max_weights)
        real(dp) :: cm, sm, cn, sn, zeta, dr, dz
        integer :: first, nw, fh, nh, m, n, k, j, mn, ng

        ng = size(self%coef)
        call hermite_weights(self%rho_full, u(1), .true., first, nw, w, dw)
        zeta = self%nfp * u(3)
        dr = 0.0_dp
        dz = 0.0_dp
        do m = 0, self%mpol - 1
            cm = cos(m * u(2))
            sm = sin(m * u(2))
            do n = 0, self%ntor
                cn = cos(n * zeta)
                sn = sin(n * zeta)
                do k = 1, nw
                    j = first + k - 1
                    dr = dr + w(k) * (cm * cn * dp_(idx(n, m, j, R_CC)) + sm * sn * dp_(idx(n, m, j, R_SS)))
                    dz = dz + w(k) * (sm * cn * dp_(idx(n, m, j, Z_SC)) + cm * sn * dp_(idx(n, m, j, Z_CS)))
                end do
            end do
        end do
        dx = [dr * cos(u(3)), dr * sin(u(3)), dz]
        call hermite_weights(self%rho_half, u(1), .true., fh, nh, wh, dwh)
        db = 0.0_dp
        do mn = 1, self%mnmax_nyq
            do k = 1, nh
                db = db + wh(k) * cos(self%xm_nyq(mn) * u(2) - self%xn_nyq(mn) * u(3)) * &
                    dp_(ng + mn + (fh + k - 2) * self%mnmax_nyq)
            end do
        end do
    contains
        pure integer function idx(n, m, j, blk)
            integer, intent(in) :: n, m, j, blk
            idx = 1 + n + (self%ntor + 1) * (m + self%mpol * (j - 1 + self%ns * (blk - 1)))
        end function idx
    end subroutine state_param_jvp

    pure subroutine state_param_vjp(self, u, xbar, bbar, pbar)
        class(vmecpp_state_t), intent(in) :: self
        real(dp), intent(in) :: u(3), xbar(3), bbar
        real(dp), intent(inout) :: pbar(:)
        real(dp) :: w(max_weights), dw(max_weights), wh(max_weights), dwh(max_weights)
        real(dp) :: cm, sm, cn, sn, zeta, rbar, zbar
        integer :: first, nw, fh, nh, m, n, k, j, mn, ng, i

        ng = size(self%coef)
        call hermite_weights(self%rho_full, u(1), .true., first, nw, w, dw)
        zeta = self%nfp * u(3)
        rbar = xbar(1) * cos(u(3)) + xbar(2) * sin(u(3))
        zbar = xbar(3)
        do m = 0, self%mpol - 1
            cm = cos(m * u(2))
            sm = sin(m * u(2))
            do n = 0, self%ntor
                cn = cos(n * zeta)
                sn = sin(n * zeta)
                do k = 1, nw
                    j = first + k - 1
                    i = idx(n, m, j, R_CC)
                    pbar(i) = pbar(i) + w(k) * cm * cn * rbar
                    i = idx(n, m, j, R_SS)
                    pbar(i) = pbar(i) + w(k) * sm * sn * rbar
                    i = idx(n, m, j, Z_SC)
                    pbar(i) = pbar(i) + w(k) * sm * cn * zbar
                    i = idx(n, m, j, Z_CS)
                    pbar(i) = pbar(i) + w(k) * cm * sn * zbar
                end do
            end do
        end do
        if (bbar == 0.0_dp) return
        call hermite_weights(self%rho_half, u(1), .true., fh, nh, wh, dwh)
        do mn = 1, self%mnmax_nyq
            do k = 1, nh
                i = ng + mn + (fh + k - 2) * self%mnmax_nyq
                pbar(i) = pbar(i) + wh(k) * cos(self%xm_nyq(mn) * u(2) - self%xn_nyq(mn) * u(3)) * bbar
            end do
        end do
    contains
        pure integer function idx(n, m, j, blk)
            integer, intent(in) :: n, m, j, blk
            idx = 1 + n + (self%ntor + 1) * (m + self%mpol * (j - 1 + self%ns * (blk - 1)))
        end function idx
    end subroutine state_param_vjp


    subroutine state_vector_jets(self, u, jet, status)
        class(vmecpp_state_t), intent(in)::self
        real(dp), intent(in)::u(3)
        type(vector_jets_t), intent(out)::jet
        integer, intent(out)::status
        type(vmec_state_t)::wout
        type(vector_jets_t)::local
        real(dp)::rc(2), zs(2), rr(2, 4), zr(2, 4)
        integer::n, m, j, k, a, b, ng, nb, nr, nparam, col, block_, offset
        integer, parameter::blocks(4) = [R_CC, R_SS, Z_SC, Z_CS]
        status = locate_unsupported
        if (.not. allocated(self%bsupumnc)) return
        wout%ns = self%ns
        wout%mnmax = 2*(self%ntor + 1)*self%mpol
        wout%mnmax_nyq = self%mnmax_nyq
        wout%nfp = self%nfp
        wout%rho_full = self%rho_full
        wout%rho_half = self%rho_half
        wout%xm_nyq = self%xm_nyq
        wout%xn_nyq = self%xn_nyq
        wout%bmnc = self%bmnc
        wout%bsupumnc = self%bsupumnc
        wout%bsupvmnc = self%bsupvmnc
        allocate (wout%xm(wout%mnmax), wout%xn(wout%mnmax), &
                  wout%rmnc(wout%mnmax, self%ns), wout%zmns(wout%mnmax, self%ns))
        do m = 0, self%mpol - 1
            do n = 0, self%ntor
                k = 1 + 2*(n + (self%ntor + 1)*m)
                wout%xm(k:k + 1) = m
                wout%xn(k:k + 1) = [n*self%nfp, -n*self%nfp]
                do j = 1, self%ns
                    call vmec_product_pair(self%coef(n, m, j, R_CC), self%coef(n, m, &
                        j, R_SS), &
                                           self%coef(n, m, j, Z_SC), self%coef(n, m, &
                                               j, Z_CS), rc, zs, rr, zr)
                    wout%rmnc(k:k + 1, j) = rc
                    wout%zmns(k:k + 1, j) = zs
                end do
            end do
        end do
        call wout%vector_jets(u, local, status)
        if (status /= locate_ok) return
        jet%x = local%x
        jet%xu = local%xu
        jet%xuu = local%xuu
        jet%b = local%b
        jet%bu = local%bu
        jet%buu = local%buu
        ng = size(self%coef)
        nb = size(self%bmnc)
        nr = wout%mnmax*self%ns
        nparam = self%n_params()
        allocate (jet%xp(3, nparam), jet%xup(3, 3, nparam), jet%bp(3, nparam), &
            jet%bup(3, 3, nparam))
        jet%xp = 0
        jet%xup = 0
        jet%bp = 0
        jet%bup = 0
        do m = 0, self%mpol - 1
            do n = 0, self%ntor
                do j = 1, self%ns
                    k = 1 + 2*(n + (self%ntor + 1)*m) + wout%mnmax*(j - 1)
                    do block_ = 1, 4
                        col = 1 + n + (self%ntor + 1)*(m + self%mpol*(j - 1 + self%ns*(blocks(block_) - 1)))
                        do a = 1, 2
                            do b = 1, 2
                                offset = 0
                                if (b == 2) offset = nr
                                if (b == 1) then
                                    call accumulate(k + a - 1 + offset, col, rr(a, &
                                        block_))
                                else
                                    call accumulate(k + a - 1 + offset, col, zr(a, &
                                        block_))
                                end if
                            end do
                        end do
                    end do
                end do
            end do
        end do
        jet%bp(:, ng + nb + 1:) = local%bp(:, 2*nr + nb + 1:)
        jet%bup(:, :, ng + nb + 1:) = local%bup(:, :, 2*nr + nb + 1:)
    contains
        subroutine accumulate(source, dest, weight)
            integer, intent(in)::source, dest
            real(dp), intent(in)::weight
            jet%xp(:, dest) = jet%xp(:, dest) + weight*local%xp(:, source)
            jet%xup(:, :, dest) = jet%xup(:, :, dest) + weight*local%xup(:, :, source)
            jet%bp(:, dest) = jet%bp(:, dest) + weight*local%bp(:, source)
            jet%bup(:, :, dest) = jet%bup(:, :, dest) + weight*local%bup(:, :, source)
        end subroutine
    end subroutine state_vector_jets

end module tiago_vmecpp_state
