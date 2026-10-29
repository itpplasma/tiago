module tiago_vmec_state
    !! VMEC nested-surface state as a real_space_state_t (labels rho = sqrt(s),
    !! VMEC theta, cylindrical phi = VMEC zeta).
    !!
    !! Every Fourier coefficient of R, Z (full radial grid) and |B| (half grid)
    !! is a C1 cubic Hermite interpolant in rho of its radial node values
    !! (tiago_radial_interp); the end cubics are extended beyond the last node,
    !! so the chart continues smoothly a little outside rho = 1. The parameters
    !! are these node values,
    !!     p = (rmnc(:, 1:ns), zmns(:, 1:ns), bmnc(:, 2:ns)),
    !! on which the map depends linearly at fixed labels. A VMEC adjoint maps a
    !! cotangent of p to its own inputs. Stellarator-symmetric (lasym = F) only.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_get_var, &
        nf90_noerr, nf90_nowrite
    use tiago_real_space_state, only: real_space_state_t, vector_jets_t, locate_ok, locate_singular
    use tiago_radial_interp, only: hermite_weights, max_weights
    use tiago_vmec_r_x, only: vmec_r_x
    use tiago_vmec_r_y, only: vmec_r_y
    use tiago_vmec_z, only: vmec_z
    use tiago_vmec_cos, only: vmec_cos
    use tiago_vmec_vector_component, only: vmec_vector_component
    implicit none
    private

    type, extends(real_space_state_t), public :: vmec_state_t
        integer :: ns = 0, mnmax = 0, mnmax_nyq = 0, nfp = 1
        real(dp), allocatable :: xm(:), xn(:), xm_nyq(:), xn_nyq(:)
        real(dp), allocatable :: rho_full(:), rho_half(:)
        real(dp), allocatable :: rmnc(:, :), zmns(:, :), bmnc(:, :) !! bmnc(:, 1:ns-1) = half grid
        real(dp), allocatable :: bsupumnc(:,:), bsupvmnc(:,:)
    contains
        procedure :: vector_jets => vmec_vector_jets
        procedure :: load => vmec_load
        procedure :: n_params => vmec_n_params
        procedure :: get_params => vmec_get_params
        procedure :: set_params => vmec_set_params
        procedure :: eval => vmec_eval
        procedure :: initial_guess => vmec_initial_guess
        procedure :: param_jvp => vmec_param_jvp
        procedure :: param_vjp => vmec_param_vjp
    end type vmec_state_t

contains

    subroutine vmec_load(self, path, ierr, message, vector)
        class(vmec_state_t), intent(inout) :: self
        character(len=*), intent(in) :: path
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        logical, intent(in), optional :: vector
        integer :: ncid, lasym, j
        real(dp), allocatable :: buf(:, :)

        ierr = 0
        message = ''
        if (nf90_open(trim(path), nf90_nowrite, ncid) /= nf90_noerr) then
            ierr = 1
            message = 'cannot open VMEC wout '//trim(path)
            return
        end if
        call get_int('ns', self%ns)
        call get_int('mnmax', self%mnmax)
        call get_int('mnmax_nyq', self%mnmax_nyq)
        call get_int('nfp', self%nfp)
        call get_int('lasym__logical__', lasym)
        if (ierr /= 0) return
        if (lasym /= 0) then
            ierr = 2
            message = 'vmec_state_t supports stellarator-symmetric equilibria only'
            return
        end if
        allocate(self%xm(self%mnmax), self%xn(self%mnmax), self%xm_nyq(self%mnmax_nyq), &
            self%xn_nyq(self%mnmax_nyq), self%rmnc(self%mnmax, self%ns), &
            self%zmns(self%mnmax, self%ns), buf(self%mnmax_nyq, self%ns))
        call get_var1('xm', self%xm)
        call get_var1('xn', self%xn)
        call get_var1('xm_nyq', self%xm_nyq)
        call get_var1('xn_nyq', self%xn_nyq)
        call get_var2('rmnc', self%rmnc)
        call get_var2('zmns', self%zmns)
        call get_var2('bmnc', buf)
        if (present(vector)) then
            if (vector) then
                self%bmnc = buf(:, 2:self%ns)
                call get_var2('bsupumnc', buf)
                self%bsupumnc = buf(:, 2:self%ns)
                call get_var2('bsupvmnc', buf)
                self%bsupvmnc = buf(:, 2:self%ns)
            end if
        end if
        if (nf90_close(ncid) /= nf90_noerr .and. ierr == 0) ierr = 1
        if (ierr /= 0) return
        if (.not. allocated(self%bmnc)) self%bmnc = buf(:, 2:self%ns)
        self%rho_full = [(sqrt(real(j - 1, dp) / real(self%ns - 1, dp)), j = 1, self%ns)]
        self%rho_half = [(sqrt((real(j, dp) - 1.5_dp) / real(self%ns - 1, dp)), &
            j = 2, self%ns)]
    contains
        subroutine get_int(name, value)
            character(len=*), intent(in) :: name
            integer, intent(out) :: value
            integer :: varid
            value = 0
            if (ierr /= 0) return
            if (nf90_inq_varid(ncid, name, varid) /= nf90_noerr) then
                ierr = 3
                message = 'wout variable missing: '//name
            else if (nf90_get_var(ncid, varid, value) /= nf90_noerr) then
                ierr = 3
                message = 'cannot read wout variable '//name
            end if
        end subroutine get_int

        subroutine get_var1(name, values)
            character(len=*), intent(in) :: name
            real(dp), intent(out) :: values(:)
            integer :: varid
            if (ierr /= 0) return
            if (nf90_inq_varid(ncid, name, varid) /= nf90_noerr) then
                ierr = 3
                message = 'wout variable missing: '//name
            else if (nf90_get_var(ncid, varid, values) /= nf90_noerr) then
                ierr = 3
                message = 'cannot read wout variable '//name
            end if
        end subroutine get_var1

        subroutine get_var2(name, values)
            character(len=*), intent(in) :: name
            real(dp), intent(out) :: values(:, :)
            integer :: varid
            if (ierr /= 0) return
            if (nf90_inq_varid(ncid, name, varid) /= nf90_noerr) then
                ierr = 3
                message = 'wout variable missing: '//name
            else if (nf90_get_var(ncid, varid, values) /= nf90_noerr) then
                ierr = 3
                message = 'cannot read wout variable '//name
            end if
        end subroutine get_var2
    end subroutine vmec_load

    pure integer function vmec_n_params(self)
        class(vmec_state_t), intent(in) :: self
        vmec_n_params = 2 * self%mnmax * self%ns + self%mnmax_nyq * (self%ns - 1)
        if (allocated(self%bsupumnc)) vmec_n_params = vmec_n_params + 2*size(self%bsupumnc)
    end function vmec_n_params

    pure subroutine vmec_get_params(self, p)
        class(vmec_state_t), intent(in) :: self
        real(dp), intent(out) :: p(:)
        integer :: nr, nb

        nr = self%mnmax * self%ns
        nb = self%mnmax_nyq * (self%ns - 1)
        p(1:nr) = reshape(self%rmnc, [nr])
        p(nr + 1:2 * nr) = reshape(self%zmns, [nr])
        p(2 * nr + 1:2 * nr + nb) = reshape(self%bmnc, [nb])
        if (allocated(self%bsupumnc)) then
            p(2*nr+nb+1:2*nr+2*nb) = reshape(self%bsupumnc, [nb])
            p(2*nr+2*nb+1:2*nr+3*nb) = reshape(self%bsupvmnc, [nb])
        end if
    end subroutine vmec_get_params

    subroutine vmec_set_params(self, p)
        class(vmec_state_t), intent(inout) :: self
        real(dp), intent(in) :: p(:)
        integer :: nr, nb

        nr = self%mnmax * self%ns
        nb = self%mnmax_nyq * (self%ns - 1)
        self%rmnc = reshape(p(1:nr), [self%mnmax, self%ns])
        self%zmns = reshape(p(nr + 1:2 * nr), [self%mnmax, self%ns])
        self%bmnc = reshape(p(2 * nr + 1:2 * nr + nb), [self%mnmax_nyq, self%ns - 1])
        if (allocated(self%bsupumnc)) then
            self%bsupumnc = reshape(p(2*nr+nb+1:2*nr+2*nb), shape(self%bmnc))
            self%bsupvmnc = reshape(p(2*nr+2*nb+1:2*nr+3*nb), shape(self%bmnc))
        end if
    end subroutine vmec_set_params

    pure subroutine vmec_eval(self, u, x, dxdu, b, dbdu)
        class(vmec_state_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3), dxdu(3, 3), b, dbdu(3)
        real(dp) :: wf(max_weights), dwf(max_weights), wh(max_weights), dwh(max_weights)
        integer :: ff, nf, fh, nh, mn
        real(dp) :: r, dr(3), z, dz(3), c, s, cr, dcr, cz, dcz, ang, cp, sp

        call hermite_weights(self%rho_full, u(1), .true., ff, nf, wf, dwf)
        call hermite_weights(self%rho_half, u(1), .true., fh, nh, wh, dwh)
        r = 0.0_dp
        z = 0.0_dp
        dr = 0.0_dp
        dz = 0.0_dp
        do mn = 1, self%mnmax
            ang = self%xm(mn) * u(2) - self%xn(mn) * u(3)
            c = cos(ang)
            s = sin(ang)
            cr = dot_product(wf(:nf), self%rmnc(mn, ff:ff + nf - 1))
            dcr = dot_product(dwf(:nf), self%rmnc(mn, ff:ff + nf - 1))
            cz = dot_product(wf(:nf), self%zmns(mn, ff:ff + nf - 1))
            dcz = dot_product(dwf(:nf), self%zmns(mn, ff:ff + nf - 1))
            r = r + cr * c
            z = z + cz * s
            dr = dr + [dcr * c, -self%xm(mn) * cr * s, self%xn(mn) * cr * s]
            dz = dz + [dcz * s, self%xm(mn) * cz * c, -self%xn(mn) * cz * c]
        end do
        b = 0.0_dp
        dbdu = 0.0_dp
        do mn = 1, self%mnmax_nyq
            ang = self%xm_nyq(mn) * u(2) - self%xn_nyq(mn) * u(3)
            c = cos(ang)
            s = sin(ang)
            cr = dot_product(wh(:nh), self%bmnc(mn, fh:fh + nh - 1))
            dcr = dot_product(dwh(:nh), self%bmnc(mn, fh:fh + nh - 1))
            b = b + cr * c
            dbdu = dbdu + [dcr * c, -self%xm_nyq(mn) * cr * s, self%xn_nyq(mn) * cr * s]
        end do
        cp = cos(u(3))
        sp = sin(u(3))
        x = [r * cp, r * sp, z]
        dxdu(1, :) = dr * cp
        dxdu(2, :) = dr * sp
        dxdu(3, :) = dz
        dxdu(1, 3) = dxdu(1, 3) - r * sp
        dxdu(2, 3) = dxdu(2, 3) + r * cp
    end subroutine vmec_eval

    pure subroutine vmec_initial_guess(self, x, u)
        !! The grid search of default_initial_guess, with the radial
        !! interpolation done once per grid radius instead of per point.
        class(vmec_state_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: u(3)
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: wf(max_weights), dwf(max_weights), cr(self%mnmax), &
            cz(self%mnmax), v(3), xv(3), r, z, ang, best, d
        integer :: ff, nf, i, j, mn

        best = huge(1.0_dp)
        u = [0.5_dp, 0.0_dp, atan2(x(2), x(1))]
        do i = 1, 4
            call hermite_weights(self%rho_full, 0.25_dp * i, .true., ff, nf, wf, dwf)
            do mn = 1, self%mnmax
                cr(mn) = dot_product(wf(:nf), self%rmnc(mn, ff:ff + nf - 1))
                cz(mn) = dot_product(wf(:nf), self%zmns(mn, ff:ff + nf - 1))
            end do
            do j = 0, 15
                v = [0.25_dp * i, 2.0_dp * pi * j / 16.0_dp, atan2(x(2), x(1))]
                r = 0.0_dp
                z = 0.0_dp
                do mn = 1, self%mnmax
                    ang = self%xm(mn) * v(2) - self%xn(mn) * v(3)
                    r = r + cr(mn) * cos(ang)
                    z = z + cz(mn) * sin(ang)
                end do
                xv = [r * cos(v(3)), r * sin(v(3)), z]
                d = sum((xv - x)**2)
                if (d < best) then
                    best = d
                    u = v
                end if
            end do
        end do
    end subroutine vmec_initial_guess

    pure subroutine vmec_param_jvp(self, u, dp_, dx, db)
        class(vmec_state_t), intent(in) :: self
        real(dp), intent(in) :: u(3), dp_(:)
        real(dp), intent(out) :: dx(3), db
        real(dp) :: wf(max_weights), dwf(max_weights), wh(max_weights), dwh(max_weights)
        integer :: ff, nf, fh, nh, mn, k, nr
        real(dp) :: dr, dz, ang

        nr = self%mnmax * self%ns
        call hermite_weights(self%rho_full, u(1), .true., ff, nf, wf, dwf)
        call hermite_weights(self%rho_half, u(1), .true., fh, nh, wh, dwh)
        dr = 0.0_dp
        dz = 0.0_dp
        do mn = 1, self%mnmax
            ang = self%xm(mn) * u(2) - self%xn(mn) * u(3)
            do k = 1, nf
                dr = dr + wf(k) * cos(ang) * dp_(mn + (ff + k - 2) * self%mnmax)
                dz = dz + wf(k) * sin(ang) * dp_(nr + mn + (ff + k - 2) * self%mnmax)
            end do
        end do
        db = 0.0_dp
        do mn = 1, self%mnmax_nyq
            ang = self%xm_nyq(mn) * u(2) - self%xn_nyq(mn) * u(3)
            do k = 1, nh
                db = db + wh(k) * cos(ang) * dp_(2 * nr + mn + (fh + k - 2) * self%mnmax_nyq)
            end do
        end do
        dx = [dr * cos(u(3)), dr * sin(u(3)), dz]
    end subroutine vmec_param_jvp

    pure subroutine vmec_param_vjp(self, u, xbar, bbar, pbar)
        class(vmec_state_t), intent(in) :: self
        real(dp), intent(in) :: u(3), xbar(3), bbar
        real(dp), intent(inout) :: pbar(:)
        real(dp) :: wf(max_weights), dwf(max_weights), wh(max_weights), dwh(max_weights)
        integer :: ff, nf, fh, nh, mn, k, nr, i
        real(dp) :: rbar, ang

        nr = self%mnmax * self%ns
        call hermite_weights(self%rho_full, u(1), .true., ff, nf, wf, dwf)
        call hermite_weights(self%rho_half, u(1), .true., fh, nh, wh, dwh)
        rbar = xbar(1) * cos(u(3)) + xbar(2) * sin(u(3))
        do mn = 1, self%mnmax
            ang = self%xm(mn) * u(2) - self%xn(mn) * u(3)
            do k = 1, nf
                i = mn + (ff + k - 2) * self%mnmax
                pbar(i) = pbar(i) + wf(k) * cos(ang) * rbar
                pbar(nr + i) = pbar(nr + i) + wf(k) * sin(ang) * xbar(3)
            end do
        end do
        do mn = 1, self%mnmax_nyq
            ang = self%xm_nyq(mn) * u(2) - self%xn_nyq(mn) * u(3)
            do k = 1, nh
                i = 2 * nr + mn + (fh + k - 2) * self%mnmax_nyq
                pbar(i) = pbar(i) + wh(k) * cos(ang) * bbar
            end do
        end do
    end subroutine vmec_param_vjp

    subroutine vmec_vector_jets(self, u, jet, status, parameters)
        class(vmec_state_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        type(vector_jets_t), intent(out) :: jet
        integer, intent(out) :: status
        logical, intent(in), optional :: parameters
        real(dp) :: wf(max_weights, 4), wh(max_weights, 4), radial(4)
        real(dp) :: val, d(3), d2(3, 3), d3(3, 3, 3), x3(3, 3, 3, 3)
        real(dp) :: v(2), vu(3, 2), vuu(3, 3, 2), vp(2), vup(3, 2)
        real(dp) :: junk2(3, 3), zero3(3, 3, 3), zero2(3, 3, 2)
        real(dp), allocatable :: x2p(:, :, :, :), cp(:, :), cup(:, :, :)
        integer :: ff, nf, fh, nh, mn, k, i, j, nr, nb, np, col
        real(dp) :: basis
        logical :: columns
        status = locate_singular
        if (.not. allocated(self%bsupumnc)) return
        if (.not. allocated(self%bsupvmnc)) return
        columns = .true.
        if (present(parameters)) columns = parameters
        np = 0
        if (columns) np = self%n_params()
        nr = self%mnmax*self%ns
        nb = self%mnmax_nyq*(self%ns - 1)
        allocate (jet%xp(3, np), jet%xup(3, 3, np), jet%bp(3, np), jet%bup(3, 3, np), &
                  x2p(3, 3, 3, np), cp(2, np), cup(3, 2, np))
        jet%x = 0
        jet%xu = 0
        jet%xuu = 0
        jet%xp = 0
        jet%xup = 0
        jet%bp = 0
        jet%bup = 0
        x2p = 0
        x3 = 0
        v = 0
        vu = 0
        vuu = 0
        cp = 0
        cup = 0
        zero3 = 0
        zero2 = 0
        call hermite_weights(self%rho_full, u(1), .true., ff, nf, &
                             wf(:, 1), wf(:, 2), wf(:, 3), wf(:, 4))
        call hermite_weights(self%rho_half, u(1), .true., fh, nh, &
                             wh(:, 1), wh(:, 2), wh(:, 3), wh(:, 4))
        do mn = 1, self%mnmax
            do k = 1, nf
                col = mn + (ff + k - 2)*self%mnmax
                do i = 1, 3
                    select case (i)
                    case (1)
                        call vmec_r_x(wf(k, :), self%xm(mn), self%xn(mn), u(2), u(3), &
                                      val, d, d2, d3)
                    case (2)
                        call vmec_r_y(wf(k, :), self%xm(mn), self%xn(mn), u(2), u(3), &
                                      val, d, d2, d3)
                    case (3)
                        call vmec_z(wf(k, :), self%xm(mn), self%xn(mn), u(2), u(3), &
                                    val, d, d2, d3)
                    end select
                    basis = val
                    if (columns) then
                        j = col
                        if (i == 3) j = nr + col
                        jet%xp(i, j) = val
                        jet%xup(i, :, j) = d
                        x2p(i, :, :, j) = d2
                    end if
                    if (i == 3) then
                        val = self%zmns(mn, ff + k - 1)
                    else
                        val = self%rmnc(mn, ff + k - 1)
                    end if
                    jet%x(i) = jet%x(i) + basis*val
                    jet%xu(i, :) = jet%xu(i, :) + d*val
                    jet%xuu(i, :, :) = jet%xuu(i, :, :) + d2*val
                    x3(i, :, :, :) = x3(i, :, :, :) + d3*val
                end do
            end do
        end do
        do mn = 1, self%mnmax_nyq
            do k = 1, nh
                call vmec_cos(wh(k, :), self%xm_nyq(mn), self%xn_nyq(mn), &
                              u(2), u(3), val, d, d2, d3)
                do i = 1, 2
                    if (columns) then
                        col = 2*nr + i*nb + mn + (fh + k - 2)*self%mnmax_nyq
                        cp(i, col) = val
                        cup(:, i, col) = d
                    end if
                    if (i == 1) then
                        radial(1) = self%bsupumnc(mn, fh + k - 1)
                    else
                        radial(1) = self%bsupvmnc(mn, fh + k - 1)
                    end if
                    v(i) = v(i) + radial(1)*val
                    vu(:, i) = vu(:, i) + radial(1)*d
                    vuu(:, :, i) = vuu(:, :, i) + radial(1)*d2
                end do
            end do
        end do
        do i = 1, 3
            call vmec_vector_component(jet%xu(i, :), jet%xuu(i, :, :), x3(i, :, :, :), &
                                       v, vu, vuu, jet%b(i), jet%bu(i, :), jet%buu(i, &
                                                                                   :, &
                                                                                       :))
            do j = 1, np
                if (any(jet%xp(:, j) /= 0) .or. any(jet%xup(:, :, j) /= 0)) then
                    call vmec_vector_component(jet%xup(i, :, j), x2p(i, :, :, j), &
                                               zero3, &
                                               v, vu, zero2, jet%bp(i, j), jet%bup(i, &
                                                                                   :, &
                                                                                       j), junk2)
                end if
                if (any(cp(:, j) /= 0) .or. any(cup(:, :, j) /= 0)) then
                    vp = cp(:, j)
                    vup = cup(:, :, j)
                    call vmec_vector_component(jet%xu(i, :), jet%xuu(i, :, :), zero3, &
                                               vp, vup, zero2, val, d, junk2)
                    jet%bp(i, j) = jet%bp(i, j) + val
                    jet%bup(i, :, j) = jet%bup(i, :, j) + d
                end if
            end do
        end do
        status = locate_ok
    end subroutine vmec_vector_jets

end module tiago_vmec_state
