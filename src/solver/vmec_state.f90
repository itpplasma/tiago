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
    use tiago_real_space_state, only: real_space_state_t
    use tiago_radial_interp, only: hermite_weights, max_weights
    implicit none
    private

    type, extends(real_space_state_t), public :: vmec_state_t
        integer :: ns = 0, mnmax = 0, mnmax_nyq = 0, nfp = 1
        real(dp), allocatable :: xm(:), xn(:), xm_nyq(:), xn_nyq(:)
        real(dp), allocatable :: rho_full(:), rho_half(:)
        real(dp), allocatable :: rmnc(:, :), zmns(:, :), bmnc(:, :) !! bmnc(:, 1:ns-1) = half grid
    contains
        procedure :: load => vmec_load
        procedure :: n_params => vmec_n_params
        procedure :: get_params => vmec_get_params
        procedure :: set_params => vmec_set_params
        procedure :: eval => vmec_eval
        procedure :: param_jvp => vmec_param_jvp
        procedure :: param_vjp => vmec_param_vjp
    end type vmec_state_t

contains

    subroutine vmec_load(self, path, ierr, message)
        class(vmec_state_t), intent(inout) :: self
        character(len=*), intent(in) :: path
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
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
        if (nf90_close(ncid) /= nf90_noerr .and. ierr == 0) ierr = 1
        if (ierr /= 0) return
        self%bmnc = buf(:, 2:self%ns)
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
end module tiago_vmec_state
