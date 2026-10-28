program test_mse_state
    !! MSE angles from the interior VMEC field (tiago_mse_state).
    !! Usage: test_mse_state <wout.nc>
    !! 1. Midpoints at X(u) on half-grid surfaces: gamma from the state equals
    !!    gamma from B = B^u dX/du + B^v dX/dv evaluated here directly from the
    !!    wout Fourier coefficients (bsupumnc, bsupvmnc on the half grid; R, Z
    !!    and their angle derivatives from rmnc, zmns averaged to the half
    !!    grid). The two differ by the radial interpolation of R and Z only.
    !! 2. d gamma / d p against central differences of the state parameters.
    !! 3. Midpoints outside the boundary and states without B are errors.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit, i32 => int32
    use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_get_var, nf90_nowrite
    use tiago_diagnostic_types, only: mse_channel_t
    use tiago_mse, only: mse_angles
    use tiago_vmec_state, only: vmec_state_t
    use tiago_mse_state, only: mse_state_angles
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: nch = 6
    character(len=1024) :: wout
    type(vmec_state_t) :: state, scalar_only
    type(mse_channel_t) :: ch(nch), outside(1)
    real(dp), allocatable :: rmnc(:, :), zmns(:, :), bu(:, :), bv(:, :), xm(:), xn(:), &
        xmn(:), xnn(:), gamma(:), ref(:), jac(:, :), p0(:), v(:), gp(:), gm(:), field(:, :)
    real(dp) :: u(3, nch), x(3), dxdu(3, 3), b, db(3), bref(3, nch), h, fd, jv, scale
    integer :: failures, j, k, ns, ierr
    integer(i32) :: code
    character(len=:), allocatable :: message

    failures = 0
    call get_command_argument(1, wout)
    call state%load(trim(wout), ierr, message, vector=.true.)
    if (ierr /= 0) error stop 'cannot load wout'
    call read_wout()

    do k = 1, nch
        j = 4 + 2 * k
        u(:, k) = [state%rho_half(j - 1), 0.9_dp * k, 0.23_dp * k]
        call state%eval(u(:, k), x, dxdu, b, db)
        ch(k) = channel(x, 0.31_dp * k)
        bref(:, k) = direct_field(j, u(2, k), u(3, k))
    end do

    ! 1. State field against the direct wout evaluation.
    call mse_state_angles(state, ch, gamma, code, message, field=field)
    call check(code == 0, 'mse_state_angles: '//message)
    call mse_angles(ch, bref, ref, code, message)
    call check(code == 0, 'reference angles: '//message)
    do k = 1, nch
        call close('|B| state vs wout', norm2(field(:, k)), norm2(bref(:, k)), 2.0e-3_dp)
        call close_abs('gamma state vs wout [rad]', gamma(k), ref(k), 2.0e-3_dp)
    end do

    ! 2. Parameter derivative at fixed midpoints.
    call mse_state_angles(state, ch, gamma, code, message, param_jacobian=jac)
    allocate(p0(state%n_params()))
    call state%get_params(p0)
    v = [(sin(1.3_dp * k + 0.4_dp), k = 1, size(p0))] * max(abs(p0), 1.0e-2_dp * maxval(abs(p0)))
    h = 1.0e-6_dp
    call state%set_params(p0 + h * v)
    call mse_state_angles(state, ch, gp, code, message)
    call state%set_params(p0 - h * v)
    call mse_state_angles(state, ch, gm, code, message)
    call state%set_params(p0)
    scale = 0.0_dp
    do k = 1, nch
        scale = max(scale, abs(dot_product(jac(k, :), v)))
    end do
    do k = 1, nch
        fd = (gp(k) - gm(k)) / (2.0_dp * h)
        jv = dot_product(jac(k, :), v)
        call close_abs('d gamma / dp . v vs central difference', jv, fd, 1.0e-6_dp * scale)
    end do

    ! 3. Errors instead of angles.
    outside(1) = channel([2.5_dp, 0.0_dp, 0.0_dp], 0.0_dp)
    call mse_state_angles(state, outside, gamma, code, message)
    call check(code /= 0 .and. index(message, 'outside the plasma') > 0, &
        'midpoint outside the plasma is rejected')
    call scalar_only%load(trim(wout), ierr, message)
    call mse_state_angles(scalar_only, ch(1:1), gamma, code, message)
    call check(code /= 0 .and. index(message, 'no vector field') > 0, &
        'a state without B is rejected')

    if (failures > 0) then
        write(error_unit, '(I0,A)') failures, ' check(s) failed'
        error stop 1
    end if
    print '(A)', 'test_mse_state passed'

contains

    type(mse_channel_t) function channel(point, angle) result(c)
        !! Beam along a tilted horizontal direction through point, sightline
        !! from a lens above and outboard; both lines meet at point.
        real(dp), intent(in) :: point(3), angle
        real(dp) :: v(3), n(3)
        v = [cos(angle + 2.0_dp), sin(angle + 2.0_dp), 0.1_dp]
        n = [cos(angle), sin(angle), 0.3_dp]
        c%label = 'synthetic'
        c%beam_point = point - 0.7_dp * v
        c%beam_direction = v
        c%sight_point = point + 0.9_dp * n
        c%sight_direction = n
        c%max_gap = 1.0e-9_dp
    end function channel

    function direct_field(j, theta, phi) result(bc)
        !! Cartesian B on half-grid surface j (wout index) at (theta, phi).
        integer, intent(in) :: j
        real(dp), intent(in) :: theta, phi
        real(dp) :: bc(3), r, rt, rp, zt, zp, but, bvt, a, c, s, rh(size(xm)), zh(size(xm))
        integer :: mn
        rh = 0.5_dp * (rmnc(:, j - 1) + rmnc(:, j))
        zh = 0.5_dp * (zmns(:, j - 1) + zmns(:, j))
        r = 0.0_dp
        rt = 0.0_dp
        rp = 0.0_dp
        zt = 0.0_dp
        zp = 0.0_dp
        do mn = 1, size(xm)
            a = xm(mn) * theta - xn(mn) * phi
            c = cos(a)
            s = sin(a)
            r = r + rh(mn) * c
            rt = rt - xm(mn) * rh(mn) * s
            rp = rp + xn(mn) * rh(mn) * s
            zt = zt + xm(mn) * zh(mn) * c
            zp = zp - xn(mn) * zh(mn) * c
        end do
        but = 0.0_dp
        bvt = 0.0_dp
        do mn = 1, size(xmn)
            a = xmn(mn) * theta - xnn(mn) * phi
            but = but + bu(mn, j) * cos(a)
            bvt = bvt + bv(mn, j) * cos(a)
        end do
        ! B = B^theta dX/dtheta + B^phi dX/dphi in (R, phi, Z), then Cartesian.
        bc = [(but * rt + bvt * rp) * cos(phi) - bvt * r * sin(phi), &
            (but * rt + bvt * rp) * sin(phi) + bvt * r * cos(phi), &
            but * zt + bvt * zp]
    end function direct_field

    subroutine read_wout()
        integer :: ncid, varid, s, mnmax, mnmax_nyq
        s = nf90_open(trim(wout), nf90_nowrite, ncid)
        s = nf90_inq_varid(ncid, 'ns', varid)
        s = nf90_get_var(ncid, varid, ns)
        s = nf90_inq_varid(ncid, 'mnmax', varid)
        s = nf90_get_var(ncid, varid, mnmax)
        s = nf90_inq_varid(ncid, 'mnmax_nyq', varid)
        s = nf90_get_var(ncid, varid, mnmax_nyq)
        allocate(xm(mnmax), xn(mnmax), xmn(mnmax_nyq), xnn(mnmax_nyq), rmnc(mnmax, ns), &
            zmns(mnmax, ns), bu(mnmax_nyq, ns), bv(mnmax_nyq, ns))
        s = nf90_inq_varid(ncid, 'xm', varid)
        s = nf90_get_var(ncid, varid, xm)
        s = nf90_inq_varid(ncid, 'xn', varid)
        s = nf90_get_var(ncid, varid, xn)
        s = nf90_inq_varid(ncid, 'xm_nyq', varid)
        s = nf90_get_var(ncid, varid, xmn)
        s = nf90_inq_varid(ncid, 'xn_nyq', varid)
        s = nf90_get_var(ncid, varid, xnn)
        s = nf90_inq_varid(ncid, 'rmnc', varid)
        s = nf90_get_var(ncid, varid, rmnc)
        s = nf90_inq_varid(ncid, 'zmns', varid)
        s = nf90_get_var(ncid, varid, zmns)
        s = nf90_inq_varid(ncid, 'bsupumnc', varid)
        s = nf90_get_var(ncid, varid, bu)
        s = nf90_inq_varid(ncid, 'bsupvmnc', varid)
        s = nf90_get_var(ncid, varid, bv)
        s = nf90_close(ncid)
    end subroutine read_wout

    subroutine check(ok, name)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: name
        if (ok) then
            print '(A)', 'ok   '//name
        else
            write(error_unit, '(A)') 'FAIL '//name
            failures = failures + 1
        end if
    end subroutine check

    subroutine close(name, a, b, rtol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: a, b, rtol
        call close_abs(name, a, b, rtol * abs(b))
    end subroutine close

    subroutine close_abs(name, a, b, tol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: a, b, tol
        if (abs(a - b) > tol .or. a /= a) then
            write(error_unit, '(A,3ES14.6)') 'FAIL '//name//': ', a, b, a - b
            failures = failures + 1
        else
            print '(A,3ES14.6)', 'ok   '//name//': ', a, b, a - b
        end if
    end subroutine close_abs
end program test_mse_state
