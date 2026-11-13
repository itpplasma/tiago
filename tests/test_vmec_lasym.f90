program test_vmec_lasym
    !! Non-symmetric VMEC states (lasym = T) in vmec_state_t.
    !! Usage: test_vmec_lasym <symmetric wout.nc> <rotated wout to write>
    !! The symmetric equilibrium turned rigidly by phi0 about the Z axis is a
    !! non-symmetric wout: with A = m theta - n zeta,
    !!     c cos(A + n phi0) = c cos(n phi0) cos A - c sin(n phi0) sin A,
    !!     s sin(A + n phi0) = s cos(n phi0) sin A + s sin(n phi0) cos A,
    !! so every coefficient acquires a sine (R, |B|, B^theta, B^phi) or cosine
    !! (Z) partner. Position, |B|, the vector field and the MSE angle of the
    !! rotated state at rotated points must equal the rotated values of the
    !! symmetric state (the MSE reference +Z is unchanged by the rotation).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit, i32 => int32
    use netcdf
    use tiago_diagnostic_types, only: mse_channel_t
    use tiago_vmec_state, only: vmec_state_t
    use tiago_mse_state, only: mse_state_angles
    implicit none

    real(dp), parameter :: phi0 = 0.3_dp
    integer, parameter :: nch = 6
    character(len=1024) :: wout, rotated
    type(vmec_state_t) :: sym, asym
    type(mse_channel_t) :: ch(nch), chr(nch)
    real(dp) :: rot(3, 3), u(3), x(3), dxdu(3, 3), b, db(3), xr(3), dxdur(3, 3), br, dbr(3)
    real(dp), allocatable :: gamma(:), gammar(:), field(:, :), fieldr(:, :), jac(:, :), &
        p0(:), v(:), gp(:), gm(:)
    real(dp) :: h, scale, fd, jv
    integer :: failures, k, ierr
    integer(i32) :: code
    character(len=:), allocatable :: message

    failures = 0
    call get_command_argument(1, wout)
    call get_command_argument(2, rotated)
    call write_rotated(trim(wout), trim(rotated))
    call sym%load(trim(wout), ierr, message, vector=.true.)
    if (ierr /= 0) error stop 'cannot load symmetric wout'
    call asym%load(trim(rotated), ierr, message, vector=.true.)
    call check(ierr == 0, 'rotated lasym wout loads: '//message)
    if (ierr /= 0) error stop 1

    rot = reshape([cos(phi0), sin(phi0), 0.0_dp, -sin(phi0), cos(phi0), 0.0_dp, &
                   0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
    do k = 1, nch
        u = [sym%rho_half(3 + 2*k), 0.9_dp*k, 0.23_dp*k]
        call sym%eval(u, x, dxdu, b, db)
        call asym%eval(u + [0.0_dp, 0.0_dp, phi0], xr, dxdur, br, dbr)
        call close_abs('rotated position [m]', norm2(xr - matmul(rot, x)), 0.0_dp, 1.0e-12_dp)
        call close_abs('rotated d position / d label [m]', &
                       maxval(abs(dxdur - matmul(rot, dxdu))), 0.0_dp, 1.0e-11_dp)
        call close_abs('|B| [T]', br, b, 1.0e-12_dp*abs(b))
        call close_abs('d|B| / d label [T]', maxval(abs(dbr - db)), 0.0_dp, 1.0e-11_dp)
        ch(k) = channel(x, 0.31_dp*k)
        chr(k) = ch(k)
        chr(k)%beam_point = matmul(rot, ch(k)%beam_point)
        chr(k)%beam_direction = matmul(rot, ch(k)%beam_direction)
        chr(k)%sight_point = matmul(rot, ch(k)%sight_point)
        chr(k)%sight_direction = matmul(rot, ch(k)%sight_direction)
    end do

    call mse_state_angles(sym, ch, gamma, code, message, field=field)
    call check(code == 0, 'symmetric angles: '//message)
    call mse_state_angles(asym, chr, gammar, code, message, field=fieldr)
    call check(code == 0, 'rotated angles: '//message)
    if (failures > 0) error stop 1
    do k = 1, nch
        call close_abs('rotated B [T]', norm2(fieldr(:, k) - matmul(rot, field(:, k))), &
                       0.0_dp, 1.0e-10_dp*norm2(field(:, k)))
        call close_abs('gamma [rad]', gammar(k), gamma(k), 1.0e-10_dp)
    end do

    ! The symmetric coefficients stay the parameters of a lasym state.
    call mse_state_angles(asym, chr, gammar, code, message, param_jacobian=jac)
    call check(code == 0, 'lasym parameter Jacobian: '//message)
    allocate (p0(asym%n_params()))
    call asym%get_params(p0)
    v = [(sin(1.3_dp*k + 0.4_dp), k=1, size(p0))]*max(abs(p0), 1.0e-2_dp*maxval(abs(p0)))
    h = 1.0e-6_dp
    call asym%set_params(p0 + h*v)
    call mse_state_angles(asym, chr, gp, code, message)
    call asym%set_params(p0 - h*v)
    call mse_state_angles(asym, chr, gm, code, message)
    call asym%set_params(p0)
    scale = maxval(abs(matmul(jac, v)))
    do k = 1, nch
        fd = (gp(k) - gm(k))/(2.0_dp*h)
        jv = dot_product(jac(k, :), v)
        call close_abs('lasym d gamma / dp . v vs central difference', jv, fd, 1.0e-6_dp*scale)
    end do

    if (failures > 0) then
        write (error_unit, '(I0,A)') failures, ' check(s) failed'
        error stop 1
    end if
    print '(A)', 'test_vmec_lasym passed'

contains

    type(mse_channel_t) function channel(point, angle) result(c)
        real(dp), intent(in) :: point(3), angle
        real(dp) :: v(3), n(3)
        v = [cos(angle + 2.0_dp), sin(angle + 2.0_dp), 0.1_dp]
        n = [cos(angle), sin(angle), 0.3_dp]
        c%label = 'synthetic'
        c%beam_point = point - 0.7_dp*v
        c%beam_direction = v
        c%sight_point = point + 0.9_dp*n
        c%sight_direction = n
        c%max_gap = 1.0e-9_dp
    end function channel

    subroutine write_rotated(source, target)
        character(len=*), intent(in) :: source, target
        integer :: ncid, out, ns, mnmax, nyq, nfp, d_r, d_mn, d_nyq, var, j
        real(dp), allocatable :: xm(:), xn(:), xmn(:), xnn(:), cs(:), sn(:), cn(:), snn(:)
        real(dp), allocatable :: a(:, :)
        character(len=8), parameter :: cosine(4) = [character(8) :: &
            'rmnc', 'bmnc', 'bsupumnc', 'bsupvmnc']
        character(len=8), parameter :: partner(4) = [character(8) :: &
            'rmns', 'bmns', 'bsupumns', 'bsupvmns']

        call nc(nf90_open(source, nf90_nowrite, ncid))
        call get_int(ncid, 'ns', ns)
        call get_int(ncid, 'mnmax', mnmax)
        call get_int(ncid, 'mnmax_nyq', nyq)
        call get_int(ncid, 'nfp', nfp)
        allocate (xm(mnmax), xn(mnmax), xmn(nyq), xnn(nyq))
        call get_vec(ncid, 'xm', xm)
        call get_vec(ncid, 'xn', xn)
        call get_vec(ncid, 'xm_nyq', xmn)
        call get_vec(ncid, 'xn_nyq', xnn)
        cs = cos(xn*phi0)
        sn = sin(xn*phi0)
        cn = cos(xnn*phi0)
        snn = sin(xnn*phi0)

        call nc(nf90_create(target, nf90_clobber, out))
        call nc(nf90_def_dim(out, 'radius', ns, d_r))
        call nc(nf90_def_dim(out, 'mn_mode', mnmax, d_mn))
        call nc(nf90_def_dim(out, 'mn_mode_nyq', nyq, d_nyq))
        call def_int(out, 'ns')
        call def_int(out, 'mnmax')
        call def_int(out, 'mnmax_nyq')
        call def_int(out, 'nfp')
        call def_int(out, 'lasym__logical__')
        call nc(nf90_def_var(out, 'xm', nf90_double, [d_mn], var))
        call nc(nf90_def_var(out, 'xn', nf90_double, [d_mn], var))
        call nc(nf90_def_var(out, 'xm_nyq', nf90_double, [d_nyq], var))
        call nc(nf90_def_var(out, 'xn_nyq', nf90_double, [d_nyq], var))
        call nc(nf90_def_var(out, 'rmnc', nf90_double, [d_mn, d_r], var))
        call nc(nf90_def_var(out, 'rmns', nf90_double, [d_mn, d_r], var))
        call nc(nf90_def_var(out, 'zmns', nf90_double, [d_mn, d_r], var))
        call nc(nf90_def_var(out, 'zmnc', nf90_double, [d_mn, d_r], var))
        do j = 2, 4
            call nc(nf90_def_var(out, trim(cosine(j)), nf90_double, [d_nyq, d_r], var))
            call nc(nf90_def_var(out, trim(partner(j)), nf90_double, [d_nyq, d_r], var))
        end do
        call nc(nf90_enddef(out))
        call put_int(out, 'ns', ns)
        call put_int(out, 'mnmax', mnmax)
        call put_int(out, 'mnmax_nyq', nyq)
        call put_int(out, 'nfp', nfp)
        call put_int(out, 'lasym__logical__', 1)
        call put_vec(out, 'xm', xm)
        call put_vec(out, 'xn', xn)
        call put_vec(out, 'xm_nyq', xmn)
        call put_vec(out, 'xn_nyq', xnn)
        allocate (a(mnmax, ns))
        call get_mat(ncid, 'rmnc', a)
        call put_mat(out, 'rmnc', spread(cs, 2, ns)*a)
        call put_mat(out, 'rmns', -spread(sn, 2, ns)*a)
        call get_mat(ncid, 'zmns', a)
        call put_mat(out, 'zmns', spread(cs, 2, ns)*a)
        call put_mat(out, 'zmnc', spread(sn, 2, ns)*a)
        deallocate (a)
        allocate (a(nyq, ns))
        do j = 2, 4
            call get_mat(ncid, trim(cosine(j)), a)
            call put_mat(out, trim(cosine(j)), spread(cn, 2, ns)*a)
            call put_mat(out, trim(partner(j)), -spread(snn, 2, ns)*a)
        end do
        call nc(nf90_close(out))
        call nc(nf90_close(ncid))
    end subroutine write_rotated

    subroutine def_int(out, name)
        integer, intent(in) :: out
        character(len=*), intent(in) :: name
        integer :: var
        call nc(nf90_def_var(out, name, nf90_int, var))
    end subroutine def_int
    subroutine put_int(out, name, value)
        integer, intent(in) :: out
        character(len=*), intent(in) :: name
        integer, intent(in) :: value
        integer :: var
        call nc(nf90_inq_varid(out, name, var))
        call nc(nf90_put_var(out, var, value))
    end subroutine put_int
    subroutine put_vec(out, name, value)
        integer, intent(in) :: out
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: value(:)
        integer :: var
        call nc(nf90_inq_varid(out, name, var))
        call nc(nf90_put_var(out, var, value))
    end subroutine put_vec
    subroutine put_mat(out, name, value)
        integer, intent(in) :: out
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: value(:, :)
        integer :: var
        call nc(nf90_inq_varid(out, name, var))
        call nc(nf90_put_var(out, var, value))
    end subroutine put_mat

    subroutine get_int(ncid, name, value)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, intent(out) :: value
        integer :: var
        call nc(nf90_inq_varid(ncid, name, var))
        call nc(nf90_get_var(ncid, var, value))
    end subroutine get_int

    subroutine get_vec(ncid, name, value)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        real(dp), intent(out) :: value(:)
        integer :: var
        call nc(nf90_inq_varid(ncid, name, var))
        call nc(nf90_get_var(ncid, var, value))
    end subroutine get_vec

    subroutine get_mat(ncid, name, value)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        real(dp), intent(out) :: value(:, :)
        integer :: var
        call nc(nf90_inq_varid(ncid, name, var))
        call nc(nf90_get_var(ncid, var, value))
    end subroutine get_mat

    subroutine nc(status)
        integer, intent(in) :: status
        if (status /= nf90_noerr) then
            write (error_unit, '(A)') 'netCDF: '//trim(nf90_strerror(status))
            error stop 1
        end if
    end subroutine nc

    subroutine check(ok, name)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: name
        if (ok) then
            print '(A)', 'ok   '//name
        else
            write (error_unit, '(A)') 'FAIL '//name
            failures = failures + 1
        end if
    end subroutine check

    subroutine close_abs(name, a, b, tol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: a, b, tol
        if (abs(a - b) > tol .or. a /= a) then
            write (error_unit, '(A,3ES14.6)') 'FAIL '//name//': ', a, b, a - b
            failures = failures + 1
        else
            print '(A,3ES14.6)', 'ok   '//name//': ', a, b, a - b
        end if
    end subroutine close_abs
end program test_vmec_lasym
