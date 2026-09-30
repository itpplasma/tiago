program test_plasma_support
    !! Plasma sheet-current model checks that need no reference code.
    !! Usage: test_plasma_support <wout.nc>
    !!  1. Ampere: a closed poloidal loop around the plasma gives mu0 * I_tor (VMEC ctor).
    !!  2. curl A = B: the flux of A around a small square equals B.n * area.
    !!  4. The boundary-shape Jacobian matches central finite differences.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_get_var, nf90_nowrite
    use tiago_plasma_support, only: plasma_support_t, vmec_boundary_t, read_vmec_boundary
    implicit none
    !!  3. The boundary-field Jacobian reproduces the plasma signals exactly.

    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 4.0e-7_dp * pi
    integer, parameter :: n = 400
    type(plasma_support_t) :: plasma
    character(len=512) :: wout
    real(dp) :: ctor, rmajor, aminor, t, circ, flux_a, flux_b, h
    real(dp) :: pts(3, n), dls(3, n), field(3, n), sq(3, 4), sq_dl(3, 4), a(3, 4), c(3, 1), bc(3, 1)
    real(dp), allocatable :: w(:, :, :), resp(:, :), gx(:, :, :), wf(:, :, :), gxf(:, :, :)
    real(dp), allocatable :: shape_pot(:, :), shape_field(:, :)
    type(vmec_boundary_t) :: vb
    integer :: k, failures, owner(n), sq_owner(4)

    failures = 0
    call get_command_argument(1, wout)
    call read_scalars(trim(wout), ctor, rmajor, aminor)
    call plasma%init_from_vmec(trim(wout), 32, 32)

    ! 1. Ampere around the plasma cross-section at phi = 0 (circle in the x-z plane).
    do k = 1, n
        t = 2.0_dp * pi * (real(k, dp) - 0.5_dp) / real(n, dp)
        pts(:, k) = [rmajor + 3.0_dp * aminor * cos(t), 0.0_dp, 3.0_dp * aminor * sin(t)]
        dls(:, k) = [-sin(t), 0.0_dp, cos(t)] * 3.0_dp * aminor * 2.0_dp * pi / real(n, dp)
    end do
    call plasma%sample_bfield(pts, field)
    circ = sum(field * dls)
    call check('Ampere: |loop integral of B| = mu0 |I_tor|', abs(circ), mu0 * abs(ctor), 2.0e-3_dp)

    ! 2. curl A = B on a 1 cm square in the x-y plane outside the plasma.
    h = 0.01_dp
    c(:, 1) = [rmajor + 3.0_dp * aminor, 0.0_dp, 0.1_dp]
    sq(:, 1) = c(:, 1) + [0.5_dp * h, 0.0_dp, 0.0_dp]
    sq(:, 2) = c(:, 1) + [0.0_dp, 0.5_dp * h, 0.0_dp]
    sq(:, 3) = c(:, 1) - [0.5_dp * h, 0.0_dp, 0.0_dp]
    sq(:, 4) = c(:, 1) - [0.0_dp, 0.5_dp * h, 0.0_dp]
    sq_dl(:, 1) = [0.0_dp, h, 0.0_dp]
    sq_dl(:, 2) = [-h, 0.0_dp, 0.0_dp]
    sq_dl(:, 3) = [0.0_dp, -h, 0.0_dp]
    sq_dl(:, 4) = [h, 0.0_dp, 0.0_dp]
    call plasma%sample_vector_potential(sq, a)
    call plasma%sample_bfield(c, bc)
    flux_a = sum(a * sq_dl)
    flux_b = bc(3, 1) * h * h
    call check('curl A = B (small loop flux)', flux_a, flux_b, 1.0e-4_dp)

    ! 3. Jacobian with respect to the VMEC boundary field coefficients.
    owner = 1
    call plasma%field_weights(pts, dls, owner, 1, w)
    call plasma%mode_response(w, resp)
    call check('Jacobian . coefficients = Ampere loop', sum(resp(1, :) * plasma%coefficients), &
        circ, 1.0e-12_dp)
    sq_owner = 1
    call plasma%potential_weights(sq, sq_dl, sq_owner, 1, w)
    call plasma%mode_response(w, resp)
    call check('Jacobian . coefficients = small-loop flux', sum(resp(1, :) * plasma%coefficients), &
        flux_a, 1.0e-12_dp)

    ! 4. Shape Jacobian vs central differences in single boundary coefficients.
    call read_vmec_boundary(trim(wout), vb)
    call plasma%init_from_boundary(vb, 32, 32)
    call plasma%potential_weights(sq, sq_dl, sq_owner, 1, w, gx)
    call plasma%shape_response(w, gx, shape_pot)
    call plasma%field_weights(pts, dls, owner, 1, wf, gxf)
    call plasma%shape_response(wf, gxf, shape_field)
    call check_shape(vb, 1)
    call check_shape(vb, 2)
    call check_shape(vb, size(vb%xm) / 2)
    call check_shape(vb, size(vb%xm))

    call plasma%finalize()
    if (failures > 0) error stop 1
    print '(A)', 'test_plasma_support passed'

contains

    subroutine check_shape(vb, mode)
        !! rmnc(mode) and zmns(mode) columns of both signals vs central differences.
        type(vmec_boundary_t), intent(in) :: vb
        integer, intent(in) :: mode
        type(vmec_boundary_t) :: pert
        character(len=64) :: name
        real(dp) :: hstep, fd(2), plus(2), minus(2)
        integer :: which, col

        do which = 1, 2
            col = mode + (which - 1) * size(vb%xm)
            hstep = 1.0e-6_dp
            pert = vb
            call shift(pert, which, mode, hstep)
            call two_signals(pert, plus)
            pert = vb
            call shift(pert, which, mode, -hstep)
            call two_signals(pert, minus)
            fd = (plus - minus) / (2.0_dp * hstep)
            write(name, '(A,I0,A,I0)') 'shape Jacobian, flux, column ', col, ' mode ', mode
            call check_abs(trim(name), shape_pot(1, col), fd(1), 1.0e-6_dp * maxval(abs(shape_pot)))
            write(name, '(A,I0,A,I0)') 'shape Jacobian, Ampere, column ', col, ' mode ', mode
            call check_abs(trim(name), shape_field(1, col), fd(2), 1.0e-6_dp * maxval(abs(shape_field)))
        end do
    end subroutine check_shape

    subroutine shift(vb, which, mode, step)
        type(vmec_boundary_t), intent(inout) :: vb
        integer, intent(in) :: which, mode
        real(dp), intent(in) :: step
        if (which == 1) vb%rmnc(mode) = vb%rmnc(mode) + step
        if (which == 2) vb%zmns(mode) = vb%zmns(mode) + step
    end subroutine shift

    subroutine two_signals(vb, values)
        !! Small-square flux (A) and Ampere loop (B) of the boundary vb.
        type(vmec_boundary_t), intent(in) :: vb
        real(dp), intent(out) :: values(2)
        type(plasma_support_t) :: p
        real(dp) :: av(3, 4), bv(3, n)

        call p%init_from_boundary(vb, 32, 32)
        call p%sample_vector_potential(sq, av)
        call p%sample_bfield(pts, bv)
        values = [sum(av * sq_dl), sum(bv * dls)]
        call p%finalize()
    end subroutine two_signals

    subroutine check_abs(name, actual, expected, atol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual, expected, atol

        if (abs(actual - expected) > atol) then
            write(error_unit, '(A,2ES24.16)') 'FAIL '//name//': ', actual, expected
            failures = failures + 1
        else
            print '(A,2ES24.16)', 'ok   '//name//': ', actual, expected
        end if
    end subroutine check_abs

    subroutine read_scalars(path, ctor, rmajor, aminor)
        character(len=*), intent(in) :: path
        real(dp), intent(out) :: ctor, rmajor, aminor
        integer :: ncid, varid, status

        status = nf90_open(path, nf90_nowrite, ncid)
        if (status /= 0) error stop 'cannot open wout file'
        status = nf90_inq_varid(ncid, 'ctor', varid)
        status = nf90_get_var(ncid, varid, ctor)
        status = nf90_inq_varid(ncid, 'Rmajor_p', varid)
        status = nf90_get_var(ncid, varid, rmajor)
        status = nf90_inq_varid(ncid, 'Aminor_p', varid)
        status = nf90_get_var(ncid, varid, aminor)
        status = nf90_close(ncid)
    end subroutine read_scalars

    subroutine check(name, actual, expected, rtol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual, expected, rtol

        if (abs(actual - expected) > rtol * abs(expected)) then
            write(error_unit, '(A,2ES24.16)') 'FAIL '//name//': ', actual, expected
            failures = failures + 1
        else
            print '(A,2ES24.16)', 'ok   '//name//': ', actual, expected
        end if
    end subroutine check
end program test_plasma_support
