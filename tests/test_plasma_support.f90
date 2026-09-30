program test_plasma_support
    !! Plasma sheet-current model checks that need no reference code.
    !! Usage: test_plasma_support <wout.nc>
    !!  1. Ampere: a closed poloidal loop around the plasma gives mu0 * I_tor (VMEC ctor).
    !!  2. curl A = B: the flux of A around a small square equals B.n * area.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_get_var, nf90_nowrite
    use tiago_plasma_support, only: plasma_support_t
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 4.0e-7_dp * pi
    integer, parameter :: n = 400
    type(plasma_support_t) :: plasma
    character(len=512) :: wout
    real(dp) :: ctor, rmajor, aminor, t, circ, flux_a, flux_b, h
    real(dp) :: pts(n, 3), dls(n, 3), field(n, 3), sq(4, 3), sq_dl(4, 3), a(4, 3), c(1, 3), bc(1, 3)
    integer :: k, failures

    failures = 0
    call get_command_argument(1, wout)
    call read_scalars(trim(wout), ctor, rmajor, aminor)
    call plasma%init_from_vmec(trim(wout), 32, 32)

    ! 1. Ampere around the plasma cross-section at phi = 0 (circle in the x-z plane).
    do k = 1, n
        t = 2.0_dp * pi * (real(k, dp) - 0.5_dp) / real(n, dp)
        pts(k, :) = [rmajor + 3.0_dp * aminor * cos(t), 0.0_dp, 3.0_dp * aminor * sin(t)]
        dls(k, :) = [-sin(t), 0.0_dp, cos(t)] * 3.0_dp * aminor * 2.0_dp * pi / real(n, dp)
    end do
    call plasma%sample_bfield(pts, field)
    circ = sum(field * dls)
    call check('Ampere: |loop integral of B| = mu0 |I_tor|', abs(circ), mu0 * abs(ctor), 2.0e-3_dp)

    ! 2. curl A = B on a 1 cm square in the x-y plane outside the plasma.
    h = 0.01_dp
    c(1, :) = [rmajor + 3.0_dp * aminor, 0.0_dp, 0.1_dp]
    sq(1, :) = c(1, :) + [0.5_dp * h, 0.0_dp, 0.0_dp]
    sq(2, :) = c(1, :) + [0.0_dp, 0.5_dp * h, 0.0_dp]
    sq(3, :) = c(1, :) - [0.5_dp * h, 0.0_dp, 0.0_dp]
    sq(4, :) = c(1, :) - [0.0_dp, 0.5_dp * h, 0.0_dp]
    sq_dl(1, :) = [0.0_dp, h, 0.0_dp]
    sq_dl(2, :) = [-h, 0.0_dp, 0.0_dp]
    sq_dl(3, :) = [0.0_dp, -h, 0.0_dp]
    sq_dl(4, :) = [h, 0.0_dp, 0.0_dp]
    call plasma%sample_vector_potential(sq, a)
    call plasma%sample_bfield(c, bc)
    flux_a = sum(a * sq_dl)
    flux_b = bc(1, 3) * h * h
    call check('curl A = B (small loop flux)', flux_a, flux_b, 1.0e-4_dp)

    call plasma%finalize()
    if (failures > 0) error stop 1
    print '(A)', 'test_plasma_support passed'

contains

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
