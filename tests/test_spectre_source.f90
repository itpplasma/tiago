program test_spectre_source
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use field_spectre, only: spectre_field_t, create_spectre_field
    use tiago_spectre_source, only: spectre_source_t
    use tiago_vacuum_forward, only: vacuum_solver_t
    implicit none
    type(spectre_source_t) :: source
    type(spectre_field_t) :: field
    type(vacuum_solver_t) :: solver
    character(len=1024) :: fixture, reffile, mode, line
    real(dp) :: row(11), u(3), x(3, 1), b(3, 1), expected(3)
    real(dp) :: ctr(3), e(3, 3), worst, h, delta(3), ap(3, 1), am(3, 1)
    real(dp) :: da(3, 3), curl(3), p(3, 1), rel
    integer :: ierr, unit, io, nrows, k

    call get_command_argument(1, fixture)
    call get_command_argument(2, reffile)
    call get_command_argument(3, mode)
    call create_spectre_field(field, trim(fixture), ierr)
    if (ierr /= 0) error stop 'fixture load failed'
    call source%init(trim(fixture), 32, 32)
    call solver%init('')
    call solver%set_plasma_source(source)
    if (len_trim(mode) > 0) then
        if (trim(mode) == 'reject_missing_flux') then
            rel = solver%plasma_flux()
            print *, rel
            error stop 'expected missing linked flux rejection'
        end if
        u = [0.3_dp, 0.7_dp, 0.2_dp]
        call field%coords%evaluate_cart(u, x(:, 1))
        if (trim(mode) == 'reject_exterior') then
            call source%sample_bfield(x, b)
        else
            x(:, 1) = [100.0_dp, 0.0_dp, 0.0_dp]
            call source%sample_interior_bfield(x, b)
        end if
        error stop 'expected domain rejection'
    end if

    ! Independent SPECTRE reference: native local-s contravariant components,
    ! exported by SPECTRE itself, rather than recomputed by the libneo reader.
    open (newunit=unit, file=trim(reffile), status='old', action='read')
    read (unit, '(A)') line
    worst = 0.0_dp
    nrows = 0
    do
        read (unit, *, iostat=io) row
        if (io /= 0) exit
        nrows = nrows + 1
        u = [row(1) - 1.0_dp + 0.5_dp*(row(2) + 1.0_dp), row(3), row(4)]
        call field%coords%evaluate_cart(u, x(:, 1))
        call solver%sample_bfield(x, b, interior=.true.)
        call field%coords%covariant_basis(u, e)
        ctr = [0.5_dp*row(7), row(8), row(9)]
        expected = matmul(e, ctr)
        worst = max(worst, norm2(b(:, 1) - expected)/norm2(expected))
    end do
    close (unit)
    if (nrows /= 90) error stop 'missing reference rows'
    call check('90 interior Cartesian B samples vs SPECTRE reference', worst, 1.0e-8_dp)

    ! curl A = B, finite differences in physical Cartesian coordinates, away
    ! from the sheet. This checks the exterior A required by flux loops.
    worst = 0.0_dp
    do k = 1, 3
        x(:, 1) = [14.0_dp + real(k, dp), 0.3_dp, 0.2_dp*real(k, dp)]
        h = 1.0e-4_dp
        do ierr = 1, 3
            delta = 0.0_dp
            delta(ierr) = h
            p(:, 1) = x(:, 1) + delta
            call source%sample_vector_potential(p, ap)
            p(:, 1) = x(:, 1) - delta
            call source%sample_vector_potential(p, am)
            da(:, ierr) = (ap(:, 1) - am(:, 1))/(2.0_dp*h)
        end do
        curl = [da(3, 2) - da(2, 3), da(1, 3) - da(3, 1), da(2, 1) - da(1, 2)]
        call solver%sample_bfield(x, b)
        rel = norm2(curl - b(:, 1))/norm2(b(:, 1))
        worst = max(worst, rel)
    end do
    call check('exterior Cartesian curl A = B', worst, 1.0e-6_dp)
contains
    subroutine check(label, error, tolerance)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: error, tolerance
        print '(A,1X,ES14.6)', label, error
        if (.not. ieee_is_finite(error)) error stop 'nonfinite error'
        if (error > tolerance) error stop 'SPECTRE source check failed'
    end subroutine check
end program test_spectre_source
