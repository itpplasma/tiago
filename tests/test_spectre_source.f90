program test_spectre_source
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use field_spectre, only: spectre_field_t, create_spectre_field
    use tiago_spectre_source, only: spectre_source_t
    use tiago_plasma_field_source, only: field_hint_t
    use tiago_plasma_support, only: plasma_support_t
    use tiago_vacuum_forward, only: vacuum_solver_t
    implicit none
    type(spectre_source_t) :: source, not_initialized
    type(plasma_support_t) :: no_interior
    type(spectre_field_t) :: field
    type(vacuum_solver_t) :: solver
    type(field_hint_t) :: hint
    character(len=1024) :: fixture, reffile, mode, line
    real(dp) :: row(11), u(3), x(3, 1), b(3, 1), expected(3)
    real(dp) :: ctr(3), e(3, 3), worst, h, delta(3), ap(3, 1), am(3, 1)
    real(dp) :: da(3, 3), curl(3), p(3, 1), rel, checked(3, 1)
    integer :: point_status(1), batch_status(3)
    real(dp) :: batch_points(3, 3), batch_values(3, 3), warm(3), warm_worst
    integer :: ierr, unit, io, nrows, k, warm_status

    call get_command_argument(1, fixture)
    call get_command_argument(2, reffile)
    call get_command_argument(3, mode)
    call create_spectre_field(field, trim(fixture), ierr)
    if (ierr /= 0) error stop 'fixture load failed'
    call source%init(trim(fixture), 32, 32)
    x(:, 1) = [100.0_dp, 0.0_dp, 0.0_dp]
    call not_initialized%sample_interior_bfield_status(x, checked, point_status)
    if (point_status(1) == 0 .or. any(checked /= 0.0_dp)) &
        error stop 'uninitialized status sampling must fail with zero field'
    hint%u = [0.3_dp, 0.7_dp, 0.2_dp]
    hint%volume = 1
    call not_initialized%sample_interior_bfield_with_hint(x(:, 1), warm, &
                                                        warm_status, hint)
    call failed_hint('uninitialized hinted sample', warm, warm_status, hint)
    call no_interior%sample_interior_bfield_status(x, checked, point_status)
    if (point_status(1) == 0 .or. any(checked /= 0.0_dp)) &
        error stop 'unsupported status sampling must fail with zero field'
    hint%volume = 1
    call no_interior%sample_interior_bfield_with_hint(x(:, 1), warm, &
                                                    warm_status, hint)
    call failed_hint('unsupported hinted sample', warm, warm_status, hint)
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
    warm_worst = 0.0_dp
    nrows = 0
    do
        read (unit, *, iostat=io) row
        if (io /= 0) exit
        nrows = nrows + 1
        u = [row(1) - 1.0_dp + 0.5_dp*(row(2) + 1.0_dp), row(3), row(4)]
        call field%coords%evaluate_cart(u, x(:, 1))
        call solver%sample_bfield(x, b, interior=.true.)
        call source%sample_interior_bfield_status(x, checked, point_status)
        if (point_status(1) /= 0) error stop 'valid interior status rejected'
        if (any(checked /= b)) error stop 'status sampling changed interior field'
        call field%coords%covariant_basis(u, e)
        ctr = [0.5_dp*row(7), row(8), row(9)]
        expected = matmul(e, ctr)
        worst = max(worst, norm2(b(:, 1) - expected)/norm2(expected))
        call source%sample_interior_bfield_with_hint(x(:, 1), warm, &
                                                    warm_status, hint)
        if (warm_status /= 0 .or. hint%volume < 1) &
            error stop 'valid hinted sample rejected'
        warm_worst = max(warm_worst, norm2(warm - expected)/norm2(expected))
        if (nrows == 1) then
            if (field%data%Mvol < 2) error stop 'wrong-volume test needs two volumes'
            hint%u = u
            hint%volume = modulo(nint(row(1)), field%data%Mvol) + 1
            call source%sample_interior_bfield_with_hint(x(:, 1), warm, &
                                                        warm_status, hint)
            if (warm_status /= 0) error stop 'wrong-volume hint blocked cold fallback'
            call check('wrong-volume fallback vs SPECTRE reference', &
                       norm2(warm - expected)/norm2(expected), 1.0e-8_dp)
        end if
    end do
    close (unit)
    if (nrows /= 90) error stop 'missing reference rows'
    call check('90 interior Cartesian B samples vs SPECTRE reference', worst, 1.0e-8_dp)
    call check('90 hinted Cartesian B samples vs SPECTRE reference', &
               warm_worst, 1.0e-8_dp)
    call boundary_hints()

    hint%u = ieee_value(0.0_dp, ieee_quiet_nan)
    hint%volume = 1
    call source%sample_interior_bfield_with_hint(x(:, 1), warm, warm_status, hint)
    if (warm_status /= 0) error stop 'nonfinite hint blocked cold fallback'
    call check('nonfinite hint recovery', &
               norm2(warm - b(:, 1))/norm2(b(:, 1)), 1.0e-8_dp)
    p(:, 1) = [ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, 0.0_dp]
    call source%sample_interior_bfield_with_hint(p(:, 1), warm, warm_status, hint)
    call failed_hint('nonfinite point invalidates hint', warm, warm_status, hint)
    call source%sample_interior_bfield_with_hint(x(:, 1), warm, warm_status, hint)
    if (warm_status /= 0) error stop 'valid point after nonfinite point rejected'
    call check('recovery after nonfinite point', &
               norm2(warm - b(:, 1))/norm2(b(:, 1)), 1.0e-8_dp)
    p(:, 1) = [100.0_dp, 0.0_dp, 0.0_dp]
    call source%sample_interior_bfield_with_hint(p(:, 1), warm, warm_status, hint)
    call failed_hint('exterior point invalidates hint', warm, warm_status, hint)
    hint%u = u
    hint%volume = field%data%Mvol + 1
    call source%sample_interior_bfield_with_hint(x(:, 1), warm, warm_status, hint)
    if (warm_status /= 0) error stop 'out-of-range hint blocked recovery'
    call check('recovery after exterior point', &
               norm2(warm - b(:, 1))/norm2(b(:, 1)), 1.0e-8_dp)

    batch_points(:, 1) = x(:, 1)
    batch_points(:, 2) = [100.0_dp, 0.0_dp, 0.0_dp]
    batch_points(:, 3) = [ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, 0.0_dp]
    call source%sample_interior_bfield_status(batch_points, batch_values, batch_status)
    if (batch_status(1) /= 0 .or. any(batch_status(2:) == 0)) &
        error stop 'mixed batch status does not identify failed points'
    if (any(batch_values(:, 1) /= b(:, 1))) error stop 'mixed batch changed valid field'
    if (any(batch_values(:, 2:) /= 0.0_dp)) error stop 'failed batch samples must be zero'

    x(:, 1) = [100.0_dp, 0.0_dp, 0.0_dp]
    call source%sample_interior_bfield_status(x, checked, point_status)
    if (point_status(1) == 0) error stop 'exterior interior-status sample accepted'
    if (any(checked /= 0.0_dp)) error stop 'failed sample must return a zero field'

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
    subroutine boundary_hints()
        real(dp) :: labels(3), point(3, 1), cold(3, 1), hinted(3)
        type(field_hint_t) :: seed
        integer :: side, volume, cold_status(1), hinted_status

        if (field%data%Nvol < 2) error stop 'interface test needs two plasma volumes'
        do side = -1, 1, 2
            labels = [1.0_dp + real(side, dp)*1.0e-10_dp, 0.7_dp, 0.2_dp]
            call field%coords%evaluate_cart(labels, point(:, 1))
            call source%sample_interior_bfield_status(point, cold, cold_status)
            if (cold_status(1) /= 0) error stop 'internal interface cold sample rejected'
            do volume = 1, 2
                ! Both neighboring seeds can snap to the same interface.
                seed%u = labels
                seed%volume = volume
                call source%sample_interior_bfield_with_hint(point(:, 1), hinted, &
                                                            hinted_status, seed)
                if (hinted_status /= 0) error stop 'internal interface hinted sample rejected'
                call check('interface sheet side vs cold sample', &
                           norm2(hinted - cold(:, 1))/norm2(cold(:, 1)), 1.0e-10_dp)
            end do
        end do

        do side = -1, 1, 2
            labels = [real(field%data%Nvol, dp) + real(side, dp)*1.0e-10_dp, &
                      0.7_dp, 0.2_dp]
            call field%coords%evaluate_cart(labels, point(:, 1))
            call source%sample_interior_bfield_status(point, cold, cold_status)
            seed%u = labels
            seed%volume = field%data%Nvol
            call source%sample_interior_bfield_with_hint(point(:, 1), hinted, &
                                                        hinted_status, seed)
            if (side < 0) then
                if (cold_status(1) /= 0 .or. hinted_status /= 0) &
                    error stop 'just-inside plasma sample rejected'
                call check('just-inside plasma vs cold sample', &
                           norm2(hinted - cold(:, 1))/norm2(cold(:, 1)), 1.0e-10_dp)
            else
                if (cold_status(1) == 0 .or. any(cold /= 0.0_dp)) &
                    error stop 'just-outside plasma cold sample accepted'
                call failed_hint('just-outside plasma clears hint', hinted, hinted_status, seed)
            end if
        end do
    end subroutine boundary_hints

    subroutine failed_hint(label, value, status, hint)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: value(3)
        integer, intent(in) :: status
        type(field_hint_t), intent(in) :: hint

        if (status /= 0 .and. all(value == 0.0_dp) .and. hint%volume == 0 &
            .and. all(hint%u == 0.0_dp)) return
        print '(A)', label
        error stop 'failed hinted sample must return zero field and clear hint'
    end subroutine failed_hint

    subroutine check(label, error, tolerance)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: error, tolerance
        print '(A,1X,ES14.6)', label, error
        if (.not. ieee_is_finite(error)) error stop 'nonfinite error'
        if (error > tolerance) error stop 'SPECTRE source check failed'
    end subroutine check
end program test_spectre_source
