program test_spectre_vmec
    !! Nested-limit cross-code comparison. The fixed VMEC interfaces are not
    !! force-relaxed; this tests conversion + diagnostics, not MRxMHD balance.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_get_var, &
                      nf90_nowrite, nf90_noerr
    use tiago_spectre_source, only: spectre_source_t
    use tiago_plasma_support, only: plasma_support_t
    use tiago_vacuum_forward, only: vacuum_solver_t, quadrature_rule_t
    use tiago_diagnostic_types, only: flux_loop_t, segmented_rogowski_t, &
                                      bprobe_t, loop_point_t
    use tiago_flux_loops, only: finalize_flux_signals
    implicit none
    type(spectre_source_t) :: source
    type(plasma_support_t) :: vmec
    type(vacuum_solver_t) :: sp_solver, vm_solver
    type(quadrature_rule_t) :: rule
    type(flux_loop_t), allocatable :: loops(:)
    type(segmented_rogowski_t), allocatable :: segs(:)
    type(bprobe_t) :: probes(24)
    real(dp) :: points(3, 24), bsp(3, 24), bvm(3, 24), asp(3, 24), avm(3, 24)
    real(dp) :: rmajor, aminor, t, phi, r
    real(dp) :: row(6), point(3, 1), b(3, 1), worst, sumsq, relative
    real(dp), allocatable :: sp(:), vm(:)
    integer :: j, k, n, unit, status, ngrid
    character(len=1024) :: h5file, wout, reference, line
    real(dp), parameter :: twopi = 2.0_dp*acos(-1.0_dp)

    call get_command_argument(1, h5file)
    call get_command_argument(2, wout)
    call get_command_argument(3, reference)
    ngrid = 64
    call get_command_argument(4, line)
    if (len_trim(line) > 0) read (line, *) ngrid
    call read_radii(trim(wout), rmajor, aminor)
    call vmec%init_from_vmec(trim(wout), ngrid, ngrid)
    call source%init(trim(h5file), ngrid, ngrid, linked_flux=vmec%diamagnetic_flux)
    ! Independent VMEC Fourier/SciPy oracle at the same physical points,
    ! including both radial sample locations in every SPECTRE volume.
    open (newunit=unit, file=trim(reference), status='old', action='read')
    read (unit, '(A)') line
    worst = 0.0_dp
    sumsq = 0.0_dp
    n = 0
    do
        read (unit, *, iostat=status) row
        if (status /= 0) exit
        n = n + 1
        point(:, 1) = row(:3)
        call source%sample_interior_bfield(point, b)
        relative = norm2(b(:, 1) - row(4:6))/norm2(row(4:6))
        worst = max(worst, relative)
        sumsq = sumsq + relative**2
    end do
    close (unit)
    if (n /= 512) error stop 'missing VMEC reference points'
    call check('nested interior B max relative difference', worst, 2.0e-3_dp)
    call check('nested interior B RMS relative difference', sqrt(sumsq/n), 6.0e-4_dp)
    do k = 1, 3
        phi = 0.37_dp*real(k - 1, dp)
        do j = 1, 8
            n = j + 8*(k - 1)
            t = twopi*real(j - 1, dp)/8.0_dp
            r = rmajor + 3.0_dp*aminor*cos(t)
            points(:, n) = [r*cos(phi), r*sin(phi), 3.0_dp*aminor*sin(t)]
        end do
    end do
    call source%sample_bfield(points, bsp)
    call vmec%sample_bfield(points, bvm)
    call source%sample_vector_potential(points, asp)
    call vmec%sample_vector_potential(points, avm)
  call check('nested max component B difference [T]', maxval(abs(bsp - bvm)), 2.0e-3_dp)
call check('nested max component A difference [T m]', maxval(abs(asp - avm)), 2.0e-3_dp)
    print '(A,2ES24.16)', 'linked flux SPECTRE / VMEC [Wb] ', &
        source%diamagnetic_flux, vmec%diamagnetic_flux
    call check('nested linked toroidal flux difference [Wb]', &
               abs(source%diamagnetic_flux - vmec%diamagnetic_flux), 1.0e-9_dp)

    call sp_solver%init('')
    call sp_solver%set_plasma_source(source)
    call vm_solver%init('')
    call vm_solver%enable_plasma_from_vmec(trim(wout), ngrid, ngrid)
    do j = 1, 24
        probes(j)%position = points(:, j)
        probes(j)%normal = [0.0_dp, 0.0_dp, 1.0_dp]
        probes(j)%eff_area = 1.0e-4_dp
    end do
    call sp_solver%bprobes(probes, sp)
    call vm_solver%bprobes(probes, vm)
    call check('nested probe difference [Wb]', maxval(abs(sp - vm)), 2.0e-7_dp)

    allocate (loops(1), segs(1))
    allocate (loops(1)%points(32), segs(1)%path(33), segs(1)%segment_area(32))
    ! A poloidal loop at three minor radii: every chord remains exterior.
    do j = 1, 32
        t = twopi*real(j - 1, dp)/32.0_dp
        loops(1)%points(j) = loop_point_t(rmajor + 3.0_dp*aminor*cos(t), &
                                          0.0_dp, 3.0_dp*aminor*sin(t))
    end do
    loops(1)%idia = 1
    segs(1)%path(:32) = loops(1)%points
    segs(1)%path(33) = loops(1)%points(1)
    segs(1)%segment_area = 1.0e-4_dp
    rule%gauss = .true.
    rule%samples_per_segment = 4
    call sp_solver%flux_loops(loops, sp, rule)
    call vm_solver%flux_loops(loops, vm, rule)
    call finalize_flux_signals(loops, sp, sp_solver%plasma_flux())
    call finalize_flux_signals(loops, vm, vm_solver%plasma_flux())
  call check('nested linked flux loop difference [Wb]', maxval(abs(sp - vm)), 5.0e-2_dp)
    call sp_solver%segrog(segs, sp, rule)
    call vm_solver%segrog(segs, vm, rule)
    call check('nested Rogowski difference [T m^3]', maxval(abs(sp - vm)), 5.0e-6_dp)
contains
    subroutine check(label, error, tolerance)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: error, tolerance
        print '(A,1X,ES14.6)', label, error
        if (.not. ieee_is_finite(error)) error stop 'nonfinite error'
        if (error > tolerance) error stop 'nested comparison failed'
    end subroutine check

    subroutine read_radii(path, rmajor, aminor)
        character(len=*), intent(in) :: path
        real(dp), intent(out) :: rmajor, aminor
        integer :: id, vid, status
        status = nf90_open(path, nf90_nowrite, id)
        if (status /= nf90_noerr) error stop 'cannot open VMEC reference'
        status = nf90_inq_varid(id, 'Rmajor_p', vid)
        if (status /= nf90_noerr) error stop 'missing Rmajor_p'
        status = nf90_get_var(id, vid, rmajor)
        if (status /= nf90_noerr) error stop 'cannot read Rmajor_p'
        status = nf90_inq_varid(id, 'Aminor_p', vid)
        if (status /= nf90_noerr) error stop 'missing Aminor_p'
        status = nf90_get_var(id, vid, aminor)
        if (status /= nf90_noerr) error stop 'cannot read Aminor_p'
        status = nf90_close(id)
        if (status /= nf90_noerr) error stop 'cannot close VMEC reference'
    end subroutine read_radii
end program test_spectre_vmec
