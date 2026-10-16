module uniform_test_source
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_field_source, only: plasma_field_source_t
    implicit none

    type, extends(plasma_field_source_t) :: uniform_source_t
    contains
        procedure :: sample_bfield => uniform_b
        procedure :: sample_vector_potential => uniform_a
        procedure :: sample_interior_bfield => interior_b
    end type uniform_source_t
contains
    subroutine uniform_b(self, points, bfield)
        class(uniform_source_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: bfield(:, :)
        if (size(points, 2) /= size(bfield, 2)) error stop 'field batch shape mismatch'
        bfield = 0.0_dp
        bfield(3, :) = 2.0_dp
    end subroutine uniform_b

    subroutine uniform_a(self, points, avec)
        class(uniform_source_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: avec(:, :)
        avec(1, :) = -points(2, :)
        avec(2, :) = points(1, :)
        avec(3, :) = 0.0_dp
    end subroutine uniform_a

    subroutine interior_b(self, points, values)
        class(uniform_source_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer :: j
        do j = 1, size(points, 2)
            values(:, j) = [7.0_dp, 8.0_dp, 9.0_dp]
        end do
    end subroutine interior_b
end module uniform_test_source

program test_field_source
    use uniform_test_source, only: uniform_source_t
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_support, only: plasma_support_t
    use tiago_vacuum_forward, only: vacuum_solver_t
    use tiago_diagnostic_types, only: flux_loop_t, segmented_rogowski_t, &
                                      bprobe_t, loop_point_t
    use tiago_flux_loops, only: finalize_flux_signals
    implicit none
    type(uniform_source_t) :: source
    type(plasma_support_t) :: direct
    type(vacuum_solver_t) :: solver, with_coils
    type(flux_loop_t), allocatable :: loops(:)
    type(segmented_rogowski_t), allocatable :: segs(:)
    type(bprobe_t) :: probes(1)
    real(dp), allocatable :: signals(:), reference(:)
    real(dp) :: points(3, 2), b(3, 2), expected(3, 2), err
    character(len=1024) :: coilfile, wout

    call get_command_argument(1, coilfile)
    call get_command_argument(2, wout)
    source%enabled = .true.
    source%diamagnetic_flux = 1.25_dp
    source%linked_flux_known = .true.
    call solver%init('')
    call solver%set_plasma_source(source)
    allocate (loops(1), segs(1))
    allocate (loops(1)%points(4), segs(1)%path(2), segs(1)%segment_area(1))
    loops(1)%points = [loop_point_t(0.0_dp, 0.0_dp, 0.0_dp), &
                       loop_point_t(1.0_dp, 0.0_dp, 0.0_dp), &
                       loop_point_t(1.0_dp, 1.0_dp, 0.0_dp), &
                       loop_point_t(0.0_dp, 1.0_dp, 0.0_dp)]
    call solver%flux_loops(loops, signals)
    call check(maxval(abs(signals - 2.0_dp)), 1.0e-14_dp, 'uniform flux = B * area')
    loops(1)%idia = 1
    call finalize_flux_signals(loops, signals, solver%plasma_flux())
    call check(abs(signals(1) - 3.25_dp), 1.0e-14_dp, 'linked flux metadata')

    segs(1)%path = [loop_point_t(0.0_dp, 0.0_dp, 0.0_dp), &
                    loop_point_t(0.0_dp, 0.0_dp, 1.0_dp)]
    segs(1)%segment_area = 0.4_dp
    call solver%segrog(segs, signals)
    call check(abs(signals(1) - 0.8_dp), 1.0e-14_dp, 'uniform Rogowski')
    probes(1)%position = [2.0_dp, 1.0_dp, 0.0_dp]
    probes(1)%normal = [0.0_dp, 0.0_dp, 1.0_dp]
    probes(1)%eff_area = 0.25_dp
    call solver%bprobes(probes, signals)
    call check(abs(signals(1) - 0.5_dp), 1.0e-14_dp, 'uniform probe')

    call with_coils%init(trim(coilfile))
    call with_coils%set_plasma_source(source)
    points(:, 1) = [2.0_dp, 1.0_dp, 0.1_dp]
    points(:, 2) = [1.0_dp, 2.0_dp, 0.3_dp]
    call with_coils%sample_bfield(points, b, interior=.true.)
    expected(:, 1) = [7.0_dp, 8.0_dp, 9.0_dp]
    expected(:, 2) = expected(:, 1)
    call check(maxval(abs(b - expected)), 0.0_dp, 'interior excludes coil addition')
    call with_coils%sample_bfield(points, b)
    call with_coils%coils%field(points, expected)
    expected(3, :) = expected(3, :) + 2.0_dp
    call check(maxval(abs(b - expected)), 0.0_dp, 'exterior adds coils once')

    ! A real VMEC source must preserve the existing field exactly, including
    ! dispatch through a copied abstract source and switching back to VMEC.
    call direct%init_from_vmec(trim(wout), 32, 32)
    call solver%enable_plasma_from_vmec(trim(wout), 32, 32)
    call direct%sample_bfield(points, expected)
    call solver%sample_bfield(points, b)
    err = maxval(abs(b - expected))
    call check(err, 0.0_dp, 'VMEC B unchanged through abstract dispatch')
    call solver%flux_loops(loops, reference)
    call solver%set_plasma_source(direct)
    call solver%flux_loops(loops, signals)
    call check(maxval(abs(signals - reference)), 0.0_dp, &
        'VMEC A unchanged through source copy')
    call direct%finalize()
    call solver%sample_bfield(points, b)
    call check(maxval(abs(b - expected)), 0.0_dp, 'source owns its field data')
    call solver%enable_plasma_from_vmec(trim(wout), 32, 32)
    call check(abs(solver%plasma_flux() - solver%plasma%diamagnetic_flux), &
               0.0_dp, 'switch back to VMEC')
    call solver%finalize()
    print '(A)', 'field source tests passed'
contains
    subroutine check(error, tolerance, label)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        real(dp), intent(in) :: error, tolerance
        character(len=*), intent(in) :: label
        print '(A,1X,ES14.6)', label, error
        if (.not. ieee_is_finite(error)) error stop 'nonfinite error'
        if (error > tolerance) error stop 'field source check failed'
    end subroutine check
end program test_field_source
