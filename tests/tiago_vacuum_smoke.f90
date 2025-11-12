program tiago_vacuum_smoke
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_diagnostic_types, only: flux_loop_t, segmented_rogowski_t
    use tiago_flux_loops, only: read_flux_loop_file
    use tiago_segmented_rogowski, only: read_segmented_rogowski_file
    use tiago_vacuum_forward, only: vacuum_solver_t, quadrature_rule_t
    implicit none

    character(len=256) :: coil_path
    character(len=256) :: flux_path
    character(len=256) :: segrog_path
    integer :: argc
    integer :: ierr
    character(len=:), allocatable :: message
    type(vacuum_solver_t) :: solver
    type(quadrature_rule_t) :: rule
    type(quadrature_rule_t) :: seg_rule
    real(dp), allocatable :: fluxes(:)
    real(dp), allocatable :: voltages(:)
    type(flux_loop_t), allocatable :: loops(:)
    type(segmented_rogowski_t), allocatable :: segs(:)

    argc = command_argument_count()
    if (argc < 3) then
        write(error_unit, '(A)') 'usage: tiago_vacuum_smoke <coil> <flux> <seg>'
        stop 1
    end if

    call get_command_argument(1, coil_path)
    call get_command_argument(2, flux_path)
    call get_command_argument(3, segrog_path)

    call read_flux_loop_file(trim(flux_path), loops, ierr, message)
    if (ierr /= 0) then
        write(error_unit, '(2A)') 'flux parse failed: ', trim(message)
        stop 2
    end if

    call read_segmented_rogowski_file(trim(segrog_path), segs, ierr, message, &
        3.40e-4_dp)
    if (ierr /= 0) then
        write(error_unit, '(2A)') 'segrog parse failed: ', trim(message)
        stop 2
    end if

    call solver%init(trim(coil_path))

    rule%samples_per_segment = 6
    call solver%flux_loops(loops, fluxes, rule)
    call solver%segrog(segs, voltages)

    call assert_close(fluxes(1), 1.0271e-11_dp, 5.0e-13_dp, 'flux sample')
    call assert_close(voltages(1), 3.2570e-14_dp, 1.0e-15_dp, &
        'segrog sample')

    call solver%finalize()
end program tiago_vacuum_smoke

subroutine assert_close(value, expected, tolerance, label)
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    real(dp), intent(in) :: value
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: label
    if (abs(value - expected) > tolerance) then
        write(error_unit, '(A,1X,ES12.4,1X,ES12.4)') trim(label)//' mismatch', &
            value, expected
        stop 3
    end if
end subroutine assert_close
