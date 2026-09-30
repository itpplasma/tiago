program test_vacuum_semantics
    !! Regression tests for DIAGNO flux-loop semantics (no external codes needed).
    !! Usage: test_vacuum_semantics <coils_sample.coils>
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_diagnostic_types, only: flux_loop_t
    use tiago_flux_loops, only: read_flux_loop_file
    use tiago_vacuum_forward, only: vacuum_solver_t, quadrature_rule_t
    implicit none

    type(vacuum_solver_t) :: solver
    type(quadrature_rule_t) :: rule
    character(len=512) :: coil_path
    integer :: failures
    character(len=*), parameter :: square_rows = &
        '0.4 0.4 0.3' // new_line('a') // '0.6 0.4 0.3' // new_line('a') // &
        '0.6 0.6 0.3' // new_line('a') // '0.4 0.6 0.3' // new_line('a')

    failures = 0
    call get_command_argument(1, coil_path)
    call solver%init(trim(coil_path))
    rule%samples_per_segment = 6

    call test_open_polygon_is_closed()

    call solver%finalize()
    if (failures > 0) then
        write(error_unit, '(I0,A)') failures, ' check(s) failed'
        error stop 1
    end if
    print '(A)', 'test_vacuum_semantics passed'

contains

    subroutine test_open_polygon_is_closed()
        !! DIAGNO closes every flux loop, whether or not the first point is repeated.
        character(len=*), parameter :: text = &
            '2' // new_line('a') // &
            '4 0 0 OPEN' // new_line('a') // square_rows // &
            '5 0 0 CLOSED' // new_line('a') // square_rows // '0.4 0.4 0.3' // new_line('a')
        real(dp), allocatable :: flux(:)

        call eval_file('semantics_open.diagno', text, flux)
        call check_close('open polygon equals closed polygon', flux(1), flux(2), 1.0e-12_dp)
    end subroutine test_open_polygon_is_closed

    subroutine eval_file(path, text, flux)
        character(len=*), intent(in) :: path, text
        real(dp), allocatable, intent(out) :: flux(:)
        type(flux_loop_t), allocatable :: loops(:)
        integer :: unit, ierr
        character(len=:), allocatable :: message

        open(newunit=unit, file=path, status='replace', action='write')
        write(unit, '(A)', advance='no') text
        close(unit)
        call read_flux_loop_file(path, loops, ierr, message)
        if (ierr /= 0) error stop 'eval_file: '//message
        call solver%flux_loops(loops, flux, rule)
    end subroutine eval_file

    subroutine check_close(name, actual, expected, rtol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual, expected, rtol

        if (abs(actual - expected) > rtol * max(abs(expected), tiny(1.0_dp)) &
                .or. actual /= actual) then
            write(error_unit, '(A,2ES24.16)') 'FAIL '//name//': ', actual, expected
            failures = failures + 1
        else
            print '(A)', 'ok   '//name
        end if
    end subroutine check_close
end program test_vacuum_semantics
