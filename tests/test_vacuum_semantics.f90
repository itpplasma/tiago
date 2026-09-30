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
    call test_one_period_loop()

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

    subroutine test_one_period_loop()
        !! iflflg=1: a loop over one field period, closed to its first point
        !! rotated by 2*pi/nfp, times nfp. For an nfp-symmetric coil set this is
        !! exactly the full toroidal polygon.
        integer, parameter :: nfp = 3, n = 12
        real(dp), parameter :: radius = 0.9_dp, pi = acos(-1.0_dp)
        type(vacuum_solver_t) :: sym
        character(len=:), allocatable :: text
        real(dp), allocatable :: flux(:)
        real(dp) :: phi
        integer :: k

        call write_symmetric_coils('semantics_sym.coils', nfp)
        call sym%init('semantics_sym.coils')
        call sym%set_nfp(nfp)

        text = '2' // new_line('a') // itoa(n) // ' 1 0 PERIOD' // new_line('a')
        do k = 0, n - 1
            phi = 2.0_dp * pi * k / real(nfp * n, dp)
            text = text // xyz(radius * cos(phi), radius * sin(phi), 0.05_dp)
        end do
        text = text // itoa(nfp * n + 1) // ' 0 0 FULL' // new_line('a')
        do k = 0, nfp * n
            phi = 2.0_dp * pi * k / real(nfp * n, dp)
            text = text // xyz(radius * cos(phi), radius * sin(phi), 0.05_dp)
        end do

        call eval_file('semantics_period.diagno', text, flux, sym)
        call check_close('iflflg=1 period loop equals full toroidal loop', flux(1), &
            flux(2), 1.0e-10_dp)
        call sym%finalize()
    end subroutine test_one_period_loop

    subroutine write_symmetric_coils(path, nfp)
        !! nfp tilted square coils, rotated copies of each other.
        character(len=*), intent(in) :: path
        integer, intent(in) :: nfp
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: base(3, 5), c, s
        integer :: unit, ip, k

        base = reshape([1.0_dp, 0.1_dp, -0.2_dp,  1.4_dp, 0.1_dp, -0.2_dp, &
                        1.4_dp, 0.3_dp,  0.2_dp,  1.0_dp, 0.3_dp,  0.2_dp, &
                        1.0_dp, 0.1_dp, -0.2_dp], [3, 5])
        open(newunit=unit, file=path, status='replace', action='write')
        write(unit, '(A)') 'periods 1', 'begin filament', 'mirror NIL'
        do ip = 0, nfp - 1
            c = cos(2.0_dp * pi * ip / nfp)
            s = sin(2.0_dp * pi * ip / nfp)
            do k = 1, 4
                write(unit, '(4ES24.15)') c * base(1, k) - s * base(2, k), &
                    s * base(1, k) + c * base(2, k), base(3, k), 1.0e3_dp
            end do
            write(unit, '(4ES24.15,A)') c * base(1, 5) - s * base(2, 5), &
                s * base(1, 5) + c * base(2, 5), base(3, 5), 0.0_dp, ' 1 COIL'
        end do
        write(unit, '(A)') 'end'
        close(unit)
    end subroutine write_symmetric_coils

    function itoa(i) result(text)
        integer, intent(in) :: i
        character(len=:), allocatable :: text
        character(len=16) :: buf
        write(buf, '(I0)') i
        text = trim(buf)
    end function itoa

    function xyz(x, y, z) result(text)
        real(dp), intent(in) :: x, y, z
        character(len=:), allocatable :: text
        character(len=80) :: buf
        write(buf, '(3ES25.16)') x, y, z
        text = trim(buf) // new_line('a')
    end function xyz

    subroutine eval_file(path, text, flux, alt_solver)
        character(len=*), intent(in) :: path, text
        real(dp), allocatable, intent(out) :: flux(:)
        type(vacuum_solver_t), intent(in), optional :: alt_solver
        type(flux_loop_t), allocatable :: loops(:)
        integer :: unit, ierr
        character(len=:), allocatable :: message

        open(newunit=unit, file=path, status='replace', action='write')
        write(unit, '(A)', advance='no') text
        close(unit)
        call read_flux_loop_file(path, loops, ierr, message)
        if (ierr /= 0) error stop 'eval_file: '//message
        if (present(alt_solver)) then
            call alt_solver%flux_loops(loops, flux, rule)
        else
            call solver%flux_loops(loops, flux, rule)
        end if
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
