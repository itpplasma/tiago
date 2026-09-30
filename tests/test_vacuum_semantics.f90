program test_vacuum_semantics
    !! Regression tests for DIAGNO flux-loop semantics (no external codes needed).
    !! Usage: test_vacuum_semantics <coils_sample.coils>
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_diagnostic_types, only: flux_loop_t
    use tiago_flux_loops, only: read_flux_loop_file, finalize_flux_signals, lint_flux_loops
    use tiago_segmented_rogowski, only: read_segmented_rogowski_file
    use tiago_diagnostic_types, only: segmented_rogowski_t, bprobe_t
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
    call test_idia()
    call test_duplicate_coil_point()
    call test_file_formats()
    call test_bprobe_square_center()
    call test_response_reconstructs_signals()
    call test_gauss_quadrature()

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

    subroutine test_idia()
        !! idia=1 (diamagnetic, adds phiedge = 0 in vacuum) leaves the signal
        !! unchanged; idia=-k subtracts the flux of loop k.
        character(len=*), parameter :: big_rows = &
            '0.3 0.3 0.3' // new_line('a') // '0.7 0.3 0.3' // new_line('a') // &
            '0.7 0.7 0.3' // new_line('a') // '0.3 0.7 0.3' // new_line('a')
        character(len=*), parameter :: text = '4' // new_line('a') // &
            '4 0 0 SMALL' // new_line('a') // square_rows // &
            '4 0 1 SMALL_DIA' // new_line('a') // square_rows // &
            '4 0 0 BIG' // new_line('a') // big_rows // &
            '4 0 -1 BIG_MINUS_SMALL' // new_line('a') // big_rows
        character(len=*), parameter :: self_ref = '1' // new_line('a') // &
            '4 0 -1 SELF' // new_line('a') // square_rows
        type(flux_loop_t), allocatable :: loops(:)
        real(dp), allocatable :: flux(:)
        integer :: unit, ierr
        character(len=:), allocatable :: message

        call eval_file('semantics_idia.diagno', text, flux)
        call check_close('idia=1 unchanged in vacuum', flux(2), flux(1), 1.0e-14_dp)
        call check_close('idia=-1 subtracts loop 1', flux(4), flux(3) - flux(1), 1.0e-14_dp)

        open(newunit=unit, file='semantics_self.diagno', status='replace', action='write')
        write(unit, '(A)', advance='no') self_ref
        close(unit)
        call read_flux_loop_file('semantics_self.diagno', loops, ierr, message)
        if (ierr == 0) then
            write(error_unit, '(A)') 'FAIL idia self-reference must be rejected'
            failures = failures + 1
        else
            print '(A)', 'ok   idia self-reference rejected'
        end if
    end subroutine test_idia

    subroutine test_duplicate_coil_point()
        !! A repeated coil node is a zero-length segment and must contribute nothing.
        character(len=*), parameter :: text = '1' // new_line('a') // &
            '4 0 0 SQ' // new_line('a') // square_rows
        type(vacuum_solver_t) :: dup
        real(dp), allocatable :: ref(:), flux(:)
        integer :: unit

        open(newunit=unit, file='semantics_dup.coils', status='replace', action='write')
        write(unit, '(A)') 'periods 1', 'begin filament', 'mirror NIL', &
            ' 0 0 0 1', ' 1 0 0 1', ' 1 0 0 1', ' 1 1 0 1', ' 0 1 0 1', ' 0 0 0 0 1 SQ', 'end'
        close(unit)
        call dup%init('semantics_dup.coils')
        call eval_file('semantics_dup.diagno', text, ref)
        call eval_file('semantics_dup.diagno', text, flux, dup)
        call check_close('duplicated coil point is harmless', flux(1), ref(1), 1.0e-14_dp)
        call dup%finalize()
    end subroutine test_duplicate_coil_point

    subroutine test_file_formats()
        !! DIAGNO (3I6,A48) headers, labels with blanks, per-point segrog areas,
        !! trailing data and duplicate labels.
        character(len=*), parameter :: fixed = '     2' // new_line('a') // &
            '     4     0     0 Loop A/upper, 1' // new_line('a') // square_rows // &
            '4 0 0 free label' // new_line('a') // square_rows
        character(len=*), parameter :: trailing = '1' // new_line('a') // &
            '4 0 0 ONE' // new_line('a') // square_rows // '4 0 0 EXTRA' // new_line('a')
        character(len=*), parameter :: dup = '2' // new_line('a') // &
            '4 0 0 SAME' // new_line('a') // square_rows // '4 0 0 SAME' // new_line('a') // square_rows
        character(len=*), parameter :: seg = '     1' // new_line('a') // &
            '     3     0     0 SEG' // new_line('a') // &
            '0 0 0 1.0e-4' // new_line('a') // '0 0 1 5.0e-4' // new_line('a') // '0 0 2 0.0'
        type(flux_loop_t), allocatable :: loops(:)
        type(segmented_rogowski_t), allocatable :: segs(:)
        integer :: ierr
        character(len=:), allocatable :: message, report

        call write_text('formats_fixed.diagno', fixed)
        call read_flux_loop_file('formats_fixed.diagno', loops, ierr, message)
        call check_true('fixed-width header keeps full label', ierr == 0 .and. &
            loops(1)%label == 'Loop A/upper, 1' .and. loops(2)%label == 'free label')

        call write_text('formats_trailing.diagno', trailing)
        call read_flux_loop_file('formats_trailing.diagno', loops, ierr, message)
        call check_true('data beyond the header count is rejected', ierr /= 0)

        call write_text('formats_dup.diagno', dup)
        call read_flux_loop_file('formats_dup.diagno', loops, ierr, message)
        call check_true('duplicate labels fail lint', .not. lint_flux_loops(loops, report))

        call write_text('formats_seg.diagno', seg)
        call read_segmented_rogowski_file('formats_seg.diagno', segs, ierr, message, 9.0_dp)
        call check_true('per-point eff_area overrides the default', ierr == 0 .and. &
            all(segs(1)%segment_area == [1.0e-4_dp, 5.0e-4_dp]))
    end subroutine test_file_formats

    subroutine write_text(path, text)
        character(len=*), intent(in) :: path, text
        integer :: unit

        open(newunit=unit, file=path, status='replace', action='write')
        write(unit, '(A)') text
        close(unit)
    end subroutine write_text

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name
        logical, intent(in) :: condition

        if (condition) then
            print '(A)', 'ok   '//name
        else
            write(error_unit, '(A)') 'FAIL '//name
            failures = failures + 1
        end if
    end subroutine check_true

    subroutine test_bprobe_square_center()
        !! Field at the centre of a square loop (side a, current I):
        !! B_z = 2 sqrt(2) mu0 I / (pi a). coils_sample.coils: a = 1 m, I = 1 A.
        type(bprobe_t) :: probes(1)
        real(dp), allocatable :: signal(:)
        real(dp), parameter :: pi = acos(-1.0_dp)

        probes(1)%label = 'CENTER'
        probes(1)%position = [0.5_dp, 0.5_dp, 0.0_dp]
        probes(1)%normal = [0.0_dp, 0.0_dp, 1.0_dp]
        call solver%bprobes(probes, signal)
        call check_close('B-probe at square-loop centre', signal(1), &
            2.0_dp * sqrt(2.0_dp) * 4.0e-7_dp * pi / pi, 1.0e-12_dp)
    end subroutine test_bprobe_square_center

    subroutine test_response_reconstructs_signals()
        !! Two groups with EXTCUR 3 and -2: sum_g M_g EXTCUR_g equals the signals.
        type(vacuum_solver_t) :: two
        type(flux_loop_t), allocatable :: loops(:)
        type(segmented_rogowski_t), allocatable :: segs(:)
        type(bprobe_t), allocatable :: probes(:)
        real(dp), allocatable :: flux(:), fr(:, :), sr(:, :), pr(:, :)
        integer :: ierr
        character(len=:), allocatable :: message

        call write_text('response.coils', 'periods 1' // new_line('a') // 'begin filament' // &
            new_line('a') // 'mirror NIL' // new_line('a') // &
            ' 0 0 0 1' // new_line('a') // ' 1 0 0 1' // new_line('a') // ' 1 1 0 1' // &
            new_line('a') // ' 0 1 0 1' // new_line('a') // ' 0 0 0 0 1 A' // new_line('a') // &
            ' 0 0 0.1 5' // new_line('a') // ' 1 0 0.1 5' // new_line('a') // ' 1 1 0.1 5' // &
            new_line('a') // ' 0 1 0.1 5' // new_line('a') // ' 0 0 0.1 0 2 B' // &
            new_line('a') // 'end')
        call write_text('response.extcur', '&INDATA EXTCUR = 3.0, -2.0 /')
        call two%init('response.coils', 'response.extcur')
        call write_text('response.diagno', '1' // new_line('a') // '4 0 0 SQ' // new_line('a') // &
            square_rows)
        call read_flux_loop_file('response.diagno', loops, ierr, message)
        call two%flux_loops(loops, flux, rule)
        call two%response(loops, segs, probes, fr, sr, pr, rule)
        call check_close('response: sum_g M_g EXTCUR_g = signal', &
            3.0_dp * fr(1, 1) - 2.0_dp * fr(1, 2), flux(1), 1.0e-12_dp)
        call two%finalize()
    end subroutine test_response_reconstructs_signals

    subroutine test_gauss_quadrature()
        !! Gauss-Legendre converges spectrally: 8 points per segment agree with 32
        !! (midpoint needs thousands of samples for the same accuracy).
        type(flux_loop_t), allocatable :: loops(:)
        type(quadrature_rule_t) :: g8, g32
        real(dp), allocatable :: f8(:), f32(:)
        integer :: ierr
        character(len=:), allocatable :: message

        call write_text('gauss.diagno', '1' // new_line('a') // '4 0 0 SQ' // new_line('a') // &
            '0.2 0.2 0.1' // new_line('a') // '0.8 0.2 0.1' // new_line('a') // &
            '0.8 0.8 0.1' // new_line('a') // '0.2 0.8 0.1')
        call read_flux_loop_file('gauss.diagno', loops, ierr, message)
        g8%samples_per_segment = 8
        g8%gauss = .true.
        g32%samples_per_segment = 32
        g32%gauss = .true.
        call solver%flux_loops(loops, f8, g8)
        call solver%flux_loops(loops, f32, g32)
        call check_close('8-point Gauss = 32-point Gauss', f8(1), f32(1), 1.0e-9_dp)
    end subroutine test_gauss_quadrature

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
        call finalize_flux_signals(loops, flux)
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
