program test_coil_loader
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic
    use neo_biotsavart_field, only: biotsavart_field_t
    use tiago_coil_loader, only: load_coils_into_field
    implicit none

    type(biotsavart_field_t) :: field
    character(len=:), allocatable :: data_dir
    character(len=512) :: arg

    ! Usage: tiago_coil_loader_tests <tests/data directory>
    call get_command_argument(1, arg)
    data_dir = trim(arg)

    call test_neo_input(field, data_dir)
    call test_stellopt_input(field, data_dir)
    call test_stellopt_with_extcur(field, data_dir)
    call test_more_points(field, data_dir)
    call test_extcur_forms(field)

    print *, 'tiago_coil_loader tests passed'
end program test_coil_loader

subroutine test_neo_input(field, data_dir)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_biotsavart_field, only: biotsavart_field_t
    use neo_biotsavart, only: coils_deinit
    use tiago_coil_loader, only: load_coils_into_field
    implicit none
    type(biotsavart_field_t), intent(inout) :: field
    character(len=*), intent(in) :: data_dir
    real(dp), parameter :: tol = 1.0e-12_dp

    call load_coils_into_field(field, data_dir//'/coil_sample.neo')
    call assert_equal(size(field%coils%x), 5, 'neo: point count')
    call assert_close(field%coils%x(1), 0.0_dp, tol, 'neo: x(1)')
    call assert_close(field%coils%y(3), 1.0_dp, tol, 'neo: y(3)')
    call assert_close(field%coils%current(2), 1.0_dp, tol, 'neo: current(2)')
    call coils_deinit(field%coils)
end subroutine test_neo_input

subroutine test_stellopt_input(field, data_dir)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_biotsavart_field, only: biotsavart_field_t
    use neo_biotsavart, only: coils_deinit
    use tiago_coil_loader, only: load_coils_into_field
    implicit none
    type(biotsavart_field_t), intent(inout) :: field
    character(len=*), intent(in) :: data_dir
    real(dp), parameter :: tol = 1.0e-12_dp

    call load_coils_into_field(field, data_dir//'/coils_sample.coils')
    call assert_equal(size(field%coils%x), 5, 'stellopt: point count')
    call assert_close(field%coils%current(1), 1.0_dp, tol, 'stellopt: current(1)')
    call assert_close(field%coils%current(5), 0.0_dp, tol, 'stellopt: closing current zero')
    call coils_deinit(field%coils)
end subroutine test_stellopt_input

subroutine test_stellopt_with_extcur(field, data_dir)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_biotsavart_field, only: biotsavart_field_t
    use neo_biotsavart, only: coils_deinit
    use tiago_coil_loader, only: load_coils_into_field
    implicit none
    type(biotsavart_field_t), intent(inout) :: field
    character(len=*), intent(in) :: data_dir
    real(dp), parameter :: tol = 1.0e-12_dp

    call load_coils_into_field(field, data_dir//'/coils_sample.coils', &
         extcur_path=data_dir//'/extcur_sample.diagno')
    call assert_close(field%coils%current(1), 2.0_dp, tol, 'stellopt extcur: current scaled')
    call assert_close(field%coils%current(5), 0.0_dp, tol, 'stellopt extcur: closing current stays zero')
    call coils_deinit(field%coils)
end subroutine test_stellopt_with_extcur

subroutine test_more_points(field, data_dir)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic
    use neo_biotsavart_field, only: biotsavart_field_t
    use neo_biotsavart, only: coils_deinit
    use tiago_coil_loader, only: load_coils_into_field
    implicit none
    type(biotsavart_field_t), intent(inout) :: field
    character(len=*), intent(in) :: data_dir
    integer :: npts

    call load_coils_into_field(field, data_dir//'/coils_more_points.coils', &
         extcur_path=data_dir//'/extcur_more_points.input')
    npts = size(field%coils%x)
    call assert_equal(npts, 10, 'more_points: node count')
    if (any(ieee_is_nan(field%coils%current))) then
        call assert_equal(0, 1, 'more_points: currents contain NaN')
    end if
    call coils_deinit(field%coils)
end subroutine test_more_points

subroutine test_extcur_forms(field)
    !! Group A starts at 1 A, group B at 2 A. current(1) belongs to A, current(6) to B.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_biotsavart_field, only: biotsavart_field_t
    use neo_biotsavart, only: coils_deinit
    use tiago_coil_loader, only: load_coils_into_field
    implicit none
    type(biotsavart_field_t), intent(inout) :: field
    real(dp), parameter :: tol = 1.0e-12_dp
    integer :: unit

    open(newunit=unit, file='extcur_two_groups.coils', status='replace', action='write')
    write(unit, '(A)') 'periods 1', 'begin filament', 'mirror NIL', &
        ' 0 0 0 1', ' 1 0 0 1', ' 1 1 0 1', ' 0 1 0 1', ' 0 0 0 0 1 A', &
        ' 0 0 1 2', ' 1 0 1 2', ' 1 1 1 2', ' 0 1 1 2', ' 0 0 1 0 2 B', 'end'
    close(unit)

    call check('&INDATA'//new_line('a')//' EXTCUR(1) = 3.0'//new_line('a')//' EXTCUR(2) = 0.0'// &
        new_line('a')//'/', 3.0_dp, 0.0_dp, 'zero switches a group off')
    call check('&INDATA  EXTCUR = 3.0, 5.0 /', 3.0_dp, 5.0_dp, 'namelist array')
    call check('&INDATA  EXTCUR(1) = 3.0  EXTCUR(2) = 5.0D0 /', 3.0_dp, 5.0_dp, 'several per line')
    call check('&INDATA  EXTCUR = 2*4.0 /', 4.0_dp, 4.0_dp, 'repeat count')
    call check('&INDATA  EXTCUR(1:2) = 6.0 7.0 /', 6.0_dp, 7.0_dp, 'slice')
    call check('&INDATA  EXTCUR(2:) = 8.0 /', 0.0_dp, 8.0_dp, 'open slice; group 1 off')
    call check('&INDATA  LEXTCUR = 9  EXTCUR(2) = 5.0 ! EXTCUR(1)=7'//new_line('a')//'/', &
        0.0_dp, 5.0_dp, 'only EXTCUR(2) given (group 1 off); comments and LEXTCUR ignored')
    call check('3.0'//new_line('a')//'5.0', 3.0_dp, 5.0_dp, 'plain list')

contains

    subroutine check(text, expect_a, expect_b, label)
        character(len=*), intent(in) :: text, label
        real(dp), intent(in) :: expect_a, expect_b

        open(newunit=unit, file='extcur_case.input', status='replace', action='write')
        write(unit, '(A)') text
        close(unit)
        call load_coils_into_field(field, 'extcur_two_groups.coils', extcur_path='extcur_case.input')
        call assert_close(field%coils%current(1), expect_a, tol, 'extcur '//label//' (A)')
        call assert_close(field%coils%current(6), expect_b, tol, 'extcur '//label//' (B)')
        call coils_deinit(field%coils)
    end subroutine check
end subroutine test_extcur_forms

subroutine assert_equal(actual, expected, label)
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none
    integer, intent(in) :: actual, expected
    character(len=*), intent(in) :: label
    if (actual /= expected) then
        write(error_unit, '(A,": expected=",I0,", got=",I0)') trim(label), expected, actual
        stop 1
    end if
end subroutine assert_equal

subroutine assert_close(actual, expected, tol, label)
    use, intrinsic :: iso_fortran_env, only: error_unit, dp => real64
    implicit none
    real(dp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > tol) then
        write(error_unit, '(A,": expected=",ES12.5,", got=",ES12.5,", tol=",ES12.5)') &
            trim(label), expected, actual, tol
        stop 1
    end if
end subroutine assert_close
