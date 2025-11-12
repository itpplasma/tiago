program test_coil_loader
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic
    use neo_biotsavart_field, only: biotsavart_field_t
    use tiago_coil_loader, only: load_coils_into_field
    implicit none

    type(biotsavart_field_t) :: field

    call test_neo_input(field)
    call test_stellopt_input(field)
    call test_stellopt_with_extcur(field)
    call test_more_points(field)

    print *, 'tiago_coil_loader tests passed'
end program test_coil_loader

subroutine test_neo_input(field)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_biotsavart_field, only: biotsavart_field_t
    use neo_biotsavart, only: coils_deinit
    use tiago_coil_loader, only: load_coils_into_field
    implicit none
    type(biotsavart_field_t), intent(inout) :: field
    real(dp), parameter :: tol = 1.0e-12_dp

    call load_coils_into_field(field, '../tests/data/coil_sample.neo')
    call assert_equal(size(field%coils%x), 5, 'neo: point count')
    call assert_close(field%coils%x(1), 0.0_dp, tol, 'neo: x(1)')
    call assert_close(field%coils%y(3), 1.0_dp, tol, 'neo: y(3)')
    call assert_close(field%coils%current(2), 1.0_dp, tol, 'neo: current(2)')
    call coils_deinit(field%coils)
end subroutine test_neo_input

subroutine test_stellopt_input(field)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_biotsavart_field, only: biotsavart_field_t
    use neo_biotsavart, only: coils_deinit
    use tiago_coil_loader, only: load_coils_into_field
    implicit none
    type(biotsavart_field_t), intent(inout) :: field
    real(dp), parameter :: tol = 1.0e-12_dp

    call load_coils_into_field(field, '../tests/data/coils_sample.stellopt')
    call assert_equal(size(field%coils%x), 5, 'stellopt: point count')
    call assert_close(field%coils%current(1), 1.0_dp, tol, 'stellopt: current(1)')
    call assert_close(field%coils%current(5), 0.0_dp, tol, 'stellopt: closing current zero')
    call coils_deinit(field%coils)
end subroutine test_stellopt_input

subroutine test_stellopt_with_extcur(field)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_biotsavart_field, only: biotsavart_field_t
    use neo_biotsavart, only: coils_deinit
    use tiago_coil_loader, only: load_coils_into_field
    implicit none
    type(biotsavart_field_t), intent(inout) :: field
    real(dp), parameter :: tol = 1.0e-12_dp

    call load_coils_into_field(field, '../tests/data/coils_sample.stellopt', &
         extcur_path='../tests/data/extcur_sample.diagno')
    call assert_close(field%coils%current(1), 2.0_dp, tol, 'stellopt extcur: current scaled')
    call assert_close(field%coils%current(5), 0.0_dp, tol, 'stellopt extcur: closing current stays zero')
    call coils_deinit(field%coils)
end subroutine test_stellopt_with_extcur

subroutine test_more_points(field)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic
    use neo_biotsavart_field, only: biotsavart_field_t
    use neo_biotsavart, only: coils_deinit
    use tiago_coil_loader, only: load_coils_into_field
    implicit none
    type(biotsavart_field_t), intent(inout) :: field
    integer :: npts

    call load_coils_into_field(field, '../tests/data/coils_more_points.coils', &
         extcur_path='../tests/data/extcur_more_points.input')
    npts = size(field%coils%x)
    call assert_equal(npts, 10, 'more_points: node count')
    if (any(ieee_is_nan(field%coils%current))) then
        call assert_equal(0, 1, 'more_points: currents contain NaN')
    end if
    call coils_deinit(field%coils)
end subroutine test_more_points

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
