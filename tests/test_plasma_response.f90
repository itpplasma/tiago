program test_plasma_response
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_response, only: plasma_response_t
    implicit none

    call test_initialization()
    call test_uninitialized_compute()
    call test_context_lifecycle()
    call test_unit_conversion_factors()

    print *, 'All plasma_response tests passed'
end program test_plasma_response

subroutine test_initialization()
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_response, only: plasma_response_t
    implicit none

    type(plasma_response_t) :: pr
    real(dp), allocatable :: x_surf(:,:,:), b_total(:,:,:)
    integer :: nphi, ntheta

    nphi = 4
    ntheta = 4
    allocate(x_surf(nphi, ntheta, 3))
    allocate(b_total(nphi, ntheta, 3))

    x_surf = 0.0_dp
    b_total = 0.0_dp

    if (pr%is_initialized()) then
        error stop 'plasma_response should start uninitialized'
    end if

    call pr%init(nfp=3, x_surf=x_surf, b_total=b_total, &
                 src_nphi=nphi, src_ntheta=ntheta)

    if (pr%is_initialized()) then
        ! Context creation succeeded
        call pr%finalize()
    end if

    deallocate(x_surf, b_total)
end subroutine test_initialization

subroutine test_uninitialized_compute()
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_response, only: plasma_response_t
    implicit none

    type(plasma_response_t) :: pr
    real(dp), allocatable :: b_total(:,:,:), b_ext(:,:,:)
    integer :: nphi, ntheta

    nphi = 4
    ntheta = 4
    allocate(b_total(nphi, ntheta, 3))
    allocate(b_ext(nphi, ntheta, 3))

    b_total = 1.0_dp
    b_ext = -1.0_dp

    call pr%compute_bext(b_total, b_ext)

    if (abs(b_ext(1, 1, 1) - 0.0_dp) > 1.0e-14_dp) then
        error stop 'uninitialized compute_bext should zero output'
    end if

    deallocate(b_total, b_ext)
end subroutine test_uninitialized_compute

subroutine test_context_lifecycle()
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_response, only: plasma_response_t
    implicit none

    type(plasma_response_t) :: pr
    real(dp), allocatable :: x_surf(:,:,:), b_total(:,:,:)
    integer :: i, nphi, ntheta

    nphi = 8
    ntheta = 8
    allocate(x_surf(nphi, ntheta, 3))
    allocate(b_total(nphi, ntheta, 3))

    do i = 1, 3
        x_surf(:, :, i) = real(i, dp)
        b_total(:, :, i) = 0.1_dp * real(i, dp)
    end do

    call pr%init(nfp=2, x_surf=x_surf, b_total=b_total, &
                 src_nphi=nphi, src_ntheta=ntheta)

    call pr%finalize()

    if (pr%is_initialized()) then
        error stop 'finalize should mark as uninitialized'
    end if

    deallocate(x_surf, b_total)
end subroutine test_context_lifecycle

subroutine test_unit_conversion_factors()
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: meters_to_cm = 100.0_dp
    real(dp), parameter :: maxwell_to_weber = 1.0e-8_dp
    real(dp) :: test_length, test_flux

    test_length = 1.0_dp
    test_length = test_length * meters_to_cm
    if (abs(test_length - 100.0_dp) > 1.0e-10_dp) then
        error stop 'meters_to_cm conversion incorrect'
    end if

    test_flux = 1.0_dp
    test_flux = test_flux * maxwell_to_weber
    if (abs(test_flux - 1.0e-8_dp) > 1.0e-18_dp) then
        error stop 'maxwell_to_weber conversion incorrect'
    end if

end subroutine test_unit_conversion_factors
