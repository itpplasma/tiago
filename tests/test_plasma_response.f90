program test_plasma_response
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_response, only: plasma_response_t
    implicit none

    call test_initialization()
    call test_uninitialized_compute()
    call test_context_lifecycle()
    call test_unit_conversion_factors()
    call test_virtual_casing_offsurface()
    call test_vector_potential_interface()

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

    call build_mock_surface(nphi, ntheta, x_surf, b_total)

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

    call build_mock_surface(nphi, ntheta, x_surf, b_total)

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

subroutine test_virtual_casing_offsurface()
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_response, only: plasma_response_t
    implicit none

    type(plasma_response_t) :: pr
    real(dp), allocatable :: x_surf(:, :, :)
    real(dp), allocatable :: b_total(:, :, :)
    real(dp), allocatable :: points(:, :)
    real(dp), allocatable :: vc_field(:, :)
    real(dp), allocatable :: bs_field(:, :)
    integer :: nphi
    integer :: ntheta
    integer :: npts
    real(dp) :: tol

    nphi = 4
    ntheta = 4
    npts = 3
    allocate(x_surf(nphi, ntheta, 3))
    allocate(b_total(nphi, ntheta, 3))
    allocate(points(npts, 3))
    allocate(vc_field(npts, 3))
    allocate(bs_field(npts, 3))

    call build_mock_surface(nphi, ntheta, x_surf, b_total)
    b_total = 0.0_dp

    points(1, :) = [2.0_dp, 0.0_dp, 0.1_dp]
    points(2, :) = [-1.5_dp, 0.2_dp, 0.2_dp]
    points(3, :) = [0.5_dp, -1.8_dp, 0.3_dp]

    call pr%init(nfp=1, x_surf=x_surf, b_total=b_total, &
        src_nphi=nphi, src_ntheta=ntheta)

    call pr%compute_bext_points(b_total, points, vc_field)
    call pr%compute_surface_bext_points(points, bs_field)
    tol = 1.0e-6_dp
    if (maxval(abs(vc_field - bs_field)) > tol) then
        error stop 'virtual casing off-surface mismatch'
    end if

    call pr%finalize()
    deallocate(x_surf, b_total, points, vc_field, bs_field)
end subroutine test_virtual_casing_offsurface

subroutine test_vector_potential_interface()
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_response, only: plasma_response_t
    implicit none

    type(plasma_response_t) :: pr
    real(dp), allocatable :: x_surf(:, :, :), b_total(:, :, :)
    real(dp), allocatable :: b_ext(:, :, :)
    real(dp), allocatable :: points(:, :)
    real(dp), allocatable :: a_pts(:, :)
    real(dp), allocatable :: b_pts(:, :)
    integer :: nphi, ntheta
    integer :: npts
    integer :: iphi
    integer :: itheta
    integer :: idx

    nphi = 4
    ntheta = 4
    allocate(x_surf(nphi, ntheta, 3))
    allocate(b_total(nphi, ntheta, 3))
    allocate(b_ext(nphi, ntheta, 3))
    npts = nphi * ntheta
    allocate(points(npts, 3))
    allocate(a_pts(npts, 3))
    allocate(b_pts(npts, 3))

    call build_mock_surface(nphi, ntheta, x_surf, b_total)
    idx = 0
    do iphi = 1, nphi
        do itheta = 1, ntheta
            idx = idx + 1
            points(idx, :) = x_surf(iphi, itheta, :)
        end do
    end do

    call pr%init(nfp=1, x_surf=x_surf, b_total=b_total, &
        src_nphi=nphi, src_ntheta=ntheta)

    call pr%compute_vector_potential_points(points, a_pts)
    if (maxval(abs(a_pts)) > 1.0e-12_dp) then
        error stop 'zero plasma should yield zero vector potential'
    end if

    call pr%compute_surface_bext_points(points, b_pts)
    if (maxval(abs(b_pts)) > 1.0e-12_dp) then
        error stop 'zero plasma should yield zero Biot-Savart field'
    end if

    b_ext = 1.0_dp
    call pr%compute_bext(b_total, b_ext)
    if (maxval(abs(b_ext)) > 1.0e-12_dp) then
        error stop 'zero plasma should yield zero B_ext'
    end if

    call pr%finalize()

    deallocate(x_surf, b_total, b_ext, points, a_pts, b_pts)
end subroutine test_vector_potential_interface

subroutine build_mock_surface(nphi, ntheta, x_surf, b_total)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    integer, intent(in) :: nphi
    integer, intent(in) :: ntheta
    real(dp), intent(out) :: x_surf(nphi, ntheta, 3)
    real(dp), intent(out) :: b_total(nphi, ntheta, 3)
    integer :: iphi
    integer :: itheta
    real(dp) :: phi
    real(dp) :: theta
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: base_r = 1.0_dp
    real(dp), parameter :: minor_r = 0.1_dp

    do iphi = 1, nphi
        phi = 2.0_dp * pi * real(iphi - 1, dp) / real(nphi, dp)
        do itheta = 1, ntheta
            theta = 2.0_dp * pi * real(itheta - 1, dp) / real(ntheta, dp)
            x_surf(iphi, itheta, 1) = (base_r + minor_r * cos(theta)) * cos(phi)
            x_surf(iphi, itheta, 2) = (base_r + minor_r * cos(theta)) * sin(phi)
            x_surf(iphi, itheta, 3) = minor_r * sin(theta)
            b_total(iphi, itheta, 1) = 0.0_dp
            b_total(iphi, itheta, 2) = 0.0_dp
            b_total(iphi, itheta, 3) = 0.0_dp
        end do
    end do
end subroutine build_mock_surface
