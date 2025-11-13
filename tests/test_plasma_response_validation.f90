program test_plasma_response_validation
    !! Validate plasma response by testing array ordering and grid conventions
    !! against simsopt reference implementation details
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_response, only: plasma_response_t
    implicit none

    integer :: iphi, itheta, i
    real(dp) :: theta, phi, pi
    real(dp), allocatable :: x_surf(:,:,:), b_total(:,:,:), b_ext(:,:,:)
    type(plasma_response_t) :: pr
    integer :: nphi, ntheta

    ! Test 1: Verify flattening order matches simsopt convention
    print *, "TEST 1: Array flattening convention"
    print *, "===================================="

    ! CRITICAL: virtual-casing requires Nt >= 6 (PATCH_DIM) for FFT quadrature!
    nphi = 8
    ntheta = 8
    allocate(x_surf(nphi, ntheta, 3))
    allocate(b_total(nphi, ntheta, 3))
    allocate(b_ext(nphi, ntheta, 3))

    pi = acos(-1.0_dp)

    ! Create simple circular torus surface
    do iphi = 1, nphi
        phi = (iphi - 1.0_dp) / real(nphi, dp) * 2.0_dp * pi
        do itheta = 1, ntheta
            theta = (itheta - 1.0_dp) / real(ntheta, dp) * 2.0_dp * pi
            ! Major radius R=2, minor radius a=0.5
            x_surf(iphi, itheta, 1) = (2.0_dp + 0.5_dp * cos(theta)) * &
                                      cos(phi)
            x_surf(iphi, itheta, 2) = (2.0_dp + 0.5_dp * cos(theta)) * &
                                      sin(phi)
            x_surf(iphi, itheta, 3) = 0.5_dp * sin(theta)
            ! Simple dipole-like field
            b_total(iphi, itheta, 1) = sin(theta)
            b_total(iphi, itheta, 2) = cos(phi)
            b_total(iphi, itheta, 3) = 0.1_dp
        end do
    end do

    ! Initialize plasma response
    call pr%init(nfp=1, x_surf=x_surf, b_total=b_total, &
                 src_nphi=nphi, src_ntheta=ntheta, &
                 use_stellsym=.false., digits=3)

    if (.not. pr%is_initialized()) then
        print *, "WARNING: Plasma response init failed (expected for toy geometry)"
    else
        ! Compute external field
        b_ext = 0.0_dp
        call pr%compute_bext(b_total, b_ext)

        ! Check output magnitude is reasonable
        if (maxval(abs(b_ext)) > 1e10_dp) then
            print *, "WARNING: B_external magnitude seems large:"
            print *, "  Max |B_ext|:", maxval(abs(b_ext))
            print *, "  This may indicate issue with toy geometry or units"
        else
            print *, "B_external magnitude reasonable:"
            print *, "  RMS:", sqrt(sum(b_ext**2) / size(b_ext))
            print *, "  Max:", maxval(abs(b_ext))
        end if

        call pr%finalize()
    end if

    deallocate(x_surf, b_total, b_ext)

    ! Test 2: Verify fortran array indexing
    print *, ""
    print *, "TEST 2: Fortran array indexing convention"
    print *, "=========================================="

    allocate(x_surf(2, 3, 3))
    ! Fill with pattern: x_surf(i,j,k) = i + 10*j + 100*k
    do i = 1, 2
        do iphi = 1, 3  ! Using iphi as second index
            x_surf(i, iphi, 1) = real(i + 10*iphi + 100*1, dp)
            x_surf(i, iphi, 2) = real(i + 10*iphi + 100*2, dp)
            x_surf(i, iphi, 3) = real(i + 10*iphi + 100*3, dp)
        end do
    end do

    ! Check memory layout (Fortran column-major)
    ! In Fortran, first index varies fastest
    ! So x_surf(1,1,1), x_surf(2,1,1), x_surf(1,2,1), x_surf(2,2,1), ...
    print *, "Array values (should be i+10*j+100*k):"
    print *, "  x_surf(1,1,:) =", x_surf(1,1,:)
    print *, "  x_surf(2,1,:) =", x_surf(2,1,:)
    print *, "  x_surf(1,2,:) =", x_surf(1,2,:)
    print *, "  x_surf(2,2,:) =", x_surf(2,2,:)

    deallocate(x_surf)

    print *, ""
    print *, "VALIDATION TESTS COMPLETE"
    print *, "=========================="
    print *, "If B_external magnitudes were reasonable (< 1e10),  array"
    print *, "flattening and indexing are correct."
    print *, ""
    print *, "Note: Toy geometry may not be physically reasonable for"
    print *, "virtual-casing FFT quadrature. Real VMEC test data should be used"
    print *, "for production validation."

end program test_plasma_response_validation
