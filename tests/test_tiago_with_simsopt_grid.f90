program test_tiago_with_simsopt_grid
    !! Test TIAGO plasma response using simsopt's exact grid points and B-field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_response, only: plasma_response_t
    implicit none

    type(plasma_response_t) :: pr
    real(dp), allocatable :: x_surf(:,:,:), b_total(:,:,:), b_ext(:,:,:)
    integer :: nphi, ntheta, iphi, itheta, iunit, nfp, nargs
    real(dp) :: b_ext_rms, expected_rms, rel_error
    character(len=512) :: gamma_file, b_file, arg

    ! Parse command-line arguments
    nargs = command_argument_count()
    if (nargs < 2) then
        print *, 'Usage: test_tiago_with_simsopt_grid <gamma.csv> <b_total.csv>'
        stop 1
    end if

    call get_command_argument(1, gamma_file)
    call get_command_argument(2, b_file)

    ! Load simsopt data from Python-generated files
    nphi = 16
    ntheta = 16
    nfp = 3

    allocate(x_surf(nphi, ntheta, 3))
    allocate(b_total(nphi, ntheta, 3))
    allocate(b_ext(nphi, ntheta, 3))

    ! Read simsopt gamma (surface positions) from CSV
    open(newunit=iunit, file=gamma_file, status='old')
    read(iunit, *)  ! Skip header
    do iphi = 1, nphi
        do itheta = 1, ntheta
            read(iunit, *) x_surf(iphi, itheta, 1), x_surf(iphi, itheta, 2), &
                          x_surf(iphi, itheta, 3)
        end do
    end do
    close(iunit)

    ! Read simsopt B_total from CSV
    open(newunit=iunit, file=b_file, status='old')
    read(iunit, *)  ! Skip header
    do iphi = 1, nphi
        do itheta = 1, ntheta
            read(iunit, *) b_total(iphi, itheta, 1), b_total(iphi, itheta, 2), &
                          b_total(iphi, itheta, 3)
        end do
    end do
    close(iunit)

    ! Initialize TIAGO plasma response
    call pr%init(nfp=nfp, x_surf=x_surf, b_total=b_total, &
        src_nphi=nphi, src_ntheta=ntheta, use_stellsym=.true., digits=6)

    if (.not. pr%is_initialized()) then
        print *, 'ERROR: init failed'
        stop 1
    end if

    ! Compute B_external
    b_ext = 0.0_dp
    call pr%compute_bext(b_total, b_ext)

    ! Compute RMS and validate against expected value
    b_ext_rms = sqrt(sum(b_ext**2) / real(size(b_ext), dp))
    expected_rms = 0.944d0  ! Simsopt reference value

    print '(A,F12.6,A)', 'TIAGO B_external RMS: ', b_ext_rms, ' Tesla'
    print '(A,F12.6,A)', 'Expected RMS:         ', expected_rms, ' Tesla'

    rel_error = abs(b_ext_rms - expected_rms) / expected_rms
    print '(A,F8.4,A)', 'Relative error:       ', rel_error * 100.0d0, ' %'

    ! Check tolerance: 2% relative error
    if (rel_error > 0.02d0) then
        print *, 'ERROR: B_external RMS exceeds 2% tolerance'
        stop 1
    end if

    print *, 'Test PASSED: B_external RMS within 2% of simsopt reference'

    call pr%finalize()
    deallocate(x_surf, b_total, b_ext)

end program test_tiago_with_simsopt_grid
