program test_tiago_with_simsopt_grid
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_response, only: plasma_response_t
    implicit none

    type(plasma_response_t) :: response
    real(dp), allocatable :: x_surf(:, :, :), b_total(:, :, :), b_ext(:, :, :)
    real(dp), allocatable :: ref_ext(:, :, :)
    integer :: nphi, ntheta, iphi, itheta, unit, nargs
    character(len=512) :: gamma_file, btotal_file, output_file, ref_file
    logical :: check_reference
    real(dp) :: rel_error, rms_diff

    nargs = command_argument_count()
    if (nargs < 3) then
        print *, 'Usage: test_tiago_with_simsopt_grid <gamma.csv> <b_total.csv> <tiago_bext.csv> [simsopt_bext.csv]'
        stop 1
    end if

    call get_command_argument(1, gamma_file)
    call get_command_argument(2, btotal_file)
    call get_command_argument(3, output_file)
    check_reference = .false.
    if (nargs >= 4) then
        call get_command_argument(4, ref_file)
        if (len_trim(ref_file) > 0) check_reference = .true.
    end if

    nphi = 16
    ntheta = 16
    allocate(x_surf(nphi, ntheta, 3))
    allocate(b_total(nphi, ntheta, 3))
    allocate(b_ext(nphi, ntheta, 3))

    open(newunit=unit, file=gamma_file, status='old')
    read(unit, *)
    do iphi = 1, nphi
        do itheta = 1, ntheta
            read(unit, *) x_surf(iphi, itheta, 1), x_surf(iphi, itheta, 2), x_surf(iphi, itheta, 3)
        end do
    end do
    close(unit)

    open(newunit=unit, file=btotal_file, status='old')
    read(unit, *)
    do iphi = 1, nphi
        do itheta = 1, ntheta
            read(unit, *) b_total(iphi, itheta, 1), b_total(iphi, itheta, 2), b_total(iphi, itheta, 3)
        end do
    end do
    close(unit)

    call response%init(nfp=3, x_surf=x_surf, b_total=b_total, src_nphi=nphi, src_ntheta=ntheta)
    if (.not. response%is_initialized()) stop 'plasma response init failed'

    call response%compute_bext(b_total, b_ext)
    call write_vector_csv(output_file, b_ext, x_surf)

    if (check_reference) then
        allocate(ref_ext(nphi, ntheta, 3))
        call read_vector_csv(ref_file, ref_ext)
        ! ref_ext is simsopt B_external (field from coils, currents outside plasma)
        ! b_ext is TIAGO B_external (should match)

        rms_diff = sqrt(sum((b_ext - ref_ext)**2) / real(size(b_ext), dp))
        rel_error = rms_diff / max(1.0e-12_dp, sqrt(sum(ref_ext**2) / real(size(ref_ext), dp)))
        print '(A,1pe12.5)', 'RMS absolute error (TIAGO vs simsopt B_external) = ', rms_diff
        print '(A,1pe12.5)', 'Relative RMS error                              = ', rel_error
        print '(A,1pe12.5)', 'TIAGO B_external RMS                            = ', sqrt(sum(b_ext**2) / real(size(b_ext), dp))
        print '(A,1pe12.5)', 'Simsopt B_external RMS                          = ', sqrt(sum(ref_ext**2) / real(size(ref_ext), dp))
        if (rel_error > 0.02_dp) then
            print *, 'ERROR: Virtual casing B_external deviates more than 2% from simsopt reference'
            print *, 'STRICT TEST FAILED - FIX THE IMPLEMENTATION!'
            error stop 1
        else
            print *, 'SUCCESS: Virtual casing B_external matches simsopt reference within 2% tolerance'
        end if
        deallocate(ref_ext)
    end if

    call response%finalize()
    deallocate(x_surf, b_total, b_ext)

contains

    subroutine read_vector_csv(path, field)
        character(len=*), intent(in) :: path
        real(dp), intent(out) :: field(:, :, :)
        integer :: unit, iphi, itheta
        integer :: idx_phi, idx_theta
        character(len=256) :: header
        logical :: has_indices

        open(newunit=unit, file=path, status='old')
        read(unit, '(A)') header
        has_indices = index(header, 'iphi') > 0

        do iphi = 1, size(field, 1)
            do itheta = 1, size(field, 2)
                if (has_indices) then
                    read(unit, *) idx_phi, idx_theta, field(iphi, itheta, 1), &
                        field(iphi, itheta, 2), field(iphi, itheta, 3)
                else
                    read(unit, *) field(iphi, itheta, 1), field(iphi, itheta, 2), &
                        field(iphi, itheta, 3)
                end if
            end do
        end do
        close(unit)
    end subroutine read_vector_csv

    subroutine write_vector_csv(path, field, points)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: field(:, :, :)
        real(dp), intent(in) :: points(:, :, :)
        integer :: unit, iphi, itheta
        open(newunit=unit, file=path, status='replace')
        write(unit, '(A)') 'iphi,itheta,X,Y,Z,Bx,By,Bz'
        do iphi = 1, size(field, 1)
            do itheta = 1, size(field, 2)
                write(unit,'(I0,",",I0,",",ES23.15E3,",",ES23.15E3,",",ES23.15E3,",",ES23.15E3,",",ES23.15E3,",",ES23.15E3)') &
                    iphi, itheta, points(iphi, itheta, 1), points(iphi, itheta, 2), points(iphi, itheta, 3), &
                    field(iphi, itheta, 1), field(iphi, itheta, 2), field(iphi, itheta, 3)
            end do
        end do
        close(unit)
    end subroutine write_vector_csv

end program test_tiago_with_simsopt_grid
