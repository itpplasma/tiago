program test_tiago_with_simsopt_grid
    !! Test TIAGO plasma response using simsopt's exact grid points and B-field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_response, only: plasma_response_t
    implicit none

    type(plasma_response_t) :: pr, pr_off
    real(dp), allocatable :: x_surf(:,:,:), b_total(:,:,:), b_ext(:,:,:)
    real(dp), allocatable :: normals(:,:,:), b_ext_offset(:,:,:), x_eval(:,:,:)
    integer :: nphi, ntheta, iphi, itheta, iunit, nfp, nargs
    integer :: phi_stride, theta_stride, eval_nphi, eval_ntheta
    integer :: iphi_idx, itheta_idx
    real(dp) :: b_ext_rms, expected_rms, rel_error, offset_distance
    character(len=512) :: gamma_file, b_file, bext_file, bext_offset_file
    character(len=512) :: arg
    logical :: write_bext, write_bext_offset

    ! Parse command-line arguments
    nargs = command_argument_count()
    if (nargs < 2) then
        print *, 'Usage: test_tiago_with_simsopt_grid <gamma.csv> <b_total.csv> [b_ext.csv [b_ext_offset.csv offset_m]]'
        stop 1
    end if

    call get_command_argument(1, gamma_file)
    call get_command_argument(2, b_file)
    write_bext = .false.
    write_bext_offset = .false.
    offset_distance = 0.0_dp
    if (nargs >= 3) then
        call get_command_argument(3, bext_file)
        if (len_trim(bext_file) > 0) then
            write_bext = .true.
        end if
    end if
    if (nargs >= 5) then
        call get_command_argument(4, bext_offset_file)
        call get_command_argument(5, arg)
        if (len_trim(bext_offset_file) > 0 .and. len_trim(arg) > 0) then
            read(arg, *) offset_distance
            write_bext_offset = .true.
        end if
    end if

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

    if (write_bext) then
        print *, 'Writing TIAGO B_external samples to ', trim(bext_file)
        call write_vector_csv(trim(bext_file), b_ext)
    end if

    if (write_bext_offset) then
        if (offset_distance <= 0.0_dp) then
            print *, 'WARNING: offset distance <= 0; skipping off-surface export'
        else
            print *, 'Computing off-surface B_external (offset=', offset_distance, 'm )'
            allocate(normals(nphi, ntheta, 3))
            call compute_surface_normals(x_surf, normals)
            phi_stride = max(1, nphi / max(1, min(nphi, 8)))
            theta_stride = max(1, ntheta / max(1, min(ntheta, 8)))
            eval_nphi = max(1, nphi / phi_stride)
            eval_ntheta = max(1, ntheta / theta_stride)
            allocate(x_eval(eval_nphi, eval_ntheta, 3))
            do iphi = 1, eval_nphi
                do itheta = 1, eval_ntheta
                    iphi_idx = 1 + (iphi - 1) * phi_stride
                    itheta_idx = 1 + (itheta - 1) * theta_stride
                    x_eval(iphi, itheta, :) = x_surf(iphi_idx, itheta_idx, :) + &
                        offset_distance * normals(iphi_idx, itheta_idx, :)
                end do
            end do
            print *, 'Off-surface grid: ', eval_nphi, 'x', eval_ntheta
            allocate(b_ext_offset(eval_nphi, eval_ntheta, 3))
            call pr_off%init(nfp=nfp, x_surf=x_surf, b_total=b_total, &
                src_nphi=nphi, src_ntheta=ntheta, use_stellsym=.true., digits=3)
            call pr_off%compute_bext_at(b_total, x_eval, b_ext_offset)
            call pr_off%finalize()
            print *, 'Writing off-surface B_external to ', trim(bext_offset_file)
            call write_vector_csv(trim(bext_offset_file), b_ext_offset)
            deallocate(normals, x_eval, b_ext_offset)
        end if
    end if

    call pr%finalize()
    deallocate(x_surf, b_total, b_ext)

contains

    subroutine write_vector_csv(path, field)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: field(:,:,:)
        integer :: unit
        integer :: i, j

        open(newunit=unit, file=path, status='replace', action='write')
        write(unit, '(A)') 'iphi,itheta,Bx,By,Bz'
        do i = 1, size(field, 1)
            do j = 1, size(field, 2)
                write(unit, '(I0,",",I0,",",ES23.15E3,",",ES23.15E3,",",ES23.15E3)') i, j, &
                    field(i, j, 1), field(i, j, 2), field(i, j, 3)
            end do
        end do
        close(unit)
    end subroutine write_vector_csv

    subroutine compute_surface_normals(points, normals)
        real(dp), intent(in) :: points(:,:,:)
        real(dp), intent(out) :: normals(:,:,:)
        integer :: nphi, ntheta, iphi, itheta
        real(dp) :: dphi(3), dtheta(3), normal(3), norm_mag
        integer :: ip_next, ip_prev, it_next, it_prev

        nphi = size(points, 1)
        ntheta = size(points, 2)

        do iphi = 1, nphi
            ip_next = wrap_index(iphi + 1, nphi)
            ip_prev = wrap_index(iphi - 1, nphi)
            do itheta = 1, ntheta
                it_next = wrap_index(itheta + 1, ntheta)
                it_prev = wrap_index(itheta - 1, ntheta)

                dphi = points(ip_next, itheta, :) - points(ip_prev, itheta, :)
                dtheta = points(iphi, it_next, :) - points(iphi, it_prev, :)
                normal = cross_product(dphi, dtheta)
                norm_mag = sqrt(sum(normal**2))
                if (norm_mag > 0.0_dp) then
                    normals(iphi, itheta, :) = normal / norm_mag
                else
                    normals(iphi, itheta, :) = 0.0_dp
                end if
            end do
        end do
    end subroutine compute_surface_normals

    pure integer function wrap_index(idx, n) result(res)
        integer, intent(in) :: idx, n
        res = idx
        do while (res < 1)
            res = res + n
        end do
        do while (res > n)
            res = res - n
        end do
    end function wrap_index

    pure function cross_product(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)

        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end function cross_product

end program test_tiago_with_simsopt_grid
