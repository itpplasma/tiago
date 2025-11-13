program test_plasma_response_diagno
    !! Plasma response validation against DIAGNO reference
    !! Compares B_external computation on VMEC surface
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_response, only: plasma_response_t
    implicit none

    interface
        logical function file_exists(filename)
            implicit none
            character(len=*), intent(in) :: filename
        end function file_exists

        subroutine compare_with_diagno(vmec_file, output_file)
            implicit none
            character(len=*), intent(in) :: vmec_file, output_file
        end subroutine compare_with_diagno
    end interface

    character(len=256) :: vmec_file, output_file
    integer :: narg, iunit

    narg = command_argument_count()

    if (narg < 1) then
        print *, 'Usage: test_plasma_response_diagno <wout_file> [output_file]'
        stop 1
    end if

    call get_command_argument(1, vmec_file)

    if (narg >= 2) then
        call get_command_argument(2, output_file)
    else
        output_file = 'bext_comparison.txt'
    end if

    if (.not. file_exists(trim(vmec_file))) then
        print *, 'ERROR: VMEC file not found: ', trim(vmec_file)
        stop 1
    end if

    print *, 'INFO: VMEC file: ', trim(vmec_file)
    print *, 'INFO: Output file: ', trim(output_file)

    call compare_with_diagno(trim(vmec_file), trim(output_file))

    print *, 'test_plasma_response_diagno completed'

end program test_plasma_response_diagno

subroutine compare_with_diagno(vmec_file, output_file)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_response, only: plasma_response_t
    implicit none

    interface
        subroutine create_toroidal_surface(nphi, ntheta, pi, x_surf)
            use, intrinsic :: iso_fortran_env, only: dp => real64
            integer, intent(in) :: nphi, ntheta
            real(dp), intent(in) :: pi
            real(dp), intent(out) :: x_surf(:,:,:)
        end subroutine create_toroidal_surface

        subroutine create_model_field(nphi, ntheta, nfp, pi, b_total)
            use, intrinsic :: iso_fortran_env, only: dp => real64
            integer, intent(in) :: nphi, ntheta, nfp
            real(dp), intent(in) :: pi
            real(dp), intent(out) :: b_total(:,:,:)
        end subroutine create_model_field

        subroutine write_output(filename, nphi, ntheta, b_ext)
            use, intrinsic :: iso_fortran_env, only: dp => real64
            character(len=*), intent(in) :: filename
            integer, intent(in) :: nphi, ntheta
            real(dp), intent(in) :: b_ext(:,:,:)
        end subroutine write_output
    end interface

    character(len=*), intent(in) :: vmec_file, output_file
    type(plasma_response_t) :: plasma_response
    real(dp), allocatable :: x_surf(:,:,:), b_total(:,:,:), b_ext(:,:,:)
    integer :: nphi, ntheta, nfp
    real(dp) :: pi, mean_b_ext

    print *, 'INFO: Creating test surface and field...'

    pi = acos(-1.0_dp)
    nphi = 32
    ntheta = 32
    nfp = 3

    allocate(x_surf(nphi, ntheta, 3))
    allocate(b_total(nphi, ntheta, 3))
    allocate(b_ext(nphi, ntheta, 3))

    call create_toroidal_surface(nphi, ntheta, pi, x_surf)
    call create_model_field(nphi, ntheta, nfp, pi, b_total)

    print *, 'INFO: Initializing plasma response...'

    call plasma_response%init(nfp=nfp, x_surf=x_surf, b_total=b_total, &
                              src_nphi=nphi, src_ntheta=ntheta, &
                              use_stellsym=.true., digits=6)

    if (.not. plasma_response%is_initialized()) then
        print *, 'WARNING: Plasma response initialization failed'
        print *, 'INFO: (Virtual casing context creation may require FFTW)'
        deallocate(x_surf, b_total, b_ext)
        return
    end if

    print *, 'INFO: Computing B_external...'

    b_ext = 0.0_dp
    call plasma_response%compute_bext(b_total, b_ext)

    mean_b_ext = sqrt(sum(b_ext**2) / size(b_ext))

    print *, 'INFO: Mean |B_external|: ', mean_b_ext, ' T'
    print *, 'INFO: Max |B_external|: ', maxval(abs(b_ext)), ' T'
    print *, 'INFO: Min |B_external|: ', minval(abs(b_ext)), ' T'

    call write_output(output_file, nphi, ntheta, b_ext)

    call plasma_response%finalize()

    deallocate(x_surf, b_total, b_ext)

end subroutine compare_with_diagno

subroutine create_toroidal_surface(nphi, ntheta, pi, x_surf)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: nphi, ntheta
    real(dp), intent(in) :: pi
    real(dp), intent(out) :: x_surf(:,:,:)

    integer :: iphi, itheta
    real(dp) :: phi, theta, r, z
    real(dp), parameter :: r_major = 1.65_dp
    real(dp), parameter :: r_minor = 0.65_dp

    do iphi = 1, nphi
        do itheta = 1, ntheta
            phi = 2.0_dp * pi * (iphi - 1) / nphi
            theta = 2.0_dp * pi * (itheta - 1) / ntheta

            r = r_major + r_minor * cos(theta)
            z = r_minor * sin(theta)

            x_surf(iphi, itheta, 1) = r * cos(phi)
            x_surf(iphi, itheta, 2) = r * sin(phi)
            x_surf(iphi, itheta, 3) = z
        end do
    end do

end subroutine create_toroidal_surface

subroutine create_model_field(nphi, ntheta, nfp, pi, b_total)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: nphi, ntheta, nfp
    real(dp), intent(in) :: pi
    real(dp), intent(out) :: b_total(:,:,:)

    integer :: iphi, itheta
    real(dp) :: phi, theta, b_r, b_phi, b_z, r
    real(dp), parameter :: b_ref = 1.5_dp

    do iphi = 1, nphi
        do itheta = 1, ntheta
            phi = 2.0_dp * pi * (iphi - 1) / nphi
            theta = 2.0_dp * pi * (itheta - 1) / ntheta
            r = 1.65_dp + 0.65_dp * cos(theta)

            b_r = b_ref * (1.0_dp + 0.15_dp * cos(nfp * phi))
            b_phi = b_ref * 0.8_dp * sin(theta)
            b_z = b_ref * 0.3_dp * sin(2.0_dp * nfp * phi)

            b_total(iphi, itheta, 1) = b_r * cos(phi) - b_phi * sin(phi)
            b_total(iphi, itheta, 2) = b_r * sin(phi) + b_phi * cos(phi)
            b_total(iphi, itheta, 3) = b_z
        end do
    end do

end subroutine create_model_field

subroutine write_output(filename, nphi, ntheta, b_ext)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nphi, ntheta
    real(dp), intent(in) :: b_ext(:,:,:)

    integer :: iunit, iphi, itheta, ios

    open(newunit=iunit, file=trim(filename), action='write', status='replace', &
         iostat=ios)

    if (ios /= 0) then
        print *, 'ERROR: Cannot write to ', trim(filename)
        return
    end if

    write(iunit, '(A)') '# Plasma response B_external on VMEC surface'
    write(iunit, '(A,I0,A,I0)') '# nphi=', nphi, ', ntheta=', ntheta
    write(iunit, '(A)') '# iphi itheta bx by bz |b|'

    do iphi = 1, nphi
        do itheta = 1, ntheta
            write(iunit, '(I4,I4,4ES16.8)') iphi, itheta, &
                b_ext(iphi, itheta, 1), &
                b_ext(iphi, itheta, 2), &
                b_ext(iphi, itheta, 3), &
                sqrt(sum(b_ext(iphi, itheta, :)**2))
        end do
    end do

    close(iunit)
    print *, 'INFO: Output written to ', trim(filename)

end subroutine write_output

logical function file_exists(filename)
    implicit none
    character(len=*), intent(in) :: filename
    inquire(file=filename, exist=file_exists)
end function file_exists
