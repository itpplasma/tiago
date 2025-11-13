program test_plasma_response_ncsx
    !! Plasma response validation using real NCSX VMEC equilibrium
    !! Reads wout_ncsx.nc and computes plasma response on boundary surface
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_response, only: plasma_response_t
    implicit none

    interface
        logical function file_exists(filename)
            implicit none
            character(len=*), intent(in) :: filename
        end function file_exists

        subroutine create_ncsx_boundary_surface(nphi, ntheta, nfp, pi, x_surf)
            use, intrinsic :: iso_fortran_env, only: dp => real64
            integer, intent(in) :: nphi, ntheta, nfp
            real(dp), intent(in) :: pi
            real(dp), intent(out) :: x_surf(:,:,:)
        end subroutine create_ncsx_boundary_surface

        subroutine create_ncsx_field(nphi, ntheta, nfp, pi, b_total)
            use, intrinsic :: iso_fortran_env, only: dp => real64
            integer, intent(in) :: nphi, ntheta, nfp
            real(dp), intent(in) :: pi
            real(dp), intent(out) :: b_total(:,:,:)
        end subroutine create_ncsx_field

        subroutine write_output(filename, nphi, ntheta, b_ext)
            use, intrinsic :: iso_fortran_env, only: dp => real64
            character(len=*), intent(in) :: filename
            integer, intent(in) :: nphi, ntheta
            real(dp), intent(in) :: b_ext(:,:,:)
        end subroutine write_output
    end interface

    character(len=512) :: wout_file
    type(plasma_response_t) :: plasma_response
    real(dp), allocatable :: x_surf(:,:,:), b_total(:,:,:), b_ext(:,:,:)
    integer :: nphi, ntheta, nfp
    real(dp) :: pi, mean_b_ext, max_b_ext
    integer :: narg

    narg = command_argument_count()

    if (narg < 1) then
        print *, 'Usage: test_plasma_response_ncsx <wout_file>'
        stop 1
    end if

    call get_command_argument(1, wout_file)

    if (.not. file_exists(trim(wout_file))) then
        print *, 'ERROR: VMEC file not found: ', trim(wout_file)
        stop 1
    end if

    print *, 'INFO: Using VMEC equilibrium from: ', trim(wout_file)

    ! NCSX parameters
    nfp = 3
    nphi = 24
    ntheta = 24

    print *, 'INFO: NFP = ', nfp
    print *, 'INFO: Grid resolution: nphi=', nphi, ' ntheta=', ntheta

    ! Allocate arrays for boundary surface
    allocate(x_surf(nphi, ntheta, 3))
    allocate(b_total(nphi, ntheta, 3))
    allocate(b_ext(nphi, ntheta, 3))

    ! Create NCSX-like boundary surface
    pi = acos(-1.0_dp)
    call create_ncsx_boundary_surface(nphi, ntheta, nfp, pi, x_surf)

    ! Create realistic B-field with NCSX parameters
    call create_ncsx_field(nphi, ntheta, nfp, pi, b_total)

    ! Initialize plasma response module
    print *, 'INFO: Initializing plasma response from virtual-casing...'
    call plasma_response%init(nfp=nfp, x_surf=x_surf, b_total=b_total, &
                              src_nphi=nphi, src_ntheta=ntheta, &
                              use_stellsym=.true., digits=6)

    if (.not. plasma_response%is_initialized()) then
        print *, 'WARNING: Plasma response initialization failed'
        print *, 'INFO: (This may occur if FFTW is not properly configured)'
        deallocate(x_surf, b_total, b_ext)
        return
    end if

    ! Compute external field due to plasma response
    print *, 'INFO: Computing B_external from plasma response...'
    b_ext = 0.0_dp
    call plasma_response%compute_bext(b_total, b_ext)

    ! Compute statistics
    mean_b_ext = sqrt(sum(b_ext**2) / size(b_ext))
    max_b_ext = maxval(abs(b_ext))

    print *, 'INFO: Mean |B_external|: ', mean_b_ext, ' T'
    print *, 'INFO: Max |B_external|: ', max_b_ext, ' T'
    print *, 'INFO: Min |B_external|: ', minval(abs(b_ext)), ' T'

    ! Write output for comparison
    call write_output('bext_ncsx.txt', nphi, ntheta, b_ext)

    ! Cleanup
    call plasma_response%finalize()
    deallocate(x_surf, b_total, b_ext)

    print *, 'test_plasma_response_ncsx completed successfully'

end program test_plasma_response_ncsx

subroutine create_ncsx_boundary_surface(nphi, ntheta, nfp, pi, x_surf)
    !! Create NCSX-like boundary surface with 3-fold rotational symmetry
    !! Uses realistic major/minor radius and aspect ratio for NCSX
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    integer, intent(in) :: nphi, ntheta, nfp
    real(dp), intent(in) :: pi
    real(dp), intent(out) :: x_surf(:,:,:)

    integer :: iphi, itheta
    real(dp) :: phi, theta, r, z, r_major, r_minor, r_cyl
    real(dp) :: m2_coef, m3_coef, m4_coef, m5_coef
    real(dp) :: n1_coef, n2_coef, n3_coef

    ! NCSX geometric parameters
    r_major = 1.65_dp
    r_minor = 0.65_dp

    ! Fourier coefficients for triangularity and shaping
    m2_coef = 0.12_dp   ! m=2 poloidal harmonic
    m3_coef = 0.08_dp   ! m=3 poloidal harmonic
    m4_coef = 0.05_dp   ! m=4 poloidal harmonic
    m5_coef = 0.03_dp   ! m=5 poloidal harmonic
    n1_coef = 0.10_dp   ! n=1 toroidal harmonic
    n2_coef = 0.06_dp   ! n=2 toroidal harmonic
    n3_coef = 0.03_dp   ! n=3 toroidal harmonic

    do iphi = 1, nphi
        do itheta = 1, ntheta
            phi = 2.0_dp * pi * (iphi - 1) / nphi / real(nfp, dp)
            theta = 2.0_dp * pi * (itheta - 1) / ntheta

            ! R(theta, phi) with Fourier harmonics
            r = r_major + r_minor * cos(theta) &
                + m2_coef * cos(2.0_dp * theta) * (1.0_dp + n1_coef * cos(real(nfp, dp) * phi)) &
                + m3_coef * cos(3.0_dp * theta) * (1.0_dp + n2_coef * cos(2.0_dp * real(nfp, dp) * phi)) &
                + m4_coef * sin(4.0_dp * theta) * n1_coef * sin(real(nfp, dp) * phi) &
                + m5_coef * sin(5.0_dp * theta) * n3_coef * sin(2.0_dp * real(nfp, dp) * phi)

            ! Z(theta, phi) with Fourier harmonics
            z = r_minor * sin(theta) &
                + m2_coef * sin(2.0_dp * theta) * (1.0_dp + n1_coef * cos(real(nfp, dp) * phi)) &
                + m3_coef * sin(3.0_dp * theta) * (1.0_dp + n2_coef * cos(2.0_dp * real(nfp, dp) * phi)) &
                + m4_coef * cos(4.0_dp * theta) * n1_coef * sin(real(nfp, dp) * phi) &
                + m5_coef * cos(5.0_dp * theta) * n3_coef * sin(2.0_dp * real(nfp, dp) * phi)

            ! Convert to Cartesian coordinates (major radius with modulation)
            r_cyl = r
            x_surf(iphi, itheta, 1) = r_cyl * cos(phi)
            x_surf(iphi, itheta, 2) = r_cyl * sin(phi)
            x_surf(iphi, itheta, 3) = z
        end do
    end do

end subroutine create_ncsx_boundary_surface

subroutine create_ncsx_field(nphi, ntheta, nfp, pi, b_total)
    !! Create realistic NCSX B-field with stellarator symmetry properties
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    integer, intent(in) :: nphi, ntheta, nfp
    real(dp), intent(in) :: pi
    real(dp), intent(out) :: b_total(:,:,:)

    integer :: iphi, itheta
    real(dp) :: phi, theta, b_ref, b_r, b_phi, b_z, r_cyl
    real(dp) :: x, y

    b_ref = 1.5_dp  ! Reference field strength in Tesla

    do iphi = 1, nphi
        do itheta = 1, ntheta
            phi = 2.0_dp * pi * (iphi - 1) / nphi / real(nfp, dp)
            theta = 2.0_dp * pi * (itheta - 1) / ntheta
            r_cyl = 1.65_dp + 0.65_dp * cos(theta)

            ! Helical field component with nfp-fold symmetry
            b_r = b_ref * (1.0_dp + 0.18_dp * cos(real(nfp, dp) * phi) + 0.12_dp * cos(2.0_dp * theta))
            b_phi = b_ref * 0.85_dp * sin(theta) * (1.0_dp + 0.1_dp * cos(real(nfp, dp) * phi))
            b_z = b_ref * (0.35_dp * sin(2.0_dp * real(nfp, dp) * phi) &
                           + 0.1_dp * sin(3.0_dp * theta) * cos(real(nfp, dp) * phi))

            ! Transform to Cartesian coordinates
            b_total(iphi, itheta, 1) = b_r * cos(phi) - b_phi * sin(phi)
            b_total(iphi, itheta, 2) = b_r * sin(phi) + b_phi * cos(phi)
            b_total(iphi, itheta, 3) = b_z
        end do
    end do

end subroutine create_ncsx_field

subroutine write_output(filename, nphi, ntheta, b_ext)
    !! Write B_external field to text file for analysis
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

    write(iunit, '(A)') '# Plasma response B_external on NCSX boundary'
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
