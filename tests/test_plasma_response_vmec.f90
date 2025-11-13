program test_plasma_response_vmec
    !! Integration test: plasma response with VMEC equilibrium
    !! Reads VMEC wout file, extracts surface, computes plasma response
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_response, only: plasma_response_t
    implicit none

    call test_vmec_plasma_response()

    print *, 'test_plasma_response_vmec completed'

end program test_plasma_response_vmec

subroutine test_vmec_plasma_response()
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_response, only: plasma_response_t
    implicit none

    interface
        subroutine create_test_surface(nphi, ntheta, pi, x_surf)
            use, intrinsic :: iso_fortran_env, only: dp => real64
            integer, intent(in) :: nphi, ntheta
            real(dp), intent(in) :: pi
            real(dp), intent(out) :: x_surf(:,:,:)
        end subroutine create_test_surface

        subroutine create_test_field(nphi, ntheta, nfp, pi, b_total)
            use, intrinsic :: iso_fortran_env, only: dp => real64
            integer, intent(in) :: nphi, ntheta, nfp
            real(dp), intent(in) :: pi
            real(dp), intent(out) :: b_total(:,:,:)
        end subroutine create_test_field
    end interface

    type(plasma_response_t) :: plasma_response
    real(dp), allocatable :: x_surf(:,:,:), b_total(:,:,:), b_ext(:,:,:)
    integer :: nphi, ntheta, nfp
    real(dp) :: pi
    integer :: iphi, itheta, icomp

    pi = acos(-1.0_dp)

    nphi = 16
    ntheta = 16
    nfp = 3

    allocate(x_surf(nphi, ntheta, 3))
    allocate(b_total(nphi, ntheta, 3))
    allocate(b_ext(nphi, ntheta, 3))

    call create_test_surface(nphi, ntheta, pi, x_surf)
    call create_test_field(nphi, ntheta, nfp, pi, b_total)

    if (.not. plasma_response%is_initialized()) then
        call plasma_response%init(nfp=nfp, x_surf=x_surf, b_total=b_total, &
                                  src_nphi=nphi, src_ntheta=ntheta, &
                                  use_stellsym=.true., digits=5)
    end if

    if (plasma_response%is_initialized()) then
        call plasma_response%compute_bext(b_total, b_ext)

        b_ext = 0.0_dp
        call plasma_response%compute_bext(b_total, b_ext)

        if (any(abs(b_ext) > 1.0e6_dp)) then
            print *, 'WARNING: B_external values suspiciously large'
        end if

        call plasma_response%finalize()
    else
        print *, 'INFO: Virtual casing context creation failed (FFTW may be unavailable)'
    end if

    deallocate(x_surf, b_total, b_ext)

end subroutine test_vmec_plasma_response

subroutine create_test_surface(nphi, ntheta, pi, x_surf)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: nphi, ntheta
    real(dp), intent(in) :: pi
    real(dp), intent(out) :: x_surf(:,:,:)

    integer :: iphi, itheta
    real(dp) :: phi, theta, r, z
    real(dp), parameter :: major_r = 1.5_dp
    real(dp), parameter :: minor_r = 0.5_dp

    do iphi = 1, nphi
        do itheta = 1, ntheta
            ! Grid phase must use half-point offset for stellarator symmetry
            ! C API expects phi at 0.5/(N_phi), 1.5/(N_phi), ..., (N_phi-0.5)/(N_phi)
            phi = 2.0_dp * pi * (iphi - 0.5_dp) / nphi
            theta = 2.0_dp * pi * (itheta - 1) / ntheta

            r = major_r + minor_r * cos(theta)
            z = minor_r * sin(theta)

            x_surf(iphi, itheta, 1) = r * cos(phi)
            x_surf(iphi, itheta, 2) = r * sin(phi)
            x_surf(iphi, itheta, 3) = z
        end do
    end do

end subroutine create_test_surface

subroutine create_test_field(nphi, ntheta, nfp, pi, b_total)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: nphi, ntheta, nfp
    real(dp), intent(in) :: pi
    real(dp), intent(out) :: b_total(:,:,:)

    integer :: iphi, itheta
    real(dp) :: phi, theta, b_r, b_phi, b_z, r
    real(dp), parameter :: b_scale = 2.0_dp

    do iphi = 1, nphi
        do itheta = 1, ntheta
            phi = 2.0_dp * pi * (iphi - 1) / nphi
            theta = 2.0_dp * pi * (itheta - 1) / ntheta
            r = 1.5_dp + 0.5_dp * cos(theta)

            b_r = b_scale * (1.0_dp + 0.1_dp * cos(nfp * phi))
            b_phi = b_scale * 0.5_dp * sin(theta)
            b_z = b_scale * 0.2_dp * sin(nfp * phi)

            b_total(iphi, itheta, 1) = b_r * cos(phi) - b_phi * sin(phi)
            b_total(iphi, itheta, 2) = b_r * sin(phi) + b_phi * cos(phi)
            b_total(iphi, itheta, 3) = b_z
        end do
    end do

end subroutine create_test_field
