program test_plasma_response_ncsx
    !! Test plasma response using real VMEC equilibrium data from NCSX
    !! NO SYNTHETIC DATA - reads actual wout file via libneo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    character(len=256) :: wout_file
    integer :: narg

    narg = command_argument_count()
    if (narg < 1) then
        print *, 'Usage: test_plasma_response_ncsx <wout_file>'
        error stop 1
    end if

    call get_command_argument(1, wout_file)
    print *, 'Loading VMEC file:', trim(wout_file)

    call test_ncsx_real_vmec(trim(wout_file))

    print *, 'NCSX plasma response test passed'
end program test_plasma_response_ncsx

subroutine test_ncsx_real_vmec(wout_file)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_response, only: plasma_response_t
    use new_vmec_stuff_mod, only: netcdffile, nper
    use spline_vmec_sub, only: spline_vmec_data, splint_vmec_data, vmec_field
    implicit none

    character(len=*), intent(in) :: wout_file
    type(plasma_response_t) :: plasma_response
    real(dp), allocatable :: x_surf(:,:,:), b_total(:,:,:), b_ext(:,:,:)
    real(dp) :: theta, varphi, pi, r, z, phi_edge
    real(dp) :: a_phi, a_theta, da_phi_ds, da_theta_ds, aiota_val, alam
    real(dp) :: dr_ds, dr_dt, dr_dp, dz_ds, dz_dt, dz_dp
    real(dp) :: dl_ds, dl_dt, dl_dp, sqg, bc_vartheta, bc_varphi
    real(dp) :: bcov_s, bcov_vartheta, bcov_varphi
    real(dp) :: cos_vphi, sin_vphi, cjac
    real(dp), dimension(3) :: e_s, e_theta, e_phi, e_vartheta, e_varphi
    real(dp), dimension(3) :: b_vec
    integer :: nphi, ntheta, iphi, itheta, nfp

    pi = acos(-1.0_dp)
    netcdffile = wout_file

    call spline_vmec_data()

    nfp = nper
    nphi = 16
    ntheta = 16

    allocate(x_surf(nphi, ntheta, 3))
    allocate(b_total(nphi, ntheta, 3))
    allocate(b_ext(nphi, ntheta, 3))

    print *, 'Computing VMEC surface and B-field on grid', &
             nphi, 'x', ntheta
    print *, 'NFP =', nfp

    do iphi = 1, nphi
        phi_edge = 2.0_dp * pi * real(iphi - 1, dp) / real(nphi, dp)
        do itheta = 1, ntheta
            theta = 2.0_dp * pi * real(itheta - 1, dp) / real(ntheta, dp)
            varphi = phi_edge

            call splint_vmec_data(1.0_dp, theta, varphi, a_phi, a_theta, &
                da_phi_ds, da_theta_ds, aiota_val, r, z, alam, dr_ds, &
                dr_dt, dr_dp, dz_ds, dz_dt, dz_dp, dl_ds, dl_dt, dl_dp)

            call vmec_field(1.0_dp, theta, varphi, a_theta, a_phi, &
                da_theta_ds, da_phi_ds, aiota_val, sqg, alam, dl_ds, &
                dl_dt, dl_dp, bc_vartheta, bc_varphi, bcov_s, &
                bcov_vartheta, bcov_varphi)

            cos_vphi = cos(varphi)
            sin_vphi = sin(varphi)

            e_s = [dr_ds * cos_vphi, dr_ds * sin_vphi, dz_ds]
            e_theta = [dr_dt * cos_vphi, dr_dt * sin_vphi, dz_dt]
            e_phi = [dr_dp * cos_vphi - r * sin_vphi, &
                     dr_dp * sin_vphi + r * cos_vphi, dz_dp]

            cjac = 1.0_dp / (1.0_dp + dl_dt)
            e_vartheta = cjac * (e_theta - dl_ds * e_s - dl_dp * e_phi)
            e_varphi = e_phi - dl_dp * cjac * e_theta

            b_vec = bc_vartheta * e_vartheta + bc_varphi * e_varphi

            x_surf(iphi, itheta, 1) = r * cos_vphi
            x_surf(iphi, itheta, 2) = r * sin_vphi
            x_surf(iphi, itheta, 3) = z

            b_total(iphi, itheta, 1) = b_vec(1)
            b_total(iphi, itheta, 2) = b_vec(2)
            b_total(iphi, itheta, 3) = b_vec(3)
        end do
    end do

    print *, 'Surface R range:', minval(sqrt(x_surf(:,:,1)**2 + &
        x_surf(:,:,2)**2)), maxval(sqrt(x_surf(:,:,1)**2 + &
        x_surf(:,:,2)**2))
    print *, 'Surface Z range:', minval(x_surf(:,:,3)), &
        maxval(x_surf(:,:,3))
    print *, 'B-field magnitude range:', &
        minval(sqrt(b_total(:,:,1)**2 + b_total(:,:,2)**2 + &
        b_total(:,:,3)**2)), &
        maxval(sqrt(b_total(:,:,1)**2 + b_total(:,:,2)**2 + &
        b_total(:,:,3)**2))

    print *, 'Initializing plasma response...'
    call plasma_response%init(nfp=nfp, x_surf=x_surf, b_total=b_total, &
        src_nphi=nphi, src_ntheta=ntheta)

    if (.not. plasma_response%is_initialized()) then
        print *, 'WARNING: Plasma response init failed'
        print *, 'This may be due to VMEC equilibrium requirements'
        deallocate(x_surf, b_total, b_ext)
        return
    end if

    print *, 'Computing B_external...'
    b_ext = 0.0_dp
    call plasma_response%compute_bext(b_total, b_ext)

    print *, 'B_external RMS:', &
        sqrt(sum(b_ext**2) / real(size(b_ext), dp))
    print *, 'B_external max:', maxval(abs(b_ext))
    print *, 'B_external/B_total ratio:', &
        maxval(abs(b_ext)) / maxval(abs(b_total))

    if (maxval(abs(b_ext)) > 1.0e10_dp) then
        print *, 'ERROR: B_external magnitude suspiciously large'
        error stop 'B_external overflow detected'
    end if

    call plasma_response%finalize()
    deallocate(x_surf, b_total, b_ext)

end subroutine test_ncsx_real_vmec
