program debug_tiago_virtual_casing
    !! Debug version that saves arrays for comparison with simsopt
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_response, only: plasma_response_t
    use new_vmec_stuff_mod, only: netcdffile, nper
    use spline_vmec_sub, only: spline_vmec_data, splint_vmec_data
    use vmec_field_tools, only: vmec_field_cylindrical
    implicit none

    character(len=256) :: wout_file
    type(plasma_response_t) :: pr
    real(dp), allocatable :: x_surf(:,:,:), b_total(:,:,:), b_ext(:,:,:)
    real(dp) :: theta, varphi, pi, r, z
    real(dp) :: a_phi, a_theta, da_phi_ds, da_theta_ds, aiota_val, alam
    real(dp) :: dr_ds, dr_dt, dr_dp, dz_ds, dz_dt, dz_dp, dl_ds, dl_dt, dl_dp
    real(dp) :: BR, Bphi, BZ, Bmag
    real(dp) :: cos_vphi, sin_vphi
    integer :: nphi, ntheta, iphi, itheta, nfp, iunit
    real(dp), parameter :: cm_to_m = 1.0e-2_dp
    real(dp), parameter :: gauss_to_tesla = 1.0e-4_dp

    call get_command_argument(1, wout_file)
    pi = acos(-1.0_dp)
    netcdffile = wout_file

    call spline_vmec_data()

    nfp = nper
    nphi = 16
    ntheta = 16

    allocate(x_surf(nphi, ntheta, 3))
    allocate(b_total(nphi, ntheta, 3))
    allocate(b_ext(nphi, ntheta, 3))

    print *, 'TIAGO DEBUG [VERSION 2.0]: Computing VMEC surface and B-field'
    print *, 'Grid:', nphi, 'x', ntheta, ', NFP=', nfp
    write(*, *) '*** ABOUT TO ENTER LOOP ***'
    call flush(6)

    do iphi = 1, nphi
        varphi = 2.0_dp * pi * real(iphi - 1, dp) / real(nphi, dp)
        do itheta = 1, ntheta
            theta = 2.0_dp * pi * real(itheta - 1, dp) / real(ntheta, dp)

            ! Get position (R, Z) in cm
            call splint_vmec_data(1.0_dp, theta, varphi, a_phi, a_theta, &
                da_phi_ds, da_theta_ds, aiota_val, r, z, alam, dr_ds, &
                dr_dt, dr_dp, dz_ds, dz_dt, dz_dp, dl_ds, dl_dt, dl_dp)

            ! Get B-field in cylindrical coordinates (CGS: Gauss)
            call vmec_field_cylindrical(1.0_dp, theta, varphi, BR, Bphi, BZ, Bmag)

            cos_vphi = cos(varphi)
            sin_vphi = sin(varphi)

            if (iphi == 1 .and. itheta == 1) then
                print *, '=== USING vmec_field_cylindrical() AT [iphi=1, itheta=1] ==='
                print *, 'Input: theta=', theta, ' varphi=', varphi
                print *, 'Position: r=', r, ' z=', z, ' (cm from libneo)'
                print *, 'B-field (CGS cylindrical): BR=', BR, ' Bphi=', Bphi, ' BZ=', BZ
                print *, '|B| (Gauss) =', Bmag
                print *, 'Converting to Cartesian and Tesla...'
                print *, '=== END DEBUG ==='
            end if

            ! Position: convert from cylindrical (R, phi, Z) to Cartesian (x, y, z) in meters
            x_surf(iphi, itheta, 1) = r * cos_vphi * cm_to_m
            x_surf(iphi, itheta, 2) = r * sin_vphi * cm_to_m
            x_surf(iphi, itheta, 3) = z * cm_to_m

            ! B-field: convert from cylindrical (BR, Bphi, BZ) to Cartesian (Bx, By, Bz)
            ! and from Gauss to Tesla
            b_total(iphi, itheta, 1) = (BR * cos_vphi - Bphi * sin_vphi) * gauss_to_tesla
            b_total(iphi, itheta, 2) = (BR * sin_vphi + Bphi * cos_vphi) * gauss_to_tesla
            b_total(iphi, itheta, 3) = BZ * gauss_to_tesla
        end do
    end do

    print *, 'x_surf[1,1,:] =', x_surf(1,1,:)
    print *, 'x_surf range:', minval(x_surf), maxval(x_surf)
    print *, 'b_total[1,1,:] =', b_total(1,1,:)
    print *, 'b_total range:', minval(b_total), maxval(b_total)

    open(newunit=iunit, file='tiago_x_surf.dat', status='replace')
    do iphi = 1, nphi
        do itheta = 1, ntheta
            write(iunit, '(2I5,3ES20.12)') iphi, itheta, x_surf(iphi,itheta,:)
        end do
    end do
    close(iunit)

    open(newunit=iunit, file='tiago_b_total.dat', status='replace')
    do iphi = 1, nphi
        do itheta = 1, ntheta
            write(iunit, '(2I5,3ES20.12)') iphi, itheta, b_total(iphi,itheta,:)
        end do
    end do
    close(iunit)

    print *, 'Calling plasma_response%init...'
    call pr%init(nfp=nfp, x_surf=x_surf, b_total=b_total, &
        src_nphi=nphi, src_ntheta=ntheta, use_stellsym=.true., digits=6)

    if (.not. pr%is_initialized()) then
        print *, 'ERROR: init failed'
        stop 1
    end if

    b_ext = 0.0_dp
    call pr%compute_bext(b_total, b_ext)

    print *, 'b_ext[1,1,:] =', b_ext(1,1,:)
    print *, 'b_ext range:', minval(b_ext), maxval(b_ext)
    print *, 'b_ext RMS:', sqrt(sum(b_ext**2) / real(size(b_ext), dp))

    open(newunit=iunit, file='tiago_b_ext.dat', status='replace')
    do iphi = 1, nphi
        do itheta = 1, ntheta
            write(iunit, '(2I5,3ES20.12)') iphi, itheta, b_ext(iphi,itheta,:)
        end do
    end do
    close(iunit)

    print *, 'Saved:'
    print *, '  tiago_x_surf.dat'
    print *, '  tiago_b_total.dat'
    print *, '  tiago_b_ext.dat'

    call pr%finalize()
    deallocate(x_surf, b_total, b_ext)

end program debug_tiago_virtual_casing
