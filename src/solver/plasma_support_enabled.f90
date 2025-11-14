module tiago_plasma_support
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    use tiago_plasma_response, only: plasma_response_t
    use new_vmec_stuff_mod, only: netcdffile, nper
    use spline_vmec_sub, only: spline_vmec_data, splint_vmec_data
    use vmec_field_tools, only: vmec_field_cylindrical
    implicit none

    logical, parameter :: tiago_plasma_available = .true.

    type :: plasma_support_t
        type(plasma_response_t) :: ctx
        real(dp), allocatable :: x_surf(:, :, :)
        real(dp), allocatable :: b_total(:, :, :)
        logical :: enabled = .false.
        integer(i32) :: nphi = 0_i32
        integer(i32) :: ntheta = 0_i32
        integer(i32) :: nfp = 0_i32
    contains
        procedure :: init_from_vmec => plasma_init_from_vmec
        procedure :: finalize => plasma_finalize
        procedure :: sample_bfield => plasma_sample_bfield
        procedure :: sample_vector_potential => plasma_sample_vector_potential
        procedure :: sample_surface_bfield => plasma_sample_surface_bfield
        procedure :: has_data => plasma_has_data
    end type plasma_support_t

contains

    subroutine plasma_init_from_vmec(self, wout_file, src_nphi, src_ntheta)
        class(plasma_support_t), intent(inout) :: self
        character(len=*), intent(in) :: wout_file
        integer(i32), intent(in) :: src_nphi
        integer(i32), intent(in) :: src_ntheta

        call self%finalize()
        call load_vmec_surface(trim(wout_file), src_nphi, src_ntheta, self%x_surf, self%b_total, self%nfp)
        call self%ctx%init(nfp=self%nfp, x_surf=self%x_surf, b_total=self%b_total, &
            src_nphi=src_nphi, src_ntheta=src_ntheta)
        self%enabled = self%ctx%is_initialized()
        if (self%enabled) then
            self%nphi = src_nphi
            self%ntheta = src_ntheta
        else
            call self%finalize()
            error stop 'plasma virtual casing initialization failed'
        end if
    end subroutine plasma_init_from_vmec

    subroutine plasma_finalize(self)
        class(plasma_support_t), intent(inout) :: self
        if (self%ctx%is_initialized()) call self%ctx%finalize()
        if (allocated(self%x_surf)) deallocate(self%x_surf)
        if (allocated(self%b_total)) deallocate(self%b_total)
        self%enabled = .false.
        self%nphi = 0_i32
        self%ntheta = 0_i32
        self%nfp = 0_i32
    end subroutine plasma_finalize

    subroutine plasma_sample_bfield(self, points, bfield)
        class(plasma_support_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: bfield(:, :)

        if (.not. self%enabled) then
            bfield = 0.0_dp
            return
        end if

        if (size(points, 2) /= 3) error stop 'sample_bfield expects Cartesian points'
        if (size(bfield, 1) /= size(points, 1) .or. size(bfield, 2) /= 3) then
            error stop 'sample_bfield: output array shape mismatch'
        end if

        ! Batch evaluation - all points at once!
        call self%ctx%compute_bext_points(self%b_total, points, bfield)
    end subroutine plasma_sample_bfield

    subroutine plasma_sample_vector_potential(self, points, avec)
        class(plasma_support_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: avec(:, :)

        ! Virtual-casing library does not provide batch vector potential evaluation.
        ! Use numerical differentiation of B-field if A-field is required off-surface.
        error stop 'plasma vector potential evaluation not supported by virtual-casing API'
    end subroutine plasma_sample_vector_potential

    subroutine plasma_sample_surface_bfield(self, points, bfield)
        class(plasma_support_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: bfield(:, :)

        ! For B-field on the plasma surface, use VMEC directly via vmec_field_tools module.
        ! Virtual-casing is only needed for B-field OUTSIDE the plasma boundary.
        error stop 'plasma on-surface B-field: use VMEC vmec_field_cylindrical() directly'
    end subroutine plasma_sample_surface_bfield

    logical function plasma_has_data(self)
        class(plasma_support_t), intent(in) :: self
        plasma_has_data = self%enabled
    end function plasma_has_data

    subroutine load_vmec_surface(wout_file, nphi, ntheta, x_surf, b_total, &
            nfp)
        character(len=*), intent(in) :: wout_file
        integer(i32), intent(in) :: nphi
        integer(i32), intent(in) :: ntheta
        real(dp), allocatable, intent(out) :: x_surf(:, :, :)
        real(dp), allocatable, intent(out) :: b_total(:, :, :)
        integer(i32), intent(out) :: nfp

        integer :: iphi
        integer :: itheta
        real(dp) :: phi_edge
        real(dp) :: theta
        real(dp) :: varphi
        real(dp) :: a_phi
        real(dp) :: a_theta
        real(dp) :: da_phi_ds
        real(dp) :: da_theta_ds
        real(dp) :: aiota_val
        real(dp) :: alam
        real(dp) :: r
        real(dp) :: z
        real(dp) :: dr_ds
        real(dp) :: dr_dt
        real(dp) :: dr_dp
        real(dp) :: dz_ds
        real(dp) :: dz_dt
        real(dp) :: dz_dp
        real(dp) :: dl_ds
        real(dp) :: dl_dt
        real(dp) :: dl_dp
        real(dp) :: br
        real(dp) :: bphi
        real(dp) :: bz
        real(dp) :: bmag
        real(dp) :: cos_vphi
        real(dp) :: sin_vphi
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp), parameter :: cm_to_m = 1.0e-2_dp
        real(dp), parameter :: gauss_to_tesla = 1.0e-4_dp

        netcdffile = trim(wout_file)
        call spline_vmec_data()
        nfp = nper

        if (allocated(x_surf)) then
            deallocate(x_surf)
        end if
        if (allocated(b_total)) then
            deallocate(b_total)
        end if
        allocate(x_surf(nphi, ntheta, 3))
        allocate(b_total(nphi, ntheta, 3))

        do iphi = 1, nphi
            phi_edge = 2.0_dp * pi * real(iphi - 1, dp) / real(nphi, dp)
            do itheta = 1, ntheta
                theta = 2.0_dp * pi * real(itheta - 1, dp) &
                    / real(ntheta, dp)
                varphi = phi_edge

                call splint_vmec_data(1.0_dp, theta, varphi, a_phi, a_theta, &
                    da_phi_ds, da_theta_ds, aiota_val, r, z, alam, dr_ds, &
                    dr_dt, dr_dp, dz_ds, dz_dt, dz_dp, dl_ds, dl_dt, dl_dp)

                call vmec_field_cylindrical(1.0_dp, theta, varphi, br, bphi, &
                    bz, bmag)

                cos_vphi = cos(varphi)
                sin_vphi = sin(varphi)

                x_surf(iphi, itheta, 1) = r * cos_vphi * cm_to_m
                x_surf(iphi, itheta, 2) = r * sin_vphi * cm_to_m
                x_surf(iphi, itheta, 3) = z * cm_to_m

                b_total(iphi, itheta, 1) = (br * cos_vphi - bphi * sin_vphi) &
                    * gauss_to_tesla
                b_total(iphi, itheta, 2) = (br * sin_vphi + bphi * cos_vphi) &
                    * gauss_to_tesla
                b_total(iphi, itheta, 3) = bz * gauss_to_tesla
            end do
        end do
    end subroutine load_vmec_surface

end module tiago_plasma_support
