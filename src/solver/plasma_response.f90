module tiago_plasma_response
    !! Plasma response contribution to magnetic field from VMEC equilibrium
    !! using virtual casing principle via the virtual_casing C library
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, &
                                            c_bool, c_long, c_null_ptr, c_associated, c_loc
    implicit none

    private
    public :: plasma_response_t

    ! C interface to virtual_casing library (double precision)
    interface
        ! Create context
        function VirtualCasingCreateContextD() bind(C, name='VirtualCasingCreateContextD')
            use iso_c_binding
            type(c_ptr) :: VirtualCasingCreateContextD
        end function

        ! Destroy context - takes pointer-to-pointer (void**)
        subroutine VirtualCasingDestroyContextD(ctx_ptr) &
                bind(C, name='VirtualCasingDestroyContextD')
            use iso_c_binding
            type(c_ptr), value :: ctx_ptr
        end subroutine

        ! Setup from surface geometry and B-field
        subroutine VirtualCasingSetupD(digits, nfp, half_period, &
                                       nt, np, x, &
                                       src_nt, src_np, trg_nt, trg_np, ctx) &
                 bind(C, name='VirtualCasingSetupD')
            use iso_c_binding
            integer(c_int), value :: digits, nfp
            logical(c_bool), value :: half_period
            integer(c_long), value :: nt, np
            real(c_double), intent(in) :: x(*)
            integer(c_long), value :: src_nt, src_np, trg_nt, trg_np
            type(c_ptr), value :: ctx
        end subroutine

        ! Compute B_external from total B field
        subroutine VirtualCasingComputeBextD(bext, b, nt, np, ctx) &
                 bind(C, name='VirtualCasingComputeBextD')
            use iso_c_binding
            real(c_double), intent(out) :: bext(*)
            real(c_double), intent(in) :: b(*)
            integer(c_long), value :: nt, np
            type(c_ptr), value :: ctx
        end subroutine
    end interface

    ! Fortran wrapper type for virtual casing plasma response
    type :: plasma_response_t
        private
        type(c_ptr) :: ctx = c_null_ptr
        logical :: initialized = .false.
        integer :: nfp = 0
        integer :: src_nphi = 0, src_ntheta = 0
        integer :: trg_nphi = 0, trg_ntheta = 0
    contains
        procedure :: init => plasma_response_init
        procedure :: compute_bext => plasma_response_compute_bext
        procedure :: finalize => plasma_response_finalize
        procedure :: is_initialized => plasma_response_is_initialized
    end type plasma_response_t

contains

    ! Initialize virtual casing from VMEC surface and B-field
    subroutine plasma_response_init(self, nfp, x_surf, b_total, &
                                     src_nphi, src_ntheta, &
                                     trg_nphi, trg_ntheta, &
                                     use_stellsym, digits)
        class(plasma_response_t), intent(inout) :: self
        integer, intent(in) :: nfp
        real(dp), intent(in) :: x_surf(:,:,:)  ! (nphi, ntheta, 3)
        real(dp), intent(in) :: b_total(:,:,:)  ! (src_nphi, src_ntheta, 3)
        integer, intent(in) :: src_nphi, src_ntheta
        integer, intent(in), optional :: trg_nphi, trg_ntheta, digits
        logical, intent(in), optional :: use_stellsym

        integer :: nphi, ntheta, trg_nphi_val, trg_ntheta_val, digits_val
        logical(c_bool) :: stellsym_c
        real(c_double), allocatable :: x_flat(:), b_flat(:)
        integer :: i, j, k, idx

        ! Surface grid dimensions
        nphi = size(x_surf, 1)
        ntheta = size(x_surf, 2)

        ! Set defaults
        trg_nphi_val = nphi
        if (present(trg_nphi)) trg_nphi_val = trg_nphi

        trg_ntheta_val = ntheta
        if (present(trg_ntheta)) trg_ntheta_val = trg_ntheta

        digits_val = 6
        if (present(digits)) digits_val = digits

        stellsym_c = .true.
        if (present(use_stellsym)) stellsym_c = use_stellsym

        ! Create context
        self%ctx = VirtualCasingCreateContextD()
        if (.not. c_associated(self%ctx)) then
            return
        end if

        ! Flatten surface coordinates to C order: {x11, x12, ..., y11, y12, ..., z11, ...}
        allocate(x_flat(nphi * ntheta * 3))
        do k = 1, 3
            do j = 1, ntheta
                do i = 1, nphi
                    idx = (k-1) * nphi * ntheta + (j-1) * nphi + i
                    x_flat(idx) = real(x_surf(i, j, k), c_double)
                end do
            end do
        end do

        ! Flatten total B-field similarly
        allocate(b_flat(src_nphi * src_ntheta * 3))
        do k = 1, 3
            do j = 1, src_ntheta
                do i = 1, src_nphi
                    idx = (k-1) * src_nphi * src_ntheta + (j-1) * src_nphi + i
                    b_flat(idx) = real(b_total(i, j, k), c_double)
                end do
            end do
        end do

        ! Call setup
        call VirtualCasingSetupD(int(digits_val, c_int), int(nfp, c_int), stellsym_c, &
                                int(nphi, c_long), int(ntheta, c_long), x_flat, &
                                int(src_nphi, c_long), int(src_ntheta, c_long), &
                                int(trg_nphi_val, c_long), int(trg_ntheta_val, c_long), &
                                self%ctx)

        self%initialized = .true.
        self%nfp = nfp
        self%src_nphi = src_nphi
        self%src_ntheta = src_ntheta
        self%trg_nphi = trg_nphi_val
        self%trg_ntheta = trg_ntheta_val

        deallocate(x_flat, b_flat)
    end subroutine plasma_response_init

    ! Compute B_external (plasma contribution) from total B-field
    subroutine plasma_response_compute_bext(self, b_total, b_external)
        class(plasma_response_t), intent(in) :: self
        real(dp), intent(in) :: b_total(:,:,:)  ! (src_nphi, src_ntheta, 3)
        real(dp), intent(out) :: b_external(:,:,:)  ! (src_nphi, src_ntheta, 3)

        real(c_double), allocatable :: b_total_flat(:), b_ext_flat(:)
        integer :: i, j, k, idx, nphi, ntheta

        if (.not. self%initialized) then
            b_external = 0.0_dp
            return
        end if

        nphi = self%src_nphi
        ntheta = self%src_ntheta

        ! Flatten input B-field
        allocate(b_total_flat(nphi * ntheta * 3))
        allocate(b_ext_flat(nphi * ntheta * 3))

        do k = 1, 3
            do j = 1, ntheta
                do i = 1, nphi
                    idx = (k-1) * nphi * ntheta + (j-1) * nphi + i
                    b_total_flat(idx) = real(b_total(i, j, k), c_double)
                end do
            end do
        end do

        ! Call virtual casing computation
        call VirtualCasingComputeBextD(b_ext_flat, b_total_flat, &
                                      int(nphi, c_long), int(ntheta, c_long), &
                                      self%ctx)

        ! Unflatten result
        do k = 1, 3
            do j = 1, ntheta
                do i = 1, nphi
                    idx = (k-1) * nphi * ntheta + (j-1) * nphi + i
                    b_external(i, j, k) = real(b_ext_flat(idx), dp)
                end do
            end do
        end do

        deallocate(b_total_flat, b_ext_flat)
    end subroutine plasma_response_compute_bext

    ! Finalize and free resources
    subroutine plasma_response_finalize(self)
        class(plasma_response_t), intent(inout) :: self
        type(c_ptr), target :: ctx_tmp
        if (c_associated(self%ctx)) then
            ctx_tmp = self%ctx
            call VirtualCasingDestroyContextD(c_loc(ctx_tmp))
            self%ctx = c_null_ptr
            self%initialized = .false.
        end if
    end subroutine plasma_response_finalize

    ! Check if initialized
    pure logical function plasma_response_is_initialized(self)
        class(plasma_response_t), intent(in) :: self
        plasma_response_is_initialized = self%initialized
    end function plasma_response_is_initialized

end module tiago_plasma_response
