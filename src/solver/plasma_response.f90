module tiago_plasma_response
    !! Plasma response contribution computed via surface Biot-Savart integrals
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_surface_biot_savart, only: surface_biot_savart_t
    implicit none

    private
    public :: plasma_response_t

    type :: plasma_response_t
        type(surface_biot_savart_t) :: evaluator
        real(dp), allocatable :: surface_nodes(:, :, :)
        integer :: src_nphi = 0
        integer :: src_ntheta = 0
        integer :: nfp = 0
        logical :: initialized = .false.
    contains
        procedure :: init => plasma_response_init
        procedure :: compute_bext => plasma_response_compute_bext
        procedure :: compute_bext_at => plasma_response_compute_bext_at
        procedure :: compute_bext_points => plasma_response_compute_bext_points
        procedure :: compute_surface_bext_points => plasma_response_compute_surface
        procedure :: compute_vector_potential_points => plasma_response_compute_vector_potential
        procedure :: finalize => plasma_response_finalize
        procedure :: is_initialized => plasma_response_is_initialized
    end type plasma_response_t

contains

    subroutine plasma_response_init(self, nfp, x_surf, b_total, src_nphi, src_ntheta)
        class(plasma_response_t), intent(inout) :: self
        integer, intent(in) :: nfp
        real(dp), intent(in) :: x_surf(:, :, :)
        real(dp), intent(in) :: b_total(:, :, :)
        integer, intent(in) :: src_nphi
        integer, intent(in) :: src_ntheta

        call self%finalize()

        if (src_nphi < 2 .or. src_ntheta < 2) then
            error stop 'plasma_response requires at least a 2x2 surface grid'
        end if
        if (size(x_surf, 1) /= src_nphi .or. size(x_surf, 2) /= src_ntheta) then
            error stop 'plasma_response: surface grid mismatch'
        end if
        if (size(b_total, 1) /= src_nphi .or. size(b_total, 2) /= src_ntheta) then
            error stop 'plasma_response: B_total grid mismatch'
        end if

        allocate(self%surface_nodes(src_nphi, src_ntheta, 3))
        self%surface_nodes = x_surf
        call self%evaluator%init(x_surf, b_total, nfp)

        self%src_nphi = src_nphi
        self%src_ntheta = src_ntheta
        self%nfp = nfp
        self%initialized = .true.
    end subroutine plasma_response_init

    subroutine plasma_response_compute_bext(self, b_total, b_external)
        class(plasma_response_t), intent(in) :: self
        real(dp), intent(in) :: b_total(:, :, :)
        real(dp), intent(out) :: b_external(:, :, :)

        if (.not. self%initialized) then
            b_external = 0.0_dp
            return
        end if
        if (size(b_total, 1) /= self%src_nphi .or. size(b_total, 2) /= self%src_ntheta) then
            error stop 'compute_bext: input grid mismatch'
        end if
        if (size(b_external, 1) /= self%src_nphi .or. size(b_external, 2) /= self%src_ntheta) then
            error stop 'compute_bext: output shape mismatch'
        end if

        call sample_on_surface(self, self%surface_nodes, b_external)
    end subroutine plasma_response_compute_bext

    subroutine plasma_response_compute_bext_at(self, b_total, x_eval, b_external)
        class(plasma_response_t), intent(in) :: self
        real(dp), intent(in) :: b_total(:, :, :)
        real(dp), intent(in) :: x_eval(:, :, :)
        real(dp), intent(out) :: b_external(:, :, :)

        if (size(b_total, 1) /= self%src_nphi .or. size(b_total, 2) /= self%src_ntheta) then
            error stop 'compute_bext_at: input grid mismatch'
        end if
        call sample_grid(self, x_eval, b_external)
    end subroutine plasma_response_compute_bext_at

    subroutine plasma_response_compute_bext_points(self, b_total, points, b_external)
        class(plasma_response_t), intent(in) :: self
        real(dp), intent(in) :: b_total(:, :, :)
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: b_external(:, :)

        if (size(b_total, 1) /= self%src_nphi .or. size(b_total, 2) /= self%src_ntheta) then
            error stop 'compute_bext_points: input grid mismatch'
        end if
        call sample_points(self, points, b_external)
    end subroutine plasma_response_compute_bext_points

    subroutine plasma_response_compute_surface(self, points, b_external)
        class(plasma_response_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: b_external(:, :)
        call sample_points(self, points, b_external)
    end subroutine plasma_response_compute_surface

    subroutine plasma_response_compute_vector_potential(self, points, avec)
        class(plasma_response_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: avec(:, :)

        if (.not. self%initialized) then
            avec = 0.0_dp
            return
        end if
        if (size(points, 2) /= 3) error stop 'vector potential expects 3D points'
        if (size(avec, 1) /= size(points, 1) .or. size(avec, 2) /= 3) then
            error stop 'vector potential output mismatch'
        end if

        call self%evaluator%sample_vector_potential(points, avec)
    end subroutine plasma_response_compute_vector_potential

    subroutine plasma_response_finalize(self)
        class(plasma_response_t), intent(inout) :: self
        if (allocated(self%surface_nodes)) deallocate(self%surface_nodes)
        call self%evaluator%finalize()
        self%initialized = .false.
        self%src_nphi = 0
        self%src_ntheta = 0
        self%nfp = 0
    end subroutine plasma_response_finalize

    pure logical function plasma_response_is_initialized(self)
        class(plasma_response_t), intent(in) :: self
        plasma_response_is_initialized = self%initialized
    end function plasma_response_is_initialized

    subroutine sample_on_surface(self, points, field)
        class(plasma_response_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :, :)
        real(dp), intent(out) :: field(:, :, :)

        real(dp), allocatable :: flat_points(:, :)
        real(dp), allocatable :: flat_field(:, :)
        integer :: npts

        npts = size(points, 1) * size(points, 2)
        allocate(flat_points(npts, 3))
        allocate(flat_field(npts, 3))
        call reshape_points(points, flat_points)
        call self%evaluator%sample_b_field(flat_points, flat_field)
        call reshape_field(flat_field, field)
        deallocate(flat_points, flat_field)
    end subroutine sample_on_surface

    subroutine sample_grid(self, x_eval, field)
        class(plasma_response_t), intent(in) :: self
        real(dp), intent(in) :: x_eval(:, :, :)
        real(dp), intent(out) :: field(:, :, :)

        real(dp), allocatable :: flat_points(:, :)
        real(dp), allocatable :: flat_field(:, :)
        integer :: npts

        if (.not. self%initialized) then
            field = 0.0_dp
            return
        end if

        npts = size(x_eval, 1) * size(x_eval, 2)
        allocate(flat_points(npts, 3))
        allocate(flat_field(npts, 3))
        call reshape_points(x_eval, flat_points)
        call self%evaluator%sample_b_field(flat_points, flat_field)
        call reshape_field(flat_field, field)
        deallocate(flat_points, flat_field)
    end subroutine sample_grid

    subroutine sample_points(self, points, field)
        class(plasma_response_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: field(:, :)

        if (.not. self%initialized) then
            field = 0.0_dp
            return
        end if
        if (size(points, 2) /= 3) error stop 'sample_points expects Cartesian inputs'
        if (size(field, 1) /= size(points, 1) .or. size(field, 2) /= 3) then
            error stop 'sample_points output mismatch'
        end if

        call self%evaluator%sample_b_field(points, field)
    end subroutine sample_points

    subroutine reshape_points(grid_points, flat_points)
        real(dp), intent(in) :: grid_points(:, :, :)
        real(dp), intent(out) :: flat_points(:, :)
        integer :: iphi, itheta, idx

        idx = 0
        do iphi = 1, size(grid_points, 1)
            do itheta = 1, size(grid_points, 2)
                idx = idx + 1
                flat_points(idx, :) = grid_points(iphi, itheta, :)
            end do
        end do
    end subroutine reshape_points

    subroutine reshape_field(flat_field, grid_field)
        real(dp), intent(in) :: flat_field(:, :)
        real(dp), intent(out) :: grid_field(:, :, :)
        integer :: iphi, itheta, idx

        idx = 0
        do iphi = 1, size(grid_field, 1)
            do itheta = 1, size(grid_field, 2)
                idx = idx + 1
                grid_field(iphi, itheta, :) = flat_field(idx, :)
            end do
        end do
    end subroutine reshape_field

end module tiago_plasma_response
