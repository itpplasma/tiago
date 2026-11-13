module tiago_plasma_field_source
    !! Diagnostic-facing fields in Cartesian SI units. Exterior methods return
    !! the plasma contribution; the interior method returns the TOTAL field.
    !! A virtual-casing source must supply the linked toroidal flux separately
    !! for DIAGNO idia=1 loops; its exterior A does not include that flux.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: plasma_field_source_t, field_hint_t

    type :: field_hint_t
        !! Caller-owned inverse seed; volume zero requests a cold inverse.
        real(dp) :: u(3) = 0.0_dp
        integer :: volume = 0
    end type field_hint_t

    type, abstract :: plasma_field_source_t
        real(dp) :: diamagnetic_flux = 0.0_dp
        logical :: linked_flux_known = .false.
        logical :: enabled = .false.
    contains
        procedure(sample_b_interface), deferred :: sample_bfield
        procedure(sample_a_interface), deferred :: sample_vector_potential
        procedure :: sample_interior_bfield => unsupported_interior
        procedure :: sample_interior_bfield_status => unsupported_interior_status
        procedure :: sample_interior_bfield_with_hint => interior_with_hint
        procedure :: has_data
    end type plasma_field_source_t

    abstract interface
        subroutine sample_b_interface(self, points, bfield)
            import :: plasma_field_source_t, dp
            class(plasma_field_source_t), intent(in) :: self
            real(dp), intent(in) :: points(:, :)
            real(dp), intent(out) :: bfield(:, :)
        end subroutine sample_b_interface
        subroutine sample_a_interface(self, points, avec)
            import :: plasma_field_source_t, dp
            class(plasma_field_source_t), intent(in) :: self
            real(dp), intent(in) :: points(:, :)
            real(dp), intent(out) :: avec(:, :)
        end subroutine sample_a_interface
    end interface

contains

    subroutine interior_with_hint(self, point, value, status, hint)
        class(plasma_field_source_t), intent(in) :: self
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: value(3)
        integer, intent(out) :: status
        type(field_hint_t), intent(inout) :: hint
        real(dp) :: points(3, 1), values(3, 1)
        integer :: statuses(1)

        points(:, 1) = point
        call self%sample_interior_bfield_status(points, values, statuses)
        value = values(:, 1)
        status = statuses(1)
        if (status /= 0) hint = field_hint_t()
    end subroutine interior_with_hint

    logical function has_data(self)
        class(plasma_field_source_t), intent(in) :: self
        has_data = self%enabled
    end function has_data

    subroutine unsupported_interior_status(self, points, values, status)
        !! Nonfatal sampling lets tracers stop a line at the represented domain.
        class(plasma_field_source_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: status(:)

        values = 0.0_dp
        status = 1
    end subroutine unsupported_interior_status

    subroutine unsupported_interior(self, points, values)
        class(plasma_field_source_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)

        values = 0.0_dp
        if (size(points, 1) /= 3) error stop 'field points must have three components'
        error stop 'plasma source does not provide a total interior field'
    end subroutine unsupported_interior
end module tiago_plasma_field_source
