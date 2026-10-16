module tiago_plasma_field_source
    !! Diagnostic-facing fields in Cartesian SI units. Exterior methods return
    !! the plasma contribution; the interior method returns the TOTAL field.
    !! A virtual-casing source must supply the linked toroidal flux separately
    !! for DIAGNO idia=1 loops; its exterior A does not include that flux.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: plasma_field_source_t

    type, abstract :: plasma_field_source_t
        real(dp) :: diamagnetic_flux = 0.0_dp
        logical :: linked_flux_known = .false.
        logical :: enabled = .false.
    contains
        procedure(sample_b_interface), deferred :: sample_bfield
        procedure(sample_a_interface), deferred :: sample_vector_potential
        procedure :: sample_interior_bfield => unsupported_interior
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

    logical function has_data(self)
        class(plasma_field_source_t), intent(in) :: self
        has_data = self%enabled
    end function has_data

    subroutine unsupported_interior(self, points, values)
        class(plasma_field_source_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)

        values = 0.0_dp
        if (size(points, 1) /= 3) error stop 'field points must have three components'
        error stop 'plasma source does not provide a total interior field'
    end subroutine unsupported_interior
end module tiago_plasma_field_source
