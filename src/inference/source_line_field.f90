module tiago_source_line_field
    !! Field-line tracing field from the interior total field of an
    !! equilibrium source (for example SPECTRE), so that field-line labels
    !! can be computed on fields that are not on an mgrid.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_field_line_label, only: line_field_t
    use tiago_plasma_field_source, only: plasma_field_source_t
    implicit none
    private

    type, extends(line_field_t), public :: source_line_field_t
        class(plasma_field_source_t), pointer :: source => null()
    contains
        procedure :: eval => source_eval
    end type source_line_field_t

contains

    subroutine source_eval(self, r, phi, z, b, inside)
        class(source_line_field_t), intent(in) :: self
        real(dp), intent(in) :: r, phi, z
        real(dp), intent(out) :: b(3)
        logical, intent(out) :: inside
        real(dp) :: point(3, 1), cart(3, 1), c, s
        integer :: status(1)

        c = cos(phi)
        s = sin(phi)
        point(:, 1) = [r * c, r * s, z]
        call self%source%sample_interior_bfield_status(point, cart, status)
        inside = status(1) == 0
        b = 0.0_dp
        if (.not. inside) return
        b = [c * cart(1, 1) + s * cart(2, 1), -s * cart(1, 1) + c * cart(2, 1), cart(3, 1)]
    end subroutine source_eval

end module tiago_source_line_field
