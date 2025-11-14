module tiago_plasma_support
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    implicit none

    logical, parameter :: tiago_plasma_available = .false.

    type :: plasma_support_t
        logical :: enabled = .false.
    contains
        procedure :: init_from_vmec => plasma_init_stub
        procedure :: finalize => plasma_finalize_stub
        procedure :: sample_bfield => plasma_sample_stub
        procedure :: sample_vector_potential => plasma_vector_potential_stub
        procedure :: sample_surface_bfield => plasma_sample_surface_stub
        procedure :: has_data => plasma_has_data_stub
    end type plasma_support_t

contains

    subroutine plasma_init_stub(self, wout_file, nphi, ntheta)
        class(plasma_support_t), intent(inout) :: self
        character(len=*), intent(in) :: wout_file
        integer(i32), intent(in) :: nphi
        integer(i32), intent(in) :: ntheta
        call plasma_disabled_error()
    end subroutine plasma_init_stub

    subroutine plasma_finalize_stub(self)
        class(plasma_support_t), intent(inout) :: self
        self%enabled = .false.
    end subroutine plasma_finalize_stub

    subroutine plasma_sample_stub(self, points, bfield)
        class(plasma_support_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: bfield(:, :)
        bfield = 0.0_dp
        call plasma_disabled_error()
    end subroutine plasma_sample_stub

    subroutine plasma_vector_potential_stub(self, points, avec)
        class(plasma_support_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: avec(:, :)
        avec = 0.0_dp
        call plasma_disabled_error()
    end subroutine plasma_vector_potential_stub

    subroutine plasma_sample_surface_stub(self, points, bfield)
        class(plasma_support_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: bfield(:, :)
        bfield = 0.0_dp
        call plasma_disabled_error()
    end subroutine plasma_sample_surface_stub

    logical function plasma_has_data_stub(self)
        class(plasma_support_t), intent(in) :: self
        plasma_has_data_stub = .false.
    end function plasma_has_data_stub

    subroutine plasma_disabled_error()
        use, intrinsic :: iso_fortran_env, only: error_unit
        write(error_unit, '(A)') 'TIAGO built without plasma response; reconfigure with TIAGO_ENABLE_PLASMA=ON'
        error stop 1
    end subroutine plasma_disabled_error

end module tiago_plasma_support
