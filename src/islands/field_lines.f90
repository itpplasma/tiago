module tiago_field_lines
    !! Total magnetic field for field-line tracing, sampled in right-handed
    !! cylindrical coordinates x = (R, phi, Z) [m, rad, m] with components
    !! b = (B_R, B_phi, B_Z) [T]. A field is any sum of filament coils and one
    !! plasma source that provides its TOTAL interior field (SPECTRE, or a
    !! nested source plus coils), so vacuum, nested and island states are
    !! traced and measured by the same code. A nested state (VMEC) adds its
    !! interior field inside rho <= 1, so a nested equilibrium plus a changed
    !! coil field is traced the same way.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_coil_kernels, only: coil_set_t, build_coil_set
    use tiago_plasma_field_source, only: plasma_field_source_t
    use tiago_real_space_state, only: real_space_state_t, located_t, &
        vector_jets_t, locate_ok
    use tiago_cylindrical_components, only: cylindrical_components
    implicit none
    private
    public :: line_field_t, cartesian_line_field_t, read_scaled_coils

    type, abstract :: line_field_t
    contains
        procedure(sample_i), deferred :: sample
    end type line_field_t

    abstract interface
        subroutine sample_i(self, x, b, status)
            !! status 0: valid; nonzero: outside the field's domain.
            import :: line_field_t, dp
            class(line_field_t), intent(in) :: self
            real(dp), intent(in) :: x(3)
            real(dp), intent(out) :: b(3)
            integer, intent(out) :: status
        end subroutine sample_i
    end interface

    type, extends(line_field_t) :: cartesian_line_field_t
        !! Coil filaments plus an optional plasma source's total interior
        !! field (SPECTRE) and an optional nested state's field (VMEC). A
        !! point outside the nested state's plasma (rho > 1) is outside the
        !! domain.
        type(coil_set_t) :: coils
        logical :: has_coils = .false.
        class(plasma_field_source_t), allocatable :: plasma
        class(real_space_state_t), allocatable :: state
        real(dp) :: r_min = 0.0_dp
    contains
        procedure :: sample => cartesian_sample
    end type cartesian_line_field_t

contains

    subroutine cartesian_sample(self, x, b, status)
        class(cartesian_line_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3)
        integer, intent(out) :: status
        real(dp) :: point(3, 1), bc(3, 1), total(3), cyl(2), jac(2, 3)
        type(located_t) :: loc
        type(vector_jets_t) :: jet
        integer :: st

        b = 0.0_dp
        status = 1
        if (.not. all(ieee_is_finite(x))) return
        if (x(1) <= self%r_min) return
        point(:, 1) = [x(1)*cos(x(2)), x(1)*sin(x(2)), x(3)]
        total = 0.0_dp
        if (self%has_coils) then
            call self%coils%field(point, bc)
            total = total + bc(:, 1)
        end if
        if (allocated(self%plasma)) then
            call self%plasma%sample_interior_bfield(point, bc)
            total = total + bc(:, 1)
        end if
        if (allocated(self%state)) then
            call self%state%locate(point(:, 1), loc)
            if (loc%status /= locate_ok .or. loc%u(1) > 1.0_dp) return
            call self%state%vector_jets(loc%u, jet, st, parameters=.false.)
            if (st /= locate_ok) return
            total = total + jet%b
        end if
        call cylindrical_components(x(2), total(1), total(2), cyl, jac)
        b = [cyl(1), cyl(2), total(3)]
        if (all(ieee_is_finite(b))) status = 0
    end subroutine cartesian_sample

    subroutine read_scaled_coils(paths, scales, coils, extcur)
        !! STELLOPT/MAKEGRID coil files, each with its file currents multiplied
        !! by one scale factor, joined into one filament set. A nonblank
        !! extcur(k) names an EXTCUR file that sets the group currents of file
        !! k before scaling. A file's closing point carries zero current, so
        !! the joints between files and coils are current-free jumpers.
        use tiago_coil_loader, only: load_coils_into_field
        use neo_biotsavart_field, only: biotsavart_field_t
        character(len=*), intent(in) :: paths(:)
        real(dp), intent(in) :: scales(:)
        type(coil_set_t), intent(out) :: coils
        character(len=*), intent(in), optional :: extcur(:)
        type(biotsavart_field_t) :: one
        logical :: with_extcur
        real(dp), allocatable :: x(:), y(:), z(:), current(:)
        integer :: k

        if (size(paths) /= size(scales)) error stop 'one scale per coil file'
        allocate (x(0), y(0), z(0), current(0))
        do k = 1, size(paths)
            with_extcur = .false.
            if (present(extcur)) with_extcur = len_trim(extcur(k)) > 0
            if (with_extcur) then
                call load_coils_into_field(one, trim(paths(k)), trim(extcur(k)))
            else
                call load_coils_into_field(one, trim(paths(k)))
            end if
            x = [x, one%coils%x]
            y = [y, one%coils%y]
            z = [z, one%coils%z]
            current = [current, scales(k)*one%coils%current]
        end do
        coils = build_coil_set(x, y, z, current)
    end subroutine read_scaled_coils
end module tiago_field_lines
