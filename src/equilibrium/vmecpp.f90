module tiago_vmecpp
    !! Fortran binding of the VMEC++ adapter (src/equilibrium/vmecpp_adapter.h):
    !! fixed-boundary equilibria with the m = 1 gauge pinned, their geometry and
    !! radial profiles, wout fields, and the implicit adjoint of the solve.
    use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_int, c_double, c_char, &
        c_null_char, c_f_pointer, c_associated
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    interface
        function c_error() bind(C, name='tiago_vmecpp_error') result(message)
            import c_ptr
            type(c_ptr) :: message
        end function c_error
        function c_create(path, handle) bind(C, name='tiago_vmecpp_create') result(status)
            import c_char, c_ptr, c_int
            character(kind=c_char), intent(in) :: path(*)
            type(c_ptr), intent(out) :: handle
            integer(c_int) :: status
        end function c_create
        subroutine c_destroy(handle) bind(C, name='tiago_vmecpp_destroy')
            import c_ptr
            type(c_ptr), value :: handle
        end subroutine c_destroy
        function c_get_input(handle, name, i, j, value) bind(C, name='tiago_vmecpp_get_input') &
                result(status)
            import c_ptr, c_char, c_int, c_double
            type(c_ptr), value :: handle
            character(kind=c_char), intent(in) :: name(*)
            integer(c_int), value :: i, j
            real(c_double), intent(out) :: value
            integer(c_int) :: status
        end function c_get_input
        function c_set_input(handle, name, i, j, value) bind(C, name='tiago_vmecpp_set_input') &
                result(status)
            import c_ptr, c_char, c_int, c_double
            type(c_ptr), value :: handle
            character(kind=c_char), intent(in) :: name(*)
            integer(c_int), value :: i, j
            real(c_double), value :: value
            integer(c_int) :: status
        end function c_set_input
        function c_input_length(handle, name, length) bind(C, name='tiago_vmecpp_input_length') &
                result(status)
            import c_ptr, c_char, c_int
            type(c_ptr), value :: handle
            character(kind=c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: length
            integer(c_int) :: status
        end function c_input_length
        function c_get_int(handle, name, value) bind(C, name='tiago_vmecpp_get_int') result(status)
            import c_ptr, c_char, c_int
            type(c_ptr), value :: handle
            character(kind=c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: value
            integer(c_int) :: status
        end function c_get_int
        function c_get_string(handle, name, buffer, length) bind(C, name='tiago_vmecpp_get_string') &
                result(status)
            import c_ptr, c_char, c_int
            type(c_ptr), value :: handle
            character(kind=c_char), intent(in) :: name(*)
            character(kind=c_char), intent(out) :: buffer(*)
            integer(c_int), value :: length
            integer(c_int) :: status
        end function c_get_string
        function c_solve(handle, converged) bind(C, name='tiago_vmecpp_solve') result(status)
            import c_ptr, c_int
            type(c_ptr), value :: handle
            integer(c_int), intent(out) :: converged
            integer(c_int) :: status
        end function c_solve
        function c_geometry(handle, coefficients, toroidal_flux, poloidal_flux) &
                bind(C, name='tiago_vmecpp_geometry') result(status)
            import c_ptr, c_double, c_int
            type(c_ptr), value :: handle
            real(c_double), intent(out) :: coefficients(*), toroidal_flux(*), poloidal_flux(*)
            integer(c_int) :: status
        end function c_geometry
        function c_radial(handle, phip_full, phip_half, iota_half, current_half, mass_half, &
                lamscale) bind(C, name='tiago_vmecpp_radial') result(status)
            import c_ptr, c_double, c_int
            type(c_ptr), value :: handle
            real(c_double), intent(out) :: phip_full(*), phip_half(*), iota_half(*)
            real(c_double), intent(out) :: current_half(*), mass_half(*), lamscale
            integer(c_int) :: status
        end function c_radial
        function c_wout(handle, name, values, length) bind(C, name='tiago_vmecpp_wout') &
                result(status)
            import c_ptr, c_char, c_double, c_int
            type(c_ptr), value :: handle
            character(kind=c_char), intent(in) :: name(*)
            real(c_double), intent(out) :: values(*)
            integer(c_int), value :: length
            integer(c_int) :: status
        end function c_wout
        function c_adjoint(handle, ncot, geometry_bar, boundary_bar, profile_bar, max_dense) &
                bind(C, name='tiago_vmecpp_adjoint') result(status)
            import c_ptr, c_int, c_double
            type(c_ptr), value :: handle
            integer(c_int), value :: ncot, max_dense
            real(c_double), intent(in) :: geometry_bar(*)
            real(c_double), intent(out) :: boundary_bar(*), profile_bar(*)
            integer(c_int) :: status
        end function c_adjoint
    end interface

    type, public :: vmecpp_t
        type(c_ptr) :: handle = c_null_ptr
    contains
        procedure :: create => vmecpp_create
        procedure :: destroy => vmecpp_destroy
        procedure :: get_input => vmecpp_get_input
        procedure :: set_input => vmecpp_set_input
        procedure :: input_length => vmecpp_input_length
        procedure :: get_int => vmecpp_get_int
        procedure :: get_string => vmecpp_get_string
        procedure :: solve => vmecpp_solve
        procedure :: geometry => vmecpp_geometry
        procedure :: radial => vmecpp_radial
        procedure :: wout => vmecpp_wout
        procedure :: adjoint => vmecpp_adjoint
    end type vmecpp_t

    public :: vmecpp_error

contains

    function vmecpp_error() result(message)
        character(len=:), allocatable :: message
        character(kind=c_char), pointer :: chars(:)
        integer :: n

        call c_f_pointer(c_error(), chars, [4096])
        n = 0
        do while (n < 4096)
            if (chars(n + 1) == c_null_char) exit
            n = n + 1
        end do
        allocate(character(len=n) :: message)
        message = transfer(chars(:n), message)
    end function vmecpp_error

    subroutine check(status, what)
        integer(c_int), intent(in) :: status
        character(len=*), intent(in) :: what
        if (status /= 0) error stop 'VMEC++ ('//what//'): '//vmecpp_error()
    end subroutine check

    function cstr(text) result(c)
        character(len=*), intent(in) :: text
        character(kind=c_char, len=len_trim(text) + 1) :: c
        c = trim(text)//c_null_char
    end function cstr

    subroutine vmecpp_create(self, json_path)
        class(vmecpp_t), intent(inout) :: self
        character(len=*), intent(in) :: json_path
        call self%destroy()
        call check(c_create(cstr(json_path), self%handle), 'reading '//json_path)
    end subroutine vmecpp_create

    subroutine vmecpp_destroy(self)
        class(vmecpp_t), intent(inout) :: self
        if (c_associated(self%handle)) call c_destroy(self%handle)
        self%handle = c_null_ptr
    end subroutine vmecpp_destroy

    real(dp) function vmecpp_get_input(self, name, i, j) result(value)
        class(vmecpp_t), intent(in) :: self
        character(len=*), intent(in) :: name
        integer, intent(in), optional :: i, j
        real(c_double) :: v
        call check(c_get_input(self%handle, cstr(name), opt(i), opt(j), v), 'get '//name)
        value = v
    end function vmecpp_get_input

    subroutine vmecpp_set_input(self, name, value, i, j)
        class(vmecpp_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: value
        integer, intent(in), optional :: i, j
        call check(c_set_input(self%handle, cstr(name), opt(i), opt(j), value), 'set '//name)
    end subroutine vmecpp_set_input

    integer function vmecpp_input_length(self, name) result(length)
        class(vmecpp_t), intent(in) :: self
        character(len=*), intent(in) :: name
        integer(c_int) :: n
        call check(c_input_length(self%handle, cstr(name), n), 'length '//name)
        length = n
    end function vmecpp_input_length

    integer function vmecpp_get_int(self, name) result(value)
        class(vmecpp_t), intent(in) :: self
        character(len=*), intent(in) :: name
        integer(c_int) :: v
        call check(c_get_int(self%handle, cstr(name), v), 'get '//name)
        value = v
    end function vmecpp_get_int

    function vmecpp_get_string(self, name) result(value)
        class(vmecpp_t), intent(in) :: self
        character(len=*), intent(in) :: name
        character(len=:), allocatable :: value
        character(kind=c_char) :: buffer(64)
        call check(c_get_string(self%handle, cstr(name), buffer, 64_c_int), 'get '//name)
        allocate(character(len=64) :: value)
        value = transfer(buffer, value)
        value = trim(value)
    end function vmecpp_get_string

    logical function vmecpp_solve(self, message) result(converged)
        !! .false. (and the reason in message) if VMEC++ found no equilibrium.
        class(vmecpp_t), intent(inout) :: self
        character(len=:), allocatable, intent(out), optional :: message
        integer(c_int) :: ok
        call check(c_solve(self%handle, ok), 'solve')
        converged = ok == 1
        if (present(message)) message = vmecpp_error()
    end function vmecpp_solve

    subroutine vmecpp_geometry(self, coefficients, toroidal_flux, poloidal_flux)
        !! coefficients(0:ntor, 0:mpol-1, ns, 12): VMEC++ product-basis blocks.
        class(vmecpp_t), intent(in) :: self
        real(dp), allocatable, intent(out) :: coefficients(:, :, :, :)
        real(dp), allocatable, intent(out) :: toroidal_flux(:), poloidal_flux(:)
        integer :: ns, mpol, ntor
        ns = self%get_int('ns')
        mpol = self%get_int('mpol')
        ntor = self%get_int('ntor')
        allocate(coefficients(0:ntor, 0:mpol - 1, ns, 12), toroidal_flux(ns), poloidal_flux(ns))
        call check(c_geometry(self%handle, coefficients, toroidal_flux, poloidal_flux), 'geometry')
    end subroutine vmecpp_geometry

    subroutine vmecpp_radial(self, phip_full, phip_half, iota_half, current_half, mass_half, &
            lamscale)
        class(vmecpp_t), intent(in) :: self
        real(dp), allocatable, intent(out) :: phip_full(:), phip_half(:), iota_half(:)
        real(dp), allocatable, intent(out) :: current_half(:), mass_half(:)
        real(dp), intent(out) :: lamscale
        integer :: ns
        ns = self%get_int('ns')
        allocate(phip_full(ns), phip_half(ns - 1), iota_half(ns - 1), current_half(ns - 1), &
            mass_half(ns - 1))
        call check(c_radial(self%handle, phip_full, phip_half, iota_half, current_half, &
            mass_half, lamscale), 'radial profiles')
    end subroutine vmecpp_radial

    function vmecpp_wout(self, name, length) result(values)
        class(vmecpp_t), intent(in) :: self
        character(len=*), intent(in) :: name
        integer, intent(in) :: length
        real(dp) :: values(length)
        call check(c_wout(self%handle, cstr(name), values, int(length, c_int)), 'wout '//name)
    end function vmecpp_wout

    subroutine vmecpp_adjoint(self, geometry_bar, boundary_bar, profile_bar, max_dense)
        !! geometry_bar(size of the 12 geometry blocks, ncot) ->
        !! boundary_bar(2 * mpol * (2 ntor + 1), ncot), profile_bar(3 * (ns - 1), ncot).
        class(vmecpp_t), intent(inout) :: self
        real(dp), intent(in) :: geometry_bar(:, :)
        real(dp), allocatable, intent(out) :: boundary_bar(:, :), profile_bar(:, :)
        integer, intent(in) :: max_dense
        integer :: ns, mpol, ntor
        ns = self%get_int('ns')
        mpol = self%get_int('mpol')
        ntor = self%get_int('ntor')
        allocate(boundary_bar(2 * mpol * (2 * ntor + 1), size(geometry_bar, 2)))
        allocate(profile_bar(3 * (ns - 1), size(geometry_bar, 2)))
        call check(c_adjoint(self%handle, int(size(geometry_bar, 2), c_int), geometry_bar, &
            boundary_bar, profile_bar, int(max_dense, c_int)), 'adjoint')
    end subroutine vmecpp_adjoint

    integer(c_int) function opt(i)
        integer, intent(in), optional :: i
        opt = 1
        if (present(i)) opt = i
    end function opt
end module tiago_vmecpp
