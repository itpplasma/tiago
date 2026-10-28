program gen_wall
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, real_expr, pi_expr, operator(+), &
        operator(-), operator(*), operator(/), operator(**), exp, sqrt
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: p(6), f, b, y
    character(16) :: names(6)
    call arena%init()
    call args(['energy0', 'p      '])
    f = p(1)*p(2)**2
    call emit('wall_legacy_energy', f, 2, [2])
    call args(['energy_ev', 'weight   ', 'fraction '])
    f = real_expr(arena, 1.602176634e-19_dp)*p(1)*p(2)*p(3)
    call emit('wall_energy', f, 3, [1])
    call args(['energy', 'area  ', 'dt    '])
    f = p(1)/(p(2)*p(3))
    call emit('wall_flux', f, 3, [1])
    call args(['flux ', 'kappa', 'cap  ', 'lag0 ', 'lag1 '])
    f = 2*p(1)*(sqrt(p(4)) - sqrt(p(5)))/sqrt(pi_expr(arena)*p(2)*p(3))
    call emit('wall_thermal', f, 5, [1])
    call args(['temp  ', 'lambda', 'eps   ', 'bg    '])
    y = exp(-real_expr(arena, 6.62607015e-34_dp*299792458.0_dp &
        /1.380649e-23_dp)/(p(2)*p(1)))
    b = real_expr(arena, 2*6.62607015e-34_dp*299792458.0_dp**2) &
        *y/(p(2)**5*(1 - y))
    f = p(3)*b + (1 - p(3))*p(4)
    call emit('wall_radiance', f, 4, [1])
contains
    subroutine args(input)
        character(*), intent(in) :: input(:)
        integer :: j
        do j = 1, size(input)
            names(j) = input(j)
            p(j) = sym(arena, trim(input(j)))
        end do
    end subroutine
    subroutine emit(name, expression, n, active)
        character(*), intent(in) :: name
        type(expr_t), intent(in) :: expression
        integer, intent(in) :: n, active(:)
        type(kernel_spec_t) :: spec
        type(expr_t) :: roots(1 + size(active))
        character(:), allocatable :: code
        character(16) :: deriv
        integer :: j, unit, last
        allocate (spec%args(n), spec%outputs(size(roots)))
        do j = 1, n
            spec%args(j) = str(trim(names(j)))
        end do
        roots(1) = expression
        spec%outputs(1) = str('value')
        do j = 1, size(active)
            roots(j + 1) = diff(expression, p(active(j)))
            write (deriv, '(A,I0)') 'd', active(j)
            spec%outputs(j + 1) = str(trim(deriv))
        end do
        spec%name = str(name)
        spec%module_name = str('tiago_'//name)
        spec%generator = str('gen/gen_wall.f90')
        spec%generator_revision = str('fortsym 21d157b394e775a999a8d9ed0be409eac1eb9574; CRC64 corrected')
        spec%regenerate_command = str('build-codegen/gen_wall (see docs/wall-synthetic.md)')
        spec%pure_procedure = .true.
        code = chars(emit_kernel(roots, spec))
        last = len(code)
        do while (last > 0)
            if (code(last:last) /= new_line('a')) exit
            last = last - 1
        end do
        open (newunit=unit, file='src/generated/'//name//'.f90', status='replace')
        write (unit, '(A)') code(:last)
        close (unit)
    end subroutine
end program
