program gen_magnetic_pickup
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
        operator(*), sin, cos
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: w, gain, delay, hr, hi, roots(6)
    type(kernel_spec_t) :: spec
    integer :: unit, last
    character(:), allocatable :: code
    call arena%init()
    w = sym(arena, 'omega')
    gain = sym(arena, 'gain')
    delay = sym(arena, 'delay')
    hr = sym(arena, 'hr')
    hi = sym(arena, 'hi')
    roots(1) = w*gain*(hi*cos(w*delay) - hr*sin(w*delay))
    roots(2) = -w*gain*(hr*cos(w*delay) + hi*sin(w*delay))
    roots(3) = diff(roots(1), gain)
    roots(4) = diff(roots(2), gain)
    roots(5) = diff(roots(1), delay)
    roots(6) = diff(roots(2), delay)
    spec%args = [str('omega'), str('gain'), str('delay'), str('hr'), str('hi')]
    spec%outputs = [str('tr'), str('ti'), str('gr'), str('gi'), str('dr'), str('di')]
    spec%name = str('magnetic_pickup_leaf')
    spec%module_name = str('tiago_magnetic_pickup_leaf')
    spec%generator = str('gen/gen_magnetic_pickup.f90')
    spec%generator_revision = str('fortsym combined cancellation native 20261028')
    spec%regenerate_command = str('build-codegen/gen_magnetic_pickup')
    spec%pure_procedure = .true.
    open (newunit=unit, file='src/generated/magnetic_pickup_leaf.f90', status='replace')
    code = chars(emit_kernel(roots, spec))
    last = len_trim(code)
    do while (last > 0)
        if (code(last:last) /= new_line('a')) exit
        last = last - 1
    end do
    write (unit, '(a)') code(:last)
    close (unit)
end program
