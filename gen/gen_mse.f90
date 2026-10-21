program gen_mse
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
        operator(*), operator(/), operator(**), atan2
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none
    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(engine_result_t) :: simplified
    type(expr_t) :: x(9), roots(6), num, den
    type(kernel_spec_t) :: spec
    character(len=2), parameter :: names(9) = &
        ['vx', 'vy', 'vz', 'hx', 'hy', 'hz', 'zx', 'zy', 'zz']
    character(len=2), parameter :: field_names(9) = &
        ['bx', 'by', 'bz', 'ax', 'ay', 'az', 'cx', 'cy', 'cz']
    integer :: i, unit

    call arena%init()
    engine = make_symengine_engine(arena)
    do i = 1, 9
        x(i) = sym(arena, names(i))
    end do
    ! h . (v x B) = (h x v) . B; likewise for zero axis z.
    roots(1) = x(5)*x(3) - x(6)*x(2)
    roots(2) = x(6)*x(1) - x(4)*x(3)
    roots(3) = x(4)*x(2) - x(5)*x(1)
    roots(4) = x(8)*x(3) - x(9)*x(2)
    roots(5) = x(9)*x(1) - x(7)*x(3)
    roots(6) = x(7)*x(2) - x(8)*x(1)
    allocate (spec%args(9), spec%outputs(6))
    do i = 1, 9
        spec%args(i) = str(names(i))
    end do
    spec%outputs = [str('ax'), str('ay'), str('az'), &
        str('cx'), str('cy'), str('cz')]
    call emit('mse_rows', roots)

    do i = 1, 9
        x(i) = sym(arena, field_names(i))
        spec%args(i) = str(field_names(i))
    end do
    num = x(4)*x(1) + x(5)*x(2) + x(6)*x(3)
    den = x(7)*x(1) + x(8)*x(2) + x(9)*x(3)
    roots(1) = atan2(num, den)
    do i = 1, 3
        roots(i + 1) = diff(roots(1), x(i))
    end do
    deallocate (spec%outputs)
    spec%outputs = [str('angle'), str('gx'), str('gy'), str('gz')]
    call emit('mse_local', roots(:4))
contains
    subroutine emit(name, expressions)
        character(len=*), intent(in) :: name
        type(expr_t), intent(inout) :: expressions(:)
        integer :: j
        character(len=:), allocatable :: code
        do j = 1, size(expressions)
            simplified = engine%simplify(expressions(j))
            if (.not. simplified%ok) error stop 'MSE simplification failed'
            expressions(j) = simplified%value
        end do
        spec%name = str(name)
        spec%module_name = str('tiago_'//name)
        spec%generator = str('gen/gen_mse.f90')
        spec%generator_revision = &
            str('fortsym 7e8ef7abdd90c49bb5f75e2c119e7d181656bcc6')
        spec%regenerate_command = str('build-codegen/gen_mse (see gen/README.md)')
        spec%pure_procedure = .true.
        open (newunit=unit, file='src/generated/'//name//'.f90', status='replace')
        code = chars(emit_kernel(expressions, spec))
        j = len(code)
        do while (j > 0)
            if (code(j:j) /= new_line('a')) exit
            j = j - 1
        end do
        write (unit, '(A)') code(:j)
        close (unit)
    end subroutine emit
end program gen_mse
