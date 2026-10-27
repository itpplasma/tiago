program gen_recorded
    !! Recorded gain and offset, with partial derivatives.
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
                            operator(*), operator(/)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none
    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(engine_result_t) :: simplified
    type(expr_t) :: raw, gain, offset, roots(4)
    type(kernel_spec_t) :: spec
    character(len=:), allocatable :: code
    integer :: i, j, unit

    call arena%init()
    engine = make_symengine_engine(arena)
    raw = sym(arena, 'raw')
    gain = sym(arena, 'gain')
    offset = sym(arena, 'offset')
    roots(1) = gain*raw + offset
    roots(2) = diff(roots(1), raw)
    roots(3) = diff(roots(1), gain)
    roots(4) = diff(roots(1), offset)
    do i = 1, 4
        simplified = engine%simplify(roots(i))
        if (.not. simplified%ok) error stop 'recorded signal simplification failed'
        roots(i) = simplified%value
    end do
    allocate (spec%args(3), spec%outputs(4))
    spec%args = [str('raw'), str('gain'), str('offset')]
    spec%outputs = [str('signal'), str('d_raw'), str('d_gain'), str('d_offset')]
    spec%name = str('recorded_affine')
    spec%module_name = str('tiago_recorded_affine')
    spec%generator = str('gen/gen_recorded.f90')
    spec%generator_revision = str('fortsym 21d157b394e775a999a8d9ed0be409eac1eb9574')
    spec%regenerate_command = str('build-codegen/gen_recorded (see gen/README.md)')
    spec%pure_procedure = .true.
    code = chars(emit_kernel(roots, spec))
    j = len(code)
    do while (j > 0)
        if (code(j:j) /= new_line('a')) exit
        j = j - 1
    end do
    open (newunit=unit, file='src/generated/recorded_affine.f90', status='replace')
    write (unit, '(A)') code(:j)
    close (unit)
end program gen_recorded
