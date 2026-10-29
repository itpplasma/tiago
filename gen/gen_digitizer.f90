program gen_digitizer
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(-), operator(*), operator(/)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none
    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(expr_t) :: volts, vres, voff
    type(expr_t) :: roots(4)
    call arena%init()
    engine = make_symengine_engine(arena)
    volts = sym(arena, 'volts')
    vres = sym(arena, 'vres')
    voff = sym(arena, 'voff')
    ! Ideal linear scale of an ADC input: counts = (volts - offset) / resolution.
    roots(1) = (volts - voff)/vres
    roots(2) = diff(roots(1), volts)
    roots(3) = diff(roots(1), vres)
    roots(4) = diff(roots(1), voff)
    call emit('digitizer_signal', [str('volts'), str('vres'), str('voff')], &
        [str('counts'), str('d_volts'), str('d_vres'), str('d_voff')], roots)
contains
    subroutine emit(name, args, outputs, expressions)
        use fortsym_string, only: str_t
        character(*), intent(in) :: name
        type(str_t), intent(in) :: args(:), outputs(:)
        type(expr_t), intent(inout) :: expressions(:)
        type(kernel_spec_t) :: spec
        type(engine_result_t) :: simplified
        character(:), allocatable :: code
        integer :: i, j, unit
        do i = 1, size(expressions)
            simplified = engine%simplify(expressions(i))
            if (.not. simplified%ok) error stop 'digitizer simplification failed'
            expressions(i) = simplified%value
        end do
        spec%args = args
        spec%outputs = outputs
        spec%name = str(name)
        spec%module_name = str('tiago_'//name)
        spec%generator = str('gen/gen_digitizer.f90')
        spec%generator_revision = str('fortsym 21d157b394e775a999a8d9ed0be409eac1eb9574')
        spec%regenerate_command = str('build-codegen/gen_digitizer')
        spec%pure_procedure = .true.
        code = chars(emit_kernel(expressions, spec))
        j = len_trim(code)
        do while (j > 0)
            if (code(j:j) /= new_line('a')) exit
            j = j - 1
        end do
        open (newunit=unit, file='src/generated/'//name//'.f90', status='replace')
        write (unit, '(A)') code(:j)
        close (unit)
    end subroutine
end program
