program gen_robust_signals
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, log, sqrt, operator(+), &
                            operator(-), operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: residual, nu, u, value
    call arena%init()
    residual = sym(arena, 'residual'); nu = sym(arena, 'nu')
    u = residual**2/nu
    value = sqrt((nu + num(arena, 1))*log(num(arena, 1) + u))
    call emit('student_residual', value)
    value = residual*sqrt((nu + num(arena, 1))/nu)* &
            (num(arena, 1) - u/num(arena, 4) + num(arena, 13)*u**2/num(arena, 96))
    call emit('student_residual_small', value)
contains
    subroutine emit(name, value)
        character(*), intent(in) :: name
        type(expr_t), intent(in) :: value
        type(expr_t) :: roots(2)
        type(kernel_spec_t) :: spec
        character(:), allocatable :: code
        integer :: unit
        roots = [value, diff(value, residual)]
        spec%name = str(name); spec%module_name = str('tiago_'//name)
        spec%args = [str('residual'), str('nu')]
        spec%outputs = [str('value'), str('d_residual')]
        spec%generator = str('gen/gen_robust_signals.f90')
        spec%generator_revision = str('fortsym 21d157b394e775a999a8d9ed0be409eac1eb9574')
        spec%regenerate_command = str('build-codegen/gen_robust_signals')
        spec%pure_procedure = .true.
        code = chars(emit_kernel(roots, spec))
        open (newunit=unit, file='src/generated/'//name//'.f90', status='replace')
        write (unit, '(A)') code(:len(code) - 1)
        close (unit)
    end subroutine
end program
