program gen_nonthermal
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, tanh, cosh, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: s, p0, zc, zk, clo, chi, klo, khi, c, k, value, roots(5)
    type(kernel_spec_t) :: spec
    character(len=:), allocatable :: code
    integer :: unit
    call arena%init()
    s = sym(arena, 's'); p0 = sym(arena, 'p0')
    zc = sym(arena, 'zc'); zk = sym(arena, 'zk')
    clo = sym(arena, 'clo'); chi = sym(arena, 'chi')
    klo = sym(arena, 'klo'); khi = sym(arena, 'khi')
    ! Equivalent logistic map without exponentials that overflow in prior tails.
    c = clo + (chi - clo)*(1 + tanh(zc/2))/2
    k = klo + (khi - klo)*(1 + tanh(zk/2))/2
    value = c*p0*(1 - s)**k
    roots = [value, diff(value, p0), diff(value, zc), diff(value, zk), diff(value, s)]
    ! sech(t)**2 = 1-tanh(t)**2 avoids overflow in the differentiated map.
    roots(3) = subs(roots(3), cosh(zc/2)**(-2), 1-tanh(zc/2)**2)
    roots(4) = subs(roots(4), cosh(zk/2)**(-2), 1-tanh(zk/2)**2)
    spec%name = str('nonthermal_pressure')
    spec%module_name = str('tiago_nonthermal_pressure')
    spec%args = [str('s'), str('p0'), str('zc'), str('zk'), &
        str('clo'), str('chi'), str('klo'), str('khi')]
    spec%outputs = [str('p'), str('dp0'), str('dzc'), str('dzk'), str('ds')]
    spec%generator = str('gen/gen_nonthermal.f90')
    spec%generator_revision = str('fortsym 21d157b394e775a999a8d9ed0be409eac1eb9574')
    spec%regenerate_command = str('build-codegen/gen_nonthermal')
    spec%pure_procedure = .true.
    code = chars(emit_kernel(roots, spec))
    open (newunit=unit, file='src/generated/nonthermal_pressure.f90', status='replace')
    write (unit, '(A)') code(:len(code) - 1)
    close (unit)
end program gen_nonthermal
