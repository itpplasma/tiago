program gen_kinetic
    !! Kinetic pressure and its partial derivatives (src/generated/kinetic_pressure.f90):
    !!   p = e (ne Te + n_i Ti) + p_fast,
    !!   n_i = n_D + n_Z, n_D = ne (Z - Zeff) / (Z - 1), n_Z = ne (Zeff - 1) / (Z (Z - 1)),
    !! from quasi-neutrality ne = n_D + Z n_Z and Zeff ne = n_D + Z^2 n_Z.
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
    type(expr_t) :: ne, te, ti, zeff, pf, z, e, nd, nz, roots(6)
    type(kernel_spec_t) :: spec
    character(len=:), allocatable :: code
    integer :: i, j, unit

    call arena%init()
    engine = make_symengine_engine(arena)
    ne = sym(arena, 'ne')
    te = sym(arena, 'te')
    ti = sym(arena, 'ti')
    zeff = sym(arena, 'zeff')
    pf = sym(arena, 'pfast')
    z = sym(arena, 'zimp')
    e = sym(arena, 'charge')
    nd = ne * (z - zeff) / (z - 1)
    nz = ne * (zeff - 1) / (z * (z - 1))
    roots(1) = e * (ne * te + (nd + nz) * ti) + pf
    roots(2) = diff(roots(1), ne)
    roots(3) = diff(roots(1), te)
    roots(4) = diff(roots(1), ti)
    roots(5) = diff(roots(1), zeff)
    roots(6) = diff(roots(1), pf)
    do i = 1, 6
        simplified = engine%simplify(roots(i))
        if (.not. simplified%ok) error stop 'kinetic pressure simplification failed'
        roots(i) = simplified%value
    end do
    allocate(spec%args(7), spec%outputs(6))
    spec%args = [str('ne'), str('te'), str('ti'), str('zeff'), str('pfast'), str('zimp'), str('charge')]
    spec%outputs = [str('p'), str('dp_dne'), str('dp_dte'), str('dp_dti'), str('dp_dzeff'), &
        str('dp_dpfast')]
    spec%name = str('kinetic_pressure')
    spec%module_name = str('tiago_kinetic_pressure')
    spec%generator = str('gen/gen_kinetic.f90')
    spec%generator_revision = str('fortsym 48dda45a34c0da0db1c78cd4130b3a97375e3c8e')
    spec%regenerate_command = str('build-codegen/gen_kinetic (see gen/README.md)')
    spec%pure_procedure = .true.
    code = chars(emit_kernel(roots, spec))
    j = len(code)
    do while (j > 0)
        if (code(j:j) /= new_line('a')) exit
        j = j - 1
    end do
    open(newunit=unit, file='src/generated/kinetic_pressure.f90', status='replace')
    write(unit, '(A)') code(:j)
    close(unit)
end program gen_kinetic
