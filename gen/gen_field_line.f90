program gen_field_line
    !! Field-line direction in cylindrical coordinates with the toroidal angle
    !! as the independent variable (src/generated/field_line_rhs.f90):
    !!   dR/dphi = R B_R / B_phi,   dZ/dphi = R B_Z / B_phi.
    use fortsym_string, only: str, chars
    use fortsym_expr, only: expr_t, sym, operator(*), operator(/)
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none
    type(symengine_engine_t) :: engine
    type(engine_result_t) :: simplified
    type(expr_t) :: r, br, bphi, bz, roots(2)
    type(kernel_spec_t) :: spec
    character(len=:), allocatable :: code
    integer :: i, j, unit

    engine = make_symengine_engine()
    r = sym('r')
    br = sym('br')
    bphi = sym('bphi')
    bz = sym('bz')
    roots(1) = r * br / bphi
    roots(2) = r * bz / bphi
    do i = 1, 2
        simplified = engine%simplify(roots(i))
        if (.not. simplified%ok) error stop 'field-line simplification failed'
        roots(i) = simplified%value
    end do
    allocate(spec%args(4), spec%outputs(2))
    spec%args = [str('r'), str('br'), str('bphi'), str('bz')]
    spec%outputs = [str('dr_dphi'), str('dz_dphi')]
    spec%name = str('field_line_rhs')
    spec%module_name = str('tiago_field_line_rhs')
    spec%generator = str('gen/gen_field_line.f90')
    spec%generator_revision = str('fortsym dbaa8e319ac5343e8c12a7938ca6a233d672a682')
    spec%regenerate_command = str('build-codegen/gen_field_line (see gen/README.md)')
    spec%pure_procedure = .true.
    code = chars(emit_kernel(roots, spec))
    j = len(code)
    do while (j > 0)
        if (code(j:j) /= new_line('a')) exit
        j = j - 1
    end do
    open(newunit=unit, file='src/generated/field_line_rhs.f90', status='replace')
    write(unit, '(A)') code(:j)
    close(unit)
end program gen_field_line
