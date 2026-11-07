program gen_cxrs
    !! Charge-exchange emission leaf of tiago_cxrs: a beam neutral density with a
    !! transverse profile g and attenuation exp(-a), times the fully stripped
    !! impurity density n_Z = ne (Zeff - 1) / (Z (Z - 1)) of a single impurity.
    !! The charge-exchange rate coefficient is constant per channel and absorbed
    !! in the brightness scale.
    use fortsym_string, only: str, chars
    use fortsym_expr, only: expr_t, sym, operator(-), operator(*), operator(/), exp
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none
    type(symengine_engine_t) :: engine
    type(expr_t) :: ne, zeff, a, g, zi, roots(4)

    engine = make_symengine_engine()
    ne = sym('ne')
    zeff = sym('zeff')
    a = sym('attenuation')
    g = sym('transverse')
    zi = sym('z_imp')
    roots(1) = g*exp(-a)*ne*(zeff - 1)/(zi*(zi - 1))
    roots(2) = diff(roots(1), ne)
    roots(3) = diff(roots(1), zeff)
    roots(4) = diff(roots(1), a)
    call emit('cxrs_emission', [str('ne'), str('zeff'), str('attenuation'), &
        str('transverse'), str('z_imp')], [str('em'), str('dem_dne'), &
        str('dem_dzeff'), str('dem_da')], roots)
contains
    subroutine emit(name, args, outputs, expressions)
        use fortsym_string, only: str_t
        character(len=*), intent(in) :: name
        type(str_t), intent(in) :: args(:), outputs(:)
        type(expr_t), intent(inout) :: expressions(:)
        type(kernel_spec_t) :: spec
        type(engine_result_t) :: simplified
        character(len=:), allocatable :: code
        integer :: i, unit
        do i = 1, size(expressions)
            simplified = engine%simplify(expressions(i))
            if (.not. simplified%ok) error stop 'cxrs simplification failed'
            expressions(i) = simplified%value
        end do
        spec%args = args
        spec%outputs = outputs
        spec%name = str(name)
        spec%module_name = str('tiago_'//name)
        spec%generator = str('gen/gen_cxrs.f90')
        spec%generator_revision = &
            str('fortsym bca0bc4e9ebe3eeadf3a8d1c2afc821f95f6b09c')
        spec%regenerate_command = str('build-codegen/gen_cxrs (see gen/README.md)')
        spec%pure_procedure = .true.
        open (newunit=unit, file='src/generated/'//name//'.f90', status='replace')
        code = chars(emit_kernel(expressions, spec))
        i = len(code)
        do while (i > 0)
            if (code(i:i) /= new_line('a')) exit
            i = i - 1
        end do
        write (unit, '(A)') code(:i)
        close (unit)
    end subroutine emit
end program gen_cxrs
