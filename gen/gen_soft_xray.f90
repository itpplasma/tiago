program gen_soft_xray
    !! Soft X-ray emissivity leaves of tiago_soft_xray:
    !!   sxr_bremsstrahlung          eps = C ne^2 Zeff sqrt(Te) F(Te), F and dF/dTe
    !!                               given, C = 1.69e-38 W m^3 eV^-1/2
    !!                               (NRL Plasma Formulary, Gaunt factor 1);
    !!   sxr_filter_segment          one interval [a, a + h] of a piecewise-linear
    !!                               filter transmission t(E):
    !!                               int t(E) exp(-E / Te) dE / Te;
    !!   sxr_filter_segment_small    the same for x = h / Te < 0.01, as a series
    !!                               free of the cancellation in 1 - (1 + x) e^-x;
    !!   sxr_filter_tail             constant transmission t above E: t exp(-E / Te).
    !! Every leaf returns its value and the derivative with respect to Te.
    use fortsym_string, only: str, str_t
    use fortsym_expr, only: expr_t, sym, real_expr, operator(+), operator(-), &
        operator(*), operator(/), operator(**), exp, sqrt
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none
    type(symengine_engine_t) :: engine
    type(expr_t) :: ne, te, zeff, f, df, a, h, ta, tb, x, e1, g, roots(4)

    engine = make_symengine_engine()
    ne = sym('ne')
    te = sym('te')
    zeff = sym('zeff')
    f = sym('filter')
    df = sym('dfilter_dte')
    roots(1) = real_expr(1.69e-38_8)*ne**2*zeff*sqrt(te)*f
    roots(2) = diff(roots(1), ne)
    roots(3) = diff(roots(1), te) + diff(roots(1), f)*df
    roots(4) = diff(roots(1), zeff)
    call emit('sxr_bremsstrahlung', [str('ne'), str('te'), str('zeff'), str('filter'), &
        str('dfilter_dte')], [str('eps'), str('deps_dne'), str('deps_dte'), &
        str('deps_dzeff')], roots)

    a = sym('lower')
    h = sym('width')
    ta = sym('t_lower')
    tb = sym('t_upper')
    x = h/te
    e1 = exp(-a/te)
    roots(1) = e1*(ta*(1 - exp(-x)) + (tb - ta)*(1 - (1 + x)*exp(-x))/x)
    roots(2) = diff(roots(1), te)
    ! Unsimplified: SymEngine would combine exp(-x) with exp(x), which overflows
    ! for x = h / Te >> 1.
    call emit('sxr_filter_segment', [str('lower'), str('width'), str('t_lower'), &
        str('t_upper'), str('te')], [str('f'), str('df_dte')], roots(:2), .false.)

    ! 1 - exp(-x) and (1 - (1 + x) exp(-x)) / x to x^6; the first omitted terms
    ! are below 2e-16 relative for x < 0.01.
    g = x - x**2/2 + x**3/6 - x**4/24 + x**5/120 - x**6/720
    roots(1) = e1*(ta*g + (tb - ta)*(x/2 - x**2/3 + x**3/8 - x**4/30 + x**5/144 &
        - x**6/840))
    roots(2) = diff(roots(1), te)
    call emit('sxr_filter_segment_small', [str('lower'), str('width'), str('t_lower'), &
        str('t_upper'), str('te')], [str('f'), str('df_dte')], roots(:2))

    roots(1) = sym('t_edge')*exp(-sym('edge')/te)
    roots(2) = diff(roots(1), te)
    call emit('sxr_filter_tail', [str('edge'), str('t_edge'), str('te')], &
        [str('f'), str('df_dte')], roots(:2))
contains
    subroutine emit(name, args, outputs, expressions, simplify)
        use fortsym_string, only: chars
        character(len=*), intent(in) :: name
        type(str_t), intent(in) :: args(:), outputs(:)
        type(expr_t), intent(inout) :: expressions(:)
        logical, intent(in), optional :: simplify
        type(kernel_spec_t) :: spec
        type(engine_result_t) :: simplified
        character(len=:), allocatable :: code
        integer :: i, unit
        do i = 1, size(expressions)
            if (present(simplify)) then
                if (.not. simplify) exit
            end if
            simplified = engine%simplify(expressions(i))
            if (.not. simplified%ok) error stop 'soft X-ray simplification failed'
            expressions(i) = simplified%value
        end do
        spec%args = args
        spec%outputs = outputs
        spec%name = str(name)
        spec%module_name = str('tiago_'//name)
        spec%generator = str('gen/gen_soft_xray.f90')
        spec%generator_revision = &
            str('fortsym f579ad118f552b09d4dd6385ebd753ec1500922a')
        spec%regenerate_command = str('build-codegen/gen_soft_xray (see gen/README.md)')
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
end program gen_soft_xray
