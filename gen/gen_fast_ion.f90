program gen_fast_ion
    !! Pointwise leaves of the fast-ion diagnostics (tiago_fusion_rates):
    !!   dd_cross_section  D(d,n)3He cross-section sigma(E) [mb] at the
    !!                     centre-of-mass energy E [keV], Bosch and Hale,
    !!                     Nucl. Fusion 32, 611 (1992), Eqs. (8), (9), Table IV;
    !!   dd_reactivity     Maxwellian D(d,n)3He reactivity <sigma v>(T) [cm^3/s]
    !!                     at the ion temperature T [keV], same paper, Eq. (12),
    !!                     Table VII;
    !!   cx_hydrogen       total charge-exchange cross-section [cm^2] of
    !!                     H+ + H(n = 1) at the collision energy E [keV/amu],
    !!                     Janev and Smith (1993) Eq. 44, Table 9, with the
    !!                     coefficients of FIDASIM's p_cx_1_janev.
    !! Every leaf returns its value and the derivative with respect to its argument.
    use fortsym_string, only: str, str_t
    use fortsym_expr, only: expr_t, sym, real_expr, rat, log, operator(+), operator(-), &
        operator(*), operator(/), operator(**), exp, sqrt
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none
    type(symengine_engine_t) :: engine
    type(expr_t) :: e, t, s, theta, xi, roots(2)
    real(8), parameter :: bg = 31.3970_8

    engine = make_symengine_engine()
    e = sym('e')
    s = real_expr(5.3701e4_8) + e*(real_expr(3.3027e2_8) + e*(real_expr(-1.2706e-1_8) &
        + e*(real_expr(2.9327e-5_8) + e*real_expr(-2.5151e-9_8))))
    roots(1) = s/(e*exp(real_expr(bg)/sqrt(e)))
    roots(2) = diff(roots(1), e)
    call emit('dd_cross_section', [str('e')], [str('sigma'), str('dsigma_de')], roots, &
        .false.)

    t = sym('t')
    theta = t/(1 - t*real_expr(5.85778e-3_8)/(1 + t*(real_expr(7.68222e-3_8) &
        + t*real_expr(-2.964e-6_8))))
    xi = (real_expr(bg**2/4)/theta)**rat(1, 3)
    roots(1) = real_expr(5.43360e-12_8)*theta*sqrt(xi/(real_expr(937814.0_8)*t**3)) &
        *exp(-3*xi)
    roots(2) = diff(roots(1), t)
    call emit('dd_reactivity', [str('t')], [str('sv'), str('dsv_dt')], roots, .false.)

    roots(1) = real_expr(3.2345e-16_8)*log(real_expr(2.3588e2_8)/e + real_expr(2.3713_8)) &
        /(1 + real_expr(3.8371e-2_8)*e + real_expr(3.8068e-6_8)*e**rat(7, 2) &
        + real_expr(1.1832e-10_8)*e**rat(27, 5))
    roots(2) = diff(roots(1), e)
    call emit('cx_hydrogen', [str('e')], [str('sigma'), str('dsigma_de')], roots, .false.)
contains
    subroutine emit(name, args, outputs, expressions, simplify)
        use fortsym_string, only: chars
        character(len=*), intent(in) :: name
        type(str_t), intent(in) :: args(:), outputs(:)
        type(expr_t), intent(inout) :: expressions(:)
        logical, intent(in) :: simplify
        type(kernel_spec_t) :: spec
        type(engine_result_t) :: simplified
        character(len=:), allocatable :: code
        integer :: i, unit
        if (simplify) then
            do i = 1, size(expressions)
                simplified = engine%simplify(expressions(i))
                if (.not. simplified%ok) error stop 'fast-ion simplification failed'
                expressions(i) = simplified%value
            end do
        end if
        spec%args = args
        spec%outputs = outputs
        spec%name = str(name)
        spec%module_name = str('tiago_'//name)
        spec%generator = str('gen/gen_fast_ion.f90')
        spec%generator_revision = &
            str('fortsym 257ea66f6071f516a02f5400eeb32e6abac6e9d4')
        spec%regenerate_command = str('build-codegen/gen_fast_ion (see gen/README.md)')
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
end program gen_fast_ion
