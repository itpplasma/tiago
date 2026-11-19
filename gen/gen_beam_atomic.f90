program gen_beam_atomic
    !! Atomic kernels of tiago_beam_neutrals.
    !!
    !! beam_stopping_suzuki: stopping cross-section of a hydrogen beam, Suzuki
    !! et al., Plasma Phys. Control. Fusion 40 (1998) 2097, equations (24),
    !! (26) and (28), with its partial derivatives in ne, Te and Zeff.
    !!
    !!   sigma_H = A1 1e-16 / E (1 + A2 e + A3 e^2)
    !!             {1 + [1 - exp(-A4 N)]^A5 (A6 + A7 e + A8 e^2)} (1 + A9 U + A10 U^2)
    !!   S_Z     = sum_ijk B_ijk e^(i-1) (ln N)^(j-1) U^(k-1)
    !!   sigma   = sigma_H [1 + (Zeff - 1) S_Z]  [cm^2, returned in m^2]
    !!
    !! with e = ln(E / keV amu^-1), N = ne / 1e19 m^-3, U = ln(Te / keV).
    !!
    !! hydrogen_1s_rates: H+ + H(1s) charge exchange (eq. 44, table 9) and
    !! ionization (eq. 40, table 8), and e + H(1s) ionization (eq. 14,
    !! table 3) of Janev, Reiter and Samm, Collision processes in
    !! low-temperature hydrogen plasmas, Report Juel-4105 (2003).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortsym_string, only: str, str_t, chars
    use fortsym_expr, only: expr_t, sym, real_expr, operator(+), operator(-), &
        operator(*), operator(/), operator(**), exp, log
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none
    type(symengine_engine_t) :: engine
    type(expr_t) :: e, ne, te, zeff, eps, lnn, nn, u, sh, sz, sigma, a(10), b(12)
    type(expr_t) :: suzuki(4), rates(3), ep, ee, x, s
    type(str_t) :: args(26)
    integer :: i, j, k, m

    engine = make_symengine_engine()
    e = sym('e_amu')
    ne = sym('ne')
    te = sym('te')
    zeff = sym('zeff')
    do i = 1, 10
        a(i) = sym('a'//digits(i))
    end do
    do i = 1, 12
        b(i) = sym('b'//digits(i))
    end do
    eps = log(e)
    nn = ne/real_expr(1.0e19_dp)
    lnn = log(nn)
    u = log(te/real_expr(1.0e3_dp))
    sh = a(1)*real_expr(1.0e-20_dp)/e*(1 + a(2)*eps + a(3)*eps**2) &
        *(1 + (1 - exp(-a(4)*nn))**a(5)*(a(6) + a(7)*eps + a(8)*eps**2)) &
        *(1 + a(9)*u + a(10)*u**2)
    sz = real_expr(0.0_dp)
    m = 0
    do i = 0, 2
        do j = 0, 1
            do k = 0, 1
                m = m + 1
                sz = sz + b(m)*eps**i*lnn**j*u**k
            end do
        end do
    end do
    sigma = sh*(1 + (zeff - 1)*sz)
    suzuki(1) = sigma
    suzuki(2) = diff(sigma, ne)
    suzuki(3) = diff(sigma, te)
    suzuki(4) = diff(sigma, zeff)
    args(1) = str('e_amu')
    args(2) = str('ne')
    args(3) = str('te')
    args(4) = str('zeff')
    do i = 1, 10
        args(4 + i) = str('a'//digits(i))
    end do
    do i = 1, 12
        args(14 + i) = str('b'//digits(i))
    end do
    call emit('beam_stopping_suzuki', args, [str('sigma'), str('dsigma_dne'), &
        str('dsigma_dte'), str('dsigma_dzeff')], suzuki)

    ! Janev, Reiter and Samm (2003): n = 1 columns of tables 9, 8 and 3.
    ep = sym('e_amu')
    rates(1) = real_expr(3.2345e-20_dp)*log(real_expr(2.3588e2_dp)/ep + real_expr(2.3713_dp)) &
        /(1 + real_expr(3.8371e-2_dp)*ep + real_expr(3.8068e-6_dp)*ep**real_expr(3.5_dp) &
        + real_expr(1.1832e-10_dp)*ep**real_expr(5.4_dp))
    rates(2) = real_expr(2.016e-3_dp*1.0e-20_dp)*(ep**real_expr(3.7154_dp) &
        *exp(-real_expr(3.9890e-2_dp)*ep)/(1 + real_expr(3.1413e-1_dp)*ep**real_expr(2.1254_dp)) &
        + real_expr(6.3990e3_dp)*exp(-real_expr(6.1897e1_dp)/ep) &
        *log(1 + real_expr(9.2731e3_dp)*ep)/ep)
    ee = sym('e_electron')
    x = 1 - real_expr(13.6_dp)/ee
    s = real_expr(-0.032226_dp)*x + real_expr(-0.034539_dp)*x**2 + real_expr(1.4003_dp)*x**3 &
        + real_expr(-2.8115_dp)*x**4 + real_expr(2.2986_dp)*x**5
    rates(3) = real_expr(1.0e-17_dp)/(real_expr(13.6_dp)*ee) &
        *(real_expr(0.18450_dp)*log(ee/real_expr(13.6_dp)) + s)
    call emit('hydrogen_1s_cross_sections', [str('e_amu'), str('e_electron')], &
        [str('sigma_cx'), str('sigma_pion'), str('sigma_eion')], rates)
contains
    pure function digits(i) result(text)
        integer, intent(in) :: i
        character(len=:), allocatable :: text
        character(len=8) :: buffer
        write (buffer, '(I0)') i
        text = trim(buffer)
    end function digits

    subroutine emit(name, args, outputs, expressions)
        character(len=*), intent(in) :: name
        type(str_t), intent(in) :: args(:), outputs(:)
        type(expr_t), intent(inout) :: expressions(:)
        type(kernel_spec_t) :: spec
        type(engine_result_t) :: simplified
        character(len=:), allocatable :: code
        integer :: i, unit
        do i = 1, size(expressions)
            simplified = engine%simplify(expressions(i))
            if (.not. simplified%ok) error stop 'beam atomic simplification failed'
            expressions(i) = simplified%value
        end do
        spec%args = args
        spec%outputs = outputs
        spec%name = str(name)
        spec%module_name = str('tiago_'//name)
        spec%generator = str('gen/gen_beam_atomic.f90')
        spec%generator_revision = &
            str('fortsym 257ea66f6071f516a02f5400eeb32e6abac6e9d4')
        spec%regenerate_command = str('build-codegen/gen_beam_atomic (see gen/README.md)')
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
end program gen_beam_atomic
