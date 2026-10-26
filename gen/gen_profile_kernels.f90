program gen_profile_kernels
    !! Pointwise kernels of the profile diagnostics:
    !!   hermite_basis   cubic Hermite basis jets on [0, 1];
    !!   circular_state  shifted-circle tokamak map, |B| = B0 R0 / R and their
    !!                   derivatives in the labels (rho, theta, phi) and the
    !!                   parameters (R0, a, Delta, B0), the analytic test state;
    !!   cold_resonance  |B| of the n-th electron-cyclotron harmonic at f;
    !!   xics_emissivity coronal line emissivity and its partial derivatives;
    !!   student_t_prior negative log density (up to a constant) of a unit
    !!                   Student-t prior and its first two derivatives.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, real_expr, pi_expr, operator(+), &
        operator(-), operator(*), operator(/), operator(**), sin, cos, exp, sqrt, log
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none
    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(engine_result_t) :: simplified
    type(expr_t) :: t, u(3), p(4), x(3), b, roots(40), f, n, ne, te, zeff, zi, el, zt, nu
    type(kernel_spec_t) :: spec
    character(len=5), parameter :: unames(3) = ['rho  ', 'theta', 'phi  ']
    character(len=5), parameter :: pnames(4) = ['r0   ', 'a    ', 'delta', 'b0   ']
    character(len=3), parameter :: hnames(8) = &
        ['h00', 'h10', 'h01', 'h11', 'd00', 'd10', 'd01', 'd11']
    character(len=16) :: ref
    type(expr_t) :: r
    integer :: i, j, k, unit

    call arena%init()
    engine = make_symengine_engine(arena)

    ! Cubic Hermite basis and its t-derivative.
    t = sym(arena, 't')
    roots(1) = 2*t**3 - 3*t**2 + 1
    roots(2) = t**3 - 2*t**2 + t
    roots(3) = -2*t**3 + 3*t**2
    roots(4) = t**3 - t**2
    do i = 1, 4
        roots(4 + i) = diff(roots(i), t)
    end do
    call reset_spec(1, 8)
    spec%args(1) = str('t')
    do i = 1, 8
        spec%outputs(i) = str(hnames(i))
    end do
    call emit('hermite_basis', roots(:8))

    ! Shifted circles: R = R0 + Delta (1 - rho^2) + a rho cos(theta),
    ! Z = a rho sin(theta), toroidal angle phi; |B| = B0 R0 / R.
    do i = 1, 3
        u(i) = sym(arena, trim(unames(i)))
    end do
    do i = 1, 4
        p(i) = sym(arena, trim(pnames(i)))
    end do
    r = p(1) + p(3)*(1 - u(1)**2) + p(2)*u(1)*cos(u(2))
    x(1) = r*cos(u(3))
    x(2) = r*sin(u(3))
    x(3) = p(2)*u(1)*sin(u(2))
    b = p(4)*p(1)/r
    call reset_spec(7, 0)
    do i = 1, 3
        spec%args(i) = str(trim(unames(i)))
    end do
    do i = 1, 4
        spec%args(3 + i) = str(trim(pnames(i)))
    end do
    spec%outputs = [str('x'), str('dxdu'), str('dxdp'), str('b'), str('dbdu'), str('dbdp')]
    spec%output_shapes = [str('(3)'), str('(3, 3)'), str('(3, 4)'), str(''), str('(3)'), &
        str('(4)')]
    allocate(spec%output_references(32))
    k = 0
    do i = 1, 3
        k = k + 1
        roots(k) = x(i)
        write(ref, '(A,I0,A)') 'x(', i, ')'
        spec%output_references(k) = str(trim(ref))
    end do
    do j = 1, 3
        do i = 1, 3
            k = k + 1
            roots(k) = diff(x(i), u(j))
            write(ref, '(A,I0,A,I0,A)') 'dxdu(', i, ', ', j, ')'
            spec%output_references(k) = str(trim(ref))
        end do
    end do
    do j = 1, 4
        do i = 1, 3
            k = k + 1
            roots(k) = diff(x(i), p(j))
            write(ref, '(A,I0,A,I0,A)') 'dxdp(', i, ', ', j, ')'
            spec%output_references(k) = str(trim(ref))
        end do
    end do
    k = k + 1
    roots(k) = b
    spec%output_references(k) = str('b')
    do j = 1, 3
        k = k + 1
        roots(k) = diff(b, u(j))
        write(ref, '(A,I0,A)') 'dbdu(', j, ')'
        spec%output_references(k) = str(trim(ref))
    end do
    do j = 1, 4
        k = k + 1
        roots(k) = diff(b, p(j))
        write(ref, '(A,I0,A)') 'dbdp(', j, ')'
        spec%output_references(k) = str(trim(ref))
    end do
    call emit('circular_state', roots(:k))

    ! Cold resonance n e B / (2 pi m_e) = f, CODATA 2018 e and m_e.
    f = sym(arena, 'f')
    n = sym(arena, 'harmonic')
    roots(1) = 2*pi_expr(arena)*real_expr(arena, 9.1093837015e-31_dp)*f / &
        (n*real_expr(arena, 1.602176634e-19_dp))
    roots(2) = diff(roots(1), f)
    call reset_spec(2, 2)
    spec%args = [str('f'), str('harmonic')]
    spec%outputs = [str('b_res'), str('db_res_df')]
    call emit('cold_resonance', roots(:2))

    ! Coronal excitation of an impurity line, up to a constant:
    ! eps = ne n_Z exp(-E / Te) / sqrt(Te),  n_Z = ne (Zeff - 1) / (Z (Z - 1)).
    ne = sym(arena, 'ne')
    te = sym(arena, 'te')
    zeff = sym(arena, 'zeff')
    zi = sym(arena, 'z_imp')
    el = sym(arena, 'e_line')
    roots(1) = ne*ne*(zeff - 1)/(zi*(zi - 1))*exp(-el/te)/sqrt(te)
    roots(2) = diff(roots(1), ne)
    roots(3) = diff(roots(1), te)
    roots(4) = diff(roots(1), zeff)
    call reset_spec(5, 4)
    spec%args = [str('ne'), str('te'), str('zeff'), str('z_imp'), str('e_line')]
    spec%outputs = [str('eps'), str('deps_dne'), str('deps_dte'), str('deps_dzeff')]
    call emit('xics_emissivity', roots(:4))

    ! Unit Student-t prior with nu degrees of freedom, -log density up to its
    ! normalization: phi = (nu + 1) / 2 log(1 + z^2 / nu).
    zt = sym(arena, 'z')
    nu = sym(arena, 'nu')
    roots(1) = (nu + 1)/2*log(1 + zt*zt/nu)
    roots(2) = diff(roots(1), zt)
    roots(3) = diff(roots(2), zt)
    call reset_spec(2, 3)
    spec%args = [str('z'), str('nu')]
    spec%outputs = [str('phi'), str('dphi_dz'), str('d2phi_dz2')]
    call emit('student_t_prior', roots(:3))
contains
    subroutine reset_spec(nargs, nout)
        integer, intent(in) :: nargs, nout
        type(kernel_spec_t) :: fresh
        spec = fresh
        allocate(spec%args(nargs), spec%outputs(nout))
    end subroutine reset_spec

    subroutine emit(name, expressions)
        character(len=*), intent(in) :: name
        type(expr_t), intent(inout) :: expressions(:)
        integer :: jj
        character(len=:), allocatable :: code
        do jj = 1, size(expressions)
            simplified = engine%simplify(expressions(jj))
            if (.not. simplified%ok) error stop 'simplification failed'
            expressions(jj) = simplified%value
        end do
        spec%name = str(name)
        spec%module_name = str('tiago_'//name)
        spec%generator = str('gen/gen_profile_kernels.f90')
        spec%generator_revision = &
            str('fortsym 48dda45a34c0da0db1c78cd4130b3a97375e3c8e')
        spec%regenerate_command = str('build-codegen/gen_profile_kernels (see gen/README.md)')
        spec%pure_procedure = .true.
        open (newunit=unit, file='src/generated/'//name//'.f90', status='replace')
        code = chars(emit_kernel(expressions, spec))
        jj = len(code)
        do while (jj > 0)
            if (code(jj:jj) /= new_line('a')) exit
            jj = jj - 1
        end do
        write (unit, '(A)') code(:jj)
        close (unit)
    end subroutine emit
end program gen_profile_kernels
