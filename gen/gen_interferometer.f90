program gen_interferometer
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, real_expr, pi_expr, operator(+), &
                            operator(-), operator(*), operator(/), sin, cos, atan2
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none
    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(expr_t) :: v(8), f(2), re, pi, electron, gas
    character(len=16) :: names(8)
    integer :: i

    call arena%init()
    engine = make_symengine_engine(arena)
    ! CODATA 2022 classical electron radius [m], fixed physical constant.
    re = real_expr(arena, 2.8179403205e-15_dp)
    pi = pi_expr(arena)
    names = [character(len=16) :: 'column', 'wavelength', 'neutral', &
             'neutral_scale', 'environment', 'env_gain', 'phase0', 'path']
    call variables(8)
    electron = 3*re*v(2)*v(1)/2
    gas = 4*pi*v(3)*v(4)/v(2)
    ! phi(2 omega) - 2 phi(omega): the common vacuum path cancels.
    f(1) = electron - gas + v(6)*v(5) + v(7) &
           + 4*pi*v(8)/v(2) - 4*pi*v(8)/v(2)
    call emit('dispersion_phase', 8, f(:1))
    names(1) = 'phase'
    call variables(7)
    ! Rebuild after changing the first symbol.
    f(1) = (v(1) + 4*pi*v(3)*v(4)/v(2) - v(6)*v(5) - v(7)) &
           /(3*re*v(2)/2)
    call emit('dispersion_column', 7, f(:1))

    names(:5) = [character(len=16) :: 'phase', 'ai', 'aq', 'oi', 'oq']
    call variables(5)
    f(1) = v(4) + v(2)*cos(v(1))
    f(2) = v(5) + v(3)*sin(v(1))
    call emit('dispersion_iq', 5, f)

    names(:6) = [character(len=16) :: 'ival', 'qval', 'ai', 'aq', 'oi', 'oq']
    call variables(6)
    f(1) = atan2((v(2) - v(6))/v(4), (v(1) - v(5))/v(3))
    call emit('dispersion_atan', 6, f(:1))
contains
    subroutine variables(n)
        integer, intent(in) :: n
        do i = 1, n
            v(i) = sym(arena, trim(names(i)))
        end do
    end subroutine

    subroutine emit(name, n, values)
        character(len=*), intent(in) :: name
        integer, intent(in) :: n
        type(expr_t), intent(in) :: values(:)
        type(kernel_spec_t) :: spec
        type(expr_t), allocatable :: roots(:)
        type(engine_result_t) :: simplified
        character(len=:), allocatable :: code
        character(len=32) :: ref, shape_value, shape_jac
        integer :: j, k, m, r, u, last, nargs
        m = size(values)
        nargs = n
        ! Common path cancels symbolically; retain its zero derivative row.
        if (name == 'dispersion_phase') nargs = n - 1
        allocate (spec%args(nargs), roots(m*(n + 1)))
        do j = 1, nargs
            spec%args(j) = str(trim(names(j)))
        end do
        spec%outputs = [str('value'), str('jac')]
        write (shape_value, '(A,I0,A)') '(', m, ')'
        write (shape_jac, '(A,I0,A,I0,A)') '(', m, ',', n, ')'
        spec%output_shapes = [str(trim(shape_value)), str(trim(shape_jac))]
        allocate (spec%output_references(size(roots)))
        do j = 1, m
            roots(j) = values(j)
            write (ref, '(A,I0,A)') 'value(', j, ')'
            spec%output_references(j) = str(trim(ref))
        end do
        r = m
        do k = 1, n
            do j = 1, m
                r = r + 1
                roots(r) = diff(values(j), v(k))
                write (ref, '(A,I0,A,I0,A)') 'jac(', j, ',', k, ')'
                spec%output_references(r) = str(trim(ref))
            end do
        end do
        do j = 1, size(roots)
            simplified = engine%simplify(roots(j))
            if (.not. simplified%ok) error stop 'interferometer simplification'
            roots(j) = simplified%value
        end do
        spec%name = str(name)
        spec%module_name = str('tiago_'//name)
        spec%generator = str('gen/gen_interferometer.f90')
        spec%generator_revision = &
            str('fortsym 21d157b394e775a999a8d9ed0be409eac1eb9574')
        spec%regenerate_command = &
            str('build-codegen/gen_interferometer (see gen/README.md)')
        spec%pure_procedure = .true.
        code = chars(emit_kernel(roots, spec))
        last = len_trim(code)
        do while (last > 0)
            if (code(last:last) /= new_line('a')) exit
            last = last - 1
        end do
        open (newunit=u, file='src/generated/'//name//'.f90', status='replace')
        write (u, '(A)') code(:last)
        close (u)
    end subroutine
end program
