program gen_magnetics
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
    type(expr_t) :: raw, gain, offset, drift, time, origin, phi, para, b, r, mu, pi
    type(expr_t) :: pressure, roots(5)
    call arena%init()
    engine = make_symengine_engine(arena)
    raw = sym(arena, 'raw')
    gain = sym(arena, 'gain')
    offset = sym(arena, 'offset')
    drift = sym(arena, 'drift')
    time = sym(arena, 'time')
    origin = sym(arena, 'origin')
    roots(1) = gain*raw + offset + drift*(time - origin)
    roots(2) = diff(roots(1), raw)
    roots(3) = diff(roots(1), gain)
    roots(4) = diff(roots(1), offset)
    roots(5) = diff(roots(1), drift)
    call emit('integrator_signal', [str('raw'), str('gain'), str('offset'), &
                                    str('drift'), str('time'), str('origin')], &
   [str('signal'), str('d_raw'), str('d_gain'), str('d_offset'), str('d_drift')], roots)
    phi = sym(arena, 'flux')
    para = sym(arena, 'paramagnetic_flux')
    b = sym(arena, 'field')
    r = sym(arena, 'radius')
    mu = sym(arena, 'mu0')
    pi = sym(arena, 'circle_pi')
    ! Circular cylinder relation; the correction is supplied, never inferred here.
    roots(1) = 3*pi*r*b*(para - phi)/mu
    roots(2) = diff(roots(1), phi)
    roots(3) = diff(roots(1), b)
    roots(4) = diff(roots(1), r)
    roots(5) = diff(roots(1), para)
    call emit('diamagnetic_energy', [str('flux'), str('field'), str('radius'), &
                              str('paramagnetic_flux'), str('mu0'), str('circle_pi')], &
  [str('energy'), str('d_flux'), str('d_field'), str('d_radius'), str('d_para')], roots)
    pressure = sym(arena, 'pressure')
    roots(1) = 3*pressure/2
    roots(2) = diff(roots(1), pressure)
    call emit('thermal_energy_density', [str('pressure')], &
              [str('energy'), str('d_pressure')], roots(:2))
contains
    subroutine emit(name, args, outputs, expressions)
        use fortsym_string, only: str_t
        character(*), intent(in) :: name
        type(str_t), intent(in) :: args(:), outputs(:)
        type(expr_t), intent(inout) :: expressions(:)
        type(kernel_spec_t) :: spec
        type(engine_result_t) :: simplified
        character(:), allocatable :: code
        integer :: i, j, unit
        do i = 1, size(expressions)
            simplified = engine%simplify(expressions(i))
            if (.not. simplified%ok) error stop 'magnetics simplification failed'
            expressions(i) = simplified%value
        end do
        spec%args = args
        spec%outputs = outputs
        spec%name = str(name)
        spec%module_name = str('tiago_'//name)
        spec%generator = str('gen/gen_magnetics.f90')
       spec%generator_revision = str('fortsym 21d157b394e775a999a8d9ed0be409eac1eb9574')
        spec%regenerate_command = str('build-codegen/gen_magnetics (see gen/README.md)')
        spec%pure_procedure = .true.
        code = chars(emit_kernel(expressions, spec))
        j = len_trim(code)
        do while (j > 0)
            if (code(j:j) /= new_line('a')) exit
            j = j - 1
        end do
        open (newunit=unit, file='src/generated/'//name//'.f90', status='replace')
        write (unit, '(A)') code(:j)
        close (unit)
    end subroutine
end program
