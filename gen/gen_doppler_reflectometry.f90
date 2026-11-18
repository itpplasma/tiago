program gen_doppler_reflectometry
    use fortsym_string, only: str, chars
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), operator(*), operator(/)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    implicit none
    type(kernel_spec_t) :: spec
    type(expr_t) :: x(5), roots(6)
    character(len=:), allocatable :: code
    integer :: i, unit

    x = [sym('radial_field'), sym('phase_velocity'), sym('wavenumber'), &
         sym('b_squared'), sym('b_projection')]
    roots(1) = x(3)*(x(1)*x(5)/x(4) + x(2))
    do i = 1, size(x)
        roots(i + 1) = diff(roots(1), x(i))
    end do
    spec%args = [str('radial_field'), str('phase_velocity'), str('wavenumber'), &
                 str('b_squared'), str('b_projection')]
    spec%outputs = [str('omega'), str('d_radial_field'), str('d_phase_velocity'), &
                    str('d_wavenumber'), str('d_b_squared'), str('d_b_projection')]
    spec%name = str('doppler_reflectometry')
    spec%module_name = str('tiago_doppler_reflectometry')
    spec%generator = str('gen/gen_doppler_reflectometry.f90')
    spec%regenerate_command = str('build-codegen/gen_doppler_reflectometry (see gen/README.md)')
    spec%pure_procedure = .true.
    code = chars(emit_kernel(roots, spec))
    open (newunit=unit, file='src/generated/doppler_reflectometry.f90', status='replace')
    i = len(code)
    do while (i > 0)
        if (code(i:i) /= new_line('a')) exit
        i = i - 1
    end do
    write (unit, '(A)') code(:i)
    close (unit)
end program gen_doppler_reflectometry
