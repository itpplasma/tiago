program gen_axisym_green
    ! Flux per radian and poloidal field of a unit toroidal filament at (rc, zc),
    ! divided by mu0, with K and E supplied at m = 4 r rc / D^2.
    use fortsym_string, only: str, str_t, chars
    use fortsym_expr, only: expr_t, sym, pi_expr, sqrt, operator(+), &
        operator(-), operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    implicit none
    type(expr_t) :: r, z, rc, zc, kk, ee, d2, m, q, g, dk, de, roots(3), modulus(2)
    r = sym('r')
    z = sym('z')
    rc = sym('rc')
    zc = sym('zc')
    kk = sym('kk')
    ee = sym('ee')
    d2 = (r + rc)**2 + (z - zc)**2
    m = 4*r*rc/d2
    ! q = 1 - m without cancellation, for K and E near the filament.
    q = ((r - rc)**2 + (z - zc)**2)/d2
    modulus(1) = m
    modulus(2) = q
    call emit('axisym_green_modulus', [str('r'), str('z'), str('rc'), str('zc')], &
        [str('m'), str('q')], modulus)
    g = sqrt(d2)/(4*pi_expr())*((2 - m)*kk - 2*ee)
    ! dK/dm and dE/dm in terms of K and E.
    dk = (ee - q*kk)/(2*m*q)
    de = (ee - kk)/(2*m)
    roots(1) = g
    roots(2) = -(diff(g, z) + (diff(g, kk)*dk + diff(g, ee)*de)*diff(m, z))/r
    roots(3) = (diff(g, r) + (diff(g, kk)*dk + diff(g, ee)*de)*diff(m, r))/r
    call emit('axisym_green', [str('r'), str('z'), str('rc'), str('zc'), str('kk'), &
        str('ee')], [str('psi'), str('br'), str('bz')], roots)
contains
    subroutine emit(name, args, outputs, expressions)
        character(*), intent(in) :: name
        type(str_t), intent(in) :: args(:), outputs(:)
        type(expr_t), intent(in) :: expressions(:)
        type(kernel_spec_t) :: spec
        character(:), allocatable :: code
        integer :: j, unit
        spec%args = args
        spec%outputs = outputs
        spec%name = str(name)
        spec%module_name = str('tiago_'//name)
        spec%generator = str('gen/gen_axisym_green.f90')
        spec%generator_revision = str('fortsym 257ea66f6071f516a02f5400eeb32e6abac6e9d4')
        spec%regenerate_command = str('build-codegen/gen_axisym_green (see gen/README.md)')
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
