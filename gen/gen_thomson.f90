program gen_thomson
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, real_expr, pi_expr, operator(+), &
                            operator(-), operator(*), operator(/), operator(**), exp, sqrt
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: p(16), f, a, chi, d, u2, y, xi, delta
    integer :: i
    character(16) :: argnames(16)
    call arena%init()

    call args(['te    ', 'ratio ', 'mu    ', 't     '])
    a = real_expr(arena, 510998.95069_dp)/p(1)
    delta = (p(2) - 1)**2/(2*(1 - p(3))*p(2))
    chi = sqrt(1 + delta)
    u2 = (1 + p(3))/(1 - p(3))
    y = 1/sqrt(chi**2 + u2)
    xi = chi + p(4)/a
    f = 2*chi/(y**3*a)*(y**2/(xi**2 + u2)**real_expr(arena, 1.5_dp) &
                        - 3/(xi**2 + u2)**real_expr(arena, 2.5_dp))
    call emit('ts_depolarization', f, [1])

    call args(['te     ', 'ratio  ', 'mu     ', 'lambda0', 'k2     ', 'depol  '])
    a = real_expr(arena, 510998.95069_dp)/p(1)
    d = 2*(1 - p(3))*p(2) + (p(2) - 1)**2
    delta = (p(2) - 1)**2/(2*(1 - p(3))*p(2))
    ! Rationalized chi-1 avoids cancellation in the cold limit.
    ! Photon spectrum = energy spectrum times lambda/lambda0.
    f = exp(-a*delta/(sqrt(1 + delta) + 1))*(1 + p(6)) &
        /(2*p(4)*p(5)*p(2)**2*sqrt(d))
    call emit('ts_spectral', f, [1, 5, 6])

    call args(['te', 'k1', 'k2', 'k3'])
    f = -real_expr(arena, 510998.95069_dp)/p(1)**2*(p(3) - (p(2) + p(4))/2)
    call emit('ts_bessel_chain', f, [integer ::])

    call args(['ne      ', 'band    ', 'energy  ', 'lambda0 ', 'omega   ', &
               'length  ', 'optics  ', 'absolute', 'outcoupl', 'spectral', 'gain    '])
    f = p(1)*p(2)*p(3)*p(4)*p(5)*p(6)*p(7)*p(8)*p(9)*p(10) &
        *real_expr(arena, (2.8179403205e-15_dp)**2/(6.62607015e-34_dp*299792458.0_dp))
    call emit('ts_photoelectrons', f, [(i, i=1, 10)])
    f = f*p(11)
    call emit('ts_recorded', f, [(i, i=1, 11)])

    call args(['pressure', 'gas_temp', 'energy  ', 'lambda0 ', 'omega   ', &
               'length  ', 'response', 'absolute', 'outcoupl', 'gain    ', 'alpha   '])
    f = p(1)/(real_expr(arena, 1.380649e-23_dp)*p(2)) &
        *16*pi_expr(arena)**4*p(11)**2/p(4)**4 &
        *p(3)*p(4)/(real_expr(arena, 6.62607015e-34_dp*299792458.0_dp)) &
        *p(5)*p(6)*p(7)*p(8)*p(9)*p(10)
    call emit('ts_rayleigh', f, [(i, i=1, 11)])

    call args(['photons ', 'bg      ', 'gain    ', 'digitvar', 'excess  ', 'bgscale '])
    f = p(3)**2*p(5)*(p(1) + (1 + p(6))*p(2)) + p(4)
    call emit('ts_noise', f, [(i, i=1, 6)])
contains
    subroutine args(names)
        character(len=*), intent(in) :: names(:)
        integer :: j
        do j = 1, size(names)
            p(j) = sym(arena, trim(names(j)))
            argnames(j) = names(j)
        end do
    end subroutine
    subroutine emit(name, expression, active)
        character(len=*), intent(in) :: name
        type(expr_t), intent(in) :: expression
        integer, intent(in) :: active(:)
        type(kernel_spec_t) :: spec
        type(expr_t) :: roots(1 + size(active))
        character(:), allocatable :: code
        character(16) :: deriv
        integer :: j, n, unit, last
        ! The argument count follows the named primitive, independently of AD rows.
        select case (name)
        case ('ts_depolarization', 'ts_bessel_chain'); n = 4
        case ('ts_spectral', 'ts_noise'); n = 6
        case ('ts_photoelectrons'); n = 10
        case default; n = 11
        end select
        allocate (spec%args(n), spec%outputs(size(roots)))
        do j = 1, n
            spec%args(j) = str(trim(argnames(j)))
        end do
        roots(1) = expression
        spec%outputs(1) = str('value')
        do j = 1, size(active)
            roots(j + 1) = diff(expression, p(active(j)))
            write (deriv, '(A,I0)') 'd', active(j)
            spec%outputs(j + 1) = str(trim(deriv))
        end do
        spec%name = str(name)
        spec%module_name = str('tiago_'//name)
        spec%generator = str('gen/gen_thomson.f90')
        spec%generator_revision = str('fortsym 21d157b394e775a999a8d9ed0be409eac1eb9574')
        spec%regenerate_command = str('build-codegen/gen_thomson (see gen/README.md)')
        spec%pure_procedure = .true.
        code = chars(emit_kernel(roots, spec))
        last = len_trim(code)
        do while (last > 0)
            if (code(last:last) /= new_line('a')) exit
            last = last - 1
        end do
        open (newunit=unit, file='src/generated/'//name//'.f90', status='replace')
        write (unit, '(A)') code(:last)
        close (unit)
    end subroutine
end program
