program gen_ece_opacity
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, real_expr, pi_expr, sqrt, exp, &
                            operator(+), operator(-), operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: p(8), theta, wc, wp2, r, z, n2, alpha, f5norm, bcor, a, b, &
                    nc2, pol, denom, f9loss
    character(16) :: names(8)
    call arena%init()
    names = [character(16) :: 'te', 'ne', 'bfield', 'frequency', 'f5re', 'f5im', 'f7re', 'f7im']
    call variables(8)
    theta = p(1)/real_expr(arena, 510998.95069_dp)
    wc = real_expr(arena, 1.602176634e-19_dp/9.1093837139e-31_dp)*p(3)
    wp2 = real_expr(arena, 1.602176634e-19_dp**2/ &
                    (8.8541878188e-12_dp*9.1093837139e-31_dp))*p(2)
    r = wp2/wc**2
    z = (1 - wc/(2*pi_expr(arena)*p(4)))/theta
    call emit('ece_opacity_coordinates', 4, [theta, wc, r, z])
    z = (1 - 2*wc/(2*pi_expr(arena)*p(4)))/theta
    call emit('ece_opacity_coordinates2', 4, [theta, wc, r, z])

    ! Bornatici et al. (1983), Eqs. (3.1.25), (3.1.40a); weak damping.
    denom = 1 + r*p(7)/2
    n2 = (1 - r)/denom
    alpha = sqrt(n2)*r*wc*(-p(8))/(2*real_expr(arena, 299792458.0_dp)*denom)
    call emit('ece_opacity_o1', 8, [alpha, n2])

    ! Eq. (3.1.41), n=2. The absorptive F_9/2 is elementary on z<0.
    n2 = 1 - r/4
    f9loss = real_expr(arena, 0.270088205852269108921625521271_dp)*(-z)** &
             real_expr(arena, 3.5_dp)*exp(z)
    alpha = n2**real_expr(arena, 1.5_dp)*r*theta*wc/ &
            real_expr(arena, 299792458.0_dp)*f9loss
    call emit('ece_opacity_o2', 8, [alpha, n2])
    call emit('ece_opacity_o2_blue', 8, [real_expr(arena, 0.0_dp), n2])

    ! Eqs. (3.1.28), (3.1.15a): fundamental X in the tenuous regime.
    alpha = r*wc*(-p(6))/(2*theta*real_expr(arena, 299792458.0_dp))
    call emit('ece_opacity_x1_tenuous', 8, [alpha, real_expr(arena, 1.0_dp)])

    ! Eqs. (3.1.30)-(3.1.32): finite-density X1, including FLR cancellations.
    z = (1 - wc/(2*pi_expr(arena)*p(4)))/theta
    f5norm = p(5)**2 + p(6)**2
    bcor = (1 - r*p(7))**2 + (r*p(8))**2 + &
           real_expr(arena, 0.8_dp)*r*(-z)*(p(5) - r*(p(5)*p(7) + p(6)*p(8))) + &
           real_expr(arena, 6.0_dp/35.0_dp)*r**2*z**2*f5norm
    call emit('ece_opacity_x1_flr', 8, [bcor])
    n2 = 2 - r
    alpha = bcor*sqrt(real_expr(arena, 2.0_dp))*(1 - r/2)** &
            real_expr(arena, 1.5_dp)/r*theta*wc/real_expr(arena, 299792458.0_dp)* &
            (-p(6))/f5norm
    call emit('ece_opacity_x1_finite', 8, [alpha, n2])

    ! Eqs. (3.1.17)-(3.1.20), (3.1.35)-(3.1.38), Hermitian weak-damping root.
    nc2 = 1 - r*(4 - r)/(4*(3 - r))
    a = -real_expr(arena, 1.5_dp)*r/(3 - r)*p(7)
    b = -2*(1 - r/6)*a
    n2 = 2*nc2/(1 + b + sqrt((1 + b)**2 + 4*a*nc2))
    pol = r*(1 + 3*n2*p(7))/(2*(3 - r*(1 + real_expr(arena, 1.5_dp)*n2*p(7))))
    denom = 1 + r*(1 + pol)**2*p(7)/2
    alpha = sqrt(n2)*(1 + pol)**2/denom*r*wc/ &
            real_expr(arena, 299792458.0_dp)*(-p(8))
    call emit('ece_opacity_x2', 8, [alpha, n2])

    names(:4) = [character(16) :: 'alpha', 'index2', 'te', 'frequency']
    call variables(4)
    alpha = p(1)*p(2)*real_expr(arena, 1.602176634e-19_dp/299792458.0_dp**2)* &
            p(3)*p(4)**2
    call emit('ece_opacity_emission', 4, [alpha])
contains
    subroutine variables(n)
        integer, intent(in) :: n
        integer :: j
        do j = 1, n
            p(j) = sym(arena, trim(names(j)))
        end do
    end subroutine
    subroutine emit(name, n, values)
        character(*), intent(in) :: name
        integer, intent(in) :: n
        type(expr_t), intent(in) :: values(:)
        type(expr_t) :: roots(size(values)*(n + 1))
        type(kernel_spec_t) :: spec
        character(:), allocatable :: code
        character(40) :: ref, shape
        integer :: i, j, k, u, last, m
        m = size(values)
        allocate (spec%args(n), spec%output_references(size(roots)))
        k = 0
        do i = 1, m
            k = k + 1
            roots(k) = values(i)
            write (ref, '(A,I0,A)') 'value(', i, ')'
            spec%output_references(k) = str(trim(ref))
            do j = 1, n
                k = k + 1
                roots(k) = diff(values(i), p(j))
                write (ref, '(A,I0,A,I0,A)') 'row(', i, ',', j, ')'
                spec%output_references(k) = str(trim(ref))
            end do
        end do
        do j = 1, n
            spec%args(j) = str(trim(names(j)))
        end do
        spec%outputs = [str('value'), str('row')]
        write (ref, '(A,I0,A)') '(', m, ')'
        write (shape, '(A,I0,A,I0,A)') '(', m, ',', n, ')'
        spec%output_shapes = [str(trim(ref)), str(trim(shape))]
        spec%name = str(name)
        spec%module_name = str('tiago_'//name)
        spec%generator = str('gen/gen_ece_opacity.f90')
        spec%generator_revision = str('fortsym 21d157b394e775a999a8d9ed0be409eac1eb9574')
        spec%regenerate_command = str('build-codegen/gen_ece_opacity (see gen/README.md)')
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
