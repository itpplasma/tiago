program gen_reflectometry
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortsym_string, only: str, chars
    use fortsym_expr, only: expr_t, sym, real_expr, pi_expr, operator(+), &
        operator(-), operator(*), operator(/), operator(**), sqrt
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    implicit none
    real(dp), parameter :: light = 299792458.0_dp
    real(dp), parameter :: charge = 1.602176634e-19_dp
    real(dp), parameter :: eps0 = 8.8541878128e-12_dp
    real(dp), parameter :: electron_mass = 9.1093837015e-31_dp
    type(expr_t) :: v(12), k0, nb, b2, nperp2, cross_n(3), projection
    type(expr_t) :: omega, wavenumber, xx, index
    type(expr_t), allocatable :: roots(:)
    type(kernel_spec_t) :: spec
    character(len=6), parameter :: names(12) = [character(len=6) :: &
        'nx', 'ny', 'nz', 'bx', 'by', 'bz', 'gx', 'gy', 'gz', 'freq', 'er', 'vph']
    character(len=40) :: ref
    character(len=:), allocatable :: code
    integer :: i, unit

    ! Turning-point Doppler shift. N is the refractive-index vector, B the
    ! field, G the density gradient. Backscattering selects q = -2 k0 N_perp
    ! and the outward normal is -G/|G|; |N_perp| cancels in the E x B term.
    do i = 1, 12
        v(i) = sym(trim(names(i)))
    end do
    k0 = 2*pi_expr()*v(10)/real_expr(light)
    nb = v(1)*v(4) + v(2)*v(5) + v(3)*v(6)
    b2 = v(4)**2 + v(5)**2 + v(6)**2
    nperp2 = v(1)**2 + v(2)**2 + v(3)**2 - nb**2/b2
    cross_n(1) = -(v(8)*v(6) - v(9)*v(5))
    cross_n(2) = -(v(9)*v(4) - v(7)*v(6))
    cross_n(3) = -(v(7)*v(5) - v(8)*v(4))
    projection = (v(1)*cross_n(1) + v(2)*cross_n(2) + v(3)*cross_n(3))/ &
        sqrt(v(7)**2 + v(8)**2 + v(9)**2)
    wavenumber = -2*k0*sqrt(nperp2)
    omega = -2*k0*(v(11)*projection/b2 + v(12)*sqrt(nperp2))
    allocate (roots(26))
    roots(1) = omega
    roots(2) = wavenumber
    do i = 1, 12
        roots(2 + i) = diff(omega, v(i))
        roots(14 + i) = diff(wavenumber, v(i))
    end do
    call reset(12)
    do i = 1, 12
        spec%args(i) = str(trim(names(i)))
    end do
    spec%outputs = [str('omega'), str('wavenumber'), str('d_omega'), &
        str('d_wavenumber')]
    spec%output_shapes = [str(''), str(''), str('(12)'), str('(12)')]
    allocate (spec%output_references(26))
    spec%output_references(1) = str('omega')
    spec%output_references(2) = str('wavenumber')
    do i = 1, 12
        write (ref, '(A,I0,A)') 'd_omega(', i, ')'
        spec%output_references(2 + i) = str(trim(ref))
        write (ref, '(A,I0,A)') 'd_wavenumber(', i, ')'
        spec%output_references(14 + i) = str(trim(ref))
    end do
    call emit('doppler_turning_shift')
    deallocate (roots)

    ! Ordinary-mode index N = sqrt(1 - X) and its inverse, the phase and
    ! group-delay integrands, with derivatives in density and frequency.
    v(1) = sym('ne')
    v(2) = sym('freq')
    xx = real_expr(charge**2/(eps0*electron_mass))*v(1)/(2*pi_expr()*v(2))**2
    index = sqrt(1 - xx)
    allocate (roots(7))
    roots(1) = xx
    roots(2) = index
    roots(3) = diff(index, v(1))
    roots(4) = diff(index, v(2))
    roots(5) = 1/index
    roots(6) = diff(1/index, v(1))
    roots(7) = real_expr(eps0*electron_mass/charge**2)*(2*pi_expr()*v(2))**2
    call reset(2)
    spec%args(1) = str('ne')
    spec%args(2) = str('freq')
    spec%outputs = [str('x_ratio'), str('n_index'), str('dn_dne'), &
        str('dn_dfreq'), str('inverse_index'), str('dinverse_dne'), &
        str('cutoff_density')]
    call emit('omode_index')
    deallocate (roots)

contains
    subroutine reset(nargs)
        integer, intent(in) :: nargs
        spec = kernel_spec_t()
        allocate (spec%args(nargs))
    end subroutine
    subroutine emit(kernel_name)
        character(len=*), intent(in) :: kernel_name
        integer :: n
        spec%name = str(kernel_name)
        spec%module_name = str('tiago_'//kernel_name)
        spec%generator = str('gen/gen_reflectometry.f90')
        spec%regenerate_command = &
            str('build-codegen/gen_reflectometry (see gen/README.md)')
        spec%pure_procedure = .true.
        code = chars(emit_kernel(roots, spec))
        n = len(code)
        do while (n > 0)
            if (code(n:n) /= new_line('a')) exit
            n = n - 1
        end do
        open (newunit=unit, file='src/generated/'//kernel_name//'.f90', &
            status='replace')
        write (unit, '(A)') code(:n)
        close (unit)
    end subroutine
end program gen_reflectometry
