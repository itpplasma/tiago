program gen_islands
    ! Pointwise kernels of the field-line and island module (src/islands).
    use fortsym_string, only: str, chars
    use fortsym_expr, only: expr_t, sym, sin, cos, atan2, operator(+), &
        operator(-), operator(*), operator(/)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    implicit none
    type(expr_t) :: r, br, bphi, bz, phi, bx, by, a, b, c, d, x, y, m, n

    r = sym('r'); br = sym('br'); bphi = sym('bphi'); bz = sym('bz')
    phi = sym('phi'); bx = sym('bx'); by = sym('by')
    a = sym('a'); b = sym('b'); c = sym('c'); d = sym('d')
    x = sym('x'); y = sym('y'); m = sym('m'); n = sym('n')
    ! dR/dphi and dZ/dphi of a field line in right-handed (R, phi, Z).
    call emit('fieldline_slope', [r, br, bphi, bz], &
        [character(8) :: 'r', 'br', 'bphi', 'bz'], [r*br/bphi, r*bz/bphi])
    ! Cartesian (Bx, By) to cylindrical (BR, Bphi) at toroidal angle phi.
    call emit('cylindrical_components', [phi, bx, by], &
        [character(8) :: 'phi', 'bx', 'by'], &
        [bx*cos(phi) + by*sin(phi), by*cos(phi) - bx*sin(phi)])
    ! Greene's residue (2 - tr J)/4 and det J of a 2x2 return-map Jacobian.
    call emit('map_invariants', [a, b, c, d], [character(8) :: 'a', 'b', 'c', 'd'], &
        [(2 - a - d)/4, a*d - b*c])
    ! Helical phase m theta - n phi of a section point (x, y) about the axis.
    call emit('helical_phase', [x, y, m, n, phi], &
        [character(8) :: 'x', 'y', 'm', 'n', 'phi'], [m*atan2(y, x) - n*phi])
contains
    subroutine emit(name, v, names, values)
        character(len=*), intent(in) :: name
        type(expr_t), intent(in) :: v(:), values(:)
        character(len=*), intent(in) :: names(:)
        type(kernel_spec_t) :: spec
        type(expr_t), allocatable :: roots(:)
        character(len=:), allocatable :: code
        character(len=32) :: ref, shape_value, shape_jac
        integer :: j, k, nv, nf, q, u, last
        nv = size(v); nf = size(values)
        allocate (spec%args(nv), roots(nf*(nv + 1)))
        do j = 1, nv
            spec%args(j) = str(trim(names(j)))
        end do
        spec%outputs = [str('value'), str('jac')]
        write (shape_value, '(A,I0,A)') '(', nf, ')'
        write (shape_jac, '(A,I0,A,I0,A)') '(', nf, ',', nv, ')'
        spec%output_shapes = [str(trim(shape_value)), str(trim(shape_jac))]
        allocate (spec%output_references(size(roots)))
        do j = 1, nf
            roots(j) = values(j)
            write (ref, '(A,I0,A)') 'value(', j, ')'
            spec%output_references(j) = str(trim(ref))
        end do
        q = nf
        do k = 1, nv
            do j = 1, nf
                q = q + 1
                roots(q) = diff(values(j), v(k))
                write (ref, '(A,I0,A,I0,A)') 'jac(', j, ',', k, ')'
                spec%output_references(q) = str(trim(ref))
            end do
        end do
        spec%name = str(name)
        spec%module_name = str('tiago_'//name)
        spec%generator = str('gen/gen_islands.f90')
        spec%generator_revision = &
            str('fortsym dbaa8e319ac5343e8c12a7938ca6a233d672a682')
        spec%regenerate_command = &
            str('build-codegen/gen_islands (see gen/README.md)')
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
