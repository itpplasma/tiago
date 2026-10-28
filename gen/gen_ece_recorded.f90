program gen_ece_recorded
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
        operator(*), operator(/), operator(**), sqrt, exp
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: v(6), f, radius
    character(16) :: names(6)

    call arena%init()
    names(1) = 'h'
    call variables(1)
    call emit('ece_launch_scale', 1, sqrt(1 - 2*v(1)))

    names(:3) = [character(16) :: 'bx', 'by', 'bz']
    call variables(3)
    call emit('ece_field_norm', 3, sqrt(v(1)**2 + v(2)**2 + v(3)**2))

    names = [character(16) :: 'x0', 'y0', 'z0', 'x1', 'y1', 'z1']
    call variables(6)
    f = sqrt((v(4) - v(1))**2 + (v(5) - v(2))**2 + (v(6) - v(3))**2)
    call emit('ece_segment_length', 6, f)

    names = [character(16) :: 'nx', 'ny', 'nz', 'bx', 'by', 'bz']
    call variables(6)
    f = (v(1)*v(4) + v(2)*v(5) + v(3)*v(6))**2/ &
        ((v(1)**2 + v(2)**2 + v(3)**2)*(v(4)**2 + v(5)**2 + v(6)**2))
    call emit('ece_parallel_fraction', 6, f)
    ! Maxwellian resonant-shell quadrature for an independent tenuous-X
    ! weakly relativistic finite-angle reference (Bornatici et al., 1983).
    names(:3) = [character(16) :: 'a', 'z', 'mu']
    call variables(3)
    radius = sqrt(v(1)**2 - v(2))
    f = exp(v(2) - 2*v(1)**2 - 2*v(1)*radius*v(3))
    call emit('ece_doppler_x1', 3, radius**3*(1 - v(3)**2)*f)
    call emit('ece_doppler_x2', 3, radius**5*(1 - v(3)**2)**2*f)
contains
    subroutine variables(n)
        integer, intent(in) :: n
        integer :: j
        do j = 1, n
            v(j) = sym(arena, trim(names(j)))
        end do
    end subroutine

    subroutine emit(name, n, value)
        character(*), intent(in) :: name
        integer, intent(in) :: n
        type(expr_t), intent(in) :: value
        type(expr_t) :: roots(n + 1)
        type(kernel_spec_t) :: spec
        character(:), allocatable :: code
        character(32) :: ref
        integer :: j, unit, last
        roots(1) = value
        allocate (spec%args(n), spec%output_references(n + 1))
        do j = 1, n
            spec%args(j) = str(trim(names(j)))
            roots(j + 1) = diff(value, v(j))
            write (ref, '(A,I0,A)') 'row(', j, ')'
            spec%output_references(j + 1) = str(trim(ref))
        end do
        spec%output_references(1) = str('value')
        spec%outputs = [str('value'), str('row')]
        write (ref, '(A,I0,A)') '(', n, ')'
        spec%output_shapes = [str(''), str(trim(ref))]
        spec%name = str(name)
        spec%module_name = str('tiago_'//name)
        spec%generator = str('gen/gen_ece_recorded.f90')
        spec%generator_revision = &
            str('fortsym 21d157b394e775a999a8d9ed0be409eac1eb9574')
        spec%regenerate_command = &
            str('build-codegen/gen_ece_recorded (see gen/README.md)')
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
