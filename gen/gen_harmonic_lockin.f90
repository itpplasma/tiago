program gen_harmonic_lockin
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(-), &
        operator(*), operator(/), sin, cos
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: v(4), f(3)
    character(len=16) :: names(4)
    integer :: i

    call arena%init()
    names(:4) = [character(len=16) :: "sample", "offset", "gain", "theta"]
    call variables(4)
    f(1) = (v(1) - v(2))/v(3)
    f(2) = cos(2*v(4))
    f(3) = sin(v(4))
    call emit("harmonic_mixer", 4, f)
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
        character(len=:), allocatable :: code
        character(len=32) :: ref, shape_value, shape_jac
        integer :: j, k, m, r, u, last
        m = size(values)
        allocate (spec%args(n), roots(m*(n + 1)))
        do j = 1, n
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
        spec%name = str(name)
        spec%module_name = str('tiago_'//name)
        spec%generator = str('gen/gen_harmonic_lockin.f90')
        spec%generator_revision = &
            str('fortsym ae635c80b42c922f5f2028c5cec67859bb8d3929')
        spec%regenerate_command = &
            str('build-codegen/gen_harmonic_lockin (see gen/README.md)')
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
