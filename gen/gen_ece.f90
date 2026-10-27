program gen_ece
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, real_expr, pi_expr, operator(+), &
                          operator(-), operator(*), operator(/), operator(**), sqrt
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: pos(3), par(8), field(4)
    type(expr_t) :: v(8), h, n2, b2, nb, xx, yy, mu2, disc, den
    type(expr_t), allocatable :: roots(:)
    type(kernel_spec_t) :: spec
    character(len=4), parameter :: names(8) = &
                        ['nx  ', 'ny  ', 'nz  ', 'bx  ', 'by  ', 'bz  ', 'ne  ', 'freq']
    character(len=40) :: ref
    character(len=:), allocatable :: code, name
    integer :: i, j, k, mode, unit

    call arena%init()
    do i = 1, 8
        v(i) = sym(arena, trim(names(i)))
    end do
    n2 = v(1)**2 + v(2)**2 + v(3)**2
    b2 = v(4)**2 + v(5)**2 + v(6)**2
    nb = v(1)*v(4) + v(2)*v(5) + v(3)*v(6)
    xx = real_expr(arena, 1.602176634e-19_dp**2/ &
                   (8.8541878128e-12_dp*9.1093837015e-31_dp))*v(7)/ &
         (2*pi_expr(arena)*v(8))**2
    yy = (real_expr(arena, 1.602176634e-19_dp/9.1093837015e-31_dp)/ &
          (2*pi_expr(arena)*v(8)))**2*b2
    mu2 = nb**2/(n2*b2)
    disc = sqrt(yy**2*(1 - mu2)**2 + 4*yy*(1 - xx)**2*mu2)
    do mode = 1, 2
        den = 2*(1 - xx) - yy*(1 - mu2)
        if (mode == 1) den = den + disc
        if (mode == 2) den = den - disc
        h = (n2 - 1 + 2*xx*(1 - xx)/den)/2
        allocate (roots(73))
        roots(1) = h
        do i = 1, 8
            roots(1 + i) = diff(h, v(i))
            do j = 1, 8
                roots(9 + (i - 1)*8 + j) = diff(diff(h, v(i)), v(j))
            end do
        end do
        call reset(8)
        do i = 1, 8
            spec%args(i) = str(trim(names(i)))
        end do
        spec%outputs = [str('h'), str('dh'), str('d2h')]
        spec%output_shapes = [str(''), str('(8)'), str('(8, 8)')]
        allocate (spec%output_references(73))
        spec%output_references(1) = str('h')
        do i = 1, 8
            write (ref, '(A,I0,A)') 'dh(', i, ')'
            spec%output_references(1 + i) = str(trim(ref))
            do j = 1, 8
                write (ref, '(A,I0,A,I0,A)') 'd2h(', i, ',', j, ')'
                spec%output_references(9 + (i - 1)*8 + j) = str(trim(ref))
            end do
        end do
        name = 'ece_hamiltonian_o'
        if (mode == 2) name = 'ece_hamiltonian_x'
        call emit(name)
        deallocate (roots)
    end do

    ! Unmagnetized limit avoids the removable 0/0 in the angular variables.
    h = (n2 - 1 + xx)/2
    allocate (roots(73))
    roots(1) = h
    do i = 1, 8
        roots(1 + i) = diff(h, v(i))
        do j = 1, 8
            roots(9 + (i - 1)*8 + j) = diff(diff(h, v(i)), v(j))
        end do
    end do
    call emit('ece_hamiltonian_unmagnetized')
    deallocate (roots)
    do i = 1, 3
        write (ref, '(A,I0)') 'x', i
        pos(i) = sym(arena, trim(ref))
    end do
    do i = 1, 8
        write (ref, '(A,I0)') 'p', i
        par(i) = sym(arena, trim(ref))
    end do
    field(1) = par(1)*(1 + pos(2)/10)
    field(2) = par(2)*(1 + pos(3)/20)
    field(3) = par(3)*(1 + pos(1)/30)
    field(4) = par(4) + par(5)*pos(1) + par(6)*pos(2) + par(7)*pos(3) &
               + par(8)*(pos(1)**2 + pos(2)**2)
    allocate (roots(4 + 12 + 36 + 32 + 96))
    call reset(11)
    do i = 1, 3
        write (ref, '(A,I0)') 'x', i
        spec%args(i) = str(trim(ref))
    end do
    do i = 1, 8
        write (ref, '(A,I0)') 'p', i
        spec%args(i + 3) = str(trim(ref))
    end do
    spec%outputs = [str('q'), str('dq'), str('d2q'), str('qp'), str('dqp')]
    spec%output_shapes = [str('(4)'), str('(3, 4)'), str('(3, 3, 4)'), &
                          str('(4, 8)'), str('(3, 4, 8)')]
    allocate (spec%output_references(size(roots)))
    k = 0
    do i = 1, 4
        k = k + 1
        roots(k) = field(i)
        write (ref, '(A,I0,A)') 'q(', i, ')'
        spec%output_references(k) = str(trim(ref))
        do j = 1, 3
            k = k + 1
            roots(k) = diff(field(i), pos(j))
            write (ref, '(A,I0,A,I0,A)') 'dq(', j, ',', i, ')'
            spec%output_references(k) = str(trim(ref))
            block
                integer :: l
                do l = 1, 3
                    k = k + 1
                    roots(k) = diff(diff(field(i), pos(j)), pos(l))
                    write (ref, '(A,I0,A,I0,A,I0,A)') 'd2q(', j, ',', l, ',', i, ')'
                    spec%output_references(k) = str(trim(ref))
                end do
            end block
        end do
        do j = 1, 8
            k = k + 1
            roots(k) = diff(field(i), par(j))
            write (ref, '(A,I0,A,I0,A)') 'qp(', i, ',', j, ')'
            spec%output_references(k) = str(trim(ref))
            block
                integer :: l
                do l = 1, 3
                    k = k + 1
                    roots(k) = diff(diff(field(i), par(j)), pos(l))
                    write (ref, '(A,I0,A,I0,A,I0,A)') 'dqp(', l, ',', i, ',', j, ')'
                    spec%output_references(k) = str(trim(ref))
                end do
            end block
        end do
    end do
    call emit('ece_test_medium')

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
        spec%generator = str('gen/gen_ece.f90')
        spec%generator_revision = &
            str('fortsym 21d157b394e775a999a8d9ed0be409eac1eb9574')
        spec%regenerate_command = str('build-codegen/gen_ece (see gen/README.md)')
        spec%pure_procedure = .true.
        code = chars(emit_kernel(roots, spec))
        n = len_trim(code)
        do while (n > 0)
            if (code(n:n) /= new_line('a')) exit
            n = n - 1
        end do
        open (newunit=unit, file='src/generated/'//kernel_name//'.f90', &
              status='replace')
        write (unit, '(A)') code(:n)
        close (unit)
    end subroutine
end program
