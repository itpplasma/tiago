program gen_ece_transfer
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
                            operator(*), operator(/), operator(**), exp
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: v(5), f, p, loss, term
    character(16) :: names(5)
    integer :: i

    call arena%init()
    names(1) = 'tau'
    call variables(1)
    f = 1 - exp(-v(1))
    call emit('ece_loss', 1, f)
    ! The alternating series avoids subtracting nearly equal exponentials.
    term = v(1)
    loss = term
    do i = 2, 10
        term = -term*v(1)/i
        loss = loss + term
    end do
    call emit('ece_loss_small', 1, loss)

    names(:3) = [character(16) :: 'incoming', 'te', 'loss']
    call variables(3)
    f = v(1) + v(3)*(v(2) - v(1))
    call emit('ece_thermal_cell', 3, f)

    names(1) = 'rho'
    call variables(1)
    call emit('ece_reflection_loss', 1, 1 - v(1))

    names = [character(16) :: 'near_emission', 'far_emission', 'loss', 'wall_loss', &
             'wall_te']
    call variables(5)
    p = 1 - v(3)
    ! Reciprocal ray between two walls with common reflectivity/temperature.
    ! Factor the denominator to retain accuracy for small loss and 1-rho.
    f = (v(1) + (1 - v(4))*p*v(2))/ &
        ((v(4) + (1 - v(4))*v(3))*(1 + (1 - v(4))*p)) + &
        v(4)*p*v(5)/(v(4) + (1 - v(4))*v(3))
    call emit('ece_two_stream_wall', 5, f)

    names = [character(16) :: 'temperature', 'receiver_te', 'gain', 'bandwidth', &
        'duration']
    call variables(5)
    f = (v(3)*(v(1) + v(2)))**2/(v(4)*v(5))
    call emit('ece_radiometer_noise', 5, f)
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
        spec%generator = str('gen/gen_ece_transfer.f90')
        spec%generator_revision = &
            str('fortsym 21d157b394e775a999a8d9ed0be409eac1eb9574')
        spec%regenerate_command = &
            str('build-codegen/gen_ece_transfer (see gen/README.md)')
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
