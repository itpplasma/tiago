program test_geqdsk_cells
    !! Plasma current and plasma part of flux loops and pickups of analytic
    !! Cerfon-Freidberg equilibria (tests/analytic/generate.py) from the cell
    !! quadrature of tiago_geqdsk_plasma, against Gauss quadrature of the exact
    !! current density. The grid-node sum misses these by about 1e-3, because
    !! the boundary cuts grid cells.
    !! Usage: test_geqdsk_cells <geqdsk> <machine> <reference.csv> <tolerance>
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_geqdsk_plasma, only: geqdsk_plasma_t, read_geqdsk_plasma
    use tiago_axisym_machine, only: axisym_machine_t, read_axisym_machine
    implicit none
    type(geqdsk_plasma_t) :: eq
    type(axisym_machine_t) :: machine
    character(len=1024) :: gpath, mpath, cpath, arg, line
    character(len=:), allocatable :: message
    real(dp), allocatable :: rp(:), zp(:), wp(:), s(:), ref(:)
    real(dp) :: tol, current, value, scale
    integer :: ierr, unit, ios, k, comma
    call get_command_argument(1, gpath)
    call get_command_argument(2, mpath)
    call get_command_argument(3, cpath)
    call get_command_argument(4, arg)
    read(arg, *) tol
    call read_geqdsk_plasma(trim(gpath), eq, ierr, message)
    if (ierr /= 0) error stop message
    call read_axisym_machine(trim(mpath), machine, ierr, message)
    if (ierr /= 0) error stop message
    eq%cell_quadrature = .true.
    call eq%filaments(rp, zp, wp)
    s = machine%sensors%signals(rp, zp, wp)
    allocate(ref(0))
    open(newunit=unit, file=trim(cpath), status='old', action='read')
    read(unit, '(A)') line
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        comma = index(line, ',')
        read(line(comma + 1:), *) value
        if (line(:comma - 1) == 'plasma_current') then
            current = value
        else
            ref = [ref, value]
        end if
    end do
    close(unit)
    call check('plasma current', sum(wp), current, tol * abs(current))
    scale = maxval(abs(ref(:4)))
    do k = 1, 4
        call check('flux loop', s(k), ref(k), tol * scale)
    end do
    scale = maxval(abs(ref(5:)))
    do k = 5, size(ref)
        call check('pickup', s(k), ref(k), tol * scale)
    end do
    print '(A)', 'tiago_geqdsk_cells: OK'
contains
    subroutine check(what, got, want, atol)
        character(len=*), intent(in) :: what
        real(dp), intent(in) :: got, want, atol
        if (abs(got - want) > atol) then
            print '(A,3ES14.6)', what//': got, want, difference ', got, want, got - want
            error stop 1
        end if
    end subroutine check
end program test_geqdsk_cells
