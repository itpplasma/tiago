program test_chease_interior
    !! Cross-code check of the interior: CHEASE re-solves the fixed-boundary
    !! equilibrium on the last closed surface, p' and F F' of the FreeGS
    !! G-EQDSK (tests/chease/chease_namelist), and TIAGO compares the plasma
    !! current, the magnetic axis, psi_N, rho_tor and the plasma part of the
    !! flux-loop and pickup signals of both files.
    !! Usage: test_chease_interior <chease> <namelist> <geqdsk> <machine> <workdir>
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_geqdsk_plasma, only: geqdsk_plasma_t, read_geqdsk_plasma
    use tiago_axisym_machine, only: axisym_machine_t, read_axisym_machine
    implicit none

    type(geqdsk_plasma_t) :: freegs, chease
    type(axisym_machine_t) :: machine
    character(len=4096) :: exe, namelist, gfile, mfile, work
    character(len=:), allocatable :: message
    real(dp), allocatable :: rf(:), zf(:), wf(:), rc(:), zc(:), wc(:), sf(:), sc(:)
    real(dp) :: r, z, pf, pc, a, b, scale, dpsi, drho
    integer :: ierr, status, failures, i, j, n

    failures = 0
    call get_command_argument(1, exe)
    call get_command_argument(2, namelist)
    call get_command_argument(3, gfile)
    call get_command_argument(4, mfile)
    call get_command_argument(5, work)
    call execute_command_line('mkdir -p '//trim(work)//' && cp '//trim(gfile)//' '//trim(work) &
        //'/EXPEQ && cp '//trim(namelist)//' '//trim(work)//'/chease_namelist && cd ' &
        //trim(work)//' && '//trim(exe)//' > chease.out 2>&1', exitstat=status)
    if (status /= 0) error stop 'CHEASE failed'
    call read_geqdsk_plasma(trim(gfile), freegs, ierr, message)
    if (ierr /= 0) error stop message
    call read_geqdsk_plasma(trim(work)//'/EQDSK_COCOS_07.OUT', chease, ierr, message)
    if (ierr /= 0) error stop message
    call read_axisym_machine(trim(mfile), machine, ierr, message)
    if (ierr /= 0) error stop message

    call freegs%filaments(rf, zf, wf)
    call chease%filaments(rc, zc, wc)
    print '(A,2ES16.8)', 'plasma current [A]: ', sum(wf), chease%current
    print '(A,4F12.6)', 'axis R, Z [m]:      ', freegs%r_axis, freegs%z_axis, &
        chease%r_axis, chease%z_axis
    call close_to(chease%current, sum(wf), 5.0e-4_dp * abs(sum(wf)), 'CHEASE plasma current')
    call close_to(chease%r_axis, freegs%r_axis, 1.0e-3_dp, 'axis R')
    call close_to(chease%z_axis, freegs%z_axis, 1.0e-3_dp, 'axis Z')

    ! psi_N on a 9 x 9 grid of points inside the plasma.
    n = 0
    dpsi = 0
    drho = 0
    do j = 1, 9
        do i = 1, 9
            r = 0.9_dp + 0.7_dp * (i - 1) / 8
            z = -0.5_dp + 1.0_dp * (j - 1) / 8
            call freegs%label(r, z, pf, a, b)
            if (pf >= 0.95_dp) cycle
            call chease%label(r, z, pc, a, b)
            call close_to(pc, pf, 2.0e-3_dp, 'psi_N')
            call close_to(chease%rho_tor(pf), freegs%rho_tor(pf), 5.0e-3_dp, 'rho_tor')
            dpsi = max(dpsi, abs(pc - pf))
            drho = max(drho, abs(chease%rho_tor(pf) - freegs%rho_tor(pf)))
            n = n + 1
        end do
    end do
    print '(A,I0,A,2ES10.3)', 'psi_N at ', n, ' points, largest psi_N and rho_tor differences: ', &
        dpsi, drho

    sf = machine%sensors%signals(rf, zf, wf)
    sc = machine%sensors%signals(rc, zc, wc)
    n = size(machine%sensors%flux)
    scale = maxval(abs(sf(:n)))
    print '(A,ES10.3)', 'flux loops, largest difference / largest signal: ', &
        maxval(abs(sc(:n) - sf(:n))) / scale
    do i = 1, n
        call close_to(sc(i), sf(i), 2.0e-3_dp * scale, 'plasma part of a flux loop')
    end do
    scale = maxval(abs(sf(n + 1:size(sf) - 1)))
    print '(A,ES10.3)', 'pickups, largest difference / largest signal:    ', &
        maxval(abs(sc(n + 1:size(sf) - 1) - sf(n + 1:size(sf) - 1))) / scale
    do i = n + 1, size(sf) - 1
        call close_to(sc(i), sf(i), 2.0e-3_dp * scale, 'plasma part of a pickup')
    end do
    if (failures > 0) then
        write(error_unit, '(I0,A)') failures, ' CHEASE interior checks failed'
        error stop 1
    end if
    print '(A)', 'CHEASE interior: all checks passed'

contains

    subroutine close_to(x, y, tolerance, what)
        real(dp), intent(in) :: x, y, tolerance
        character(*), intent(in) :: what
        if (abs(x - y) > tolerance) then
            write(error_unit, '(A,": ",ES23.15," vs ",ES23.15)') what, x, y
            failures = failures + 1
        end if
    end subroutine close_to

end program test_chease_interior
