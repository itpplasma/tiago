program test_geqdsk_plasma
    !! G-EQDSK plasma source and labels against FreeGS on its TestTokamakSensor
    !! machine (tests/freegs/generate_test_tokamak.py): plasma current, the
    !! flux-loop, pickup and Rogowski signals of coils plus plasma, psi_N at
    !! points above and below the midplane of the up-down asymmetric plasma,
    !! rho_tor and the field B_R, B_phi, B_Z.
    !! Usage: test_geqdsk_plasma <geqdsk> <machine> <freegs.csv>
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_geqdsk_plasma, only: geqdsk_plasma_t, read_geqdsk_plasma
    use tiago_axisym_machine, only: axisym_machine_t, read_axisym_machine
    implicit none

    type(geqdsk_plasma_t) :: eq
    type(axisym_machine_t) :: machine
    character(len=1024) :: gpath, mpath, cpath
    character(len=:), allocatable :: message
    character(len=256) :: line
    character(len=32) :: kind, label
    real(dp), allocatable :: rp(:), zp(:), wp(:), rf(:), zf(:), wf(:), s(:), current(:)
    real(dp) :: r, z, value, psin, a, b, scale(3), field(3)
    integer :: ierr, unit, ios, failures, k, c, nflux, npick
    real(dp), allocatable :: ref(:)

    failures = 0
    call get_command_argument(1, gpath)
    call get_command_argument(2, mpath)
    call get_command_argument(3, cpath)
    call read_geqdsk_plasma(trim(gpath), eq, ierr, message)
    if (ierr /= 0) error stop message
    call read_axisym_machine(trim(mpath), machine, ierr, message)
    if (ierr /= 0) error stop message
    call check(eq%exponent == 0 .and. eq%sigma > 0, 'FreeGS psi per radian, same sign')

    call eq%filaments(rp, zp, wp)
    nflux = size(machine%sensors%flux)
    npick = size(machine%sensors%pickup)
    allocate(current(machine%ngroups), ref(machine%sensors%count()))
    open(newunit=unit, file=trim(cpath), status='old', action='read')
    read(unit, '(A)') line
    k = 0
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        call split(line, kind, label, r, z, value)
        select case (trim(kind))
        case ('current')
            if (trim(label) == 'plasma') then
                call close_to(sum(wp), value, 1.0e-5_dp * abs(value), 'plasma current')
            else
                k = k + 1
                current(k) = value
            end if
        case ('flux')
            ref(index_of(label, 'flux')) = value
        case ('pickup')
            ref(nflux + index_of(label, 'pickup')) = value
        case ('rogowski')
            ref(nflux + npick + 1) = value
        case ('psin')
            call eq%label(r, z, psin, a, b)
            call close_to(psin, value, 1.0e-6_dp, 'psi_N at '//trim(label))
        case ('br', 'bphi', 'bz')
            call eq%field(r, z, field)
            c = merge(1, merge(2, 3, trim(kind) == 'bphi'), trim(kind) == 'br')
            ! measured deviation below 1e-6 T for B_p ~ 0.1 T and B_phi ~ 2 T
            call close_to(field(c), value, 5.0e-6_dp, &
                trim(kind)//' at '//trim(label))
        case ('rhotor')
            call close_to(eq%rho_tor(r), value, 3.0e-4_dp, 'rho_tor at '//trim(label))
        end select
    end do
    close(unit)

    wf = [machine%turns * current(machine%group), wp]
    rf = [machine%r, rp]
    zf = [machine%z, zp]
    s = machine%sensors%signals(rf, zf, wf)
    scale = [maxval(abs(ref(:nflux))), maxval(abs(ref(nflux + 1:nflux + npick))), &
        abs(ref(size(ref)))]
    do k = 1, size(s)
        if (k <= nflux) then
            call close_to(s(k), ref(k), 2.0e-5_dp * scale(1), 'flux loop')
        else if (k <= nflux + npick) then
            call close_to(s(k), ref(k), 5.0e-5_dp * scale(2), 'pickup')
        else
            call close_to(s(k), ref(k), 1.0e-5_dp * scale(3), 'Rogowski')
        end if
    end do
    if (failures > 0) then
        write(error_unit, '(I0,A)') failures, ' G-EQDSK plasma checks failed'
        error stop 1
    end if
    print '(A)', 'G-EQDSK plasma: all checks passed'

contains

    integer function index_of(name, what) result(i)
        character(*), intent(in) :: name, what
        if (what == 'flux') then
            do i = 1, size(machine%sensors%flux)
                if (machine%sensors%flux(i)%label == trim(name)) return
            end do
        else
            do i = 1, size(machine%sensors%pickup)
                if (machine%sensors%pickup(i)%label == trim(name)) return
            end do
        end if
        error stop 'sensor of the reference not in the machine file'
    end function index_of

    subroutine split(text, kind, label, r, z, value)
        character(*), intent(in) :: text
        character(*), intent(out) :: kind, label
        real(dp), intent(out) :: r, z, value
        character(len=len(text)) :: t
        integer :: i
        t = text
        do i = 1, len_trim(t)
            if (t(i:i) == ',') t(i:i) = ' '
        end do
        read(t, *) kind, label, r, z, value
    end subroutine split

    subroutine close_to(x, y, tolerance, what)
        real(dp), intent(in) :: x, y, tolerance
        character(*), intent(in) :: what
        if (abs(x - y) > tolerance) then
            write(error_unit, '(A,": ",ES23.15," vs ",ES23.15)') what, x, y
            failures = failures + 1
        end if
    end subroutine close_to

    subroutine check(condition, what)
        logical, intent(in) :: condition
        character(*), intent(in) :: what
        if (.not. condition) then
            write(error_unit, '(A)') what
            failures = failures + 1
        end if
    end subroutine check

end program test_geqdsk_plasma
