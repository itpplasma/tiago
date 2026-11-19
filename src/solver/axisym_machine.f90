module tiago_axisym_machine
    !! Text description of an axisymmetric machine, one item per line:
    !!     coil <label> <group> <R> <Z> <dR> <dZ> <turns>     (dR = dZ = 0: filament)
    !!     flux <label> <R> <Z>
    !!     pickup <label> <R> <Z> <angle> <length> <area>
    !!     rogowski <label> <R1> <Z1> <R2> <Z2> ...
    !! in m, degrees and m^2. The pickup axis is (cos angle, sin angle) in the
    !! (R, Z) plane and its length runs along the axis. Lines starting with #
    !! are comments.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_axisym_magnetics, only: axisym_sensors_t, axisym_flux_loop_t, &
        axisym_pickup_t, axisym_rogowski_t, axisym_rectangle
    implicit none
    private

    real(dp), parameter :: deg = acos(-1.0_dp) / 180
    integer, parameter :: coil_points = 4

    !> Coil filaments with signed turns per unit group current.
    type, public :: axisym_machine_t
        type(axisym_sensors_t) :: sensors
        real(dp), allocatable :: r(:), z(:), turns(:)
        integer, allocatable :: group(:)
        integer :: ngroups = 0
    end type axisym_machine_t

    public :: read_axisym_machine

contains

    subroutine read_axisym_machine(path, machine, ierr, message)
        character(*), intent(in) :: path
        type(axisym_machine_t), intent(out) :: machine
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        character(len=4096) :: line
        character(len=64) :: kind, label
        real(dp) :: v(6), buffer(512)
        real(dp), allocatable :: rf(:), zf(:), wf(:)
        integer :: unit, ios, g, n, lineno
        type(axisym_flux_loop_t) :: loop
        type(axisym_pickup_t) :: pickup
        type(axisym_rogowski_t) :: rogowski
        ierr = 0
        message = ''
        allocate(machine%r(0), machine%z(0), machine%turns(0), machine%group(0))
        allocate(machine%sensors%flux(0), machine%sensors%pickup(0), machine%sensors%rogowski(0))
        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1
            message = 'cannot open machine file '//trim(path)
            return
        end if
        lineno = 0
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            lineno = lineno + 1
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            read(line, *, iostat=ios) kind, label
            if (ios /= 0) call fail()
            select case (trim(kind))
            case ('coil')
                read(line, *, iostat=ios) kind, label, g, v(1:5)
                if (ios /= 0 .or. g < 1) call fail()
                if (v(3) > 0 .and. v(4) > 0) then
                    call axisym_rectangle(v(1), v(2), v(3), v(4), v(5), coil_points, rf, zf, wf)
                else
                    rf = [v(1)]
                    zf = [v(2)]
                    wf = [v(5)]
                end if
                machine%r = [machine%r, rf]
                machine%z = [machine%z, zf]
                machine%turns = [machine%turns, wf]
                machine%group = [machine%group, spread(g, 1, size(rf))]
                machine%ngroups = max(machine%ngroups, g)
            case ('flux')
                read(line, *, iostat=ios) kind, label, v(1:2)
                if (ios /= 0) call fail()
                loop%label = trim(label)
                loop%r = v(1)
                loop%z = v(2)
                machine%sensors%flux = [machine%sensors%flux, loop]
            case ('pickup')
                read(line, *, iostat=ios) kind, label, v(1:5)
                if (ios /= 0) call fail()
                pickup%label = trim(label)
                pickup%r = v(1)
                pickup%z = v(2)
                pickup%nr = cos(v(3) * deg)
                pickup%nz = sin(v(3) * deg)
                pickup%tr = pickup%nr
                pickup%tz = pickup%nz
                pickup%length = v(4)
                pickup%area = v(5)
                pickup%points = merge(4, 1, v(4) > 0)
                machine%sensors%pickup = [machine%sensors%pickup, pickup]
            case ('rogowski')
                do n = 2 * size(buffer) / 2, 2, -2
                    read(line, *, iostat=ios) kind, label, buffer(:n)
                    if (ios == 0) exit
                end do
                if (ios /= 0 .or. n < 6) call fail()
                rogowski%label = trim(label)
                rogowski%r = buffer(1:n:2)
                rogowski%z = buffer(2:n:2)
                machine%sensors%rogowski = [machine%sensors%rogowski, rogowski]
            case default
                call fail()
            end select
            if (ierr /= 0) exit
        end do
        close(unit)
    contains
        subroutine fail()
            character(len=16) :: number
            write(number, '(I0)') lineno
            ierr = 2
            message = 'malformed line '//trim(number)//' in '//trim(path)
        end subroutine fail
    end subroutine read_axisym_machine

end module tiago_axisym_machine
