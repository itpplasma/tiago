module tiago_bprobes
    !! DIAGNO magnetic probe files (bprobes_file): a count, then per probe
    !!     x y z theta_inc phi_inc eff_area        (or R phi[rad] z ... with rphiz)
    !! with the probe normal (sin phi_inc cos theta_inc, sin phi_inc sin theta_inc,
    !! cos phi_inc) and the angles in degrees. Probes are labelled PROBE_0001, ...
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    use tiago_diagnostic_types, only: bprobe_t
    implicit none
    private

    real(dp), parameter :: deg = acos(-1.0_dp) / 180.0_dp

    public :: read_bprobe_file

contains

    subroutine read_bprobe_file(path, probes, ierr, message, rphiz)
        character(len=*), intent(in) :: path
        type(bprobe_t), allocatable, intent(out) :: probes(:)
        integer(i32), intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        logical, intent(in), optional :: rphiz

        integer :: unit, ios, i, count
        real(dp) :: row(6), th, ph
        character(len=16) :: label

        ierr = 0_i32
        message = ''
        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1_i32
            message = 'unable to open B-probe file: '//trim(path)
            return
        end if
        read(unit, *, iostat=ios) count
        if (ios /= 0 .or. count <= 0) then
            ierr = 2_i32
            message = 'invalid probe count in '//trim(path)
            close(unit)
            return
        end if
        allocate(probes(count))
        do i = 1, count
            read(unit, *, iostat=ios) row
            if (ios /= 0) then
                ierr = 5_i32
                write(label, '(I0)') i
                message = 'missing or malformed row for probe '//trim(label)//' in '//trim(path)
                close(unit)
                return
            end if
            if (present(rphiz)) then
                if (rphiz) row(1:2) = [row(1) * cos(row(2)), row(1) * sin(row(2))]
            end if
            th = row(4) * deg
            ph = row(5) * deg
            write(label, '(A,I4.4)') 'PROBE_', i
            probes(i)%label = trim(label)
            probes(i)%position = row(1:3)
            probes(i)%normal = [sin(ph) * cos(th), sin(ph) * sin(th), cos(ph)]
            probes(i)%eff_area = row(6)
        end do
        close(unit)
    end subroutine read_bprobe_file
end module tiago_bprobes
