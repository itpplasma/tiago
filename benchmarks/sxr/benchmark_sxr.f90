program benchmark_sxr
    !! Soft X-ray chords on a VMEC equilibrium for the comparison with
    !! STELLOPT's chisq_sxr (see README.md). One geometry and one set of
    !! profiles serve both codes:
    !!   benchmark_sxr stellopt <out>          STELLOPT &OPTIMUM lines for the
    !!                                         profiles and chords
    !!   benchmark_sxr tiago <wout.nc> <csv>   TIAGO signals for npanel = 50 ...
    !!                                         1600 panels of 4 Gauss points
    !! Profiles (s the normalized toroidal flux, as STELLOPT power_series):
    !! ne = 1e20 (1 - s) m^-3, Te = 2 keV (1 - s), Zeff = 1.5, so that the
    !! STELLOPT emissivity Zeff ne^2 Te^(1/2) and the unfiltered TIAGO
    !! bremsstrahlung vanish smoothly at the boundary; both are zero outside.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_vmec_state, only: vmec_state_t
    use tiago_real_space_state, only: located_t
    use tiago_state_sampler, only: sample_state
    use tiago_profile_diagnostic, only: plasma_sample_t, Q_S, Q_NE, Q_TE, Q_ZEFF
    use tiago_sight_geometry, only: chord_t
    use tiago_soft_xray, only: soft_xray_t, pinhole_camera_t, camera_chords, &
        SXR_BREMSSTRAHLUNG_MODEL
    implicit none

    real(dp), parameter :: n0 = 1.0e20_dp, t0 = 2.0e3_dp, zeff0 = 1.5_dp
    integer, parameter :: panels(6) = [50, 100, 200, 400, 800, 1600]
    character(len=1024) :: mode, a, b
    type(chord_t), allocatable :: chords(:)

    call get_command_argument(1, mode)
    call get_command_argument(2, a)
    call get_command_argument(3, b)
    chords = geometry(panels(1))
    select case (trim(mode))
    case ('stellopt')
        call write_stellopt(trim(a))
    case ('tiago')
        call evaluate(trim(a), trim(b))
    case default
        write(error_unit, '(A)') 'usage: benchmark_sxr stellopt <out> | tiago <wout.nc> <csv>'
        error stop 1
    end select

contains

    function geometry(npanel) result(all)
        !! Three cameras, 12 channels each: above the plasma and on the
        !! outboard side in the symmetry plane phi = 0, and an oblique fan at
        !! phi = 18 degrees (between symmetry planes of the three-period LI383).
        integer, intent(in) :: npanel
        type(chord_t), allocatable :: all(:), part(:)
        type(pinhole_camera_t) :: cam
        real(dp), allocatable :: g(:)
        real(dp) :: phi, er(3), ez(3), ep(3)
        integer :: k, c

        allocate(all(0))
        do c = 1, 3
            phi = merge(0.0_dp, acos(-1.0_dp) / 10.0_dp, c < 3)
            er = [cos(phi), sin(phi), 0.0_dp]
            ep = [-sin(phi), cos(phi), 0.0_dp]
            ez = [0.0_dp, 0.0_dp, 1.0_dp]
            allocate(cam%detector(3, 12))
            select case (c)
            case (1)
                cam%name = 'top'
                cam%aperture = 1.45_dp * er + 0.75_dp * ez
                do k = 1, 12
                    cam%detector(:, k) = cam%aperture + 0.05_dp * ez &
                        + 0.004_dp * (k - 6.5_dp) * er
                end do
            case (2)
                cam%name = 'outboard'
                cam%aperture = 2.25_dp * er
                do k = 1, 12
                    cam%detector(:, k) = cam%aperture + 0.05_dp * er &
                        + 0.0045_dp * (k - 6.5_dp) * ez
                end do
            case default
                cam%name = 'oblique'
                cam%aperture = 2.2_dp * er + 0.45_dp * ez + 0.1_dp * ep
                do k = 1, 12
                    cam%detector(:, k) = cam%aperture + 0.035_dp * er + 0.035_dp * ez &
                        + 0.0022_dp * (k - 3) * (er - ez) / sqrt(2.0_dp)
                end do
            end select
            cam%length = 1.6_dp
            call camera_chords(cam, npanel, part, g)
            all = [all, part]
            deallocate(cam%detector)
        end do
    end function geometry

    subroutine write_stellopt(path)
        !! Chords as STELLOPT (R, phi, Z) end points; sigma = 1 selects them.
        character(len=*), intent(in) :: path
        real(dp) :: x0(3), x1(3)
        integer :: u, c

        open(newunit=u, file=path, status='replace', action='write')
        write(u, '(A)') "  NE_TYPE = 'power_series'"
        write(u, '(A,ES23.15,A,ES23.15)') '  NE_OPT(0:1) = ', n0, ', ', -n0
        write(u, '(A)') "  TE_TYPE = 'power_series'"
        write(u, '(A,ES23.15,A,ES23.15)') '  TE_OPT(0:1) = ', t0, ', ', -t0
        write(u, '(A)') "  ZEFF_TYPE = 'power_series'"
        write(u, '(A,ES23.15)') '  ZEFF_OPT(0) = ', zeff0
        do c = 1, size(chords)
            x0 = chords(c)%origin
            x1 = x0 + chords(c)%length * chords(c)%direction / norm2(chords(c)%direction)
            write(u, '(A,I0,A,ES23.15,A,I0,A)') '  TARGET_SXR(', c, ') = ', 1.0_dp, &
                '  SIGMA_SXR(', c, ') = 1.0'
            write(u, '(3(A,I0,A,ES23.15))') '  R0_SXR(', c, ') = ', hypot(x0(1), x0(2)), &
                '  PHI0_SXR(', c, ') = ', atan2(x0(2), x0(1)), '  Z0_SXR(', c, ') = ', x0(3)
            write(u, '(3(A,I0,A,ES23.15))') '  R1_SXR(', c, ') = ', hypot(x1(1), x1(2)), &
                '  PHI1_SXR(', c, ') = ', atan2(x1(2), x1(1)), '  Z1_SXR(', c, ') = ', x1(3)
        end do
        close(u)
    end subroutine write_stellopt

    subroutine evaluate(wout, csv)
        !! STELLOPT's integrand without its constant: Zeff ne^2 Te^(1/2)
        !! [1e40 m^-6 eV^(1/2)], i.e. the TIAGO bremsstrahlung divided by 1.69e-38.
        character(len=*), intent(in) :: wout, csv
        type(vmec_state_t) :: vmec
        type(soft_xray_t) :: sxr
        type(plasma_sample_t) :: smp
        type(located_t), allocatable :: loc(:)
        real(dp), allocatable :: x(:, :), dx(:, :, :), y(:, :)
        real(dp) :: t_start, t_end
        integer :: ierr, k, p, u, c
        character(len=:), allocatable :: message

        call vmec%load(wout, ierr, message)
        if (ierr /= 0) then
            write(error_unit, '(A)') message
            error stop 1
        end if
        allocate(y(size(chords), size(panels)))
        do p = 1, size(panels)
            call sxr%setup(geometry(panels(p)))
            sxr%model = SXR_BREMSSTRAHLUNG_MODEL
            call cpu_time(t_start)
            call sxr%sample_points([real(dp) ::], x, dx)
            call sample_state(vmec, x, smp, loc)
            do k = 1, smp%n
                smp%q(Q_ZEFF, k) = zeff0
                if (.not. smp%inside(k)) cycle
                smp%q(Q_NE, k) = n0 * (1.0_dp - smp%q(Q_S, k))
                smp%q(Q_TE, k) = t0 * (1.0_dp - smp%q(Q_S, k))
            end do
            call sxr%predict(smp, [real(dp) ::], y(:, p))
            call cpu_time(t_end)
            y(:, p) = y(:, p) / 1.69e-38_dp / 1.0e40_dp
            print '(A,I5,A,F8.3,A)', 'npanel ', panels(p), ': ', t_end - t_start, ' s'
        end do
        open(newunit=u, file=csv, status='replace', action='write')
        write(u, '(A,*(A,I0))') 'channel,label', (',npanel_', panels(p), p = 1, size(panels))
        do c = 1, size(chords)
            write(u, '(I0,A,A,*(:,",",ES23.15))') c, ',', chords(c)%label, y(c, :)
        end do
        close(u)
    end subroutine evaluate
end program benchmark_sxr
