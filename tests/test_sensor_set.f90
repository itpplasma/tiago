program test_sensor_set
    !! Sensor-set ingestion (tests/data/sensors/synthetic.sensors), signals and
    !! Jacobians of saddle, flux, diamagnetic and compensation loops, probes and
    !! linear combinations.
    !! Usage: test_sensor_set <manifest> [<wout.nc>]
    !! Oracles: Maxwell's mutual inductance of coaxial circles (elliptic
    !! integrals), the on-axis field of a circular coil, zero vacuum flux through
    !! a poloidal-plane loop, linearity in the coil currents, and with a VMEC
    !! plasma the exact linearity of the sheet field in its coefficients.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32, error_unit
    use tiago_sensor_set, only: sensor_set_t, read_sensor_set, kind_saddle_coil, &
        kind_diamagnetic_loop, kind_compensation_loop, kind_combination, kind_bprobe
    use tiago_sensor_signals, only: sensor_set_signals, sensor_set_coil_jacobian, &
        sensor_set_plasma_jacobian, sensor_set_linked_flux_row
    use tiago_vacuum_forward, only: vacuum_solver_t, quadrature_rule_t
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 4.0e-7_dp * pi
    real(dp), parameter :: current = 1000.0_dp, linked = 0.37_dp
    integer, parameter :: ncoil = 4096, nloop = 1024
    character(len=1024) :: manifest, dir, wout
    type(sensor_set_t) :: set
    type(vacuum_solver_t) :: solver
    type(quadrature_rule_t) :: rule
    integer(i32) :: ierr
    character(len=:), allocatable :: message
    real(dp), allocatable :: y(:), jac(:, :), row(:), ybar(:), rawbar(:), raw(:), yc(:)
    integer :: failures, k, status

    failures = 0
    call get_command_argument(1, manifest)
    call get_environment_variable('TIAGO_TEST_SENSOR_GEOMETRY', dir, status=status)
    if (status /= 0) error stop 'set TIAGO_TEST_SENSOR_GEOMETRY'
    call write_geometry(trim(dir))
    call read_sensor_set(trim(manifest), set, ierr, message)
    if (ierr /= 0) then
        write(error_unit, '(A)') message
        error stop 1
    end if
    call check(set%n_signals() == 6, 'six signals')
    call check(set%kinds(set%index_of('SADDLE')) == kind_saddle_coil, 'saddle kind')
    call check(set%kinds(set%index_of('DIA')) == kind_diamagnetic_loop, 'diamagnetic kind')
    call check(set%kinds(set%index_of('COMP')) == kind_compensation_loop, 'compensation kind')
    call check(set%kinds(set%index_of('PROBE_0001')) == kind_bprobe, 'probe kind')
    call check(set%kinds(6) == kind_combination, 'combination kind')

    call solver%init(trim(dir)//'/coil.coils')
    rule%samples_per_segment = 4
    rule%gauss = .true.
    call sensor_set_signals(solver, set, y, rule, linked_flux=linked)
    call close('saddle coil = M(1, 0.5, 0.3) I', y(set%index_of('SADDLE')), &
        mutual(1.0_dp, 0.5_dp, 0.3_dp) * current, 2.0e-5_dp)
    call close('flux loop = M(1, 1.5, -0.2) I', y(set%index_of('FLUX')), &
        mutual(1.0_dp, 1.5_dp, -0.2_dp) * current, 2.0e-5_dp)
    call close('diamagnetic loop = turns * linked flux', y(set%index_of('DIA')), &
        2.0_dp * linked, 1.0e-9_dp)
    call close('compensation loop = M(1, 0.8, 0) I', y(set%index_of('COMP')), &
        mutual(1.0_dp, 0.8_dp, 0.0_dp) * current, 2.0e-5_dp)
    call close('probe = mu0 I / (2 a) * area', y(set%index_of('PROBE_0001')), &
        0.5_dp * mu0 * current * 2.0e-3_dp, 2.0e-5_dp)
    call close('compensated diamagnetic combination', y(6), &
        y(set%index_of('DIA')) - 0.5_dp * y(set%index_of('COMP')), 1.0e-14_dp)

    ! Coil Jacobian: the signals are linear in EXTCUR (here one group, 1000 A).
    call sensor_set_coil_jacobian(solver, set, jac, rule)
    call sensor_set_linked_flux_row(set, row)
    call check(size(jac, 2) == 1, 'one coil group')
    do k = 1, set%n_signals()
        call close('coil Jacobian * EXTCUR + linked row', jac(k, 1) * solver%extcur(1) &
            + row(k) * linked, y(k), 1.0e-11_dp)
    end do

    ! VJP of the combination map: <ybar, C raw> = <C^T ybar, raw>.
    allocate(raw(set%n_raw()), yc(set%n_signals()), ybar(set%n_signals()), &
        rawbar(set%n_raw()))
    raw = [(sin(1.3_dp * k), k = 1, set%n_raw())]
    ybar = [(cos(0.7_dp * k), k = 1, set%n_signals())]
    call set%combine(raw, yc)
    call set%combine_vjp(ybar, rawbar)
    call close('combination VJP transposition', dot_product(ybar, yc), &
        dot_product(rawbar, raw), 1.0e-14_dp)

    if (command_argument_count() >= 2) then
        call get_command_argument(2, wout)
        call check_plasma(trim(wout))
    end if
    call check_rejections(trim(dir))
    if (failures > 0) then
        write(error_unit, '(I0,A)') failures, ' check(s) failed'
        error stop 1
    end if
    print '(A)', 'test_sensor_set passed'

contains

    subroutine check_plasma(path)
        !! Plasma part = mode Jacobian . coefficients + linked row . phiedge.
        character(len=*), intent(in) :: path
        real(dp), allocatable :: yv(:), yp(:), mode(:, :), shape(:, :)
        real(dp) :: phi

        call sensor_set_signals(solver, set, yv, rule, linked_flux=0.0_dp)
        call solver%enable_plasma_from_vmec(path, 32_i32, 32_i32)
        call sensor_set_signals(solver, set, yp, rule)
        phi = solver%plasma_flux()
        call sensor_set_plasma_jacobian(solver, set, mode, shape, rule)
        call check(size(shape, 1) == set%n_signals(), 'shape Jacobian rows')
        do k = 1, set%n_signals()
            call close('plasma signal = mode Jacobian . coefficients', yp(k) - yv(k), &
                dot_product(mode(k, :), solver%plasma%coefficients) + row(k) * phi, &
                1.0e-10_dp)
        end do
    end subroutine check_plasma

    subroutine check_rejections(d)
        character(len=*), intent(in) :: d
        call expect_error(d//'/bad_kind.sensors', 'tiago-sensors 1'//new_line('a')// &
            'group rogowski_coil loops.diagno', 'unknown sensor kind')
        call expect_error(d//'/bad_dia.sensors', 'tiago-sensors 1'//new_line('a')// &
            'group flux_loop loops.diagno select=DIA', 'declare it as diamagnetic_loop')
        call expect_error(d//'/bad_env.sensors', 'tiago-sensors 1'//new_line('a')// &
            'root ${TIAGO_SURELY_UNSET_VARIABLE}', 'is not set')
        call expect_error(d//'/bad_combo.sensors', 'tiago-sensors 1'//new_line('a')// &
            'group flux_loop loops.diagno select=FLUX'//new_line('a')// &
            'combine X 1.0 NOPE', 'unknown sensor NOPE')
        call expect_error(d//'/bad_header.sensors', 'group flux_loop loops.diagno', &
            'first line must be')
    end subroutine check_rejections

    subroutine expect_error(path, text, expected)
        character(len=*), intent(in) :: path, text, expected
        type(sensor_set_t) :: bad
        integer :: unit
        open(newunit=unit, file=path, status='replace', action='write')
        write(unit, '(A)') text
        close(unit)
        call read_sensor_set(path, bad, ierr, message)
        if (ierr == 0 .or. index(message, expected) == 0) then
            write(error_unit, '(A)') 'FAIL rejection '//path//': '//message
            failures = failures + 1
        end if
    end subroutine expect_error

    real(dp) function mutual(a, b, h)
        !! Mutual inductance of coaxial circles of radii a, b at distance h.
        real(dp), intent(in) :: a, b, h
        real(dp) :: k2, k, ek, ee
        k2 = 4.0_dp * a * b / ((a + b)**2 + h**2)
        k = sqrt(k2)
        call elliptic(k, ek, ee)
        mutual = mu0 * sqrt(a * b) * ((2.0_dp / k - k) * ek - 2.0_dp / k * ee)
    end function mutual

    subroutine elliptic(k, ek, ee)
        !! Complete elliptic integrals K(k), E(k) by the AGM.
        real(dp), intent(in) :: k
        real(dp), intent(out) :: ek, ee
        real(dp) :: a, b, c, an, s, p
        integer :: n
        a = 1.0_dp
        b = sqrt(1.0_dp - k * k)
        c = k
        s = 0.5_dp * c * c
        p = 0.5_dp
        do n = 1, 40
            an = 0.5_dp * (a + b)
            c = 0.5_dp * (a - b)
            b = sqrt(a * b)
            a = an
            p = 2.0_dp * p
            s = s + p * c * c
            if (abs(c) < 1.0e-17_dp) exit
        end do
        ek = pi / (2.0_dp * a)
        ee = ek * (1.0_dp - s)
    end subroutine elliptic

    subroutine write_geometry(d)
        character(len=*), intent(in) :: d
        integer :: unit, j
        real(dp) :: t

        open(newunit=unit, file=d//'/coil.coils', status='replace', action='write')
        write(unit, '(A)') 'periods 1', 'begin filament', 'mirror NIL'
        do j = 0, ncoil - 1
            t = 2.0_dp * pi * j / ncoil
            write(unit, '(4ES25.16)') cos(t), sin(t), 0.0_dp, current
        end do
        write(unit, '(4ES25.16,A)') 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, ' 1 CIRCLE'
        write(unit, '(A)') 'end'
        close(unit)

        open(newunit=unit, file=d//'/loops.diagno', status='replace', action='write')
        write(unit, '(I0)') 4
        call circle(unit, 'SADDLE', 0, 0.5_dp, 0.3_dp, .false.)
        call circle(unit, 'FLUX', 0, 1.5_dp, -0.2_dp, .false.)
        call circle(unit, 'DIA', 1, 0.5_dp, 3.0_dp, .true.)
        call circle(unit, 'COMP', 0, 0.8_dp, 0.0_dp, .false.)
        close(unit)

        open(newunit=unit, file=d//'/probes.txt', status='replace', action='write')
        write(unit, '(I0)') 1
        write(unit, '(6ES25.16)') 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.0e-3_dp
        close(unit)
    end subroutine write_geometry

    subroutine circle(unit, label, idia, radius, offset, poloidal)
        !! Horizontal circle at height offset, or a circle in the poloidal plane
        !! y = 0 centred on R = offset.
        integer, intent(in) :: unit, idia
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: radius, offset
        logical, intent(in) :: poloidal
        integer :: j
        real(dp) :: t

        write(unit, '(I0,1X,I0,1X,I0,1X,A)') nloop, 0, idia, label
        do j = 0, nloop - 1
            t = 2.0_dp * pi * j / nloop
            if (poloidal) then
                write(unit, '(3ES25.16)') offset + radius * cos(t), 0.0_dp, radius * sin(t)
            else
                write(unit, '(3ES25.16)') radius * cos(t), radius * sin(t), offset
            end if
        end do
    end subroutine circle

    subroutine close(name, a, b, tol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: a, b, tol
        if (.not. (abs(a - b) <= tol * max(abs(b), 1.0e-12_dp))) then
            write(error_unit, '(A,2ES24.16)') 'FAIL '//name//': ', a, b
            failures = failures + 1
        end if
    end subroutine close

    subroutine check(ok, name)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: name
        if (.not. ok) then
            write(error_unit, '(A)') 'FAIL '//name
            failures = failures + 1
        end if
    end subroutine check
end program test_sensor_set
