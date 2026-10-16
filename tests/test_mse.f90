program test_mse
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use tiago_diagnostic_types, only: mse_channel_t
    use tiago_mse, only: mse_geometry, mse_angles, read_mse_file
    use tiago_vacuum_forward, only: vacuum_solver_t
    use tiago_coil_kernels, only: build_coil_set
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), q = sqrt(0.5_dp)
    type(mse_channel_t) :: channels(3), ch(1), bad
    type(mse_channel_t), allocatable :: parsed(:)
    type(vacuum_solver_t) :: solver
    real(dp) :: field(3, 3), f(3, 1), p(3), a(6), gap, phi, er(3), ep(3)
    real(dp) :: radii(3), expected(3), reference, step, old, perturbed(3, 1)
    real(dp), allocatable :: angles(:), grad(:, :), plus(:), minus(:), jac(:, :)
    real(dp) :: v(3), n(3), z(3), h(3), e(3), num, den
    character(len=*), parameter :: row = &
        'GOOD 3 -1 0 1 1 0 6 0 0 1 0 0 0 0 1 0.001'

    integer :: j, k, unit
    integer(i32) :: ierr
    character(len=512) :: data_dir, example
    character(len=:), allocatable :: message

    call get_command_argument(1, data_dir)
    call get_command_argument(2, example)
    radii = [2.0_dp, 4.0_dp, 5.0_dp]
    ! Bphi = 2/R [T], Bz = 0.2 T. Hand-calculated atan(0.1 R).
    expected = [0.19739555984988076_dp, 0.38050637711236490_dp, &
        0.46364760900080612_dp]
    do j = 1, 3
        phi = (j - 1)*pi/3
        er = [cos(phi), sin(phi), 0.0_dp]
        ep = [-sin(phi), cos(phi), 0.0_dp]
        p = radii(j)*er
        channels(j) = channel_at(p, er + ep, er)
        field(:, j) = (2/radii(j))*ep + [0.0_dp, 0.0_dp, 0.2_dp]
        call mse_geometry(channels(j), p, a, ierr, message, gap)
        call success()
        call close_enough('midpoint', norm2(p - radii(j)*er), 0.0_dp, 2.0e-14_dp)
        call close_enough('intersection gap', gap, 0.0_dp, 2.0e-14_dp)
        call close_enough('A1', a(1), q, 2.0e-14_dp)
        call close_enough('A2,A3,A4 zero', norm2(a(2:4)), 0.0_dp, 2.0e-14_dp)
        call close_enough('A5', a(5), -q, 2.0e-14_dp)
        call close_enough('A6', a(6), q, 2.0e-14_dp)
    end do
    call mse_angles(channels, field, angles, ierr, message, grad)
    call success()
    do j = 1, 3
        call close_enough('toroidal + vertical angle', angles(j), &
            expected(j), 2.0e-14_dp)
    end do
    print '(A,3F19.15)', 'analytic angles [rad]: ', angles

    ch = channels(1:1)
    f = field(:, 1:1)
    call mse_angles(ch, -f, plus, ierr, message)
    call success()
    call close_enough('B reversal modulo pi', plus(1), angles(1), 2.0e-14_dp)
    f = 0.0_dp
    f(2, 1) = 1.0_dp
    call mse_angles(ch, f, plus, ierr, message)
    call success()
    call close_enough('pure toroidal', plus(1), 0.0_dp, 2.0e-14_dp)
    f(:, 1) = [0.0_dp, 0.0_dp, 1.0_dp]
    call mse_angles(ch, f, plus, ierr, message)
    call success()
    call close_enough('denominator zero', plus(1), -pi/2, 2.0e-14_dp)
    f(:, 1) = [0.0_dp, -1.0_dp, 0.2_dp]
    call mse_angles(ch, f, plus, ierr, message)
    call success()
    call close_enough('negative pitch', plus(1), -expected(1), 2.0e-14_dp)
    f = 0.0_dp
    call mse_angles(ch, f, plus, ierr, message)
    call rejected('zero field')
    f(:, 1) = [1.0_dp, 1.0_dp, 0.0_dp]
    call mse_angles(ch, f, plus, ierr, message)
    call rejected('field along beam')
    f(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call mse_angles(ch, f, plus, ierr, message)
    call rejected('nonfinite field')

    ! A fully 3-D oracle using the electric vector and Stokes Q,U directly.
    v = [1.0_dp, 2.0_dp, 3.0_dp]
    n = [-2.0_dp, 1.0_dp, 0.7_dp]
    ch(1) = channel_at([3.0_dp, -2.0_dp, 0.5_dp], v, n)
    n = n/norm2(n)
    z = [0.0_dp, 0.0_dp, 1.0_dp] - n(3)*n
    z = z/norm2(z)
    h = cross(n, z)
    f(:, 1) = [-0.3_dp, 1.2_dp, 0.15_dp]
    e = cross(v, f(:, 1))
    num = dot_product(e, h)
    den = dot_product(e, z)
    reference = 0.5_dp*atan2(2*num*den, den*den - num*num)
    call mse_angles(ch, f, plus, ierr, message, grad)
    call success()
    call close_enough('3-D Stokes oracle', plus(1), reference, 2.0e-14_dp)
    step = 1.0e-6_dp
    do k = 1, 3
        perturbed = f
        perturbed(k, 1) = perturbed(k, 1) + step
        call mse_angles(ch, perturbed, plus, ierr, message)
        call success()
        perturbed(k, 1) = perturbed(k, 1) - 2*step
        call mse_angles(ch, perturbed, minus, ierr, message)
        call success()
        call close_enough('Cartesian dgamma/dB FD', &
            (plus(1) - minus(1))/(2*step), grad(k, 1), 2.0e-9_dp)
    end do
    ch(1)%beam_direction = -10*ch(1)%beam_direction
    call mse_angles(ch, f, plus, ierr, message)
    call success()
    call close_enough('beam reversal and rescaling', plus(1), reference, 2.0e-14_dp)
    ch(1)%sight_direction = -3*ch(1)%sight_direction
    call mse_angles(ch, f, plus, ierr, message)
    call success()
    call close_enough('sightline reversal sign', plus(1), -reference, 2.0e-14_dp)

    ! Skew lines differ only vertically: closest-approach midpoint z=0.01.
    bad = channels(1)
    bad%sight_point(3) = bad%sight_point(3) + 0.02_dp
    bad%max_gap = 0.03_dp
    call mse_geometry(bad, p, a, ierr, message, gap)
    call success()
    call close_enough('skew midpoint z', p(3), 0.01_dp, 2.0e-14_dp)
    call close_enough('skew gap', gap, 0.02_dp, 2.0e-14_dp)
    bad%max_gap = 0.001_dp
    call mse_geometry(bad, p, a, ierr, message)
    call rejected('excess gap')
    bad = channels(1)
    bad%sight_direction = bad%beam_direction
    call mse_geometry(bad, p, a, ierr, message)
    call rejected('parallel lines')
    bad%sight_direction(3) = 1.0e-8_dp
    call mse_geometry(bad, p, a, ierr, message)
    call rejected('nearly parallel lines')
    bad = channels(1)
    bad%reference = bad%sight_direction
    call mse_geometry(bad, p, a, ierr, message)
    call rejected('reference along sightline')
    bad = channels(1)
    bad%beam_direction = 0.0_dp
    call mse_geometry(bad, p, a, ierr, message)
    call rejected('zero beam direction')
    bad = channel_at([0.0_dp, 0.0_dp, 0.1_dp], &
        [1.0_dp, 1.0_dp, 0.0_dp], [1.0_dp, 0.0_dp, 0.0_dp])
    call mse_geometry(bad, p, a, ierr, message)
    call rejected('cylindrical axis')

    call read_mse_file(trim(example), parsed, ierr, message)
    call success()
    if (size(parsed) /= 2) error stop 'MSE example count'
    call mse_geometry(parsed(1), p, a, ierr, message)
    call success()
    call close_enough('example midpoint', norm2(p - [4.0_dp, 0.0_dp, 0.0_dp]), &
        0.0_dp, 2.0e-14_dp)
    open (newunit=unit, file='invalid_mse.txt', status='replace')
    write (unit, '(A)') '1', 'BAD 1 2 3'
    close (unit)
    call read_mse_file('invalid_mse.txt', parsed, ierr, message)
    call rejected('truncated geometry file')
    if (allocated(parsed)) error stop 'failed parse published partial geometry'

    call parser_rejects('0', '', '')
    call parser_rejects('1', row//' 99', '')
    call parser_rejects('1', row, row)
    call parser_rejects('2', row, row)
    call parser_rejects('1', 'BAD 3 -1 0 0 0 0 6 0 0 1 0 0 0 0 1 0.001', '')
    call parser_rejects('1', 'BAD 3 -1 0 NaN 1 0 6 0 0 1 0 0 0 0 1 0.001', '')
    call parser_rejects('1', 'BAD 3 -1 0 1 1 0 6 0 0 1 0 0 0 0 1 -1', '')

    ! Two orthogonal square coils: Bz=C I1, Bphi=C I2 at their common center.
    ! Thus gamma=atan(I1/I2), J=(I2,-I1)/(I1**2+I2**2).
    ch(1) = channel_at([4.0_dp, 0.0_dp, 0.0_dp], &
        [1.0_dp, 1.0_dp, 0.0_dp], [1.0_dp, 0.0_dp, 0.0_dp])
    call solver%init(trim(data_dir)//'/mse_coils.coils')
    call solver%mse(ch, angles, ierr, message, jac)
    call success()
    call close_enough('square-coil analytic angle', angles(1), expected(3), 2.0e-14_dp)
    call close_enough('coil J1', jac(1, 1), 0.4_dp, 2.0e-14_dp)
    call close_enough('coil J2', jac(1, 2), -0.2_dp, 2.0e-14_dp)
    do j = 1, 2
        old = solver%extcur(j)
        solver%extcur(j) = old + step
        call rebuild_coils()
        call solver%mse(ch, plus, ierr, message)
        call success()
        solver%extcur(j) = old - step
        call rebuild_coils()
        call solver%mse(ch, minus, ierr, message)
        call success()
        solver%extcur(j) = old
        call rebuild_coils()
        call close_enough('coil EXTCUR Jacobian FD', &
            (plus(1) - minus(1))/(2*step), jac(1, j), 2.0e-9_dp)
    end do
    print '(A,F19.15,A,2F19.15)', 'square-coil angle: ', angles(1), ' J: ', jac
    call solver%finalize()
    print '(A)', 'test_mse passed'
contains
    subroutine parser_rejects(count, first, second)
        character(len=*), intent(in) :: count, first, second
        open (newunit=unit, file='invalid_mse.txt', status='replace')
        write (unit, '(A)') count
        if (len(first) > 0) write (unit, '(A)') first
        if (len(second) > 0) write (unit, '(A)') second
        close (unit)
        call read_mse_file('invalid_mse.txt', parsed, ierr, message)
        call rejected('malformed channel file')
        if (allocated(parsed)) error stop 'partial geometry after failed parse'
    end subroutine parser_rejects

    function channel_at(point, beam, sight) result(c)
        real(dp), intent(in) :: point(3), beam(3), sight(3)
        type(mse_channel_t) :: c
        c%label = 'TEST'
        c%beam_point = point - beam
        c%beam_direction = beam
        c%sight_point = point + 2*sight
        c%sight_direction = sight
    end function channel_at

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross

    subroutine rebuild_coils()
        solver%coils = build_coil_set(solver%x, solver%y, solver%z, &
            solver%unit_current*solver%extcur(solver%group))
    end subroutine rebuild_coils

    subroutine success()
        if (ierr /= 0) then
            print '(A)', message
            error stop 'unexpected MSE error'
        end if
    end subroutine success

    subroutine rejected(name)
        character(len=*), intent(in) :: name
        if (ierr == 0) then
            print '(A)', name
            error stop 'invalid MSE input accepted'
        end if
    end subroutine rejected

    subroutine close_enough(name, actual, expected, atol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual, expected, atol
        if (abs(actual - expected) <= atol) return
        print '(A,2ES24.16)', name//': ', actual, expected
        error stop 1
    end subroutine close_enough
end program test_mse
