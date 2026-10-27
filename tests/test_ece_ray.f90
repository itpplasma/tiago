
program test_ece_ray
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use ece_ray_fixture, only: polynomial_medium_t
    use tiago_ece_ray, only: cold_ray_t, cold_hamiltonian, trace_cold_ray, &
                     cold_ray_jvp, cold_ray_vjp, trace_cold_rays, ECE_O_MODE, ECE_X_MODE
    use fortnum_status, only: fortnum_status_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), frequency = 9.0e10_dp
    real(dp), parameter :: electron_mass = 9.1093837015e-31_dp
    real(dp), parameter :: charge = 1.602176634e-19_dp, &
                           eps0 = 8.8541878128e-12_dp
    real(dp), &
        parameter :: critical = (2*pi*frequency)**2*eps0*electron_mass/charge**2
    type(polynomial_medium_t) :: medium, plus, minus
    type(cold_ray_t) :: ray, rp, rm
    type(fortnum_status_t) :: status
    real(dp) :: initial(6), direction(3), expected(6), error, h, dh(8), &
                d2h(8, 8)
    real(dp) :: xx, yy, n2, t, duration, turning_time, norm, step
    real(dp) :: q(4), dq(3, 4), d2q(3, 3, 4), qp(4, 8), dqp(3, 4, 8)
    real(dp) :: initial_dot(6), dp_(8), dfreq, final_dot(6), fd(6), &
                final_bar(6)
    real(dp) :: initial_bar(6), pbar(8), fbar, qdot(4), lhs, rhs, ip(6), &
                im(6)
    integer :: mode, k

    ! Appleton-Hartree, collisionless electron plasma, ions fixed.
    ! Stix, Waves in Plasmas (1992), Ch. 1: O: P; X: RL/S at k perpendicular B.
    xx = 0.15_dp
    medium%p = 0.0_dp
    medium%p(3) = 1.0_dp
    medium%p(4) = xx*critical
    yy = (charge/(electron_mass*2*pi*frequency))**2
    direction = [1.0_dp, 0.0_dp, 0.0_dp]
    do mode = ECE_O_MODE, ECE_X_MODE
        n2 = 1 - xx
        if (mode == ECE_X_MODE) n2 = 1 - xx*(1 - xx)/(1 - xx - yy)
        call cold_hamiltonian(sqrt(n2)*direction, medium%p(:3), &
                              medium%p(4), frequency, mode, h, dh, d2h)
        call check(abs(h) < 2.0e-14_dp, 'perpendicular dispersion')
    end do
    print '(A,ES12.4)', 'perpendicular dispersion |H|: ', abs(h)

    medium%p(4:) = 0.0_dp
    initial = [0.2_dp, -0.1_dp, 0.3_dp, 0.8_dp, 0.6_dp, 0.0_dp]
    do mode = ECE_O_MODE, ECE_X_MODE
        call trace_cold_ray(medium, frequency, mode, initial, 2.0_dp, ray, &
                            status)
        call check(status%code == 0, 'vacuum trace '//status%msg)
        expected = initial
        expected(:3) = initial(:3) + 2*initial(4:)
        error = maxval(abs(ray%trace%y(:, size(ray%trace%t)) - expected))
        call check(error < 3.0e-12_dp, 'straight vacuum ray')
        print '(A,I0,A,ES12.4)', 'vacuum mode ', mode, &
            ' position/index error [m/1]: ', error
    end do

    call check_batch()

    ! O-mode in a linear-density slab, B perpendicular to the ray plane:
    ! X=g*x, x=x0+Nx0*t-g*t^2/4, Nx=Nx0-g*t/2. The turning point has
    ! X=1-Ny^2; the normal-incidence cutoff is x=1/g.
    medium%p = 0.0_dp
    medium%p(3) = 1.0_dp
    medium%p(5) = 0.6_dp*critical
    initial = [0.2_dp, 0.0_dp, 0.0_dp, sqrt(0.72_dp), 0.4_dp, 0.0_dp]
    turning_time = 2*initial(4)/0.6_dp
    call trace_cold_ray(medium, frequency, ECE_O_MODE, initial, turning_time, &
                        ray, status)
    call check(status%code == 0, 'slab turning trace '//status%msg)
    error = abs(ray%trace%y(1, size(ray%trace%t)) - 1.4_dp)
    call check(error < 2.0e-9_dp, 'slab analytic turning point')
    call check(abs(ray%trace%y(4, size(ray%trace%t))) < 2.0e-9_dp, &
               'slab zero normal index')
    print '(A,ES12.4)', 'slab turning-point error [m]: ', error
    duration = 1.5_dp*turning_time
    call trace_cold_ray(medium, frequency, ECE_O_MODE, initial, duration, &
                        ray, status)
    call check(status%code == 0, 'slab reflected trace '//status%msg)
    error = 0.0_dp
    do k = 1, size(ray%trace%t)
        t = ray%trace%t(k)
        expected = initial
        expected(1) = initial(1) + initial(4)*t - 0.6_dp*t*t/4
        expected(2) = initial(2) + initial(5)*t
        expected(4) = initial(4) - 0.6_dp*t/2
        error = max(error, maxval(abs(expected - ray%trace%y(:, k))))
    end do
    call check(error < 3.0e-9_dp, 'slab analytic reflected ray')
    call check(ray%trace%y(4, size(ray%trace%t)) < 0.0_dp, &
               'slab reflection')
    print '(A,ES12.4)', 'slab whole-ray position/index error [m/1]: ', error

    medium%p(3) = 0.0_dp
    initial = [0.2_dp, 0.0_dp, 0.0_dp, sqrt(0.88_dp), 0.0_dp, 0.0_dp]
    turning_time = 2*initial(4)/0.6_dp
    call trace_cold_ray(medium, frequency, ECE_O_MODE, initial, turning_time, &
                        ray, status)
    call check(status%code == 0, 'normal cutoff trace '//status%msg)
    error = abs(ray%trace%y(1, size(ray%trace%t)) - 1/0.6_dp)
    call check(error < 2.0e-9_dp, 'normal-incidence analytic cutoff')
    print '(A,ES12.4)', 'normal-incidence cutoff error [m]: ', error
    call cold_ray_jvp(medium, ray, 0*initial, 0*medium%p, &
                      0.0_dp, final_dot, status)
    call check(status%code /= 0, 'active field-null derivative rejected')
    medium%p(3) = 1.0_dp
    ! A cylindrical parabolic density gives a harmonic oscillator with unit
    ! scale length: x=x0*cos(t)+N0*sin(t), N=N0*cos(t)-x0*sin(t).
    medium%p(5) = 0.0_dp
    medium%p(8) = critical
    initial = [-0.7_dp, 0.1_dp, 0.0_dp, 0.7_dp, 0.1_dp, 0.0_dp]
    call trace_cold_ray(medium, frequency, ECE_O_MODE, initial, 2*pi, ray, &
                        status, &
                        rtol=1.0e-11_dp, atol=1.0e-13_dp)
    call check(status%code == 0, 'cylinder trace '//status%msg)
    error = 0.0_dp
    do k = 1, size(ray%trace%t)
        t = ray%trace%t(k)
        expected(:3) = initial(:3)*cos(t) + initial(4:)*sin(t)
        expected(4:) = initial(4:)*cos(t) - initial(:3)*sin(t)
        error = max(error, maxval(abs(expected - ray%trace%y(:, k))))
    end do
    call check(error < 4.0e-10_dp, 'cylinder analytic orbit')
    print '(A,ES12.4)', 'parabolic cylinder whole-ray error [m/1]: ', error

    ! A fully oblique ray in a nonuniform synthetic medium checks derivatives
    ! in the field, density, initial position/index and observing frequency.
    medium%p = [0.2_dp, 0.3_dp, 1.2_dp, 0.08_dp*critical, &
                0.01_dp*critical, -0.005_dp*critical, 0.002_dp*critical, &
                0.02_dp*critical]
    direction = [0.8_dp, 0.55_dp, 0.23_dp]
    direction = direction/norm2(direction)
    dp_ = [0.01_dp, -0.02_dp, 0.03_dp, 0.004_dp*critical, &
           0.001_dp*critical, -0.002_dp*critical, 0.0003_dp*critical, &
           0.003_dp*critical]
    dfreq = 0.005_dp*frequency
    initial(:3) = [0.1_dp, -0.2_dp, 0.3_dp]
    initial_dot(:3) = [0.02_dp, 0.01_dp, -0.03_dp]
    step = 1.0e-4_dp
    final_bar = [0.7_dp, -0.2_dp, 0.4_dp, -0.3_dp, 0.6_dp, 0.5_dp]
    do mode = ECE_O_MODE, ECE_X_MODE
        call launch(medium, frequency, mode, initial)
        call medium%sample(initial(:3), q, dq, d2q, qp, dqp)
        call cold_hamiltonian(initial(4:), q(:3), q(4), frequency, mode, h, &
                              dh, d2h)
        qdot = matmul(transpose(dq), initial_dot(:3)) + matmul(qp, dp_)
        norm = -(dot_product(dh(4:7), &
                             qdot) + dh(8)*dfreq)/dot_product(dh(:3), direction)
        initial_dot(4:) = norm*direction
        call trace_cold_ray(medium, frequency, mode, initial, 1.5_dp, ray, &
                            status, &
                            rtol=1.0e-11_dp, atol=1.0e-13_dp)
        call check(status%code == 0, 'oblique trace '//status%msg)
        call cold_ray_jvp(medium, ray, initial_dot, dp_, dfreq, final_dot, &
                          status)
        call check(status%code == 0, 'ray JVP '//status%msg)
        plus%p = medium%p + step*dp_
        minus%p = medium%p - step*dp_
        ip(:3) = initial(:3) + step*initial_dot(:3)
        im(:3) = initial(:3) - step*initial_dot(:3)
        call launch(plus, frequency + step*dfreq, mode, ip)
        call launch(minus, frequency - step*dfreq, mode, im)
        call trace_cold_ray(plus, frequency + step*dfreq, mode, ip, 1.5_dp, &
                            rp, status, &
                            rtol=1.0e-11_dp, atol=1.0e-13_dp)
        call check(status%code == 0, 'positive FD trace '//status%msg)
        call trace_cold_ray(minus, frequency - step*dfreq, mode, im, 1.5_dp, &
                            rm, status, &
                            rtol=1.0e-11_dp, atol=1.0e-13_dp)
        call check(status%code == 0, 'negative FD trace '//status%msg)
        fd = (rp%trace%y(:, size(rp%trace%t)) - rm%trace%y(:, &
                                                           size(rm%trace%t)))/(2*step)
        error = maxval(abs(fd - final_dot))/maxval(abs(fd))
        call check(error < 2.0e-7_dp, 'ray finite differences')
        print '(A,I0,A,ES12.4)', 'oblique mode ', mode, &
            ' JVP relative FD error: ', error
        call cold_ray_vjp(medium, ray, final_bar, initial_bar, pbar, fbar, &
                          status)
        call check(status%code == 0, 'ray VJP '//status%msg)
        lhs = dot_product(final_dot, final_bar)
        rhs = dot_product(initial_dot, initial_bar) + dot_product(dp_, &
                                                                  pbar) + dfreq*fbar
        error = abs(lhs - rhs)/max(abs(lhs), abs(rhs))
        call check(error < 2.0e-12_dp, 'ray transpose identity')
        print '(A,I0,A,ES12.4)', 'oblique mode ', mode, &
            ' adjoint relative error: ', error
    end do

    medium%p = 0.0_dp
    medium%p(3) = sqrt(2.0_dp)*electron_mass*2*pi*frequency/charge
    medium%p(4) = 0.1_dp*critical
    call launch(medium, frequency, ECE_X_MODE, initial)
    call trace_cold_ray(medium, frequency, ECE_X_MODE, initial, 0.2_dp, ray, &
                        status)
    call check(status%code /= 0, 'slow branch |N| > 1 rejected')
    initial(4:) = [2.0_dp, 0.0_dp, 0.0_dp]
    call trace_cold_ray(medium, frequency, ECE_O_MODE, initial, 1.0_dp, ray, &
                        status)
    call check(status%code /= 0, 'off-shell launch rejected')
    call trace_cold_ray(medium, -frequency, ECE_O_MODE, initial, 1.0_dp, &
                        ray, status)
    call check(status%code /= 0, 'negative frequency rejected')
contains
    subroutine check_batch()
        type(cold_ray_t), allocatable :: rays(:)
        type(fortnum_status_t), allocatable :: statuses(:)
        type(cold_ray_t) :: reference
        type(fortnum_status_t) :: result, single_status
        real(dp) :: starts(6, 4), frequencies(4), durations(4), difference
        integer :: modes(4), i
        do i = 1, 4
            starts(:, i) = initial
        end do
        frequencies = frequency
        durations = 0.2_dp
        modes = [ECE_O_MODE, ECE_X_MODE, ECE_O_MODE, ECE_X_MODE]
        call trace_cold_rays(medium, frequencies, modes, starts, durations, &
                             rays, statuses, result)
        call check(result%code == 0, 'batch status')
        difference = 0.0_dp
        do i = 1, 4
            call check(statuses(i)%code == 0, 'per-ray batch status')
            call trace_cold_ray(medium, frequency, modes(i), starts(:, i), &
                                durations(i), reference, single_status)
            call check(single_status%code == 0, 'serial batch reference')
            difference = max(difference, maxval(abs( &
                                         reference%trace%y(:, size(reference%trace%t)) &
                                          - rays(i)%trace%y(:, size(rays(i)%trace%t)))))
        end do
        call check(difference < 1.0e-14_dp, 'batch identity')
        frequencies(2) = -frequency
        call trace_cold_rays(medium, frequencies, modes, starts, durations, &
                             rays, statuses, result)
        call check(result%code /= 0, 'batch propagates failure')
        call check(statuses(2)%code /= 0, 'failed ray status')
        call check(statuses(1)%code == 0, 'independent ray still succeeds')
        call trace_cold_rays(medium, frequencies, modes, starts(:, :3), &
                             durations, rays, statuses, result)
        call check(result%code /= 0, 'batch shape validation')
        call trace_cold_rays(medium, frequencies(:0), modes(:0), starts(:, &
                                                                        :0), &
                             durations(:0), rays, statuses, result)
        call check(result%code == 0, 'empty batch')
    end subroutine
    subroutine launch(m, freq, mode, y)
        type(polynomial_medium_t), intent(in) :: m
        real(dp), intent(in) :: freq
        integer, intent(in) :: mode
        real(dp), intent(inout) :: y(6)
        real(dp) :: local_q(4), grad(3, 4), hess(3, 3, 4), p(4, 8), gp(3, 4, 8)
        real(dp) :: h0, g(8), hh(8, 8)
        call m%sample(y(:3), local_q, grad, hess, p, gp)
        call cold_hamiltonian(direction, local_q(:3), local_q(4), freq, &
                              mode, h0, g, hh)
        y(4:) = sqrt(1 - 2*h0)*direction
    end subroutine
    subroutine check(ok, message)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: message
        if (.not. ok) then
            print *, 'FAIL: '//message
            error stop 1
        end if
    end subroutine
end program
