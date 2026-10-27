program test_ece_geometry
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use ece_geometry_fixture, only: setup, setup_vmecpp
    use tiago_ece_state_medium, only: state_medium_t
    use tiago_vmec_state, only: vmec_state_t
    use tiago_state_sampler, only: medium_sample_jvp, medium_sample_vjp
    use tiago_ece_ray, only: cold_ray_t, trace_cold_ray, cold_hamiltonian, &
                             cold_ray_jvp, cold_ray_vjp, ECE_O_MODE, ECE_X_MODE
    use fortnum_status, only: fortnum_status_t
    use tiago_vmecpp_state, only: vmecpp_state_t
    implicit none
    type(state_medium_t) :: medium, shifted
    type(cold_ray_t) :: ray, plus, minus
    type(fortnum_status_t) :: status
    real(dp), parameter :: freq = 9.0e10_dp, step = 1.0e-6_dp
    real(dp) :: x(3), q(4), dq(3, 4), d2q(3, 3, 4), qa(4), da(3, 4), ha(3, 3, 4)
    real(dp) :: qb(4), db(3, 4), hb(3, 3, 4), expected(4), eg(3, 4), eh(3, 3, 4)
    real(dp) :: r, s, sg(3), sh(3, 3), er(3), rr(3, 3), ne, scale(4), eye(3, 3)
    real(dp) :: err, emax, xdot(3), qdot(4), gdot(3, 4), bar(4), gbar(3, 4), xbar(3)
    real(dp) :: initial(6), direction(3), dh(8), d2h(8, 8), h, fd(6), ydot(6)
    real(dp) :: ybar(6), ibar(6), fbar, expect6(6), duration
    real(dp), allocatable :: qp(:, :), dqp(:, :, :), ap(:, :), dap(:, :, :)
    real(dp), allocatable :: p0(:), v(:), pbar(:), zero(:)
    integer :: np, ns, i, j, mode
    character(1024) :: path
    logical :: product_state
    product_state = .false.
    if (command_argument_count() > 0) then
        call get_command_argument(1, path)
        product_state = trim(path) == 'vmecpp'
    end if
    call setup(medium)
    if (product_state) call setup_vmecpp(medium)
    ns = medium%state%n_params()
    np = medium%n_params()
    allocate (qp(4, np), dqp(3, 4, np), ap(4, np), dap(3, 4, np), p0(np), v(np), &
              pbar(np), zero(np))
    zero = 0
    call medium%state%get_params(p0(:ns))
    p0(ns + 1:) = medium%z
    scale = [1.0_dp, 1.0_dp, 1.0_dp, 1.0e-19_dp]
    x = [3.4_dp*cos(0.3_dp), 3.4_dp*sin(0.3_dp), 0.2_dp]
    call medium%sample(x, q, dq, d2q, qp, dqp)
    r = norm2(x(:2))
    s = (r - 3)**2 + x(3)**2
    ne = 1.0e19_dp*exp(-0.4_dp*s)
    expected = [-0.3_dp*x(2), 0.3_dp*x(1), 0.0_dp, ne]
    eg = 0
    eg(2, 1) = -0.3_dp
    eg(1, 2) = 0.3_dp
    er = [x(1)/r, x(2)/r, 0.0_dp]
    eye = 0
    do i = 1, 3
        eye(i, i) = 1
    end do
    rr = eye
    rr(3, 3) = 0
    do i = 1, 3
        rr(:, i) = (rr(:, i) - er*er(i))/r
    end do
    sg = 2*(r - 3)*er + [0.0_dp, 0.0_dp, 2*x(3)]
    sh = 2*(r - 3)*rr
    do i = 1, 3
        sh(:, i) = sh(:, i) + 2*er*er(i)
    end do
    sh(3, 3) = 2
    eh = 0
    eg(:, 4) = -0.4_dp*ne*sg
    do i = 1, 3
        eh(:, i, 4) = ne*(0.16_dp*sg*sg(i) - 0.4_dp*sh(:, i))
    end do
    err = maxval(abs((q - expected)*scale))
    call check(err < 2.0e-12_dp, 'analytic circular field/density value', err)
    err = maxval(abs((dq - eg)*spread(scale, 1, 3)))
    call check(err < 2.0e-11_dp, 'analytic first spatial derivative [scaled/m]', err)
    err = maxval(abs((d2q - eh)*spread(spread(scale, 1, 3), 1, 3)))
    call check(err < 2.0e-10_dp, 'analytic second spatial derivative [scaled/m2]', err)
    emax = 0
    do i = 1, 3
        call medium%sample(x + step*eye(:, i), qa, da, ha, ap, dap)
        call medium%sample(x - step*eye(:, i), qb, db, hb, ap, dap)
        err = maxval(abs(((qa - qb)/(2*step) - dq(i, :))*scale))
        emax = max(emax, err)
        err = maxval(abs(((da - db)/(2*step) - d2q(:, i, :))*spread(scale, 1, 3)))
        emax = max(emax, err)
    end do
    call check(emax < 1.0e-7_dp, 'spatial jets versus central differences', emax)
    ! Nonzero poloidal and toroidal Fourier perturbations exercise mixed jets.
    select type (st => medium%state)
    type is (vmec_state_t)
        st%bsupumnc(1, :) = 0.08_dp
        st%bsupumnc(2, :) = 0.02_dp*st%rho_half
        st%xn(2) = 1
        st%xn_nyq(2) = 1
    type is (vmecpp_state_t)
        st%bsupumnc(1, :) = 0.08_dp
        st%bsupumnc(2, :) = 0.02_dp*st%rho_half
        st%coef(1, 1, :, 1) = 0.02_dp*st%rho_full
        st%coef(1, 1, :, 2) = -0.01_dp*st%rho_full
        st%coef(1, 1, :, 5) = 0.015_dp*st%rho_full
        st%coef(1, 1, :, 6) = 0.01_dp*st%rho_full
        st%xn_nyq(2) = 1
    end select
    call medium%state%get_params(p0(:ns))
    call medium%sample(x, q, dq, d2q, qp, dqp)
    emax = 0
    do j = 1, np
        shifted = medium
        v = 0
        v(j) = max(0.03_dp, abs(p0(j)))
        call change(shifted, p0 + step*v)
        call shifted%sample(x, qa, da, ha, ap, dap)
        call change(shifted, p0 - step*v)
        call shifted%sample(x, qb, db, hb, ap, dap)
        err = maxval(abs(((qa - qb)/(2*step) - qp(:, j)*v(j))*scale))
        emax = max(emax, err)
        err = maxval(abs(((da - db)/(2*step) - dqp(:, :, j)*v(j))*spread(scale, 1, 3)))
        emax = max(emax, err)
    end do
    call check(emax < 2.0e-7_dp, &
               'all state/profile parameter columns versus finite differences', emax)
    v = [(0.03_dp*sin(real(i, dp)), i=1, np)]
    xdot = [0.02_dp, -0.03_dp, 0.01_dp]
    bar = [0.3_dp, -0.2_dp, 0.1_dp, 0.4e-19_dp]
    gbar = spread(bar, 1, 3)
    call medium_sample_jvp(dq, d2q, qp, dqp, xdot, v, qdot, gdot)
    call medium_sample_vjp(dq, d2q, qp, dqp, bar, gbar, xbar, pbar)
    err = abs(dot_product(bar, qdot) + sum(gbar*gdot) - dot_product(xbar, &
                                                                    xdot) - dot_product(pbar, v))
    call check(err < 2.0e-12_dp, 'sample value/gradient JVP/VJP identity', err)
    ! Moving paths, field geometry, density latents and frequency, O and X modes.
    initial(:3) = x
    direction = [0.8_dp, 0.3_dp, 0.2_dp]
    direction = direction/norm2(direction)
    duration = 0.16_dp
    do mode = ECE_O_MODE, ECE_X_MODE
        call medium%sample(x, q, dq, d2q, qp, dqp)
        call cold_hamiltonian(direction, q(:3), q(4), freq, mode, h, dh, d2h)
        initial(4:) = sqrt(1 - 2*h)*direction
        call trace_cold_ray(medium, freq, mode, initial, duration, ray, status, &
                            rtol=1.0e-11_dp, atol=1.0e-13_dp, hmax=0.02_dp)
        call check(status%code == 0, 'state-medium cold ray '//status%msg, &
                   real(status%code, dp))
        call cold_ray_jvp(medium, ray, 0*initial, v, 0.001_dp*freq, ydot, status)
        call check(status%code == 0, 'path JVP status', real(status%code, dp))
        shifted = medium
        call change(shifted, p0 + step*v)
        call trace_cold_ray(shifted, freq*(1 + step*0.001_dp), mode, initial, &
                            duration, &
                            plus, status, rtol=1.0e-11_dp, atol=1.0e-13_dp, &
                            hmax=0.02_dp)
        call check(status%code == 0, 'positive path difference '//status%msg, &
                   real(status%code, dp))
        call change(shifted, p0 - step*v)
        call trace_cold_ray(shifted, freq*(1 - step*0.001_dp), mode, initial, &
                            duration, &
                            minus, status, rtol=1.0e-11_dp, atol=1.0e-13_dp, &
                            hmax=0.02_dp)
        call check(status%code == 0, 'negative path difference '//status%msg, &
                   real(status%code, dp))
        fd = (plus%trace%y(:, size(plus%trace%t)) - minus%trace%y(:, &
                                                                  size(minus%trace%t)))/(2*step)
        err = maxval(abs(fd - ydot))/max(1.0e-6_dp, maxval(abs(ydot)))
        call check(err < 3.0e-5_dp, &
                   'O/X composed moving-path finite difference relative error', err)
        call medium%sample(ray%trace%y(:3, size(ray%trace%t)), q, dq, d2q, qp, dqp)
        call medium_sample_jvp(dq, d2q, qp, dqp, ydot(:3), v, qdot, gdot)
        call change(shifted, p0 + step*v)
        call shifted%sample(plus%trace%y(:3, size(plus%trace%t)), qa, da, ha, ap, dap)
        call change(shifted, p0 - step*v)
        call shifted%sample(minus%trace%y(:3, size(minus%trace%t)), qb, db, hb, ap, dap)
        err = maxval(abs(((qa - qb)/(2*step) - qdot)*scale))
        call check(err < 2.0e-8_dp, 'profiles sampled on moving O/X endpoints', err)
        ybar = [0.3_dp, 0.4_dp, -0.2_dp, 0.1_dp, -0.3_dp, 0.2_dp]
        call cold_ray_vjp(medium, ray, ybar, ibar, pbar, fbar, status)
        call check(status%code == 0, 'path VJP status', real(status%code, dp))
        err = abs(dot_product(ybar, ydot) - dot_product(pbar, v) - fbar*0.001_dp*freq)
        call check(err < 1.0e-10_dp, 'O/X composed path transpose identity', err)
    end do
    call setup(medium)
    medium%profiles%p(1)%unit = 0
    initial = [3.4_dp, 0.0_dp, 0.2_dp, 0.8_dp, 0.0_dp, 0.6_dp]
    call trace_cold_ray(medium, freq, ECE_O_MODE, initial, 0.15_dp, ray, status)
    call check(status%code == 0, 'vacuum state ray', real(status%code, dp))
    expect6 = initial
    expect6(:3) = initial(:3) + 0.15_dp*initial(4:)
    err = maxval(abs(ray%trace%y(:, size(ray%trace%t)) - expect6))
    call check(err < 2.0e-12_dp, 'vacuum straight ray [m/1]', err)
    call medium%sample([30.0_dp, 0.0_dp, 0.0_dp], q, dq, d2q, qp, dqp)
    call check(all(ieee_is_nan(q)), 'outside chart rejected', 0.0_dp)
    if (command_argument_count() > 0 .and. .not. product_state) then
        call get_command_argument(1, path)
        select type (st => medium%state)
        type is (vmec_state_t)
            block
                type(vmec_state_t) :: loaded
                integer :: ierr
                character(:), allocatable :: msg
                call loaded%load(trim(path), ierr, msg, vector=.true.)
                call check(ierr == 0, 'load published VMEC wout '//msg, real(ierr, dp))
                medium%state = loaded
            end block
        end select
        call medium%state%eval([0.47_dp, 0.6_dp, 0.2_dp], x, rr, ne, er)
        deallocate (qp, dqp, ap, dap)
        np = medium%n_params()
        allocate (qp(4, np), dqp(3, 4, np), ap(4, np), dap(3, 4, np))
        medium%profiles%p(1)%unit = 1.0e19_dp
        call medium%sample(x, q, dq, d2q, qp, dqp)
        emax = 0
        do i = 1, 3
            call medium%sample(x + step*eye(:, i), qa, da, ha, ap, dap)
            call medium%sample(x - step*eye(:, i), qb, db, hb, ap, dap)
            err = maxval(abs(((qa - qb)/(2*step) - dq(i, :))*scale))
            emax = max(emax, err)
            err = maxval(abs(((da - db)/(2*step) - d2q(:, i, :))*spread(scale, 1, 3)))
            emax = max(emax, err)
        end do
        call check(emax < 1.0e-5_dp, 'published VMEC spatial jet finite differences', &
                   emax)
        initial(:3) = x
        direction = [0.8_dp, 0.3_dp, 0.2_dp]
        direction = direction/norm2(direction)
        do mode = ECE_O_MODE, ECE_X_MODE
            call cold_hamiltonian(direction, q(:3), q(4), freq, mode, h, dh, d2h)
            initial(4:) = sqrt(1 - 2*h)*direction
            call trace_cold_ray(medium, freq, mode, initial, 0.002_dp, ray, status)
            call check(status%code == 0, 'published VMEC O/X ray '//status%msg, &
                       real(status%code, dp))
            call check(ray%max_hamiltonian_error < 2.0e-8_dp, &
                       'VMEC Hamiltonian drift', &
                       ray%max_hamiltonian_error)
        end do
    end if
contains
    subroutine check(ok, label, error)
        logical, intent(in)::ok
        character(*), intent(in)::label
        real(dp), intent(in)::error
        print '(A,1X,ES14.6)', trim(label), error
        if (.not. ok) error stop 'ECE geometry test failed'
    end subroutine
    subroutine change(m, p)
        type(state_medium_t), intent(inout)::m
        real(dp), intent(in)::p(:)
        integer::n
        n = m%state%n_params()
        call m%state%set_params(p(:n))
        m%z = p(n + 1:)
    end subroutine
end program
