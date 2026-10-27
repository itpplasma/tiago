program test_ece_transfer
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use tiago_ece_transfer, only: ece_transfer_t, thermal_transfer, &
        thermal_transfer_batch
    use tiago_ece_receiver, only: ece_receiver_t, ece_response_channel_t
    use tiago_calibration, only: calibration_prior_t
    use tiago_ece_radiometer_noise, only: ece_radiometer_noise
    implicit none
    integer :: failures = 0
    call limits()
    call layered_wall()
    call transfer_derivatives()
    call convergence()
    call receiver()
    call composed()
    call noise_kernel()
    call invalid()
    if (failures /= 0) error stop 'ECE transfer failed'
    print '(A)', 'ECE transport and receiver: PASS'
contains
    subroutine check(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        print '(A)', 'FAIL: '//message
        failures = failures + 1
    end subroutine

    subroutine limits()
        ! Formal constant-source transfer solution: Rybicki & Lightman,
        ! Radiative Processes in Astrophysics (1979), Eq. (1.30).
        type(ece_transfer_t) :: t
        real(dp) :: reference, thin_error, thick_error, vacuum_error, equilibrium_error
        real(dp) :: tb(1), ab(1), lb(1), rb, wb, dt
        integer :: status
        call thermal_transfer([2100.0_dp], [.7_dp], [1.4_dp], 0.0_dp, 0.0_dp, t, status)
        reference = 2100*(1 - exp(-.98_dp))
        call check(status == 0 .and. abs(t%temperature - reference) < 1.0e-11_dp, &
                   'uniform slab')
        call thermal_transfer([2100.0_dp], [1.0e-12_dp], [1.0_dp], 0.0_dp, 0.0_dp, t, &
            status)
        thin_error = abs(t%temperature/(2100.0_dp*1.0e-12_dp) - 1)
        call check(status == 0 .and. thin_error < 6.0e-13_dp, &
                   'thin emissivity integral alpha Te L')
        call t%jvp([0.0_dp], [1.0_dp], [0.0_dp], 0.0_dp, 0.0_dp, dt)
        call check(abs(dt/2100.0_dp - 1) < 1.1e-12_dp, 'thin absorption tangent')
        call thermal_transfer([2100.0_dp], [100.0_dp], [1.0_dp], .8_dp, 7.0_dp, t, &
            status)
        thick_error = abs(t%temperature - 2100)
        call check(status == 0 .and. thick_error < 1.0e-11_dp, &
            'thick thermal equilibrium')
        call thermal_transfer([2100.0_dp], [0.0_dp], [1.0_dp], .8_dp, 7.0_dp, t, status)
        vacuum_error = abs(t%temperature - 7)
        call check(status == 0 .and. vacuum_error < 1.0e-13_dp, &
                   'vacuum with reflecting thermal walls')
        call thermal_transfer([380.0_dp], [1.0e-12_dp], [1.0_dp], &
                              1 - 1.0e-12_dp, 380.0_dp, t, status)
        equilibrium_error = abs(t%temperature - 380)
        call check(status == 0 .and. equilibrium_error < 1.0e-10_dp, &
                   'nearly closed cavity equilibrium')
        call thermal_transfer([2100.0_dp], [0.0_dp], [1.0_dp], 0.0_dp, 7.0_dp, t, &
            status)
        call t%vjp(1.0_dp, tb, ab, lb, rb, wb)
        call check(abs(ab(1) - 2093) < 1.0e-11_dp, 'zero opacity derivative')
        call check(abs(tb(1)) + abs(lb(1)) + abs(rb) < 1.0e-13_dp, &
            'zero opacity inactive rows')
        call thermal_transfer([real(dp) ::], [real(dp) ::], [real(dp) ::], &
                              .5_dp, 7.0_dp, t, status)
        call check(status == 0 .and. abs(t%temperature - 7) < 1.0e-13_dp, &
            'empty vacuum ray')
        print '(A,ES13.5)', 'thin_relative_error=', thin_error
        print '(A,ES13.5)', 'thick_error_eV=', thick_error
        print '(A,ES13.5)', 'vacuum_error_eV=', vacuum_error
        print '(A,ES13.5)', 'cavity_equilibrium_error_eV=', equilibrium_error
    end subroutine

    subroutine layered_wall()
        type(ece_transfer_t) :: t
        real(dp), parameter :: te(3) = [200.0_dp, 900.0_dp, 500.0_dp], &
                               alpha(3) = [.2_dp, 3.0_dp, .7_dp], length(3) = [.4_dp, &
                                   .2_dp, .3_dp]
        real(dp) :: incoming, near, far, p, rho, tw, error, reference
        integer :: i, k, status
        rho = .82_dp
        tw = .03_dp
        incoming = tw
        ! Independent repeated boundary reflections, rather than the closed form.
        do k = 1, 200
            near = incoming
            do i = 3, 1, -1
                p = exp(-alpha(i)*length(i))
                near = p*near + (1 - p)*te(i)
            end do
            far = (1 - rho)*tw + rho*near
            do i = 1, 3
                p = exp(-alpha(i)*length(i))
                far = p*far + (1 - p)*te(i)
            end do
            incoming = (1 - rho)*tw + rho*far
        end do
        reference = near
        call thermal_transfer(te, alpha, length, rho, tw, t, status)
        error = abs(t%temperature - reference)
        call check(status == 0 .and. error < 2.0e-11_dp, &
                   'multiple reflections of nonuniform source')
        print '(A,ES13.5)', 'wall_series_error_eV=', error
    end subroutine

    subroutine transfer_derivatives()
        type(ece_transfer_t) :: t, plus, minus
        real(dp) :: x(11), dx(11), xp(11), xm(11), gradient(11), jac, fd, h, &
                    max_error, transpose_error
        integer :: i, status
        x = [200.0_dp, 900.0_dp, 500.0_dp, .2_dp, 3.0_dp, .7_dp, .4_dp, .2_dp, .3_dp, &
             .82_dp, .03_dp]
        call thermal_transfer(x(:3), x(4:6), x(7:9), x(10), x(11), t, status)
        call check(status == 0, 'derivative trace')
        max_error = 0.0_dp
        do i = 1, size(x)
            dx = 0.0_dp
            dx(i) = 1.0_dp
            h = 1.0e-5_dp*max(abs(x(i)), .01_dp)
            ! This input is affine; a larger step avoids subtraction noise.
            if (i == 11) h = 1.0e-3_dp
            xp = x + h*dx
            xm = x - h*dx
            call thermal_transfer(xp(:3), xp(4:6), xp(7:9), xp(10), xp(11), plus, &
                status)
            call check(status == 0, 'positive derivative evaluation')
            call thermal_transfer(xm(:3), xm(4:6), xm(7:9), xm(10), xm(11), minus, &
                status)
            call check(status == 0, 'negative derivative evaluation')
            call t%jvp(dx(:3), dx(4:6), dx(7:9), dx(10), dx(11), jac)
            fd = (plus%temperature - minus%temperature)/(2*h)
            max_error = max(max_error, abs(jac - fd)/max(abs(jac), 1.0e-4_dp))
        end do
        call check(max_error < 2.0e-8_dp, &
            'every thermal, opacity, length and wall derivative')
        dx = [.3_dp, -.2_dp, .1_dp, .04_dp, -.02_dp, .05_dp, .01_dp, .02_dp, -.01_dp, &
              .03_dp, .001_dp]
        call t%jvp(dx(:3), dx(4:6), dx(7:9), dx(10), dx(11), jac)
        call t%vjp(.7_dp, gradient(:3), gradient(4:6), gradient(7:9), gradient(10), &
                   gradient(11))
        transpose_error = abs(.7_dp*jac - dot_product(gradient, dx))/max(1.0_dp, &
            abs(.7_dp*jac))
        call check(transpose_error < 1.0e-13_dp, 'transfer transpose identity')
        print '(A,ES13.5)', 'transfer_max_column_fd_relative_error=', max_error
        print '(A,ES13.5)', 'transfer_transpose_relative_error=', transpose_error
    end subroutine

    subroutine convergence()
        type(ece_transfer_t) :: t
        real(dp), allocatable :: te(:), alpha(:), length(:)
        real(dp) :: reference, error, previous
        integer :: i, k, n, status
        ! Te(x)=200+100*x eV, alpha=2/m, 0<x<1m, no incident radiation.
        reference = 200*(1 - exp(-2.0_dp)) + 100*((1 - exp(-2.0_dp))/2 - exp(-2.0_dp))
        previous = huge(1.0_dp)
        do k = 0, 4
            n = 16*4**k
            allocate (te(n), alpha(n), length(n))
            te = [(200 + 100*(i - .5_dp)/n, i=1, n)]
            alpha = 2.0_dp
            length = 1.0_dp/n
            call thermal_transfer(te, alpha, length, 0.0_dp, 0.0_dp, t, status)
            error = abs(t%temperature - reference)
            call check(status == 0 .and. error < previous/10, &
                'nonuniform slab refinement')
            print '(A,I0,A,ES13.5)', 'slab_cells=', n, ',error_eV=', error
            previous = error
            deallocate (te, alpha, length)
        end do
        call check(error < 1.0e-6_dp, 'linear-source analytic integral')
    end subroutine

    subroutine receiver()
        type(ece_response_channel_t) :: channels(2)
        type(ece_receiver_t) :: model
        type(calibration_prior_t) :: prior
        real(dp) :: mean(6), factor(6, 6), covariance(6, 6), eta(6), te(4), dt(4), de(6)
        real(dp) :: y(2), var(2), yp(2), ym(2), vp(2), vm(2), dy(2), dv(2), tb(4), eb(6)
        real(dp) :: yb(2), vb(2), reference, err, varerr, transpose_error, h, &
                    bandwidth(2), tf(2)
        integer :: i, j, c, status
        do c = 1, 2
            channels(c)%rays = reshape([1, 2, 3, 4], [2, 2])
            channels(c)%frequency = [1.0e11_dp, 1.2e11_dp]
            channels(c)%weights = [2.0e9_dp, 1.0e9_dp]
            channels(c)%response = [1.0_dp, .5_dp]
            channels(c)%antenna = [1.0_dp, 3.0_dp]
            channels(c)%duration = .001_dp*c
        end do
        ! Shared ray indices exercise the parallel adjoint reduction.
        mean = [.003_dp, .01_dp, 10.0_dp, .004_dp, .02_dp, 20.0_dp]
        factor = 0.0_dp
        do i = 1, 6
            factor(i, i) = max(abs(mean(i))*.1_dp, .001_dp)
        end do
        factor(2, 1) = .003_dp
        factor(5, 3) = -.002_dp
        covariance = matmul(factor, transpose(factor))
        call prior%init(mean, covariance, status)
        call check(status == 0, 'correlated calibration prior')
        call model%init(channels, 4, prior, status)
        call check(status == 0 .and. model%n_nuisance() == 6, 'receiver init')
        bandwidth = model%effective_bandwidth()
        call check(maxval(abs(bandwidth/((2.5e9_dp)**2/2.25e9_dp) - 1)) < 1.0e-14_dp, &
                   'effective noise bandwidth')
        te = [200.0_dp, 400.0_dp, 800.0_dp, 1200.0_dp]
        eta = 0.0_dp
        call model%evaluate(te, eta, y, var, status)
        tf = [350.0_dp, 1100.0_dp]
        reference = (2*tf(1) + .5_dp*tf(2))/2.5_dp
        call check(status == 0 .and. &
            abs(y(1) - (mean(1)*reference + mean(2))) < 1.0e-14_dp, &
                   'bandwidth and antenna average')
        reference = mean(1)**2* &
                    (2.0e9_dp*(tf(1) + mean(3))**2 + .25e9_dp*(tf(2) + mean(3))**2)/ &
                    ((2.5e9_dp)**2*channels(1)%duration)
        call check(abs(var(1)/reference - 1) < 1.0e-14_dp, &
            'nonuniform spectrum thermal noise')
        call model%evaluate([500.0_dp, 500.0_dp, 500.0_dp, 500.0_dp], eta, y, var, &
            status)
        reference = (mean(1)*(500 + mean(3)))**2/(bandwidth(1)*channels(1)%duration)
        call check(status == 0 .and. abs(var(1)/reference - 1) < 1.0e-14_dp, &
                   'flat spectrum radiometer equation')
        err = 0.0_dp
        varerr = 0.0_dp
        h = 1.0e-4_dp
        do i = 1, 10
            dt = [(merge(100.0_dp, 0.0_dp, i == j), j = 1, 4)]
            de = [(merge(.5_dp, 0.0_dp, i == 4 + j), j = 1, 6)]
            call model%evaluate(te + h*dt, eta + h*de, yp, vp, status)
            call check(status == 0, 'receiver FD positive')
            call model%evaluate(te - h*dt, eta - h*de, ym, vm, status)
            call check(status == 0, 'receiver FD negative')
            call model%jvp(te, eta, dt, de, dy, dv)
            err = max(err, maxval(abs(dy - (yp - ym)/(2*h)))/max(maxval(abs(dy)), &
                1.0e-12_dp))
            varerr = max(varerr, &
                maxval(abs(dv - (vp - vm)/(2*h)))/max(maxval(abs(dv)), 1.0e-15_dp))
        end do
        call check(err < 2.0e-8_dp .and. varerr < 2.0e-8_dp, &
                   'all receiver and calibration columns')
        dt = [20.0_dp, -30.0_dp, 70.0_dp, -40.0_dp]
        de = [.1_dp, -.2_dp, .3_dp, -.2_dp, .1_dp, .2_dp]
        yb = [.7_dp, -.2_dp]
        vb = [1.0e6_dp, -3.0e6_dp]
        call model%jvp(te, eta, dt, de, dy, dv)
        call model%vjp(te, eta, yb, vb, tb, eb)
        reference = dot_product(yb, dy) + dot_product(vb, dv)
        transpose_error = abs(reference - dot_product(tb, dt) - &
            dot_product(eb, de))/max(abs(reference), 1.0_dp)
        call check(transpose_error < 1.0e-13_dp, &
                   'receiver transpose with shared rays and correlated prior')
        print '(A,ES13.5)', 'receiver_voltage_max_column_fd_relative_error=', err
        print '(A,ES13.5)', 'receiver_variance_max_column_fd_relative_error=', varerr
        print '(A,ES13.5)', 'receiver_transpose_relative_error=', transpose_error
        te = 0.0_dp
        call model%evaluate(te, eta, y, var, status)
        call check(status == 0 .and. maxval(abs(y - mean([2, 5]))) < 1.0e-14_dp, &
                   'receiver zero source')
        eta(1) = -100.0_dp
        call model%evaluate(te, eta, y, var, status)
        call check(status /= 0, 'reject nonphysical calibration latent')
        channels(1)%weights = 0.0_dp
        call model%init(channels, 4, prior, status)
        call check(status /= 0, 'reject zero response area')
    end subroutine

    subroutine noise_kernel()
        ! Condon & Ransom, Essential Radio Astronomy (2016), Eq. (3.154):
        ! sigma_T = Tsys/sqrt(B*t), one total-power polarization.
        ! https://www.cv.nrao.edu/~sransom/web/Ch3.html
        real(dp) :: x(5), xp(5), xm(5), value, row(5), vp, vm, scratch(5), h, error
        integer :: i
        x = [500.0_dp, 10.0_dp, .003_dp, 2.0e9_dp, .001_dp]
        call ece_radiometer_noise(x(1), x(2), x(3), x(4), x(5), value, row)
        call check(abs(value/((.003_dp*510)**2/(2.0e9_dp*.001_dp)) - 1) < 1.0e-14_dp, &
                   'published ideal radiometer noise')
        error = 0.0_dp
        do i = 1, 5
            h = 1.0e-5_dp*x(i)
            xp = x
            xm = x
            xp(i) = xp(i) + h
            xm(i) = xm(i) - h
            call ece_radiometer_noise(xp(1), xp(2), xp(3), xp(4), xp(5), vp, scratch)
            call ece_radiometer_noise(xm(1), xm(2), xm(3), xm(4), xm(5), vm, scratch)
            error = max(error, abs((vp - vm)/(2*h)/row(i) - 1))
        end do
        call check(error < 2.0e-9_dp, 'noise kernel all derivative rows')
        print '(A,ES13.5)', 'noise_max_column_fd_relative_error=', error
    end subroutine

    subroutine composed()
        type(ece_response_channel_t) :: channel(1)
        type(ece_receiver_t) :: model
        type(calibration_prior_t) :: prior
        type(ece_transfer_t) :: trace(2), perturbed(2)
        real(dp) :: te(2, 2), alpha(2, 2), length(2, 2), dte(2, 2), dalpha(2, 2), &
            dl(2, 2)
        real(dp) :: rho(2), tw(2), drho(2), dtw(2), eta(3), deta(3), cov(3, 3), h
        real(dp) :: t(2), dt(2), tb(2), eb(3), etb(2), ab(2), lb(2), rb, wb
        real(dp) :: y(1), var(1), yp(1), vp(1), ym(1), vm(1), dy(1), dv(1), error, &
            lhs, rhs
        integer :: status, codes(2), k
        channel(1)%rays = reshape([1, 2], [1, 2])
        channel(1)%frequency = [1.0e11_dp, 1.2e11_dp]
        channel(1)%weights = [1.0e9_dp, 1.0e9_dp]
        channel(1)%response = [1.0_dp, .7_dp]
        channel(1)%antenna = [1.0_dp]
        channel(1)%duration = .001_dp
        cov = 0.0_dp
        cov(1, 1) = 1.0e-6_dp
        cov(2, 2) = .01_dp
        cov(3, 3) = 1.0_dp
        call prior%init([.003_dp, .01_dp, 10.0_dp], cov, status)
        call check(status == 0, 'composed prior')
        call model%init(channel, 2, prior, status)
        call check(status == 0, 'composed receiver')
        te = reshape([200.0_dp, 800.0_dp, 700.0_dp, 300.0_dp], [2, 2])
        alpha = reshape([.7_dp, .2_dp, .8_dp, .4_dp], [2, 2])
        length = .3_dp
        rho = .4_dp
        tw = .03_dp
        dte = reshape([30.0_dp, -70.0_dp, 20.0_dp, 50.0_dp], [2, 2])
        dalpha = .03_dp
        dl = .02_dp
        drho = .05_dp
        dtw = .002_dp
        eta = 0.0_dp
        deta = [.1_dp, -.2_dp, .3_dp]
        call thermal_transfer_batch(te, alpha, length, rho, tw, trace, codes)
        call check(all(codes == 0), 'composed transfer')
        t = trace%temperature
        do k = 1, 2
            call trace(k)%jvp(dte(:, k), dalpha(:, k), dl(:, k), drho(k), dtw(k), dt(k))
        end do
        call model%evaluate(t, eta, y, var, status)
        call check(status == 0, 'composed voltage')
        call model%jvp(t, eta, dt, deta, dy, dv)
        h = 1.0e-4_dp
        call thermal_transfer_batch(te + h*dte, alpha + h*dalpha, length + h*dl, &
                                    rho + h*drho, tw + h*dtw, perturbed, codes)
        call check(all(codes == 0), 'composed FD positive transfer')
        call model%evaluate(perturbed%temperature, eta + h*deta, yp, vp, status)
        call check(status == 0, 'composed FD positive voltage')
        call thermal_transfer_batch(te - h*dte, alpha - h*dalpha, length - h*dl, &
                                    rho - h*drho, tw - h*dtw, perturbed, codes)
        call check(all(codes == 0), 'composed FD negative transfer')
        call model%evaluate(perturbed%temperature, eta - h*deta, ym, vm, status)
        call check(status == 0, 'composed FD negative voltage')
        error = max(abs((yp(1) - ym(1))/(2*h)/dy(1) - 1), &
            abs((vp(1) - vm(1))/(2*h)/dv(1) - 1))
        call check(error < 2.0e-8_dp, &
            'complete prescribed-ray to voltage/noise tangent')
        call model%vjp(t, eta, [.7_dp], [1.0e6_dp], tb, eb)
        lhs = .7_dp*dy(1) + 1.0e6_dp*dv(1)
        rhs = dot_product(eb, deta)
        do k = 1, 2
            call trace(k)%vjp(tb(k), etb, ab, lb, rb, wb)
            rhs = rhs + dot_product(etb, dte(:, k)) + dot_product(ab, dalpha(:, k)) + &
                  dot_product(lb, dl(:, k)) + rb*drho(k) + wb*dtw(k)
        end do
        call check(abs(lhs - rhs) < 1.0e-12_dp*max(1.0_dp, abs(lhs)), &
            'composed transpose')
        print '(A,ES13.5)', 'composed_voltage_noise_fd_relative_error=', error
    end subroutine

    subroutine invalid()
        type(ece_transfer_t) :: t, batch(2)
        integer :: status, codes(2)
        real(dp) :: nan
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        call thermal_transfer([1.0_dp], [-1.0_dp], [1.0_dp], 0.0_dp, 0.0_dp, t, status)
        call check(status /= 0 .and. .not. t%valid, 'reject negative opacity')
        call thermal_transfer([1.0_dp], [1.0_dp], [1.0_dp], 1.0_dp, 0.0_dp, t, status)
        call check(status /= 0 .and. .not. t%valid, 'reject perfectly closed wall')
        call thermal_transfer([nan], [1.0_dp], [1.0_dp], 0.0_dp, 0.0_dp, t, status)
        call check(status /= 0 .and. .not. t%valid, 'reject NaN')
        call thermal_transfer([1.0_dp], [real(dp) ::], [1.0_dp], 0.0_dp, 0.0_dp, t, &
            status)
        call check(status /= 0 .and. .not. t%valid, 'reject input shape')
        call thermal_transfer_batch(reshape([100.0_dp, 100.0_dp], [1, 2]), &
                                    reshape([1.0_dp, -1.0_dp], [1, 2]), &
                                        reshape([1.0_dp, 1.0_dp], [1, 2]), &
                                    [0.0_dp, 0.0_dp], [0.0_dp, 0.0_dp], batch, codes)
        call check(codes(1) == 0 .and. codes(2) /= 0, 'independent batch statuses')
        call check(batch(1)%valid .and. .not. batch(2)%valid, &
            'independent batch validity')
    end subroutine
end program
