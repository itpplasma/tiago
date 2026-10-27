program test_thomson
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit, int64
    use fortnum_quadrature, only: gauss_legendre_ab
    use tiago_thomson_spectrum, only: thomson_spectrum_t
    use tiago_thomson_scattering, only: thomson_volume_t
    use tiago_thomson_polychromator, only: recorded_thomson_t, polychromator_channel_t
    use tiago_calibration, only: calibration_prior_t
    use tiago_profile_diagnostic, only: plasma_sample_t, N_QUANTITIES, Q_NE, Q_TE, &
                                        check_diagnostic_derivatives
    use tiago_thomson_reduction, only: pulse_integral, thomson_noise
    use tiago_signal_reduction, only: linear_reduction_t
    use tiago_signal_source, only: signal_window_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), lambda0 = 532.0e-9_dp
    integer :: fails = 0
    character(32) :: argument
    call get_command_argument(1, argument)
    if (trim(argument) == 'benchmark') then
        call benchmark()
        stop
    end if
    if (trim(argument) == 'spectra') then
        call spectrum_csv()
        stop
    end if
    call spectra()
    call instrument()
    call reduction()
    call posterior()
    if (fails > 0) error stop 1
    print '(A)', 'Thomson recorded signals: PASS'
contains
    subroutine check(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        write (error_unit, '(A)') 'FAIL: '//message
        fails = fails + 1
    end subroutine
    subroutine spectra()
        type(thomson_spectrum_t) :: model
        ! Independent velocity-space integration of Mirnov & Den Hartog,
        ! PPCF 59 (2017) 063001, Eq. (60), DOI 10.1088/1361-6587/aa61e2.
        ! Integrate the delta function analytically; 96 Laguerre gamma nodes
        ! and 512 periodic azimuth nodes. These are calculated reference values,
        ! not experimental data or a claim of published tabulated spectra.
        real(dp), parameter :: temp(7) = [10.0_dp, 1000.0_dp, 1000.0_dp, &
                                          20000.0_dp, 20000.0_dp, 20000.0_dp, 20000.0_dp]
        real(dp), parameter :: ratio(7) = [1.0_dp, .9_dp, 1.1_dp, .6_dp, .8_dp, 1.0_dp, 1.4_dp]
        real(dp), parameter :: angle(7) = [pi/2, pi/3, pi/2, pi, pi/2, pi/2, pi]
        real(dp), parameter :: reference(7) = [63.763498826447091_dp, .68417375672270264_dp, &
                                               1.5597409983497064_dp, 1.2348698846135278_dp, 1.5502001045355356_dp, &
                                               1.2340140027862891_dp, .24275169448784631_dp]
        real(dp) :: s(1), ds(1), sp(1), sm(1), wave(1), err, fd_error, h, area
        real(dp) :: z(800), w(800), l(800), values(800), derivs(800), integral, integral_error
        real(dp) :: sigma, gaussian, cold_error, second, center, normalized_fd, tt
        real(dp), parameter :: test_angle(3) = [.1_dp, pi/2, pi]
        real(dp), parameter :: norm_temp(3) = [10.0_dp, 1000.0_dp, 20000.0_dp]
        integer :: i, j
        call model%init()
        err = 0.0_dp
        fd_error = 0.0_dp
        do i = 1, 7
            wave = ratio(i)*lambda0
            call model%cross_section(temp(i), lambda0, angle(i), wave, s, ds)
            err = max(err, abs(s(1)*lambda0/reference(i) - 1))
            if (temp(i) <= 10.0_dp .or. temp(i) >= 20000.0_dp) cycle
            h = temp(i)*1.0e-4_dp
            ! Separate buffers: the derivative output must not alias a value.
            call raw_value(temp(i) + h, angle(i), wave, sp)
            call raw_value(temp(i) - h, angle(i), wave, sm)
            fd_error = max(fd_error, abs((sp(1) - sm(1))/(2*h) - ds(1))/max(abs(ds(1)), tiny(1.0_dp)))
        end do
        call check(err < 2.0e-11_dp, 'independent relativistic velocity-space reference')
        call check(fd_error < 2.0e-7_dp, 'spectrum temperature derivative')
        integral_error = 0.0_dp
        do i = 1, 3
            do j = 1, 3
                ! A different wavelength rule checks the composite normalization.
                sigma = sqrt(2*norm_temp(i)*(1 - cos(test_angle(j)))/510998.95069_dp)
                call gauss_legendre_ab(800, -12*sigma, 12*sigma, z, w)
                l = lambda0*exp(z)
                call model%normalized(norm_temp(i), lambda0, test_angle(j), l, values, derivs, area)
                integral = dot_product(w*l, values)
                integral_error = max(integral_error, abs(integral - 1))
                call check(abs(dot_product(w*l, derivs)) < 1.0e-10_dp, 'normalized derivative integral zero')
            end do
        end do
        normalized_fd = 0.0_dp
        do i = 1, 3
            tt = max(10.01_dp, min(19990.0_dp, norm_temp(i)))
            h = tt*1.0e-4_dp
            wave = 1.005_dp*lambda0
            call model%normalized(tt, lambda0, pi/2, wave, s, ds)
            call normalized_value(tt + h, wave, sp)
            call normalized_value(tt - h, wave, sm)
            normalized_fd = max(normalized_fd, abs((sp(1) - sm(1))/(2*h)/ds(1) - 1))
        end do
        call check(normalized_fd < 2.0e-7_dp, 'normalized spectrum finite difference')
        call check(integral_error < 2.0e-10_dp, 'unit photon spectrum normalization 10 eV--20 keV')
        sigma = lambda0*sqrt(2*10.0_dp/510998.95069_dp)
        wave = lambda0
        call model%normalized(10.0_dp, lambda0, pi/2, wave, s, ds)
        gaussian = 1/(sqrt(2*pi)*sigma)
        cold_error = abs(s(1)/gaussian - 1)
        call check(cold_error < 2.0e-4_dp, '10 eV nonrelativistic Gaussian peak')
        call gauss_legendre_ab(800, lambda0 - 12*sigma, lambda0 + 12*sigma, l, w)
        call model%normalized(10.0_dp, lambda0, pi/2, l, values, derivs)
        center = dot_product(w, values*l)
        second = dot_product(w, values*(l - center)**2)
        call check(abs(second/sigma**2 - 1) < 3.0e-4_dp, '10 eV Gaussian variance')
        print '(A,ES12.4)', 'relativistic reference max relative error: ', err
        print '(A,ES12.4)', 'spectrum Te finite difference relative error: ', fd_error
        print '(A,ES12.4)', 'normalization max absolute error: ', integral_error
        print '(A,ES12.4)', 'normalized spectrum FD relative error: ', normalized_fd
        print '(A,ES12.4)', '10 eV Gaussian peak relative error: ', cold_error
        print '(A,ES12.4)', '10 eV Gaussian variance relative error: ', abs(second/sigma**2 - 1)
    end subroutine

    subroutine normalized_value(te, wave, s)
        real(dp), intent(in) :: te, wave(:)
        real(dp), intent(out) :: s(:)
        type(thomson_spectrum_t) :: model
        real(dp) :: ds(size(s))
        call model%init()
        call model%normalized(te, lambda0, pi/2, wave, s, ds)
    end subroutine
    subroutine spectrum_csv()
        type(thomson_spectrum_t) :: model
        real(dp) :: wave(501), s(501), ds(501), t(3)
        integer :: i, j
        call model%init()
        t = [10.0_dp, 1000.0_dp, 20000.0_dp]
        print '(A)', 'temperature_eV,wavelength_nm,density_per_nm'
        do j = 1, 3
            if (j == 1) then
                wave = [(lambda0*(.97_dp + .06_dp*(i - 1)/500), i=1, 501)]
            else
                wave = [(lambda0*(.4_dp + 1.2_dp*(i - 1)/500), i=1, 501)]
            end if
            call model%normalized(t(j), lambda0, pi/2, wave, s, ds)
            do i = 1, 501
                print '(F10.2,A,F12.5,A,ES20.12)', t(j), ',', wave(i)*1e9_dp, ',', s(i)*1e-9_dp
            end do
        end do
    end subroutine

    subroutine raw_value(te, theta, wave, s)
        real(dp), intent(in) :: te, theta, wave(:)
        real(dp), intent(out) :: s(:)
        type(thomson_spectrum_t) :: model
        real(dp) :: ds(size(s))
        call model%init()
        call model%cross_section(te, lambda0, theta, wave, s, ds)
    end subroutine

    subroutine setup(model, sample, n, rayleigh, aligned)
        type(recorded_thomson_t), intent(out) :: model
        type(plasma_sample_t), intent(out) :: sample
        integer, intent(in) :: n
        logical, optional, intent(in) :: rayleigh, aligned
        type(thomson_volume_t) :: volumes(1)
        type(polychromator_channel_t) :: channels(n)
        type(calibration_prior_t) :: prior, gas
        real(dp) :: mean(3 + 2*n), cov(3 + 2*n, 3 + 2*n), gm(3), gc(3, 3), sigma(3 + 2*n)
        real(dp) :: x(96), w(96), lo, hi
        integer :: i, j, status, band
        logical :: use_gas, align
        use_gas = .false.
        align = .true.
        if (present(rayleigh)) use_gas = rayleigh
        if (present(aligned)) align = aligned
        volumes(1)%label = 'synthetic volume'
        volumes(1)%center = [2.0_dp, 0.0_dp, 0.0_dp]
        volumes(1)%length = .02_dp
        volumes(1)%npoints = 3
        mean = 1.0_dp
        mean(3) = .8_dp
        sigma = .02_dp
        sigma(1) = .25_dp
        sigma(2) = .03_dp
        do i = 1, n
            band = mod(i - 1, 3) + 1
            select case (band)
            case (1); lo = .75_dp; hi = .88_dp
            case (2); lo = .88_dp; hi = .96_dp
            case default; lo = .96_dp; hi = 1.10_dp
            end select
            call gauss_legendre_ab(96, lo*lambda0, hi*lambda0, x, w)
            channels(i)%name = 'synthetic band'
            channels(i)%wavelength = x
            channels(i)%weights = w
            allocate (channels(i)%transmission(96), channels(i)%qe(96))
            channels(i)%transmission = 1.0_dp
            channels(i)%qe = .6_dp
            channels(i)%lambda0 = lambda0
            channels(i)%theta = pi/2
            channels(i)%energy = .5_dp
            channels(i)%solid_angle = .01_dp
            if (band == 3) channels(i)%elastic_response = .6_dp
            j = 4 + 2*(i - 1)
            mean(j + 1) = 1.0e-10_dp
            sigma(j + 1) = .02e-10_dp
        end do
        cov = 0.0_dp
        do i = 1, size(mean)
            cov(i, i) = sigma(i)**2
        end do
        call prior%init(mean, cov, status)
        if (status /= 0) error stop 'synthetic calibration prior'
        gm = [100.0_dp, 300.0_dp, 1.64e-30_dp]
        gc = 0.0_dp
        gc(1, 1) = 1.0_dp
        gc(2, 2) = 1.0_dp
        gc(3, 3) = (.01e-30_dp)**2
        call gas%init(gm, gc, status)
        if (status /= 0) error stop 'synthetic gas prior'
        if (use_gas) then
            call model%init(volumes, channels, prior, status, align, gas)
        else
            call model%init(volumes, channels, prior, status, align)
        end if
        if (status /= 0) error stop 'synthetic polychromator setup'
        sample%n = size(model%sampling%w)
        allocate (sample%q(N_QUANTITIES, sample%n))
        sample%q = 0.0_dp
        sample%q(Q_NE, :) = 2.0e19_dp
        sample%q(Q_TE, :) = 1000.0_dp
    end subroutine
    subroutine instrument()
        type(recorded_thomson_t) :: model, gas_model
        type(plasma_sample_t) :: sample, gas_sample
        real(dp), allocatable :: eta(:), x(:, :), dx(:, :, :)
        real(dp) :: y(3), pe(3), yp(3), ym(3), err, adj, ray_err, ray_adj, ratio, lo, hi, mid
        real(dp) :: truth(3), relerr, qd(N_QUANTITIES, 3), ed(11), dy(3), h
        integer :: i
        call setup(model, sample, 3, aligned=.false.)
        allocate (eta(model%n_nuisance()))
        eta = .15_dp
        call model%sample_points(eta, x, dx)
        call check(abs(x(1, 1) - 2.15_dp) < 1.0e-14_dp, 'existing volume alignment')
        call check(all(abs(dx(:, :, 3:)) <= 0.0_dp), 'calibration fixes sample points')
        call model%predict(sample, eta, y)
        call model%photoelectrons(sample, eta, pe)
        call photon_checks(model, sample, eta)
        call check(maxval(abs(y/pe - (1.0e-10_dp + .02e-10_dp*.15_dp))) < 1.0e-24_dp, &
                   'photoelectron to V s gain')
        call check_diagnostic_derivatives(model, sample, eta, 1.0e-5_dp, err, adj)
        call check(err < 2.0e-8_dp, 'recorded Thomson FD in ne Te and calibration')
        call check(adj < 2.0e-12_dp, 'recorded Thomson adjoint')
        ! Every calibration column is exercised individually, including gain,
        ! whose dimensional scale is too small for a generic random direction.
        qd = 0.0_dp
        relerr = 0.0_dp
        h = 1.0e-4_dp
        do i = 3, size(eta)
            ed = 0.0_dp
            ed(i) = 1.0_dp
            call model%jvp(sample, eta, qd, ed, dy)
            eta(i) = eta(i) + h
            call model%predict(sample, eta, yp)
            eta(i) = eta(i) - 2*h
            call model%predict(sample, eta, ym)
            eta(i) = eta(i) + h
            relerr = max(relerr, maxval(abs((yp - ym)/(2*h) - dy))/max(maxval(abs(dy)), tiny(1.0_dp)))
        end do
        call check(relerr < 2.0e-9_dp, 'each physical calibration chain FD')
        ! Ratio removes density and a shared absolute calibration.
        eta = 0.0_dp
        call model%predict(sample, eta, truth)
        ratio = truth(1)/truth(2)
        lo = 100.0_dp
        hi = 5000.0_dp
        do i = 1, 45
            mid = (lo + hi)/2
            sample%q(Q_TE, :) = mid
            call model%predict(sample, eta, y)
            if (y(1)/y(2) < ratio) then
                lo = mid
            else
                hi = mid
            end if
        end do
        call check(abs(mid/1000.0_dp - 1) < 1.0e-10_dp, 'ideal rectangular filter ratio recovers Te')
        print '(A,ES12.4)', 'ratio recovered Te relative error: ', abs(mid/1000.0_dp - 1)
        sample%q(Q_NE, :) = 0.0_dp
        sample%q(Q_TE, :) = 0.0_dp
        call model%predict(sample, eta, y)
        call check(all(abs(y) <= 0.0_dp), 'vacuum gives zero scattered signal')
        sample%q(Q_NE, :) = 2.0e19_dp
        sample%q(Q_TE, :) = [400.0_dp, 1000.0_dp, 1800.0_dp]
        call check_diagnostic_derivatives(model, sample, eta, 1.0e-5_dp, ray_err, ray_adj)
        call check(ray_err < 2.0e-8_dp .and. ray_adj < 2.0e-12_dp, 'nonuniform volume spectrum derivative')
        call setup(gas_model, gas_sample, 3, rayleigh=.true.)
        deallocate (eta)
        allocate (eta(gas_model%n_nuisance()))
        eta = .1_dp
        call check_diagnostic_derivatives(gas_model, gas_sample, eta, 1.0e-5_dp, ray_err, ray_adj)
        call check(ray_err < 2.0e-8_dp, 'Rayleigh calibration gas and optical FD')
        call check(ray_adj < 2.0e-12_dp, 'Rayleigh adjoint')
        call photon_checks(gas_model, gas_sample, eta)
        call gas_model%predict(gas_sample, eta, y)
        call check(all(abs(y(:2)) <= 0.0_dp) .and. y(3) > 0.0_dp, 'elastic line respects filter transmission')
        ! Synthetic polarizability with an independent classical dipole law.
        ! Classical dipole law: Sneep & Ubachs, JQSRT 92 (2005) 293--310,
        ! DOI 10.1016/j.jqsrt.2004.07.025; exact supplied alpha is synthetic here.
        eta = 0.0_dp
        call gas_model%predict(gas_sample, eta, y)
        mid = 100.0_dp/(1.380649e-23_dp*300.0_dp)*16*pi**4*(1.64e-30_dp)**2/lambda0**4 &
              *.5_dp*lambda0/(6.62607015e-34_dp*299792458.0_dp)*.01_dp*.02_dp*.6_dp*.8_dp*1.0e-10_dp
        call check(abs(y(3)/mid - 1) < 3.0e-14_dp, 'Rayleigh independent dipole cross-section')
        print '(A,ES12.4)', 'Thomson JVP finite difference relative error: ', err
        print '(A,ES12.4)', 'Thomson adjoint relative error: ', adj
        print '(A,ES12.4)', 'individual calibration FD relative error: ', relerr
        print '(A,ES12.4)', 'Rayleigh JVP finite difference relative error: ', ray_err
        print '(A,ES12.4)', 'Rayleigh adjoint relative error: ', ray_adj
        print '(A,ES12.4)', 'Rayleigh dipole relative error: ', abs(y(3)/mid - 1)
        ! Sneep & Ubachs (2005), Table 2: Ar refractivity 268e-6 at
        ! density 2.547e25 m^-3, wavelength 532.2 nm. Lorentz--Lorenz gives
        ! alpha_v; compare total cross section against their n-based value
        ! 4.56e-31 m^2 and measured (4.45 +/- 0.30)e-31 m^2 (one sigma).
        gas_model%gas%mean(3) = 3/(4*pi*2.547e25_dp) &
                                *((1 + 268e-6_dp)**2 - 1)/((1 + 268e-6_dp)**2 + 2)
        gas_model%channels(3)%lambda0 = 532.2e-9_dp
        call gas_model%predict(gas_sample, eta, y)
        mid = y(3)/(100.0_dp/(1.380649e-23_dp*300.0_dp) &
                    *.5_dp*532.2e-9_dp/(6.62607015e-34_dp*299792458.0_dp) &
                    *.01_dp*.02_dp*.6_dp*.8_dp*1.0e-10_dp)*8*pi/3
        call check(abs(mid/4.56e-31_dp - 1) < .003_dp, 'published argon n-based cross section')
        call check(abs(mid - 4.45e-31_dp) < .30e-31_dp, 'published argon measured uncertainty')
        print '(A,ES12.4)', 'published Rayleigh n-based relative error: ', abs(mid/4.56e-31_dp - 1)
        print '(A,ES12.4)', 'published Rayleigh measured sigma discrepancy: ', abs(mid - 4.45e-31_dp)/.30e-31_dp
    end subroutine
    subroutine photon_checks(model, sample, eta)
        type(recorded_thomson_t), intent(in) :: model
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        type(plasma_sample_t) :: shifted
        real(dp) :: ed(size(eta)), eb(size(eta)), qd(N_QUANTITIES, sample%n), qb(N_QUANTITIES, sample%n)
        real(dp) :: yp(3), ym(3), dy(3), seed(3), lhs, rhs, h
        ed = .2_dp
        qd = 0.0_dp
        qd(Q_NE, :) = sample%q(Q_NE, :)*.03_dp
        qd(Q_TE, :) = sample%q(Q_TE, :)*.02_dp
        h = 1e-5_dp
        call model%photon_jvp(sample, eta, qd, ed, dy)
        shifted = sample
        shifted%q = sample%q + h*qd
        call model%photoelectrons(shifted, eta + h*ed, yp)
        shifted%q = sample%q - h*qd
        call model%photoelectrons(shifted, eta - h*ed, ym)
        call check(maxval(abs((yp - ym)/(2*h) - dy))/maxval(abs(dy)) < 2e-8_dp, 'photoelectron JVP FD')
        seed = [.3_dp, -.7_dp, .8_dp]
        call model%photon_vjp(sample, eta, seed, qb, eb)
        lhs = dot_product(seed, dy)
        rhs = sum(qb*qd) + dot_product(eb, ed)
        call check(abs(lhs - rhs)/max(abs(lhs), abs(rhs)) < 2e-12_dp, 'photoelectron VJP adjoint')
    end subroutine

    subroutine reduction()
        type(linear_reduction_t) :: step
        type(signal_window_t) :: input, output, shifted
        real(dp) :: t(7), variance, row(6), d(7), jv(1), bar(7), weights(7), expected
        real(dp) :: p(6), pp(6), pm(6), vp, vm, r(6), error
        integer :: status, i
        t = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]*1.0e-9_dp
        call pulse_integral(step, t, [4, 7], [1, 3], status)
        call check(status == 0, 'pulse reduction construction')
        input%time = t
        input%samples = [2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 2.0_dp]
        input%variance = spread(.04_dp, 1, 7)
        allocate (input%covariance(7, 7))
        input%covariance = .01_dp
        input%units = 'V'
        input%provenance%source_id = 'synthetic digitiser'
        input%provenance%version = '1'
        input%provenance%record_id = 'synthetic pulse'
        input%provenance%request_key = 'synthetic trace'
        input%provenance%metadata = 'analytic triangle plus constant baseline'
        call input%seal()
        call step%apply(input, output, status)
        call check(status == 0, 'recorded window reduction')
        call check(abs(output%samples(1) - 3.0e-9_dp) < 1.0e-23_dp, 'analytic pulse area after baseline')
        weights = [-.75_dp, -1.5_dp, -.75_dp, .5_dp, 1.0_dp, 1.0_dp, .5_dp]*1.0e-9_dp
        expected = .04_dp*sum(weights**2)
        call check(abs(output%covariance(1, 1)/expected - 1) < 2.0e-14_dp, 'correlated baseline covariance')
        d = [(sin(real(i, dp)), i=1, 7)]
        call step%jvp(input, d, jv)
        call step%vjp(input, [.7_dp], bar)
        call check(abs(.7_dp*jv(1) - dot_product(bar, d)) < 1.0e-23_dp, 'pulse adjoint')
        shifted = input
        shifted%samples = input%samples + 1.0e-4_dp*d
        call shifted%seal()
        call step%apply(shifted, output, status)
        call check(abs((output%samples(1) - 3.0e-9_dp)/1.0e-4_dp - jv(1)) < 1.0e-19_dp, 'pulse FD')
        p = [100.0_dp, 20.0_dp, 2.0_dp, 3.0_dp, 1.2_dp, .5_dp]
        call thomson_noise(p(1), p(2), p(3), p(4), p(5), p(6), variance, row)
        call check(abs(variance - 627.0_dp) < 1.0e-12_dp, 'Poisson background digitiser noise')
        error = 0.0_dp
        do i = 1, 6
            pp = p; pm = p
            pp(i) = pp(i) + 1.0e-4_dp
            pm(i) = pm(i) - 1.0e-4_dp
            call thomson_noise(pp(1), pp(2), pp(3), pp(4), pp(5), pp(6), vp, r)
            call thomson_noise(pm(1), pm(2), pm(3), pm(4), pm(5), pm(6), vm, r)
            error = max(error, abs((vp - vm)/2.0e-4_dp - row(i))/max(1.0_dp, abs(row(i))))
        end do
        call check(error < 1.0e-8_dp, 'noise derivative rows')
        print '(A,ES12.4)', 'pulse covariance relative error: ', abs(output%covariance(1, 1)/expected - 1)
        print '(A,ES12.4)', 'noise derivative relative error: ', error
    end subroutine

    subroutine posterior()
        type(recorded_thomson_t) :: model, gas_model
        type(plasma_sample_t) :: sample, gas_sample
        integer, parameter :: n = 14, m = 18
        real(dp) :: truth(n), theta(n), data(6), sig(6), res(m), jac(m, n), hess(n, n), step(n)
        real(dp) :: cov(n, n), sigma(n), pull(n)
        integer :: i, it, info
        call setup(model, sample, 3)
        call setup(gas_model, gas_sample, 3, rayleigh=.true.)
        ! Synthetic spectral calibration records narrow the relative response.
        ! The absolute scale starts biased; Rayleigh observations constrain it.
        model%calibration%mean(4:8:2) = [1.10_dp, .91_dp, 1.09_dp]
        gas_model%calibration = model%calibration
        truth = 0.0_dp
        truth(1) = log(2.2_dp)
        truth(2) = log(1.4_dp)
        truth(3) = (1.18_dp - 1)/.25_dp
        truth(6) = (1.12_dp - 1.10_dp)/.02_dp
        truth(8) = (.90_dp - .91_dp)/.02_dp
        truth(10) = (1.08_dp - 1.09_dp)/.02_dp
        call posterior_predictions(model, gas_model, sample, gas_sample, truth, data)
        ! Deterministic half-sigma perturbations, not a coverage claim.
        sig = max(abs(data)*.01_dp, 1.0e-14_dp)
        data = data + sig*[.3_dp, -.4_dp, .2_dp, 0.0_dp, 0.0_dp, -.3_dp]
        theta = 0.0_dp
        theta(1) = log(1.8_dp)
        theta(2) = log(.9_dp)
        do it = 1, 15
            call posterior_residual(model, gas_model, sample, gas_sample, theta, data, sig, res, jac)
            hess = matmul(transpose(jac), jac)
            step = -matmul(transpose(jac), res)
            call dposv('U', n, 1, hess, n, step, n, info)
            if (info /= 0) error stop 'synthetic posterior normal solve'
            theta = theta + step
            if (maxval(abs(step)) < 1.0e-8_dp) exit
        end do
        call posterior_residual(model, gas_model, sample, gas_sample, theta, data, sig, res, jac)
        cov = matmul(transpose(jac), jac)
        call dpotrf('U', n, cov, n, info)
        if (info /= 0) error stop 'synthetic posterior factor'
        call dpotri('U', n, cov, n, info)
        if (info /= 0) error stop 'synthetic posterior covariance'
        do i = 1, n
            sigma(i) = sqrt(cov(i, i))
        end do
        pull = abs(theta - truth)/sigma
        call check(maxval(pull) < 2.0_dp, 'ne Te absolute spectral and gas truth within 2 posterior sd')
        call check(sum(res**2) < 8.0_dp, 'synthetic posterior residual and priors')
        print '(A,I0)', 'synthetic posterior Gauss-Newton iterations: ', it
        print '(A,ES12.4)', 'synthetic posterior max marginal truth pull: ', maxval(pull)
        print '(A,ES12.4)', 'synthetic posterior chi squared including priors: ', sum(res**2)
        print '(A,2ES16.7)', 'recovered ne [m^-3], Te [eV]: ', 1.0e19_dp*exp(theta(1)), 1000*exp(theta(2))
        print '(A,3ES16.7)', 'absolute/spectral1/spectral2 posterior pulls: ', pull(3), pull(6), pull(8)
    end subroutine

    subroutine posterior_predictions(model, gas_model, sample, gas_sample, theta, y)
        type(recorded_thomson_t), intent(in) :: model, gas_model
        type(plasma_sample_t), intent(in) :: sample, gas_sample
        real(dp), intent(in) :: theta(14)
        real(dp), intent(out) :: y(6)
        type(plasma_sample_t) :: smp
        smp = sample
        smp%q(Q_NE, :) = 1.0e19_dp*exp(theta(1))
        smp%q(Q_TE, :) = 1000.0_dp*exp(theta(2))
        call model%predict(smp, theta(3:11), y(:3))
        call gas_model%predict(gas_sample, theta(3:14), y(4:))
    end subroutine
    subroutine posterior_residual(model, gas_model, sample, gas_sample, theta, data, sig, res, jac)
        type(recorded_thomson_t), intent(in) :: model, gas_model
        type(plasma_sample_t), intent(in) :: sample, gas_sample
        real(dp), intent(in) :: theta(14), data(6), sig(6)
        real(dp), intent(out) :: res(18), jac(18, 14)
        type(plasma_sample_t) :: smp
        real(dp) :: y(6), qd(N_QUANTITIES, sample%n), zero(N_QUANTITIES, sample%n)
        real(dp) :: ed(12), dy(6)
        integer :: j, k
        smp = sample
        smp%q(Q_NE, :) = 1.0e19_dp*exp(theta(1))
        smp%q(Q_TE, :) = 1000.0_dp*exp(theta(2))
        call posterior_predictions(model, gas_model, sample, gas_sample, theta, y)
        res(:6) = (y - data)/sig
        res(7:) = theta(3:)
        jac = 0.0_dp
        zero = 0.0_dp
        do j = 1, 14
            qd = 0.0_dp
            ed = 0.0_dp
            if (j == 1) qd(Q_NE, :) = smp%q(Q_NE, :)
            if (j == 2) qd(Q_TE, :) = smp%q(Q_TE, :)
            do k = 1, 12
                if (j == k + 2) ed(k) = 1.0_dp
            end do
            call model%jvp(smp, theta(3:11), qd, ed(:9), dy(:3))
            call gas_model%jvp(gas_sample, theta(3:14), zero, ed, dy(4:))
            jac(:6, j) = dy/sig
            if (j > 2) jac(j + 4, j) = 1.0_dp
        end do
    end subroutine

    subroutine benchmark()
        type(recorded_thomson_t) :: model
        type(plasma_sample_t) :: sample
        real(dp), allocatable :: eta(:), y(:), dy(:), qd(:, :), ed(:), qb(:, :), eb(:)
        integer(int64) :: begin, finish, rate
        real(dp) :: elapsed, checksum
        integer :: i
        call setup(model, sample, 768)
        allocate (eta(model%n_nuisance()), ed(model%n_nuisance()), eb(model%n_nuisance()), &
                  y(768), dy(768), qd(N_QUANTITIES, sample%n), qb(N_QUANTITIES, sample%n))
        eta = 0.0_dp; ed = .1_dp; qd = 0.0_dp
        qd(Q_NE, :) = 1.0e18_dp; qd(Q_TE, :) = 30.0_dp
        call model%predict(sample, eta, y)
        call system_clock(begin, rate)
        do i = 1, 6
            call model%predict(sample, eta, y)
            call model%jvp(sample, eta, qd, ed, dy)
            call model%vjp(sample, eta, y, qb, eb)
        end do
        call system_clock(finish)
        elapsed = real(finish - begin, dp)/real(rate, dp)
        checksum = sum(y) + sum(dy) + sum(qb) + sum(eb)
        print '(A)', 'channels,volume_points,wavelength_nodes,repetitions,seconds,checksum'
        print '(A,F12.6,A,ES24.16)', '768,3,96,6,', elapsed, ',', checksum
    end subroutine
end program
