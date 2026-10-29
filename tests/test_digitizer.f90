program test_digitizer
    use, intrinsic :: iso_fortran_env, only: dp => real64, int32, error_unit
    use tiago_profile_diagnostic, only: plasma_sample_t, N_QUANTITIES, Q_TE, Q_NE, &
        check_diagnostic_derivatives
    use tiago_point_probe, only: point_probe_t
    use tiago_integrator_diagnostic, only: integrator_diagnostic_t
    use tiago_digitizer_diagnostic, only: digitizer_diagnostic_t
    use tiago_calibration, only: calibration_prior_t
    use tiago_adc, only: adc_decode, adc_quantize, adc_quantization_variance
    implicit none
    integer, parameter :: n = 100000, bits = 16
    real(dp), parameter :: vres = 1.0e-3_dp, voff = 0.5_dp
    real(dp) :: err, adj, variance, mean, stride, u
    integer :: i, status
    call check_adc()
    call check_chain()
    call check_variance()
    print '(A)', 'digitiser stage: PASS'
contains
    subroutine check(ok, text)
        logical, intent(in) :: ok
        character(*), intent(in) :: text
        if (ok) return
        write (error_unit, '(A)') 'FAIL: '//text
        error stop 1
    end subroutine

    subroutine check_adc()
        integer, parameter :: m = 64
        real(dp) :: v(m), round_trip
        integer(int32) :: c(m)
        integer(int32), parameter :: lo = shiftl(-1_int32, bits - 1)
        integer(int32), parameter :: hi = shiftl(1_int32, bits - 1) - 1_int32
        logical :: sat(m)
        integer :: j
        ! In-band round trip: decode(quantize(v)) is within half a code.
        do j = 1, m
            v(j) = voff + sin(0.7_dp*real(j, dp))*30000.0_dp*vres
        end do
        do j = 1, m
            call adc_quantize(v(j), vres, voff, bits, c(j), sat(j))
        end do
        call check(.not. any(sat), 'in-band samples are not saturated')
        round_trip = maxval(abs(adc_decode(c, spread(vres, 1, m), &
            spread(voff, 1, m)) - v))
        call check(round_trip <= 0.5_dp*vres, 'round trip error is at most half a code')
        ! Exactly at the code limits: no saturation; beyond: clipped and flagged.
        call adc_quantize(voff + real(hi, dp)*vres, vres, voff, bits, c(1), sat(1))
        call check(c(1) == hi .and. .not. sat(1), 'positive code limit is exact')
        call adc_quantize(voff + real(hi + 1, dp)*vres, vres, voff, bits, c(2), sat(2))
        call check(c(2) == hi .and. sat(2), 'above full scale saturates high')
        call adc_quantize(voff + real(lo - 1, dp)*vres, vres, voff, bits, c(3), sat(3))
        call check(c(3) == lo .and. sat(3), 'below full scale saturates low')
        call adc_quantize(voff + real(lo, dp)*vres, vres, voff, bits, c(4), sat(4))
        call check(c(4) == lo .and. .not. sat(4), 'negative code limit is exact')
        call check(abs(adc_quantization_variance(vres) - vres*vres/12.0_dp) &
            < tiny(1.0_dp), 'quantisation variance is vres^2/12')
        print '(A,ES12.4)', 'adc round trip max error [V]: ', round_trip
    end subroutine

    subroutine check_chain()
        type(point_probe_t) :: probe
        type(integrator_diagnostic_t) :: recorder
        type(digitizer_diagnostic_t) :: digitizer
        type(calibration_prior_t) :: iprior, dprior
        type(plasma_sample_t) :: sample
        real(dp) :: eta(7), y(2), expected(2), analytic
        real(dp), allocatable :: x(:, :), dx(:, :, :)
        ! Toy probe from tests/test_recorded_magnetics.f90: two fixed channels.
        probe%x0 = reshape([2.0_dp, 0.0_dp, 0.0_dp, 2.1_dp, 0.0_dp, 0.1_dp], [3, 2])
        probe%quantity = [Q_TE, Q_NE]
        call iprior%init([2.0_dp, 0.3_dp, 0.01_dp], &
            reshape([0.04_dp, 0.01_dp, 0.0_dp, &
            0.01_dp, 0.09_dp, 0.003_dp, &
            0.0_dp, 0.003_dp, 0.01_dp], [3, 3]), status)
        call check(status == 0, 'integrator calibration prior')
        call recorder%init(probe, iprior, [0.0_dp, 0.4_dp], 0.0_dp, status)
        call check(status == 0, 'integrator stage')
        ! Per-channel [vres_1, vres_2, voff_1, voff_2], deliberately asymmetric.
        call dprior%init([1.0e-3_dp, 2.0e-3_dp, 0.1_dp, -0.2_dp], &
            reshape([4.0e-8_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 9.0e-8_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 1.0e-4_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 1.6e-4_dp], [4, 4]), status)
        call check(status == 0, 'digitiser calibration prior')
        call digitizer%init(recorder, dprior, bits, status)
        call check(status == 0, 'digitiser stage')
        call check(digitizer%n_signals() == 2, 'signal count')
        call check(digitizer%n_nuisance() == 7, 'integrator and digitiser latents')
        sample%n = 2
        allocate (sample%q(N_QUANTITIES, 2))
        sample%q = 0.0_dp
        sample%q(Q_TE, 1) = 100.0_dp
        sample%q(Q_NE, 2) = 1.0e19_dp
        ! Zero standardized latents: physical parameters equal the prior means.
        eta = 0.0_dp
        call digitizer%sample_points(eta, x, dx)
        call check(maxval(abs(dx(:, :, 3:))) == 0.0_dp, 'calibration fixes geometry')
        call digitizer%predict(sample, eta, y)
        ! Expected counts by hand: (gain*q + offset + drift*t - voff)/vres.
        expected(1) = ((2.0_dp*100.0_dp + 0.3_dp) - 0.1_dp)/1.0e-3_dp
        expected(2) = ((2.0_dp*1.0e19_dp + 0.3_dp + 0.01_dp*0.4_dp) + 0.2_dp)/2.0e-3_dp
        analytic = maxval(abs(y - expected)/abs(expected))
        call check(analytic < 1.0e-12_dp, 'analytic counts and parameter order')
        eta(4:) = [0.1_dp, -0.2_dp, 0.3_dp, -0.1_dp]
        call check_diagnostic_derivatives(digitizer, sample, eta, 1.0e-5_dp, err, adj)
        call check(err < 1.0e-7_dp, 'digitiser chain finite differences')
        call check(adj < 1.0e-12_dp, 'digitiser chain transpose')
        print '(A,2ES12.4)', 'digitiser FD / transpose relative errors: ', err, adj
        call reject_bad_priors()
    end subroutine

    subroutine reject_bad_priors()
        type(point_probe_t) :: probe
        type(digitizer_diagnostic_t) :: digitizer
        type(calibration_prior_t) :: prior
        probe%x0 = reshape([2.0_dp, 0.0_dp, 0.0_dp], [3, 1])
        probe%quantity = [Q_TE]
        call prior%init([1.0e-3_dp, 0.1_dp], &
            reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), status)
        call digitizer%init(probe, prior, bits, status)
        call check(status == 0, 'positive volts-per-count accepted')
        call prior%init([-1.0e-3_dp, 0.1_dp], &
            reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), status)
        call digitizer%init(probe, prior, bits, status)
        call check(status /= 0, 'negative volts-per-count rejected')
        call prior%init([0.0_dp, 0.1_dp], &
            reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), status)
        call digitizer%init(probe, prior, bits, status)
        call check(status /= 0, 'zero volts-per-count rejected')
    end subroutine

    subroutine check_variance()
        ! Uniform quantisation error in code units has variance 1/12 when the
        ! signal spans many codes (Widrow and Kollar 2008, ch. 4). A golden
        ! irrational stride spreads the sweep over 20000 codes without lock-in.
        integer(int32) :: c(n)
        real(dp) :: volts(n), res(n), ofs(n)
        logical :: sat(n)
        real(dp) :: golden
        integer :: j
        golden = sqrt(2.0_dp) - 1.0_dp
        do j = 1, n
            u = real(j, dp)*golden
            u = u - real(int(u), dp)
            volts(j) = vres*((u - 0.5_dp)*20000.0_dp)
        end do
        res = vres
        ofs = 0.0_dp
        call adc_quantize(volts, res, ofs, bits, c, sat)
        call check(.not. any(sat), 'noise sweep stays inside the code range')
        ! Quantisation error in code units: nint(u) - u, variance 1/12.
        volts = real(c, dp) - volts/vres
        mean = sum(volts)/n
        variance = sum((volts - mean)**2)/n
        err = abs(variance - 1.0_dp/12.0_dp)/(1.0_dp/12.0_dp)
        call check(err < 2.0e-2_dp, 'empirical quantisation variance equals 1/12')
        print '(A,ES12.4)', 'empirical quantisation variance / (1/12): ', &
            variance*12.0_dp
        print '(A,ES12.4)', 'variance relative error: ', err
    end subroutine
end program
