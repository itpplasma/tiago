program test_magnetic_pickup
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use tiago_magnetic_pickup, only: magnetic_pickup_t
    use tiago_calibration, only: calibration_prior_t
    implicit none
    integer, parameter :: n = 729
    real(dp), parameter :: pi = acos(-1.0_dp), dt = 1.e-5_dp
    integer(int64), parameter :: origin = 1700000000000000000_int64
    type(magnetic_pickup_t) :: model, bad
    type(calibration_prior_t) :: prior
    integer(int64) :: clock(n), broken(n)
    complex(dp) :: h(n), hbad(n)
    real(dp) :: flux(n), y(n), expected(n), direction(n), tangent(n), &
        bar(n), fluxbar(n), zp(2), zm(2), yplus(n), yminus(n), zbar(2), &
        zdot(2), z(2), cov(2, 2), omega, time, lhs, rhs, eps, nan
    integer :: j, status
    cov = reshape([.04_dp, 2.e-7_dp, 2.e-7_dp, 4.e-10_dp], [2, 2])
    call prior%init([1.7_dp, 3.e-5_dp], cov, status)
    call check(status == 0, 'prior covariance')
    prior%record%source_id = 'synthetic public fixture'
    prior%record%record_id = 'full-chain-voltage-gain-delay'
    prior%record%metadata = 'H: synthetic full chain; V/V; 0 to 1 MHz'
    prior%valid_from = 0
    prior%valid_to = n*dt
    omega = 2*pi*7/(n*dt)
    h = cmplx(0.8_dp, -0.3_dp, dp)
    h(1) = 1
    do j = 2, n/2 + 1
        h(n - j + 2) = conjg(h(j))
    end do
    do j = 1, n
        clock(j) = origin + int(j - 1, int64)*10000_int64
        time = (j - 1)*dt
        flux(j) = 2.e-5_dp*cos(omega*time)
        expected(j) = 2.e-5_dp*omega*1.7_dp* &
            (0.8_dp*sin(omega*(time - 3.e-5_dp)) - &
            0.3_dp*cos(omega*(time - 3.e-5_dp)))
        direction(j) = 1.e-5_dp*sin(2*pi*13*(j - 1)/n)
        bar(j) = cos(2*pi*7*(j - 1)/n) + .2_dp*sin(2*pi*13*(j - 1)/n)
    end do
    z = 0
    call initialize(model, clock, prior, h, status)
    call check(status == 0, 'valid complete chain')
    call model%predict(flux, z, y, status)
    call check(status == 0, 'sinusoid accepted')
    call check(maxval(abs(y - expected)) < 2.e-12_dp, 'independent sine phase delay')
    call model%predict(flux + 1.e-3_dp, z, yplus, status)
    call check(status == 0, 'constant offset accepted')
    call check(maxval(abs(yplus - y)) < 2.e-11_dp, 'Faraday DC null')
    zdot = [.3_dp, -.7_dp]
    call model%jvp(flux, z, direction, zdot, tangent, status)
    call check(status == 0, 'JVP')
    eps = 1.e-5_dp
    zp = z + eps*zdot
    zm = z - eps*zdot
    call model%predict(flux + eps*direction, zp, yplus, status)
    call check(status == 0, 'positive difference')
    call model%predict(flux - eps*direction, zm, yminus, status)
    call check(status == 0, 'negative difference')
    call check(maxval(abs((yplus - yminus)/(2*eps) - tangent)) < 2.e-9_dp, &
        'independent finite difference flux gain delay')
    call model%vjp(flux, z, bar, fluxbar, zbar, status)
    call check(status == 0, 'VJP')
    lhs = dot_product(bar, tangent)
    rhs = dot_product(fluxbar, direction) + dot_product(zbar, zdot)
    call check(abs(lhs - rhs) < 1.e-10_dp*max(1.0_dp, abs(lhs)), 'adjoint duality')
    hbad = 1
    call initialize(model, clock, prior, hbad, status)
    call model%predict(flux, z, y, status)
    expected = 2.e-5_dp*omega*1.7_dp*sin(omega* &
        ([(real(j - 1, dp)*dt, j=1, n)] - 3.e-5_dp))
    call check(status == 0, 'identity electronics')
    call check(maxval(abs(y - expected)) < 2.e-12_dp, 'identity Faraday derivative')
    broken = clock
    broken(4) = broken(4) + 1
    call initialize(bad, broken, prior, h, status)
    call check(status /= 0, 'nonuniform nanosecond clock rejected')
    broken = clock
    broken(4) = broken(3)
    call initialize(bad, broken, prior, h, status)
    call check(status /= 0, 'nonpositive clock rejected')
    broken = clock
    broken(1) = -1
    call initialize(bad, broken, prior, h, status)
    call check(status /= 0, 'negative UTC rejected')
    call bad%init(clock, 0_int64, prior, h, 1.e6_dp, 'Wb', 'V', 'V/V', status)
    call check(status /= 0, 'inexact large reference differences rejected')
    call bad%init(clock(:n - 1), origin, prior, h(:n - 1), 1.e6_dp, &
        'Wb', 'V', 'V/V', status)
    call check(status /= 0, 'even N Nyquist rejected')
    call bad%init(clock, origin, prior, h, 1.e6_dp, 'T', 'V', 'V/V', status)
    call check(status /= 0, 'field units rejected')
    call bad%init(clock, origin, prior, h, 1.e6_dp, 'Wb', 'V', 'A/V', status)
    call check(status /= 0, 'electronics units rejected')
    call bad%init(clock, origin, prior, h, 1.0_dp, 'Wb', 'V', 'V/V', status)
    call check(status /= 0, 'unsupported grid rejected')
    prior%valid_to = dt
    call initialize(bad, clock, prior, h, status)
    call check(status /= 0, 'calibration interval rejected')
    prior%valid_to = n*dt
    prior%record%record_id = ''
    call initialize(bad, clock, prior, h, status)
    call check(status /= 0, 'missing provenance rejected')
    prior%record%record_id = 'full-chain-voltage-gain-delay'
    hbad = h
    hbad(2) = 2
    call initialize(bad, clock, prior, hbad, status)
    call check(status /= 0, 'non-Hermitian response rejected')
    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    hbad = h
    hbad(2) = cmplx(nan, 0.0_dp, dp)
    call initialize(bad, clock, prior, hbad, status)
    call check(status /= 0, 'nonfinite response rejected')
    flux(2) = nan
    call model%predict(flux, z, y, status)
    call check(status /= 0, 'nonfinite flux rejected')
    flux(2) = 0
    z(1) = nan
    call model%predict(flux, z, y, status)
    call check(status /= 0, 'nonfinite calibration rejected')
    z = [0.0_dp, 1.e6_dp]
    call model%predict(flux, z, y, status)
    call check(status /= 0, 'wrapped delay rejected')
    print *, 'magnetic pickup: analytic, derivative and rejection checks passed'
contains
    subroutine initialize(m, t, p, response, code)
        type(magnetic_pickup_t), intent(inout) :: m
        integer(int64), intent(in) :: t(:)
        type(calibration_prior_t), intent(in) :: p
        complex(dp), intent(in) :: response(:)
        integer, intent(out) :: code
        call m%init(t, origin, p, response, 1.e6_dp, 'Wb', 'V', 'V/V', code)
    end subroutine
    subroutine check(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine
end program
