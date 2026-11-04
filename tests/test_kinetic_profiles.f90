program test_kinetic_profiles
    !! Kinetic profiles: the generated pressure kernel against quasi-neutrality
    !! solved independently, the profile and pressure reverse modes against
    !! central differences, interpolation of fitted latents, and a flat-top
    !! profile (plateau) being representable with a finite prior cost.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_kinetic_profiles, only: kinetic_profiles_t, P_NE, P_TE, P_TI, P_ZEFF, P_FAST, &
        N_PROFILES, elementary_charge, KERNEL_MINERVA
    use tiago_kinetic_pressure, only: kinetic_pressure
    implicit none
    type(kinetic_profiles_t) :: kp
    real(dp), allocatable :: z(:), zbar(:), s(:), f(:, :), df(:, :), fbar(:, :), sbar(:)
    real(dp), allocatable :: p(:), dps(:), pbar(:), zp(:), zm(:), fp(:, :), fm(:, :), pp(:), pm(:)
    real(dp) :: lhs, rhs, h, err, g(5), nd, nz, a(2, 2), rhs2(2), det
    integer :: i, j, fails

    fails = 0
    ! --- generated pressure against quasi-neutrality solved as a 2x2 system
    call kinetic_pressure(3.0e19_dp, 1.5e3_dp, 1.0e3_dp, 2.2_dp, 500.0_dp, 6.0_dp, &
        elementary_charge, lhs, g(1), g(2), g(3), g(4), g(5))
    a = reshape([1.0_dp, 1.0_dp, 6.0_dp, 36.0_dp], [2, 2])   ! [n_D + 6 n_Z; n_D + 36 n_Z]
    rhs2 = [3.0e19_dp, 2.2_dp * 3.0e19_dp]
    det = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)
    nd = (rhs2(1) * a(2, 2) - a(1, 2) * rhs2(2)) / det
    nz = (a(1, 1) * rhs2(2) - a(2, 1) * rhs2(1)) / det
    rhs = elementary_charge * (3.0e19_dp * 1.5e3_dp + (nd + nz) * 1.0e3_dp) + 500.0_dp
    call expect(abs(lhs - rhs) < 1.0e-12_dp * rhs, 'pressure kernel vs quasi-neutrality')

    ! --- profiles with every component latent
    call kp%init([0.0_dp, 0.15_dp, 0.35_dp, 0.6_dp, 0.8_dp, 1.0_dp])
    call kp%set_profile(P_NE, [(5.0e19_dp * (1.0_dp - 0.8_dp * kp%greville(i)), i = 1, kp%ncoef)], &
        .true., 0.5_dp, 0.8_dp, 0.1_dp, 0.6_dp, 0.4_dp)
    call kp%set_profile(P_TE, [(2.0e3_dp * (1.0_dp - 0.9_dp * kp%greville(i)), i = 1, kp%ncoef)], &
        .true., 0.5_dp, 0.8_dp, 0.1_dp, 0.6_dp, 0.4_dp)
    call kp%set_profile(P_TI, [(1.5e3_dp * (1.0_dp - 0.9_dp * kp%greville(i)), i = 1, kp%ncoef)], &
        .true., 0.5_dp, 0.8_dp, 0.1_dp, 0.6_dp, 0.4_dp)
    call kp%set_profile(P_ZEFF, [(1.8_dp, i = 1, kp%ncoef)], .true., 0.5_dp, 1.0_dp, 1.0_dp, 0.5_dp, 0.5_dp)
    call kp%set_profile(P_FAST, [(2.0e3_dp * (1.0_dp - kp%greville(i))**2 + 1.0_dp, i = 1, kp%ncoef)], &
        .true., 0.5_dp, 0.8_dp, 0.1_dp, 0.6_dp, 0.4_dp)
    call expect(kp%n_latent == 5 * kp%ncoef, 'latent count')
    allocate(z(kp%n_latent), zbar(kp%n_latent))
    z = [(0.3_dp * sin(1.3_dp * i), i = 1, kp%n_latent)]
    s = [0.0_dp, 0.05_dp, 0.2_dp, 0.47_dp, 0.81_dp, 0.99_dp, 1.08_dp]   ! last: scrape-off
    allocate(f(N_PROFILES, size(s)), df(N_PROFILES, size(s)), fbar(N_PROFILES, size(s)), &
        sbar(size(s)), fp(N_PROFILES, size(s)), fm(N_PROFILES, size(s)))
    call kp%evaluate(z, s, f, df)
    do j = 1, size(s)
        do i = 1, N_PROFILES
            fbar(i, j) = cos(0.7_dp * i + 0.3_dp * j) / max(abs(f(i, j)), 1.0e-30_dp)
        end do
    end do
    zbar = 0.0_dp
    call kp%vjp(z, s, fbar, zbar, sbar)
    h = 1.0e-6_dp
    zp = z + h * [(cos(0.41_dp * i), i = 1, kp%n_latent)]
    zm = z - h * [(cos(0.41_dp * i), i = 1, kp%n_latent)]
    call kp%evaluate(zp, s, fp, df)
    call kp%evaluate(zm, s, fm, df)
    lhs = sum(fbar * (fp - fm)) / (2.0_dp * h)
    rhs = dot_product(zbar, [(cos(0.41_dp * i), i = 1, kp%n_latent)])
    call expect(abs(lhs - rhs) < 1.0e-7_dp * abs(rhs), 'profile vjp in z vs central differences')
    call kp%evaluate(z, s, f, df)
    call kp%evaluate(z, s + h, fp, df)
    call kp%evaluate(z, s - h, fm, df)
    call kp%evaluate(z, s, f, df)
    err = 0.0_dp
    do j = 2, size(s)
        err = max(err, abs(sum(fbar(:, j) * (fp(:, j) - fm(:, j))) / (2.0_dp * h) - sbar(j)) &
            / max(abs(sbar(j)), 1.0_dp))
    end do
    call expect(err < 1.0e-6_dp, 'profile vjp in s vs central differences')

    ! --- pressure reverse mode and s derivative
    allocate(p(size(s)), dps(size(s)), pp(size(s)), pm(size(s)))
    pbar = [(1.0_dp + 0.1_dp * j, j = 1, size(s))]
    call kp%pressure(z, s, p, dps)
    zbar = 0.0_dp
    call kp%pressure_vjp(z, s, pbar, zbar)
    call kp%pressure(zp, s, pp)
    call kp%pressure(zm, s, pm)
    lhs = dot_product(pbar, pp - pm) / (2.0_dp * h)
    rhs = dot_product(zbar, [(cos(0.41_dp * i), i = 1, kp%n_latent)])
    call expect(abs(lhs - rhs) < 1.0e-7_dp * abs(rhs), 'pressure vjp vs central differences')
    call kp%pressure(z, s + h, pp)
    call kp%pressure(z, s - h, pm)
    call expect(maxval(abs((pp(2:) - pm(2:)) / (2.0_dp * h) - dps(2:))) < 1.0e-6_dp * maxval(abs(dps)), &
        'pressure s derivative')

    ! --- fitted latents interpolate given values, a plateau has finite cost
    z = 0.0_dp
    call kp%fit_latent(P_TE, kp%greville, [(min(2.0e3_dp, 2.6e3_dp * (1.0_dp - kp%greville(i))) + 10.0_dp, &
        i = 1, kp%ncoef)], z)
    call kp%evaluate(z, kp%greville, f, df)
    call expect(maxval(abs(f(P_TE, :) - [(min(2.0e3_dp, 2.6e3_dp * (1.0_dp - kp%greville(i))) + 10.0_dp, &
        i = 1, kp%ncoef)])) < 1.0e-9_dp * 2.0e3_dp, 'fit_latent interpolates')
    ! A local flattening with zero gradient at s = 0.45 (width 0.3, as from a
    ! finite-transport island) under a prior centred on the profile without it
    ! (kernel scales 0.3 core / 0.15 edge, sigma_f = 1): its prior cost |z|^2
    ! lies below the 99.9 % chi^2 quantile for 8 degrees of freedom (26.1).
    ! An exp-gradient prior would give this profile zero density. Sharp-cornered
    ! The squared-exponential (Minerva) kernel puts it far outside (|z|^2 ~ 130).
    call kp%set_profile(P_TE, [(te_ramp(kp%greville(i), .false.), i = 1, kp%ncoef)], &
        .true., 1.0_dp, 0.3_dp, 0.15_dp, 0.6_dp, 0.2_dp)
    z = 0.0_dp
    call kp%fit_latent(P_TE, kp%greville, [(te_ramp(kp%greville(i), .true.), i = 1, kp%ncoef)], z)
    lhs = sum(z(kp%p(P_TE)%offset + 1:kp%p(P_TE)%offset + kp%ncoef)**2)
    call expect(kp%ncoef == 8, 'eight coefficients for six breakpoints')
    call expect(lhs < 26.1_dp, 'flattened profile inside the Matern-5/2 prior bulk')
    call kp%set_profile(P_TE, [(te_ramp(kp%greville(i), .false.), i = 1, kp%ncoef)], &
        .true., 1.0_dp, 0.3_dp, 0.15_dp, 0.6_dp, 0.2_dp, KERNEL_MINERVA)
    z = 0.0_dp
    call kp%fit_latent(P_TE, kp%greville, [(te_ramp(kp%greville(i), .true.), i = 1, kp%ncoef)], z)
    rhs = sum(z(kp%p(P_TE)%offset + 1:kp%p(P_TE)%offset + kp%ncoef)**2)
    call expect(rhs > 50.0_dp, 'squared-exponential kernel penalizes the flattening')
    ! The optional peaked component includes the thermal on-axis pressure chain.
    block
        type(kinetic_profiles_t) :: peaked
        real(dp), allocatable :: zz(:), zb(:), za(:), base(:)
        real(dp) :: ss(4), val(4), plus(4), minus(4), slope(4), central(1), v0(1)
        real(dp) :: thermal_plus(4), thermal_minus(4), thermal_slope(4)
        real(dp) :: c, exponent, direction
        integer :: col
        ss = [0.0_dp, 0.3_dp, 0.8_dp, 1.0_dp]
        call peaked%init([0.0_dp, 0.5_dp, 1.0_dp])
        call peaked%set_profile(P_NE, [(1.0e19_dp, i=1, peaked%ncoef)], &
            .true., 0.2_dp, 0.8_dp, 0.2_dp, 0.5_dp, 0.3_dp)
        call peaked%set_profile(P_TE, [(800.0_dp, i=1, peaked%ncoef)], .false.)
        call peaked%set_profile(P_TI, [(600.0_dp, i=1, peaked%ncoef)], .false.)
        allocate(zz(peaked%n_latent))
        zz = 0
        allocate(base(4))
        call peaked%pressure(zz, ss, base, thermal_slope)
        call peaked%pressure(zz, [0.0_dp], central)
        call peaked%set_nonthermal([0.0_dp, 2.0_dp], [2.0_dp, 8.0_dp])
        zz = [zz, 0.3_dp, -0.4_dp]
        allocate(zb(size(zz)))
        call peaked%pressure(zz, ss, val, slope)
        c = 2.0_dp / (1.0_dp + exp(-zz(size(zz)-1)))
        exponent = 2.0_dp + 6.0_dp / (1.0_dp + exp(-zz(size(zz))))
        call expect(maxval(abs(val - base - c*central(1)*(1-ss)**exponent)) < 1.0e-10_dp, &
            'bounded peaked pressure vs independent formula')
        zb = 0
        call peaked%pressure_vjp(zz, ss, [1.0_dp, -0.3_dp, 0.7_dp, 0.2_dp], zb)
        do col = 1, size(zz)
            za = zz
            za(col) = zz(col) + h
            call peaked%pressure(za, ss, plus)
            za(col) = zz(col) - h
            call peaked%pressure(za, ss, minus)
            direction = dot_product([1.0_dp, -0.3_dp, 0.7_dp, 0.2_dp], (plus-minus)/(2*h))
            call expect(abs(zb(col)-direction) < 1.0e-6_dp*max(abs(direction), 1.0_dp), &
                'peaked pressure VJP includes amplitude, exponent and thermal p0')
        end do
        call peaked%pressure(zz, ss+h, plus)
        call peaked%pressure(zz, ss-h, minus)
        call peaked%set_nonthermal([0.0_dp, 2.0_dp], [2.0_dp, 8.0_dp], enabled=.false.)
        call peaked%pressure(zz, ss+h, thermal_plus)
        call peaked%pressure(zz, ss-h, thermal_minus)
        call expect(maxval(abs(slope-thermal_slope &
            -(plus-minus-thermal_plus+thermal_minus)/(2*h))) < 1.0e-4_dp, &
            'peaked component flux derivative and zero edge continuation')
        call peaked%set_nonthermal([0.0_dp, 2.0_dp], [2.0_dp, 8.0_dp])
        zz(size(zz)-1:) = [30.0_dp, -30.0_dp]
        call peaked%pressure(zz, [0.0_dp], v0)
        call expect(v0(1) <= 3.0_dp*central(1), 'nonthermal amplitude obeys prior bound')
        zz(size(zz)-1:) = [-1000.0_dp, 1000.0_dp]
        call peaked%pressure(zz, ss, val)
        zb = 0
        call peaked%pressure_vjp(zz, ss, [1.0_dp, -0.3_dp, 0.7_dp, 0.2_dp], zb)
        call expect(all(ieee_is_finite(val)) .and. all(ieee_is_finite(zb)), &
            'nonthermal transform and derivatives remain finite in prior tails')
        call peaked%set_nonthermal([0.0_dp, 2.0_dp], [2.0_dp, 8.0_dp], enabled=.false.)
        call peaked%pressure(zz(:peaked%n_latent), ss, val)
        call expect(maxval(abs(val-base)) < 1.0e-12_dp, 'thermal-only model remains selectable')
    end block
    if (fails > 0) error stop 1
    print '(A)', 'kinetic profiles: all checks passed'
contains
    real(dp) function te_ramp(x, flat)
        real(dp), intent(in) :: x
        logical, intent(in) :: flat
        real(dp), parameter :: sc = 0.45_dp, w = 0.3_dp
        real(dp) :: g
        g = x
        if (flat .and. abs(x - sc) < 0.5_dp * w) g = x - w / (2.0_dp * acos(-1.0_dp)) * &
            sin(2.0_dp * acos(-1.0_dp) * (x - sc) / w)
        te_ramp = 2.6e3_dp * (1.0_dp - g) + 10.0_dp
    end function te_ramp
    subroutine expect(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            fails = fails + 1
            write(error_unit, '(A)') 'FAIL: '//what
        end if
    end subroutine expect
end program test_kinetic_profiles
