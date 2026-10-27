program test_kinetic_profiles
    !! Kinetic profiles: the generated pressure kernel against quasi-neutrality
    !! solved independently, the profile and pressure reverse modes against
    !! central differences, interpolation of fitted latents, and a flat-top
    !! profile (plateau) being representable with a finite prior cost.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
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
