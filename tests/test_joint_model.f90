program test_joint_model
    !! Joint equilibrium/profile chain on pinned LI383 (fixed boundary):
    !!  1. a line-segment kinetic pressure equal to the input's power-series
    !!     pressure on the half grid reproduces the power-series equilibrium;
    !!  2. the flux map inverts the geometry (round trip of (s, theta, phi));
    !!  3. the joint Jacobian (magnetics, profile probes, |B| probes at fixed
    !!     points, alignment; implicit VMEC++ adjoint including |B|, kinetic
    !!     pressure, flux-map VJP) agrees with central
    !!     differences of re-solved equilibria for every latent kind;
    !!  4. MAP and Laplace run on noise-free synthetic data and the evidence is
    !!     finite.
    !! Usage: test_joint_model <VMEC input> <sensors.awk> <workdir>
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, error_unit
    use tiago_reconstruction, only: rec, jacobian_at
    use tiago_joint_model, only: joint_model_t, diagnostic_slot_t
    use tiago_synthetic_case, only: setup_magnetics, place_probes, synthesize_magnetics, &
        synthesize_profiles
    use tiago_kinetic_profiles, only: P_NE, P_TE, P_TI, P_ZEFF
    use tiago_real_space_state, only: located_t, locate_ok
    use tiago_point_probe, only: point_probe_t
    use tiago_profile_diagnostic, only: Q_B
    use tiago_inference_engines, only: map_fit, laplace, laplace_t
    use fortnum_rng, only: rng_t, rng_seed
    use fortnum_status, only: fortnum_status_t
    implicit none
    character(len=512) :: input, awk, work
    type(joint_model_t), target :: model
    type(rng_t) :: rng
    type(fortnum_status_t) :: st
    type(laplace_t) :: lap
    type(located_t) :: loc(3)
    real(dp), allocatable :: theta0(:), theta_true(:), theta(:), r(:), rp(:), rm(:), jac(:, :)
    real(dp), allocatable :: y_series(:), p(:), s(:), am(:), colfd(:)
    real(dp) :: h, err, rr, zz, d(2, 3), pres_scale, chi2_start
    integer :: i, k, fails, steps, col
    integer, allocatable :: cols(:)
    character(len=:), allocatable :: message
    logical :: ok

    fails = 0
    call get_command_argument(1, input)
    call get_command_argument(2, awk)
    call get_command_argument(3, work)
    call execute_command_line('mkdir -p "'//trim(work)//'" && awk -v R0=1.42 -v rs=0.85 -v nphi=4 '// &
        '-v nth=4 -v nprobe=12 -v dir="'//trim(work)//'" -f "'//trim(awk)//'"')
    call setup_magnetics(trim(input), trim(work)//'/work', trim(work), &
        [character(len=64) :: 'CURTOR', 'PHIEDGE'], [-174250.0_dp, 0.514386_dp], &
        [2.0e4_dp, 0.02_dp], 1.0e-14_dp)

    ! --- 1. line-segment pressure reproduces the power-series equilibrium
    call rec%eq%solve([-174250.0_dp, 0.514386_dp], ok)
    if (.not. ok) error stop 'no reference equilibrium'
    y_series = rec%eq%y
    s = rec%eq%s_half()
    pres_scale = rec%eq%vmec%get_input('pres_scale')
    am = [(rec%eq%vmec%get_input('am', k), k = 1, rec%eq%vmec%input_length('am'))]
    p = [(pres_scale * sum([(am(k) * s(i)**(k - 1), k = 1, size(am))]), i = 1, size(s))]
    call rec%eq%set_pressure(p)
    call rec%eq%solve([-174250.0_dp, 0.514386_dp], ok)
    err = maxval(abs(rec%eq%y - y_series)) / maxval(abs(y_series))
    print '(A,ES10.2)', 'line-segment vs power-series pressure, max rel. edge difference: ', err
    call expect(ok .and. err < 1.0e-8_dp, 'line-segment pressure reproduces the power series')

    ! --- kinetic profiles: ne, Te, Ti latent, Zeff fixed
    call model%prof%init([0.0_dp, 0.2_dp, 0.45_dp, 0.7_dp, 1.0_dp])
    associate (g => model%prof%greville)
        call model%prof%set_profile(P_NE, 1.0e20_dp * (1.0_dp - 0.6_dp * g), .true., &
            0.3_dp, 0.4_dp, 0.2_dp, 0.6_dp, 0.2_dp)
        call model%prof%set_profile(P_TE, 2.5e3_dp * (1.0_dp - g) + 50.0_dp, .true., &
            0.3_dp, 0.4_dp, 0.2_dp, 0.6_dp, 0.2_dp)
        call model%prof%set_profile(P_TI, 2.0e3_dp * (1.0_dp - g) + 50.0_dp, .true., &
            0.3_dp, 0.4_dp, 0.2_dp, 0.6_dp, 0.2_dp)
        call model%prof%set_profile(P_ZEFF, 1.5_dp + 0.0_dp * g, .false.)
    end associate
    allocate(model%diags(0))
    call model%setup()
    theta0 = [rec%prior / rec%scale, [(0.0_dp, i = 1, model%n_prof)]]
    allocate(r(model%n_residuals()))
    call model%residual(theta0, r, ok)
    if (.not. ok) error stop 'no equilibrium at the prior mean'

    ! --- 2. point location round trip
    err = 0.0_dp
    do i = 1, 3
        call model%state%point(0.3_dp * i, 0.7_dp * i, 0.2_dp * i, rr, zz, d)
        call model%state%locate([rr * cos(0.2_dp * i), rr * sin(0.2_dp * i), zz], loc(i))
        err = max(err, maxval(abs(loc(i)%u - [0.3_dp * i, 0.7_dp * i, 0.2_dp * i])))
    end do
    call expect(err < 1.0e-10_dp .and. all(loc%status == locate_ok), 'point location round trip')
    ! state contract: label derivatives of X and |B| against central differences,
    ! param_vjp as the adjoint of param_jvp
    block
        real(dp) :: u(3), xu(3), dxdu(3, 3), b0, dbdu(3), xp(3), xm(3), bp, bm, jt(3, 3), dbt(3)
        real(dp) :: dx(3), db, xbar(3), lhs, rhs
        real(dp), allocatable :: pdot(:), pbar(:)
        u = [0.55_dp, 1.1_dp, 0.3_dp]
        call model%state%eval(u, xu, dxdu, b0, dbdu)
        err = 0.0_dp
        do k = 1, 3
            u(k) = u(k) + 1.0e-6_dp
            call model%state%eval(u, xp, jt, bp, dbt)
            u(k) = u(k) - 2.0e-6_dp
            call model%state%eval(u, xm, jt, bm, dbt)
            u(k) = u(k) + 1.0e-6_dp
            err = max(err, maxval(abs((xp - xm) / 2.0e-6_dp - dxdu(:, k))) / maxval(abs(dxdu)), &
                abs((bp - bm) / 2.0e-6_dp - dbdu(k)) / max(maxval(abs(dbdu)), 1.0e-12_dp))
        end do
        call expect(err < 1.0e-7_dp .and. b0 > 0.0_dp, 'state label derivatives of X and |B|')
        allocate(pdot(model%state%n_params()), pbar(model%state%n_params()))
        pdot = [(sin(0.37_dp * i), i = 1, size(pdot))]
        xbar = [0.3_dp, -0.7_dp, 0.5_dp]
        call model%state%param_jvp(u, pdot, dx, db)
        pbar = 0.0_dp
        call model%state%param_vjp(u, xbar, 0.9_dp, pbar)
        lhs = dot_product(xbar, dx) + 0.9_dp * db
        rhs = dot_product(pbar, pdot)
        call expect(abs(lhs - rhs) < 1.0e-12_dp * max(abs(lhs), 1.0_dp), 'state param_vjp adjoint of param_jvp')
    end block

    ! --- probes and noise-free synthetic data at a perturbed truth
    call place_probes(model, [0.25_dp, 0.45_dp, 0.65_dp, 0.85_dp], 0.3_dp, 0.0_dp, &
        acos(-1.0_dp) / 6.0_dp, 0.005_dp)
    call add_b_probes()
    call model%setup()
    deallocate(r)
    allocate(r(model%n_residuals()), rp(model%n_residuals()), rm(model%n_residuals()))
    theta0 = [rec%prior / rec%scale, [(0.0_dp, i = 1, model%n - model%n_rec)]]
    theta_true = theta0 + [0.03_dp, -0.02_dp, [(0.6_dp * sin(1.7_dp * i), i = 1, model%n_prof)], &
        0.5_dp, -0.4_dp]
    call rng_seed(rng, 1_int64, st)
    call model%residual(theta_true, r, ok)
    call synthesize_magnetics(model, theta_true(1:2) * rec%scale, 0.01_dp, 1.0e-4_dp, 1.0e-8_dp, &
        1.0e-7_dp, rng, .false.)
    call synthesize_profiles(model, theta_true, 0.05_dp, [0.0_dp, 0.0_dp, 0.0_dp], rng, .false.)

    ! --- 3. Jacobian against central differences of re-solved equilibria;
    !     cold solves, since hot restarts from different start states leave
    !     non-smooth residuals at FTOL (as in test_vmecpp_equilibrium)
    call rec%eq%vmec%set_input('hot_restart', 0.0_dp)
    theta = 0.5_dp * (theta0 + theta_true)
    call model%residual(theta, r, ok)
    call model%jacobian(theta, jac)
    ! magnetic Fisher information of each rec parameter: the joint magnetic
    ! rows must carry the synthesized measurement sigmas, i.e. agree with the
    ! reconstruction rows whitened by rec%meas_sigma (regression: a covariance
    ! factored before synthesize_magnetics kept unit sigmas)
    block
        real(dp), allocatable :: jrec(:, :)
        real(dp) :: f_joint, f_rec, prior_sd
        call jacobian_at(theta(1:model%n_rec) * rec%scale, jrec)
        do k = 1, model%n_rec
            f_joint = sum(jac(1:model%n_mag, k)**2)
            f_rec = sum((jrec(1:model%n_mag, k) * rec%scale(k))**2)
            prior_sd = rec%prior_sigma(k) / rec%scale(k)
            print '(2A,3ES11.3)', trim(rec%names(k)), ': magnetic Fisher joint/rec, '// &
                'post/prior sd from magnetics + prior: ', f_joint, f_rec, &
                1.0_dp / sqrt(1.0_dp + f_joint * prior_sd**2)
            call expect(abs(f_joint - f_rec) <= 1.0e-10_dp * f_rec, &
                'joint magnetic rows whitened by the measurement sigmas')
        end do
    end block
    cols = [1, 2, model%n_rec + 3, model%n_rec + 7 + 4, model%n_rec + 14 + 2, model%n - 1, model%n]
    h = 1.0e-4_dp
    block
        character(len=32) :: hs
        integer :: stat
        call get_environment_variable('TIAGO_FD_STEP', hs, status=stat)
        if (stat == 0) read(hs, *) h
    end block
    do k = 1, size(cols)
        col = cols(k)
        theta(col) = theta(col) + h
        call model%residual(theta, rp, ok)
        theta(col) = theta(col) - 2.0_dp * h
        call model%residual(theta, rm, ok)
        theta(col) = theta(col) + h
        colfd = (rp - rm) / (2.0_dp * h)
        err = maxval(abs(colfd - jac(:, col))) / max(maxval(abs(colfd)), 1.0e-12_dp)
        print '(A,I3,A,ES10.2,A,ES10.2)', 'column ', col, ': max |J - FD| / max |FD| = ', err, &
            '   max |FD| = ', maxval(abs(colfd))
        print '(A,3ES10.2)', '   block max|FD|: ', maxval(abs(colfd(1:model%n_mag))), &
            maxval(abs(colfd(model%n_mag + 1:model%n_mag + model%diags(1)%d%n_signals()))), &
            maxval(abs(colfd(model%n_mag + model%diags(1)%d%n_signals() + 1:model%n_mag + model%n_diag_rows)))
        print '(A,3ES10.2)', '   blocks mag/probe-a/probe-b: ', &
            blockerr(1, model%n_mag), blockerr(model%n_mag + 1, model%n_mag + model%diags(1)%d%n_signals()), &
            blockerr(model%n_mag + model%diags(1)%d%n_signals() + 1, model%n_mag + model%n_diag_rows)
        call expect(err < 2.0e-4_dp, 'Jacobian column vs central differences')
    end do

    call rec%eq%vmec%set_input('hot_restart', 1.0_dp)
    ! --- 4. MAP and Laplace on noise-free data from the prior mean
    theta = theta0
    call model%residual(theta, r, ok)
    chi2_start = sum(r**2)
    call map_fit(model, theta, 30, steps, message)
    call model%residual(theta, r, ok)
    print '(A,I3,A,ES11.3,A,ES11.3,2A)', 'MAP: ', steps, ' steps, chi^2 ', chi2_start, ' -> ', &
        sum(r**2), ', ', message
    call expect(sum(r**2) < 0.1_dp * chi2_start, 'MAP lowers chi^2')
    call laplace(model, theta, lap)
    call expect(lap%info == 0, 'Laplace Hessian positive definite')
    call expect(abs(lap%log_evidence) < huge(1.0_dp) .and. lap%proper, 'finite evidence')
    print '(A,F12.3,A,F10.3)', 'log evidence ', lap%log_evidence, ', chi^2 ', lap%chi2
    print '(A,*(F8.3))', 'pull of truth (theta* - truth)/sd: ', &
        [((theta(i) - theta_true(i)) / sqrt(lap%cov(i, i)), i = 1, model%n)]
    if (fails > 0) error stop 1
    print '(A)', 'joint model: all checks passed'
contains
    subroutine add_b_probes()
        !! Third diagnostic: |B| at the Ti points (fixed in space), so the
        !! Jacobian check covers |B| at fixed points through the equilibrium.
        type(point_probe_t) :: bp
        type(point_probe_t), allocatable :: ti
        type(diagnostic_slot_t), allocatable :: old(:)
        integer :: j
        select type (q => model%diags(2)%d)
        type is (point_probe_t)
            allocate(ti, source=q)
        end select
        bp%x0 = ti%x0
        bp%quantity = [(Q_B, j = 1, size(ti%quantity))]
        bp%aligned = .true.
        call move_alloc(model%diags, old)
        allocate(model%diags(3))
        model%diags(1:2) = old
        allocate(model%diags(3)%d, source=bp)
        allocate(model%diags(3)%eta_sigma(0))
        model%diags(3)%gain_sigma = 0.0_dp
        model%diags(3)%data = [(0.0_dp, j = 1, size(bp%quantity))]
        model%diags(3)%covariance = reshape([(merge(1.0_dp, 0.0_dp, mod(j, size(bp%quantity) + 1) == 0), &
            j = 0, size(bp%quantity)**2 - 1)], [size(bp%quantity), size(bp%quantity)])
        model%diags(3)%model_sigma = [(0.0_dp, j = 1, size(bp%quantity))]
    end subroutine add_b_probes

    real(dp) function blockerr(a, b)
        integer, intent(in) :: a, b
        blockerr = maxval(abs(colfd(a:b) - jac(a:b, col))) / max(maxval(abs(colfd(a:b))), 1.0e-12_dp)
    end function blockerr
    subroutine expect(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            fails = fails + 1
            write(error_unit, '(A)') 'FAIL: '//what
        end if
    end subroutine expect
end program test_joint_model
