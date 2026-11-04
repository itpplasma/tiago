module mixture_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_local_posterior, only: local_posterior_t
    use fortnum_rng, only: rng_t, rng_normal
    implicit none
    type, extends(local_posterior_t) :: fixture_t
        logical :: folded = .true., proper = .true., bounded = .false.
        real(dp) :: observed = 6, sigma = 0.25_dp, prior_mean = 0.5_dp, prior_sd = 3
        real(dp) :: current = 0
    contains
        procedure :: residual
        procedure :: jacobian
        procedure :: n_residuals
        procedure :: log_normalization
        procedure :: proper_prior
        procedure :: draw_prior
        procedure :: save_equilibrium
    end type
contains
    subroutine residual(self, theta, r, ok)
        class(fixture_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), intent(out) :: r(:)
        logical, intent(out) :: ok
        real(dp) :: prediction
        self%current = theta(1)
        prediction = theta(1)
        if (self%folded) prediction = abs(prediction)
        r = [(prediction - self%observed)/self%sigma, &
            (theta(1) - self%prior_mean)/self%prior_sd]
        ok = theta(1) /= 0
        if (self%bounded) ok = ok .and. theta(1) <= 4
    end subroutine
    subroutine jacobian(self, theta, jac)
        class(fixture_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), allocatable, intent(out) :: jac(:, :)
        allocate (jac(2, 1))
        jac(:, 1) = [1/self%sigma, 1/self%prior_sd]
        if (self%folded) jac(1, 1) = sign(jac(1, 1), theta(1))
    end subroutine
    integer function n_residuals(self) result(n)
        class(fixture_t), intent(in) :: self
        n = 2
    end function
    real(dp) function log_normalization(self) result(c)
        class(fixture_t), intent(in) :: self
        c = -log(2*acos(-1.0_dp)*self%sigma*self%prior_sd)
    end function
    logical function proper_prior(self) result(proper)
        class(fixture_t), intent(in) :: self
        proper = self%proper
    end function
    subroutine draw_prior(self, rng, theta, ok)
        class(fixture_t), intent(in) :: self
        type(rng_t), intent(inout) :: rng
        real(dp), intent(out) :: theta(:)
        logical, intent(out) :: ok
        call rng_normal(rng, theta(1))
        theta = self%prior_mean + self%prior_sd*theta
        ok = self%proper
    end subroutine
    subroutine save_equilibrium(self, prefix, path)
        class(fixture_t), intent(inout) :: self
        character(len=*), intent(in) :: prefix
        character(len=:), allocatable, intent(out) :: path
        integer :: u
        path = prefix//'.state'
        open (newunit=u, file=path, status='replace')
        write (u, *) self%current
        close (u)
    end subroutine
end module

program test_local_mixture
    use fortnum_rng, only: rng_seed
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use mixture_fixture
    use tiago_local_mixture, only: mode_start_t, mixture_options_t, local_mixture_t, &
        multi_start_map, &
        sample_local_mixture
    use fortnum_status, only: fortnum_status_t
    implicit none
    type(fixture_t) :: model
    type(mode_start_t) :: starts(36), family(2), bad(1)
    type(mixture_options_t) :: opts
    type(local_mixture_t) :: mixture
    type(fortnum_status_t) :: status
    type(rng_t) :: rng
    real(dp), allocatable :: sample(:)
    character(len=512) :: directory
    real(dp) :: variance, expected_mean, expected_z, sign_mode, total_variance, &
        ratio, saved
    integer :: i, j, positive, negative, u, fails
    fails = 0
    call get_command_argument(1, directory)
    opts%decrease_tolerance = 1.0e-12_dp
    opts%psis_draws = 200
    opts%seed = 817
    starts(1)%theta = [-4.0_dp]
    starts(2)%theta = [4.0_dp]
    starts(3)%theta = [-5.0_dp]
    starts(4)%theta = [5.0_dp]
    variance = 1/(1/model%sigma**2 + 1/model%prior_sd**2)
    total_variance = model%sigma**2 + model%prior_sd**2
    call multi_start_map(model, starts, 1, trim(directory)//'/two_mode', opts, &
        mixture, status)
    call expect(status%code == 0, 'two-mode fit status')
    call expect(size(mixture%modes) == 2, &
        'both planted modes; duplicate and prior starts collapse')
    if (size(mixture%modes) /= 2) stop 1
    positive = 0
    negative = 0
    do i = 1, 2
        associate (m => mixture%modes(i), g => mixture%modes(i)%gaussian)
            sign_mode = sign(1.0_dp, g%mean(1))
            if (sign_mode > 0) positive = i
            if (sign_mode < 0) negative = i
            expected_mean = variance*(sign_mode*model%observed/model%sigma**2 + &
                model%prior_mean/model%prior_sd**2)
            expected_z = -0.5_dp*(log(2*acos(-1.0_dp)*total_variance) + &
                (model%observed - sign_mode*model%prior_mean)**2/ &
                total_variance)
            call expect(abs(g%mean(1) - expected_mean) < 1.0e-7_dp, &
                'analytic mode location')
            call expect(abs(g%cov(1, 1) - variance) < 1.0e-12_dp, &
                'analytic local covariance')
            call expect(abs(g%log_det_hessian + log(variance)) < 1.0e-12_dp, &
                'Hessian determinant')
            call expect(abs(g%log_evidence - expected_z) < 1.0e-10_dp, &
                'normalized analytic local evidence')
            call expect(m%psis%ess > 199 .and. m%psis%failures == 0, &
                'existing PSIS: Gaussian proposal exact')
            call expect(ieee_is_finite(m%psis%khat), 'PSIS k-hat reported')
            open (newunit=u, file=m%equilibrium_file, status='old')
            read (u, *) saved
            close (u)
            call expect(abs(saved - g%mean(1)) < 1.0e-12_dp, &
                'mode state saved at MAP after PSIS')
        end associate
    end do
    ratio = mixture%modes(positive)%weight/mixture%modes(negative)%weight
    call expect(abs(ratio - exp(2*model%observed*model%prior_mean/total_variance)) < &
        1.0e-9_dp, 'known relative evidence recovered')
    call expect(abs(sum(mixture%modes%weight) - 1) < 1.0e-14_dp, &
        'normalized finite mixture')
    do i = 5, size(starts)
        call expect(mixture%attempts(i)%mode > 0, 'prior draw converges to a mode')
    end do
    print '(A,F12.8)', 'two-mode evidence ratio: ', ratio
    call rng_seed(rng, 29_int64, status)
    do j = 1, 100
        call sample_local_mixture(mixture, rng, i, sample, status)
        call expect(status%code == 0, 'finite mixture sampling')
        call expect(abs(sample(1)) > 4 .and. &
            sample(1)*mixture%modes(i)%gaussian%mean(1) > 0, &
            'sample keeps its mode identity')
    end do

    ! Prior draws alone must also discover both signs (fixed seed).
    call multi_start_map(model, starts(5:), 1, trim(directory), opts, mixture, status)
    call expect(status%code == 0 .and. size(mixture%modes) == 2, &
        'prior-only mode discovery')

    block
        type(fixture_t) :: scaled
        type(mode_start_t) :: scaled_starts(2)
        scaled%observed = 0.06_dp
        scaled%sigma = 0.0025_dp
        scaled%prior_mean = 0.005_dp
        scaled%prior_sd = 0.03_dp
        scaled_starts(1)%theta = [-0.04_dp]
        scaled_starts(2)%theta = [0.04_dp]
        call multi_start_map(scaled, scaled_starts, 1, trim(directory), opts, &
            mixture, status)
        call expect(status%code == 0 .and. size(mixture%modes) == 2, &
            'Mahalanobis deduplication invariant under coordinate rescaling')
    end block

    model%folded = .false.
    call multi_start_map(model, starts, 1, trim(directory), opts, mixture, status)
    call expect(status%code == 0 .and. size(mixture%modes) == 1, &
        'unimodal model returns one mode')
    call expect(all(mixture%attempts%mode == 1), &
        'all starts deduplicated in unimodal case')
    expected_z = -0.5_dp*(log(2*acos(-1.0_dp)*total_variance) + &
        (model%observed - model%prior_mean)**2/total_variance)
    call expect(abs(mixture%log_evidence - expected_z) < 1.0e-10_dp, &
        'unimodal absolute evidence')

    ! Distinct model/branch identities are never merged, even with equal parameters.
    family(1)%theta = [4.0_dp]
    family(2)%theta = [5.0_dp]
    family%family = ['A', 'B']
    family%log_family_prior = log([0.25_dp, 0.75_dp])
    call multi_start_map(model, family, 1, trim(directory), opts, mixture, status, &
        prepare)
    call expect(status%code == 0 .and. size(mixture%modes) == 2, &
        'separate branch identities retained')
    ratio = exp(((model%observed - 0.5_dp)**2 - &
        (model%observed - 1.5_dp)**2)/(2*total_variance))*3
    call expect(abs(mixture%modes(2)%weight/mixture%modes(1)%weight - ratio) < &
        1.0e-9_dp, 'activated family likelihood and discrete prior')
    model%prior_mean = 0.5_dp

    model%proper = .false.
    call multi_start_map(model, starts(:4), 1, trim(directory), opts, mixture, status)
    call expect(status%code == 0 .and. .not. mixture%evidence_valid, &
        'improper priors flagged')
    call expect(.not. ieee_is_finite(mixture%modes(1)%weight), &
        'no evidence weights for improper priors')
    model%proper = .true.
    bad(1)%theta = [0.0_dp]
    call multi_start_map(model, bad, 1, trim(directory), opts, mixture, status)
    call expect(status%code /= 0 .and. size(mixture%modes) == 0, &
        'failed start excluded')
    call expect(.not. mixture%attempts(1)%converged, 'failed start reported')
    opts%max_iterations = 1
    call multi_start_map(model, starts(:1), 1, trim(directory), opts, mixture, status)
    call expect(status%code /= 0 .and. size(mixture%modes) == 0, &
        'iteration limit excluded')
    opts%max_iterations = 40
    model%bounded = .true.
    opts%decrease_tolerance = 1.0e-3_dp
    bad(1)%theta = [4.0_dp]
    call multi_start_map(model, bad, 1, trim(directory), opts, mixture, status)
    call expect(status%code /= 0 .and. size(mixture%modes) == 0, &
        'infeasible optimum is not convergence from large damping')
    if (fails /= 0) stop 1
contains
    subroutine expect(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (condition) return
        print '(A)', 'FAIL: '//label
        fails = fails + 1
    end subroutine
    subroutine prepare(index, posterior, ok)
        use tiago_local_posterior, only: local_posterior_t
        integer, intent(in) :: index
        class(local_posterior_t), intent(inout) :: posterior
        logical, intent(out) :: ok
        ok = .false.
        select type (posterior)
            type is (fixture_t)
            posterior%folded = .false.
            posterior%prior_mean = 0.5_dp + (index - 1)
            ok = index >= 1 .and. index <= 2
        end select
    end subroutine
end program
