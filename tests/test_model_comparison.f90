module comparison_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_local_posterior, only: local_posterior_t
    use fortnum_rng, only: rng_t, rng_normal
    implicit none
    type, extends(local_posterior_t) :: signature_t
        logical :: alternative = .false.
        real(dp) :: y(2) = [0.0_dp, 3.0_dp], sigma = 0.15_dp
        real(dp) :: profile_sd = 1, signature_sd = 1, signature_mean = 0.3_dp
        real(dp) :: sensitivity = 1, constant_residual = 0
    contains
        procedure :: residual
        procedure :: jacobian
        procedure :: n_residuals
        procedure :: log_normalization
        procedure :: proper_prior
        procedure :: draw_prior
    end type
contains
    subroutine residual(self, theta, r, ok)
        class(signature_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), intent(out) :: r(:)
        logical, intent(out) :: ok
        real(dp) :: signal
        signal = 0
        if (self%alternative) signal = self%sensitivity*abs(theta(2))
        r(1) = (theta(1) - self%y(1))/self%sigma
        r(2) = (theta(1) + signal - self%y(2))/self%sigma
        r(3) = theta(1)/self%profile_sd
        if (self%alternative) &
            r(4) = (theta(2) - self%signature_mean)/self%signature_sd
        if (self%constant_residual /= 0) r(size(r)) = self%constant_residual
        ok = .true.
    end subroutine
    subroutine jacobian(self, theta, jac)
        class(signature_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), allocatable, intent(out) :: jac(:, :)
        allocate (jac(self%n_residuals(), size(theta)))
        jac = 0
        jac(1:2, 1) = 1/self%sigma
        jac(3, 1) = 1/self%profile_sd
        if (self%alternative) then
            jac(2, 2) = self%sensitivity*sign(1.0_dp, theta(2))/self%sigma
            jac(4, 2) = 1/self%signature_sd
        end if
    end subroutine
    integer function n_residuals(self) result(n)
        class(signature_t), intent(in) :: self
        n = 3
        if (self%alternative) n = n + 1
        if (self%constant_residual /= 0) n = n + 1
    end function
    real(dp) function log_normalization(self) result(c)
        class(signature_t), intent(in) :: self
        c = -1.5_dp*log(2*acos(-1.0_dp)) &
            - 2*log(self%sigma) - log(self%profile_sd)
        if (self%alternative) c = c - 0.5_dp*log(2*acos(-1.0_dp)) &
            - log(self%signature_sd)
    end function
    logical function proper_prior(self) result(proper)
        class(signature_t), intent(in) :: self
        proper = .true.
    end function
    subroutine draw_prior(self, rng, theta, ok)
        class(signature_t), intent(in) :: self
        type(rng_t), intent(inout) :: rng
        real(dp), intent(out) :: theta(:)
        logical, intent(out) :: ok
        call rng_normal(rng, theta(1))
        theta(1) = self%profile_sd*theta(1)
        if (self%alternative) then
            call rng_normal(rng, theta(2))
            theta(2) = self%signature_mean + self%signature_sd*theta(2)
        end if
        ok = .true.
    end subroutine

    real(dp) function marginal_likelihood(model, amplitude) result(value)
        type(signature_t), intent(in) :: model
        real(dp), intent(in) :: amplitude
        real(dp) :: noise, prior, determinant, y0, y1, quadratic
        noise = model%sigma**2
        prior = model%profile_sd**2
        determinant = noise*(noise + 2*prior)
        y0 = model%y(1)
        y1 = model%y(2) - model%sensitivity*abs(amplitude)
        quadratic = ((noise + prior)*(y0**2 + y1**2) &
            - 2*prior*y0*y1)/determinant
        value = exp(-0.5_dp*quadratic)/(2*acos(-1.0_dp)*sqrt(determinant))
    end function

    subroutine integrate(model, panels, evidence, exceedance)
        type(signature_t), intent(in) :: model
        integer, intent(in) :: panels
        real(dp), intent(out) :: evidence, exceedance
        real(dp) :: a, h, density
        integer :: i
        evidence = 0
        exceedance = 0
        h = 20*model%signature_sd/panels
        do i = 1, panels
            a = model%signature_mean - 10*model%signature_sd + (i - 0.5_dp)*h
            density = marginal_likelihood(model, a) &
                *exp(-0.5_dp*((a - model%signature_mean)/model%signature_sd)**2) &
                /(sqrt(2*acos(-1.0_dp))*model%signature_sd)
            evidence = evidence + h*density
            if (abs(a) > 0.5_dp) exceedance = exceedance + h*density
        end do
        exceedance = exceedance/evidence
    end subroutine
end module

program test_model_comparison
    use comparison_fixture
    use tiago_local_mixture, only: mode_start_t, mixture_options_t, &
        local_mixture_t, multi_start_map, &
        compare_local_mixtures
    use tiago_inference_engines, only: map_fit
    use fortnum_status, only: fortnum_status_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    type(signature_t) :: model
    type(local_mixture_t) :: mixtures(2)
    type(mode_start_t) :: null_start(1), starts(4)
    type(mixture_options_t) :: options
    type(fortnum_status_t) :: status
    real(dp), allocatable :: probabilities(:)
    real(dp) :: z, reference, refined, exceedance, refined_exceedance, null_z
    real(dp) :: theta(1), expected, probability, local_exceedance, sd, mean
    character(len=1024) :: directory
    character(len=:), allocatable :: message
    logical :: converged
    integer :: i, j, steps

    call get_command_argument(1, directory)
    options%decrease_tolerance = 1.0e-12_dp
    options%psis_draws = 100
    options%seed = 1821
    null_start(1)%theta = [0.0_dp]
    do i = 1, 4
        allocate (starts(i)%theta(2))
        starts(i)%theta(1) = 0
        starts(i)%theta(2) = (-1.0_dp)**i*(1 + i*0.1_dp)
    end do
    do j = 1, 4
        model%sensitivity = 1
        select case (j)
        case (1)
            model%y(2) = 3
        case (2)
            model%y(2) = 0.7_dp
        case (3)
            model%y(2) = 0.05_dp
        case (4)
            model%sensitivity = 0
            model%y(2) = 0
        end select
        model%alternative = .false.
        call multi_start_map(model, null_start, 1, trim(directory)//'/null', &
            options, mixtures(1), status)
        call require(status%code == 0, 'null fit')
        null_z = marginal_likelihood(model, 0.0_dp)
        call require(abs(mixtures(1)%log_evidence - log(null_z)) < 1.e-10_dp, &
            'independent marginalized null evidence')
        model%alternative = .true.
        call multi_start_map(model, starts, 2, trim(directory)//'/alternative', &
            options, mixtures(2), status)
        call require(status%code == 0, 'alternative fit')
        call compare_local_mixtures(mixtures, log([0.7_dp, 0.3_dp]), &
            probabilities, z, status)
        call require(status%code == 0, 'different-dimensional comparison')
        call integrate(model, 40000, reference, exceedance)
        call integrate(model, 80000, refined, refined_exceedance)
        call require(abs(reference/refined - 1) < 1.e-6_dp, 'integration refinement')
        call require(abs(exceedance - refined_exceedance) < 1.e-5_dp, &
            'signature exceedance refinement')
        probability = 0.3_dp*refined/(0.7_dp*null_z + 0.3_dp*refined)
        local_exceedance = 0
        do i = 1, size(mixtures(2)%modes)
            mean = mixtures(2)%modes(i)%gaussian%mean(2)
            sd = sqrt(mixtures(2)%modes(i)%gaussian%cov(2, 2))
            local_exceedance = local_exceedance + mixtures(2)%modes(i)%weight &
                *0.5_dp*(erfc((0.5_dp - mean)/(sqrt(2.0_dp)*sd)) &
                + erfc((0.5_dp + mean)/(sqrt(2.0_dp)*sd)))
        end do
        print '(I2,1X,I2,6(1X,ES15.7))', j, size(mixtures(2)%modes), &
            mixtures(2)%log_evidence, log(refined), probabilities(2), probability, &
            probabilities(2)*local_exceedance, probability*refined_exceedance
        if (j == 1 .or. j == 4) then
            call require(abs(mixtures(2)%log_evidence - log(refined)) < 1.e-8_dp, &
                'Gaussian regime absolute evidence')
            call require(abs(probabilities(2) - probability) < 1.e-8_dp, &
                'Gaussian regime model probability')
            call require(abs(z - log(0.7_dp*null_z + 0.3_dp*refined)) < 1.e-8_dp, &
                'model-averaged evidence')
            call require(abs(probabilities(2)*local_exceedance &
                - probability*refined_exceedance) < 1.e-7_dp, &
                'model-averaged exceedance')
        end if
        if (j == 1 .or. j == 2) &
            call require(size(mixtures(2)%modes) == 2, 'both alias modes retained')
        if (j == 3) &
            call require(abs(mixtures(2)%log_evidence - log(refined)) > 0.1_dp, &
            'reference integration exposes boundary error')
        if (j == 4) &
            call require(abs(probabilities(2) - 0.3_dp) < 1.e-10_dp, &
            'insensitive observation retains model prior')
    end do
    call compare_local_mixtures(mixtures, [0.0_dp, 0.0_dp], probabilities, z, status)
    call require(status%code /= 0, 'unnormalized model priors rejected')
    mixtures(1)%evidence_valid = .false.
    call compare_local_mixtures(mixtures, log([0.7_dp, 0.3_dp]), &
        probabilities, z, status)
    call require(status%code /= 0, 'unqualified input evidence rejected')
    call require(.not. any(ieee_is_finite(probabilities)), 'no partial probabilities')

    model%alternative = .false.
    model%constant_residual = 1.e8_dp
    model%y = [2.0_dp, 2.0_dp]
    model%sigma = 1
    theta = 0
    expected = 4.0_dp/3
    call map_fit(model, theta, 100, steps, message, converged, 1.e-14_dp)
    call require(converged, 'constant residual does not stall MAP')
    call require(abs(theta(1) - expected) < 1.e-8_dp, 'analytic linear MAP')
contains
    subroutine require(ok, description)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: description
        if (.not. ok) then
            print *, 'FAIL: ', description
            stop 1
        end if
    end subroutine
end program
