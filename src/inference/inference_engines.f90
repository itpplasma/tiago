module tiago_inference_engines
    !! Inference engines on a joint_model_t:
    !!   map_fit         Levenberg-Marquardt on the whitened residual (the
    !!                   algorithm of tiago_reconstruction's fit; VMEC++
    !!                   hot-restarts from the last equilibrium);
    !!   laplace         local covariance and normalized Laplace evidence;
    !!                   observed curvature for hierarchical AR(1) profiles,
    !!                   Gauss-Newton curvature otherwise;
    !!   sample_nuts     fortnum NUTS with the Laplace covariance as metric and
    !!                   gradients -J^T r from the implicit-adjoint chain;
    !!   psis_validate   importance sampling from the Laplace Gaussian with
    !!                   Pareto smoothing (k-hat) and the importance-weighted
    !!                   posterior moments.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    use fortnum_rng, only: rng_t, rng_seed, rng_normal
    use fortnum_nuts, only: nuts_options_t, nuts_result_t, nuts_sample
    use fortnum_psis, only: psis_smooth
    use tiago_joint_model, only: joint_model_t
    implicit none
    private

    type :: context_t
        type(joint_model_t), pointer :: model => null()
    end type context_t

    type, public :: laplace_t
        real(dp), allocatable :: mean(:), cov(:, :), chol(:, :)
        real(dp) :: log_evidence = 0.0_dp, log_det_hessian = 0.0_dp, chi2 = 0.0_dp
        logical :: proper = .false.
        integer :: info = 0
    end type laplace_t

    type, public :: psis_result_t
        real(dp) :: khat = huge(1.0_dp)
        integer :: draws = 0, failures = 0
        real(dp), allocatable :: mean(:), sd(:)
        real(dp) :: ess = 0.0_dp
    end type psis_result_t

    public :: map_fit, laplace, sample_nuts, psis_validate, log_posterior

contains

    subroutine map_fit(model, theta, max_iterations, steps, message)
        type(joint_model_t), intent(inout) :: model
        real(dp), intent(inout) :: theta(:)
        integer, intent(in) :: max_iterations
        integer, intent(out) :: steps
        character(len=:), allocatable, intent(out) :: message
        real(dp) :: floor
        real(dp), allocatable :: r(:), r_trial(:), jac(:, :), a(:, :), m(:, :), g(:), h(:), d(:)
        real(dp), allocatable :: theta_trial(:)
        real(dp) :: chi2, chi2_trial, predicted, rho, mu, nu
        integer :: n, i, info
        logical :: ok
        external :: dposv

        n = size(theta)
        floor = 1.e-3_dp
        if (model%prof%te_ar1%active) floor = 1.e-10_dp
        allocate(r(model%n_residuals()), r_trial(model%n_residuals()))
        call model%residual(theta, r, ok)
        steps = 0
        if (.not. ok) then
            message = 'no equilibrium at the start values'
            return
        end if
        chi2 = sum(r**2)
        call model%jacobian(theta, jac)
        a = matmul(transpose(jac), jac)
        g = matmul(transpose(jac), r)
        mu = 1.0e-3_dp
        nu = 2.0_dp
        message = 'iteration limit reached'
        do while (steps < max_iterations)
            d = [(max(a(i, i), tiny(1.0_dp)), i = 1, n)]
            m = a
            do i = 1, n
                m(i, i) = m(i, i) + mu * d(i)
            end do
            h = -g
            call dposv('U', n, 1, m, n, h, n, info)
            if (info /= 0) then
                mu = mu * nu
                nu = 2.0_dp * nu
                cycle
            end if
            predicted = -2.0_dp * dot_product(h, g) - dot_product(h, matmul(a, h))
            if (predicted <= floor) then
                message = 'converged (predicted chi^2 decrease below the solver noise floor)'
                exit
            end if
            theta_trial = theta + h
            call model%residual(theta_trial, r_trial, ok)
            chi2_trial = huge(1.0_dp)
            if (ok) chi2_trial = sum(r_trial**2)
            rho = -1.0_dp
            if (ok) rho = (chi2 - chi2_trial) / predicted
            if (rho > 0.0_dp) then
                theta = theta_trial
                r = r_trial
                chi2 = chi2_trial
                steps = steps + 1
                call model%jacobian(theta, jac)
                a = matmul(transpose(jac), jac)
                g = matmul(transpose(jac), r)
                mu = mu * max(1.0_dp / 3.0_dp, 1.0_dp - (2.0_dp * rho - 1.0_dp)**3)
                nu = 2.0_dp
            else
                mu = mu * nu
                nu = 2.0_dp * nu
                if (mu > 1.0e12_dp) then
                    message = 'converged (no further decrease of chi^2)'
                    exit
                end if
            end if
        end do
    end subroutine map_fit

    subroutine laplace(model, theta, lap)
        type(joint_model_t), intent(inout) :: model
        real(dp), intent(in) :: theta(:)
        type(laplace_t), intent(out) :: lap
        real(dp), allocatable :: r(:), jac(:, :), h(:, :)
        integer :: n, i, j
        logical :: ok
        external :: dpotrf, dpotri

        n = size(theta)
        allocate(r(model%n_residuals()))
        call model%residual(theta, r, ok)
        if (.not. ok) then
            lap%info = -1
            return
        end if
        call model%jacobian(theta, jac)
        h = matmul(transpose(jac), jac)
        if (model%prof%te_ar1%active) then
            call observed_hessian(model,theta,h,ok)
            if (.not. ok) then
                lap%info = -2
                return
            end if
        end if
        call dpotrf('U', n, h, n, lap%info)
        if (lap%info /= 0) return
        lap%log_det_hessian = 2.0_dp * sum([(log(h(i, i)), i = 1, n)])
        lap%cov = h
        call dpotri('U', n, lap%cov, n, lap%info)
        do j = 1, n
            do i = j + 1, n
                lap%cov(i, j) = lap%cov(j, i)
            end do
        end do
        lap%chol = lap%cov
        call dpotrf('L', n, lap%chol, n, lap%info)
        do j = 2, n
            lap%chol(1:j - 1, j) = 0.0_dp
        end do
        lap%mean = theta
        lap%chi2 = sum(r**2)
        lap%proper = model%proper_prior()
        lap%log_evidence = -0.5_dp * lap%chi2 + model%log_normalization() + &
            0.5_dp * n * log(2.0_dp * acos(-1.0_dp)) - 0.5_dp * lap%log_det_hessian
    end subroutine laplace

    subroutine observed_hessian(model, theta, h, ok)
        !! Hyperparameters require residual curvature beyond Gauss-Newton.
        type(joint_model_t), intent(inout) :: model
        real(dp), intent(in) :: theta(:)
        real(dp), intent(out) :: h(:, :)
        logical, intent(out) :: ok
        real(dp) :: shifted(size(theta)), plus(size(theta)), minus(size(theta)), step
        real(dp), allocatable :: r(:), jac(:, :)
        integer :: i
        allocate (r(model%n_residuals()))
        do i = 1, size(theta)
            step = 2.e-5_dp*max(1._dp, abs(theta(i)))
            shifted = theta; shifted(i) = shifted(i) + step
            call model%residual(shifted, r, ok)
            if (.not. ok) return
            call model%jacobian(shifted, jac)
            plus = matmul(r, jac)
            shifted(i) = shifted(i) - 2*step
            call model%residual(shifted, r, ok)
            if (.not. ok) return
            call model%jacobian(shifted, jac)
            minus = matmul(r, jac)
            h(:, i) = (plus - minus)/(2*step)
        end do
        h = (h + transpose(h))/2
        call model%residual(theta, r, ok)
    end subroutine

    subroutine log_posterior(n, x, f, g, context, status)
        !! fortnum NUTS callback: f = -|r|^2 / 2, g = -J^T r.
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(out) :: f, g(n)
        class(*), intent(inout) :: context
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: r(:), jac(:, :)
        logical :: ok

        call status_set(status, FORTNUM_OK, '')
        f = 0.0_dp
        g = 0.0_dp
        select type (c => context)
            type is (context_t)
            allocate(r(c%model%n_residuals()))
            call c%model%residual(x, r, ok)
            if (.not. ok) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, 'no equilibrium')
                return
            end if
            call c%model%jacobian(x, jac)
            f = -0.5_dp * sum(r**2)
            g = -matmul(r, jac)
        end select
    end subroutine log_posterior

    subroutine sample_nuts(model, lap, opts, seed, result, status)
        type(joint_model_t), intent(inout), target :: model
        type(laplace_t), intent(in) :: lap
        type(nuts_options_t), intent(in) :: opts
        integer, intent(in) :: seed
        type(nuts_result_t), intent(out) :: result
        type(fortnum_status_t), intent(out) :: status
        type(context_t) :: context
        type(rng_t) :: rng

        context%model => model
        call rng_seed(rng, int(seed, int64), status)
        if (status%code /= FORTNUM_OK) return
        call nuts_sample(log_posterior, size(lap%mean), lap%mean, lap%chol, opts, rng, context, &
            result, status)
    end subroutine sample_nuts

    subroutine psis_validate(model, lap, draws, seed, result, samples)
        !! Importance sampling with the Laplace Gaussian q as proposal:
        !! log w = log p(theta | d) - log q(theta) up to a constant. Draws for
        !! which VMEC++ finds no equilibrium get the smallest finite weight
        !! minus 50 and are counted in failures.
        type(joint_model_t), intent(inout) :: model
        type(laplace_t), intent(in) :: lap
        integer, intent(in) :: draws, seed
        type(psis_result_t), intent(out) :: result
        real(dp), allocatable, intent(out), optional :: samples(:, :)
        type(rng_t) :: rng
        type(fortnum_status_t) :: st
        real(dp), allocatable :: zeta(:), theta(:, :), lw(:), sm(:), w(:), r(:)
        logical, allocatable :: good(:)
        integer :: k, i, n
        logical :: ok

        n = size(lap%mean)
        allocate(zeta(n), theta(n, draws), lw(draws), sm(draws), good(draws), &
            r(model%n_residuals()))
        call rng_seed(rng, int(seed, int64), st)
        do k = 1, draws
            do i = 1, n
                call rng_normal(rng, zeta(i))
            end do
            theta(:, k) = lap%mean + matmul(lap%chol, zeta)
            call model%residual(theta(:, k), r, ok)
            good(k) = ok
            lw(k) = 0.0_dp
            if (ok) lw(k) = -0.5_dp * sum(r**2) + 0.5_dp * dot_product(zeta, zeta)
        end do
        result%draws = draws
        result%failures = count(.not. good)
        if (result%failures == draws) return
        where (.not. good) lw = minval(lw, mask=good) - 50.0_dp
        call psis_smooth(lw, sm, result%khat, st)
        if (st%code /= FORTNUM_OK) return
        w = exp(sm - maxval(sm))
        w = w / sum(w)
        result%ess = 1.0_dp / sum(w**2)
        result%mean = matmul(theta, w)
        allocate(result%sd(n))
        do i = 1, n
            result%sd(i) = sqrt(max(sum(w * (theta(i, :) - result%mean(i))**2), 0.0_dp))
        end do
        if (present(samples)) samples = theta
    end subroutine psis_validate

end module tiago_inference_engines
