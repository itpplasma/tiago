module tiago_local_mixture
    !! Finite, evidence-weighted local posteriors; no equilibrium averaging.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use fortnum_rng, only: rng_t, rng_seed, rng_normal, rng_uniform
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR
    use tiago_local_posterior, only: local_posterior_t
    use tiago_inference_engines, only: map_fit, laplace, laplace_t, psis_validate, &
        psis_result_t
    implicit none
    private

    type, public :: mode_start_t
        real(dp), allocatable :: theta(:) !! absent: draw from the active model prior
        character(len=64) :: family = 'default' !! model and equilibrium branch
        real(dp) :: log_family_prior = 0.0_dp !! discrete model/branch probability
    end type

    type, public :: mixture_options_t
        integer :: max_iterations = 40, psis_draws = 400, seed = 1
        real(dp) :: dedup_distance = 1.0_dp, decrease_tolerance = 1.0e-3_dp
    end type

    type, public :: local_mode_t
        type(laplace_t) :: gaussian
        type(psis_result_t) :: psis
        character(len=64) :: family = ''
        character(len=:), allocatable :: equilibrium_file, message
        real(dp) :: weight = 0, log_family_prior = 0, log_normalization = 0
        real(dp), allocatable :: residual(:)
        real(dp) :: gradient_max = 0
        integer :: start = 0, steps = 0, residuals = 0
    end type

    type, public :: mode_attempt_t
        integer :: mode = 0, steps = 0
        logical :: converged = .false.
        character(len=:), allocatable :: message
    end type

    type, public :: local_mixture_t
        type(local_mode_t), allocatable :: modes(:)
        type(mode_attempt_t), allocatable :: attempts(:)
        logical :: evidence_valid = .false.
        real(dp) :: log_evidence = 0
    end type

    abstract interface
        subroutine prepare_interface(index, model, ok)
            import local_posterior_t
            integer, intent(in) :: index
            class(local_posterior_t), intent(inout) :: model
            logical, intent(out) :: ok
        end subroutine
    end interface
    public :: multi_start_map, sample_local_mixture
contains
    subroutine multi_start_map(model, starts, n, directory, opts, mixture, status, &
            prepare)
        !! prepare activates the supplied model/branch BEFORE every start.
        !! directory must exist; its mode files are overwritten on a repeated run.
        class(local_posterior_t), intent(inout) :: model
        type(mode_start_t), intent(in) :: starts(:)
        integer, intent(in) :: n
        character(len=*), intent(in) :: directory
        type(mixture_options_t), intent(in) :: opts
        type(local_mixture_t), intent(out) :: mixture
        type(fortnum_status_t), intent(out) :: status
        procedure(prepare_interface), optional :: prepare
        type(local_mode_t) :: trial
        type(rng_t) :: rng
        real(dp), allocatable :: theta(:), r(:), jac(:, :), logw(:)
        real(dp) :: nan, peak
        integer :: i, j, match, m, unit, ios
        logical :: ok, converged
        character(len=32) :: number
        character(len=:), allocatable :: prefix, message

        call status_set(status, FORTNUM_DOMAIN_ERROR, &
            'invalid mixture options or starts')
        allocate (mixture%modes(0), mixture%attempts(size(starts)))
        nan = ieee_value(1.0_dp, ieee_quiet_nan)
        mixture%log_evidence = nan
        if (n < 1 .or. size(starts) < 1 .or. opts%max_iterations < 1 .or. &
            opts%psis_draws < 20 .or. opts%dedup_distance < 0 .or. &
            opts%decrease_tolerance <= 0) return
        if (.not. ieee_is_finite(opts%dedup_distance) .or. &
            .not. ieee_is_finite(opts%decrease_tolerance)) return
        do i = 1, size(starts)
            if (len_trim(starts(i)%family) == 0) return
            if (scan(starts(i)%family, ',"'//achar(10)//achar(13)) /= 0) return
            if (.not. ieee_is_finite(starts(i)%log_family_prior)) return
            if (allocated(starts(i)%theta)) then
                if (size(starts(i)%theta) /= n) return
                if (.not. all(ieee_is_finite(starts(i)%theta))) return
            end if
            do j = 1, i - 1
                if (starts(j)%family /= starts(i)%family) cycle
                if (starts(j)%log_family_prior /= starts(i)%log_family_prior) return
            end do
        end do
        ! Discrete family probabilities are counted once, never once per start.
        peak = maxval(starts%log_family_prior)
        m = 0
        allocate (logw(size(starts)))
        do i = 1, size(starts)
            if (any(starts(:i - 1)%family == starts(i)%family)) cycle
            m = m + 1
            logw(m) = starts(i)%log_family_prior
        end do
        if (abs(peak + log(sum(exp(logw(:m) - peak)))) > 1.0e-10_dp) return
        deallocate (logw)
        open (newunit=unit, file=directory//'/modes.csv', status='replace', iostat=ios)
        if (ios /= 0) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'cannot write mixture directory')
            return
        end if
        close (unit)
        call rng_seed(rng, int(opts%seed, int64), status)
        if (status%code /= FORTNUM_OK) return
        allocate (theta(n))
        do i = 1, size(starts)
            mixture%attempts(i)%message = 'candidate activation failed'
            ok = .true.
            if (present(prepare)) call prepare(i, model, ok)
            if (.not. ok) cycle
            call model%exact()
            if (allocated(starts(i)%theta)) then
                theta = starts(i)%theta
            else
                call model%draw_prior(rng, theta, ok)
                mixture%attempts(i)%message = 'prior draw requires a proper prior'
                if (.not. ok) cycle
            end if
            call map_fit(model, theta, opts%max_iterations, mixture%attempts(i)%steps, &
                message, converged, opts%decrease_tolerance)
            mixture%attempts(i)%message = message
            mixture%attempts(i)%converged = converged
            if (.not. converged) cycle
            trial = local_mode_t()
            call laplace(model, theta, trial%gaussian)
            if (trial%gaussian%info /= 0) then
                mixture%attempts(i)%message = 'Laplace factorization failed'
                cycle
            end if
            if (.not. ieee_is_finite(trial%gaussian%log_evidence)) then
                mixture%attempts(i)%message = 'nonfinite log evidence'
                cycle
            end if
            trial%family = starts(i)%family
            trial%log_family_prior = starts(i)%log_family_prior
            match = 0
            do j = 1, size(mixture%modes)
                if (trial%family /= mixture%modes(j)%family) cycle
                if (.not. same_mode(trial%gaussian, mixture%modes(j)%gaussian, &
                    opts%dedup_distance)) cycle
                match = j
                exit
            end do
            if (match /= 0) then
                mixture%attempts(i)%mode = match
                if (trial%gaussian%chi2 >= mixture%modes(match)%gaussian%chi2) cycle
            else
                mixture%modes = [mixture%modes, trial]
                match = size(mixture%modes)
                mixture%attempts(i)%mode = match
            end if
            trial%start = i
            trial%steps = mixture%attempts(i)%steps
            trial%message = message
            trial%log_normalization = model%log_normalization()
            trial%residuals = model%n_residuals()
            call model%jacobian(theta, jac)
            if (allocated(r)) deallocate (r)
            allocate (r(trial%residuals))
            call model%residual(theta, r, ok)
            trial%gradient_max = maxval(abs(matmul(r, jac)))
            call psis_validate(model, trial%gaussian, opts%psis_draws, &
                opts%seed + i, trial%psis)
            ! PSIS leaves the provider at its last draw. Re-solve the MAP to save it.
            call model%exact()
            call model%residual(theta, r, ok)
            if (.not. ok) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    'could not restore mode equilibrium')
                return
            end if
            trial%residual = r
            write (number, '(I0)') match
            prefix = directory//'/mode_'//trim(number)
            call model%save_equilibrium(prefix, trial%equilibrium_file)
            mixture%modes(match) = trial
        end do
        mixture%evidence_valid = size(mixture%modes) > 0
        mixture%evidence_valid = mixture%evidence_valid .and. &
            all(mixture%modes%gaussian%proper)
        if (mixture%evidence_valid) then
            logw = [(mixture%modes(j)%gaussian%log_evidence + &
                mixture%modes(j)%log_family_prior, j=1, size(mixture%modes))]
            peak = maxval(logw)
            mixture%log_evidence = peak + log(sum(exp(logw - peak)))
            mixture%modes%weight = exp(logw - mixture%log_evidence)
        else
            mixture%modes%weight = nan
        end if
        call write_reports(directory, mixture, ios)
        if (ios /= 0) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'cannot write mode reports')
        else if (size(mixture%modes) == 0) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'no converged modes')
        else
            call status_set(status, FORTNUM_OK, '')
        end if
    end subroutine

    subroutine sample_local_mixture(mixture, rng, mode, theta, status)
        !! The returned mode selects the equilibrium/model to evaluate at theta.
        type(local_mixture_t), intent(in) :: mixture
        type(rng_t), intent(inout) :: rng
        integer, intent(out) :: mode
        real(dp), allocatable, intent(out) :: theta(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: z(:)
        real(dp) :: uniform, cumulative
        integer :: i
        mode = 0
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'mixture has no evidence weights')
        if (.not. mixture%evidence_valid) return
        call rng_uniform(rng, uniform)
        cumulative = 0
        do i = 1, size(mixture%modes)
            cumulative = cumulative + mixture%modes(i)%weight
            mode = i
            if (uniform < cumulative) exit
        end do
        associate (g => mixture%modes(mode)%gaussian)
            allocate (z(size(g%mean)))
            do i = 1, size(z)
                call rng_normal(rng, z(i))
            end do
            theta = g%mean + matmul(g%chol, z)
        end associate
        call status_set(status, FORTNUM_OK, '')
    end subroutine

    logical function same_mode(a, b, distance) result(same)
        type(laplace_t), intent(in) :: a, b
        real(dp), intent(in) :: distance
        real(dp) :: delta(size(a%mean)), za(size(a%mean)), zb(size(a%mean))
        integer :: i
        delta = a%mean - b%mean
        do i = 1, size(delta)
            za(i) = (delta(i) - dot_product(a%chol(i, :i - 1), za(:i - 1)))/a%chol(i, i)
            zb(i) = (delta(i) - dot_product(b%chol(i, :i - 1), zb(:i - 1)))/b%chol(i, i)
        end do
        same = max(sum(za**2), sum(zb**2)) <= distance**2
    end function

    subroutine write_reports(directory, mixture, ios)
        character(len=*), intent(in) :: directory
        type(local_mixture_t), intent(in) :: mixture
        integer, intent(out) :: ios
        integer :: u, i, j, k
        character(len=32) :: number
        character(len=:), allocatable :: prefix
        open (newunit=u, file=directory//'/modes.csv', status='replace', iostat=ios)
        if (ios /= 0) return
        write (u, '(A)') &
            'mode,family,start,weight,log_evidence,log_family_prior,'// &
            'log_det_hessian,'// &
            'log_normalization,proper,chi2,residuals,gradient_max,steps,khat,'// &
            'psis_ess,psis_failures,equilibrium'
        do i = 1, size(mixture%modes)
            associate (m => mixture%modes(i), g => mixture%modes(i)%gaussian)
                write (u, '(*(G0,:,","))') i, trim(m%family), m%start, m%weight, &
                    g%log_evidence, &
                    m%log_family_prior, g%log_det_hessian, m%log_normalization, &
                    g%proper, &
                    g%chi2, m%residuals, m%gradient_max, m%steps, m%psis%khat, &
                    m%psis%ess, &
                    m%psis%failures, '"'//m%equilibrium_file//'"'
            end associate
        end do
        close (u, iostat=ios)
        if (ios /= 0) return
        open (newunit=u, file=directory//'/starts.csv', status='replace', iostat=ios)
        if (ios /= 0) return
        write (u, '(A)') 'start,mode,steps,converged,message'
        do i = 1, size(mixture%attempts)
            associate (a => mixture%attempts(i))
                write (u, '(*(G0,:,","))') i, a%mode, a%steps, a%converged, a%message
            end associate
        end do
        close (u, iostat=ios)
        if (ios /= 0) return
        do i = 1, size(mixture%modes)
            write (number, '(I0)') i
            prefix = directory//'/mode_'//trim(number)
            open (newunit=u, file=prefix//'_parameters.csv', status='replace', &
                iostat=ios)
            if (ios /= 0) return
            write (u, '(A)') 'parameter,map,sd'
            associate (g => mixture%modes(i)%gaussian)
                do j = 1, size(g%mean)
                    write (u, '(*(G0,:,","))') j, g%mean(j), sqrt(g%cov(j, j))
                end do
                close (u, iostat=ios)
                if (ios /= 0) return
                open (newunit=u, file=prefix//'_covariance.csv', status='replace', &
                    iostat=ios)
                if (ios /= 0) return
                do j = 1, size(g%mean)
                    write (u, '(*(G0,:,","))') (g%cov(j, k), k=1, size(g%mean))
                end do
                close (u, iostat=ios)
                if (ios /= 0) return
            end associate
            open (newunit=u, file=prefix//'_residuals.csv', status='replace', &
                iostat=ios)
            if (ios /= 0) return
            write (u, '(A)') 'row,whitened_residual'
            do j = 1, size(mixture%modes(i)%residual)
                write (u, '(*(G0,:,","))') j, mixture%modes(i)%residual(j)
            end do
            close (u, iostat=ios)
            if (ios /= 0) return
        end do
    end subroutine
end module
