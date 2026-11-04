program tiago_infer
    !! Synthetic nested joint-inference study on a fixed-boundary VMEC case:
    !! truth (drawn from the prior, or the prior mean) -> synthetic magnetics
    !! and profile probes with noise -> MAP -> Laplace covariance and evidence
    !! -> optional PSIS check and NUTS -> CSV outputs for coverage statistics
    !! and figures (tools/inference_plots.py, tools/inference_coverage.py).
    !! Usage: tiago_infer input.nml
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use tiago_reconstruction, only: rec
    use tiago_joint_model, only: joint_model_t
    use tiago_synthetic_case, only: setup_magnetics, place_probes, synthesize_magnetics, &
        synthesize_profiles
    use tiago_kinetic_profiles, only: P_NE, P_TE, P_TI, P_ZEFF, N_PROFILES
    use tiago_inference_engines, only: map_fit, laplace, laplace_t, sample_nuts, psis_validate, &
        psis_result_t
    use tiago_local_mixture, only: mode_start_t, mixture_options_t, local_mixture_t, multi_start_map
    use fortnum_rng, only: rng_t, rng_seed, rng_normal, rng_split
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_nuts, only: nuts_options_t, nuts_result_t
    implicit none

    integer, parameter :: maxk = 32, ns_out = 41
    character(len=512) :: vmec_input = '', sensors_awk = '', output_dir = 'infer'
    integer :: seed = 1, max_iterations = 30, n_warmup = 200, n_samples = 300, max_depth = 6
    integer :: psis_draws = 0, posterior_solves = 40, laplace_draws = 400
    integer :: prior_starts = 0, candidate_count = 0
    character(len=512) :: candidate_file = ''
    real(dp) :: mode_distance = 1.0_dp
    logical :: noise = .true., truth_from_prior = .true.
    real(dp) :: breakpoints(maxk) = -1.0_dp
    real(dp) :: sigma_f = 0.3_dp, l_core = 0.4_dp, l_edge = 0.2_dp, s0 = 0.6_dp, sw = 0.2_dp
    real(dp) :: ne0 = 1.0e20_dp, te0 = 2.5e3_dp, ti0 = 2.0e3_dp, zeff = 1.5_dp
    real(dp) :: probe_rho(maxk) = -1.0_dp, probe_theta = 0.3_dp, align_sigma = 0.005_dp
    real(dp) :: sigma_profile_rel = 0.05_dp, vmec_ftol = 1.0e-14_dp
    real(dp) :: curtor_prior(2) = [-174250.0_dp, 2.0e4_dp], phiedge_prior(2) = [0.514386_dp, 0.02_dp]
    namelist /infer/ vmec_input, sensors_awk, output_dir, seed, max_iterations, n_warmup, n_samples, &
        max_depth, psis_draws, posterior_solves, laplace_draws, noise, truth_from_prior, breakpoints, &
        sigma_f, l_core, l_edge, s0, sw, ne0, te0, ti0, zeff, probe_rho, probe_theta, align_sigma, &
        sigma_profile_rel, vmec_ftol, curtor_prior, phiedge_prior, prior_starts, candidate_count, &
        candidate_file, mode_distance

    type(joint_model_t), target :: model
    type(local_mixture_t) :: mixture
    type(mixture_options_t) :: mixture_opts
    type(mode_start_t), allocatable :: starts(:)
    integer :: selected_mode
    type(laplace_t) :: lap
    type(psis_result_t) :: ps
    type(nuts_result_t) :: chain
    type(nuts_options_t) :: opts
    type(rng_t) :: rng, rng_truth, rng_noise
    type(fortnum_status_t) :: st
    character(len=512) :: path
    character(len=:), allocatable :: message, out
    character(len=32), allocatable :: names(:)
    real(dp), allocatable :: theta0(:), theta_true(:), theta(:), r(:), s_out(:), zeta(:)
    real(dp) :: t0, t1, t_map, t_lap, t_psis, t_nuts, chi2_start
    integer :: unit, ios, i, k, steps, tries
    logical :: ok

    call get_command_argument(1, path)
    open(newunit=unit, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'cannot open the namelist file'
    read(unit, nml=infer, iostat=ios)
    if (ios /= 0) error stop 'error in namelist &infer'
    close(unit)
    out = trim(output_dir)
    call execute_command_line('mkdir -p "'//out//'/sensors"')
    call execute_command_line('awk -v R0=1.42 -v rs=0.85 -v nphi=4 -v nth=4 -v nprobe=12 -v dir="'// &
        out//'/sensors" -f "'//trim(sensors_awk)//'"')
    call setup_magnetics(trim(vmec_input), out//'/work', out//'/sensors', &
        [character(len=64) :: 'CURTOR', 'PHIEDGE'], [curtor_prior(1), phiedge_prior(1)], &
        [curtor_prior(2), phiedge_prior(2)], vmec_ftol)
    if (count(breakpoints >= 0.0_dp) < 2) breakpoints(1:5) = [0.0_dp, 0.2_dp, 0.45_dp, 0.7_dp, 1.0_dp]
    call model%prof%init(pack(breakpoints, breakpoints >= 0.0_dp))
    associate (g => model%prof%greville)
        call model%prof%set_profile(P_NE, ne0 * (1.0_dp - 0.6_dp * g), .true., sigma_f, l_core, l_edge, s0, sw)
        call model%prof%set_profile(P_TE, te0 * (1.0_dp - g) + 50.0_dp, .true., sigma_f, l_core, l_edge, s0, sw)
        call model%prof%set_profile(P_TI, ti0 * (1.0_dp - g) + 50.0_dp, .true., sigma_f, l_core, l_edge, s0, sw)
        call model%prof%set_profile(P_ZEFF, zeff + 0.0_dp * g, .false.)
    end associate
    allocate(model%diags(0))
    call model%setup()
    theta0 = [rec%prior / rec%scale, [(0.0_dp, i = 1, model%n_prof)]]
    allocate(r(model%n_residuals()))
    call model%residual(theta0, r, ok)
    if (.not. ok) error stop 'no equilibrium at the prior mean'
    if (count(probe_rho > 0.0_dp) == 0) probe_rho(1:4) = [0.25_dp, 0.45_dp, 0.65_dp, 0.85_dp]
    call place_probes(model, pack(probe_rho, probe_rho > 0.0_dp), probe_theta, 0.0_dp, &
        acos(-1.0_dp) / 6.0_dp, align_sigma)
    call model%setup()
    deallocate(r)
    allocate(r(model%n_residuals()), zeta(model%n))
    names = parameter_names()
    theta0 = [rec%prior / rec%scale, [(0.0_dp, i = 1, model%n - model%n_rec)]]

    ! truth: a prior draw (rec parameters from their Gaussian priors, latents
    ! standard normal), redrawn if VMEC++ finds no equilibrium
    call rng_seed(rng, int(seed, int64), st)
    call rng_split(rng, 1_int64, rng_truth, st)
    call rng_split(rng, 2_int64, rng_noise, st)
    theta_true = theta0
    do tries = 1, 20
        if (truth_from_prior) then
            do i = 1, model%n
                call rng_normal(rng_truth, zeta(i))
            end do
            theta_true = theta0
            do k = 1, rec%n
                theta_true(k) = theta0(k) + zeta(k) * rec%prior_sigma(k) / rec%scale(k)
            end do
            theta_true(model%n_rec + 1:) = zeta(model%n_rec + 1:)
        end if
        call model%residual(theta_true, r, ok)
        if (ok) exit
    end do
    if (.not. ok) error stop 'no equilibrium for 20 prior draws'
    call synthesize_magnetics(model, theta_true(1:rec%n) * rec%scale, 0.01_dp, 1.0e-4_dp, 1.0e-8_dp, &
        1.0e-7_dp, rng_noise, noise)
    call synthesize_profiles(model, theta_true, sigma_profile_rel, [0.0_dp, 0.0_dp], rng_noise, noise)
    call rec%eq%write_wout(out//'/wout_truth.nc')
    s_out = [((i - 1.0_dp) / (ns_out - 1), i = 1, ns_out)]
    call write_profile_truth(out//'/truth_profiles.csv')
    call write_surfaces(out//'/surfaces_truth.csv', 0)

    ! MAP from the prior mean
    theta = theta0
    call model%residual(theta, r, ok)
    chi2_start = sum(r**2)
    t0 = wall()
    if (prior_starts < 0 .or. candidate_count < 0) error stop 'negative mixture start count'
    if (prior_starts + candidate_count > 0) then
        call execute_command_line('mkdir -p "'//out//'/modes"')
        allocate(starts(1 + candidate_count + prior_starts))
        starts(1)%theta = theta0
        if (candidate_count > 0) then
            open(newunit=unit, file=trim(candidate_file), status='old', iostat=ios)
            if (ios /= 0) error stop 'cannot open candidate_file'
            do i = 1, candidate_count
                allocate(starts(i + 1)%theta(model%n))
                read(unit, *, iostat=ios) starts(i + 1)%theta
                if (ios /= 0) error stop 'candidate_file requires one complete latent vector per line'
            end do
            close(unit)
        end if
        mixture_opts%max_iterations = max_iterations
        mixture_opts%psis_draws = max(20, psis_draws)
        mixture_opts%seed = seed
        mixture_opts%dedup_distance = mode_distance
        call multi_start_map(model, starts, model%n, out//'/modes', mixture_opts, mixture, st)
        if (st%code /= FORTNUM_OK) error stop 'multi-start MAP failed: '//trim(st%msg)
        selected_mode = minloc(mixture%modes%gaussian%chi2, dim=1)
        if (mixture%evidence_valid) selected_mode = maxloc(mixture%modes%weight, dim=1)
        theta = mixture%modes(selected_mode)%gaussian%mean
        steps = mixture%modes(selected_mode)%steps
        message = mixture%modes(selected_mode)%message
        call model%exact()
    else
        call map_fit(model, theta, max_iterations, steps, message)
    end if
    t1 = wall()
    t_map = t1 - t0
    call laplace(model, theta, lap)
    t0 = wall()
    t_lap = t0 - t1
    if (lap%info /= 0) error stop 'Laplace Hessian not positive definite'
    call rec%eq%write_wout(out//'/wout_map.nc')

    t_psis = 0.0_dp
    if (psis_draws > 0) then
        t0 = wall()
        call psis_validate(model, lap, psis_draws, seed + 7919, ps)
        t1 = wall()
        t_psis = t1 - t0
    end if
    t_nuts = 0.0_dp
    if (n_samples > 0) then
        opts%n_warmup = n_warmup
        opts%n_samples = n_samples
        opts%max_depth = max_depth
        t0 = wall()
        call sample_nuts(model, lap, opts, seed + 104729, chain, st)
        t1 = wall()
        t_nuts = t1 - t0
        if (st%code /= FORTNUM_OK) error stop 'NUTS failed: '//trim(st%msg)
        call write_draws(out//'/nuts_draws.csv')
    end if
    call write_parameters(out//'/parameters.csv')
    call write_profiles(out//'/profiles.csv')
    call write_posterior_solves()
    call write_summary(out//'/summary.txt')
    print '(A)', 'inference written to '//out

contains

    real(dp) function wall()
        integer(int64) :: count, rate
        call system_clock(count, rate)
        wall = real(count, dp) / real(rate, dp)
    end function wall

    function parameter_names() result(nm)
        character(len=32), allocatable :: nm(:)
        character(len=2), parameter :: tags(3) = ['ne', 'Te', 'Ti']
        integer :: j, q
        allocate(nm(model%n))
        do j = 1, rec%n
            nm(j) = rec%names(j)
        end do
        do q = 1, 3
            do j = 1, model%prof%ncoef
                write(nm(rec%n + (q - 1) * model%prof%ncoef + j), '(A,"_z",I0)') tags(q), j
            end do
        end do
        do j = model%n_rec + model%n_prof + 1, model%n
            write(nm(j), '("align_",I0)') j - model%n_rec - model%n_prof
        end do
    end function parameter_names

    subroutine profiles_of(th, f)
        !! ne, Te, Ti (rows 1-3) and p (row 4) on s_out for latent vector th.
        real(dp), intent(in) :: th(:)
        real(dp), intent(out) :: f(4, ns_out)
        real(dp) :: g(N_PROFILES, ns_out), dg(N_PROFILES, ns_out), z(model%n_prof), pp(ns_out)
        z = th(model%n_rec + 1:model%n_rec + model%n_prof)
        call model%prof%evaluate(z, s_out, g, dg)
        f(1:3, :) = g([P_NE, P_TE, P_TI], :)
        call model%prof%pressure(z, s_out, pp)
        f(4, :) = pp
    end subroutine profiles_of

    subroutine write_profile_truth(file)
        character(len=*), intent(in) :: file
        real(dp) :: f(4, ns_out)
        integer :: u, j
        call profiles_of(theta_true, f)
        open(newunit=u, file=file, status='replace')
        write(u, '(A)') 's,ne,Te,Ti,p'
        do j = 1, ns_out
            write(u, '(*(ES15.7,:,","))') s_out(j), f(:, j)
        end do
        close(u)
    end subroutine write_profile_truth

    subroutine write_surfaces(file, tag)
        !! R, Z of surfaces rho = 0.25 ... 1 at three toroidal cuts of the
        !! current equilibrium (tag identifies the draw; 0 = truth, -1 = MAP).
        character(len=*), intent(in) :: file
        integer, intent(in) :: tag
        real(dp) :: rr, zz, d(2, 3), phi
        integer :: u, ic, ir, it
        logical :: exists
        inquire(file=file, exist=exists)
        if (exists .and. tag > 0) then
            open(newunit=u, file=file, status='old', position='append')
        else
            open(newunit=u, file=file, status='replace')
            write(u, '(A)') 'draw,cut,rho,theta,R,Z'
        end if
        do ic = 0, 2
            phi = ic * acos(-1.0_dp) / (3.0_dp * rec%eq%nfp)
            do ir = 1, 4
                do it = 0, 64
                    call model%state%point(0.25_dp * ir, it * 2.0_dp * acos(-1.0_dp) / 64, phi, rr, zz, d)
                    write(u, '(I0,",",I0,*(",",ES15.7))') tag, ic, 0.25_dp * ir, &
                        it * 2.0_dp * acos(-1.0_dp) / 64, rr, zz
                end do
            end do
        end do
        close(u)
    end subroutine write_surfaces

    subroutine write_draws(file)
        character(len=*), intent(in) :: file
        integer :: u, j
        open(newunit=u, file=file, status='replace')
        write(u, '(*(A,:,","))') (trim(names(j)), j = 1, model%n), 'log_density', 'accept', 'depth', &
            'divergent'
        do j = 1, size(chain%x, 2)
            write(u, '(*(G0,:,","))') chain%x(:, j), chain%log_density(j), &
                chain%accept_stat(j), chain%depth(j), merge(1, 0, chain%divergent(j))
        end do
        close(u)
    end subroutine write_draws

    subroutine write_parameters(file)
        !! Per latent: truth, MAP, Laplace sd, PSIS mean/sd, NUTS mean/sd and
        !! central 68 % and 95 % intervals.
        character(len=*), intent(in) :: file
        real(dp) :: q(4), m, sd, pm, psd
        integer :: u, j
        open(newunit=u, file=file, status='replace')
        write(u, '(A)') 'name,truth,map,laplace_sd,psis_mean,psis_sd,nuts_mean,nuts_sd,'// &
            'nuts_q025,nuts_q16,nuts_q84,nuts_q975'
        do j = 1, model%n
            q = 0.0_dp
            m = 0.0_dp
            sd = 0.0_dp
            if (n_samples > 0) then
                m = sum(chain%x(j, :)) / n_samples
                sd = sqrt(sum((chain%x(j, :) - m)**2) / max(n_samples - 1, 1))
                q = quantiles(chain%x(j, :), [0.025_dp, 0.16_dp, 0.84_dp, 0.975_dp])
            end if
            pm = 0.0_dp
            psd = 0.0_dp
            if (allocated(ps%mean)) then
                pm = ps%mean(j)
                psd = ps%sd(j)
            end if
            write(u, '(A,*(",",ES15.7))') trim(names(j)), theta_true(j), lap%mean(j), &
                sqrt(lap%cov(j, j)), pm, psd, m, sd, q
        end do
        close(u)
    end subroutine write_parameters

    subroutine write_profiles(file)
        !! Posterior bands of ne, Te, Ti, p on s_out from Laplace draws (no
        !! solve needed) and from the NUTS draws: quantiles 2.5, 16, 50, 84, 97.5 %.
        character(len=*), intent(in) :: file
        real(dp), allocatable :: lap_f(:, :, :), nuts_f(:, :, :)
        real(dp) :: z(model%n), q(5)
        type(rng_t) :: g
        integer :: u, j, k2, c
        character(len=2), parameter :: qn(4) = ['ne', 'Te', 'Ti', 'p ']
        call rng_seed(g, int(seed + 31, int64), st)
        allocate(lap_f(4, ns_out, laplace_draws))
        do k2 = 1, laplace_draws
            do j = 1, model%n
                call rng_normal(g, z(j))
            end do
            call profiles_of(lap%mean + matmul(lap%chol, z), lap_f(:, :, k2))
        end do
        if (n_samples > 0) then
            allocate(nuts_f(4, ns_out, n_samples))
            do k2 = 1, n_samples
                call profiles_of(chain%x(:, k2), nuts_f(:, :, k2))
            end do
        end if
        open(newunit=u, file=file, status='replace')
        write(u, '(A)') 'method,quantity,s,q025,q16,q50,q84,q975'
        do c = 1, 4
            do j = 1, ns_out
                q = quantiles(lap_f(c, j, :), [0.025_dp, 0.16_dp, 0.5_dp, 0.84_dp, 0.975_dp])
                write(u, '(A,",",A,*(",",ES15.7))') 'laplace', trim(qn(c)), s_out(j), q
                if (n_samples > 0) then
                    q = quantiles(nuts_f(c, j, :), [0.025_dp, 0.16_dp, 0.5_dp, 0.84_dp, 0.975_dp])
                    write(u, '(A,",",A,*(",",ES15.7))') 'nuts', trim(qn(c)), s_out(j), q
                end if
            end do
        end do
        close(u)
    end subroutine write_profiles

    subroutine write_posterior_solves()
        !! iota profiles and flux surfaces for the truth, the MAP and up to
        !! posterior_solves draws (NUTS draws if available, else Laplace).
        real(dp), allocatable :: draw(:)
        real(dp) :: z(model%n)
        type(rng_t) :: g
        integer :: u, k2, nd, every
        open(newunit=u, file=out//'/iota.csv', status='replace')
        write(u, '(A)') 'draw,s,iota'
        call model%residual(theta_true, r, ok)
        call write_iota(u, 0)
        call model%residual(lap%mean, r, ok)
        call write_iota(u, -1)
        call write_surfaces(out//'/surfaces_map.csv', -1)
        call rng_seed(g, int(seed + 57, int64), st)
        nd = posterior_solves
        every = 1
        if (n_samples > 0) every = max(1, n_samples / max(nd, 1))
        do k2 = 1, nd
            if (n_samples > 0) then
                if (k2 * every > n_samples) exit
                draw = chain%x(:, k2 * every)
            else
                do i = 1, model%n
                    call rng_normal(g, z(i))
                end do
                draw = lap%mean + matmul(lap%chol, z)
            end if
            call model%residual(draw, r, ok)
            if (.not. ok) cycle
            call write_iota(u, k2)
            call write_surfaces(out//'/surfaces_draws.csv', k2)
        end do
        close(u)
    end subroutine write_posterior_solves

    subroutine write_iota(u, tag)
        !! Rotational transform of the current equilibrium (wout iotas, half grid).
        integer, intent(in) :: u, tag
        real(dp), allocatable :: sh(:), iotas(:)
        integer :: j
        sh = rec%eq%s_half()
        iotas = rec%eq%vmec%wout('iotas', rec%eq%ns)
        do j = 1, size(sh)
            write(u, '(I0,2(",",ES15.7))') tag, sh(j), iotas(j + 1)
        end do
    end subroutine write_iota

    subroutine write_summary(file)
        character(len=*), intent(in) :: file
        integer :: u
        open(newunit=u, file=file, status='replace')
        write(u, '(A,I0)') 'seed = ', seed
        write(u, '(A,I0)') 'latents = ', model%n
        write(u, '(A,I0)') 'residual_rows = ', model%n_residuals()
        write(u, '(A,ES15.7)') 'chi2_start = ', chi2_start
        write(u, '(A,ES15.7)') 'chi2_map = ', lap%chi2
        write(u, '(A,I0)') 'map_steps = ', steps
        write(u, '(A,A)') 'map_message = ', message
        write(u, '(A,ES15.7)') 'log_evidence = ', lap%log_evidence
        write(u, '(A,L1)') 'evidence_proper = ', lap%proper
        write(u, '(A,ES15.7)') 'psis_khat = ', ps%khat
        write(u, '(A,I0)') 'psis_failures = ', ps%failures
        write(u, '(A,ES15.7)') 'psis_ess = ', ps%ess
        if (n_samples > 0) then
            write(u, '(A,ES15.7)') 'nuts_step_size = ', chain%step_size
            write(u, '(A,I0)') 'nuts_divergent = ', chain%n_divergent
            write(u, '(A,I0)') 'nuts_evaluations = ', chain%n_evaluations
            write(u, '(A,ES15.7)') 'nuts_mean_accept = ', sum(chain%accept_stat) / n_samples
            write(u, '(A,ES15.7)') 'nuts_mean_depth = ', real(sum(chain%depth), dp) / n_samples
            write(u, '(A,ES15.7)') 'nuts_ebfmi = ', ebfmi(chain%energy)
        end if
        write(u, '(A,I0)') 'equilibrium_solves = ', rec%eq%solves
        write(u, '(A,ES15.7)') 'time_map_s = ', t_map
        write(u, '(A,ES15.7)') 'time_laplace_s = ', t_lap
        write(u, '(A,ES15.7)') 'time_psis_s = ', t_psis
        write(u, '(A,ES15.7)') 'time_nuts_s = ', t_nuts
        close(u)
    end subroutine write_summary

    real(dp) function ebfmi(e)
        !! Energy Bayesian fraction of missing information (Betancourt 2016).
        real(dp), intent(in) :: e(:)
        integer :: j
        ebfmi = sum([((e(j) - e(j - 1))**2, j = 2, size(e))]) / &
            max(sum((e - sum(e) / size(e))**2), tiny(1.0_dp))
    end function ebfmi

    function quantiles(v, p) result(q)
        real(dp), intent(in) :: v(:), p(:)
        real(dp) :: q(size(p)), w(size(v)), t, h
        integer :: a, b, j, lo
        w = v
        do a = 2, size(w)
            t = w(a)
            b = a - 1
            do while (b >= 1)
                if (w(b) <= t) exit
                w(b + 1) = w(b)
                b = b - 1
            end do
            w(b + 1) = t
        end do
        do j = 1, size(p)
            h = (size(w) - 1) * p(j) + 1
            lo = min(max(int(h), 1), size(w) - 1)
            q(j) = w(lo) + (h - lo) * (w(lo + 1) - w(lo))
        end do
    end function quantiles

end program tiago_infer
