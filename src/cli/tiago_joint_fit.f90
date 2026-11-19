program tiago_joint_fit
    !! Joint fit of a fixed-boundary VMEC++ equilibrium and kinetic profiles to
    !! measured point values (ne, Te, Ti) and interferometer line integrals
    !! read from files: MAP, Laplace covariance, profile bands, posterior draws
    !! with re-solved equilibria. Device geometry and data stay in the input
    !! files (docs/joint-fit.md).
    !! Usage: tiago_joint_fit input.nml
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use tiago_reconstruction, only: rec
    use tiago_joint_model, only: joint_model_t
    use tiago_kinetic_profiles, only: P_NE, P_TE, P_TI, P_ZEFF, N_PROFILES
    use tiago_profile_diagnostic, only: Q_NE, Q_TE, Q_TI, plasma_sample_t
    use tiago_real_space_state, only: located_t
    use tiago_point_probe, only: point_probe_t
    use tiago_line_integrals, only: interferometer_t
    use tiago_sight_geometry, only: chord_t
    use tiago_inference_engines, only: map_fit, laplace, laplace_t
    use fortnum_rng, only: rng_t, rng_seed, rng_normal
    use fortnum_status, only: fortnum_status_t
    implicit none

    integer, parameter :: maxk = 64, maxp = 32, maxg = 16, ns_out = 101
    character(len=512) :: vmec_input = '', output_dir = 'joint_fit', probes_file = '', &
        chords_file = '', prior_table = ''
    character(len=64) :: names(maxp) = ''
    real(dp) :: prior(maxp) = 0.0_dp, prior_sigma(maxp) = 0.0_dp
    real(dp) :: breakpoints(maxk) = -1.0_dp, z_imp = 5.0_dp, vmec_ftol = 1.0e-12_dp
    real(dp) :: sigma_f(3) = 1.0_dp, l_core = 0.3_dp, l_edge = 0.05_dp, s0 = 0.85_dp, sw = 0.05_dp
    real(dp) :: align_sigma(maxg) = 0.0_dp, cal_sigma(maxg) = 0.0_dp, model_rel(maxg) = 0.0_dp
    real(dp) :: student_nu(maxg) = 0.0_dp
    real(dp) :: chord_cal_sigma = 0.0_dp, chord_model_rel = 0.0_dp
    integer :: max_iterations = 40, laplace_draws = 400, posterior_solves = 0, seed = 1
    logical :: predict_only = .false.
    namelist /joint_fit/ vmec_input, output_dir, probes_file, chords_file, prior_table, names, &
        prior, prior_sigma, breakpoints, z_imp, vmec_ftol, sigma_f, l_core, l_edge, s0, sw, &
        align_sigma, cal_sigma, model_rel, student_nu, chord_cal_sigma, chord_model_rel, &
        max_iterations, laplace_draws, posterior_solves, seed, predict_only

    type(joint_model_t), target :: model
    type(laplace_t) :: lap
    type(fortnum_status_t) :: st
    character(len=512) :: path
    character(len=:), allocatable :: out, message
    character(len=48), allocatable :: pnames(:)
    real(dp), allocatable :: theta(:), r(:), s_out(:)
    integer, allocatable :: group(:), quantity(:)
    real(dp), allocatable :: xp(:, :), val(:), sig(:)
    integer :: unit, ios, i, steps, n_groups, n_chords
    logical :: ok

    call get_command_argument(1, path)
    open(newunit=unit, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'cannot open the namelist file'
    read(unit, nml=joint_fit, iostat=ios)
    if (ios /= 0) error stop 'error in namelist &joint_fit'
    close(unit)
    out = trim(output_dir)
    call execute_command_line('mkdir -p "'//out//'/draws"')

    call setup_equilibrium()
    if (count(breakpoints >= 0.0_dp) < 2) error stop 'breakpoints: at least two knots'
    call model%prof%init(pack(breakpoints, breakpoints >= 0.0_dp), z_imp)
    call set_profile_priors()
    call read_probes()
    n_chords = count_rows(chords_file)
    allocate(model%diags(n_groups + merge(1, 0, n_chords > 0)))
    do i = 1, n_groups
        call fill_probe_slot(i)
    end do
    if (n_chords > 0) call fill_chord_slot(n_groups + 1)
    call model%setup()
    allocate(r(model%n_residuals()))
    theta = [rec%prior / rec%scale, [(0.0_dp, i = 1, model%n - model%n_rec)]]
    call model%residual(theta, r, ok)
    if (.not. ok) error stop 'no equilibrium at the prior mean'
    pnames = parameter_names()
    if (predict_only) then
        call write_predictions(out//'/predictions.csv', theta)
        print '(A)', 'prior-mean predictions written to '//out
        stop
    end if

    call map_fit(model, theta, max_iterations, steps, message)
    call laplace(model, theta, lap)
    if (lap%info /= 0) error stop 'Laplace Hessian not positive definite'
    call model%residual(lap%mean, r, ok)
    call rec%eq%write_wout(out//'/wout_map.nc')
    s_out = [((i - 1.0_dp) / (ns_out - 1), i = 1, ns_out)]
    call write_parameters(out//'/parameters.csv')
    call write_covariance(out//'/covariance.csv')
    call write_predictions(out//'/predictions.csv', lap%mean)
    call write_draws()
    call write_summary(out//'/summary.txt')
    print '(A)', 'joint fit written to '//out

contains

    subroutine setup_equilibrium()
        !! Equilibrium parameters with Gaussian priors and no magnetic sensors.
        integer :: k, n
        n = count(len_trim(names) > 0)
        if (n < 1) error stop 'names: at least one equilibrium parameter'
        rec%n = n
        allocate(rec%names(n), rec%is_coil(n), rec%group(n), rec%eq_index(n), rec%x0(n), &
            rec%scale(n), rec%prior(n), rec%prior_sigma(n), rec%history(0))
        do k = 1, n
            rec%names(k) = names(k)
            rec%is_coil(k) = .false.
            rec%group(k) = 0
            rec%eq_index(k) = k
        end do
        call execute_command_line('mkdir -p "'//out//'/work"')
        call rec%eq%init(trim(vmec_input), out//'/work', rec%names)
        call rec%eq%vmec%set_input('ftol', vmec_ftol)
        call rec%eq%vmec%set_input('niter', 20000.0_dp)
        rec%rule%samples_per_segment = 4
        rec%rule%gauss = .true.
        rec%nphi = 32
        rec%ntheta = 32
        allocate(rec%loops(0), rec%segs(0), rec%probes(0))
        rec%n_diag_probes = 0
        rec%n_fb = 0
        call rec%coils%init('')
        allocate(rec%extcur(0), rec%flux_m(0, 0), rec%seg_m(0, 0), rec%probe_m(0, 0))
        call rec%plasma%init('')
        rec%n_meas = 0
        allocate(rec%meas_kind(0), rec%meas_index(0), rec%meas_value(0), rec%meas_sigma(0), &
            rec%meas_label(0))
        rec%prior = prior(1:n)
        rec%prior_sigma = prior_sigma(1:n)
        rec%scale = max(abs(rec%prior), 1.0e-3_dp)
        rec%x0 = rec%prior
        rec%prior_param = pack([(k, k = 1, n)], rec%prior_sigma > 0.0_dp)
        rec%n_prior = size(rec%prior_param)
    end subroutine setup_equilibrium

    subroutine set_profile_priors()
        !! Prior means of ne, Te, Ti (latent) and Zeff (fixed) at the Greville
        !! abscissae, interpolated linearly in s from prior_table (s, ne, Te, Ti, Zeff).
        real(dp), allocatable :: t(:, :)
        real(dp) :: m(4, model%prof%ncoef)
        integer :: u, n, j, c
        n = count_rows(prior_table)
        if (n < 2) error stop 'prior_table: at least two rows'
        allocate(t(5, n))
        open(newunit=u, file=trim(prior_table), status='old', action='read')
        call skip_header(u)
        do j = 1, n
            read(u, *) t(:, j)
        end do
        close(u)
        do c = 1, 4
            do j = 1, model%prof%ncoef
                m(c, j) = interp(t(1, :), t(c + 1, :), model%prof%greville(j))
            end do
        end do
        call model%prof%set_profile(P_NE, m(1, :), .true., sigma_f(1), l_core, l_edge, s0, sw)
        call model%prof%set_profile(P_TE, m(2, :), .true., sigma_f(2), l_core, l_edge, s0, sw)
        call model%prof%set_profile(P_TI, m(3, :), .true., sigma_f(3), l_core, l_edge, s0, sw)
        call model%prof%set_profile(P_ZEFF, m(4, :), .false.)
    end subroutine set_profile_priors

    subroutine read_probes()
        !! Rows: group, quantity (1 ne [m^-3], 2 Te [eV], 3 Ti [eV]), R [m], phi [rad], Z [m],
        !! value, sigma. Groups are numbered 1..G, each with its own alignment and calibration.
        integer :: u, n, j
        n = count_rows(probes_file)
        if (n < 1) error stop 'probes_file: no rows'
        allocate(group(n), quantity(n), xp(3, n), val(n), sig(n))
        open(newunit=u, file=trim(probes_file), status='old', action='read')
        call skip_header(u)
        do j = 1, n
            read(u, *) group(j), quantity(j), xp(:, j), val(j), sig(j)
        end do
        close(u)
        if (any(quantity < 1 .or. quantity > 3)) error stop 'probes_file: quantity must be 1, 2 or 3'
        if (any(sig <= 0.0_dp)) error stop 'probes_file: sigma must be positive'
        n_groups = maxval(group)
        if (minval(group) < 1 .or. n_groups > maxg) error stop 'probes_file: groups 1..16'
    end subroutine read_probes

    subroutine fill_probe_slot(g)
        integer, intent(in) :: g
        type(point_probe_t) :: p
        logical, allocatable :: mask(:)
        integer, parameter :: qmap(3) = [Q_NE, Q_TE, Q_TI]
        integer :: m
        mask = group == g
        m = count(mask)
        if (m == 0) error stop 'probes_file: empty group'
        p%x0 = xp(:, pack([(i, i = 1, size(group))], mask))
        p%quantity = qmap(pack(quantity, mask))
        p%aligned = align_sigma(g) <= 0.0_dp
        associate (slot => model%diags(g))
            allocate(slot%d, source=p)
            if (align_sigma(g) > 0.0_dp) then
                slot%eta_sigma = [align_sigma(g), align_sigma(g)]
            else
                allocate(slot%eta_sigma(0))
            end if
            slot%cal_sigma = cal_sigma(g)
            slot%student_nu = student_nu(g)
            slot%data = pack(val, mask)
            slot%covariance = diag(pack(sig, mask)**2)
            slot%model_sigma = model_rel(g) * abs(slot%data)
        end associate
    end subroutine fill_probe_slot

    subroutine fill_chord_slot(k)
        !! Rows: x0, y0, z0, dx, dy, dz, length [m], value [m^-2], sigma [m^-2].
        integer, intent(in) :: k
        type(interferometer_t) :: interf
        type(chord_t), allocatable :: chords(:)
        real(dp) :: row(9)
        real(dp), allocatable :: v(:), s(:)
        integer :: u, j
        character(len=32) :: tag
        allocate(chords(n_chords), v(n_chords), s(n_chords))
        open(newunit=u, file=trim(chords_file), status='old', action='read')
        call skip_header(u)
        do j = 1, n_chords
            read(u, *) row
            write(tag, '("chord_",I0)') j
            chords(j)%label = trim(tag)
            chords(j)%origin = row(1:3)
            chords(j)%direction = row(4:6) / norm2(row(4:6))
            chords(j)%length = row(7)
            chords(j)%npanel = 80
            chords(j)%ngauss = 4
            v(j) = row(8)
            s(j) = row(9)
        end do
        close(u)
        call interf%init(chords)
        associate (slot => model%diags(k))
            allocate(slot%d, source=interf)
            allocate(slot%eta_sigma(0))
            slot%cal_sigma = chord_cal_sigma
            slot%data = v
            slot%covariance = diag(s**2)
            slot%model_sigma = chord_model_rel * abs(v)
        end associate
    end subroutine fill_chord_slot

    pure function diag(v) result(a)
        real(dp), intent(in) :: v(:)
        real(dp) :: a(size(v), size(v))
        integer :: j
        a = 0.0_dp
        do j = 1, size(v)
            a(j, j) = v(j)
        end do
    end function diag

    integer function count_rows(file)
        !! Data rows of a file whose first line is a header; zero if no file is named.
        character(len=*), intent(in) :: file
        character(len=1024) :: line
        integer :: u, io
        count_rows = 0
        if (len_trim(file) == 0) return
        open(newunit=u, file=trim(file), status='old', action='read', iostat=io)
        if (io /= 0) error stop 'cannot open '//trim(file)
        read(u, '(A)', iostat=io) line
        do
            read(u, '(A)', iostat=io) line
            if (io /= 0) exit
            if (len_trim(line) > 0) count_rows = count_rows + 1
        end do
        close(u)
    end function count_rows

    subroutine skip_header(u)
        integer, intent(in) :: u
        character(len=1) :: c
        read(u, '(A)') c
    end subroutine skip_header

    pure real(dp) function interp(x, y, x0)
        real(dp), intent(in) :: x(:), y(:), x0
        integer :: j
        if (x0 <= x(1)) then
            interp = y(1)
            return
        end if
        do j = 2, size(x)
            if (x0 <= x(j)) then
                interp = y(j - 1) + (y(j) - y(j - 1)) * (x0 - x(j - 1)) / (x(j) - x(j - 1))
                return
            end if
        end do
        interp = y(size(y))
    end function interp

    function parameter_names() result(nm)
        character(len=48), allocatable :: nm(:)
        character(len=2), parameter :: tags(3) = ['ne', 'Te', 'Ti']
        integer :: j, q, a
        allocate(nm(model%n))
        do j = 1, model%n
            write(nm(j), '("latent_",I0)') j
        end do
        do j = 1, rec%n
            nm(j) = rec%names(j)
        end do
        do q = 1, 3
            do j = 1, model%prof%ncoef
                write(nm(rec%n + (q - 1) * model%prof%ncoef + j), '(A,"_z",I0)') tags(q), j
            end do
        end do
        do q = 1, size(model%diags)
            associate (slot => model%diags(q))
                do a = 1, size(slot%eta_sigma)
                    write(nm(slot%eta_offset + a), '("group",I0,"_",A)') q, merge('dR', 'dZ', a == 1)
                end do
                if (slot%cal_index > 0) write(nm(slot%cal_index), '("group",I0,"_logC")') q
            end associate
        end do
    end function parameter_names

    subroutine profiles_of(th, f)
        !! ne, Te, Ti and pressure on s_out for latent vector th.
        real(dp), intent(in) :: th(:)
        real(dp), intent(out) :: f(4, ns_out)
        real(dp) :: g(N_PROFILES, ns_out), dg(N_PROFILES, ns_out), pp(ns_out)
        real(dp), allocatable :: z(:)
        z = th(model%n_rec + 1:model%n_rec + model%n_prof)
        call model%prof%evaluate(z, s_out, g, dg)
        f(1:3, :) = g([P_NE, P_TE, P_TI], :)
        call model%prof%pressure(z, s_out, pp)
        f(4, :) = pp
    end subroutine profiles_of

    subroutine write_parameters(file)
        character(len=*), intent(in) :: file
        integer :: u, j, d, a
        real(dp) :: f(model%n), c(model%n)
        ! physical value = c + f * latent: equilibrium parameters in their units,
        ! alignment shifts in m, calibrations as log C, profile latents standard normal
        f = 1.0_dp
        c = 0.0_dp
        f(1:rec%n) = rec%scale
        do d = 1, size(model%diags)
            associate (slot => model%diags(d))
                do a = 1, size(slot%eta_sigma)
                    f(slot%eta_offset + a) = slot%eta_sigma(a)
                end do
                if (slot%cal_index > 0) then
                    f(slot%cal_index) = slot%cal_sigma
                    c(slot%cal_index) = slot%cal_mean
                end if
            end associate
        end do
        open(newunit=u, file=file, status='replace')
        write(u, '(A)') 'name,map_latent,sd_latent,map_physical,sd_physical'
        do j = 1, model%n
            write(u, '(A,4(",",ES16.8))') trim(pnames(j)), lap%mean(j), sqrt(lap%cov(j, j)), &
                c(j) + f(j) * lap%mean(j), f(j) * sqrt(lap%cov(j, j))
        end do
        close(u)
    end subroutine write_parameters

    subroutine write_covariance(file)
        character(len=*), intent(in) :: file
        integer :: u, j
        open(newunit=u, file=file, status='replace')
        write(u, '(*(A,:,","))') (trim(pnames(j)), j = 1, model%n)
        do j = 1, model%n
            write(u, '(*(ES16.8,:,","))') lap%cov(j, :)
        end do
        close(u)
    end subroutine write_covariance

    subroutine write_predictions(file, th)
        !! Data, prediction at th and data sigma per signal.
        character(len=*), intent(in) :: file
        real(dp), intent(in) :: th(:)
        type(plasma_sample_t) :: sample
        type(located_t), allocatable :: loc(:)
        real(dp), allocatable :: y(:), eta(:), dxe(:, :, :), scale(:)
        integer :: u, d, j
        open(newunit=u, file=file, status='replace')
        write(u, '(A)') 'group,signal,label,data,model,sigma'
        do d = 1, size(model%diags)
            associate (slot => model%diags(d))
                allocate(y(size(slot%data)))
                call model%predict_diagnostic(th, d, y, sample, loc, eta, dxe)
                scale = model%signal_scale(th, d)
                do j = 1, size(y)
                    write(u, '(I0,",",I0,",",A,3(",",ES16.8))') d, j, '"'//trim(slot%d%label(j))//'"', &
                        slot%data(j), y(j) * scale(j), &
                        sqrt(slot%covariance(j, j) + slot%model_sigma(j)**2)
                end do
                deallocate(y)
            end associate
        end do
        close(u)
    end subroutine write_predictions

    subroutine write_draws()
        !! Profile quantiles from laplace_draws Laplace draws, and for the first
        !! posterior_solves draws the latent vector, profiles and re-solved wout.
        real(dp), allocatable :: f(:, :, :), z(:), th(:)
        real(dp) :: q(5)
        type(rng_t) :: g
        integer :: u, v, k, j, c
        character(len=2), parameter :: qn(4) = ['ne', 'Te', 'Ti', 'p ']
        character(len=16) :: tag
        call rng_seed(g, int(seed, int64), st)
        allocate(f(4, ns_out, laplace_draws), z(model%n))
        open(newunit=v, file=out//'/draws/latents.csv', status='replace')
        write(v, '(*(A,:,","))') 'draw', 'solved', (trim(pnames(j)), j = 1, model%n)
        do k = 1, laplace_draws
            do j = 1, model%n
                call rng_normal(g, z(j))
            end do
            th = lap%mean + matmul(lap%chol, z)
            call profiles_of(th, f(:, :, k))
            ok = .false.
            if (k <= posterior_solves) then
                call model%residual(th, r, ok)
                if (ok) then
                    write(tag, '(I4.4)') k
                    call rec%eq%write_wout(out//'/draws/wout_'//trim(tag)//'.nc')
                end if
            end if
            write(v, '(I0,",",I0,*(:,",",ES16.8))') k, merge(1, 0, ok), th
        end do
        close(v)
        open(newunit=u, file=out//'/profiles.csv', status='replace')
        write(u, '(A)') 'quantity,s,map,q025,q16,q50,q84,q975'
        call profiles_of(lap%mean, f(:, :, 1))
        do c = 1, 4
            do j = 1, ns_out
                q = quantiles(f(c, j, 2:), [0.025_dp, 0.16_dp, 0.5_dp, 0.84_dp, 0.975_dp])
                write(u, '(A,*(:,",",ES16.8))') trim(qn(c)), s_out(j), f(c, j, 1), q
            end do
        end do
        close(u)
        open(newunit=u, file=out//'/draws/profiles.csv', status='replace')
        write(u, '(A)') 'draw,s,ne,Te,Ti,p'
        do k = 2, laplace_draws
            do j = 1, ns_out
                write(u, '(I0,*(:,",",ES16.8))') k, s_out(j), f(:, j, k)
            end do
        end do
        close(u)
    end subroutine write_draws

    subroutine write_summary(file)
        character(len=*), intent(in) :: file
        integer :: u
        open(newunit=u, file=file, status='replace')
        write(u, '(A,I0)') 'latents = ', model%n
        write(u, '(A,I0)') 'residual_rows = ', model%n_residuals()
        write(u, '(A,ES15.7)') 'chi2_map = ', lap%chi2
        write(u, '(A,I0)') 'map_steps = ', steps
        write(u, '(A,A)') 'map_message = ', message
        write(u, '(A,ES15.7)') 'log_evidence = ', lap%log_evidence
        write(u, '(A,I0)') 'equilibrium_solves = ', rec%eq%solves
        close(u)
    end subroutine write_summary

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

end program tiago_joint_fit
