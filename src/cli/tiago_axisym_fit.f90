program tiago_axisym_fit
    !! Fit an axisymmetric free-boundary equilibrium to magnetic sensors, coil
    !! currents and profile rows through an equilibrium provider process
    !! (docs/axisymmetric-free-boundary.md). Usage: tiago_axisym_fit input.nml
    !!
    !! &axisym_fit
    !!   machine, provider, workdir      machine file, provider command, scratch directory
    !!   measurements                     CSV label,value,sigma: every sensor label of
    !!                                    the machine, and rows pressure@<psi_N> [Pa]
    !!                                    jphi@<psi_N> [A/m^2] or jdotb@<psi_N> [A T/m^2]
    !!                                    (docs/bootstrap.md), and every MSE channel
    !!                                    label [rad]
    !!   mse_channels                     MSE channel file (docs/mse.md), optional
    !!   nx, ny, rmin, rmax, zmin, zmax, tolerance, fvac, alpha_m, alpha_n
    !!   coil_current(g), plasma_current, pressure_axis
    !!   profile ('paxis' or 'knots'), knot_psin(k), pprime_knot(k), ffprime_knot(k)
    !!   parameters(i), start_values(i), parameter_scale(i), prior_values(i), prior_sigma(i)
    !!   max_iterations, fd_relative, draws, seed
    !!   edge_psin(j)                     labels of the reported outboard j_phi
    !!   output                           prefix of the result files
    !!   predict_geqdsk                   only predict the sensors of this G-EQDSK
    !!   initial                          plasma flux (.npy) every provider solve starts from
    !!   probes_file, chords_file         profile diagnostics in the format of docs/joint-fit.md;
    !!                                    either one makes the fit joint with the kinetic profiles
    !!   prior_table, breakpoints(k), z_imp, sigma_f(3), l_core, l_edge, s0, sw
    !!   align_sigma(g), cal_sigma(g), model_rel(g), chord_cal_sigma, chord_model_rel
    !!   pressure_psin(j), pressure_sigma(j)  labels and sd [Pa] of the rows p_eq - p_kin
    !! /
    !!
    !! Writes <output>_parameters.csv (MAP, local standard deviation, start),
    !! <output>_signals.csv (measured, predicted, sigma), <output>_map.geqdsk and,
    !! with draws > 0, <output>_draw_NNNN.geqdsk of Laplace draws re-solved by the
    !! provider and <output>_edge.csv with j_phi at edge_psin for the MAP and
    !! every draw. A joint fit also writes the profile latents LATENT(i) into
    !! <output>_parameters.csv, <output>_profile_signals.csv (slot, measured,
    !! predicted with calibration, sigma) and <output>_profiles.csv (n_e, T_e,
    !! T_i, p and their Laplace sd on 101 values of s); its predict mode writes
    !! the profile signals at the prior mean.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, error_unit
    use fortnum_status, only: fortnum_status_t
    use fortnum_rng, only: rng_t, rng_seed, rng_normal
    use tiago_axisym_free_boundary, only: axisym_free_boundary_t
    use tiago_axisym_machine, only: read_axisym_machine
    use tiago_geqdsk_plasma, only: geqdsk_plasma_t, read_geqdsk_plasma
    use tiago_inference_engines, only: map_fit, laplace, laplace_t
    use tiago_mse, only: read_mse_file
    use tiago_axisym_joint, only: axisym_joint_t
    use tiago_profile_inputs, only: set_profile_priors, read_probe_slots, read_chord_slot
    use tiago_kinetic_profiles, only: P_NE, P_TE, P_TI, N_PROFILES
    use tiago_profile_diagnostic, only: plasma_sample_t
    implicit none

    integer, parameter :: mp = 64, mk = 32, mg = 16
    character(len=1024) :: machine, provider, workdir, measurements, output, nml, predict_geqdsk, initial, &
        mse_channels, probes_file, chords_file, prior_table
    character(len=16) :: profile
    character(len=32) :: parameters(mp)
    integer :: nx, ny, max_iterations, draws, seed, nknots
    real(dp) :: rmin, rmax, zmin, zmax, tolerance, fvac, alpha_m, alpha_n
    real(dp) :: coil_current(mp), plasma_current, pressure_axis, fd_relative
    real(dp) :: knot_psin(mk), pprime_knot(mk), ffprime_knot(mk), edge_psin(16)
    real(dp) :: start_values(mp), parameter_scale(mp), prior_values(mp), prior_sigma(mp)
    real(dp) :: breakpoints(mk), z_imp, sigma_f(3), l_core, l_edge, s0, sw, align_sigma(mg), &
        cal_sigma(mg), model_rel(mg), chord_cal_sigma, chord_model_rel, pressure_psin(mk), &
        pressure_sigma(mk)
    namelist /axisym_fit/ machine, provider, workdir, measurements, output, nx, ny, &
        rmin, rmax, zmin, zmax, tolerance, fvac, alpha_m, alpha_n, coil_current, &
        plasma_current, pressure_axis, profile, knot_psin, pprime_knot, ffprime_knot, &
        parameters, start_values, parameter_scale, prior_values, prior_sigma, &
        max_iterations, fd_relative, draws, seed, edge_psin, predict_geqdsk, initial, mse_channels, &
        probes_file, chords_file, prior_table, breakpoints, z_imp, sigma_f, l_core, l_edge, s0, sw, &
        align_sigma, cal_sigma, model_rel, chord_cal_sigma, chord_model_rel, pressure_psin, &
        pressure_sigma

    type(axisym_free_boundary_t), target :: model
    type(axisym_joint_t), target :: joint
    type(axisym_free_boundary_t), pointer :: fb
    type(laplace_t) :: lap
    type(geqdsk_plasma_t) :: eq
    type(rng_t) :: rng
    type(fortnum_status_t) :: status
    character(len=:), allocatable :: message
    character(len=32) :: suffix
    real(dp), allocatable :: theta(:), start(:), s(:), z(:), td(:)
    integer :: unit, ierr, n, ng, ne, k, d, steps
    logical :: ok, conv, joint_mode

    machine = ''
    provider = ''
    workdir = 'axisym_work'
    measurements = ''
    output = 'axisym_fit'
    nx = 33
    ny = 33
    rmin = 0.1_dp
    rmax = 2.0_dp
    zmin = -1.0_dp
    zmax = 1.0_dp
    tolerance = 1.0e-6_dp
    fvac = 1.0_dp
    alpha_m = 1.0_dp
    alpha_n = 2.0_dp
    coil_current = 0
    plasma_current = 0
    pressure_axis = 0
    profile = 'paxis'
    knot_psin = -1
    pprime_knot = 0
    ffprime_knot = 0
    parameters = ''
    start_values = huge(1.0_dp)
    parameter_scale = 1
    prior_values = 0
    prior_sigma = 0
    max_iterations = 20
    fd_relative = 1.0e-3_dp
    draws = 0
    seed = 1
    edge_psin = -1
    predict_geqdsk = ''
    initial = ''
    mse_channels = ''
    probes_file = ''
    chords_file = ''
    prior_table = ''
    breakpoints = -1
    z_imp = 5.0_dp
    sigma_f = 1.0_dp
    l_core = 0.3_dp
    l_edge = 0.05_dp
    s0 = 0.85_dp
    sw = 0.05_dp
    align_sigma = 0
    cal_sigma = 0
    model_rel = 0
    chord_cal_sigma = 0
    chord_model_rel = 0
    pressure_psin = -1
    pressure_sigma = 0
    call get_command_argument(1, nml)
    if (len_trim(nml) == 0) error stop 'usage: tiago_axisym_fit input.nml'
    open(newunit=unit, file=trim(nml), status='old', action='read')
    read(unit, nml=axisym_fit)
    close(unit)

    model%machine_path = trim(machine)
    call read_axisym_machine(model%machine_path, model%machine, ierr, message)
    if (ierr /= 0) error stop message
    model%command = trim(provider)
    model%workdir = trim(workdir)
    call execute_command_line('mkdir -p '//model%workdir)
    model%nx = nx
    model%ny = ny
    model%rmin = rmin
    model%rmax = rmax
    model%zmin = zmin
    model%zmax = zmax
    model%tolerance = tolerance
    model%fvac = fvac
    model%alpha_m = alpha_m
    model%alpha_n = alpha_n
    ng = model%machine%ngroups
    model%coil_current = coil_current(:ng)
    model%plasma_current = plasma_current
    model%pressure_axis = pressure_axis
    model%profile = trim(profile)
    nknots = count(knot_psin >= 0)
    model%knot_psin = knot_psin(:nknots)
    model%pprime_knot = pprime_knot(:nknots)
    model%ffprime_knot = ffprime_knot(:nknots)
    model%fd_relative = fd_relative
    model%initial = trim(initial)
    n = count(len_trim(parameters) > 0)
    model%names = parameters(:n)
    model%scale = parameter_scale(:n)
    model%prior_mean = prior_values(:n)
    model%prior_sigma = prior_sigma(:n)
    if (len_trim(mse_channels) > 0) then
        call read_mse_file(trim(mse_channels), model%mse, ierr, message)
        if (ierr /= 0) error stop message
    end if
    call read_measurements(trim(measurements))
    joint_mode = len_trim(probes_file) > 0 .or. len_trim(chords_file) > 0
    fb => model
    if (joint_mode) then
        call setup_joint()
        fb => joint%fb
    end if
    if (len_trim(predict_geqdsk) > 0) then
        call predict(trim(predict_geqdsk))
        stop
    end if

    allocate(theta(n))
    do k = 1, n
        if (start_values(k) == huge(1.0_dp)) then
            theta(k) = default_value(model%names(k)) / model%scale(k)
        else
            theta(k) = start_values(k) / model%scale(k)
        end if
    end do
    if (joint_mode) then
        ! profile latents start at their prior mean
        theta = [theta, [(0.0_dp, k = 1, joint%pp%n)]]
        start = theta
        call map_fit(joint, theta, max_iterations, steps, message, conv)
    else
        start = theta
        call map_fit(model, theta, max_iterations, steps, message, conv)
    end if
    write(*, '(A,I0,A,I0,A,L1)') 'fit: ', steps, ' steps, ', fb%solves, &
        ' provider solves, converged ', conv
    if (.not. conv) write(error_unit, '(A)') 'fit: '//message
    if (joint_mode) then
        call laplace(joint, theta, lap)
    else
        call laplace(model, theta, lap)
    end if
    if (lap%info /= 0) error stop 'Laplace covariance failed'

    open(newunit=unit, file=trim(output)//'_parameters.csv', status='replace')
    write(unit, '(A)') 'name,value,sd,start'
    do k = 1, n
        write(unit, '(A,3(",",ES24.16))') trim(fb%names(k)), theta(k) * fb%scale(k), &
            sqrt(lap%cov(k, k)) * fb%scale(k), start(k) * fb%scale(k)
    end do
    do k = n + 1, size(theta)
        write(suffix, '(A,I0,A)') 'LATENT(', k - n, ')'
        write(unit, '(A,3(",",ES24.16))') trim(suffix), theta(k), sqrt(lap%cov(k, k)), start(k)
    end do
    close(unit)

    call fb%signals(theta(:n), s, ok, eq)
    if (.not. ok) error stop 'no equilibrium at the MAP'
    call execute_command_line('cp '//fb%workdir//'/equilibrium.geqdsk '//trim(output) &
        //'_map.geqdsk')
    if (joint_mode) then
        joint%pp%eq = eq
        call write_profile_signals(theta(n + 1:))
        call write_profiles()
    end if
    open(newunit=unit, file=trim(output)//'_signals.csv', status='replace')
    write(unit, '(A)') 'label,measured,predicted,sigma'
    do k = 1, size(s)
        write(unit, '(A,3(",",ES24.16))') trim(row_label(k)), model%data(k), s(k), model%sigma(k)
    end do
    close(unit)

    ne = count(edge_psin >= 0)
    if (ne > 0) then
        open(newunit=unit, file=trim(output)//'_edge.csv', status='replace')
        write(unit, '(A)') 'sample,psin,r,jphi'
        call write_edge(unit, 0, eq)
        call rng_seed(rng, int(seed, int64), status)
        allocate(z(size(theta)))
        do d = 1, draws
            do k = 1, size(theta)
                call rng_normal(rng, z(k))
            end do
            td = theta + matmul(lap%chol, z)
            call fb%signals(td(:n), s, ok, eq)
            if (.not. ok) then
                write(error_unit, '(A,I0,A)') 'draw ', d, ' did not converge, skipped'
                cycle
            end if
            write(suffix, '(A,I4.4,A)') '_draw_', d, '.geqdsk'
            call execute_command_line('cp '//fb%workdir//'/equilibrium.geqdsk '//trim(output) &
                //trim(suffix))
            call write_edge(unit, d, eq)
        end do
        close(unit)
    end if

contains

    subroutine setup_joint()
        !! Joint model of the configured free-boundary model and the profile
        !! diagnostics, profile prior and pressure rows of the namelist.
        joint%fb = model
        call joint%pp%prof%init(pack(breakpoints, breakpoints >= 0), z_imp)
        if (len_trim(prior_table) == 0) error stop 'a joint fit needs prior_table'
        call set_profile_priors(joint%pp%prof, trim(prior_table), sigma_f, l_core, l_edge, s0, sw)
        if (len_trim(probes_file) > 0) &
            call read_probe_slots(trim(probes_file), align_sigma, cal_sigma, model_rel, joint%pp%slots)
        if (len_trim(chords_file) > 0) &
            call read_chord_slot(trim(chords_file), chord_cal_sigma, chord_model_rel, joint%pp%slots)
        call joint%pp%setup()
        joint%p_psin = pack(pressure_psin, pressure_psin >= 0)
        joint%p_sigma = pressure_sigma(:size(joint%p_psin))
        if (any(joint%p_sigma <= 0)) error stop 'pressure_sigma must be positive'
    end subroutine setup_joint

    subroutine write_profile_signals(zp)
        !! Calibrated predictions of every profile slot at the latents zp in
        !! the equilibrium joint%pp%eq.
        real(dp), intent(in) :: zp(:)
        real(dp), allocatable :: y(:), eta(:), dx(:, :, :), sd(:)
        type(plasma_sample_t) :: sample
        real(dp) :: c
        integer :: u, j, i
        open(newunit=u, file=trim(output)//'_profile_signals.csv', status='replace')
        write(u, '(A)') 'slot,index,measured,predicted,sigma'
        do j = 1, size(joint%pp%slots)
            call joint%pp%predict(zp, j, y, sample, eta, dx, c)
            sd = joint%pp%sigma_of(j)
            do i = 1, size(y)
                write(u, '(I0,",",I0,3(",",ES24.16))') j, i, joint%pp%slots(j)%data(i), c * y(i), sd(i)
            end do
        end do
        close(u)
    end subroutine write_profile_signals

    subroutine write_profiles()
        !! MAP profiles and pressure on 101 values of s with the sd of 2000
        !! Laplace draws of the latents (the profiles do not depend on the
        !! equilibrium parameters).
        integer, parameter :: ns = 101, nd = 2000
        real(dp) :: sg(ns), f(N_PROFILES, ns), df(N_PROFILES, ns), p(ns), m(4, ns), m2(4, ns), &
            fm(4, ns), w(size(theta))
        integer :: u, i, j, np
        np = joint%pp%n_prof
        sg = [((i - 1) / real(ns - 1, dp), i = 1, ns)]
        call joint%pp%prof%evaluate(theta(n + 1:n + np), sg, f, df)
        call joint%pp%prof%pressure(theta(n + 1:n + np), sg, p)
        fm = transpose(reshape([f(P_NE, :), f(P_TE, :), f(P_TI, :), p], [ns, 4]))
        m = 0
        m2 = 0
        call rng_seed(rng, int(seed, int64) + 7, status)
        do j = 1, nd
            do i = 1, size(w)
                call rng_normal(rng, w(i))
            end do
            w = theta + matmul(lap%chol, w)
            call joint%pp%prof%evaluate(w(n + 1:n + np), sg, f, df)
            call joint%pp%prof%pressure(w(n + 1:n + np), sg, p)
            m = m + transpose(reshape([f(P_NE, :), f(P_TE, :), f(P_TI, :), p], [ns, 4]))
            m2 = m2 + transpose(reshape([f(P_NE, :), f(P_TE, :), f(P_TI, :), p], [ns, 4]))**2
        end do
        m = m / nd
        m2 = sqrt(max(m2 / nd - m**2, 0.0_dp))
        open(newunit=u, file=trim(output)//'_profiles.csv', status='replace')
        write(u, '(A)') 's,ne,ne_sd,te,te_sd,ti,ti_sd,p,p_sd'
        do i = 1, ns
            write(u, '(F0.4,8(",",ES24.16))') sg(i), (fm(j, i), m2(j, i), j = 1, 4)
        end do
        close(u)
    end subroutine write_profiles

    subroutine predict(path)
        !! Signals of a given G-EQDSK and the coil currents of the input,
        !! written to <output>_signals.csv, without a fit.
        character(*), intent(in) :: path
        type(geqdsk_plasma_t) :: given
        real(dp), allocatable :: rp(:), zp(:), wp(:), sig(:)
        integer :: unit, k
        call read_geqdsk_plasma(path, given, ierr, message)
        if (ierr /= 0) error stop message
        call given%filaments(rp, zp, wp)
        sig = model%predict(given, model%coil_current)
        if (joint_mode) then
            joint%pp%eq = given
            call write_profile_signals([(0.0_dp, k = 1, joint%pp%n)])
        end if
        write(*, '(A,ES14.6,A)') 'plasma current of the G-EQDSK: ', sum(wp), ' A'
        open(newunit=unit, file=trim(output)//'_signals.csv', status='replace')
        write(unit, '(A)') 'label,measured,predicted,sigma'
        do k = 1, size(sig)
            write(unit, '(A,3(",",ES24.16))') trim(row_label(k)), model%data(k), sig(k), model%sigma(k)
        end do
        close(unit)
    end subroutine predict

    subroutine write_edge(unit, sample, eq)
        integer, intent(in) :: unit, sample
        type(geqdsk_plasma_t), intent(in) :: eq
        real(dp) :: j, r
        integer :: i
        do i = 1, ne
            call eq%midplane_jphi(edge_psin(i), j, r)
            write(unit, '(I0,3(",",ES24.16))') sample, edge_psin(i), r, j
        end do
    end subroutine write_edge

    real(dp) function default_value(name) result(v)
        character(*), intent(in) :: name
        integer :: g
        if (name == 'PLASMA_CURRENT') then
            v = plasma_current
        else if (name == 'PRESSURE_AXIS') then
            v = pressure_axis
        else if (index(name, 'CURRENT(') == 1) then
            read(name(9:index(name, ')') - 1), *) g
            v = coil_current(g)
        else if (index(name, 'PPRIME(') == 1) then
            read(name(8:index(name, ')') - 1), *) g
            v = pprime_knot(g)
        else if (index(name, 'FFPRIME(') == 1) then
            read(name(9:index(name, ')') - 1), *) g
            v = ffprime_knot(g)
        else
            error stop 'unknown parameter '//trim(name)
        end if
    end function default_value

    function row_label(k) result(label)
        integer, intent(in) :: k
        character(len=:), allocatable :: label
        integer :: nf, np, nr, nq
        character(len=32) :: buf
        nf = size(model%machine%sensors%flux)
        np = size(model%machine%sensors%pickup)
        nr = size(model%machine%sensors%rogowski)
        nq = 0
        if (allocated(model%row_kind)) nq = size(model%row_kind)
        if (k > nf + np + nr + nq) then
            label = model%mse(k - nf - np - nr - nq)%label
        else if (k <= nf) then
            label = model%machine%sensors%flux(k)%label
        else if (k <= nf + np) then
            label = model%machine%sensors%pickup(k - nf)%label
        else if (k <= nf + np + nr) then
            label = model%machine%sensors%rogowski(k - nf - np)%label
        else
            write(buf, '(F0.4)') model%row_psin(k - nf - np - nr)
            select case (model%row_kind(k - nf - np - nr))
            case (1)
                label = 'pressure@'//trim(buf)
            case (2)
                label = 'jphi@'//trim(buf)
            case default
                label = 'jdotb@'//trim(buf)
            end select
        end if
    end function row_label

    subroutine read_measurements(path)
        !! Sensor rows in machine order, the profile rows in file order, then
        !! the MSE rows in channel order.
        character(*), intent(in) :: path
        character(len=256) :: line, label
        real(dp) :: value, sigma
        real(dp), allocatable :: rv(:), rs(:), rp(:), mv(:), ms(:)
        integer, allocatable :: rk(:)
        logical, allocatable :: seen(:), mseen(:)
        integer :: unit, ios, i, nsens, at, nm
        nsens = model%machine%sensors%count()
        nm = 0
        if (allocated(model%mse)) nm = size(model%mse)
        allocate(model%data(nsens), model%sigma(nsens), seen(nsens), rv(0), rs(0), rp(0), rk(0), &
            mv(nm), ms(nm), mseen(nm))
        seen = .false.
        mseen = .false.
        open(newunit=unit, file=path, status='old', action='read')
        read(unit, '(A)') line
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            do i = 1, len_trim(line)
                if (line(i:i) == ',') line(i:i) = ' '
            end do
            read(line, *) label, value, sigma
            at = index(label, '@')
            if (at > 0) then
                select case (label(:at - 1))
                case ('pressure')
                    rk = [rk, 1]
                case ('jphi')
                    rk = [rk, 2]
                case ('jdotb')
                    rk = [rk, 3]
                case default
                    error stop 'unknown profile row '//trim(label)
                end select
                read(label(at + 1:), *) value
                rp = [rp, value]
                read(line, *) label, value, sigma
                rv = [rv, value]
                rs = [rs, sigma]
                cycle
            end if
            do i = 1, nm
                if (model%mse(i)%label == trim(label)) exit
            end do
            if (i <= nm) then
                mv(i) = value
                ms(i) = sigma
                mseen(i) = .true.
                cycle
            end if
            i = sensor_index(trim(label))
            if (i == 0) error stop 'measurement of an unknown sensor: '//trim(label)
            model%data(i) = value
            model%sigma(i) = sigma
            seen(i) = .true.
        end do
        close(unit)
        if (.not. all(seen)) error stop 'every sensor of the machine needs a measurement'
        if (size(rk) > 0) then
            model%row_kind = rk
            model%row_psin = rp
            model%data = [model%data, rv]
            model%sigma = [model%sigma, rs]
        end if
        if (.not. all(mseen)) error stop 'every MSE channel needs a measurement'
        model%data = [model%data, mv]
        model%sigma = [model%sigma, ms]
    end subroutine read_measurements

    integer function sensor_index(label) result(i)
        character(*), intent(in) :: label
        integer :: k, nf, np
        nf = size(model%machine%sensors%flux)
        np = size(model%machine%sensors%pickup)
        do k = 1, nf
            i = k
            if (model%machine%sensors%flux(k)%label == label) return
        end do
        do k = 1, np
            i = nf + k
            if (model%machine%sensors%pickup(k)%label == label) return
        end do
        do k = 1, size(model%machine%sensors%rogowski)
            i = nf + np + k
            if (model%machine%sensors%rogowski(k)%label == label) return
        end do
        i = 0
    end function sensor_index

end program tiago_axisym_fit
