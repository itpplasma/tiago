module tiago_joint_model
    !! Joint posterior of equilibrium, kinetic-profile and diagnostic latents
    !! (docs/inference.md). The latent vector is
    !!
    !!     theta = [ x / scale  (reconstruction parameters: equilibrium inputs
    !!                           and coil currents, tiago_reconstruction),
    !!               z_prof     (kinetic-profile latents, N(0, I) prior),
    !!               z_eta      (diagnostic alignment, eta = eta_sigma z_eta),
    !!               z_gain     (calibration gains, g = 1 + gain_sigma z_gain),
    !!               z_cal      (log-scale calibrations, see below) ]
    !!
    !! and the whitened residual r(theta) stacks the magnetic measurements
    !! (whitened with a full covariance), the free-boundary consistency rows and
    !! parameter priors of tiago_reconstruction, the profile diagnostics
    !! L^-1 (g c y - d) with their full covariance, and the latent priors. The
    !! unnormalized log posterior is -|r|^2 / 2.
    !!
    !! Log-scale calibration of a diagnostic slot (cal_sigma > 0): signal i with
    !! cal_group(i) = l > 0 is multiplied by c_i = exp(log C + b_l), signals
    !! with cal_group(i) = 0 are not. The campaign scale is
    !! log C = cal_mean + cal_sigma z_C with z_C ~ N(0, 1), or a unit Student-t
    !! with cal_nu degrees of freedom (heavy tails, prior row
    !! sign(z) sqrt(2 phi(z)) with phi = -log density, FortSym kernel
    !! student_t_prior). With several groups (lasers) and laser_sigma > 0 the
    !! per-group log offsets are b = laser_sigma H z_b, H the orthonormal
    !! Helmert contrasts (L x (L - 1), columns orthogonal to ones), so that
    !! sum_l b_l = 0, each b_l has prior sd laser_sigma sqrt(1 - 1/L), and log C
    !! is the mean log scale of all groups. A slot whose groups are all 1 (e.g.
    !! every interferometer chord) has one shared multiplicative systematic.
    !!
    !! The pressure given to VMEC++ is the kinetic pressure of z_prof on the
    !! half grid (tiago_equilibrium%set_pressure); the magnetic Jacobian with
    !! respect to z_prof follows from VMEC++'s pressure cotangents, and the
    !! profile diagnostics' dependence on the equilibrium from the implicit VJP
    !! of point location in the VMEC++ state (tiago_vmecpp_state, a
    !! real_space_state_t) composed with the same adjoint. The |B| cotangent on
    !! the state's bmnc nodes enters the same adjoint through Tiago's
    !! differentiable copy of VMEC++'s |B| output (vmec_edge%modb). Covariances include
    !! an additive model-error term (emulator or convergence error).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_reconstruction, only: rec, residual_at, jacobian_at, invalidate_caches, n_rows
    use tiago_kinetic_profiles, only: kinetic_profiles_t, N_PROFILES, P_NE, P_TE, P_TI, P_ZEFF
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t, N_QUANTITIES, &
        Q_S, Q_B, Q_NE, Q_TE, Q_TI, Q_ZEFF
    use tiago_vmecpp_state, only: vmecpp_state_t
    use tiago_real_space_state, only: located_t
    use tiago_state_sampler, only: sample_state, sample_state_vjp
    use tiago_student_t_prior, only: student_t_prior
    implicit none
    private

    real(dp), parameter :: mu0 = 4.0e-7_dp * acos(-1.0_dp)
    real(dp), parameter :: log_2pi = log(2.0_dp * acos(-1.0_dp))

    type, public :: diagnostic_slot_t
        class(profile_diagnostic_t), allocatable :: d
        real(dp), allocatable :: data(:)                !! measured signals
        real(dp), allocatable :: covariance(:, :)       !! measurement covariance
        real(dp), allocatable :: model_sigma(:)         !! model-error standard deviation
        real(dp), allocatable :: chol(:, :)             !! of covariance + diag(model_sigma^2)
        real(dp), allocatable :: eta_sigma(:)           !! prior sd of the alignment nuisances
        real(dp) :: gain_sigma = 0.0_dp                 !! > 0: calibration gain latent
        real(dp) :: cal_sigma = 0.0_dp                  !! > 0: log-scale calibration latent
        real(dp) :: cal_mean = 0.0_dp                   !! prior mean of log C
        real(dp) :: cal_nu = 0.0_dp                     !! > 0: Student-t prior on z_C
        real(dp) :: laser_sigma = 0.0_dp                !! > 0: per-group log offsets
        integer, allocatable :: cal_group(:)            !! per signal: 0 or group (laser) 1..L
        integer :: eta_offset = 0, gain_index = 0, row_offset = 0
        integer :: cal_index = 0, laser_offset = 0, n_laser = 0
    end type diagnostic_slot_t

    type, public :: joint_model_t
        type(kinetic_profiles_t) :: prof
        type(diagnostic_slot_t), allocatable :: diags(:)
        real(dp), allocatable :: mag_covariance(:, :)   !! of the magnetic signals (n_meas)
        real(dp), allocatable :: mag_model_sigma(:)
        real(dp), allocatable :: mag_chol(:, :)
        integer :: n_rec = 0, n_prof = 0, n_eta = 0, n_gain = 0, n_cal = 0, n = 0
        integer :: n_mag = 0, n_fb = 0, n_rec_prior = 0, n_diag_rows = 0
        ! evaluation state
        real(dp), allocatable :: s_half(:)
        type(vmecpp_state_t) :: state
        integer :: state_solve = -1
        integer :: evaluations = 0, jacobians = 0
    contains
        procedure :: setup => joint_setup
        procedure :: n_residuals => joint_n_residuals
        procedure :: update_covariances => joint_update_covariances
        procedure :: residual => joint_residual
        procedure :: jacobian => joint_jacobian
        procedure :: log_normalization => joint_log_normalization
        procedure :: proper_prior => joint_proper_prior
        procedure :: unpack => joint_unpack
        procedure :: calibration => joint_calibration
        procedure :: signal_scale => joint_signal_scale
        procedure :: predict_diagnostic => joint_predict_diagnostic
        procedure :: sample_at => joint_sample_at
    end type joint_model_t

contains

    subroutine joint_setup(self)
        !! Latent layout and covariance factors; call after rec, prof and the
        !! diagnostic slots (with data and covariances) are set.
        class(joint_model_t), intent(inout) :: self
        integer :: d, rows

        self%n_rec = rec%n
        self%n_prof = self%prof%n_latent
        self%n_eta = 0
        self%n_gain = 0
        rows = 0
        if (.not. allocated(self%diags)) allocate(self%diags(0))
        do d = 1, size(self%diags)
            associate (slot => self%diags(d))
                slot%eta_offset = self%n_rec + self%n_prof + self%n_eta
                if (.not. allocated(slot%eta_sigma)) allocate(slot%eta_sigma(slot%d%n_nuisance()))
                if (size(slot%eta_sigma) /= slot%d%n_nuisance()) error stop 'joint: eta_sigma size'
                self%n_eta = self%n_eta + size(slot%eta_sigma)
                slot%row_offset = rows
                rows = rows + slot%d%n_signals()
                if (.not. allocated(slot%model_sigma)) then
                    allocate(slot%model_sigma(slot%d%n_signals()))
                    slot%model_sigma = 0.0_dp
                end if
            end associate
        end do
        do d = 1, size(self%diags)
            self%diags(d)%gain_index = 0
            if (self%diags(d)%gain_sigma > 0.0_dp) then
                self%n_gain = self%n_gain + 1
                self%diags(d)%gain_index = self%n_rec + self%n_prof + self%n_eta + self%n_gain
            end if
        end do
        self%n_cal = 0
        do d = 1, size(self%diags)
            associate (slot => self%diags(d))
                slot%cal_index = 0
                slot%laser_offset = 0
                slot%n_laser = 0
                if (slot%cal_sigma > 0.0_dp) then
                    if (.not. allocated(slot%cal_group)) then
                        allocate(slot%cal_group(slot%d%n_signals()))
                        slot%cal_group = 1
                    end if
                    if (size(slot%cal_group) /= slot%d%n_signals()) error stop 'joint: cal_group size'
                    if (any(slot%cal_group < 0)) error stop 'joint: negative cal_group'
                    slot%n_laser = 0
                    if (size(slot%cal_group) > 0) slot%n_laser = maxval(slot%cal_group)
                    if (slot%n_laser == 0) error stop 'joint: calibration without calibrated signals'
                    self%n_cal = self%n_cal + 1
                    slot%cal_index = self%n_rec + self%n_prof + self%n_eta + self%n_gain + self%n_cal
                    if (slot%laser_sigma > 0.0_dp .and. slot%n_laser > 1) then
                        slot%laser_offset = slot%cal_index
                        self%n_cal = self%n_cal + slot%n_laser - 1
                    end if
                end if
            end associate
        end do
        self%n = self%n_rec + self%n_prof + self%n_eta + self%n_gain + self%n_cal
        self%n_mag = rec%n_meas
        self%n_fb = rec%n_fb
        self%n_rec_prior = rec%n_prior
        self%n_diag_rows = rows
        if (.not. allocated(self%mag_covariance)) then
            allocate(self%mag_covariance(self%n_mag, self%n_mag))
            self%mag_covariance = 0.0_dp
            do d = 1, self%n_mag
                self%mag_covariance(d, d) = rec%meas_sigma(d)**2
            end do
        end if
        if (.not. allocated(self%mag_model_sigma)) then
            allocate(self%mag_model_sigma(self%n_mag))
            self%mag_model_sigma = 0.0_dp
        end if
        call self%update_covariances()
        self%s_half = rec%eq%s_half()
        rec%quiet = .true.
    end subroutine joint_setup

    subroutine joint_update_covariances(self)
        !! Cholesky factors of measurement plus model-error covariances.
        class(joint_model_t), intent(inout) :: self
        integer :: d
        self%mag_chol = chol_with_model(self%mag_covariance, self%mag_model_sigma)
        do d = 1, size(self%diags)
            self%diags(d)%chol = chol_with_model(self%diags(d)%covariance, self%diags(d)%model_sigma)
        end do
    end subroutine joint_update_covariances

    function chol_with_model(c, sigma) result(l)
        real(dp), intent(in) :: c(:, :), sigma(:)
        real(dp) :: l(size(c, 1), size(c, 1))
        integer :: i, info
        external :: dpotrf
        l = c
        do i = 1, size(c, 1)
            l(i, i) = l(i, i) + sigma(i)**2
        end do
        if (size(c, 1) == 0) return
        call dpotrf('L', size(c, 1), l, size(c, 1), info)
        if (info /= 0) error stop 'joint: covariance not positive definite'
        do i = 2, size(c, 1)
            l(1:i - 1, i) = 0.0_dp
        end do
    end function chol_with_model

    integer function joint_n_residuals(self)
        class(joint_model_t), intent(in) :: self
        joint_n_residuals = self%n_mag + self%n_fb + self%n_rec_prior + self%n_diag_rows + &
            self%n_prof + self%n_eta + self%n_gain + self%n_cal
    end function joint_n_residuals

    subroutine joint_calibration(self, theta, d, log_c, b)
        !! Campaign log scale log C and per-group log offsets b(n_laser) of
        !! diagnostic d (log C = 0, b = 0 without calibration latent).
        class(joint_model_t), intent(in) :: self
        real(dp), intent(in) :: theta(:)
        integer, intent(in) :: d
        real(dp), intent(out) :: log_c
        real(dp), allocatable, intent(out) :: b(:)
        integer :: l, m
        associate (slot => self%diags(d))
            allocate(b(slot%n_laser))
            log_c = 0.0_dp
            b = 0.0_dp
            if (slot%cal_index == 0) return
            log_c = slot%cal_mean + slot%cal_sigma * theta(slot%cal_index)
            if (slot%laser_offset == 0) return
            do l = 1, slot%n_laser
                do m = 1, slot%n_laser - 1
                    b(l) = b(l) + slot%laser_sigma * helmert(l, m) * theta(slot%laser_offset + m)
                end do
            end do
        end associate
    end subroutine joint_calibration

    function joint_signal_scale(self, theta, d) result(c)
        !! Calibration factor c_i = exp(log C + b_l) of each signal of
        !! diagnostic d (1 for uncalibrated signals).
        class(joint_model_t), intent(in) :: self
        real(dp), intent(in) :: theta(:)
        integer, intent(in) :: d
        real(dp), allocatable :: c(:), b(:)
        real(dp) :: log_c
        integer :: i
        allocate(c(self%diags(d)%d%n_signals()))
        c = 1.0_dp
        if (self%diags(d)%cal_index == 0) return
        call self%calibration(theta, d, log_c, b)
        do i = 1, size(c)
            associate (l => self%diags(d)%cal_group(i))
                if (l > 0) c(i) = exp(log_c + b(l))
            end associate
        end do
    end function joint_signal_scale

    pure real(dp) function helmert(l, m) result(h)
        !! Entry (l, m) of the orthonormal Helmert contrasts: column m is
        !! (1, .., 1, -m, 0, ..) / sqrt(m (m + 1)) with m ones.
        integer, intent(in) :: l, m
        h = 0.0_dp
        if (l <= m) then
            h = 1.0_dp / sqrt(real(m * (m + 1), dp))
        else if (l == m + 1) then
            h = -m / sqrt(real(m * (m + 1), dp))
        end if
    end function helmert

    pure subroutine prior_row(z, nu, r, dr)
        !! Whitened prior row of a latent z: z for N(0, 1) (nu <= 0), else
        !! sign(z) sqrt(2 phi) for the unit Student-t with nu degrees of
        !! freedom, so that r^2 / 2 = phi = -log density + const; near z = 0
        !! the series r = z k (1 - x / 4), x = z^2 / nu, k^2 = (nu + 1) / nu,
        !! replaces the cancellation in log(1 + x).
        real(dp), intent(in) :: z, nu
        real(dp), intent(out) :: r, dr
        real(dp) :: phi, dphi, d2phi, x
        if (nu <= 0.0_dp) then
            r = z
            dr = 1.0_dp
            return
        end if
        x = z * z / nu
        if (x < 1.0e-5_dp) then
            r = z * sqrt((nu + 1.0_dp) / nu) * (1.0_dp - 0.25_dp * x)
            dr = sqrt((nu + 1.0_dp) / nu) * (1.0_dp - 0.75_dp * x)
            return
        end if
        call student_t_prior(z, nu, phi, dphi, d2phi)
        r = sign(sqrt(2.0_dp * phi), z)
        dr = dphi / r
    end subroutine prior_row

    subroutine joint_unpack(self, theta, x, z, eta, gain, d)
        !! Physical reconstruction parameters x, profile latents z, and for
        !! diagnostic d its alignment eta and gain.
        class(joint_model_t), intent(in) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), allocatable, intent(out), optional :: x(:), z(:), eta(:)
        real(dp), intent(out), optional :: gain
        integer, intent(in), optional :: d
        if (present(x)) x = theta(1:self%n_rec) * rec%scale
        if (present(z)) z = theta(self%n_rec + 1:self%n_rec + self%n_prof)
        if (present(d)) then
            associate (slot => self%diags(d))
                if (present(eta)) eta = slot%eta_sigma * &
                    theta(slot%eta_offset + 1:slot%eta_offset + size(slot%eta_sigma))
                if (present(gain)) then
                    gain = 1.0_dp
                    if (slot%gain_index > 0) gain = 1.0_dp + slot%gain_sigma * theta(slot%gain_index)
                end if
            end associate
        end if
    end subroutine joint_unpack

    subroutine prepare(self, theta, r_rec, ok)
        !! Kinetic pressure into VMEC++, magnetic residual rows (solving the
        !! equilibrium if needed) and the flux map of that equilibrium.
        type(joint_model_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), intent(out) :: r_rec(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: x(:), z(:), p(:), wout_b(:), xm(:), xn(:)
        integer :: nyq

        call self%unpack(theta, x=x, z=z)
        allocate(p(size(self%s_half)))
        call self%prof%pressure(z, self%s_half, p)
        if (.not. allocated(rec%eq%pressure_half)) then
            call rec%eq%set_pressure(p)
            call invalidate_caches()
        else if (any(rec%eq%pressure_half /= p)) then
            call rec%eq%set_pressure(p)
            call invalidate_caches()
        end if
        call residual_at(x, r_rec)
        ok = all(ieee_is_finite(r_rec))
        if (.not. ok) return
        if (self%state_solve /= rec%eq%solves) then
            nyq = rec%eq%vmec%get_int('mnmax_nyq')
            xm = rec%eq%vmec%wout('xm_nyq', nyq)
            xn = rec%eq%vmec%wout('xn_nyq', nyq)
            wout_b = rec%eq%vmec%wout('bmnc', nyq * rec%eq%ns)
            call self%state%init(rec%eq%coef, rec%eq%nfp, xm, xn, reshape(wout_b, [nyq, rec%eq%ns]))
            self%state_solve = rec%eq%solves
        end if
    end subroutine prepare

    subroutine joint_sample_at(self, z, x, sample, loc)
        !! Provider fields (tiago_state_sampler) and kinetic profiles at points
        !! x(3, n) = (R, phi, Z) for the current equilibrium.
        class(joint_model_t), intent(in) :: self
        real(dp), intent(in) :: z(:), x(:, :)
        type(plasma_sample_t), intent(out) :: sample
        type(located_t), allocatable, intent(out) :: loc(:)
        real(dp), allocatable :: f(:, :), df(:, :), s(:)
        integer :: j, i
        integer, parameter :: rows(4) = [Q_NE, Q_TE, Q_TI, Q_ZEFF]
        integer, parameter :: profiles(4) = [P_NE, P_TE, P_TI, P_ZEFF]

        call sample_state(self%state, x, sample, loc)
        s = sample%q(Q_S, :)
        allocate(f(N_PROFILES, size(s)), df(N_PROFILES, size(s)))
        call self%prof%evaluate(z, s, f, df)
        do j = 1, size(s)
            do i = 1, 4
                sample%q(rows(i), j) = f(profiles(i), j)
                sample%dq_dx(:, rows(i), j) = df(profiles(i), j) * sample%dq_dx(:, Q_S, j)
            end do
        end do
    end subroutine joint_sample_at

    subroutine joint_predict_diagnostic(self, theta, d, y, sample, loc, eta, dx_deta)
        !! Gain-free prediction y of diagnostic d at the current equilibrium.
        class(joint_model_t), intent(in) :: self
        real(dp), intent(in) :: theta(:)
        integer, intent(in) :: d
        real(dp), intent(out) :: y(:)
        type(plasma_sample_t), intent(out) :: sample
        type(located_t), allocatable, intent(out) :: loc(:)
        real(dp), allocatable, intent(out) :: eta(:), dx_deta(:, :, :)
        real(dp), allocatable :: z(:), x(:, :)

        call self%unpack(theta, z=z, eta=eta, d=d)
        call self%diags(d)%d%sample_points(eta, x, dx_deta)
        call self%sample_at(z, x, sample, loc)
        call self%diags(d)%d%predict(sample, eta, y)
    end subroutine joint_predict_diagnostic

    subroutine joint_residual(self, theta, r, ok)
        !! Whitened residual; ok = .false. if VMEC++ found no equilibrium.
        class(joint_model_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), intent(out) :: r(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: r_rec(:), y(:), eta(:), dx(:, :, :)
        type(plasma_sample_t) :: sample
        type(located_t), allocatable :: loc(:)
        real(dp) :: gain, dr
        integer :: d, row, nm

        self%evaluations = self%evaluations + 1
        allocate(r_rec(n_rows()))
        call prepare(self, theta, r_rec, ok)
        if (.not. ok) then
            r = r_rec(1)
            return
        end if
        nm = self%n_mag
        r(1:nm) = whiten(self%mag_chol, r_rec(1:nm) * rec%meas_sigma(1:nm))
        r(nm + 1:nm + self%n_fb + self%n_rec_prior) = r_rec(nm + 1:)
        row = nm + self%n_fb + self%n_rec_prior
        do d = 1, size(self%diags)
            associate (slot => self%diags(d))
                allocate(y(slot%d%n_signals()))
                call self%predict_diagnostic(theta, d, y, sample, loc, eta, dx)
                call self%unpack(theta, gain=gain, d=d)
                r(row + 1:row + size(y)) = whiten(slot%chol, &
                    gain * self%signal_scale(theta, d) * y - slot%data)
                row = row + size(y)
                deallocate(y)
            end associate
        end do
        r(row + 1:) = theta(self%n_rec + 1:)
        do d = 1, size(self%diags)
            associate (k => self%diags(d)%cal_index)
                if (k > 0) call prior_row(theta(k), self%diags(d)%cal_nu, r(row + k - self%n_rec), dr)
            end associate
        end do
        ok = all(ieee_is_finite(r))
    end subroutine joint_residual

    subroutine joint_jacobian(self, theta, jac)
        !! d r / d theta at theta (after residual at the same theta or not).
        class(joint_model_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), allocatable, intent(out) :: jac(:, :)
        real(dp), allocatable :: r_rec(:), jrec(:, :), z(:), zbar(:), y(:), eta(:), dx(:, :, :)
        real(dp), allocatable :: q_bar(:, :), eta_bar(:), f_bar(:, :), s_bar(:), p_bar(:), qs_bar(:, :)
        real(dp), allocatable :: coef_cot(:, :), bmnc_cot(:, :), cot(:, :), jeq(:, :), jd(:, :), e(:)
        real(dp), allocatable :: f(:), c(:)
        type(plasma_sample_t) :: sample
        type(located_t), allocatable :: loc(:)
        real(dp) :: gain, rk, drk
        integer :: d, i, j, k, q, row, nm, nrec_rows, ndr, col, ng, l, m
        logical :: ok
        integer, parameter :: rows(4) = [Q_NE, Q_TE, Q_TI, Q_ZEFF]
        integer, parameter :: profiles(4) = [P_NE, P_TE, P_TI, P_ZEFF]

        self%jacobians = self%jacobians + 1
        allocate(jac(self%n_residuals(), self%n))
        jac = 0.0_dp
        allocate(r_rec(n_rows()))
        call prepare(self, theta, r_rec, ok)
        if (.not. ok) error stop 'joint_jacobian: no equilibrium at theta'
        call self%unpack(theta, z=z)
        nm = self%n_mag
        nrec_rows = nm + self%n_fb + self%n_rec_prior

        ! reconstruction rows: parameters (scaled), kinetic pressure through mu0 p_h
        call jacobian_at(theta(1:self%n_rec) * rec%scale, jrec)
        jac(1:nrec_rows, 1:self%n_rec) = jrec * spread(rec%scale, 1, nrec_rows)
        if (.not. allocated(rec%eq%pressure_bar)) error stop 'joint: needs an equilibrium parameter'
        if (size(rec%eq%pressure_bar, 2) /= nm + self%n_fb) error stop 'joint: pressure cotangent rows'
        allocate(zbar(self%n_prof))
        do i = 1, nm + self%n_fb
            zbar = 0.0_dp
            call self%prof%pressure_vjp(z, self%s_half, mu0 * rec%eq%pressure_bar(:, i), zbar)
            jac(i, self%n_rec + 1:self%n_rec + self%n_prof) = zbar
        end do
        ! magnetic rows: undo the diagonal whitening, apply the full covariance
        do col = 1, self%n_rec + self%n_prof
            jac(1:nm, col) = whiten(self%mag_chol, jac(1:nm, col) * rec%meas_sigma(1:nm))
        end do

        ! profile-diagnostic rows
        ndr = self%n_diag_rows
        if (ndr > 0) then
            ng = self%state%n_geometry()
            allocate(coef_cot(ng, ndr), p_bar(self%state%n_params()))
            allocate(bmnc_cot(self%state%n_params() - ng, ndr))
            row = 0
            do d = 1, size(self%diags)
                associate (slot => self%diags(d))
                    allocate(y(slot%d%n_signals()), e(slot%d%n_signals()))
                    call self%predict_diagnostic(theta, d, y, sample, loc, eta, dx)
                    call self%unpack(theta, gain=gain, d=d)
                    ! total factor of each signal: gain times calibration scale
                    c = self%signal_scale(theta, d)
                    f = gain * c
                    allocate(jd(size(y), self%n), q_bar(N_QUANTITIES, sample%n), &
                        eta_bar(size(eta)), f_bar(N_PROFILES, sample%n), s_bar(sample%n), &
                        qs_bar(N_QUANTITIES, sample%n))
                    jd = 0.0_dp
                    do i = 1, size(y)
                        e = 0.0_dp
                        e(i) = 1.0_dp
                        call slot%d%vjp(sample, eta, e, q_bar, eta_bar)
                        ! profiles: latents and s through the profiles
                        f_bar = 0.0_dp
                        do k = 1, 4
                            f_bar(profiles(k), :) = q_bar(rows(k), :)
                        end do
                        zbar = 0.0_dp
                        call self%prof%vjp(z, sample%q(Q_S, :), f_bar, zbar, s_bar)
                        jd(i, self%n_rec + 1:self%n_rec + self%n_prof) = f(i) * zbar
                        s_bar = s_bar + q_bar(Q_S, :)
                        ! alignment: direct part plus motion of the points
                        do k = 1, size(eta)
                            do j = 1, sample%n
                                do q = 1, N_QUANTITIES
                                    eta_bar(k) = eta_bar(k) + q_bar(q, j) * &
                                        dot_product(sample%dq_dx(:, q, j), dx(:, j, k))
                                end do
                            end do
                        end do
                        jd(i, slot%eta_offset + 1:slot%eta_offset + size(eta)) = &
                            f(i) * eta_bar * slot%eta_sigma
                        if (slot%gain_index > 0) jd(i, slot%gain_index) = &
                            slot%gain_sigma * c(i) * y(i)
                        l = 0
                        if (slot%cal_index > 0) l = slot%cal_group(i)
                        if (l > 0) then
                            jd(i, slot%cal_index) = slot%cal_sigma * f(i) * y(i)
                            if (slot%laser_offset > 0) then
                                do m = 1, slot%n_laser - 1
                                    jd(i, slot%laser_offset + m) = &
                                        slot%laser_sigma * helmert(l, m) * f(i) * y(i)
                                end do
                            end if
                        end if
                        ! geometry and |B| (bmnc nodes) at fixed points through
                        ! the state's VJP, then the equilibrium adjoint
                        qs_bar = 0.0_dp
                        qs_bar(Q_S, :) = f(i) * s_bar
                        qs_bar(Q_B, :) = f(i) * q_bar(Q_B, :)
                        p_bar = 0.0_dp
                        call sample_state_vjp(self%state, loc, qs_bar, p_bar)
                        coef_cot(:, row + i) = p_bar(1:ng)
                        bmnc_cot(:, row + i) = p_bar(ng + 1:)
                    end do
                    jac(nrec_rows + row + 1:nrec_rows + row + size(y), :) = jd
                    row = row + size(y)
                    deallocate(y, e, jd, q_bar, eta_bar, f_bar, s_bar, qs_bar)
                end associate
            end do
            ! equilibrium adjoint for all diagnostic rows at once
            allocate(cot(ndr, rec%eq%n_y()))
            cot = 0.0_dp
            call rec%eq%jacobian(rec%eq_x, cot, jeq, coef_cot, bmnc_cot)
            do k = 1, self%n_rec
                if (rec%is_coil(k)) cycle
                jac(nrec_rows + 1:nrec_rows + ndr, k) = jac(nrec_rows + 1:nrec_rows + ndr, k) + &
                    jeq(:, rec%eq_index(k)) * rec%scale(k)
            end do
            do i = 1, ndr
                zbar = 0.0_dp
                call self%prof%pressure_vjp(z, self%s_half, mu0 * rec%eq%pressure_bar(:, i), zbar)
                jac(nrec_rows + i, self%n_rec + 1:self%n_rec + self%n_prof) = &
                    jac(nrec_rows + i, self%n_rec + 1:self%n_rec + self%n_prof) + zbar
            end do
            ! whiten each diagnostic block
            do d = 1, size(self%diags)
                associate (slot => self%diags(d))
                    row = nrec_rows + slot%row_offset
                    do col = 1, self%n
                        jac(row + 1:row + slot%d%n_signals(), col) = &
                            whiten(slot%chol, jac(row + 1:row + slot%d%n_signals(), col))
                    end do
                end associate
            end do
        end if
        ! latent priors
        row = nrec_rows + ndr
        do k = self%n_rec + 1, self%n
            jac(row + k - self%n_rec, k) = 1.0_dp
        end do
        do d = 1, size(self%diags)
            k = self%diags(d)%cal_index
            if (k == 0) cycle
            call prior_row(theta(k), self%diags(d)%cal_nu, rk, drk)
            jac(row + k - self%n_rec, k) = drk
        end do
    end subroutine joint_jacobian

    real(dp) function joint_log_normalization(self) result(c)
        !! log p(d, theta) = -|r|^2 / 2 + c for the latent coordinates theta:
        !! Gaussian normalizations of all data rows (magnetic, consistency,
        !! profile) and of all priors (reconstruction priors in scaled units).
        class(joint_model_t), intent(in) :: self
        integer :: d, i, k
        c = -0.5_dp * self%n_residuals() * log_2pi
        do i = 1, self%n_mag
            c = c - log(self%mag_chol(i, i))
        end do
        c = c - self%n_fb * log(rec%fb_sigma)
        do i = 1, rec%n_prior
            k = rec%prior_param(i)
            c = c - log(rec%prior_sigma(k) / rec%scale(k))
        end do
        do d = 1, size(self%diags)
            do i = 1, size(self%diags(d)%chol, 1)
                c = c - log(self%diags(d)%chol(i, i))
            end do
            ! Student-t calibration prior instead of the N(0, 1) normalization
            associate (nu => self%diags(d)%cal_nu)
                if (self%diags(d)%cal_index > 0 .and. nu > 0.0_dp) c = c + 0.5_dp * log_2pi + &
                    log_gamma(0.5_dp * (nu + 1.0_dp)) - log_gamma(0.5_dp * nu) - &
                    0.5_dp * log(nu * acos(-1.0_dp))
            end associate
        end do
    end function joint_log_normalization

    logical function joint_proper_prior(self)
        !! Every reconstruction parameter has a prior (evidence well defined).
        class(joint_model_t), intent(in) :: self
        joint_proper_prior = rec%n_prior == rec%n
    end function joint_proper_prior

    function whiten(l, v) result(w)
        !! w = L^-1 v for lower-triangular L.
        real(dp), intent(in) :: l(:, :), v(:)
        real(dp) :: w(size(v))
        integer :: i
        do i = 1, size(v)
            w(i) = (v(i) - dot_product(l(i, 1:i - 1), w(1:i - 1))) / l(i, i)
        end do
    end function whiten

end module tiago_joint_model
