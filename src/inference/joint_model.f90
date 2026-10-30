module tiago_joint_model
    !! Joint posterior of equilibrium, kinetic-profile and diagnostic latents
    !! (docs/inference.md). The latent vector is
    !!
    !!     theta = [ x / scale  (reconstruction parameters: equilibrium inputs
    !!                           and coil currents, tiago_reconstruction),
    !!               z_prof     (kinetic-profile latents, N(0, I) prior),
    !!               z_eta      (diagnostic alignment, eta = eta_factor z_eta),
    !!               z_gain     (calibration gains and magnetic integrator),
    !!               z_delta    (constrained instrument discrepancy) ]
    !!
    !! and the whitened residual r(theta) stacks the magnetic measurements
    !! (whitened with a full covariance), the free-boundary consistency rows and
    !! parameter priors of tiago_reconstruction, the profile diagnostics
    !! L^-1 (g y - d) with their full covariance, and the latent priors. The
    !! unnormalized log posterior is -|r|^2 / 2.
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
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use tiago_reconstruction, only: rec, residual_at, jacobian_at, invalidate_caches, n_rows
    use tiago_kinetic_profiles, only: kinetic_profiles_t, N_PROFILES, P_NE, P_TE, P_TI, P_ZEFF
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t, N_QUANTITIES, &
        Q_S, Q_B, Q_NE, Q_TE, Q_TI, Q_ZEFF
    use tiago_vmecpp_state, only: vmecpp_state_t
    use tiago_ece_radiometer, only: ece_t
    use tiago_ece_recorded, only: recorded_ece_t, record_ece
    use tiago_ece_state_medium, only: state_medium_t
    use fortnum_status, only: fortnum_status_t
    use fortnum_multifidelity, only: autoregressive_gp_t
    use tiago_calibration_campaign, only: constrain_discrepancy
    use tiago_calibration, only: calibration_prior_t
    use tiago_integrator_signal, only: integrator_signal
    use tiago_robust_signals, only: robust_rows, robust_normalization
    use tiago_real_space_state, only: located_t
    use tiago_state_sampler, only: sample_state, sample_state_vjp
    implicit none
    private

    real(dp), parameter :: mu0 = 4.0e-7_dp * acos(-1.0_dp)
    real(dp), parameter :: log_2pi = log(2.0_dp * acos(-1.0_dp))

    type, public :: diagnostic_slot_t
        class(profile_diagnostic_t), allocatable :: d
        real(dp), allocatable :: data(:) !! measured signals
        real(dp), allocatable :: covariance(:, :) !! measurement covariance
        real(dp), allocatable :: model_sigma(:) !! model-error standard deviation
        real(dp), allocatable :: chol(:, :) !! of covariance + diag(model_sigma^2)
        real(dp), allocatable :: eta_sigma(:) !! prior sd of the alignment nuisances
        real(dp), allocatable :: eta_mean(:)
        real(dp), allocatable :: discrepancy_basis(:, :), protected_jacobian(:, :)
        integer :: discrepancy_offset = 0
        logical :: low_fidelity_te = .false.
        real(dp) :: student_nu = 0.0_dp
        real(dp) :: gain_sigma = 0.0_dp !! > 0: calibration gain latent
        integer :: eta_offset = 0, gain_index = 0, row_offset = 0
    end type diagnostic_slot_t

    type, public :: joint_model_t
        type(kinetic_profiles_t) :: prof
        type(diagnostic_slot_t), allocatable :: diags(:)
        real(dp), allocatable :: mag_covariance(:, :) !! of the magnetic signals (n_meas)
        real(dp), allocatable :: mag_model_sigma(:)
        real(dp), allocatable :: mag_chol(:, :)
        real(dp), allocatable :: eta_factor(:, :)
        real(dp) :: mag_student_nu = 0.0_dp
        type(calibration_prior_t) :: mag_calibration
        real(dp), allocatable :: mag_time(:)
        real(dp) :: mag_time_origin = 0
        integer :: mag_calibration_index = 0
        type(autoregressive_gp_t) :: surrogate
        real(dp), allocatable :: surrogate_origin(:), surrogate_base(:), surrogate_jac(:, :)
        real(dp), allocatable :: surrogate_lower(:), surrogate_upper(:)
        logical :: use_surrogate = .false.
        integer :: n_discrepancy = 0
        integer :: n_rec = 0, n_prof = 0, n_eta = 0, n_gain = 0, n = 0
        integer :: n_mag = 0, n_fb = 0, n_rec_prior = 0, n_diag_rows = 0
        ! evaluation state
        real(dp), allocatable :: s_half(:)
        type(vmecpp_state_t) :: state
        integer :: state_solve = -1
        integer :: evaluations = 0, jacobians = 0
    contains
        procedure :: setup => joint_setup
        procedure :: train_surrogate => joint_train_surrogate
        procedure :: jvp => joint_jvp
        procedure :: vjp => joint_vjp
        procedure :: n_residuals => joint_n_residuals
        procedure :: update_covariances => joint_update_covariances
        procedure :: residual => joint_residual
        procedure :: jacobian => joint_jacobian
        procedure :: log_normalization => joint_log_normalization
        procedure :: proper_prior => joint_proper_prior
        procedure :: unpack => joint_unpack
        procedure :: predict_diagnostic => joint_predict_diagnostic
        procedure :: sample_at => joint_sample_at
    end type joint_model_t

contains

    subroutine joint_setup(self)
        !! Latent layout and covariance factors; call after rec, prof and the
        !! diagnostic slots (with data and covariances) are set.
        class(joint_model_t), intent(inout) :: self
        integer :: d, rows, k

        self%state_solve = -1
        self%use_surrogate = .false.
        if (.not. ieee_is_finite(self%mag_student_nu) .or. self%mag_student_nu < 0) &
            error stop 'joint: invalid magnetic Student degrees of freedom'
        self%n_rec = rec%n
        self%n_prof = self%prof%n_latent
        self%n_eta = 0
        self%n_gain = 0
        rows = 0
        if (.not. allocated(self%diags)) allocate(self%diags(0))
        do d = 1, size(self%diags)
            associate (slot => self%diags(d))
                if (self%prof%te_ar1%active) then
                    select type (instrument => slot%d)
                        type is (ece_t)
                        slot%low_fidelity_te = .true.
                    end select
                    if (slot%low_fidelity_te .and. slot%gain_sigma > 0) &
                        error stop 'joint: AR1 scale replaces the ECE gain nuisance'
                else if (slot%low_fidelity_te) then
                    error stop 'joint: low-fidelity Te requires an AR1 profile'
                end if
                if (.not. ieee_is_finite(slot%student_nu) .or. slot%student_nu < 0) &
                    error stop 'joint: invalid Student degrees of freedom'
                slot%eta_offset = self%n_rec + self%n_prof + self%n_eta
                if (.not. allocated(slot%eta_sigma)) then
                    allocate (slot%eta_sigma(slot%d%n_nuisance()))
                    slot%eta_sigma = 1
                end if
                if (size(slot%eta_sigma) /= slot%d%n_nuisance()) error stop 'joint: eta_sigma size'
                if (.not. allocated(slot%eta_mean)) then
                    allocate (slot%eta_mean(size(slot%eta_sigma)))
                    slot%eta_mean = 0
                end if
                if (size(slot%eta_mean) /= size(slot%eta_sigma)) error stop 'joint: eta_mean size'
                if (.not. all(ieee_is_finite(slot%eta_mean))) error stop 'joint: nonfinite eta_mean'
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
        self%mag_calibration_index = 0
        if (self%mag_calibration%n_params() > 0) then
            if (self%mag_calibration%n_params() /= 3) error stop 'joint: magnetic calibration count'
            if (rec%n_meas < 1) error stop 'joint: magnetic calibration requires signals'
            if (.not. allocated(self%mag_calibration%factor)) &
                error stop 'joint: magnetic calibration requires a prior factor'
            if (any(rec%meas_kind /= rec%meas_kind(1))) &
                error stop 'joint: common integrator calibration requires a common signal unit'
            self%mag_calibration_index = self%n_rec + self%n_prof + self%n_eta + self%n_gain
            self%n_gain = self%n_gain + 3
            if (.not. allocated(self%mag_time)) error stop 'joint: magnetic recording times required'
            if (size(self%mag_time) /= rec%n_meas) error stop 'joint: magnetic time count'
            if (.not. all(ieee_is_finite(self%mag_time)) .or. &
                .not. ieee_is_finite(self%mag_time_origin)) error stop 'joint: magnetic clock'
            if (any(self%mag_time < self%mag_calibration%valid_from) .or. &
                any(self%mag_time > self%mag_calibration%valid_to)) &
                error stop 'joint: magnetic calibration validity interval'
        end if
        self%n = self%n_rec + self%n_prof + self%n_eta + self%n_gain
        if (allocated(self%eta_factor)) then
            if (any(shape(self%eta_factor) /= [self%n_eta, self%n_eta])) deallocate (self%eta_factor)
        end if
        if (.not. allocated(self%eta_factor)) then
            allocate (self%eta_factor(self%n_eta, self%n_eta))
            self%eta_factor = 0
            do d = 1, size(self%diags)
                associate (slot => self%diags(d))
                    do k = 1, size(slot%eta_sigma)
                        self%eta_factor(slot%eta_offset - self%n_rec - self%n_prof + k, &
                            slot%eta_offset - self%n_rec - self%n_prof + k) = slot%eta_sigma(k)
                    end do
                end associate
            end do
        end if
        if (any(shape(self%eta_factor) /= [self%n_eta, self%n_eta])) &
            error stop 'joint: eta_factor shape'
        if (.not. all(ieee_is_finite(self%eta_factor))) error stop 'joint: nonfinite eta_factor'

        self%n_discrepancy = 0
        do d = 1, size(self%diags)
            associate (slot => self%diags(d))
                slot%discrepancy_offset = self%n + self%n_discrepancy
                if (.not. allocated(slot%discrepancy_basis)) &
                    allocate (slot%discrepancy_basis(slot%d%n_signals(), 0))
                self%n_discrepancy = self%n_discrepancy + size(slot%discrepancy_basis, 2)
            end associate
        end do
        self%n = self%n + self%n_discrepancy
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
        integer :: d, code
        real(dp), allocatable :: constrained(:, :)
        self%use_surrogate = .false.
        if (allocated(self%surrogate_origin)) then
            deallocate(self%surrogate_origin,self%surrogate_base,self%surrogate_jac, &
                self%surrogate_lower,self%surrogate_upper)
        end if
        self%mag_chol = chol_with_model(self%mag_covariance, self%mag_model_sigma)
        do d = 1, size(self%diags)
            self%diags(d)%chol = chol_with_model(self%diags(d)%covariance, self%diags(d)%model_sigma)
        end do
        do d = 1, size(self%diags)
            associate (slot => self%diags(d))
                if (size(slot%discrepancy_basis, 2) == 0) cycle
                if (.not. allocated(slot%protected_jacobian)) &
                    error stop 'joint: discrepancy requires protected profile/calibration rows'
                call constrain_discrepancy(slot%chol, slot%protected_jacobian, &
                    slot%discrepancy_basis, constrained, code)
                if (code /= 0) error stop 'joint: discrepancy constraint failed'
                slot%discrepancy_basis = constrained
            end associate
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
            self%n_prof + self%n_eta + self%n_gain + self%n_discrepancy
    end function joint_n_residuals

    subroutine joint_unpack(self, theta, x, z, eta, gain, d)
        !! Physical reconstruction parameters x, profile latents z, and for
        !! diagnostic d its alignment eta and gain.
        class(joint_model_t), intent(in) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), allocatable, intent(out), optional :: x(:), z(:), eta(:)
        real(dp), intent(out), optional :: gain
        integer, intent(in), optional :: d
        integer :: first,last
        if (present(x)) x = theta(1:self%n_rec) * rec%scale
        if (present(z)) z = theta(self%n_rec + 1:self%n_rec + self%n_prof)
        if (present(d)) then
            associate (slot => self%diags(d))
                first=slot%eta_offset-self%n_rec-self%n_prof+1
                last=first+size(slot%eta_sigma)-1
                if(present(eta)) eta=slot%eta_mean+matmul(self%eta_factor(first:last,:), &
                    theta(self%n_rec+self%n_prof+1:self%n_rec+self%n_prof+self%n_eta))
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
        real(dp), allocatable :: x(:), z(:), p(:), wout_b(:), xm(:), xn(:), bu(:), bv(:)
        integer :: nyq

        call self%unpack(theta, x=x, z=z)
        allocate(p(size(self%s_half)))
        call self%prof%pressure(z, self%s_half, p)
        ok = all(ieee_is_finite(p))
        if (.not. ok) then
            r_rec = huge(1._dp)
            return
        end if
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
            bu = rec%eq%vmec%wout('bsupumnc', nyq*rec%eq%ns)
            bv = rec%eq%vmec%wout('bsupvmnc', nyq*rec%eq%ns)
            call self%state%init(rec%eq%coef, rec%eq%nfp, xm, xn, &
                reshape(wout_b, [nyq, rec%eq%ns]), reshape(bu, [nyq, rec%eq%ns]), &
                reshape(bv, [nyq, rec%eq%ns]))
            self%state_solve = rec%eq%solves
        end if
    end subroutine prepare

    subroutine joint_sample_at(self, z, x, sample, loc, low_fidelity_te)
        !! Provider fields (tiago_state_sampler) and kinetic profiles at points
        !! x(3, n) = (R, phi, Z) for the current equilibrium.
        class(joint_model_t), intent(in) :: self
        real(dp), intent(in) :: z(:), x(:, :)
        type(plasma_sample_t), intent(out) :: sample
        type(located_t), allocatable, intent(out) :: loc(:)
        real(dp), allocatable :: f(:, :), df(:, :), s(:)
        integer :: j, i
        logical, intent(in), optional :: low_fidelity_te
        type(kinetic_profiles_t) :: profile_view
        integer, parameter :: rows(4) = [Q_NE, Q_TE, Q_TI, Q_ZEFF]
        integer, parameter :: profiles(4) = [P_NE, P_TE, P_TI, P_ZEFF]

        call sample_state(self%state, x, sample, loc)
        s = sample%q(Q_S, :)
        allocate(f(N_PROFILES, size(s)), df(N_PROFILES, size(s)))
        profile_view = self%prof
        if (present(low_fidelity_te)) profile_view%low_fidelity_te = low_fidelity_te
        call profile_view%evaluate(z, s, f, df)
        do j = 1, size(s)
            do i = 1, 4
                sample%q(rows(i), j) = f(profiles(i), j)
                sample%dq_dx(:, rows(i), j) = df(profiles(i), j) * sample%dq_dx(:, Q_S, j)
            end do
        end do
    end subroutine joint_sample_at

    subroutine joint_predict_diagnostic(self, theta, d, y, sample, loc, eta, dx_deta, tape)
        !! Gain-free prediction y of diagnostic d at the current equilibrium.
        class(joint_model_t), intent(in) :: self
        real(dp), intent(in) :: theta(:)
        integer, intent(in) :: d
        real(dp), intent(out) :: y(:)
        type(plasma_sample_t), intent(out) :: sample
        type(located_t), allocatable, intent(out) :: loc(:)
        real(dp), allocatable, intent(out) :: eta(:), dx_deta(:, :, :)
        real(dp), allocatable :: z(:), x(:, :)
        type(recorded_ece_t), intent(out), optional :: tape
        type(recorded_ece_t) :: local_tape
        type(state_medium_t) :: medium
        type(ece_t) :: prepared
        type(fortnum_status_t) :: status

        call self%unpack(theta, z=z, eta=eta, d=d)
        select type (instrument => self%diags(d)%d)
            type is (ece_t)
            allocate (medium%state, source=self%state)
            medium%profiles = self%prof
            medium%profiles%low_fidelity_te = self%diags(d)%low_fidelity_te
            medium%z = z
            call record_ece(instrument, medium, eta, local_tape, status)
            if (status%code /= 0) then
                y = ieee_value(0._dp, ieee_quiet_nan)
                allocate (loc(0), dx_deta(3, 0, size(eta)))
                return
            end if
            y = local_tape%voltage
            call local_tape%sampled(prepared, sample)
            allocate (loc(0), dx_deta(3, sample%n, size(eta)))
            dx_deta = 0
            if (present(tape)) tape = local_tape
            return
        end select
        call self%diags(d)%d%sample_points(eta, x, dx_deta)
        call self%sample_at(z, x, sample, loc, self%diags(d)%low_fidelity_te)
        call self%diags(d)%d%predict(sample, eta, y)
    end subroutine joint_predict_diagnostic

    subroutine joint_residual(self, theta, r, ok)
        !! Whitened residual; ok = .false. if VMEC++ found no equilibrium.
        class(joint_model_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), intent(out) :: r(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: r_rec(:), y(:), eta(:), dx(:, :, :), low(:), sj(:, :)
        real(dp) :: variance
        type(fortnum_status_t) :: status
        type(plasma_sample_t) :: sample
        type(located_t), allocatable :: loc(:)
        real(dp) :: gain
        integer :: d, row, nm, i, k

        if (self%use_surrogate) then
            ok = .false.
            if (.not. allocated(self%surrogate_origin)) return
            if (any(theta < self%surrogate_lower) .or. any(theta > self%surrogate_upper)) return
            low = self%surrogate_base + matmul(self%surrogate_jac, theta - self%surrogate_origin)
            allocate (sj(self%n_residuals(), self%n))
            call self%surrogate%predict(theta, low, self%surrogate_jac, r, sj, variance, status)
            ok = status%code == 0
            do i = 1, self%n_rec_prior
                k = rec%prior_param(i)
                r(self%n_mag + self%n_fb + i) = &
                    (theta(k)*rec%scale(k) - rec%prior(k))/rec%prior_sigma(k)
            end do
            row = self%n_residuals() - (self%n - self%n_rec)
            r(row + 1:) = theta(self%n_rec + 1:)
            return
        end if
        self%evaluations = self%evaluations + 1
        allocate(r_rec(n_rows()))
        call prepare(self, theta, r_rec, ok)
        if (.not. ok) then
            r = r_rec(1)
            return
        end if
        nm = self%n_mag
        r(1:nm) = magnetic_residual(self, theta, r_rec(:nm))
        r(nm + 1:nm + self%n_fb + self%n_rec_prior) = r_rec(nm + 1:)
        row = nm + self%n_fb + self%n_rec_prior
        do d = 1, size(self%diags)
            associate (slot => self%diags(d))
                allocate(y(slot%d%n_signals()))
                call self%predict_diagnostic(theta, d, y, sample, loc, eta, dx)
                call self%unpack(theta, gain=gain, d=d)
                y = gain*y - slot%data + matmul(slot%discrepancy_basis, &
                    theta(slot%discrepancy_offset+1: &
                    slot%discrepancy_offset+size(slot%discrepancy_basis,2)))
                r(row + 1:row + size(y)) = whiten(slot%chol, y)
                row = row + size(y)
                deallocate(y)
            end associate
        end do
        call robust_rows(r(:nm), self%mag_student_nu)
        do d = 1, size(self%diags)
            associate (slot => self%diags(d))
                row = nm + self%n_fb + self%n_rec_prior + slot%row_offset
                call robust_rows(r(row + 1:row + slot%d%n_signals()), slot%student_nu)
            end associate
        end do
        row = nm + self%n_fb + self%n_rec_prior + self%n_diag_rows
        r(row + 1:) = theta(self%n_rec + 1:)
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
        type(plasma_sample_t) :: sample
        type(kinetic_profiles_t) :: observation_profiles
        type(recorded_ece_t) :: tape
        type(fortnum_status_t) :: status
        type(located_t), allocatable :: loc(:)
        real(dp), allocatable :: ub_cot(:, :), vb_cot(:, :), pb_all(:, :), eb_all(:, :), r0(:), slope(:)
        real(dp) :: gain
        integer :: d, i, j, k, q, row, nm, nrec_rows, ndr, col, ng, nb
        logical :: ok
        integer, parameter :: rows(4) = [Q_NE, Q_TE, Q_TI, Q_ZEFF]
        integer, parameter :: profiles(4) = [P_NE, P_TE, P_TI, P_ZEFF]

        if (self%use_surrogate) then
            if (.not. allocated(self%surrogate_origin)) error stop 'joint: surrogate not trained'
            if (any(theta < self%surrogate_lower) .or. any(theta > self%surrogate_upper)) &
                error stop 'joint: surrogate extrapolation'
            allocate (jac(self%n_residuals(), self%n), r0(self%n_residuals()))
            r0 = self%surrogate_base + matmul(self%surrogate_jac, theta - self%surrogate_origin)
            y = r0
            call self%surrogate%predict(theta, y, self%surrogate_jac, r0, jac, gain, status)
            if (status%code /= 0) error stop 'joint: surrogate prediction failed'
            do i = 1, self%n_rec_prior
                k = rec%prior_param(i)
                row = self%n_mag + self%n_fb + i
                jac(row, :) = 0
                jac(row, k) = rec%scale(k)/rec%prior_sigma(k)
            end do
            row = self%n_residuals() - (self%n - self%n_rec)
            jac(row + 1:, :) = 0
            do k = self%n_rec + 1, self%n
                jac(row + k - self%n_rec, k) = 1
            end do
            return
        end if
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
        if (self%mag_calibration_index > 0) then
            block
                real(dp):: c(3), raw(self%n_mag), v, dr, dg, db, dd
                real(dp):: cb(3), cjac(self%n_mag, 3)
                integer:: mi, mk
                mi = self%mag_calibration_index
                c = self%mag_calibration%physical(theta(mi + 1:mi + 3))
                raw = r_rec(:nm)*rec%meas_sigma(:nm) + rec%meas_value(:nm)
                do mk = 1, nm
                    call integrator_signal(raw(mk), c(1), c(2), c(3), self%mag_time(mk), &
                        self%mag_time_origin, v, dr, dg, db, dd)
                    jac(mk, :) = dr*jac(mk, :)
                    cb = [dg, db, dd]
                    cjac(mk, :) = matmul(cb, self%mag_calibration%factor)/rec%meas_sigma(mk)
                end do
                jac(:nm, mi + 1:mi + 3) = cjac
            end block
        end if
        ! magnetic rows: undo the diagonal whitening, apply the full covariance
        do col = 1, self%n
            jac(1:nm, col) = whiten(self%mag_chol, jac(1:nm, col) * rec%meas_sigma(1:nm))
        end do

        ! profile-diagnostic rows
        ndr = self%n_diag_rows
        if (ndr > 0) then
            ng = self%state%n_geometry()
            allocate(coef_cot(ng, ndr), p_bar(self%state%n_params()))
            nb = size(self%state%bmnc)
            allocate (bmnc_cot(nb, ndr), ub_cot(nb, ndr), vb_cot(nb, ndr))
            ub_cot = 0; vb_cot = 0
            row = 0
            do d = 1, size(self%diags)
                associate (slot => self%diags(d))
                    observation_profiles = self%prof
                    observation_profiles%low_fidelity_te = slot%low_fidelity_te
                    allocate(y(slot%d%n_signals()), e(slot%d%n_signals()))
                    call self%predict_diagnostic(theta, d, y, sample, loc, eta, dx, tape)
                    call self%unpack(theta, gain=gain, d=d)
                    allocate(jd(size(y), self%n), q_bar(N_QUANTITIES, sample%n), &
                        eta_bar(size(eta)), f_bar(N_PROFILES, sample%n), s_bar(sample%n), &
                        qs_bar(N_QUANTITIES, sample%n))
                    jd = 0.0_dp
                    select type (instrument => slot%d)
                        type is (ece_t)
                        allocate (pb_all(self%state%n_params() + self%n_prof, size(y)), &
                            eb_all(size(eta), size(y)))
                        call ece_rows(tape, pb_all, eb_all)
                    end select
                    do i = 1, size(y)
                        e = 0.0_dp
                        e(i) = 1.0_dp
                        select type (instrument => slot%d)
                            type is (ece_t)
                            p_bar = pb_all(:self%state%n_params(), i)
                            zbar = pb_all(self%state%n_params() + 1:, i)
                            eta_bar = eb_all(:, i)
                            jd(i, self%n_rec + 1:self%n_rec + self%n_prof) = gain*zbar
                            jd(i, self%n_rec + self%n_prof + 1:self%n_rec + self%n_prof + self%n_eta) = &
                                gain*matmul(eta_bar, self%eta_factor( &
                                slot%eta_offset - self%n_rec - self%n_prof + 1: &
                                slot%eta_offset - self%n_rec - self%n_prof + size(eta), :))
                            if (slot%gain_index > 0) jd(i, slot%gain_index) = slot%gain_sigma*y(i)
                            coef_cot(:, row + i) = gain*p_bar(:ng)
                            bmnc_cot(:, row + i) = gain*p_bar(ng + 1:ng + nb)
                            ub_cot(:, row + i) = gain*p_bar(ng + nb + 1:ng + 2*nb)
                            vb_cot(:, row + i) = gain*p_bar(ng + 2*nb + 1:ng + 3*nb)
                            cycle
                        end select
                        call slot%d%vjp(sample, eta, e, q_bar, eta_bar)
                        ! profiles: latents and s through the profiles
                        f_bar = 0.0_dp
                        do k = 1, 4
                            f_bar(profiles(k), :) = q_bar(rows(k), :)
                        end do
                        zbar = 0.0_dp
                        call observation_profiles%vjp(z, sample%q(Q_S, :), f_bar, zbar, s_bar)
                        jd(i, self%n_rec + 1:self%n_rec + self%n_prof) = gain * zbar
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
                        jd(i, self%n_rec + self%n_prof + 1:self%n_rec + self%n_prof + self%n_eta) = &
                            gain*matmul(eta_bar, self%eta_factor( &
                            slot%eta_offset - self%n_rec - self%n_prof + 1: &
                            slot%eta_offset - self%n_rec - self%n_prof + size(eta), :))
                        if (slot%gain_index > 0) jd(i, slot%gain_index) = slot%gain_sigma * y(i)
                        ! geometry and |B| (bmnc nodes) at fixed points through
                        ! the state's VJP, then the equilibrium adjoint
                        qs_bar = 0.0_dp
                        qs_bar(Q_S, :) = gain * s_bar
                        qs_bar(Q_B, :) = gain * q_bar(Q_B, :)
                        p_bar = 0.0_dp
                        call sample_state_vjp(self%state, loc, qs_bar, p_bar)
                        coef_cot(:, row + i) = p_bar(1:ng)
                        bmnc_cot(:, row + i) = p_bar(ng + 1:ng + nb)
                    end do
                    if (allocated(pb_all)) deallocate (pb_all, eb_all)
                    jd(:, slot%discrepancy_offset + 1: &
                        slot%discrepancy_offset + size(slot%discrepancy_basis, 2)) = slot%discrepancy_basis
                    jac(nrec_rows + row + 1:nrec_rows + row + size(y), :) = jd
                    row = row + size(y)
                    deallocate(y, e, jd, q_bar, eta_bar, f_bar, s_bar, qs_bar)
                end associate
            end do
            ! equilibrium adjoint for all diagnostic rows at once
            allocate(cot(ndr, rec%eq%n_y()))
            cot = 0.0_dp
            call rec%eq%jacobian(rec%eq_x, cot, jeq, coef_cot, bmnc_cot, ub_cot, vb_cot)
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
        allocate (r0(self%n_residuals()))
        ! robust_rows inverses are unnecessary: recompute untransformed block values.
        r0(:nm) = magnetic_residual(self, theta, r_rec(:nm))
        allocate (slope(nm))
        call robust_rows(r0(:nm), self%mag_student_nu, slope)
        do col = 1, self%n
            jac(:nm, col) = slope*jac(:nm, col)
        end do
        deallocate (slope)
        do d = 1, size(self%diags)
            associate (slot => self%diags(d))
                if (slot%student_nu == 0) cycle
                allocate (y(slot%d%n_signals()), slope(slot%d%n_signals()))
                call self%predict_diagnostic(theta, d, y, sample, loc, eta, dx)
                call self%unpack(theta, gain=gain, d=d)
                y = whiten(slot%chol, gain*y - slot%data + matmul(slot%discrepancy_basis, &
                    theta(slot%discrepancy_offset+1: &
                    slot%discrepancy_offset+size(slot%discrepancy_basis,2))))
                call robust_rows(y, slot%student_nu, slope)
                row = nrec_rows + slot%row_offset
                do col = 1, self%n
                    jac(row + 1:row + size(y), col) = slope*jac(row + 1:row + size(y), col)
                end do
                deallocate (y, slope)
            end associate
        end do
        ! latent priors
        row = nrec_rows + ndr
        do k = self%n_rec + 1, self%n
            jac(row + k - self%n_rec, k) = 1.0_dp
        end do
    end subroutine joint_jacobian

    subroutine ece_rows(tape, parameter_rows, calibration_rows)
        type(recorded_ece_t), intent(in) :: tape
        real(dp), intent(out) :: parameter_rows(:, :), calibration_rows(:, :)
        integer :: i
        !$omp parallel do default(shared) private(i) schedule(dynamic,1)
        do i = 1, size(tape%voltage)
            call ece_row(tape, i, parameter_rows(:, i), calibration_rows(:, i))
        end do
        !$omp end parallel do
    end subroutine

    subroutine ece_row(tape, i, parameter_bar, eta_bar)
        type(recorded_ece_t), intent(in) :: tape
        integer, intent(in) :: i
        real(dp), intent(out) :: parameter_bar(:), eta_bar(:)
        real(dp) :: voltage_bar(size(tape%voltage)), variance_bar(size(tape%voltage))
        type(fortnum_status_t) :: status
        voltage_bar = 0; variance_bar = 0; voltage_bar(i) = 1
        call tape%vjp(voltage_bar, variance_bar, parameter_bar, eta_bar, status)
        if (status%code /= 0) error stop 'joint: ECE adjoint failed'
    end subroutine

    subroutine joint_jvp(self, theta, direction, product)
        class(joint_model_t), intent(inout)::self
        real(dp), intent(in)::theta(:), direction(:)
        real(dp), intent(out)::product(:)
        real(dp), allocatable::jac(:, :)
        if (size(direction) /= self%n .or. size(product) /= self%n_residuals()) &
            error stop 'joint: tangent shape'
        call self%jacobian(theta, jac)
        product = matmul(jac, direction)
    end subroutine

    subroutine joint_vjp(self, theta, weight, bar)
        class(joint_model_t), intent(inout)::self
        real(dp), intent(in)::theta(:), weight(:)
        real(dp), intent(out)::bar(:)
        real(dp), allocatable::jac(:, :)
        if (size(weight) /= self%n_residuals() .or. size(bar) /= self%n) &
            error stop 'joint: cotangent shape'
        call self%jacobian(theta, jac)
        bar = matmul(weight, jac)
    end subroutine

    function magnetic_residual(self, theta, raw_residual) result(r)
        class(joint_model_t), intent(in)::self
        real(dp), intent(in)::theta(:), raw_residual(:)
        real(dp)::r(size(raw_residual)), c(3), v, dr, dg, db, dd, raw
        integer::i, k
        r = raw_residual*rec%meas_sigma(:self%n_mag)
        k = self%mag_calibration_index
        if (k > 0) then
            c = self%mag_calibration%physical(theta(k + 1:k + 3))
            do i = 1, self%n_mag
                raw = r(i) + rec%meas_value(i)
                call integrator_signal(raw, c(1), c(2), c(3), self%mag_time(i), &
                    self%mag_time_origin, v, dr, dg, db, dd)
                r(i) = v - rec%meas_value(i)
            end do
        end if
        r = whiten(self%mag_chol, r)
    end function

    subroutine joint_train_surrogate(self, design, length, amplitude, rho, nugget, status)
        class(joint_model_t), intent(inout)::self
        real(dp), intent(in)::design(:, :), length(:), amplitude, rho, nugget
        type(fortnum_status_t), intent(out)::status
        real(dp), allocatable::low(:, :), high(:, :), base(:), jac(:, :)
        integer::i
        logical::ok
        status%code = 1
        self%use_surrogate = .false.
        if (allocated(self%surrogate_origin)) then
            deallocate(self%surrogate_origin,self%surrogate_base,self%surrogate_jac, &
                self%surrogate_lower,self%surrogate_upper)
        end if
        status%msg = 'joint: invalid surrogate design or hyperparameters'
        if (size(design, 1) /= self%n .or. size(design, 2) < 2) return
        if (size(length) /= self%n) return
        if (.not. all(ieee_is_finite(design)) .or. .not. all(ieee_is_finite(length))) return
        if (any(length <= 0)) return
        if (.not. all(ieee_is_finite([amplitude, rho, nugget]))) return
        if (amplitude <= 0 .or. nugget < 0) return
        allocate (base(self%n_residuals()), low(self%n_residuals(), size(design, 2)), &
            high(self%n_residuals(), size(design, 2)))
        call self%residual(design(:, 1), base, ok)
        if (.not. ok) return
        call self%jacobian(design(:, 1), jac)
        do i = 1, size(design, 2)
            low(:, i) = base + matmul(jac, design(:, i) - design(:, 1))
            call self%residual(design(:, i), high(:, i), ok)
            if (.not. ok) return
        end do
        call self%surrogate%fit(design, low, high, length, amplitude, rho, nugget, status)
        if (status%code /= 0) return
        self%surrogate_origin = design(:, 1)
        self%surrogate_base = base
        self%surrogate_jac = jac
        self%surrogate_lower = minval(design, dim=2)
        self%surrogate_upper = maxval(design, dim=2)
    end subroutine

    real(dp) function joint_log_normalization(self) result(c)
        !! log p(d, theta) = -|r|^2 / 2 + c for the latent coordinates theta:
        !! Gaussian normalizations of all data rows (magnetic, consistency,
        !! profile) and of all priors (reconstruction priors in scaled units).
        class(joint_model_t), intent(in) :: self
        integer :: d, i, k
        c = -0.5_dp*self%n_residuals()*log_2pi + &
            self%n_mag*robust_normalization(self%mag_student_nu)
        do i = 1, self%n_mag
            c = c - log(self%mag_chol(i, i))
        end do
        c = c - self%n_fb * log(rec%fb_sigma)
        do i = 1, rec%n_prior
            k = rec%prior_param(i)
            c = c - log(rec%prior_sigma(k) / rec%scale(k))
        end do
        do d = 1, size(self%diags)
            c = c + self%diags(d)%d%n_signals()*robust_normalization(self%diags(d)%student_nu)
            do i = 1, size(self%diags(d)%chol, 1)
                c = c - log(self%diags(d)%chol(i, i))
            end do
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
