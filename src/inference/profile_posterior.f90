module tiago_profile_posterior
    !! Posterior of kinetic profiles and diagnostic nuisances in a fixed
    !! axisymmetric equilibrium from a G-EQDSK file (tiago_geqdsk_state). The
    !! latent vector is
    !!
    !!     theta = [ z_prof  (kinetic-profile latents, N(0, I) prior),
    !!               z_eta   (alignment of each slot, eta = eta_sigma z_eta),
    !!               z_cal   (log calibration of each slot with cal_sigma > 0) ]
    !!
    !! and the whitened residual stacks, for every slot, (c y - d) / sigma with
    !! c = exp(cal_sigma z_cal) and sigma^2 = sigma_meas^2 + (model_rel d)^2,
    !! followed by the prior rows theta. The log posterior is -|r|^2 / 2 plus
    !! log_normalization. The Jacobian is exact: diagnostic VJPs at fixed
    !! points, the motion of the points with eta and the profile VJP.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_rng, only: rng_t, rng_normal
    use tiago_local_posterior, only: local_posterior_t
    use tiago_geqdsk_plasma, only: geqdsk_plasma_t
    use tiago_geqdsk_state, only: sample_geqdsk
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t, &
        Q_S, Q_NE, Q_TE, Q_TI, Q_ZEFF, N_QUANTITIES
    use tiago_kinetic_profiles, only: kinetic_profiles_t, P_NE, P_TE, P_TI, P_ZEFF, &
        N_PROFILES
    implicit none
    private

    real(dp), parameter :: log_2pi = log(2 * acos(-1.0_dp))
    integer, parameter :: rows(4) = [Q_NE, Q_TE, Q_TI, Q_ZEFF]
    integer, parameter :: profiles(4) = [P_NE, P_TE, P_TI, P_ZEFF]

    type, public :: profile_slot_t
        class(profile_diagnostic_t), allocatable :: d
        real(dp), allocatable :: data(:), sigma(:) !! measured signals, their sd
        real(dp) :: model_rel = 0 !! relative model-error sd
        real(dp), allocatable :: eta_sigma(:) !! prior sd of the alignment
        real(dp) :: cal_sigma = 0 !! prior sd of log c; 0: none
        integer :: eta_offset = 0, cal_index = 0, row_offset = 0
    end type profile_slot_t

    type, extends(local_posterior_t), public :: profile_posterior_t
        type(geqdsk_plasma_t) :: eq
        type(kinetic_profiles_t) :: prof
        type(profile_slot_t), allocatable :: slots(:)
        integer :: n_prof = 0, n_eta = 0, n_cal = 0, n = 0, n_data = 0
    contains
        procedure :: setup => pp_setup
        procedure :: residual => pp_residual
        procedure :: jacobian => pp_jacobian
        procedure :: n_residuals => pp_n_residuals
        procedure :: log_normalization => pp_log_normalization
        procedure :: proper_prior => pp_proper_prior
        procedure :: draw_prior => pp_draw_prior
        procedure :: predict => pp_predict
        procedure :: sigma_of => pp_sigma_of
    end type profile_posterior_t

contains

    subroutine pp_setup(self)
        !! Latent and row layout, after state, prof and slots are set.
        class(profile_posterior_t), intent(inout) :: self
        integer :: k, eta, rowo
        self%n_prof = self%prof%n_latent
        eta = self%n_prof
        rowo = 0
        do k = 1, size(self%slots)
            associate (s => self%slots(k))
                if (.not. allocated(s%eta_sigma)) allocate(s%eta_sigma(s%d%n_nuisance()), &
                    source=0.0_dp)
                if (size(s%eta_sigma) /= s%d%n_nuisance()) error stop 'profile_posterior: eta_sigma'
                if (size(s%data) /= s%d%n_signals()) error stop 'profile_posterior: data size'
                s%eta_offset = eta
                eta = eta + size(s%eta_sigma)
                s%row_offset = rowo
                rowo = rowo + size(s%data)
            end associate
        end do
        self%n_eta = eta - self%n_prof
        self%n_data = rowo
        self%n_cal = 0
        do k = 1, size(self%slots)
            if (self%slots(k)%cal_sigma > 0) then
                self%n_cal = self%n_cal + 1
                self%slots(k)%cal_index = eta + self%n_cal
            else
                self%slots(k)%cal_index = 0
            end if
        end do
        self%n = self%n_prof + self%n_eta + self%n_cal
    end subroutine pp_setup

    integer function pp_n_residuals(self)
        class(profile_posterior_t), intent(in) :: self
        pp_n_residuals = self%n_data + self%n
    end function pp_n_residuals

    logical function pp_proper_prior(self)
        class(profile_posterior_t), intent(in) :: self
        pp_proper_prior = .true.
    end function pp_proper_prior

    subroutine pp_draw_prior(self, rng, theta, ok)
        class(profile_posterior_t), intent(in) :: self
        type(rng_t), intent(inout) :: rng
        real(dp), intent(out) :: theta(:)
        logical, intent(out) :: ok
        integer :: i
        do i = 1, size(theta)
            call rng_normal(rng, theta(i))
        end do
        ok = .true.
    end subroutine pp_draw_prior

    function pp_sigma_of(self, k) result(sd)
        class(profile_posterior_t), intent(in) :: self
        integer, intent(in) :: k
        real(dp), allocatable :: sd(:)
        associate (s => self%slots(k))
            sd = sqrt(s%sigma**2 + (s%model_rel * s%data)**2)
        end associate
    end function pp_sigma_of

    real(dp) function pp_log_normalization(self) result(c)
        class(profile_posterior_t), intent(in) :: self
        integer :: k
        c = -0.5_dp * self%n_residuals() * log_2pi
        do k = 1, size(self%slots)
            c = c - sum(log(self%sigma_of(k)))
        end do
    end function pp_log_normalization

    subroutine sample_at(self, z, x, sample)
        class(profile_posterior_t), intent(in) :: self
        real(dp), intent(in) :: z(:), x(:, :)
        type(plasma_sample_t), intent(out) :: sample
        real(dp), allocatable :: f(:, :), df(:, :), s(:)
        integer :: i, j
        call sample_geqdsk(self%eq, x, sample)
        s = sample%q(Q_S, :)
        allocate(f(N_PROFILES, size(s)), df(N_PROFILES, size(s)))
        call self%prof%evaluate(z, s, f, df)
        do j = 1, size(s)
            do i = 1, 4
                sample%q(rows(i), j) = f(profiles(i), j)
                sample%dq_dx(:, rows(i), j) = df(profiles(i), j) * sample%dq_dx(:, Q_S, j)
            end do
        end do
    end subroutine sample_at

    subroutine pp_predict(self, theta, k, y, sample, eta, dx_deta, c)
        !! Uncalibrated prediction y of slot k, its sample, alignment and the
        !! calibration factor c.
        class(profile_posterior_t), intent(in) :: self
        real(dp), intent(in) :: theta(:)
        integer, intent(in) :: k
        real(dp), allocatable, intent(out) :: y(:), eta(:), dx_deta(:, :, :)
        type(plasma_sample_t), intent(out) :: sample
        real(dp), intent(out) :: c
        real(dp), allocatable :: x(:, :)
        associate (s => self%slots(k))
            eta = s%eta_sigma * theta(s%eta_offset + 1:s%eta_offset + size(s%eta_sigma))
            call s%d%sample_points(eta, x, dx_deta)
            call sample_at(self, theta(1:self%n_prof), x, sample)
            allocate(y(size(s%data)))
            call s%d%predict(sample, eta, y)
            c = 1
            if (s%cal_index > 0) c = exp(s%cal_sigma * theta(s%cal_index))
        end associate
    end subroutine pp_predict

    subroutine pp_residual(self, theta, r, ok)
        class(profile_posterior_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), intent(out) :: r(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: y(:), eta(:), dx(:, :, :)
        type(plasma_sample_t) :: sample
        real(dp) :: c
        integer :: k
        do k = 1, size(self%slots)
            call self%predict(theta, k, y, sample, eta, dx, c)
            associate (s => self%slots(k))
                r(s%row_offset + 1:s%row_offset + size(y)) = (c * y - s%data) / self%sigma_of(k)
            end associate
        end do
        r(self%n_data + 1:) = theta
        ok = all(abs(r) < huge(1.0_dp))
    end subroutine pp_residual

    subroutine pp_jacobian(self, theta, jac)
        class(profile_posterior_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), allocatable, intent(out) :: jac(:, :)
        real(dp), allocatable :: y(:), eta(:), dx(:, :, :), sd(:), y_bar(:), q_bar(:, :), &
            eta_bar(:), f_bar(:, :), s_bar(:), z_bar(:)
        type(plasma_sample_t) :: sample
        real(dp) :: c
        integer :: k, i, j, q, row
        allocate(jac(self%n_residuals(), self%n), source=0.0_dp)
        do k = 1, size(self%slots)
            call self%predict(theta, k, y, sample, eta, dx, c)
            sd = self%sigma_of(k)
            associate (s => self%slots(k))
                allocate(y_bar(size(y)), q_bar(N_QUANTITIES, sample%n), eta_bar(size(eta)), &
                    f_bar(N_PROFILES, sample%n), s_bar(sample%n), z_bar(self%n_prof))
                do i = 1, size(y)
                    row = s%row_offset + i
                    y_bar = 0
                    y_bar(i) = c / sd(i)
                    call s%d%vjp(sample, eta, y_bar, q_bar, eta_bar)
                    ! points move with eta: q changes by dq/dx dx/deta
                    do j = 1, sample%n
                        do q = 1, N_QUANTITIES
                            eta_bar = eta_bar + q_bar(q, j) * matmul(sample%dq_dx(:, q, j), dx(:, j, :))
                        end do
                    end do
                    jac(row, s%eta_offset + 1:s%eta_offset + size(eta)) = eta_bar * s%eta_sigma
                    f_bar = 0
                    do q = 1, 4
                        f_bar(profiles(q), :) = q_bar(rows(q), :)
                    end do
                    z_bar = 0
                    call self%prof%vjp(theta(1:self%n_prof), sample%q(Q_S, :), f_bar, z_bar, s_bar)
                    jac(row, 1:self%n_prof) = z_bar
                    if (s%cal_index > 0) jac(row, s%cal_index) = s%cal_sigma * c * y(i) / sd(i)
                end do
                deallocate(y_bar, q_bar, eta_bar, f_bar, s_bar, z_bar)
            end associate
        end do
        do i = 1, self%n
            jac(self%n_data + i, i) = 1
        end do
    end subroutine pp_jacobian

end module tiago_profile_posterior
