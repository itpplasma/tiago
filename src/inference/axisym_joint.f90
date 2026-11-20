module tiago_axisym_joint
    !! Joint posterior of an axisymmetric free-boundary equilibrium and the
    !! kinetic profiles (docs/axisymmetric-free-boundary.md). The latent vector
    !! is theta = [theta_fb, theta_prof], the parameters of
    !! tiago_axisym_free_boundary followed by the latents of
    !! tiago_profile_posterior. The residual stacks
    !!
    !!     free-boundary rows of the solved equilibrium (sensors, MSE, profile
    !!         rows, Gaussian priors),
    !!     profile-diagnostic rows sampled in the same equilibrium, with their
    !!         alignment, calibration and latent prior rows,
    !!     pressure rows (p_eq(psi_N) - p_kin(s(psi_N))) / sigma_p,
    !!
    !! where p_kin = e (n_e T_e + n_i T_i) of the profiles at the toroidal-flux
    !! label s of the solved equilibrium. Free-boundary columns of the Jacobian
    !! are central differences of re-solved equilibria, with the step control
    !! of tiago_axisym_free_boundary. Profile columns are exact at the
    !! equilibrium of theta: the diagnostic VJPs of tiago_profile_posterior and
    !! the pressure VJP of the profiles.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_rng, only: rng_t
    use tiago_local_posterior, only: local_posterior_t
    use tiago_geqdsk_plasma, only: geqdsk_plasma_t
    use tiago_axisym_free_boundary, only: axisym_free_boundary_t
    use tiago_profile_posterior, only: profile_posterior_t
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    type, extends(local_posterior_t), public :: axisym_joint_t
        type(axisym_free_boundary_t) :: fb
        type(profile_posterior_t) :: pp
        !> Labels psi_N of the pressure rows and their sd [Pa].
        real(dp), allocatable :: p_psin(:), p_sigma(:)
    contains
        procedure :: residual => joint_residual
        procedure :: jacobian => joint_jacobian
        procedure :: n_residuals => joint_n_residuals
        procedure :: log_normalization => joint_log_normalization
        procedure :: proper_prior => joint_proper_prior
        procedure :: draw_prior => joint_draw_prior
        procedure :: n_fb => joint_n_fb
        procedure :: kinetic_pressure => joint_kinetic_pressure
    end type axisym_joint_t

contains

    integer function joint_n_fb(self)
        class(axisym_joint_t), intent(in) :: self
        joint_n_fb = size(self%fb%names)
    end function joint_n_fb

    integer function joint_n_residuals(self)
        class(axisym_joint_t), intent(in) :: self
        joint_n_residuals = self%fb%n_residuals() + self%pp%n_residuals() + size(self%p_psin)
    end function joint_n_residuals

    subroutine joint_kinetic_pressure(self, z, eq, s, p)
        !! Toroidal-flux labels s of p_psin in eq and the kinetic pressure there.
        class(axisym_joint_t), intent(in) :: self
        real(dp), intent(in) :: z(:)
        type(geqdsk_plasma_t), intent(in) :: eq
        real(dp), allocatable, intent(out) :: s(:), p(:)
        real(dp) :: ds
        integer :: k
        allocate(s(size(self%p_psin)), p(size(self%p_psin)))
        do k = 1, size(s)
            call eq%toroidal_flux(self%p_psin(k), s(k), ds)
        end do
        call self%pp%prof%pressure(z(:self%pp%n_prof), s, p)
    end subroutine joint_kinetic_pressure

    subroutine joint_residual(self, theta, r, ok)
        class(axisym_joint_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), intent(out) :: r(:)
        logical, intent(out) :: ok
        type(geqdsk_plasma_t) :: eq
        real(dp), allocatable :: sig(:), s(:), p(:)
        integer :: nf, nr, np, k
        nf = self%n_fb()
        nr = self%fb%n_residuals()
        np = self%pp%n_residuals()
        call self%fb%signals(theta(:nf), sig, ok, eq)
        if (.not. ok) return
        call self%fb%rows(theta(:nf), sig, r(:nr))
        self%pp%eq = eq
        call self%pp%residual(theta(nf + 1:), r(nr + 1:nr + np), ok)
        if (.not. ok) return
        call self%kinetic_pressure(theta(nf + 1:), eq, s, p)
        do k = 1, size(p)
            r(nr + np + k) = (eq%pressure(self%p_psin(k)) - p(k)) / self%p_sigma(k)
        end do
    end subroutine joint_residual

    subroutine joint_jacobian(self, theta, jac)
        class(axisym_joint_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), allocatable, intent(out) :: jac(:, :)
        type(geqdsk_plasma_t) :: eq0
        real(dp), allocatable :: r0(:), rp(:), rm(:), t(:), jpp(:, :), s(:), p(:), z_bar(:)
        real(dp) :: h
        logical :: ok, okp, okm
        integer :: nf, nr, np, n, k, attempt
        nf = self%n_fb()
        nr = self%fb%n_residuals()
        np = self%pp%n_residuals()
        n = self%n_residuals()
        allocate(jac(n, size(theta)), r0(n), rp(n), rm(n))
        jac = 0
        call self%residual(theta, r0, ok)
        if (.not. ok) error stop 'axisym joint: no equilibrium at the point'
        eq0 = self%pp%eq
        ! Profile columns at the equilibrium of theta.
        call self%pp%jacobian(theta(nf + 1:), jpp)
        jac(nr + 1:nr + np, nf + 1:) = jpp
        call self%kinetic_pressure(theta(nf + 1:), eq0, s, p)
        allocate(z_bar(self%pp%n_prof))
        do k = 1, size(s)
            z_bar = 0
            call self%pp%prof%pressure_vjp(theta(nf + 1:nf + self%pp%n_prof), s(k:k), [1.0_dp], z_bar)
            jac(nr + np + k, nf + 1:nf + self%pp%n_prof) = -z_bar / self%p_sigma(k)
        end do
        ! Free-boundary columns by re-solved central differences. Where a side
        ! does not converge, the step shrinks twice by 4, then one side is used.
        do k = 1, nf
            h = self%fb%fd_relative * max(abs(theta(k)), 1.0_dp)
            do attempt = 1, 3
                t = theta
                t(k) = theta(k) + h
                call self%residual(t, rp, okp)
                t(k) = theta(k) - h
                call self%residual(t, rm, okm)
                if (okp .and. okm) exit
                h = h / 4
            end do
            if (okp .and. okm) then
                jac(:, k) = (rp - rm) / (2 * h)
            else if (okp) then
                jac(:, k) = (rp - r0) / h
            else if (okm) then
                jac(:, k) = (r0 - rm) / h
            else
                write(error_unit, '(A)') 'axisym joint: no derivative for '//trim(self%fb%names(k))
                error stop 'axisym joint: derivative solve failed'
            end if
        end do
        self%pp%eq = eq0
    end subroutine joint_jacobian

    real(dp) function joint_log_normalization(self) result(v)
        class(axisym_joint_t), intent(in) :: self
        v = self%fb%log_normalization() + self%pp%log_normalization() &
            - sum(log(sqrt(2 * pi) * self%p_sigma))
    end function joint_log_normalization

    logical function joint_proper_prior(self)
        class(axisym_joint_t), intent(in) :: self
        joint_proper_prior = self%fb%proper_prior()
    end function joint_proper_prior

    subroutine joint_draw_prior(self, rng, theta, ok)
        class(axisym_joint_t), intent(in) :: self
        type(rng_t), intent(inout) :: rng
        real(dp), intent(out) :: theta(:)
        logical, intent(out) :: ok
        integer :: nf
        nf = self%n_fb()
        call self%fb%draw_prior(rng, theta(:nf), ok)
        if (ok) call self%pp%draw_prior(rng, theta(nf + 1:), ok)
    end subroutine joint_draw_prior

end module tiago_axisym_joint
