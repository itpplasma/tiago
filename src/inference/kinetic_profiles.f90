module tiago_kinetic_profiles
    !! Kinetic profiles ne, Te, Ti, Zeff and the fast-particle pressure p_fast
    !! as functions of the provider's flux label s, the derived pressure
    !!
    !!     p = e (ne Te + n_i Ti) + p_fast,
    !!     n_i = n_D + n_Z = ne (Z (Z - Zeff) + Zeff - 1) / (Z (Z - 1))
    !!
    !! (deuterium and one impurity of charge Z = z_imp), and their reverse-mode
    !! derivatives.
    !!
    !! Each profile is u(s) = sum_i c_i B_i(s), a clamped cubic B-spline on
    !! [0, 1], mapped to the physical value by f = unit exp(u) (positive
    !! profiles) or Zeff = 1 + (z_imp - 1) / (1 + exp(-u)). The coefficients
    !! have a Gaussian-process prior evaluated at the Greville abscissae,
    !! c = m + L z with z ~ N(0, I) the latent variables and L L^T the kernel
    !! matrix of the Minerva/surroBIER non-stationary Gibbs kernel
    !!
    !!     k(s, s') = sigma_f^2 sqrt(2 l l' / (l^2 + l'^2)) exp(-(s - s')^2 / (l^2 + l'^2)),
    !!     l(s) = (l_core + l_edge) / 2 - (l_core - l_edge) / 2 tanh((s - s0) / sw).
    !!
    !! (kernel = KERNEL_MINERVA), or its Matern-nu analogue (Paciorek and
    !! Schervish 2004) with the same prefactor and length scale and
    !! r = |s - s'| / sqrt((l^2 + l'^2) / 2):
    !!     nu = 3/2: (1 + sqrt(3) r) exp(-sqrt(3) r),
    !!     nu = 5/2: (1 + sqrt(5) r + 5 r^2 / 3) exp(-sqrt(5) r).
    !! The squared-exponential kernel makes the prior effectively
    !! low-dimensional: a local flattening costs |z|^2 ~ 1e2 under it, while the
    !! Matern kernels (default 5/2) keep such profiles in the prior bulk.
    !!
    !! The prior acts on the log value, not on the log gradient as in
    !! surroBIER's PCA prior, so flat regions (plateaus from islands) have
    !! non-zero prior density. A profile that is not latent keeps fixed
    !! coefficients c = m.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_bspline, only: bspline_workspace_t, bspline_init, bspline_set_knots, &
        bspline_eval_deriv
    use tiago_kinetic_pressure, only: kinetic_pressure
    implicit none
    private

    integer, parameter, public :: P_NE = 1, P_TE = 2, P_TI = 3, P_ZEFF = 4, P_FAST = 5
    integer, parameter, public :: N_PROFILES = 5
    integer, parameter, public :: KERNEL_MINERVA = 0, KERNEL_MATERN32 = 1, KERNEL_MATERN52 = 2
    real(dp), parameter, public :: elementary_charge = 1.602176634e-19_dp
    character(len=5), parameter, public :: profile_names(N_PROFILES) = &
        ['ne   ', 'Te   ', 'Ti   ', 'Zeff ', 'pfast']

    type, public :: profile_t
        logical :: latent = .false.
        real(dp) :: unit = 1.0_dp                   !! physical value for u = 0 (log profiles)
        real(dp), allocatable :: mean(:)            !! prior mean of the coefficients
        real(dp), allocatable :: chol(:, :)         !! lower Cholesky factor of the prior covariance
        integer :: offset = 0                       !! first latent index - 1
    end type profile_t

    type, public :: kinetic_profiles_t
        type(bspline_workspace_t) :: basis
        integer :: ncoef = 0
        real(dp), allocatable :: greville(:)
        real(dp) :: z_imp = 6.0_dp
        type(profile_t) :: p(N_PROFILES)
        integer :: n_latent = 0
    contains
        procedure :: init => profiles_init
        procedure :: set_profile => profiles_set
        procedure :: coefficients => profiles_coefficients
        procedure :: evaluate => profiles_evaluate
        procedure :: vjp => profiles_vjp
        procedure :: pressure => profiles_pressure
        procedure :: pressure_vjp => profiles_pressure_vjp
        procedure :: fit_latent => profiles_fit_latent
    end type kinetic_profiles_t

    public :: gibbs_kernel

contains

    subroutine profiles_init(self, breakpoints, z_imp)
        !! Cubic B-spline basis on the given breakpoints in s (0 ... 1).
        class(kinetic_profiles_t), intent(inout) :: self
        real(dp), intent(in) :: breakpoints(:)
        real(dp), intent(in), optional :: z_imp
        type(fortnum_status_t) :: st
        integer :: i, k

        if (present(z_imp)) self%z_imp = z_imp
        call bspline_init(self%basis, 4, size(breakpoints), st)
        if (st%code /= FORTNUM_OK) error stop 'kinetic profiles: '//trim(st%msg)
        call bspline_set_knots(self%basis, breakpoints, st)
        if (st%code /= FORTNUM_OK) error stop 'kinetic profiles: '//trim(st%msg)
        self%ncoef = self%basis%ncoef
        allocate(self%greville(self%ncoef))
        k = self%basis%order
        do i = 1, self%ncoef
            self%greville(i) = sum(self%basis%knots(i + 1:i + k - 1)) / (k - 1)
        end do
        do i = 1, N_PROFILES
            self%p(i)%latent = .false.
            allocate(self%p(i)%mean(self%ncoef))
            self%p(i)%mean = 0.0_dp
        end do
        self%p(P_NE)%unit = 1.0e19_dp
        self%p(P_TE)%unit = 1.0e3_dp
        self%p(P_TI)%unit = 1.0e3_dp
        self%p(P_FAST)%unit = 1.0e3_dp
        self%p(P_FAST)%mean = -50.0_dp      ! p_fast ~ 0 unless set
        self%p(P_ZEFF)%mean = -50.0_dp      ! Zeff ~ 1 unless set
        call layout(self)
    end subroutine profiles_init

    subroutine profiles_set(self, which, mean_values, latent, sigma_f, l_core, l_edge, s0, sw, kernel)
        !! Prior of one profile: mean_values are physical values at the
        !! Greville abscissae (the prior mean of u is their transform); the
        !! Gibbs kernel parameters apply to u. Without latent = .true. the
        !! profile is fixed at mean_values.
        class(kinetic_profiles_t), intent(inout) :: self
        integer, intent(in) :: which
        real(dp), intent(in) :: mean_values(:)
        logical, intent(in) :: latent
        real(dp), intent(in), optional :: sigma_f, l_core, l_edge, s0, sw
        integer, intent(in), optional :: kernel
        real(dp), allocatable :: k(:, :)
        integer :: i, j, kind

        if (size(mean_values) /= self%ncoef) error stop 'set_profile: one value per coefficient'
        do i = 1, self%ncoef
            self%p(which)%mean(i) = inverse_transform(self, which, mean_values(i))
        end do
        self%p(which)%latent = latent
        if (latent) then
            if (.not. (present(sigma_f) .and. present(l_core) .and. present(l_edge) .and. &
                present(s0) .and. present(sw))) error stop 'set_profile: kernel parameters required'
            kind = KERNEL_MATERN52
            if (present(kernel)) kind = kernel
            allocate(k(self%ncoef, self%ncoef))
            do j = 1, self%ncoef
                do i = 1, self%ncoef
                    k(i, j) = gibbs_kernel(self%greville(i), self%greville(j), sigma_f, l_core, &
                        l_edge, s0, sw, kind)
                end do
                k(j, j) = k(j, j) * (1.0_dp + 1.0e-8_dp)
            end do
            self%p(which)%chol = cholesky(k)
        end if
        call layout(self)
    end subroutine profiles_set

    subroutine layout(self)
        type(kinetic_profiles_t), intent(inout) :: self
        integer :: i
        self%n_latent = 0
        do i = 1, N_PROFILES
            self%p(i)%offset = self%n_latent
            if (self%p(i)%latent) self%n_latent = self%n_latent + self%ncoef
        end do
    end subroutine layout

    pure real(dp) function gibbs_kernel(a, b, sigma_f, l_core, l_edge, s0, sw, kind) result(k)
        real(dp), intent(in) :: a, b, sigma_f, l_core, l_edge, s0, sw
        integer, intent(in) :: kind
        real(dp) :: la, lb, r
        la = 0.5_dp * (l_core + l_edge) - 0.5_dp * (l_core - l_edge) * tanh((a - s0) / sw)
        lb = 0.5_dp * (l_core + l_edge) - 0.5_dp * (l_core - l_edge) * tanh((b - s0) / sw)
        k = sigma_f**2 * sqrt(2.0_dp * la * lb / (la**2 + lb**2))
        r = abs(a - b) / sqrt(0.5_dp * (la**2 + lb**2))
        select case (kind)
        case (KERNEL_MINERVA)
            k = k * exp(-0.5_dp * r**2)
        case (KERNEL_MATERN32)
            k = k * (1.0_dp + sqrt(3.0_dp) * r) * exp(-sqrt(3.0_dp) * r)
        case default
            k = k * (1.0_dp + sqrt(5.0_dp) * r + 5.0_dp * r**2 / 3.0_dp) * exp(-sqrt(5.0_dp) * r)
        end select
    end function gibbs_kernel

    function profiles_coefficients(self, z, which) result(c)
        class(kinetic_profiles_t), intent(in) :: self
        real(dp), intent(in) :: z(:)
        integer, intent(in) :: which
        real(dp) :: c(self%ncoef)
        associate (p => self%p(which))
            c = p%mean
            if (p%latent) c = c + matmul(p%chol, z(p%offset + 1:p%offset + self%ncoef))
        end associate
    end function profiles_coefficients

    subroutine profiles_evaluate(self, z, s, f, df_ds)
        !! Profile values f(N_PROFILES, n) and d f / d s at s(n) (clamped to
        !! [0, 1]; the caller decides what applies outside the plasma).
        class(kinetic_profiles_t), intent(in) :: self
        real(dp), intent(in) :: z(:), s(:)
        real(dp), intent(out) :: f(:, :), df_ds(:, :)
        real(dp) :: c(self%ncoef, N_PROFILES), b(0:1, self%ncoef), u, du, g, dg
        integer :: i, j

        do i = 1, N_PROFILES
            c(:, i) = self%coefficients(z, i)
        end do
        !$omp parallel do private(b, i, u, du, g, dg) if (size(s) > 256)
        do j = 1, size(s)
            call basis_at(self, s(j), b)
            do i = 1, N_PROFILES
                u = dot_product(b(0, :), c(:, i))
                du = dot_product(b(1, :), c(:, i))
                call transform(self, i, u, g, dg)
                f(i, j) = g
                df_ds(i, j) = dg * du
            end do
        end do
        !$omp end parallel do
    end subroutine profiles_evaluate

    subroutine profiles_vjp(self, z, s, f_bar, z_bar, s_bar)
        !! z_bar += sum f_bar df/dz; s_bar = sum_i f_bar(i, j) df_i/ds (j).
        class(kinetic_profiles_t), intent(in) :: self
        real(dp), intent(in) :: z(:), s(:), f_bar(:, :)
        real(dp), intent(inout) :: z_bar(:)
        real(dp), intent(out) :: s_bar(:)
        real(dp) :: c(self%ncoef, N_PROFILES), c_bar(self%ncoef, N_PROFILES)
        real(dp) :: b(0:1, self%ncoef), u, du, g, dg
        integer :: i, j

        do i = 1, N_PROFILES
            c(:, i) = self%coefficients(z, i)
        end do
        c_bar = 0.0_dp
        !$omp parallel do private(b, i, u, du, g, dg) reduction(+:c_bar) if (size(s) > 256)
        do j = 1, size(s)
            call basis_at(self, s(j), b)
            s_bar(j) = 0.0_dp
            do i = 1, N_PROFILES
                u = dot_product(b(0, :), c(:, i))
                du = dot_product(b(1, :), c(:, i))
                call transform(self, i, u, g, dg)
                c_bar(:, i) = c_bar(:, i) + f_bar(i, j) * dg * b(0, :)
                s_bar(j) = s_bar(j) + f_bar(i, j) * dg * du
            end do
        end do
        !$omp end parallel do
        do i = 1, N_PROFILES
            associate (p => self%p(i))
                if (p%latent) z_bar(p%offset + 1:p%offset + self%ncoef) = &
                    z_bar(p%offset + 1:p%offset + self%ncoef) + matmul(c_bar(:, i), p%chol)
            end associate
        end do
    end subroutine profiles_vjp

    subroutine profiles_pressure(self, z, s, p, dp_ds)
        !! Pressure [Pa] at s(n) and its s derivative.
        class(kinetic_profiles_t), intent(in) :: self
        real(dp), intent(in) :: z(:), s(:)
        real(dp), intent(out) :: p(:)
        real(dp), intent(out), optional :: dp_ds(:)
        real(dp), allocatable :: f(:, :), df(:, :)
        real(dp) :: g(N_PROFILES)
        integer :: j

        allocate(f(N_PROFILES, size(s)), df(N_PROFILES, size(s)))
        call self%evaluate(z, s, f, df)
        do j = 1, size(s)
            call kinetic_pressure(f(P_NE, j), f(P_TE, j), f(P_TI, j), f(P_ZEFF, j), f(P_FAST, j), &
                self%z_imp, elementary_charge, p(j), g(P_NE), g(P_TE), g(P_TI), g(P_ZEFF), g(P_FAST))
            if (present(dp_ds)) dp_ds(j) = dot_product(g, df(:, j))
        end do
    end subroutine profiles_pressure

    subroutine profiles_pressure_vjp(self, z, s, p_bar, z_bar)
        !! z_bar += (dp/dz)^T p_bar at fixed s.
        class(kinetic_profiles_t), intent(in) :: self
        real(dp), intent(in) :: z(:), s(:), p_bar(:)
        real(dp), intent(inout) :: z_bar(:)
        real(dp), allocatable :: f(:, :), df(:, :), f_bar(:, :), s_bar(:)
        real(dp) :: p, g(N_PROFILES)
        integer :: j

        allocate(f(N_PROFILES, size(s)), df(N_PROFILES, size(s)), f_bar(N_PROFILES, size(s)), &
            s_bar(size(s)))
        call self%evaluate(z, s, f, df)
        do j = 1, size(s)
            call kinetic_pressure(f(P_NE, j), f(P_TE, j), f(P_TI, j), f(P_ZEFF, j), f(P_FAST, j), &
                self%z_imp, elementary_charge, p, g(P_NE), g(P_TE), g(P_TI), g(P_ZEFF), g(P_FAST))
            f_bar(:, j) = p_bar(j) * g
        end do
        call self%vjp(z, s, f_bar, z_bar, s_bar)
    end subroutine profiles_pressure_vjp

    subroutine profiles_fit_latent(self, which, s, values, z)
        !! Latent z of one profile whose spline interpolates the physical
        !! values at the Greville abscissae (s, values given there); for
        !! synthetic truths built from analytic profiles.
        class(kinetic_profiles_t), intent(in) :: self
        integer, intent(in) :: which
        real(dp), intent(in) :: s(:), values(:)
        real(dp), intent(inout) :: z(:)
        real(dp), allocatable :: a(:, :), u(:), c(:)
        real(dp) :: b(0:1, self%ncoef)
        integer :: j, info, ipiv(self%ncoef)
        external :: dgesv

        if (size(s) /= self%ncoef) error stop 'fit_latent: one value per coefficient'
        allocate(a(self%ncoef, self%ncoef), u(self%ncoef))
        do j = 1, self%ncoef
            call basis_at(self, s(j), b)
            a(j, :) = b(0, :)
            u(j) = inverse_transform(self, which, values(j))
        end do
        call dgesv(self%ncoef, 1, a, self%ncoef, ipiv, u, self%ncoef, info)
        if (info /= 0) error stop 'fit_latent: singular collocation'
        associate (p => self%p(which))
            if (.not. p%latent) error stop 'fit_latent: profile is not latent'
            c = u - p%mean
            ! z = L^-1 c (forward substitution)
            do j = 1, self%ncoef
                z(p%offset + j) = (c(j) - dot_product(p%chol(j, 1:j - 1), &
                    z(p%offset + 1:p%offset + j - 1))) / p%chol(j, j)
            end do
        end associate
    end subroutine profiles_fit_latent

    subroutine basis_at(self, s, b)
        type(kinetic_profiles_t), intent(in) :: self
        real(dp), intent(in) :: s
        real(dp), intent(out) :: b(0:, :)
        type(fortnum_status_t) :: st
        call bspline_eval_deriv(self%basis, min(max(s, 0.0_dp), 1.0_dp), 1, b, st)
        if (s > 1.0_dp .or. s < 0.0_dp) b(1, :) = 0.0_dp
    end subroutine basis_at

    subroutine transform(self, which, u, g, dg)
        type(kinetic_profiles_t), intent(in) :: self
        integer, intent(in) :: which
        real(dp), intent(in) :: u
        real(dp), intent(out) :: g, dg
        real(dp) :: sig
        if (which == P_ZEFF) then
            sig = 1.0_dp / (1.0_dp + exp(-u))
            g = 1.0_dp + (self%z_imp - 1.0_dp) * sig
            dg = (self%z_imp - 1.0_dp) * sig * (1.0_dp - sig)
        else
            g = self%p(which)%unit * exp(u)
            dg = g
        end if
    end subroutine transform

    real(dp) function inverse_transform(self, which, value) result(u)
        type(kinetic_profiles_t), intent(in) :: self
        integer, intent(in) :: which
        real(dp), intent(in) :: value
        real(dp) :: t
        if (which == P_ZEFF) then
            t = (value - 1.0_dp) / (self%z_imp - 1.0_dp)
            t = min(max(t, 1.0e-12_dp), 1.0_dp - 1.0e-12_dp)
            u = log(t / (1.0_dp - t))
        else
            u = log(max(value, tiny(1.0_dp)) / self%p(which)%unit)
        end if
    end function inverse_transform

    function cholesky(a) result(l)
        real(dp), intent(in) :: a(:, :)
        real(dp) :: l(size(a, 1), size(a, 1))
        integer :: info, j
        external :: dpotrf
        l = a
        call dpotrf('L', size(a, 1), l, size(a, 1), info)
        if (info /= 0) error stop 'kinetic profiles: prior covariance not positive definite'
        do j = 2, size(a, 1)
            l(1:j - 1, j) = 0.0_dp
        end do
    end function cholesky

end module tiago_kinetic_profiles
