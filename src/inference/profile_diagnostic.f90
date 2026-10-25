module tiago_profile_diagnostic
    !! Contract between TIAGO's inference core and profile diagnostics
    !! (Thomson scattering, ECE, interferometry, XICS, ...). See
    !! docs/inference.md, section "Profile-diagnostic interface".
    !!
    !! A diagnostic declares real-space sample points x(eta), which may move
    !! with its own alignment nuisance parameters eta. The core evaluates the
    !! provider fields and the kinetic profiles there (plasma_sample_t) and the
    !! diagnostic maps them to predicted signals y(sample, eta). The diagnostic
    !! supplies the forward map and its tangent (jvp) and adjoint (vjp) with
    !! respect to the sampled quantities q and to eta at fixed points; the core
    !! adds the motion of the points (dpoints/deta with dq/dx) and propagates
    !! q_bar to the latent parameters through the profiles and the
    !! equilibrium provider (flux label s and |B| through its implicit adjoint).
    !!
    !! Diagnostics never see flux coordinates of the provider, its state layout
    !! or the latent parameterization. Calibration gains/offsets common to all
    !! diagnostics and the measurement covariance belong to the likelihood
    !! (tiago_likelihood), not to the diagnostic.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    !! Rows of plasma_sample_t%q.
    integer, parameter, public :: Q_S = 1      !! flux label s (normalized toroidal flux), > 1 outside
    integer, parameter, public :: Q_B = 2      !! |B| [T]
    integer, parameter, public :: Q_NE = 3     !! electron density [m^-3]
    integer, parameter, public :: Q_TE = 4     !! electron temperature [eV]
    integer, parameter, public :: Q_TI = 5     !! ion temperature [eV]
    integer, parameter, public :: Q_ZEFF = 6   !! effective charge
    integer, parameter, public :: N_QUANTITIES = 6

    type, public :: plasma_sample_t
        !! Provider fields and kinetic profiles at n real-space points.
        integer :: n = 0
        real(dp), allocatable :: x(:, :)          !! (3, n): R [m], phi [rad], Z [m]
        logical, allocatable :: inside(:)         !! inside the last closed flux surface
        real(dp), allocatable :: q(:, :)          !! (N_QUANTITIES, n)
        real(dp), allocatable :: dq_dx(:, :, :)   !! (3, N_QUANTITIES, n): d q / d (R, phi, Z)
    end type plasma_sample_t
    !! Outside the plasma (inside = .false.) ne = Te = Ti = 0, Zeff = 1, s is
    !! the provider's extrapolation and |B| remains defined; the profile rows of
    !! dq_dx are zero there.

    type, abstract, public :: profile_diagnostic_t
    contains
        procedure(count_i), deferred :: n_signals
        procedure(count_i), deferred :: n_nuisance
        procedure(label_i), deferred :: label
        procedure(points_i), deferred :: sample_points
        procedure(predict_i), deferred :: predict
        procedure(jvp_i), deferred :: jvp
        procedure(vjp_i), deferred :: vjp
    end type profile_diagnostic_t

    abstract interface
        integer function count_i(self)
            import :: profile_diagnostic_t
            class(profile_diagnostic_t), intent(in) :: self
        end function count_i

        function label_i(self, i) result(text)
            !! Name of signal i (for reports and plots).
            import :: profile_diagnostic_t
            class(profile_diagnostic_t), intent(in) :: self
            integer, intent(in) :: i
            character(len=:), allocatable :: text
        end function label_i

        subroutine points_i(self, eta, x, dx_deta)
            !! Sample points at nuisance eta(n_nuisance) and their derivative.
            import :: profile_diagnostic_t, dp
            class(profile_diagnostic_t), intent(in) :: self
            real(dp), intent(in) :: eta(:)
            real(dp), allocatable, intent(out) :: x(:, :)            !! (3, n_points)
            real(dp), allocatable, intent(out) :: dx_deta(:, :, :)   !! (3, n_points, n_nuisance)
        end subroutine points_i

        subroutine predict_i(self, sample, eta, y)
            !! Predicted signals y(n_signals) at the sampled quantities.
            import :: profile_diagnostic_t, plasma_sample_t, dp
            class(profile_diagnostic_t), intent(in) :: self
            type(plasma_sample_t), intent(in) :: sample
            real(dp), intent(in) :: eta(:)
            real(dp), intent(out) :: y(:)
        end subroutine predict_i

        subroutine jvp_i(self, sample, eta, q_dot, eta_dot, y_dot)
            !! y_dot = dy/dq q_dot + dy/deta|_x eta_dot (points held fixed).
            import :: profile_diagnostic_t, plasma_sample_t, dp
            class(profile_diagnostic_t), intent(in) :: self
            type(plasma_sample_t), intent(in) :: sample
            real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
            real(dp), intent(out) :: y_dot(:)
        end subroutine jvp_i

        subroutine vjp_i(self, sample, eta, y_bar, q_bar, eta_bar)
            !! q_bar = (dy/dq)^T y_bar, eta_bar = (dy/deta|_x)^T y_bar.
            import :: profile_diagnostic_t, plasma_sample_t, dp
            class(profile_diagnostic_t), intent(in) :: self
            type(plasma_sample_t), intent(in) :: sample
            real(dp), intent(in) :: eta(:), y_bar(:)
            real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        end subroutine vjp_i
    end interface

    public :: check_diagnostic_derivatives

contains

    subroutine check_diagnostic_derivatives(diag, sample, eta, h, jvp_error, adjoint_error)
        !! Test helper for diagnostic authors. jvp_error: maximum relative
        !! difference between jvp and central differences of predict along a
        !! fixed pseudo-random direction (q, eta), step h relative to |q|;
        !! adjoint_error: relative mismatch of <y_bar, J v> and <J^T y_bar, v>.
        !! Inside points perturb all rows, outside points only s and |B|.
        class(profile_diagnostic_t), intent(in) :: diag
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), h
        real(dp), intent(out) :: jvp_error, adjoint_error
        type(plasma_sample_t) :: plus, minus
        real(dp), allocatable :: q_dot(:, :), eta_dot(:), y_dot(:), yp(:), ym(:), fd(:)
        real(dp), allocatable :: y_bar(:), q_bar(:, :), eta_bar(:)
        real(dp) :: lhs, rhs, scale
        integer :: i, j, m

        m = diag%n_signals()
        allocate(q_dot(N_QUANTITIES, sample%n), eta_dot(size(eta)), y_dot(m), yp(m), ym(m), &
            y_bar(m), q_bar(N_QUANTITIES, sample%n), eta_bar(size(eta)))
        do j = 1, sample%n
            do i = 1, N_QUANTITIES
                q_dot(i, j) = sin(1.7_dp * i + 0.31_dp * j) * max(abs(sample%q(i, j)), 1.0e-3_dp)
                if (.not. sample%inside(j) .and. i > Q_B) q_dot(i, j) = 0.0_dp
            end do
        end do
        eta_dot = [(cos(0.9_dp * i) * max(abs(eta(i)), 1.0e-3_dp), i = 1, size(eta))]
        y_bar = [(cos(0.37_dp * i + 0.2_dp), i = 1, m)]
        call diag%jvp(sample, eta, q_dot, eta_dot, y_dot)
        plus = sample
        minus = sample
        plus%q = sample%q + h * q_dot
        minus%q = sample%q - h * q_dot
        call diag%predict(plus, eta + h * eta_dot, yp)
        call diag%predict(minus, eta - h * eta_dot, ym)
        fd = (yp - ym) / (2.0_dp * h)
        scale = max(maxval(abs(fd)), maxval(abs(y_dot)), tiny(1.0_dp))
        jvp_error = maxval(abs(fd - y_dot)) / scale
        call diag%vjp(sample, eta, y_bar, q_bar, eta_bar)
        lhs = dot_product(y_bar, y_dot)
        rhs = sum(q_bar * q_dot) + dot_product(eta_bar, eta_dot)
        adjoint_error = abs(lhs - rhs) / max(abs(lhs), abs(rhs), tiny(1.0_dp))
    end subroutine check_diagnostic_derivatives

end module tiago_profile_diagnostic
