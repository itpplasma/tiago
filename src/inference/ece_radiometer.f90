module tiago_ece_radiometer
    !! Electron-cyclotron emission, cold-resonance model, as a profile
    !! diagnostic.
    !!
    !! Channel i (frequency f_i, harmonic n_i) samples its line of sight from
    !! the antenna at t_j = j L / nscan. Its cold resonance is the first
    !! interval [t_{j-1}, t_j] inside the plasma, counted from the antenna,
    !! where |B| - B_res changes sign, B_res = 2 pi m_e f_i / (n_i e). Linear
    !! interpolation in t locates the resonance,
    !!     a = (B_res - B_{j-1}) / (B_j - B_{j-1}),
    !! and the signal is the radiation temperature of an optically thick,
    !! unbroadened layer, Te_res = (1 - a) Te_{j-1} + a Te_j. The derivatives
    !! include the motion of the resonance through q_bar(Q_B). The
    !! interpolation error is O(h^2) in the sample spacing h = L / nscan.
    !!
    !! Extension point (documented, not implemented): a radiation-transport
    !! model (TRAVIS-like: relativistic, Doppler-broadened absorption and
    !! emission, finite optical depth, refraction) replaces the single resonance
    !! by the weighted ray integral T_rad = int Te alpha exp(-tau) ds, a
    !! function of ne, Te and |B| at the same samples; see
    !! docs/profile-diagnostics.md.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t, Q_B, Q_TE
    use tiago_sight_geometry, only: to_rpz
    use tiago_cold_resonance, only: cold_resonance
    implicit none
    private

    type, public :: ece_channel_t
        character(len=:), allocatable :: label
        real(dp) :: frequency = 0.0_dp !! [Hz]
        integer :: harmonic = 2
        real(dp) :: origin(3) = 0.0_dp !! antenna, Cartesian [m]
        real(dp) :: direction(3) = [-1.0_dp, 0.0_dp, 0.0_dp]
        real(dp) :: length = 1.0_dp
        integer :: nscan = 1000
    end type ece_channel_t

    type, public, extends(profile_diagnostic_t) :: ece_t
        type(ece_channel_t), allocatable :: channels(:)
        integer, allocatable :: first(:) !! first sample of each channel
    contains
        procedure :: init => ece_init
        procedure :: n_signals => ece_n_signals
        procedure :: n_nuisance => ece_n_nuisance
        procedure :: label => ece_label
        procedure :: sample_points => ece_points
        procedure :: predict => ece_predict
        procedure :: jvp => ece_jvp
        procedure :: vjp => ece_vjp
        procedure :: resonance => ece_resonance
    end type ece_t

contains

    subroutine ece_init(self, channels)
        class(ece_t), intent(inout) :: self
        type(ece_channel_t), intent(in) :: channels(:)
        integer :: c

        self%channels = channels
        if (allocated(self%first)) deallocate(self%first)
        allocate(self%first(size(channels) + 1))
        self%first(1) = 1
        do c = 1, size(channels)
            self%first(c + 1) = self%first(c) + channels(c)%nscan + 1
        end do
    end subroutine ece_init

    integer function ece_n_signals(self)
        class(ece_t), intent(in) :: self
        ece_n_signals = size(self%channels)
    end function ece_n_signals

    integer function ece_n_nuisance(self)
        class(ece_t), intent(in) :: self
        ece_n_nuisance = 0
    end function ece_n_nuisance

    function ece_label(self, i) result(text)
        class(ece_t), intent(in) :: self
        integer, intent(in) :: i
        character(len=:), allocatable :: text
        text = 'Te_ECE '//self%channels(i)%label
    end function ece_label

    subroutine ece_points(self, eta, x, dx_deta)
        class(ece_t), intent(in) :: self
        real(dp), intent(in) :: eta(:)
        real(dp), allocatable, intent(out) :: x(:, :), dx_deta(:, :, :)
        real(dp) :: d(3)
        integer :: c, j

        allocate(x(3, self%first(size(self%channels) + 1) - 1))
        allocate(dx_deta(3, size(x, 2), 0))
        do c = 1, size(self%channels)
            d = self%channels(c)%direction / norm2(self%channels(c)%direction)
            do j = 0, self%channels(c)%nscan
                x(:, self%first(c) + j) = to_rpz(self%channels(c)%origin + &
                    self%channels(c)%length * real(j, dp) / self%channels(c)%nscan * d)
            end do
        end do
    end subroutine ece_points

    subroutine ece_resonance(self, sample, c, k, a, b_res)
        !! Bracketing samples k - 1, k and fraction a of channel c; k = 0 if
        !! the channel has no resonance inside the plasma.
        class(ece_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        integer, intent(in) :: c
        integer, intent(out) :: k
        real(dp), intent(out) :: a, b_res
        real(dp) :: db_res_df, g0, g1
        integer :: j

        call cold_resonance(self%channels(c)%frequency, real(self%channels(c)%harmonic, dp), &
            b_res, db_res_df)
        k = 0
        a = 0.0_dp
        do j = self%first(c) + 1, self%first(c + 1) - 1
            if (.not. sample%inside(j - 1)) cycle
            if (.not. sample%inside(j)) cycle
            g0 = sample%q(Q_B, j - 1) - b_res
            g1 = sample%q(Q_B, j) - b_res
            if (g0 == 0.0_dp .or. g0 * g1 < 0.0_dp) then
                k = j
                a = g0 / (g0 - g1)
                return
            end if
        end do
    end subroutine ece_resonance

    subroutine ece_predict(self, sample, eta, y)
        class(ece_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: a, b_res
        integer :: c, k

        do c = 1, size(self%channels)
            call self%resonance(sample, c, k, a, b_res)
            y(c) = 0.0_dp
            if (k > 0) y(c) = (1.0_dp - a) * sample%q(Q_TE, k - 1) + a * sample%q(Q_TE, k)
        end do
    end subroutine ece_predict

    subroutine ece_jvp(self, sample, eta, q_dot, eta_dot, y_dot)
        class(ece_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        real(dp) :: a, b_res, da
        integer :: c, k

        do c = 1, size(self%channels)
            call self%resonance(sample, c, k, a, b_res)
            y_dot(c) = 0.0_dp
            if (k == 0) cycle
            da = -((1.0_dp - a) * q_dot(Q_B, k - 1) + a * q_dot(Q_B, k)) / &
                (sample%q(Q_B, k) - sample%q(Q_B, k - 1))
            y_dot(c) = (1.0_dp - a) * q_dot(Q_TE, k - 1) + a * q_dot(Q_TE, k) &
                + (sample%q(Q_TE, k) - sample%q(Q_TE, k - 1)) * da
        end do
    end subroutine ece_jvp

    subroutine ece_vjp(self, sample, eta, y_bar, q_bar, eta_bar)
        class(ece_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        real(dp) :: a, b_res, abar
        integer :: c, k

        q_bar = 0.0_dp
        eta_bar = 0.0_dp
        do c = 1, size(self%channels)
            call self%resonance(sample, c, k, a, b_res)
            if (k == 0) cycle
            q_bar(Q_TE, k - 1) = q_bar(Q_TE, k - 1) + (1.0_dp - a) * y_bar(c)
            q_bar(Q_TE, k) = q_bar(Q_TE, k) + a * y_bar(c)
            abar = (sample%q(Q_TE, k) - sample%q(Q_TE, k - 1)) * y_bar(c) / &
                (sample%q(Q_B, k) - sample%q(Q_B, k - 1))
            q_bar(Q_B, k - 1) = q_bar(Q_B, k - 1) - (1.0_dp - a) * abar
            q_bar(Q_B, k) = q_bar(Q_B, k) - a * abar
        end do
    end subroutine ece_vjp
end module tiago_ece_radiometer
