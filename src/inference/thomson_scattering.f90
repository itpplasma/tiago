module tiago_thomson_scattering
    !! Thomson scattering as a profile diagnostic. Channel i is a scattering
    !! volume of length L_i along the laser, sampled by Gauss-Legendre points
    !! with normalized weights w_k:
    !!     ne_i = sum_k w_k ne_k,    Te_i = sum_k w_k ne_k Te_k / ne_i,
    !! since each part of the volume scatters in proportion to its density.
    !! Signals are (ne_1 .. ne_n, Te_1 .. Te_n). With aligned = .false. a rigid
    !! shift (dR, dZ) of all volumes is the nuisance eta, as for the point probe.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t, Q_NE, Q_TE
    use tiago_sight_geometry, only: gauss_legendre01, to_rpz
    implicit none
    private

    type, public :: thomson_volume_t
        character(len=:), allocatable :: label
        real(dp) :: center(3) = 0.0_dp !! Cartesian [m]
        real(dp) :: direction(3) = [0.0_dp, 0.0_dp, 1.0_dp] !! laser direction
        real(dp) :: length = 0.0_dp !! [m]
        integer :: npoints = 3
    end type thomson_volume_t

    type, public, extends(profile_diagnostic_t) :: thomson_t
        type(thomson_volume_t), allocatable :: volumes(:)
        logical :: aligned = .true.
        real(dp), allocatable :: x0(:, :), w(:) !! sample points (R, phi, Z), weights
        integer, allocatable :: owner(:)
    contains
        procedure :: init => thomson_init
        procedure :: n_signals => thomson_n_signals
        procedure :: n_nuisance => thomson_n_nuisance
        procedure :: label => thomson_label
        procedure :: sample_points => thomson_points
        procedure :: predict => thomson_predict
        procedure :: jvp => thomson_jvp
        procedure :: vjp => thomson_vjp
    end type thomson_t

contains

    subroutine thomson_init(self, volumes, aligned)
        class(thomson_t), intent(inout) :: self
        type(thomson_volume_t), intent(in) :: volumes(:)
        logical, intent(in), optional :: aligned
        real(dp), allocatable :: xg(:), wg(:)
        real(dp) :: d(3)
        integer :: c, k, m, n, np

        self%volumes = volumes
        if (present(aligned)) self%aligned = aligned
        n = 0
        do c = 1, size(volumes)
            n = n + npts(volumes(c))
        end do
        if (allocated(self%x0)) deallocate(self%x0, self%w, self%owner)
        allocate(self%x0(3, n), self%w(n), self%owner(n))
        m = 0
        do c = 1, size(volumes)
            np = npts(volumes(c))
            call gauss_legendre01(np, xg, wg)
            d = volumes(c)%direction / norm2(volumes(c)%direction)
            do k = 1, np
                m = m + 1
                self%x0(:, m) = to_rpz(volumes(c)%center + &
                    (xg(k) - 0.5_dp) * volumes(c)%length * d)
                self%w(m) = wg(k)
                self%owner(m) = c
            end do
        end do
    contains
        pure integer function npts(v)
            type(thomson_volume_t), intent(in) :: v
            npts = merge(max(v%npoints, 1), 1, v%length > 0.0_dp)
        end function npts
    end subroutine thomson_init

    integer function thomson_n_signals(self)
        class(thomson_t), intent(in) :: self
        thomson_n_signals = 2 * size(self%volumes)
    end function thomson_n_signals

    integer function thomson_n_nuisance(self)
        class(thomson_t), intent(in) :: self
        thomson_n_nuisance = merge(0, 2, self%aligned)
    end function thomson_n_nuisance

    function thomson_label(self, i) result(text)
        class(thomson_t), intent(in) :: self
        integer, intent(in) :: i
        character(len=:), allocatable :: text
        integer :: n
        n = size(self%volumes)
        if (i <= n) then
            text = 'ne '//self%volumes(i)%label
        else
            text = 'Te '//self%volumes(i - n)%label
        end if
    end function thomson_label

    subroutine thomson_points(self, eta, x, dx_deta)
        class(thomson_t), intent(in) :: self
        real(dp), intent(in) :: eta(:)
        real(dp), allocatable, intent(out) :: x(:, :), dx_deta(:, :, :)

        x = self%x0
        allocate(dx_deta(3, size(x, 2), self%n_nuisance()))
        if (self%aligned) return
        x(1, :) = x(1, :) + eta(1)
        x(3, :) = x(3, :) + eta(2)
        dx_deta = 0.0_dp
        dx_deta(1, :, 1) = 1.0_dp
        dx_deta(3, :, 2) = 1.0_dp
    end subroutine thomson_points

    subroutine thomson_predict(self, sample, eta, y)
        class(thomson_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: nt(size(self%volumes))
        integer :: k, c, n

        n = size(self%volumes)
        y(:2 * n) = 0.0_dp
        nt = 0.0_dp
        do k = 1, size(self%w)
            c = self%owner(k)
            y(c) = y(c) + self%w(k) * sample%q(Q_NE, k)
            nt(c) = nt(c) + self%w(k) * sample%q(Q_NE, k) * sample%q(Q_TE, k)
        end do
        where (y(:n) > 0.0_dp) y(n + 1:2 * n) = nt / y(:n)
    end subroutine thomson_predict

    subroutine thomson_jvp(self, sample, eta, q_dot, eta_dot, y_dot)
        class(thomson_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        real(dp) :: y(2 * size(self%volumes)), dnt(size(self%volumes))
        integer :: k, c, n

        n = size(self%volumes)
        call self%predict(sample, eta, y)
        y_dot(:2 * n) = 0.0_dp
        dnt = 0.0_dp
        do k = 1, size(self%w)
            c = self%owner(k)
            y_dot(c) = y_dot(c) + self%w(k) * q_dot(Q_NE, k)
            dnt(c) = dnt(c) + self%w(k) * (q_dot(Q_NE, k) * sample%q(Q_TE, k) &
                + sample%q(Q_NE, k) * q_dot(Q_TE, k))
        end do
        do c = 1, n
            if (y(c) > 0.0_dp) y_dot(n + c) = (dnt(c) - y(n + c) * y_dot(c)) / y(c)
        end do
    end subroutine thomson_jvp

    subroutine thomson_vjp(self, sample, eta, y_bar, q_bar, eta_bar)
        class(thomson_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        real(dp) :: y(2 * size(self%volumes)), nbar(size(self%volumes)), ntbar(size(self%volumes))
        integer :: k, c, n

        n = size(self%volumes)
        call self%predict(sample, eta, y)
        nbar = y_bar(:n)
        ntbar = 0.0_dp
        do c = 1, n
            if (y(c) > 0.0_dp) then
                nbar(c) = nbar(c) - y_bar(n + c) * y(n + c) / y(c)
                ntbar(c) = y_bar(n + c) / y(c)
            end if
        end do
        q_bar = 0.0_dp
        eta_bar = 0.0_dp
        do k = 1, size(self%w)
            c = self%owner(k)
            q_bar(Q_NE, k) = self%w(k) * (nbar(c) + ntbar(c) * sample%q(Q_TE, k))
            q_bar(Q_TE, k) = self%w(k) * ntbar(c) * sample%q(Q_NE, k)
        end do
    end subroutine thomson_vjp
end module tiago_thomson_scattering
