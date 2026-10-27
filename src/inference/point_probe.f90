module tiago_point_probe
    !! Reference profile diagnostic: ideal local measurements of ne, Te or Ti
    !! at fixed points, with a rigid alignment shift (dR, dZ) of all points as
    !! nuisance parameters. It exercises the profile-diagnostic contract in
    !! synthetic tests; physical instruments (Thomson scattering volumes and
    !! spectral response, ECE, line integrals) are separate implementations.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t, Q_NE, Q_TE, Q_TI
    implicit none
    private

    type, public, extends(profile_diagnostic_t) :: point_probe_t
        real(dp), allocatable :: x0(:, :)        !! (3, n) nominal points R, phi, Z
        integer, allocatable :: quantity(:)      !! Q_NE, Q_TE or Q_TI per signal
        logical :: aligned = .true.              !! .false.: eta = (dR, dZ) nuisance
    contains
        procedure :: n_signals => probe_n_signals
        procedure :: n_nuisance => probe_n_nuisance
        procedure :: label => probe_label
        procedure :: sample_points => probe_points
        procedure :: predict => probe_predict
        procedure :: jvp => probe_jvp
        procedure :: vjp => probe_vjp
    end type point_probe_t

contains

    integer function probe_n_signals(self)
        class(point_probe_t), intent(in) :: self
        probe_n_signals = size(self%quantity)
    end function probe_n_signals

    integer function probe_n_nuisance(self)
        class(point_probe_t), intent(in) :: self
        probe_n_nuisance = merge(0, 2, self%aligned)
    end function probe_n_nuisance

    function probe_label(self, i) result(text)
        class(point_probe_t), intent(in) :: self
        integer, intent(in) :: i
        character(len=:), allocatable :: text
        character(len=48) :: buffer
        character(len=2), parameter :: names(3) = ['ne', 'Te', 'Ti']
        write(buffer, '(A,"(R=",F6.3,",Z=",F6.3,")")') names(self%quantity(i) - Q_NE + 1), &
            self%x0(1, i), self%x0(3, i)
        text = trim(buffer)
    end function probe_label

    subroutine probe_points(self, eta, x, dx_deta)
        class(point_probe_t), intent(in) :: self
        real(dp), intent(in) :: eta(:)
        real(dp), allocatable, intent(out) :: x(:, :), dx_deta(:, :, :)
        integer :: n

        n = size(self%quantity)
        x = self%x0
        allocate(dx_deta(3, n, self%n_nuisance()))
        if (self%aligned) return
        x(1, :) = x(1, :) + eta(1)
        x(3, :) = x(3, :) + eta(2)
        dx_deta = 0.0_dp
        dx_deta(1, :, 1) = 1.0_dp
        dx_deta(3, :, 2) = 1.0_dp
    end subroutine probe_points

    subroutine probe_predict(self, sample, eta, y)
        class(point_probe_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        integer :: i
        do i = 1, size(self%quantity)
            y(i) = sample%q(self%quantity(i), i)
        end do
    end subroutine probe_predict

    subroutine probe_jvp(self, sample, eta, q_dot, eta_dot, y_dot)
        class(point_probe_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        integer :: i
        do i = 1, size(self%quantity)
            y_dot(i) = q_dot(self%quantity(i), i)
        end do
    end subroutine probe_jvp

    subroutine probe_vjp(self, sample, eta, y_bar, q_bar, eta_bar)
        class(point_probe_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        integer :: i
        q_bar = 0.0_dp
        eta_bar = 0.0_dp
        do i = 1, size(self%quantity)
            q_bar(self%quantity(i), i) = y_bar(i)
        end do
    end subroutine probe_vjp

end module tiago_point_probe
