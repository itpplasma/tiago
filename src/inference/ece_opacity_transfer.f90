module tiago_ece_opacity_transfer
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_ece_opacity, only: ece_opacity_batch, ECE_X1_FINITE
    use tiago_ece_transfer, only: thermal_transfer, ece_transfer_t
    implicit none
    private
    type, public :: ece_opacity_transfer_t
        private
        type(ece_transfer_t) :: transfer
        real(dp), allocatable :: alpha_row(:, :)
        logical, public :: valid = .false.
        real(dp), public :: temperature = 0.0_dp, optical_depth = 0.0_dp
    contains
        procedure :: jvp => opacity_transfer_jvp
        procedure :: vjp => opacity_transfer_vjp
    end type
    public :: thermal_ece_transfer
contains
    subroutine thermal_ece_transfer(te, ne, bfield, frequency, mode, harmonic, length, &
                                    rho, wall_te, trace, status, x1_regime)
        real(dp), intent(in) :: te(:), ne(:), bfield(:), frequency, length(:), rho, wall_te
        integer, intent(in) :: mode, harmonic
        type(ece_opacity_transfer_t), intent(out) :: trace
        integer, intent(out) :: status
        integer, intent(in), optional :: x1_regime
        real(dp) :: alpha(size(te)), emission(size(te)), n2(size(te)), rows(3, 4, size(te))
        integer :: codes(size(te)), regime, n
        status = 1
        n = size(te)
        if (size(ne) /= n .or. size(bfield) /= n .or. size(length) /= n) return
        regime = ECE_X1_FINITE
        if (present(x1_regime)) regime = x1_regime
        call ece_opacity_batch(te, ne, bfield, spread(frequency, 1, n), spread(mode, 1, n), &
                               spread(harmonic, 1, n), alpha, emission, n2, rows, codes, regime)
        if (any(codes /= 0)) then
            status = maxval(codes)
            return
        end if
        call thermal_transfer(te, alpha, length, rho, wall_te, trace%transfer, status)
        if (status /= 0) return
        trace%alpha_row = rows(1, :, :)
        trace%temperature = trace%transfer%temperature
        trace%optical_depth = trace%transfer%optical_depth
        trace%valid = .true.
    end subroutine
    subroutine opacity_transfer_jvp(self, te_dot, ne_dot, b_dot, frequency_dot, &
                                    length_dot, rho_dot, wall_dot, result)
        class(ece_opacity_transfer_t), intent(in) :: self
        real(dp), intent(in) :: te_dot(:), ne_dot(:), b_dot(:), frequency_dot, &
                                length_dot(:), rho_dot, wall_dot
        real(dp), intent(out) :: result
        real(dp) :: alpha_dot(size(te_dot))
        integer :: n
        if (.not. self%valid) error stop 'ECE opacity: invalid transfer tape'
        n = size(self%alpha_row, 2)
        if (size(te_dot) /= n .or. size(ne_dot) /= n .or. size(b_dot) /= n) &
            error stop 'ECE opacity: tangent shape'
        alpha_dot = self%alpha_row(1, :)*te_dot + self%alpha_row(2, :)*ne_dot + &
                    self%alpha_row(3, :)*b_dot + self%alpha_row(4, :)*frequency_dot
        call self%transfer%jvp(te_dot, alpha_dot, length_dot, rho_dot, wall_dot, result)
    end subroutine
    subroutine opacity_transfer_vjp(self, seed, te_bar, ne_bar, b_bar, frequency_bar, &
                                    length_bar, rho_bar, wall_bar)
        class(ece_opacity_transfer_t), intent(in) :: self
        real(dp), intent(in) :: seed
        real(dp), intent(out) :: te_bar(:), ne_bar(:), b_bar(:), frequency_bar, &
                                 length_bar(:), rho_bar, wall_bar
        real(dp) :: alpha_bar(size(te_bar))
        integer :: n
        if (.not. self%valid) error stop 'ECE opacity: invalid transfer tape'
        n = size(self%alpha_row, 2)
        if (size(te_bar) /= n .or. size(ne_bar) /= n .or. size(b_bar) /= n) &
            error stop 'ECE opacity: adjoint shape'
        call self%transfer%vjp(seed, te_bar, alpha_bar, length_bar, rho_bar, wall_bar)
        te_bar = te_bar + alpha_bar*self%alpha_row(1, :)
        ne_bar = alpha_bar*self%alpha_row(2, :)
        b_bar = alpha_bar*self%alpha_row(3, :)
        frequency_bar = sum(alpha_bar*self%alpha_row(4, :))
    end subroutine
end module
