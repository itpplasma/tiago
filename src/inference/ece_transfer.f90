module tiago_ece_transfer
    !! Thermal transfer on a supplied reciprocal ray; cells run observer first.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_ece_loss, only: ece_loss
    use tiago_ece_loss_small, only: ece_loss_small
    use tiago_ece_thermal_cell, only: ece_thermal_cell
    use tiago_ece_reflection_loss, only: ece_reflection_loss
    use tiago_ece_two_stream_wall, only: ece_two_stream_wall
    implicit none
    private

    type, public :: ece_transfer_t
        private
        logical, public :: valid = .false.
        real(dp), public :: temperature = 0.0_dp, optical_depth = 0.0_dp
        real(dp), allocatable :: alpha(:), length(:), dloss(:)
        real(dp), allocatable :: near_row(:, :), far_row(:, :)
        real(dp) :: wall_row(5) = 0.0_dp, total_dloss = 0.0_dp
    contains
        procedure :: jvp => transfer_jvp
        procedure :: vjp => transfer_vjp
    end type

    public :: thermal_transfer, thermal_transfer_batch
contains
    pure subroutine attenuation(tau, loss, derivative)
        real(dp), intent(in) :: tau
        real(dp), intent(out) :: loss, derivative
        real(dp) :: row(1)
        if (tau < 1.0e-3_dp) then
            call ece_loss_small(tau, loss, row)
        else
            call ece_loss(tau, loss, row)
        end if
        derivative = row(1)
    end subroutine

    subroutine thermal_transfer(te, alpha, length, rho, wall_te, trace, status)
        real(dp), intent(in) :: te(:), alpha(:), length(:), rho, wall_te
        type(ece_transfer_t), intent(out) :: trace
        integer, intent(out) :: status
        real(dp) :: loss(size(te)), tau(size(te)), total_loss, en, ef, value, &
                    wall_loss, reflection_row(1)
        integer :: j, n

        status = 1
        n = size(te)
        if (size(alpha) /= n .or. size(length) /= n) return
        if (.not. all(ieee_is_finite(te))) return
        if (.not. all(ieee_is_finite(alpha))) return
        if (.not. all(ieee_is_finite(length))) return
        if (.not. all(ieee_is_finite([rho, wall_te]))) return
        if (any(te < 0.0_dp) .or. any(alpha < 0.0_dp)) return
        if (any(length < 0.0_dp)) return
        if (rho < 0.0_dp .or. rho >= 1.0_dp .or. wall_te < 0.0_dp) return
        tau = alpha*length
        if (.not. all(ieee_is_finite(tau))) return
        trace%optical_depth = sum(tau)
        if (.not. ieee_is_finite(trace%optical_depth)) return
        trace%alpha = alpha
        trace%length = length
        allocate (trace%dloss(n), trace%near_row(3, n), trace%far_row(3, n))
        do j = 1, n
            call attenuation(tau(j), loss(j), trace%dloss(j))
        end do
        en = 0.0_dp
        do j = n, 1, -1
            call ece_thermal_cell(en, te(j), loss(j), value, trace%near_row(:, j))
            en = value
        end do
        ef = 0.0_dp
        do j = 1, n
            call ece_thermal_cell(ef, te(j), loss(j), value, trace%far_row(:, j))
            ef = value
        end do
        call attenuation(trace%optical_depth, total_loss, trace%total_dloss)
        call ece_reflection_loss(rho, wall_loss, reflection_row)
        call ece_two_stream_wall(en, ef, total_loss, wall_loss, wall_te, &
                                 trace%temperature, trace%wall_row)
        trace%wall_row(4) = trace%wall_row(4)*reflection_row(1)
        if (.not. ieee_is_finite(trace%temperature)) return
        if (.not. all(ieee_is_finite(trace%wall_row))) return
        trace%valid = .true.
        status = 0
    end subroutine

    subroutine thermal_transfer_batch(te, alpha, length, rho, wall_te, trace, status)
        real(dp), intent(in) :: te(:, :), alpha(:, :), length(:, :), rho(:), wall_te(:)
        type(ece_transfer_t), intent(out) :: trace(:)
        integer, intent(out) :: status(:)
        integer :: j, n
        status = 1
        n = size(te, 2)
        if (any(shape(alpha) /= shape(te))) return
        if (any(shape(length) /= shape(te))) return
        if (size(rho) /= n .or. size(wall_te) /= n) return
        if (size(trace) /= n .or. size(status) /= n) return
        !$omp parallel do schedule(static)
        do j = 1, n
            call thermal_transfer(te(:, j), alpha(:, j), length(:, j), &
                                  rho(j), wall_te(j), trace(j), status(j))
        end do
        !$omp end parallel do
    end subroutine

    subroutine transfer_jvp(self, te_dot, alpha_dot, length_dot, rho_dot, wall_dot, &
        result)
        class(ece_transfer_t), intent(in) :: self
        real(dp), intent(in) :: te_dot(:), alpha_dot(:), length_dot(:), rho_dot, &
            wall_dot
        real(dp), intent(out) :: result
        real(dp) :: dtau(size(te_dot)), dl(size(te_dot)), den, def
        integer :: j, n
        if (.not. self%valid) error stop 'ECE transfer: invalid derivative tape'
        n = size(self%alpha)
        if (size(te_dot) /= n .or. size(alpha_dot) /= n .or. size(length_dot) /= n) &
            error stop 'ECE transfer: tangent shape'
        dtau = self%length*alpha_dot + self%alpha*length_dot
        dl = self%dloss*dtau
        den = 0.0_dp
        do j = n, 1, -1
            den = dot_product(self%near_row(:, j), [den, te_dot(j), dl(j)])
        end do
        def = 0.0_dp
        do j = 1, n
            def = dot_product(self%far_row(:, j), [def, te_dot(j), dl(j)])
        end do
        result = dot_product(self%wall_row, &
                             [den, def, self%total_dloss*sum(dtau), rho_dot, wall_dot])
    end subroutine

    subroutine transfer_vjp(self, seed, te_bar, alpha_bar, length_bar, rho_bar, &
        wall_bar)
        class(ece_transfer_t), intent(in) :: self
        real(dp), intent(in) :: seed
        real(dp), intent(out) :: te_bar(:), alpha_bar(:), length_bar(:), rho_bar, &
            wall_bar
        real(dp) :: lb(size(te_bar)), tb(size(te_bar)), enb, efb
        integer :: j, n
        if (.not. self%valid) error stop 'ECE transfer: invalid derivative tape'
        n = size(self%alpha)
        if (size(te_bar) /= n .or. size(alpha_bar) /= n .or. size(length_bar) /= n) &
            error stop 'ECE transfer: adjoint shape'
        te_bar = 0.0_dp
        lb = 0.0_dp
        enb = seed*self%wall_row(1)
        do j = 1, n
            te_bar(j) = te_bar(j) + enb*self%near_row(2, j)
            lb(j) = lb(j) + enb*self%near_row(3, j)
            enb = enb*self%near_row(1, j)
        end do
        efb = seed*self%wall_row(2)
        do j = n, 1, -1
            te_bar(j) = te_bar(j) + efb*self%far_row(2, j)
            lb(j) = lb(j) + efb*self%far_row(3, j)
            efb = efb*self%far_row(1, j)
        end do
        tb = lb*self%dloss + seed*self%wall_row(3)*self%total_dloss
        alpha_bar = tb*self%length
        length_bar = tb*self%alpha
        rho_bar = seed*self%wall_row(4)
        wall_bar = seed*self%wall_row(5)
    end subroutine
end module
