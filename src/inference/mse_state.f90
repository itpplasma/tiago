module tiago_mse_state
    !! MSE pitch angles from the interior field of a real-space state.
    !!
    !! The channel midpoint (docs/mse.md) is located in the state's labels u,
    !! B is the state's vector field there (real_space_state_t%vector_jets),
    !! and tiago_mse turns it into the pi-polarization angle. The derivative
    !! with respect to the state parameters p holds the midpoint fixed in
    !! real space, so the labels move with p (implicit-function theorem):
    !!     dB/dp|_x = dB/dp|_u - dB/du J^{-1} dX/dp|_u,    J = dX/du.
    !! Unlike the virtual-casing source, this field is valid inside the plasma
    !! only; midpoints outside the boundary are errors.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    use tiago_diagnostic_types, only: mse_channel_t
    use tiago_mse, only: mse_geometry, mse_angles
    use tiago_real_space_state, only: real_space_state_t, located_t, vector_jets_t, locate_ok
    implicit none
    private

    public :: mse_state_angles

contains

    subroutine mse_state_angles(state, channels, angles, ierr, message, param_jacobian, &
            field, rho)
        class(real_space_state_t), intent(in) :: state
        type(mse_channel_t), intent(in) :: channels(:)
        real(dp), allocatable, intent(out) :: angles(:)
        integer(i32), intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        !> d gamma / d p, (nchannel, n_params) [rad per parameter]
        real(dp), allocatable, intent(out), optional :: param_jacobian(:, :)
        !> Cartesian B at the midpoints, (3, nchannel) [T], and their rho
        real(dp), allocatable, intent(out), optional :: field(:, :), rho(:)
        type(located_t), allocatable :: loc(:)
        type(vector_jets_t), allocatable :: jets(:)
        real(dp) :: point(3), coefficients(6)
        real(dp), allocatable :: b(:, :), g(:, :), r(:)
        integer :: j, status
        character(len=16) :: num

        allocate(b(3, size(channels)), r(size(channels)), jets(size(channels)), loc(size(channels)))
        do j = 1, size(channels)
            write(num, '(I0)') j
            call mse_geometry(channels(j), point, coefficients, ierr, message)
            if (ierr /= 0) then
                message = 'MSE channel '//trim(num)//': '//message
                return
            end if
            ierr = 1_i32
            call state%locate(point, loc(j))
            if (loc(j)%status /= locate_ok .or. loc(j)%u(1) > 1.0_dp) then
                message = 'MSE channel '//trim(num)//': midpoint outside the plasma'
                return
            end if
            call state%vector_jets(loc(j)%u, jets(j), status)
            if (status /= locate_ok) then
                message = 'MSE channel '//trim(num)//': the state supplies no vector field'
                return
            end if
            b(:, j) = jets(j)%b
            r(j) = loc(j)%u(1)
        end do
        call mse_angles(channels, b, angles, ierr, message, field_jacobian=g)
        if (ierr /= 0) return
        if (present(field)) field = b
        if (present(rho)) rho = r
        if (.not. present(param_jacobian)) return
        allocate(param_jacobian(size(channels), state%n_params()))
        do j = 1, size(channels)
            associate (jet => jets(j))
                param_jacobian(j, :) = matmul(g(:, j), &
                    jet%bp - matmul(jet%bu, matmul(loc(j)%jinv, jet%xp)))
            end associate
        end do
    end subroutine mse_state_angles
end module tiago_mse_state
