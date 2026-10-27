module tiago_magnetic_response
    !! Exact affine coefficient response at fixed boundary/coil geometry.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t
    use tiago_sensor_set, only: sensor_set_t
    use tiago_vacuum_forward, only: vacuum_solver_t, quadrature_rule_t
    use tiago_sensor_signals, only: sensor_set_signals, sensor_set_coil_jacobian, &
                                  sensor_set_plasma_jacobian, sensor_set_linked_flux_row
    implicit none
    private
    public :: magnetic_response_t
    type, extends(profile_diagnostic_t) :: magnetic_response_t
        private
        type(sensor_set_t) :: sensors
        real(dp), allocatable :: base(:), response(:, :)
    contains
        procedure :: init
        procedure :: n_signals
        procedure :: n_nuisance
        procedure :: label
        procedure :: sample_points
        procedure :: predict
        procedure :: jvp
        procedure :: vjp
    end type
contains
    subroutine init(self, solver, sensors, rule)
        class(magnetic_response_t), intent(inout) :: self
        type(vacuum_solver_t), intent(in) :: solver
        type(sensor_set_t), intent(in) :: sensors
        type(quadrature_rule_t), optional, intent(in) :: rule
        real(dp), allocatable :: coil(:, :), modes(:, :), flux(:)
        integer :: n, nc, nm
        real(dp) :: linked
        self%sensors = sensors
        linked = solver%plasma_flux(required=any(sensors%loops%idia == 1))
        call sensor_set_signals(solver, sensors, self%base, rule, linked)
        call sensor_set_coil_jacobian(solver, sensors, coil, rule)
        call sensor_set_linked_flux_row(sensors, flux)
        n = size(self%base)
        ! External field sources have no coefficient derivative contract.
        allocate (modes(n, 0))
        if (.not. allocated(solver%field_source)) then
            if (solver%plasma%enabled) then
                call sensor_set_plasma_jacobian(solver, sensors, modes, rule=rule)
            end if
        end if
        nc = size(coil, 2)
        nm = size(modes, 2)
        if (allocated(self%response)) deallocate (self%response)
        allocate (self%response(n, nc + nm + 1))
        self%response(:, :nc) = coil
        self%response(:, nc + 1:nc + nm) = modes
        self%response(:, nc + nm + 1) = flux
    end subroutine

    integer function n_signals(self)
        class(magnetic_response_t), intent(in) :: self
        n_signals = size(self%base)
    end function
    integer function n_nuisance(self)
        class(magnetic_response_t), intent(in) :: self
        n_nuisance = size(self%response, 2)
    end function
    function label(self, i) result(text)
        class(magnetic_response_t), intent(in) :: self
        integer, intent(in) :: i
        character(:), allocatable :: text
        text = self%sensors%labels(i)%s
    end function
    subroutine sample_points(self, eta, x, dx_deta)
        class(magnetic_response_t), intent(in) :: self
        real(dp), intent(in) :: eta(:)
        real(dp), allocatable, intent(out) :: x(:, :), dx_deta(:, :, :)
        allocate (x(3, 0), dx_deta(3, 0, self%n_nuisance()))
    end subroutine
    subroutine predict(self, sample, eta, y)
        class(magnetic_response_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        y = self%base + matmul(self%response, eta)
    end subroutine
    subroutine jvp(self, sample, eta, q_dot, eta_dot, y_dot)
        class(magnetic_response_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        y_dot = matmul(self%response, eta_dot)
    end subroutine
    subroutine vjp(self, sample, eta, y_bar, q_bar, eta_bar)
        class(magnetic_response_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        q_bar = 0.0_dp
        eta_bar = matmul(transpose(self%response), y_bar)
    end subroutine
end module
