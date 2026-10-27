module recorded_magnetic_test_source
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_field_source, only: plasma_field_source_t
    implicit none
    real(dp), parameter :: mu0 = 1.25663706212e-6_dp
    type, extends(plasma_field_source_t) :: analytic_source_t
        real(dp) :: field = 2.0_dp, pressure = 0.0_dp, radius = 0.5_dp
    contains
        procedure :: sample_bfield => field
        procedure :: sample_vector_potential => potential
    end type
contains
    subroutine field(self, points, bfield)
        class(analytic_source_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: bfield(:, :)
        real(dp) :: r2
        integer :: i
        bfield = 0.0_dp
        do i = 1, size(points, 2)
            r2 = sum(points(:2, i)**2)
            if (self%pressure == 0.0_dp) then
                bfield(3, i) = self%field
            else
                bfield(3, i) = sqrt(self%field**2 - &
                               2*mu0*self%pressure*max(0.0_dp, 1 - r2/self%radius**2)) &
                               - self%field
            end if
        end do
    end subroutine
    subroutine potential(self, points, avec)
        class(analytic_source_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: avec(:, :)
        real(dp) :: r2, limited, u, v, integral, a
        integer :: i
        do i = 1, size(points, 2)
            r2 = sum(points(:2, i)**2)
            if (self%pressure == 0.0_dp) then
                a = self%field/2
            else
                limited = min(r2, self%radius**2)
                v = sqrt(self%field**2 - 2*mu0*self%pressure)
                u = sqrt(v*v + 2*mu0*self%pressure*limited/self%radius**2)
                integral = limited*((u*u + u*v + v*v)/(3*(u + v)) &
                                    - self%field/2)
                a = integral/r2
            end if
            avec(:, i) = [-a*points(2, i), a*points(1, i), 0.0_dp]
        end do
    end subroutine
end module

program test_recorded_magnetics
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use recorded_magnetic_test_source, only: analytic_source_t, mu0
    use tiago_diagnostic_types, only: loop_point_t
    use tiago_sensor_set, only: sensor_set_t
    use tiago_vacuum_forward, only: vacuum_solver_t, quadrature_rule_t
    use tiago_sensor_signals, only: sensor_set_plasma_jacobian, sensor_set_signals
    use tiago_plasma_support, only: vmec_boundary_t, read_vmec_boundary
    use tiago_profile_diagnostic, only: plasma_sample_t, check_diagnostic_derivatives, &
                                        N_QUANTITIES, Q_NE, Q_TE, Q_TI, Q_ZEFF
    use tiago_calibration, only: calibration_prior_t
    use tiago_magnetic_response, only: magnetic_response_t
    use tiago_integrator_diagnostic, only: integrator_diagnostic_t
    use tiago_energy_consistency, only: energy_consistency_t
    implicit none
    integer, parameter :: n = 2048, nv = 64
    real(dp), parameter :: pi = acos(-1.0_dp), a = 0.5_dp, b0 = 2.0_dp, &
                           major = 10.0_dp, charge = 1.602176634e-19_dp
    type(analytic_source_t) :: source
    type(vmec_boundary_t) :: boundary
    type(vacuum_solver_t) :: solver, coil_solver
    type(sensor_set_t) :: sensors
    type(quadrature_rule_t) :: rule
    type(magnetic_response_t) :: magnetic
    type(integrator_diagnostic_t) :: recorded
    type(energy_consistency_t) :: energy
    type(calibration_prior_t) :: prior
    type(plasma_sample_t) :: empty, sample
    real(dp) :: angle, err, adj, area, reference(4), y(4), eta(4), &
                wdia, wkin, eeta(4), p0, exact_energy, covariance(3, 3), radial
    real(dp), allocatable :: coil_eta(:), state_eta(:), jac(:, :), actual(:)
    character(1024) :: coilfile, wout
    integer :: i, status
    call get_command_argument(1, coilfile)
    call get_command_argument(2, wout)
    allocate (sensors%loops(2), sensors%segs(1), sensors%probes(1), &
              sensors%labels(4), sensors%row_start(5), sensors%col(4), sensors%val(4))
    allocate (sensors%loops(1)%points(n))
    do i = 1, n
        angle = 2*pi*(i - 1)/n
        sensors%loops(1)%points(i) = loop_point_t(a*cos(angle), a*sin(angle), 0.0_dp)
    end do
    sensors%loops(2) = sensors%loops(1)
    sensors%loops(2)%idia = 1
    allocate (sensors%segs(1)%path(2), sensors%segs(1)%segment_area(1))
    sensors%segs(1)%path = [loop_point_t(0.0_dp, 0.0_dp, 0.0_dp), &
                            loop_point_t(0.0_dp, 0.0_dp, 1.0_dp)]
    sensors%segs(1)%segment_area = 0.4_dp
    sensors%probes(1)%position = [0.0_dp, 0.0_dp, 0.0_dp]
    sensors%probes(1)%normal = [0.0_dp, 0.0_dp, 1.0_dp]
    sensors%probes(1)%eff_area = 0.25_dp
    sensors%labels(1)%s = 'circle'
    sensors%labels(2)%s = 'diamagnetic'
    sensors%labels(3)%s = 'rogowski'
    sensors%labels(4)%s = 'probe'
    sensors%row_start = [1, 2, 3, 4, 5]
    sensors%col = [1, 2, 3, 4]
    sensors%val = 1.0_dp
    rule%gauss = .true.
    rule%samples_per_segment = 4
    source%enabled = .true.
    source%linked_flux_known = .true.
    source%diamagnetic_flux = 0.3_dp
    call solver%init('')
    call solver%set_plasma_source(source)
    call magnetic%init(solver, sensors, rule)
    empty%n = 0
    allocate (empty%q(N_QUANTITIES, 0))
    call magnetic%predict(empty, [0.0_dp], y)
    area = pi*a*a
    reference = [b0*area, b0*area + 0.3_dp, 0.8_dp, 0.5_dp]
    err = maxval(abs(y - reference))/maxval(abs(reference))
    call check(err < 2.0e-6_dp, 'circular loop / uniform Rogowski and B-probe')
    print '(A,ES14.6)', 'uniform field maximum relative analytic error: ', err
    covariance = reshape([0.04_dp, 0.01_dp, 0.0_dp, &
                         0.01_dp, 0.09_dp, 0.003_dp, 0.0_dp, 0.003_dp, 0.01_dp], [3, 3])
    call prior%init([1.7_dp, -0.2_dp, 0.08_dp], covariance, status)
    call check(status == 0, 'calibration prior')
    call recorded%init(magnetic, prior, [0.0_dp, 0.4_dp, 1.3_dp, 2.0_dp], &
                       0.4_dp, status)
    call check(status == 0, 'recorded magnetic adapter')
    eta = [0.02_dp, 0.3_dp, -0.4_dp, 0.2_dp]
    call check_diagnostic_derivatives(recorded, empty, eta, 1.0e-4_dp, err, adj)
    call check(err < 1.0e-7_dp, 'recorded magnetic finite differences')
    call check(adj < 1.0e-12_dp, 'recorded magnetic adjoint')
    print '(A,2ES14.6)', 'recorded magnetic FD / transpose relative errors: ', err, adj
    call recover_calibration()
    ! A different real solver is the consumer oracle for the coil response.
    call coil_solver%init(trim(coilfile))
    call magnetic%init(coil_solver, sensors, rule)
    allocate (coil_eta(magnetic%n_nuisance()))
    coil_eta = 0.0_dp
    coil_eta(:size(coil_eta) - 1) = -coil_solver%extcur
    call magnetic%predict(empty, coil_eta, y)
    call check(maxval(abs(y)) < 1.0e-14_dp, 'zero-current vacuum limit')
    coil_eta = 0.1_dp
    call check_diagnostic_derivatives(magnetic, empty, coil_eta, 1.0e-5_dp, err, adj)
    call check(err < 1.0e-8_dp .and. adj < 1.0e-12_dp, 'coil response derivatives')
    ! Fixed-boundary mode increments must reproduce a re-evaluated VMEC source.
    call solver%enable_plasma_from_vmec(trim(wout), 16, 16)
    call magnetic%init(solver, sensors, rule)
    allocate (state_eta(magnetic%n_nuisance()))
    state_eta = 0.0_dp
    call sensor_set_plasma_jacobian(solver, sensors, jac, rule=rule)
    call magnetic%predict(empty, state_eta, reference)
    state_eta(:size(jac, 2)) = 0.001_dp
    call magnetic%predict(empty, state_eta, y)
    call check(maxval(abs(y - reference - matmul(jac, state_eta(:size(jac, 2))))) &
               < 1.0e-12_dp, 'VMEC exact mode response composition')
    call read_vmec_boundary(trim(wout), boundary)
    boundary%bumnc = boundary%bumnc + 0.001_dp
    boundary%bvmnc = boundary%bvmnc + 0.001_dp
    if (boundary%lasym) then
        boundary%bumns = boundary%bumns + 0.001_dp
        boundary%bvmns = boundary%bvmns + 0.001_dp
    end if
    call solver%enable_plasma_from_boundary(boundary, 16, 16)
    call sensor_set_signals(solver, sensors, actual, rule)
    err = maxval(abs(y - actual))/maxval(abs(actual))
    call check(err < 1.0e-12_dp, 'coefficient increments vs rebuilt VMEC source')
    print '(A,ES14.6)', 'VMEC response vs rebuilt source relative error: ', err
    call check_diagnostic_derivatives(magnetic, empty, state_eta, 1.0e-4_dp, err, adj)
    call check(err < 1.0e-7_dp .and. adj < 1.0e-12_dp, 'VMEC response derivatives')
    print '(A,2ES14.6)', 'VMEC coefficient FD / transpose relative errors: ', err, adj
    ! Independent pressure-balance cylinder, Bz^2 + 2 mu0 p = B0^2.
    ! Pustovitov, NIFS annual report 2010, p64, Eq.(1):
    ! https://www.nifs.ac.jp/report/annrep10/pdf/064.pdf
    p0 = 1.0e-4_dp*b0*b0/(2*mu0)
    source%pressure = p0
    source%diamagnetic_flux = 0.0_dp
    call solver%set_plasma_source(source)
    call magnetic%init(solver, sensors, rule)
    call magnetic%predict(empty, [0.0_dp], y)
    eeta = [y(2), b0, major, 0.0_dp]
    sample%n = nv
  allocate (sample%q(N_QUANTITIES, nv), energy%points(3, nv), energy%volume_weights(nv))
    sample%q = 0.0_dp
    energy%points = 0.0_dp
    energy%volume_weights = 2*pi*pi*major*a*a/nv
    do i = 1, nv
        radial = (i - 0.5_dp)/nv
        sample%q(Q_NE, i) = p0*(1 - radial)/(charge*1000.0_dp)
        sample%q(Q_TE, i) = 600.0_dp
        sample%q(Q_TI, i) = 400.0_dp
        sample%q(Q_ZEFF, i) = 1.0_dp
    end do
    call energy%energies(sample, eeta, wdia, wkin)
    exact_energy = 1.5_dp*pi*pi*major*a*a*p0
    err = abs(wdia - wkin)/wkin
    call check(abs(wkin - exact_energy)/exact_energy < 1.0e-14_dp, &
               'thermal cylinder integral')
    call check(err < 2.0e-5_dp, 'large aspect ratio / low beta diamagnetic energy')
    print '(A,3ES14.6)', 'cylinder Wdia [J], Wkin [J], relative difference: ', &
        wdia, wkin, err
    call check_diagnostic_derivatives(energy, sample, eeta, 1.0e-5_dp, err, adj)
    call check(err < 1.0e-8_dp .and. adj < 1.0e-12_dp, 'energy derivatives')
    print '(A,2ES14.6)', 'energy FD / transpose relative errors: ', err, adj
    eeta(4) = 0.001_dp
    call energy%energies(sample, eeta, wdia, wkin)
    call check(wdia > wkin, 'explicit paramagnetic correction sign')
    sample%q = 0.0_dp
    sample%q(Q_ZEFF, :) = 1.0_dp
    eeta(1) = 0.0_dp
    eeta(4) = 0.0_dp
    call energy%energies(sample, eeta, wdia, wkin)
    call check(wdia == 0.0_dp .and. wkin == 0.0_dp, 'zero-pressure energy limit')
contains
    subroutine check(ok, text)
        logical, intent(in) :: ok
        character(*), intent(in) :: text
        if (ok) return
        print '(A)', 'FAIL: '//text
        error stop 1
    end subroutine
    subroutine recover_calibration()
        real(dp) :: truth(4), prediction(4), zero(4), direction(4), derivative(4)
        real(dp) :: design(4, 3), rhs(4, 1), work(128), recovered(3), physical(3)
        real(dp) :: qdot(N_QUANTITIES, 0)
        integer :: j, info
        truth = [0.0_dp, 0.7_dp, -0.4_dp, 0.3_dp]
        zero = 0.0_dp
        call recorded%predict(empty, truth, prediction)
        call recorded%predict(empty, zero, y)
        rhs(:, 1) = prediction - y
        do j = 1, 3
            direction = 0.0_dp
            direction(j + 1) = 1.0_dp
            call recorded%jvp(empty, zero, qdot, direction, derivative)
            design(:, j) = derivative
        end do
        call dgels('N', 4, 3, 1, design, 4, rhs, 4, work, size(work), info)
        call check(info == 0, 'synthetic calibration solve')
        recovered = prior%physical(rhs(:3, 1))
        physical = prior%physical(truth(2:))
        err = maxval(abs(recovered - physical))
        call check(err < 1.0e-12_dp, 'recover gain, offset, drift')
print '(A,ES14.6)', 'synthetic gain / offset / drift max absolute recovery error: ', err
    end subroutine
end program
