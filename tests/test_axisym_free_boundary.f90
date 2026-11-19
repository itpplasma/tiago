program test_axisym_free_boundary
    !! Recovery on FreeGS's TestTokamakSensor machine with the FreeGSNKE
    !! provider (tools/freegsnke_provider.py): noise-free signals of a planted
    !! equilibrium, coil currents measured to 1 %, a start off in the coil
    !! currents (0.3 %), 2 % in the plasma current and 20 % in the axis pressure. Checks
    !! that the re-solved derivatives agree at two steps, recovery within two local
    !! standard deviations and to 1e-3 relative, and the Laplace covariance.
    !! Usage: test_axisym_free_boundary <machine> <provider command> <workdir>
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_axisym_free_boundary, only: axisym_free_boundary_t
    use tiago_axisym_machine, only: read_axisym_machine
    use tiago_inference_engines, only: map_fit, laplace, laplace_t
    implicit none

    type(axisym_free_boundary_t) :: model
    type(laplace_t) :: lap
    character(len=4096) :: arg
    character(len=:), allocatable :: message
    real(dp), parameter :: truth(6) = [1.53020064498019026e+05_dp, 6.27160501916148351e+04_dp, &
        -9.80609577221965155e+04_dp, -5.75578364420235448e+04_dp, 2.0e5_dp, 1.0e3_dp]
    real(dp), allocatable :: s(:), theta(:), j1(:, :), j2(:, :)
    real(dp) :: sd, err
    integer :: ierr, k, steps, failures, nflux, npick
    logical :: ok, conv

    failures = 0
    call get_command_argument(1, arg)
    model%machine_path = trim(arg)
    call read_axisym_machine(model%machine_path, model%machine, ierr, message)
    if (ierr /= 0) error stop message
    call get_command_argument(2, arg)
    model%command = trim(arg)
    call get_command_argument(3, arg)
    model%workdir = trim(arg)
    call execute_command_line('mkdir -p '//model%workdir)
    model%fvac = 2.0_dp
    model%coil_current = truth(1:4)
    model%plasma_current = truth(5)
    model%pressure_axis = truth(6)
    model%names = [character(len=32) :: 'CURRENT(1)', 'CURRENT(2)', 'CURRENT(3)', &
        'CURRENT(4)', 'PLASMA_CURRENT', 'PRESSURE_AXIS']
    model%scale = [1.0e5_dp, 1.0e5_dp, 1.0e5_dp, 1.0e5_dp, 1.0e5_dp, 1.0e3_dp]
    model%prior_mean = [truth(1:4), 0.0_dp, 0.0_dp]
    model%prior_sigma = [0.01_dp * abs(truth(1:4)), 0.0_dp, 0.0_dp]

    theta = truth / model%scale
    call model%signals(theta, s, ok)
    if (.not. ok) error stop 'planted equilibrium did not converge'
    nflux = size(model%machine%sensors%flux)
    npick = size(model%machine%sensors%pickup)
    model%data = s
    allocate(model%sigma(size(s)))
    model%sigma(:nflux) = 1.0e-4_dp
    model%sigma(nflux + 1:nflux + npick) = 1.0e-4_dp
    model%sigma(nflux + npick + 1:) = 50.0_dp

    model%fd_relative = 3.0e-3_dp
    call model%jacobian(theta, j1)
    model%fd_relative = 1.0e-3_dp
    call model%jacobian(theta, j2)
    do k = 1, size(theta)
        print '(A,I0,ES10.3)', 'column ', k, maxval(abs(j1(:, k) - j2(:, k))) / maxval(abs(j1(:, k)))
    end do
    err = maxval(abs(j1 - j2)) / maxval(abs(j1))
    print '(A,ES10.3)', 'derivative step change (3e-3 against 1e-3): ', err
    call check(err < 2.0e-2_dp, 're-solved derivatives depend on the step')

    theta = [1.003_dp * truth(1:4), 0.98_dp * truth(5), 1.2_dp * truth(6)] / model%scale
    call map_fit(model, theta, 30, steps, message, conv)
    print '(A,I0,A,I0,A)', 'fit: ', steps, ' steps, ', model%solves, ' provider solves'
    call check(conv, 'fit did not converge: '//message)
    call laplace(model, theta, lap)
    call check(lap%proper .or. lap%info == 0, 'Laplace covariance failed')
    do k = 1, size(theta)
        sd = sqrt(lap%cov(k, k)) * model%scale(k)
        err = theta(k) * model%scale(k) - truth(k)
        print '(A,3ES14.6)', model%names(k)(1:16), theta(k) * model%scale(k), truth(k), sd
        call check(abs(err) <= max(2 * sd, 1.0e-12_dp) .or. abs(err) <= 1.0e-6_dp * abs(truth(k)), &
            trim(model%names(k))//' outside two standard deviations')
        call check(abs(err) <= 1.0e-3_dp * abs(truth(k)), trim(model%names(k))//' not recovered')
    end do
    if (failures > 0) then
        write(error_unit, '(I0,A)') failures, ' free-boundary recovery checks failed'
        error stop 1
    end if
    print '(A)', 'axisym free boundary: all checks passed'

contains

    subroutine check(condition, what)
        logical, intent(in) :: condition
        character(*), intent(in) :: what
        if (.not. condition) then
            write(error_unit, '(A)') what
            failures = failures + 1
        end if
    end subroutine check

end program test_axisym_free_boundary
