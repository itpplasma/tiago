program test_axisym_joint
    !! Joint free-boundary and profile posterior on the FreeGS TestTokamakSensor
    !! G-EQDSK with a provider stand-in that returns that file for every
    !! request, so the free-boundary columns vanish. Checks the residual layout,
    !! the pressure rows against the G-EQDSK pressure and the kinetic pressure
    !! of the profiles, and the exact profile columns (diagnostic and pressure
    !! VJPs) against central differences.
    !! Usage: test_axisym_joint <geqdsk> <machine> <stub provider> <workdir>
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_axisym_joint, only: axisym_joint_t
    use tiago_axisym_machine, only: read_axisym_machine
    use tiago_geqdsk_plasma, only: read_geqdsk_plasma
    use tiago_point_probe, only: point_probe_t
    use tiago_profile_diagnostic, only: Q_NE, Q_TE
    use tiago_kinetic_profiles, only: P_NE, P_TE, P_TI, P_ZEFF
    implicit none
    type(axisym_joint_t) :: m
    type(point_probe_t) :: probe
    character(len=1024) :: gpath, mpath, stub, work
    character(len=:), allocatable :: message
    real(dp), allocatable :: theta(:), r(:), rp(:), rm(:), jac(:, :), g(:), s(:), p(:)
    real(dp) :: h, err, scale
    integer :: ierr, i, k, nf, nr, np, failures
    logical :: ok
    real(dp), parameter :: breaks(6) = [0.0_dp, 0.3_dp, 0.6_dp, 0.85_dp, 0.95_dp, 1.0_dp]

    failures = 0
    call get_command_argument(1, gpath)
    call get_command_argument(2, mpath)
    call get_command_argument(3, stub)
    call get_command_argument(4, work)
    m%fb%machine_path = trim(mpath)
    call read_axisym_machine(m%fb%machine_path, m%fb%machine, ierr, message)
    if (ierr /= 0) error stop message
    m%fb%command = 'sh '//trim(stub)//' '//trim(gpath)
    m%fb%workdir = trim(work)
    call execute_command_line('mkdir -p '//trim(work))
    m%fb%coil_current = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    m%fb%plasma_current = 2.0e5_dp
    m%fb%profile = 'paxis'
    m%fb%names = [character(len=32) :: 'PLASMA_CURRENT']
    m%fb%scale = [2.0e5_dp]
    m%fb%prior_mean = [2.0e5_dp]
    m%fb%prior_sigma = [2.0e3_dp]
    allocate(m%fb%data(m%fb%machine%sensors%count()), source=0.0_dp)
    allocate(m%fb%sigma(size(m%fb%data)), source=1.0_dp)

    call read_geqdsk_plasma(trim(gpath), m%pp%eq, ierr, message)
    if (ierr /= 0) error stop message
    call m%pp%prof%init(breaks, 6.0_dp)
    g = m%pp%prof%greville
    call m%pp%prof%set_profile(P_NE, 1.0e19_dp * (1 - 0.6_dp * g**2), .true., &
        sigma_f=0.5_dp, l_core=0.3_dp, l_edge=0.05_dp, s0=0.85_dp, sw=0.05_dp)
    call m%pp%prof%set_profile(P_TE, 3.0e2_dp * (1 - 0.9_dp * g**2), .true., &
        sigma_f=0.5_dp, l_core=0.3_dp, l_edge=0.05_dp, s0=0.85_dp, sw=0.05_dp)
    call m%pp%prof%set_profile(P_TI, 3.0e2_dp * (1 - 0.9_dp * g**2), .true., &
        sigma_f=0.5_dp, l_core=0.3_dp, l_edge=0.05_dp, s0=0.85_dp, sw=0.05_dp)
    call m%pp%prof%set_profile(P_ZEFF, [(1.5_dp, i = 1, size(g))], .false.)
    ! n_e and T_e at 8 points on the outboard midplane of the FreeGS plasma.
    allocate(probe%x0(3, 16), probe%quantity(16))
    do i = 1, 8
        probe%x0(:, i) = [1.20_dp + 0.30_dp * (i - 1) / 7, 0.0_dp, 0.05_dp]
        probe%x0(:, 8 + i) = probe%x0(:, i)
    end do
    probe%quantity = [(Q_NE, i = 1, 8), (Q_TE, i = 1, 8)]
    probe%aligned = .false.
    allocate(m%pp%slots(1))
    allocate(m%pp%slots(1)%d, source=probe)
    m%pp%slots(1)%eta_sigma = [0.01_dp, 0.01_dp]
    m%pp%slots(1)%cal_sigma = 0.05_dp
    allocate(m%pp%slots(1)%data(16), source=1.0_dp)
    allocate(m%pp%slots(1)%sigma(16), source=1.0_dp)
    call m%pp%setup()
    m%p_psin = [0.2_dp, 0.5_dp, 0.8_dp]
    m%p_sigma = [50.0_dp, 50.0_dp, 50.0_dp]

    nf = m%n_fb()
    nr = m%fb%n_residuals()
    np = m%pp%n_residuals()
    allocate(theta(nf + m%pp%n), r(m%n_residuals()), rp(m%n_residuals()), rm(m%n_residuals()))
    theta(1) = 1.0_dp
    do i = 2, size(theta)
        theta(i) = 0.3_dp * sin(1.7_dp * i)
    end do
    call check(m%n_residuals() == size(m%fb%data) + 1 + 16 + m%pp%n + 3, 'residual count')
    call m%residual(theta, r, ok)
    call check(ok, 'residual at the point')
    call m%kinetic_pressure(theta(nf + 1:), m%pp%eq, s, p)
    do k = 1, 3
        call close_to(r(nr + np + k), (m%pp%eq%pressure(m%p_psin(k)) - p(k)) / 50.0_dp, &
            1.0e-12_dp, 'pressure row')
    end do
    call m%jacobian(theta, jac)
    ! profile columns against central differences (fixed equilibrium)
    err = 0
    scale = maxval(abs(jac(:, nf + 1:)))
    do k = nf + 1, size(theta)
        h = 1.0e-6_dp
        theta(k) = theta(k) + h
        call m%residual(theta, rp, ok)
        theta(k) = theta(k) - 2 * h
        call m%residual(theta, rm, ok)
        theta(k) = theta(k) + h
        err = max(err, maxval(abs((rp - rm) / (2 * h) - jac(:, k))) / scale)
    end do
    print '(A,ES10.3)', 'profile columns, largest difference / largest entry: ', err
    call check(err < 1.0e-6_dp, 'profile columns of the Jacobian')
    ! with a fixed G-EQDSK only the I_P prior row depends on theta(1): scale / sigma = 100
    call check(all(jac(:nr - 1, 1) == 0) .and. all(jac(nr + 1:, 1) == 0), &
        'free-boundary column with a fixed G-EQDSK')
    call close_to(jac(nr, 1), 100.0_dp, 1.0e-6_dp, 'prior row of I_P')
    call check(any(jac(nr + np + 1:, nf + 1:) /= 0), 'pressure rows depend on the profiles')
    if (failures > 0) error stop 1
    print '(A)', 'axisym joint: all checks passed'

contains

    subroutine check(condition, what)
        logical, intent(in) :: condition
        character(*), intent(in) :: what
        if (.not. condition) then
            print '(A)', 'FAILED: '//what
            failures = failures + 1
        end if
    end subroutine check

    subroutine close_to(a, b, tol, what)
        real(dp), intent(in) :: a, b, tol
        character(*), intent(in) :: what
        call check(abs(a - b) <= tol * max(1.0_dp, abs(b)), what)
    end subroutine close_to

end program test_axisym_joint
