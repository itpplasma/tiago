program test_island_vmec
    !! A VMEC field B = B^theta e_theta + B^zeta e_zeta is tangent to the
    !! interpolated flux surfaces, so a traced line keeps its label rho.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_field_lines, only: cartesian_line_field_t
    use tiago_vmec_state, only: vmec_state_t
    use tiago_real_space_state, only: located_t, locate_ok
    use tiago_island_measure, only: trace_options_t, poincare_section
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    type(cartesian_line_field_t) :: field
    type(vmec_state_t), allocatable :: state
    type(trace_options_t) :: opt
    type(located_t) :: loc
    character(len=512) :: path
    character(len=:), allocatable :: message
    real(dp) :: x(3), dx(3, 3), b, db(3), starts(2, 2), pts(2, 10, 2), drift
    integer :: ierr, st(2), j, k, failures

    failures = 0
    call get_command_argument(1, path)
    allocate (state)
    call state%load(trim(path), ierr, message, vector=.true.)
    if (ierr /= 0) error stop 'cannot load wout'
    do j = 1, 2
        ! Outboard midplane (theta = 0) at rho = 0.4 and 0.8 on phi = 0.
        call state%eval([0.4_dp*j, 0.0_dp, 0.0_dp], x, dx, b, db)
        starts(:, j) = [x(1), x(3)]
    end do
    call move_alloc(state, field%state)
    opt%rtol = 1.0e-10_dp
    opt%atol = 1.0e-12_dp
    opt%max_step = 2*pi/256
    call poincare_section(field, starts, 0.0_dp, 1.0_dp, 10, pts, st, opt)
    drift = 0
    do j = 1, 2
        if (st(j) /= 0) failures = failures + 1
        do k = 1, 10
            call field%state%locate([pts(1, k, j), 0.0_dp, pts(2, k, j)], loc)
            if (loc%status /= locate_ok) failures = failures + 1
            drift = max(drift, abs(loc%u(1) - 0.4_dp*j))
        end do
    end do
    print '(a,es10.2)', 'largest change of rho over 10 turns: ', drift
    if (drift > 1.0e-6_dp) failures = failures + 1
    if (failures > 0) error stop 'VMEC field line left its flux surface'
    print '(a)', 'VMEC field-line checks passed'
end program test_island_vmec
