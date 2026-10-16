program test_mse_plasma
    !! Finite differences of the actual field source, at fixed channel geometry.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    use tiago_diagnostic_types, only: mse_channel_t
    use tiago_vacuum_forward, only: vacuum_solver_t
    use tiago_plasma_support, only: vmec_boundary_t, read_vmec_boundary
    implicit none
    type(vacuum_solver_t) :: solver
    type(vmec_boundary_t) :: vb, shifted
    type(mse_channel_t) :: ch(1)
    character(len=512) :: wout, mode, coils
    character(len=:), allocatable :: message
    real(dp), allocatable :: angle(:), plus(:), minus(:), jac(:, :), shape(:, :)
    real(dp) :: step, fd, err_field, err_shape, point(3), b(3, 1), bc(3, 1)
    real(dp) :: num, den, expected
    integer :: j, col, nm, ns, idx(1), k
    integer(i32) :: ierr

    call get_command_argument(1, wout)
    call get_command_argument(2, coils)
    call get_command_argument(3, mode)
    call solver%init(trim(coils))
    call read_vmec_boundary(trim(wout), vb, trim(mode) /= 'contravariant', &
        trim(mode) == 'conservative')
    call solver%enable_plasma_from_boundary(vb, 16, 16)
    ! Outside the LI383 plasma; field representation is valid here.
    point = [4.0_dp, 0.3_dp, 0.2_dp]
    ch(1)%beam_direction = [1.0_dp, 1.0_dp, 0.4_dp]
    ch(1)%sight_direction = [1.0_dp, -0.2_dp, 0.1_dp]
    ch(1)%beam_point = point - ch(1)%beam_direction
    ch(1)%sight_point = point + ch(1)%sight_direction
    ch(1)%label = 'OUTSIDE'
    call solver%mse(ch, angle, ierr, message, plasma_jacobian=jac, shape_jacobian=shape)
    call success()
    ! Explicit independent electric-vector projection of coil + sheet field.
    call solver%coils%field(reshape(point, [3, 1]), bc)
    call solver%plasma%sample_bfield(reshape(point, [3, 1]), b)
    b = b + bc
    expected = electric_angle(b(:, 1))
    call check('combined-source angle', angle(1), expected, 1.0e-13_dp)
    nm = size(vb%xm_nyq)
    ns = size(vb%xm)
    err_field = 0.0_dp
    err_shape = 0.0_dp
    ! Include the strongest nonzero column in each field/shape family.
    do k = 1, 2
        idx = maxloc(abs(jac(1, (k - 1)*nm + 1:k*nm)))
        j = idx(1)
        col = (k - 1)*nm + j
        step = 1.0e-6_dp
        shifted = vb
        if (k == 1) shifted%bumnc(j) = shifted%bumnc(j) + step
        if (k == 2) shifted%bvmnc(j) = shifted%bvmnc(j) + step
        call eval(shifted, plus)
        shifted = vb
        if (k == 1) shifted%bumnc(j) = shifted%bumnc(j) - step
        if (k == 2) shifted%bvmnc(j) = shifted%bvmnc(j) - step
        call eval(shifted, minus)
        fd = angle_difference(plus(1), minus(1))/(2*step)
        call check('plasma field Jacobian FD', jac(1, col), fd, &
            2.0e-6_dp*max(1.0_dp, abs(fd)))
        err_field = max(err_field, abs(jac(1, col) - fd))

        idx = maxloc(abs(shape(1, (k - 1)*ns + 1:k*ns)))
        j = idx(1)
        col = (k - 1)*ns + j
        shifted = vb
        if (k == 1) shifted%rmnc(j) = shifted%rmnc(j) + step
        if (k == 2) shifted%zmns(j) = shifted%zmns(j) + step
        call eval(shifted, plus)
        shifted = vb
        if (k == 1) shifted%rmnc(j) = shifted%rmnc(j) - step
        if (k == 2) shifted%zmns(j) = shifted%zmns(j) - step
        call eval(shifted, minus)
        fd = angle_difference(plus(1), minus(1))/(2*step)
        call check('plasma shape Jacobian FD', shape(1, col), fd, &
            2.0e-6_dp*max(1.0_dp, abs(fd)))
        err_shape = max(err_shape, abs(shape(1, col) - fd))
    end do
    print '(A,2ES14.5)', 'maximum absolute FD errors (field, shape): ', &
        err_field, err_shape
    call solver%finalize()
    print '(A)', 'test_mse_plasma passed: '//trim(mode)
contains
    subroutine eval(boundary, values)
        type(vmec_boundary_t), intent(in) :: boundary
        real(dp), allocatable, intent(out) :: values(:)
        call solver%enable_plasma_from_boundary(boundary, 16, 16)
        call solver%mse(ch, values, ierr, message)
        call success()
    end subroutine eval

    function electric_angle(field) result(gamma)
        real(dp), intent(in) :: field(3)
        real(dp) :: gamma, n(3), z(3), h(3), e(3)
        n = ch(1)%sight_direction/norm2(ch(1)%sight_direction)
        z = [0.0_dp, 0.0_dp, 1.0_dp] - n(3)*n
        z = z/norm2(z)
        h = cross(n, z)
        e = cross(ch(1)%beam_direction, field)
        num = dot_product(h, e)
        den = dot_product(z, e)
        gamma = 0.5_dp*atan2(2*num*den, den*den - num*num)
    end function electric_angle

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
            a(1)*b(2) - a(2)*b(1)]
    end function cross

    real(dp) function angle_difference(a, b) result(d)
        real(dp), intent(in) :: a, b
        d = 0.5_dp*atan2(sin(2*(a - b)), cos(2*(a - b)))
    end function angle_difference

    subroutine success()
        if (ierr == 0) return
        print '(A)', message
        error stop 'MSE plasma evaluation failed'
    end subroutine success

    subroutine check(name, actual, expected, atol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual, expected, atol
        if (abs(actual - expected) <= atol) return
        print '(A,2ES24.16)', name//': ', actual, expected
        error stop 1
    end subroutine check
end program test_mse_plasma
