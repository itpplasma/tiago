program test_plasma_projection
    !! Closed analytic one-form on a circular torus, including asymmetric modes.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_plasma_support, only: plasma_support_t, vmec_boundary_t
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp), major = 3.0_dp, minor = 0.7_dp
    integer, parameter :: nt = 16, np = 8
    type(vmec_boundary_t) :: vb, pert
    type(plasma_support_t) :: plasma, shifted
    real(dp) :: points(3, 2), dl(3, 2), avec(3, 2), plus(2), minus(2), fd(2)
    real(dp) :: theta, phi, r, bu, bv, xt(3), xp(3), expected(3), scale, step
    real(dp), allocatable :: w(:, :, :), gx(:, :, :), resp(:, :), shape(:, :)
    integer :: owner(2), k, col, mode, block_index, failures

    failures = 0
    vb%nfp = 3
    vb%lasym = .true.
    vb%covariant = .true.
    vb%conservative = .true.
    vb%xm = [0.0_dp, 1.0_dp]
    vb%xn = [0.0_dp, 0.0_dp]
    vb%rmnc = [major, minor]
    vb%zmns = [0.0_dp, minor]
    vb%rmns = [0.0_dp, 0.0_dp]
    vb%zmnc = [0.0_dp, 0.0_dp]
    vb%xm_nyq = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    vb%xn_nyq = [0.0_dp, 3.0_dp, 3.0_dp, -3.0_dp]
    vb%bumnc = [0.02_dp, 0.001_dp, 0.002_dp, 0.003_dp]
    vb%bvmnc = [0.3_dp, 0.004_dp, 0.0_dp, 0.0_dp]
    vb%bumns = [0.0_dp, 0.001_dp, 0.0_dp, 0.002_dp]
    vb%bvmns = [0.0_dp, 0.0_dp, 0.001_dp, 0.0_dp]
    call plasma%init_from_boundary(vb, np, nt)

    ! These coefficients are the derivatives of a scalar potential plus the
    ! two constant circulation forms. Check the physical sheet at every point,
    ! independently of the implementation's projection matrix and diagnostics.
    scale = -(2.0_dp * pi / nt) * (2.0_dp * pi / (np * vb%nfp)) / (4.0_dp * pi)
    do k = 1, size(plasma%sheet, 2)
        theta = plasma%theta(k)
        phi = plasma%phi(k)
        r = major + minor * cos(theta)
        xt = [-minor * sin(theta) * cos(phi), &
            -minor * sin(theta) * sin(phi), minor * cos(theta)]
        xp = [-r * sin(phi), r * cos(phi), 0.0_dp]
        bu = 0.02_dp + 0.0002_dp * cos(theta - 3.0_dp * phi) &
            + 0.0003_dp * cos(theta + 3.0_dp * phi) &
            - 0.0003_dp * sin(theta - 3.0_dp * phi) &
            + 0.0002_dp * sin(theta + 3.0_dp * phi)
        bv = 0.3_dp + 0.004_dp * cos(3.0_dp * phi) &
            - 0.0006_dp * cos(theta - 3.0_dp * phi) &
            + 0.0009_dp * cos(theta + 3.0_dp * phi) &
            + 0.0009_dp * sin(theta - 3.0_dp * phi) &
            + 0.0006_dp * sin(theta + 3.0_dp * phi)
        expected = scale * (bu * xp - bv * xt)
        if (maxval(abs(plasma%sheet(:, k) - expected)) > 2.0e-16_dp) then
            failures = failures + 1
        end if
    end do
    if (.not. all(ieee_is_finite(plasma%sheet))) failures = failures + 1
    if (.not. ieee_is_finite(plasma%curl_norm)) failures = failures + 1
    if (.not. ieee_is_finite(plasma%current_ripple)) failures = failures + 1
    if (plasma%curl_norm > 1.0e-16_dp) failures = failures + 1
    if (plasma%current_ripple > 1.0e-12_dp) failures = failures + 1

    points(:, 1) = [4.5_dp, 0.2_dp, 0.3_dp]
    points(:, 2) = [0.3_dp, 4.5_dp, -0.2_dp]
    dl(:, 1) = [0.3_dp, 0.7_dp, 0.2_dp]
    dl(:, 2) = [0.4_dp, -0.2_dp, 0.1_dp]
    owner = [1, 2]
    call plasma%potential_weights(points, dl, owner, 2, w, gx)
    call plasma%mode_response(w, resp)
    call plasma%shape_response(w, gx, shape)
    call plasma%sample_vector_potential(points, avec)
    call check('raw-coefficient response after projection', &
        matmul(resp, plasma%coefficients), sum(avec * dl, dim=1), 1.0e-14_dp)

    ! Central differences perturb the original, unprojected input coefficients.
    step = 1.0e-6_dp
    do col = 1, size(resp, 2)
        mode = mod(col - 1, 4) + 1
        block_index = (col - 1) / 4 + 1
        pert = vb
        call shift_field(pert, block_index, mode, step)
        call signals(pert, plus)
        pert = vb
        call shift_field(pert, block_index, mode, -step)
        call signals(pert, minus)
        fd = (plus - minus) / (2.0_dp * step)
        call check('projected field derivative vs FD', resp(:, col), fd, &
            1.0e-8_dp * maxval(abs(resp)))
    end do
    do col = 1, size(shape, 2)
        mode = mod(col - 1, 2) + 1
        block_index = (col - 1) / 2 + 1
        pert = vb
        call shift_shape(pert, block_index, mode, step)
        call signals(pert, plus)
        pert = vb
        call shift_shape(pert, block_index, mode, -step)
        call signals(pert, minus)
        fd = (plus - minus) / (2.0_dp * step)
        call check('projected shape derivative vs FD', shape(:, col), fd, &
            1.0e-8_dp * maxval(abs(shape)))
    end do
    if (failures /= 0) then
        print '(A,I0)', 'plasma projection failures: ', failures
        error stop 1
    end if
    print '(A)', 'analytic closed sheet, field and shape derivatives passed'

contains

    subroutine signals(boundary, values)
        type(vmec_boundary_t), intent(in) :: boundary
        real(dp), intent(out) :: values(2)
        real(dp) :: a(3, 2)
        call shifted%init_from_boundary(boundary, np, nt)
        call shifted%sample_vector_potential(points, a)
        values = sum(a * dl, dim=1)
    end subroutine signals

    subroutine shift_field(boundary, which, mode, delta)
        type(vmec_boundary_t), intent(inout) :: boundary
        integer, intent(in) :: which, mode
        real(dp), intent(in) :: delta
        select case (which)
        case (1)
            boundary%bumnc(mode) = boundary%bumnc(mode) + delta
        case (2)
            boundary%bvmnc(mode) = boundary%bvmnc(mode) + delta
        case (3)
            boundary%bumns(mode) = boundary%bumns(mode) + delta
        case (4)
            boundary%bvmns(mode) = boundary%bvmns(mode) + delta
        end select
    end subroutine shift_field

    subroutine shift_shape(boundary, which, mode, delta)
        type(vmec_boundary_t), intent(inout) :: boundary
        integer, intent(in) :: which, mode
        real(dp), intent(in) :: delta
        select case (which)
        case (1)
            boundary%rmnc(mode) = boundary%rmnc(mode) + delta
        case (2)
            boundary%zmns(mode) = boundary%zmns(mode) + delta
        case (3)
            boundary%rmns(mode) = boundary%rmns(mode) + delta
        case (4)
            boundary%zmnc(mode) = boundary%zmnc(mode) + delta
        end select
    end subroutine shift_shape

    subroutine check(label, actual, reference, tolerance)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual(:), reference(:), tolerance
        if (.not. all(ieee_is_finite(actual)) &
            .or. .not. all(ieee_is_finite(reference)) &
            .or. .not. ieee_is_finite(tolerance)) then
            print '(A)', 'FAIL nonfinite '//label
            failures = failures + 1
        else if (maxval(abs(actual - reference)) > tolerance) then
            print '(A,ES12.4)', 'FAIL '//label//': ', maxval(abs(actual - reference))
            failures = failures + 1
        end if
    end subroutine check
end program test_plasma_projection
