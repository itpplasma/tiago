program test_coil_potential
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_coil_kernels, only: coil_set_t, build_coil_set
    implicit none

    type(coil_set_t) :: coil
    real(dp), parameter :: lengths(3) = [1.0e-12_dp, 1.0_dp, 1.0e12_dp]
    real(dp), parameter :: axial(3) = [-3.0_dp, 0.37_dp, 3.0_dp]
    real(dp), parameter :: radial(5) = [1.0e-4_dp, 0.1_dp, 1.0_dp, 1.0e4_dp, 1.0e12_dp]
    real(dp), parameter :: ratios(3) = 0.01_dp*[1.0_dp - 1.0e-8_dp, 1.0_dp, &
        1.0_dp + 1.0e-8_dp]
    real(dp) :: x(4), y(4), z(4), current(4)
    real(dp) :: point(3, 5), potential(3, 5), expected(5)
    real(dp) :: threshold_point(3, 3), threshold_potential(3, 3)
    real(dp) :: length, rho, u, tolerance
    integer :: l, r, a, polarity

    ! Integrate 1/sqrt((u-s)**2+rho**2) directly along the wire. Its asinh
    ! primitive is independent of the kernel's endpoint-distance/log formula.
    ! A repeated node contributes no current. Five points exercise worker batches.
    do polarity = -1, 1, 2
        do l = 1, size(lengths)
            length = lengths(l)
            x(1) = 0.0_dp
            x(2) = 0.5_dp*length
            x(3) = x(2)
            x(4) = length
            y = 0.0_dp
            z = 0.0_dp
            current(1) = 2.0_dp*real(polarity, dp)
            current(2) = 99.0_dp*real(polarity, dp)
            current(3) = current(1)
            current(4) = 0.0_dp
            coil = build_coil_set(x, y, z, current)
            if (coil%n_segments() /= 2) error stop 'repeated node became active'
            do a = 1, size(axial)
                u = axial(a)*length
                do r = 1, size(radial)
                    rho = radial(r)*length
                    point(1, r) = u
                    point(2, r) = rho
                    point(3, r) = 0.0_dp
                    expected(r) = 2.0e-7_dp*real(polarity, dp) &
                        *(asinh(u/rho) - asinh((u - length)/rho))
                end do
                call coil%vector_potential(point, potential)
                do r = 1, size(radial)
                    tolerance = 1.0e-8_dp
                    if (radial(r) >= 1.0e4_dp) tolerance = 1.0e-12_dp
                    call check_potential(potential(:, r), expected(r), tolerance)
                end do
            end do
        end do
    end do

    ! Public current and geometry arrays remain effective on the next call.
    coil%ix = -0.5_dp*coil%ix
    expected = -0.5_dp*expected
    call coil%vector_potential(point, potential)
    do r = 1, size(radial)
        call check_potential(potential(:, r), expected(r), 1.0e-8_dp)
    end do
    coil%xn = 2.0_dp*coil%xn
    coil%yn = 2.0_dp*coil%yn
    coil%zn = 2.0_dp*coil%zn
    coil%ix = 2.0_dp*coil%ix
    coil%iy = 2.0_dp*coil%iy
    coil%iz = 2.0_dp*coil%iz
    coil%len = 2.0_dp*coil%len
    point = 2.0_dp*point
    call coil%vector_potential(point, potential)
    do r = 1, size(radial)
        call check_potential(potential(:, r), expected(r), 1.0e-8_dp)
    end do

    ! For a centred observation point, eps = L/(2*sqrt((L/2)**2+rho**2)).
    ! Place samples on both sides of the polynomial/log boundary at every scale.
    do l = 1, size(lengths)
        length = lengths(l)
        x(1) = 0.0_dp
        x(2) = length
        y = 0.0_dp
        z = 0.0_dp
        current(1) = 2.0_dp
        current(2) = 0.0_dp
        coil = build_coil_set(x(:2), y(:2), z(:2), current(:2))
        do r = 1, size(ratios)
            rho = 0.5_dp*length*sqrt(1.0_dp/ratios(r)**2 - 1.0_dp)
            threshold_point(1, r) = 0.5_dp*length
            threshold_point(2, r) = rho
            threshold_point(3, r) = 0.0_dp
            expected(r) = 4.0e-7_dp*asinh(0.5_dp*length/rho)
        end do
        call coil%vector_potential(threshold_point, threshold_potential)
        do r = 1, size(ratios)
            call check_potential(threshold_potential(:, r), expected(r), 1.0e-12_dp)
        end do
    end do

contains

    subroutine check_potential(actual, target, relative_tolerance)
        real(dp), intent(in) :: actual(3), target, relative_tolerance

        if (.not. all(ieee_is_finite(actual))) &
            error stop 'straight-wire potential is not finite'
        if (.not. ieee_is_finite(target)) error stop 'asinh oracle is not finite'
        if (abs(actual(1) - target) > relative_tolerance*abs(target)) &
            error stop 'straight-wire potential differs from asinh integral'
        if (maxval(abs(actual(2:3))) > 0.0_dp) &
            error stop 'straight-wire potential has transverse components'
    end subroutine check_potential
end program test_coil_potential
