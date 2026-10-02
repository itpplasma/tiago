program test_inactive_coil_segments
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_coil_kernels, only: coil_set_t, build_coil_set
    implicit none

    type(coil_set_t) :: coil, zero_coil
    real(dp) :: x(4), y(4), z(4), current(4), point(3, 4), a(3, 4), b(3, 4)
    real(dp) :: expected

    x = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    y = 0.0_dp
    z = 0.0_dp
    current = [1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
    point(1, :) = 1.5_dp
    point(2:3, :) = 0.0_dp
    coil = build_coil_set(x, y, z, current)
    call coil%vector_potential(point, a)
    call coil%field(point, b)
    ! At this valid point on the inactive connector, each active wire has
    ! integral_0^1 1/(1.5-t) dt = log(3); both fields vanish on their extensions.
    expected = 2.0e-7_dp * log(3.0_dp)
    if (.not. all(ieee_is_finite(a))) error stop 'inactive connector produced nonfinite A'
    if (.not. all(ieee_is_finite(b))) error stop 'inactive connector produced nonfinite B'
    if (maxval(abs(a(1, :) - expected)) > 1.0e-14_dp * expected) &
        error stop 'potential differs from independent wire integral'
    if (maxval(abs(a(2:3, :))) > 0.0_dp) error stop 'transverse potential differs'
    if (maxval(abs(b)) > 0.0_dp) error stop 'extension field differs'

    ! Vanishing actual current does not remove the geometry needed for a
    ! later unit-current evaluation (the basis used by group responses).
    current = 0.0_dp
    zero_coil = build_coil_set(x, y, z, current)
    call zero_coil%vector_potential(point, a)
    call zero_coil%field(point, b)
    if (.not. all(ieee_is_finite(a))) error stop 'zero-current potential is not finite'
    if (.not. all(ieee_is_finite(b))) error stop 'zero-current field is not finite'
    if (maxval(abs(a)) > 0.0_dp .or. maxval(abs(b)) > 0.0_dp) &
        error stop 'zero current produced nonzero signal'
    zero_coil%ix = coil%ix
    zero_coil%iy = coil%iy
    zero_coil%iz = coil%iz
    zero_coil%active = coil%active
    call zero_coil%vector_potential(point, a)
    if (.not. all(ieee_is_finite(a))) error stop 'unit-current response is not finite'
    if (maxval(abs(a(1, :) - expected)) > 1.0e-14_dp * expected) &
        error stop 'zero actual current lost unit-current geometry'

    point(1, :) = 0.5_dp
    call coil%vector_potential(point, a)
    call coil%field(point, b)
    if (all(ieee_is_finite(a))) error stop 'active-wire potential singularity was hidden'
    if (all(ieee_is_finite(b))) error stop 'active-wire field singularity was hidden'
end program test_inactive_coil_segments
