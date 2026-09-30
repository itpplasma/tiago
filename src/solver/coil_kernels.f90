module tiago_coil_kernels
    !! Vector potential and field of piecewise-straight current filaments
    !! (Hanson & Hirshman, Phys. Plasmas 9, 4410 (2002)), in SI units.
    !!
    !! The filaments are consecutive nodes x_k with current I_k from node k to
    !! k+1 (zero current marks jumpers between coils). For a field point x and
    !! segment k with R_i = |x - x_k|, R_f = |x - x_{k+1}|, L = |x_{k+1} - x_k|,
    !! eps = L / (R_i + R_f):
    !!     A = mu0/(4 pi) I ehat log((1 + eps)/(1 - eps))
    !!     B = mu0/(4 pi) I (dl x (x - x_k)) 2 (R_i + R_f) / (R_i R_f ((R_i + R_f)^2 - L^2))
    !! Distances to the nodes are computed once per field point and shared by the
    !! two segments meeting at a node. Segment k always joins nodes k and k+1, so
    !! the segment loop is contiguous and vectorises; segments without current
    !! (jumpers, repeated points) are stored with I dl = 0 and contribute exactly
    !! zero. Their true length (at least tiny) keeps eps < 1 and all terms finite.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    real(dp), parameter :: mu0_over_4pi = 1.0e-7_dp

    type, public :: coil_set_t
        real(dp), allocatable :: xn(:), yn(:), zn(:)   !! nodes [m]
        real(dp), allocatable :: ix(:), iy(:), iz(:)   !! I * dl of segment k -> k+1 [A m]
        real(dp), allocatable :: len(:)                !! max(|dl|, tiny) [m]
        integer :: active = 0                          !! number of segments with current
    contains
        procedure :: vector_potential => coil_set_vector_potential
        procedure :: field => coil_set_field
        procedure :: n_segments => coil_set_n_segments
    end type coil_set_t

    public :: build_coil_set

contains

    function build_coil_set(x, y, z, current) result(cs)
        real(dp), intent(in) :: x(:), y(:), z(:), current(:)
        type(coil_set_t) :: cs
        real(dp) :: dx, dy, dz, l
        integer :: k, n

        n = max(size(x) - 1, 0)
        cs%xn = x
        cs%yn = y
        cs%zn = z
        allocate(cs%ix(n), cs%iy(n), cs%iz(n), cs%len(n))
        cs%ix = 0.0_dp
        cs%iy = 0.0_dp
        cs%iz = 0.0_dp
        do k = 1, n
            dx = x(k + 1) - x(k)
            dy = y(k + 1) - y(k)
            dz = z(k + 1) - z(k)
            l = sqrt(dx * dx + dy * dy + dz * dz)
            cs%len(k) = max(l, tiny(1.0_dp))
            if (current(k) == 0.0_dp .or. l == 0.0_dp) cycle
            cs%ix(k) = current(k) * dx
            cs%iy(k) = current(k) * dy
            cs%iz(k) = current(k) * dz
            cs%active = cs%active + 1
        end do
    end function build_coil_set

    integer function coil_set_n_segments(self)
        class(coil_set_t), intent(in) :: self
        coil_set_n_segments = self%active
    end function coil_set_n_segments

    subroutine coil_set_vector_potential(self, points, a)
        !! points(3, n) [m] -> a(3, n) [T m]
        class(coil_set_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: a(:, :)

        a = 0.0_dp
        if (self%n_segments() == 0) return
        call potential_kernel(size(self%xn), self%xn, self%yn, self%zn, self%ix, self%iy, &
            self%iz, self%len, size(points, 2), points, a)
    end subroutine coil_set_vector_potential

    subroutine coil_set_field(self, points, b)
        !! points(3, n) [m] -> b(3, n) [T]
        class(coil_set_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: b(:, :)

        b = 0.0_dp
        if (self%n_segments() == 0) return
        call field_kernel(size(self%xn), self%xn, self%yn, self%zn, self%ix, self%iy, &
            self%iz, self%len, size(points, 2), points, b)
    end subroutine coil_set_field

    ! Explicit-shape arrays keep the inner loops contiguous and vectorisable.

    subroutine potential_kernel(nn, xn, yn, zn, ix, iy, iz, len, np, points, a)
        integer, intent(in) :: nn, np
        real(dp), intent(in) :: xn(nn), yn(nn), zn(nn)
        real(dp), intent(in) :: ix(nn - 1), iy(nn - 1), iz(nn - 1), len(nn - 1)
        real(dp), intent(in) :: points(3, np)
        real(dp), intent(out) :: a(3, np)
        real(dp) :: r(nn), ax, ay, az, eps, f
        integer :: i, s

!$omp parallel do default(shared) private(r, i, s, ax, ay, az, eps, f) schedule(static)
        do i = 1, np
            do s = 1, nn
                r(s) = sqrt((points(1, i) - xn(s))**2 + (points(2, i) - yn(s))**2 + &
                    (points(3, i) - zn(s))**2)
            end do
            ax = 0.0_dp
            ay = 0.0_dp
            az = 0.0_dp
            !$omp simd private(eps, f) reduction(+:ax, ay, az)
            do s = 1, nn - 1
                eps = len(s) / (r(s) + r(s + 1))
                f = log((1.0_dp + eps) / (1.0_dp - eps)) / len(s)
                ax = ax + ix(s) * f
                ay = ay + iy(s) * f
                az = az + iz(s) * f
            end do
            a(:, i) = mu0_over_4pi * [ax, ay, az]
        end do
!$omp end parallel do
    end subroutine potential_kernel

    subroutine field_kernel(nn, xn, yn, zn, ix, iy, iz, len, np, points, b)
        integer, intent(in) :: nn, np
        real(dp), intent(in) :: xn(nn), yn(nn), zn(nn)
        real(dp), intent(in) :: ix(nn - 1), iy(nn - 1), iz(nn - 1), len(nn - 1)
        real(dp), intent(in) :: points(3, np)
        real(dp), intent(out) :: b(3, np)
        real(dp) :: r(nn), bx, by, bz, rx, ry, rz, rs, f
        integer :: i, s

!$omp parallel do default(shared) private(r, i, s, bx, by, bz, rx, ry, rz, rs, f) &
!$omp schedule(static)
        do i = 1, np
            do s = 1, nn
                r(s) = sqrt((points(1, i) - xn(s))**2 + (points(2, i) - yn(s))**2 + &
                    (points(3, i) - zn(s))**2)
            end do
            bx = 0.0_dp
            by = 0.0_dp
            bz = 0.0_dp
            !$omp simd private(rx, ry, rz, rs, f) reduction(+:bx, by, bz)
            do s = 1, nn - 1
                rx = points(1, i) - xn(s)
                ry = points(2, i) - yn(s)
                rz = points(3, i) - zn(s)
                rs = r(s) + r(s + 1)
                f = 2.0_dp * rs / (r(s) * r(s + 1) * (rs * rs - len(s)**2))
                bx = bx + (iy(s) * rz - iz(s) * ry) * f
                by = by + (iz(s) * rx - ix(s) * rz) * f
                bz = bz + (ix(s) * ry - iy(s) * rx) * f
            end do
            b(:, i) = mu0_over_4pi * [bx, by, bz]
        end do
!$omp end parallel do
    end subroutine field_kernel
end module tiago_coil_kernels
