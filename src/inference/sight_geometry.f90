module tiago_sight_geometry
    !! Small geometry helpers of the profile diagnostics: Gauss-Legendre rules,
    !! Cartesian to (R, phi, Z) points, and composite rules along a chord.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    type, public :: chord_t
        !! Straight chord origin + t d, t in [0, length], Cartesian [m].
        character(len=:), allocatable :: label
        real(dp) :: origin(3) = 0.0_dp
        real(dp) :: direction(3) = [1.0_dp, 0.0_dp, 0.0_dp] !! normalized on use
        real(dp) :: length = 1.0_dp
        integer :: npanel = 200 !! panels of the composite rule
        integer :: ngauss = 4 !! Gauss points per panel
    end type chord_t

    public :: gauss_legendre01, to_rpz, chord_rule

contains

    pure function to_rpz(x) result(p)
        real(dp), intent(in) :: x(3)
        real(dp) :: p(3)
        p = [hypot(x(1), x(2)), atan2(x(2), x(1)), x(3)]
    end function to_rpz

    subroutine chord_rule(chord, points, weights)
        !! Composite Gauss-Legendre points (R, phi, Z) and weights [m].
        type(chord_t), intent(in) :: chord
        real(dp), allocatable, intent(out) :: points(:, :), weights(:)
        real(dp), allocatable :: xg(:), wg(:)
        real(dp) :: d(3), h, t
        integer :: i, k, m

        call gauss_legendre01(chord%ngauss, xg, wg)
        d = chord%direction / norm2(chord%direction)
        h = chord%length / chord%npanel
        allocate(points(3, chord%npanel * chord%ngauss), weights(chord%npanel * chord%ngauss))
        m = 0
        do i = 1, chord%npanel
            do k = 1, chord%ngauss
                m = m + 1
                t = (i - 1 + xg(k)) * h
                points(:, m) = to_rpz(chord%origin + t * d)
                weights(m) = h * wg(k)
            end do
        end do
    end subroutine chord_rule

    pure subroutine gauss_legendre01(n, x, w)
        !! Gauss-Legendre nodes and weights on [0, 1].
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: x(:), w(:)
        real(dp) :: z, z1, p1, p2, p3, pp
        integer :: i, j, it

        allocate(x(n), w(n))
        do i = 1, (n + 1) / 2
            z = cos(pi * (i - 0.25_dp) / (n + 0.5_dp))
            pp = 1.0_dp
            do it = 1, 100
                p1 = 1.0_dp
                p2 = 0.0_dp
                do j = 1, n
                    p3 = p2
                    p2 = p1
                    p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j
                end do
                pp = n * (z * p1 - p2) / (z * z - 1.0_dp)
                z1 = z
                z = z1 - p1 / pp
                if (abs(z - z1) <= 1.0e-15_dp) exit
            end do
            x(i) = 0.5_dp * (1.0_dp - z)
            x(n + 1 - i) = 0.5_dp * (1.0_dp + z)
            w(i) = 1.0_dp / ((1.0_dp - z * z) * pp * pp)
            w(n + 1 - i) = w(i)
        end do
    end subroutine gauss_legendre01
end module tiago_sight_geometry
