module tiago_radial_interp
    !! C1 cubic Hermite interpolation on a non-uniform node grid, written as
    !! weights on the node values: f(r) = sum_k w(k) f(first + k - 1). Node
    !! slopes are three-point (one-sided second-order at the ends) differences,
    !! so the interpolant reproduces quadratics exactly and its derivatives with
    !! respect to the node values are the weights themselves.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_hermite_basis, only: hermite_basis
    use tiago_hermite_higher, only: hermite_higher
    implicit none
    private

    integer, parameter, public :: max_weights = 6

    public :: hermite_weights

contains

    pure subroutine hermite_weights(nodes, r, extrapolate, first, nw, w, dw, d2w, d3w)
        !! Weights w and d/dr weights dw on nodes first .. first + nw - 1.
        !! Outside the node range the end cubic is extended (extrapolate) or
        !! the end value is held constant (dw = 0).
        real(dp), intent(in) :: nodes(:), r
        logical, intent(in) :: extrapolate
        integer, intent(out) :: first, nw
        real(dp), intent(out) :: w(max_weights), dw(max_weights)
        real(dp), intent(out), optional :: d2w(max_weights), d3w(max_weights)
        real(dp) :: w2(max_weights), w3(max_weights), h2(4), h3(4)
        integer :: n, i, k, j
        real(dp) :: h, t, rr, h00, h10, h01, h11, d00, d10, d01, d11
        real(dp) :: s(3)
        integer :: sidx(3)

        w2 = 0.0_dp
        w3 = 0.0_dp
        if (present(d2w)) d2w = 0.0_dp
        if (present(d3w)) d3w = 0.0_dp
        n = size(nodes)
        w = 0.0_dp
        dw = 0.0_dp
        rr = r
        if (.not. extrapolate) rr = min(max(r, nodes(1)), nodes(n))
        if (n == 1) then
            first = 1
            nw = 1
            w(1) = 1.0_dp
            return
        end if
        i = 1
        do k = 1, n - 1
            if (rr >= nodes(k)) i = k
        end do
        first = max(1, i - 2)
        nw = min(n, i + 3) - first + 1
        h = nodes(i + 1) - nodes(i)
        t = (rr - nodes(i)) / h
        call hermite_basis(t, h00, h10, h01, h11, d00, d10, d01, d11)
        call hermite_higher(t, h2, h3)
        w2(i - first + 1) = h2(1) / h**2
        w3(i - first + 1) = h3(1) / h**3
        w2(i - first + 2) = h2(3) / h**2
        w3(i - first + 2) = h3(3) / h**3
        call add(w, dw, i - first + 1, h00, d00 / h)
        call add(w, dw, i + 1 - first + 1, h01, d01 / h)
        call slope(i, sidx, s)
        do j = 1, 3
            if (sidx(j) <= 0) cycle
            call add(w, dw, sidx(j) - first + 1, h * h10 * s(j), d10 * s(j))
            w2(sidx(j) - first + 1) = w2(sidx(j) - first + 1) + h2(2)*s(j)/h
            w3(sidx(j) - first + 1) = w3(sidx(j) - first + 1) + h3(2)*s(j)/h**2
        end do
        call slope(i + 1, sidx, s)
        do j = 1, 3
            if (sidx(j) <= 0) cycle
            call add(w, dw, sidx(j) - first + 1, h * h11 * s(j), d11 * s(j))
            w2(sidx(j) - first + 1) = w2(sidx(j) - first + 1) + h2(4)*s(j)/h
            w3(sidx(j) - first + 1) = w3(sidx(j) - first + 1) + h3(4)*s(j)/h**2
        end do
        if (.not. extrapolate .and. (r < nodes(1) .or. r > nodes(n))) then
            dw = 0.0_dp
            w2 = 0.0_dp
            w3 = 0.0_dp
        end if
        if (present(d2w)) d2w = w2
        if (present(d3w)) d3w = w3
    contains
        pure subroutine slope(j0, idx, c)
            !! Node slope as c . f(idx).
            integer, intent(in) :: j0
            integer, intent(out) :: idx(3)
            real(dp), intent(out) :: c(3)
            real(dp) :: hl, hr

            idx = 0
            c = 0.0_dp
            if (n == 2) then
                idx(1:2) = [1, 2]
                c(1:2) = [-1.0_dp, 1.0_dp] / (nodes(2) - nodes(1))
            else if (j0 == 1) then
                hl = nodes(2) - nodes(1)
                hr = nodes(3) - nodes(2)
                idx = [1, 2, 3]
                c = [-(2.0_dp * hl + hr) / (hl * (hl + hr)), (hl + hr) / (hl * hr), &
                    -hl / (hr * (hl + hr))]
            else if (j0 == n) then
                hl = nodes(n - 1) - nodes(n - 2)
                hr = nodes(n) - nodes(n - 1)
                idx = [n - 2, n - 1, n]
                c = [hr / (hl * (hl + hr)), -(hl + hr) / (hl * hr), &
                    (2.0_dp * hr + hl) / (hr * (hl + hr))]
            else
                hl = nodes(j0) - nodes(j0 - 1)
                hr = nodes(j0 + 1) - nodes(j0)
                idx = [j0 - 1, j0, j0 + 1]
                c = [-hr / (hl * (hl + hr)), (hr - hl) / (hl * hr), hl / (hr * (hl + hr))]
            end if
        end subroutine slope
    end subroutine hermite_weights

    pure subroutine add(w, dw, k, wv, dv)
        real(dp), intent(inout) :: w(:), dw(:)
        integer, intent(in) :: k
        real(dp), intent(in) :: wv, dv
        w(k) = w(k) + wv
        dw(k) = dw(k) + dv
    end subroutine add
end module tiago_radial_interp
