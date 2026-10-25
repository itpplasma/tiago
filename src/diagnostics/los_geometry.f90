module tiago_los_geometry
    !! Located sample points of profile diagnostics and their derivatives.
    !!
    !! A sample set holds real-space points x_k with weights w_k, located in a
    !! real_space_state_t. Fixed points (Thomson scattering volumes) only move in
    !! label space when the state changes. Line-of-sight samples additionally
    !! move along the chord: the plasma segments [t_a, t_b] end where rho = 1,
    !! and the Gauss points t_k = t_a + (t_b - t_a) xi_k with weights
    !! (t_b - t_a) omega_k follow these crossings. The derivatives below are
    !! those of this discrete model, so finite differences of the forward model
    !! reproduce them up to the crossing tolerance.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_real_space_state, only: real_space_state_t, located_t, locate_ok, &
        located_jvp, located_vjp
    implicit none
    private

    type, public :: chord_t
        character(len=:), allocatable :: label
        real(dp) :: origin(3) = 0.0_dp
        real(dp) :: direction(3) = [1.0_dp, 0.0_dp, 0.0_dp] !! normalized on use
        real(dp) :: length = 1.0_dp
        integer :: nscan = 400 !! scan points for the plasma crossings
        integer :: ngauss = 16 !! Gauss points per plasma segment
    end type chord_t

    type, public :: crossing_t
        real(dp) :: t = 0.0_dp
        logical :: fixed = .true. !! chord end inside the plasma
        type(located_t) :: loc
        real(dp) :: drho_dt = 0.0_dp
    end type crossing_t

    type, public :: sample_set_t
        !! Sample k belongs to output owner(k). For chords, sample k lies in
        !! segment seg(k) at fraction xi(k) with weight factor omega(k).
        integer :: n = 0
        real(dp), allocatable :: x(:, :), w(:), xi(:), omega(:), drho_dt(:)
        integer, allocatable :: owner(:), seg(:)
        type(located_t), allocatable :: loc(:)
        logical, allocatable :: inside(:)
        !> Plasma segments: crossings seg_a(j), seg_b(j) in cross(:)
        type(crossing_t), allocatable :: cross(:)
        integer, allocatable :: seg_a(:), seg_b(:)
    end type sample_set_t

    public :: gauss_legendre01
    public :: locate_samples
    public :: chord_samples
    public :: samples_jvp
    public :: samples_vjp
    public :: inside_plasma

contains

    pure logical function inside_plasma(loc)
        type(located_t), intent(in) :: loc
        inside_plasma = loc%status == locate_ok .and. loc%u(1) <= 1.0_dp
    end function inside_plasma

    subroutine locate_samples(state, set)
        !! Locate every point of a fixed-point set (x, w, owner set by caller).
        class(real_space_state_t), intent(in) :: state
        type(sample_set_t), intent(inout) :: set
        integer :: k

        if (allocated(set%loc)) deallocate(set%loc, set%inside)
        allocate(set%loc(set%n), set%inside(set%n))
        !$omp parallel do schedule(dynamic, 4)
        do k = 1, set%n
            call state%locate(set%x(:, k), set%loc(k))
            set%inside(k) = inside_plasma(set%loc(k))
        end do
        !$omp end parallel do
    end subroutine locate_samples

    subroutine chord_samples(state, chords, set)
        !! Plasma segments and Gauss samples of every chord.
        class(real_space_state_t), intent(in) :: state
        type(chord_t), intent(in) :: chords(:)
        type(sample_set_t), intent(out) :: set
        type(crossing_t), allocatable :: cross(:)
        integer, allocatable :: seg_a(:), seg_b(:), seg_owner(:)
        type(crossing_t), allocatable :: cc(:)
        real(dp), allocatable :: xg(:), wg(:)
        real(dp) :: d(3), ta, tb
        integer :: c, j, k, m, nseg, ng

        allocate(cross(0), seg_a(0), seg_b(0), seg_owner(0))
        do c = 1, size(chords)
            call chord_crossings(state, chords(c), cc)
            do j = 1, size(cc), 2
                cross = [cross, cc(j), cc(j + 1)]
                seg_a = [seg_a, size(cross) - 1]
                seg_b = [seg_b, size(cross)]
                seg_owner = [seg_owner, c]
            end do
        end do
        nseg = size(seg_a)
        set%n = 0
        do j = 1, nseg
            set%n = set%n + chords(seg_owner(j))%ngauss
        end do
        allocate(set%x(3, set%n), set%w(set%n), set%xi(set%n), set%omega(set%n), &
            set%owner(set%n), set%seg(set%n), set%drho_dt(set%n))
        m = 0
        do j = 1, nseg
            c = seg_owner(j)
            ng = chords(c)%ngauss
            call gauss_legendre01(ng, xg, wg)
            d = chords(c)%direction / norm2(chords(c)%direction)
            ta = cross(seg_a(j))%t
            tb = cross(seg_b(j))%t
            do k = 1, ng
                m = m + 1
                set%x(:, m) = chords(c)%origin + (ta + (tb - ta) * xg(k)) * d
                set%w(m) = (tb - ta) * wg(k)
                set%xi(m) = xg(k)
                set%omega(m) = wg(k)
                set%owner(m) = c
                set%seg(m) = j
            end do
        end do
        call move_alloc(cross, set%cross)
        call move_alloc(seg_a, set%seg_a)
        call move_alloc(seg_b, set%seg_b)
        call locate_samples(state, set)
        do m = 1, set%n
            d = chords(set%owner(m))%direction / norm2(chords(set%owner(m))%direction)
            set%drho_dt(m) = dot_product(set%loc(m)%jinv(1, :), d)
            ! Gauss points of a segment lie inside up to the crossing tolerance.
            set%inside(m) = set%loc(m)%status == locate_ok
        end do
    end subroutine chord_samples

    subroutine chord_crossings(state, chord, cross)
        !! Pairs (entry, exit) of rho = 1 crossings along the chord.
        class(real_space_state_t), intent(in) :: state
        type(chord_t), intent(in) :: chord
        type(crossing_t), allocatable, intent(out) :: cross(:)
        type(located_t), allocatable :: scan(:)
        logical, allocatable :: inside(:)
        real(dp) :: d(3), t
        integer :: j, n

        n = chord%nscan
        d = chord%direction / norm2(chord%direction)
        allocate(scan(0:n), inside(0:n), cross(0))
        !$omp parallel do private(t) schedule(dynamic, 8)
        do j = 0, n
            t = chord%length * real(j, dp) / real(n, dp)
            call state%locate(chord%origin + t * d, scan(j))
            inside(j) = inside_plasma(scan(j))
        end do
        !$omp end parallel do
        if (inside(0)) cross = [cross, end_point(0.0_dp, scan(0))]
        do j = 1, n
            if (inside(j) .neqv. inside(j - 1)) then
                if (inside(j)) then
                    cross = [cross, refine(chord%length * (j - 1) / real(n, dp), &
                        chord%length * j / real(n, dp))]
                else
                    cross = [cross, refine(chord%length * j / real(n, dp), &
                        chord%length * (j - 1) / real(n, dp))]
                end if
            end if
        end do
        if (inside(n)) cross = [cross, end_point(chord%length, scan(n))]
    contains
        function end_point(tt, loc) result(c)
            real(dp), intent(in) :: tt
            type(located_t), intent(in) :: loc
            type(crossing_t) :: c
            c%t = tt
            c%fixed = .true.
            c%loc = loc
            c%drho_dt = dot_product(loc%jinv(1, :), d)
        end function end_point

        function refine(t_out, t_in) result(c)
            !! Bisection for rho = 1 between an outside and an inside point.
            real(dp), intent(in) :: t_out, t_in
            type(crossing_t) :: c
            real(dp) :: a, b, tm
            type(located_t) :: loc
            integer :: it

            a = t_in
            b = t_out
            do it = 1, 200
                tm = 0.5_dp * (a + b)
                if (abs(b - a) <= 4.0_dp * epsilon(1.0_dp) * max(chord%length, abs(tm))) exit
                call state%locate(chord%origin + tm * d, loc)
                if (inside_plasma(loc)) then
                    a = tm
                else
                    b = tm
                end if
            end do
            c%t = a
            c%fixed = .false.
            call state%locate(chord%origin + a * d, c%loc)
            c%drho_dt = dot_product(c%loc%jinv(1, :), d)
        end function refine
    end subroutine chord_crossings

    subroutine samples_jvp(state, set, dp_, drho, dw)
        !! Changes of sample labels rho_k and weights w_k for a state tangent.
        class(real_space_state_t), intent(in) :: state
        type(sample_set_t), intent(in) :: set
        real(dp), intent(in) :: dp_(:)
        real(dp), intent(out) :: drho(:), dw(:)
        real(dp), allocatable :: dt(:)
        real(dp) :: du(3), dbx, dta, dtb
        integer :: k, j

        dw = 0.0_dp
        allocate(dt(0))
        if (allocated(set%cross)) then
            deallocate(dt)
            allocate(dt(size(set%cross)))
            do j = 1, size(set%cross)
                dt(j) = 0.0_dp
                if (set%cross(j)%fixed) cycle
                call located_jvp(state, set%cross(j)%loc, dp_, du, dbx)
                dt(j) = -du(1) / set%cross(j)%drho_dt
            end do
        end if
        !$omp parallel do private(du, dbx, dta, dtb, j)
        do k = 1, set%n
            drho(k) = 0.0_dp
            if (.not. set%inside(k)) cycle
            call located_jvp(state, set%loc(k), dp_, du, dbx)
            drho(k) = du(1)
            if (allocated(set%seg)) then
                j = set%seg(k)
                dta = dt(set%seg_a(j))
                dtb = dt(set%seg_b(j))
                drho(k) = drho(k) + set%drho_dt(k) * &
                    ((1.0_dp - set%xi(k)) * dta + set%xi(k) * dtb)
                dw(k) = (dtb - dta) * set%omega(k)
            end if
        end do
        !$omp end parallel do
    end subroutine samples_jvp

    subroutine samples_vjp(state, set, rhobar, wbar, pbar)
        !! pbar += adjoint of samples_jvp for cotangents of rho_k and w_k.
        class(real_space_state_t), intent(in) :: state
        type(sample_set_t), intent(in) :: set
        real(dp), intent(in) :: rhobar(:), wbar(:)
        real(dp), intent(inout) :: pbar(:)
        real(dp), allocatable :: tbar(:), acc(:)
        real(dp) :: tk
        integer :: k, j, np

        np = size(pbar)
        allocate(acc(np))
        acc = 0.0_dp
        if (allocated(set%cross)) then
            allocate(tbar(size(set%cross)))
        else
            allocate(tbar(0))
        end if
        tbar = 0.0_dp
        !$omp parallel do private(j, tk) reduction(+:acc, tbar)
        do k = 1, set%n
            if (.not. set%inside(k)) cycle
            call located_vjp(state, set%loc(k), [rhobar(k), 0.0_dp, 0.0_dp], 0.0_dp, acc)
            if (allocated(set%seg)) then
                j = set%seg(k)
                tk = rhobar(k) * set%drho_dt(k)
                tbar(set%seg_a(j)) = tbar(set%seg_a(j)) + (1.0_dp - set%xi(k)) * tk &
                    - wbar(k) * set%omega(k)
                tbar(set%seg_b(j)) = tbar(set%seg_b(j)) + set%xi(k) * tk &
                    + wbar(k) * set%omega(k)
            end if
        end do
        !$omp end parallel do
        if (allocated(set%cross)) then
            do j = 1, size(set%cross)
                if (set%cross(j)%fixed) cycle
                call located_vjp(state, set%cross(j)%loc, &
                    [-tbar(j) / set%cross(j)%drho_dt, 0.0_dp, 0.0_dp], 0.0_dp, acc)
            end do
        end if
        pbar = pbar + acc
    end subroutine samples_vjp

    subroutine gauss_legendre01(n, x, w)
        !! Gauss-Legendre nodes and weights on [0, 1].
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: x(:), w(:)
        real(dp) :: z, z1, p1, p2, p3, pp
        real(dp), parameter :: pi = acos(-1.0_dp)
        integer :: i, j, it

        allocate(x(n), w(n))
        do i = 1, (n + 1) / 2
            z = cos(pi * (i - 0.25_dp) / (n + 0.5_dp))
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
end module tiago_los_geometry
