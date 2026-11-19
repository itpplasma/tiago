module tiago_profile_reflectometry
    !! Ordinary-mode profile reflectometry at normal incidence along a sight
    !! line: cutoff, round-trip phase and group delay with profile derivatives,
    !! and the Bottollier-Curtet group-delay inversion.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_DOMAIN_ERROR
    use fortnum_quadrature, only: gauss_legendre
    use tiago_omode_index, only: omode_index
    implicit none
    private
    real(dp), parameter :: pi = acos(-1.0_dp), light = 299792458.0_dp

    type, abstract, public :: sight_density_t
        !! Electron density in m^-3 along the sight line, s in metres from the
        !! plasma edge, where the density vanishes, towards the plasma.
    contains
        procedure(count_i), deferred :: n_params
        procedure(sample_i), deferred :: sample
    end type
    abstract interface
        pure integer function count_i(self)
            import :: sight_density_t
            class(sight_density_t), intent(in) :: self
        end function
        subroutine sample_i(self, s, ne, dne_ds, ne_p, dne_ds_p)
            import :: sight_density_t, dp
            class(sight_density_t), intent(in) :: self
            real(dp), intent(in) :: s
            real(dp), intent(out) :: ne, dne_ds
            real(dp), intent(out), optional :: ne_p(:), dne_ds_p(:)
        end subroutine
    end interface

    type, public :: omode_reflection_t
        real(dp) :: cutoff = 0.0_dp
        !! Round-trip phase -pi/2 + (4 pi f/c)(d + int N ds) in rad.
        real(dp) :: phase = 0.0_dp
        !! Round-trip group delay (2/c)(d + int ds/N) in s.
        real(dp) :: group_delay = 0.0_dp
        real(dp), allocatable :: cutoff_p(:), phase_p(:), group_delay_p(:)
    end type

    public :: omode_reflection, invert_omode_group_delay
contains

    subroutine omode_reflection(profile, frequency, vacuum_distance, range, &
            result, status, nodes, search_points)
        class(sight_density_t), intent(in) :: profile
        real(dp), intent(in) :: frequency, vacuum_distance, range
        type(omode_reflection_t), intent(out) :: result
        type(fortnum_status_t), intent(out) :: status
        integer, intent(in), optional :: nodes, search_points
        integer :: m, np, i
        real(dp) :: critical, sc, slope, upper, u, s, weight
        real(dp), allocatable :: x(:), w(:), ne_p(:), dne_p(:), cut_p(:)
        real(dp) :: ne, dne, xr, n_index, dn_dne, dn_df, inv, dinv, dummy
        np = profile%n_params()
        m = 128
        if (present(nodes)) m = nodes
        allocate (result%cutoff_p(np), result%phase_p(np), result%group_delay_p(np))
        result%cutoff_p = 0.0_dp
        result%phase_p = 0.0_dp
        result%group_delay_p = 0.0_dp
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid reflectometry input')
        if (.not. (frequency > 0.0_dp .and. range > 0.0_dp)) return
        if (.not. (vacuum_distance >= 0.0_dp) .or. m < 2) return
        call omode_index(0.0_dp, frequency, xr, n_index, dn_dne, dn_df, inv, dinv, &
            critical)
        call find_cutoff(profile, critical, range, sc, slope, status, search_points)
        if (status%code /= 0) return
        allocate (x(m), w(m), ne_p(np), dne_p(np), cut_p(np))
        call profile%sample(sc, ne, dne, ne_p, dne_p)
        cut_p = -ne_p/slope
        result%cutoff = sc
        result%cutoff_p = cut_p
        ! s = sc - u^2 removes the inverse square-root singularity at the cutoff.
        upper = sqrt(sc)
        call gauss_legendre(m, x, w)
        result%phase = vacuum_distance
        result%group_delay = vacuum_distance
        do i = 1, m
            u = 0.5_dp*upper*(x(i) + 1)
            weight = 0.5_dp*upper*w(i)*2*u
            s = sc - u*u
            call profile%sample(s, ne, dne, ne_p, dne_p)
            call omode_index(ne, frequency, xr, n_index, dn_dne, dn_df, inv, dinv, &
                dummy)
            if (.not. (xr < 1.0_dp .and. ieee_is_finite(inv))) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    'density exceeds the cutoff before the located cutoff')
                return
            end if
            result%phase = result%phase + weight*n_index
            result%group_delay = result%group_delay + weight*inv
            result%phase_p = result%phase_p + weight*dn_dne*ne_p
            result%group_delay_p = result%group_delay_p &
                + weight*dinv*(ne_p + dne*cut_p)
        end do
        ! Edge term of the moving upper limit u = sqrt(sc).
        call profile%sample(0.0_dp, ne, dne)
        call omode_index(ne, frequency, xr, n_index, dn_dne, dn_df, inv, dinv, dummy)
        result%group_delay_p = result%group_delay_p + inv*cut_p
        result%phase = 4*pi*frequency/light*result%phase - 0.5_dp*pi
        result%phase_p = 4*pi*frequency/light*result%phase_p
        result%group_delay = 2/light*result%group_delay
        result%group_delay_p = 2/light*result%group_delay_p
        call status_set(status, 0, '')
    end subroutine

    subroutine find_cutoff(profile, critical, range, sc, slope, status, &
            search_points)
        class(sight_density_t), intent(in) :: profile
        real(dp), intent(in) :: critical, range
        real(dp), intent(out) :: sc, slope
        type(fortnum_status_t), intent(out) :: status
        integer, intent(in), optional :: search_points
        integer :: n, k, iteration
        real(dp) :: lo, hi, ne, dne, step, s, f_lo
        n = 400
        if (present(search_points)) n = search_points
        sc = 0.0_dp
        slope = 0.0_dp
        call status_set(status, FORTNUM_DOMAIN_ERROR, &
            'no ordinary-mode cutoff within the sight range')
        call profile%sample(0.0_dp, ne, dne)
        f_lo = ne - critical
        if (f_lo >= 0.0_dp) return
        lo = 0.0_dp
        hi = -1.0_dp
        do k = 1, n
            s = range*k/n
            call profile%sample(s, ne, dne)
            if (ne - critical >= 0.0_dp) then
                hi = s
                exit
            end if
            lo = s
        end do
        if (hi < 0.0_dp) return
        s = 0.5_dp*(lo + hi)
        do iteration = 1, 100
            call profile%sample(s, ne, dne)
            if (ne < critical) lo = s
            if (ne >= critical) hi = s
            step = 0.5_dp*(lo + hi) - s
            if (dne > 0.0_dp) step = -(ne - critical)/dne
            if (s + step <= lo .or. s + step >= hi) step = 0.5_dp*(lo + hi) - s
            s = s + step
            if (abs(step) <= 4*epsilon(1.0_dp)*max(s, range)) exit
        end do
        call profile%sample(s, ne, dne)
        if (.not. (dne > 0.0_dp)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'density gradient must be positive at the cutoff')
            return
        end if
        sc = s
        slope = dne
        call status_set(status, 0, '')
    end subroutine

    subroutine invert_omode_group_delay(frequency, group_delay, first_cutoff, &
            cutoff, status)
        !! Cutoff positions from plasma group delays (vacuum delay removed) at
        !! ascending frequencies, recursion of Bottollier-Curtet as written by
        !! Morales et al. (2017): trapezoids over reconstructed steps and the
        !! weight 2/3 in the last step. The density vanishes at s = 0.
        real(dp), intent(in) :: frequency(:), group_delay(:), first_cutoff
        real(dp), intent(out) :: cutoff(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), parameter :: last_weight = 2.0_dp/3.0_dp
        real(dp) :: phase, area, n_lo, n_hi, n_last
        integer :: n, i, m
        m = size(frequency)
        cutoff = 0.0_dp
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid inversion input')
        if (size(group_delay) /= m .or. size(cutoff) /= m .or. m < 1) return
        if (.not. (first_cutoff > 0.0_dp .and. frequency(1) > 0.0_dp)) return
        do i = 2, m
            if (.not. (frequency(i) > frequency(i - 1))) return
        end do
        cutoff(1) = first_cutoff
        phase = 4*pi*frequency(1)/light*last_weight*first_cutoff
        do n = 2, m
            phase = phase + pi*(group_delay(n - 1) + group_delay(n)) &
                *(frequency(n) - frequency(n - 1))
            area = 0.0_dp
            n_lo = 1.0_dp
            do i = 1, n - 1
                n_hi = sqrt(1 - (frequency(i)/frequency(n))**2)
                area = area + 0.5_dp*(n_lo + n_hi)* &
                    (cutoff(i) - merge(0.0_dp, cutoff(max(i - 1, 1)), i == 1))
                n_lo = n_hi
            end do
            n_last = n_lo
            cutoff(n) = cutoff(n - 1) &
                + (light*phase/(4*pi*frequency(n)) - area)/(last_weight*n_last)
            if (.not. (cutoff(n) >= cutoff(n - 1))) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    'group delays give a non-monotone cutoff sequence')
                return
            end if
        end do
        call status_set(status, 0, '')
    end subroutine
end module tiago_profile_reflectometry
