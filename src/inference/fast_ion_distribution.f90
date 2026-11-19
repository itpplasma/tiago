module tiago_fast_ion_distribution
    !! Gyrotropic fast-ion distribution f(s, E, p) on flux-label nodes s_k and
    !! energy-pitch cells (E_i, p_j) of widths dE_i, dp_j, with p = v_par / v
    !! relative to B. f is the number density per energy and unit pitch
    !! [m^-3 eV^-1], so the density in cell (i, j) is f dE_i dp_j. Between nodes
    !! f is linear in s, below s_1 it equals f(s_1) and above s_ns it vanishes,
    !! so the last node should carry f = 0. The distribution is an input from an
    !! orbit or Fokker-Planck code and carries no derivatives with respect to
    !! its own parameters.
    !!
    !! Text format (list-directed, # comment lines only before the data):
    !!   mass_amu
    !!   ns ne np
    !!   s(1:ns)
    !!   E(1:ne) [eV]  dE(1:ne) [eV]
    !!   p(1:np)       dp(1:np)
    !!   f(1:ne, 1:np, 1:ns) in array element order [m^-3 eV^-1]
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_fusion_rates, only: ATOMIC_MASS, ELEMENTARY_CHARGE
    implicit none
    private

    type, public :: fast_ion_distribution_t
        real(dp) :: mass_amu = 2.01410177812_dp
        real(dp), allocatable :: s(:), energy(:), denergy(:), pitch(:), dpitch(:)
        real(dp), allocatable :: f(:, :, :) !! (ne, np, ns)
    contains
        procedure :: locate => distribution_locate
        procedure :: speed => distribution_speed
        procedure :: is_deuterium => distribution_is_deuterium
        procedure :: check => distribution_check
    end type fast_ion_distribution_t

    public :: read_fast_ion_distribution

contains

    pure subroutine distribution_locate(self, s, k, t, dt_ds)
        !! f(s) = (1 - t) f(:, :, k) + t f(:, :, k + 1); k = 0 outside the support.
        class(fast_ion_distribution_t), intent(in) :: self
        real(dp), intent(in) :: s
        integer, intent(out) :: k
        real(dp), intent(out) :: t, dt_ds
        integer :: n

        n = size(self%s)
        t = 0.0_dp
        dt_ds = 0.0_dp
        k = 0
        if (s >= self%s(n) .or. n < 2) return
        if (s <= self%s(1)) then
            k = 1
            return
        end if
        k = 1
        do while (self%s(k + 1) <= s)
            k = k + 1
        end do
        dt_ds = 1.0_dp / (self%s(k + 1) - self%s(k))
        t = (s - self%s(k)) * dt_ds
    end subroutine distribution_locate

    elemental real(dp) function distribution_speed(self, energy) result(v)
        !! Speed [m/s] at the energy [eV], non-relativistic.
        class(fast_ion_distribution_t), intent(in) :: self
        real(dp), intent(in) :: energy
        v = sqrt(2.0_dp * energy * ELEMENTARY_CHARGE / (self%mass_amu * ATOMIC_MASS))
    end function distribution_speed

    logical function distribution_is_deuterium(self)
        class(fast_ion_distribution_t), intent(in) :: self
        distribution_is_deuterium = abs(self%mass_amu - 2.014_dp) < 0.05_dp
    end function distribution_is_deuterium

    integer function distribution_check(self) result(status)
        !! 0: valid; 1: shapes inconsistent; 2: s not increasing or a non-positive
        !! energy or width; 3: pitch outside [-1, 1]; 4: negative density.
        class(fast_ion_distribution_t), intent(in) :: self
        status = 1
        if (.not. (allocated(self%s) .and. allocated(self%energy) .and. allocated(self%pitch) &
            .and. allocated(self%f))) return
        if (size(self%denergy) /= size(self%energy) .or. size(self%dpitch) /= size(self%pitch)) &
            return
        if (any(shape(self%f) /= [size(self%energy), size(self%pitch), size(self%s)])) return
        if (size(self%s) < 2) return
        status = 2
        if (any(self%s(2:) <= self%s(:size(self%s) - 1))) return
        if (any(self%energy <= 0.0_dp) .or. any(self%denergy <= 0.0_dp)) return
        if (any(self%dpitch <= 0.0_dp)) return
        status = 3
        if (any(abs(self%pitch) > 1.0_dp)) return
        status = 4
        if (any(self%f < 0.0_dp)) return
        status = 0
    end function distribution_check

    subroutine read_fast_ion_distribution(path, dist, status)
        !! status: 0 read and valid; -1 unreadable; > 0 from dist%check.
        character(len=*), intent(in) :: path
        type(fast_ion_distribution_t), intent(out) :: dist
        integer, intent(out) :: status
        integer :: unit, ios, ns, ne, np
        character(len=512) :: line

        status = -1
        open (newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        do
            read (unit, '(A)', iostat=ios) line
            if (ios /= 0) then
                close (unit)
                return
            end if
            if (len_trim(line) == 0) cycle
            if (adjustl(line(1:1)) /= '#') exit
        end do
        backspace (unit)
        read (unit, *, iostat=ios) dist%mass_amu
        if (ios == 0) read (unit, *, iostat=ios) ns, ne, np
        if (ios /= 0 .or. ns < 2 .or. ne < 1 .or. np < 1) then
            close (unit)
            return
        end if
        allocate(dist%s(ns), dist%energy(ne), dist%denergy(ne), dist%pitch(np), &
            dist%dpitch(np), dist%f(ne, np, ns))
        read (unit, *, iostat=ios) dist%s
        if (ios == 0) read (unit, *, iostat=ios) dist%energy, dist%denergy
        if (ios == 0) read (unit, *, iostat=ios) dist%pitch, dist%dpitch
        if (ios == 0) read (unit, *, iostat=ios) dist%f
        close (unit)
        if (ios /= 0) return
        status = dist%check()
    end subroutine read_fast_ion_distribution

end module tiago_fast_ion_distribution
