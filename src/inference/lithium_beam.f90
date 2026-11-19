module tiago_lithium_beam
    !! Collisional-radiative populations of a neutral alkali beam along its
    !! path and the emission of one line, with derivatives in the electron
    !! density and temperature at the path nodes.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_DOMAIN_ERROR
    implicit none
    private

    type, public :: beam_atomic_data_t
        !! Rate coefficients <sigma v> in m^3/s, indexed (to, from). Electron
        !! rates are tabulated on ascending temperatures in eV and interpolated
        !! linearly in log temperature. Ion rates belong to the beam energy.
        real(dp), allocatable :: temperature(:)
        real(dp), allocatable :: electron_rate(:, :, :)
        real(dp), allocatable :: electron_loss(:, :)
        real(dp), allocatable :: ion_rate(:, :)
        real(dp), allocatable :: ion_loss(:)
        !! Spontaneous emission rates in 1/s, indexed (to, from).
        real(dp), allocatable :: einstein(:, :)
        real(dp) :: speed = 0.0_dp
        integer :: upper = 2, lower = 1
    end type

    type, public :: beam_emission_t
        !! Populations per injected atom, indexed (level, node).
        real(dp), allocatable :: population(:, :)
        !! Line emission A N_upper / v in photons per metre per injected atom.
        real(dp), allocatable :: emission(:)
        !! Derivatives of emission(k) in ne(j) and te(j), indexed (k, j).
        real(dp), allocatable :: emission_ne(:, :), emission_te(:, :)
    end type

    public :: lithium_beam_emission
contains

    subroutine lithium_beam_emission(atomic, z, ne, te, ion_ratio, result, status, &
            initial)
        !! Crank-Nicolson integration of dN/dz = M(z) N / v with M from the
        !! node-averaged density and temperature of each cell. Ion density is
        !! ion_ratio * ne. The beam enters at z(1) in the first level unless an
        !! initial population is given.
        type(beam_atomic_data_t), intent(in) :: atomic
        real(dp), intent(in) :: z(:), ne(:), te(:), ion_ratio
        type(beam_emission_t), intent(out) :: result
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(in), optional :: initial(:)
        integer :: nl, nn, k, j, info, side
        real(dp), allocatable :: lu(:, :, :), rhs_matrix(:, :, :), m_ne(:, :, :)
        real(dp), allocatable :: m_te(:, :, :), tangent(:), forcing(:)
        integer, allocatable :: pivots(:, :)
        real(dp) :: h, ne_mid, te_mid
        real(dp), allocatable :: rate(:, :), rate_t(:, :)
        external :: dgetrf, dgetrs

        nl = size(atomic%einstein, 1)
        nn = size(z)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid beam input')
        if (nn < 2 .or. size(ne) /= nn .or. size(te) /= nn) return
        if (.not. (atomic%speed > 0.0_dp .and. ion_ratio >= 0.0_dp)) return
        if (atomic%upper < 1 .or. atomic%upper > nl) return
        if (any(z(2:) <= z(:nn - 1)) .or. any(ne < 0.0_dp)) return
        if (any(te < atomic%temperature(1)) .or. &
            any(te > atomic%temperature(size(atomic%temperature)))) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'electron temperature outside the rate table')
            return
        end if
        allocate (result%population(nl, nn), result%emission(nn))
        allocate (result%emission_ne(nn, nn), result%emission_te(nn, nn))
        allocate (lu(nl, nl, nn - 1), rhs_matrix(nl, nl, nn - 1))
        allocate (m_ne(nl, nl, nn - 1), m_te(nl, nl, nn - 1), pivots(nl, nn - 1))
        allocate (rate(nl, nl), rate_t(nl, nl), tangent(nl), forcing(nl))
        result%population = 0.0_dp
        result%population(1, 1) = 1.0_dp
        if (present(initial)) then
            if (size(initial) /= nl) return
            result%population(:, 1) = initial
        end if
        do k = 1, nn - 1
            h = (z(k + 1) - z(k))/(2*atomic%speed)
            ne_mid = 0.5_dp*(ne(k) + ne(k + 1))
            te_mid = 0.5_dp*(te(k) + te(k + 1))
            call rate_matrix(atomic, ne_mid, te_mid, ion_ratio, rate, m_ne(:, :, k), &
                m_te(:, :, k))
            lu(:, :, k) = -h*rate
            rhs_matrix(:, :, k) = h*rate
            do j = 1, nl
                lu(j, j, k) = lu(j, j, k) + 1
                rhs_matrix(j, j, k) = rhs_matrix(j, j, k) + 1
            end do
            call dgetrf(nl, nl, lu(:, :, k), nl, pivots(:, k), info)
            if (info /= 0) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, 'singular beam step')
                return
            end if
            result%population(:, k + 1) = matmul(rhs_matrix(:, :, k), &
                result%population(:, k))
            call dgetrs('N', nl, 1, lu(:, :, k), nl, pivots(:, k), &
                result%population(:, k + 1), nl, info)
        end do
        result%emission = atomic%einstein(atomic%lower, atomic%upper)* &
            result%population(atomic%upper, :)/atomic%speed
        ! Node j enters the two adjacent cells with weight 1/2 each.
        result%emission_ne = 0.0_dp
        result%emission_te = 0.0_dp
        do j = 1, nn
            do side = 1, 2
                tangent = 0.0_dp
                do k = max(j - 1, 1), nn - 1
                    h = (z(k + 1) - z(k))/(2*atomic%speed)
                    forcing = matmul(rhs_matrix(:, :, k), tangent)
                    if (k == j - 1 .or. k == j) then
                        if (side == 1) then
                            forcing = forcing + 0.5_dp*h*matmul(m_ne(:, :, k), &
                                result%population(:, k) + result%population(:, k + 1))
                        else
                            forcing = forcing + 0.5_dp*h*matmul(m_te(:, :, k), &
                                result%population(:, k) + result%population(:, k + 1))
                        end if
                    end if
                    call dgetrs('N', nl, 1, lu(:, :, k), nl, pivots(:, k), forcing, &
                        nl, info)
                    tangent = forcing
                    if (side == 1) then
                        result%emission_ne(k + 1, j) = tangent(atomic%upper)
                    else
                        result%emission_te(k + 1, j) = tangent(atomic%upper)
                    end if
                end do
            end do
        end do
        result%emission_ne = atomic%einstein(atomic%lower, atomic%upper)* &
            result%emission_ne/atomic%speed
        result%emission_te = atomic%einstein(atomic%lower, atomic%upper)* &
            result%emission_te/atomic%speed
        if (.not. all(ieee_is_finite(result%population))) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'nonfinite beam population')
            return
        end if
        call status_set(status, 0, '')
    end subroutine

    subroutine rate_matrix(atomic, ne, te, ion_ratio, rate, rate_ne, rate_te)
        !! M = ne E(te) + ion_ratio ne I + A with gains off the diagonal and
        !! all losses, including ionisation, on it.
        type(beam_atomic_data_t), intent(in) :: atomic
        real(dp), intent(in) :: ne, te, ion_ratio
        real(dp), intent(out) :: rate(:, :), rate_ne(:, :), rate_te(:, :)
        real(dp) :: electron(size(rate, 1), size(rate, 1))
        real(dp) :: electron_t(size(rate, 1), size(rate, 1))
        real(dp) :: loss(size(rate, 1)), loss_t(size(rate, 1)), w, dw
        integer :: i, j, lo, nt
        nt = size(atomic%temperature)
        lo = 1
        do i = 1, nt - 1
            if (te >= atomic%temperature(i)) lo = i
        end do
        lo = min(lo, nt - 1)
        w = log(te/atomic%temperature(lo))/ &
            log(atomic%temperature(lo + 1)/atomic%temperature(lo))
        dw = 1/(te*log(atomic%temperature(lo + 1)/atomic%temperature(lo)))
        electron = (1 - w)*atomic%electron_rate(:, :, lo) &
            + w*atomic%electron_rate(:, :, lo + 1)
        electron_t = dw*(atomic%electron_rate(:, :, lo + 1) &
            - atomic%electron_rate(:, :, lo))
        loss = (1 - w)*atomic%electron_loss(:, lo) + w*atomic%electron_loss(:, lo + 1)
        loss_t = dw*(atomic%electron_loss(:, lo + 1) - atomic%electron_loss(:, lo))
        rate_ne = electron + ion_ratio*atomic%ion_rate
        rate_te = ne*electron_t
        rate = atomic%einstein
        do j = 1, size(rate, 1)
            rate_ne(j, j) = 0.0_dp
            rate_te(j, j) = 0.0_dp
            rate(j, j) = 0.0_dp
            rate_ne(j, j) = -(sum(electron(:, j)) - electron(j, j)) - loss(j) &
                - ion_ratio*(sum(atomic%ion_rate(:, j)) - atomic%ion_rate(j, j) &
                + atomic%ion_loss(j))
            rate_te(j, j) = -ne*(sum(electron_t(:, j)) - electron_t(j, j) + loss_t(j))
            rate(j, j) = -(sum(atomic%einstein(:, j)) - atomic%einstein(j, j))
        end do
        rate = rate + ne*rate_ne
    end subroutine
end module tiago_lithium_beam
