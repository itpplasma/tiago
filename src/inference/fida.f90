module tiago_fida
    !! Fast-ion D-alpha (FIDA) spectra on straight sight lines. A fast ion of
    !! velocity v (energy E_i, pitch p_j, gyro-phase gamma about the unit field
    !! b) charge-exchanges with beam neutrals of component l and with halo
    !! neutrals at the rate
    !!
    !!   nu = sum_l n_l sigma_cx(|v - v_l|) |v - v_l| + n_h <sigma_cx u>(|v|, T_i),
    !!
    !! (beam densities from tiago_beam_neutrals, rates from tiago_fusion_rates)
    !! and the new fast neutral emits
    !! Y photons of D-alpha (Y = photons_per_cx, collisional-radiative
    !! population of n = 3 not modelled) at the Doppler-shifted wavelength
    !!
    !!   lambda = lambda_0 (1 + v . d / c),
    !!
    !! with d the unit sight direction pointing away from the detector. Channel
    !! c, wavelength bin k of width dlambda_k returns the spectral radiance
    !! [photons s^-1 m^-2 sr^-1 nm^-1]
    !!
    !!   S_ck = scale_c / (4 pi dlambda_k) int_c dl sum_ij f_ij(s) dE_i dp_j
    !!          < Y nu chi_k(lambda) >_gamma,
    !!
    !! where chi_k is one inside bin k and the gyro-phase average uses ngyro
    !! midpoints. Stark splitting, fine structure, the instrument function and
    !! the beam-emission and bremsstrahlung background are not modelled. Signals
    !! are ordered by channel, then bin. The plasma enters through s (the
    !! distribution) and T_i (the halo rate). The field direction and the
    !! neutral densities are fixed at the sample points (fida_set_beams or
    !! fida_set_neutrals) and not differentiated.
    !!
    !! Weight functions (fida_weights) give S_ck = sum_ij W_ck(E_i, p_j) n_ij for
    !! a fast-ion density n_ij [m^-3] per cell that is uniform along the sight
    !! line.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t, Q_S, Q_TI
    use tiago_sight_geometry, only: chord_t, chord_rule
    use tiago_fast_ion_distribution, only: fast_ion_distribution_t
    use tiago_beam_neutrals, only: neutral_beam_t, beam_plasma_t, beam_densities, &
        beam_component_speed, BEAM_COMPONENTS
    use tiago_fusion_rates, only: thermal_rule_t, thermal_rule, maxwellian_rate, &
        reaction_sigma, REACTION_CX_H, SPEED_OF_LIGHT, DEUTERON_MASS
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter, public :: D_ALPHA_NM = 656.1_dp !! unshifted D-alpha [nm], air

    abstract interface
        pure function field_direction_i(x) result(b)
            !! Unit magnetic field direction at the Cartesian point x [m].
            import :: dp
            real(dp), intent(in) :: x(3)
            real(dp) :: b(3)
        end function field_direction_i
    end interface
    public :: field_direction_i

    type, public, extends(profile_diagnostic_t) :: fida_t
        type(fast_ion_distribution_t) :: dist
        type(thermal_rule_t) :: rule
        real(dp) :: lambda0 = D_ALPHA_NM
        real(dp) :: photons_per_cx = 1.0_dp !! Y
        real(dp) :: halo_mass = DEUTERON_MASS !! [kg]
        integer :: ngyro = 36
        real(dp), allocatable :: edges(:) !! wavelength bin edges [nm], increasing
        real(dp), allocatable :: speed(:) !! v_i [m/s]
        real(dp), allocatable :: x0(:, :), w(:) !! (R, phi, Z) and dl [m]
        real(dp), allocatable :: bhat(:, :), sight(:, :) !! (3, n) Cartesian unit vectors
        real(dp), allocatable :: nb(:, :), vb(:, :, :), nh(:) !! neutrals at the points
        integer, allocatable :: first(:)
        real(dp), allocatable :: scale(:)
    contains
        procedure :: init => fida_init
        procedure :: set_beams => fida_set_beams
        procedure :: set_neutrals => fida_set_neutrals
        procedure :: n_bins => fida_n_bins
        procedure :: n_signals => fida_n_signals
        procedure :: n_nuisance => fida_n_nuisance
        procedure :: label => fida_label
        procedure :: sample_points => fida_points
        procedure :: predict => fida_predict
        procedure :: jvp => fida_jvp
        procedure :: vjp => fida_vjp
        procedure :: weights => fida_weights
    end type fida_t

contains

    subroutine fida_init(self, chords, edges, dist, field_direction, status)
        !! Geometry and distribution, no neutrals yet. status 1: fewer than two
        !! or non-increasing bin edges; 2: invalid distribution.
        class(fida_t), intent(inout) :: self
        type(chord_t), intent(in) :: chords(:)
        real(dp), intent(in) :: edges(:)
        type(fast_ion_distribution_t), intent(in) :: dist
        procedure(field_direction_i) :: field_direction
        integer, intent(out) :: status
        real(dp), allocatable :: x(:, :), wl(:)
        real(dp) :: xc(3), d(3)
        integer :: c, m, n

        status = 1
        if (size(edges) < 2) return
        if (any(edges(2:) <= edges(:size(edges) - 1))) return
        status = 2
        if (dist%check() /= 0) return
        status = 0
        self%edges = edges
        self%dist = dist
        self%rule = thermal_rule(64)
        self%speed = dist%speed(dist%energy)

        if (allocated(self%x0)) deallocate(self%x0, self%w, self%first, self%bhat, self%sight, &
            self%nb, self%vb, self%nh)
        allocate(self%x0(3, 0), self%w(0), self%first(size(chords) + 1))
        self%first(1) = 1
        do c = 1, size(chords)
            call chord_rule(chords(c), x, wl)
            n = size(wl)
            self%x0 = reshape([self%x0, x], [3, size(self%w) + n])
            self%w = [self%w, wl]
            self%first(c + 1) = self%first(c) + n
        end do
        n = size(self%w)
        allocate(self%bhat(3, n), self%sight(3, n), self%nb(0, n), self%vb(3, 0, n), self%nh(n))
        self%nh = 0.0_dp
        do c = 1, size(chords)
            d = chords(c)%direction / norm2(chords(c)%direction)
            do m = self%first(c), self%first(c + 1) - 1
                xc = [self%x0(1, m) * cos(self%x0(2, m)), self%x0(1, m) * sin(self%x0(2, m)), &
                    self%x0(3, m)]
                self%sight(:, m) = d
                self%bhat(:, m) = field_direction(xc)
                self%bhat(:, m) = self%bhat(:, m) / norm2(self%bhat(:, m))
            end do
        end do
        if (allocated(self%scale)) deallocate(self%scale)
        allocate(self%scale(size(chords)))
        self%scale = 1.0_dp
    end subroutine fida_init

    subroutine fida_set_beams(self, beams, plasma, status)
        !! Component and halo densities of every beam at the sample points
        !! (tiago_beam_neutrals) for the given plasma. status: that of
        !! beam_densities for the first failing beam, else 0.
        class(fida_t), intent(inout) :: self
        type(neutral_beam_t), intent(in) :: beams(:)
        class(beam_plasma_t), intent(in) :: plasma
        integer, intent(out) :: status
        real(dp), allocatable :: xc(:, :), nb(:, :), vb(:, :, :), nh(:), nbi(:, :), nhi(:)
        real(dp) :: d(3)
        integer :: i, j, l, n

        n = size(self%w)
        xc = reshape([self%x0(1, :) * cos(self%x0(2, :)), self%x0(1, :) * sin(self%x0(2, :)), &
            self%x0(3, :)], [n, 3])
        xc = transpose(xc)
        allocate(nb(BEAM_COMPONENTS * size(beams), n), vb(3, BEAM_COMPONENTS * size(beams), n), &
            nh(n), nbi(BEAM_COMPONENTS, n), nhi(n))
        nh = 0.0_dp
        do i = 1, size(beams)
            call beam_densities(beams(i), plasma, xc, nbi, nhi, status)
            if (status /= 0) return
            d = beams(i)%direction / norm2(beams(i)%direction)
            do j = 1, BEAM_COMPONENTS
                l = (i - 1) * BEAM_COMPONENTS + j
                nb(l, :) = nbi(j, :)
                vb(:, l, :) = spread(beam_component_speed(beams(i), j) * d, 2, n)
            end do
            nh = nh + nhi
        end do
        call self%set_neutrals(nb, vb, nh)
    end subroutine fida_set_beams

    subroutine fida_set_neutrals(self, nb, vb, nh)
        !! Neutral densities nb(l, m) [m^-3] with velocities vb(3, l, m) [m/s]
        !! and halo densities nh(m) [m^-3] at the sample points m.
        class(fida_t), intent(inout) :: self
        real(dp), intent(in) :: nb(:, :), vb(:, :, :), nh(:)
        self%nb = nb
        self%vb = vb
        self%nh = nh
    end subroutine fida_set_neutrals

    pure integer function fida_n_bins(self)
        class(fida_t), intent(in) :: self
        fida_n_bins = size(self%edges) - 1
    end function fida_n_bins

    integer function fida_n_signals(self)
        class(fida_t), intent(in) :: self
        fida_n_signals = (size(self%first) - 1) * self%n_bins()
    end function fida_n_signals

    integer function fida_n_nuisance(self)
        class(fida_t), intent(in) :: self
        associate (unused => self)
        end associate
        fida_n_nuisance = 0
    end function fida_n_nuisance

    function fida_label(self, i) result(text)
        class(fida_t), intent(in) :: self
        integer, intent(in) :: i
        character(len=:), allocatable :: text
        character(len=48) :: buffer
        integer :: c, k
        c = (i - 1) / self%n_bins() + 1
        k = i - (c - 1) * self%n_bins()
        write (buffer, '(A,I0,A,F0.3)') 'fida_', c, '_', 0.5_dp * (self%edges(k) + self%edges(k + 1))
        text = trim(buffer)
    end function fida_label

    subroutine fida_points(self, eta, x, dx_deta)
        class(fida_t), intent(in) :: self
        real(dp), intent(in) :: eta(:)
        real(dp), allocatable, intent(out) :: x(:, :), dx_deta(:, :, :)
        associate (unused => eta)
        end associate
        x = self%x0
        allocate(dx_deta(3, size(x, 2), 0))
    end subroutine fida_points

    pure integer function bin_of(edges, lambda) result(k)
        real(dp), intent(in) :: edges(:), lambda
        integer :: lo, hi, mid
        k = 0
        if (lambda < edges(1) .or. lambda >= edges(size(edges))) return
        lo = 1
        hi = size(edges)
        do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (edges(mid) <= lambda) then
                lo = mid
            else
                hi = mid
            end if
        end do
        k = lo
    end function bin_of

    pure subroutine cell_spectrum(self, m, i, j, rh, drh, spec, dspec)
        !! Adds Y nu chi_k / (4 pi), gyro-averaged, of one fast ion in cell
        !! (i, j) at point m to spec(k) [photons s^-1 sr^-1], and its T_i
        !! derivative to dspec.
        class(fida_t), intent(in) :: self
        integer, intent(in) :: m, i, j
        real(dp), intent(in) :: rh, drh
        real(dp), intent(inout) :: spec(:), dspec(:)
        real(dp) :: b(3), e1(3), e2(3), v(3), rel(3), p, qp, u, sg, dsg, nu, gam, c0
        integer :: g, l, k

        b = self%bhat(:, m)
        if (abs(b(3)) < 0.9_dp) then
            e1 = [b(2), -b(1), 0.0_dp]
        else
            e1 = [0.0_dp, b(3), -b(2)]
        end if
        e1 = e1 / norm2(e1)
        e2 = [b(2) * e1(3) - b(3) * e1(2), b(3) * e1(1) - b(1) * e1(3), b(1) * e1(2) - b(2) * e1(1)]
        p = self%dist%pitch(j)
        qp = sqrt(max(0.0_dp, 1.0_dp - p**2))
        c0 = self%photons_per_cx / (4.0_dp * pi * self%ngyro)
        do g = 1, self%ngyro
            gam = 2.0_dp * pi * (g - 0.5_dp) / self%ngyro
            v = self%speed(i) * (p * b + qp * (cos(gam) * e1 + sin(gam) * e2))
            k = bin_of(self%edges, self%lambda0 * (1.0_dp + dot_product(v, self%sight(:, m)) &
                / SPEED_OF_LIGHT))
            if (k == 0) cycle
            nu = self%nh(m) * rh
            do l = 1, size(self%nb, 1)
                if (self%nb(l, m) == 0.0_dp) cycle
                rel = v - self%vb(:, l, m)
                u = norm2(rel)
                call reaction_sigma(REACTION_CX_H, u, sg, dsg)
                nu = nu + self%nb(l, m) * sg * u
            end do
            spec(k) = spec(k) + c0 * nu
            dspec(k) = dspec(k) + c0 * self%nh(m) * drh
        end do
    end subroutine cell_spectrum

    subroutine fida_point(self, m, q, val, d_s, d_ti)
        !! Spectral emission per bin at point m, per unit length and bin width
        !! [photons s^-1 m^-3 sr^-1 nm^-1], and its s and T_i derivatives.
        class(fida_t), intent(in) :: self
        integer, intent(in) :: m
        real(dp), intent(in) :: q(:)
        real(dp), intent(out) :: val(:), d_s(:), d_ti(:)
        real(dp) :: spec(size(val)), dspec(size(val)), t, dt_ds, cell, dcell
        real(dp), allocatable :: rh(:), drh(:)
        integer :: k, i, j, ne

        val = 0.0_dp
        d_s = 0.0_dp
        d_ti = 0.0_dp
        call self%dist%locate(q(Q_S), k, t, dt_ds)
        if (k == 0) return
        ne = size(self%speed)
        allocate(rh(ne), drh(ne))
        do i = 1, ne
            call maxwellian_rate(self%rule, REACTION_CX_H, self%speed(i), q(Q_TI), &
                self%halo_mass, rh(i), drh(i))
        end do
        do j = 1, size(self%dist%pitch)
            do i = 1, ne
                if (k < size(self%dist%s)) then
                    cell = (1.0_dp - t) * self%dist%f(i, j, k) + t * self%dist%f(i, j, k + 1)
                    dcell = (self%dist%f(i, j, k + 1) - self%dist%f(i, j, k)) * dt_ds
                else
                    cell = self%dist%f(i, j, k)
                    dcell = 0.0_dp
                end if
                if (cell == 0.0_dp .and. dcell == 0.0_dp) cycle
                cell = cell * self%dist%denergy(i) * self%dist%dpitch(j)
                dcell = dcell * self%dist%denergy(i) * self%dist%dpitch(j)
                spec = 0.0_dp
                dspec = 0.0_dp
                call cell_spectrum(self, m, i, j, rh(i), drh(i), spec, dspec)
                val = val + cell * spec
                d_s = d_s + dcell * spec
                d_ti = d_ti + cell * dspec
            end do
        end do
        val = val / (self%edges(2:) - self%edges(:size(self%edges) - 1))
        d_s = d_s / (self%edges(2:) - self%edges(:size(self%edges) - 1))
        d_ti = d_ti / (self%edges(2:) - self%edges(:size(self%edges) - 1))
    end subroutine fida_point

    subroutine fida_predict(self, sample, eta, y)
        class(fida_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: val(self%n_bins()), d_s(self%n_bins()), d_ti(self%n_bins())
        integer :: c, m, nb
        associate (unused => eta)
        end associate
        nb = self%n_bins()
        y = 0.0_dp
        do c = 1, size(self%first) - 1
            do m = self%first(c), self%first(c + 1) - 1
                call fida_point(self, m, sample%q(:, m), val, d_s, d_ti)
                y((c - 1) * nb + 1:c * nb) = y((c - 1) * nb + 1:c * nb) &
                    + self%scale(c) * self%w(m) * val
            end do
        end do
    end subroutine fida_predict

    subroutine fida_jvp(self, sample, eta, q_dot, eta_dot, y_dot)
        class(fida_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        real(dp) :: val(self%n_bins()), d_s(self%n_bins()), d_ti(self%n_bins())
        integer :: c, m, nb
        associate (unused => eta, unused2 => eta_dot)
        end associate
        nb = self%n_bins()
        y_dot = 0.0_dp
        do c = 1, size(self%first) - 1
            do m = self%first(c), self%first(c + 1) - 1
                call fida_point(self, m, sample%q(:, m), val, d_s, d_ti)
                y_dot((c - 1) * nb + 1:c * nb) = y_dot((c - 1) * nb + 1:c * nb) &
                    + self%scale(c) * self%w(m) * (d_s * q_dot(Q_S, m) + d_ti * q_dot(Q_TI, m))
            end do
        end do
    end subroutine fida_jvp

    subroutine fida_vjp(self, sample, eta, y_bar, q_bar, eta_bar)
        class(fida_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        real(dp) :: val(self%n_bins()), d_s(self%n_bins()), d_ti(self%n_bins())
        integer :: c, m, nb
        associate (unused => eta)
        end associate
        nb = self%n_bins()
        q_bar = 0.0_dp
        eta_bar = 0.0_dp
        do c = 1, size(self%first) - 1
            do m = self%first(c), self%first(c + 1) - 1
                call fida_point(self, m, sample%q(:, m), val, d_s, d_ti)
                q_bar(Q_S, m) = self%scale(c) * self%w(m) &
                    * dot_product(d_s, y_bar((c - 1) * nb + 1:c * nb))
                q_bar(Q_TI, m) = self%scale(c) * self%w(m) &
                    * dot_product(d_ti, y_bar((c - 1) * nb + 1:c * nb))
            end do
        end do
    end subroutine fida_vjp

    subroutine fida_weights(self, sample, c, weights)
        !! weights(k, i, j) = W_ck(E_i, p_j) [photons s^-1 m^-2 sr^-1 nm^-1 per
        !! m^-3 of fast ions in cell (i, j), uniform along the sight line].
        class(fida_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        integer, intent(in) :: c
        real(dp), intent(out) :: weights(:, :, :)
        real(dp) :: spec(self%n_bins()), dspec(self%n_bins()), rh, drh
        integer :: m, i, j

        weights = 0.0_dp
        do m = self%first(c), self%first(c + 1) - 1
            do i = 1, size(self%speed)
                call maxwellian_rate(self%rule, REACTION_CX_H, self%speed(i), sample%q(Q_TI, m), &
                    self%halo_mass, rh, drh)
                do j = 1, size(self%dist%pitch)
                    spec = 0.0_dp
                    dspec = 0.0_dp
                    call cell_spectrum(self, m, i, j, rh, drh, spec, dspec)
                    weights(:, i, j) = weights(:, i, j) + self%scale(c) * self%w(m) * spec
                end do
            end do
        end do
        do j = 1, size(self%dist%pitch)
            do i = 1, size(self%speed)
                weights(:, i, j) = weights(:, i, j) / (self%edges(2:) - self%edges(:size(self%edges) - 1))
            end do
        end do
    end subroutine fida_weights

end module tiago_fida
