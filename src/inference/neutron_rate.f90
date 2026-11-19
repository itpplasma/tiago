module tiago_neutron_rate
    !! D-D neutron emission from thermal deuterium and a fast-deuteron
    !! distribution f(s, E, p) (tiago_fast_ion_distribution). The emissivity
    !! [m^-3 s^-1] at a sample point is the sum of
    !!
    !!   thermonuclear  eps_tn = n_D^2 <sigma v>_BH(T_i) / 2,
    !!   beam-target    eps_bt = n_D sum_i F_i(s) <sigma u>(v_i, T_i),
    !!   beam-beam      eps_bb = 1/2 sum_a sum_b f_a f_b <sigma(|v_a - v_b|) |v_a - v_b|>_gyro,
    !!
    !! with n_D = n_e (Z_imp - Z_eff) / (Z_imp - 1) (one impurity of charge
    !! Z_imp), <sigma v>_BH the Bosch-Hale reactivity, F_i = sum_j f_ij dE_i dp_j
    !! the fast density at energy E_i and <sigma u> the Maxwellian-target rate
    !! of tiago_fusion_rates (thermal target at rest). The beam-beam gyro
    !! average is over the gyro-phase difference of the two cells, with
    !! |v_a - v_b|^2 = v_a^2 + v_b^2 - 2 v_a v_b (p_a p_b + q_a q_b cos psi),
    !! q = sqrt(1 - p^2). f is linear in s between nodes, so eps_bb is the
    !! quadratic form of the node distributions, precomputed per node pair.
    !! Emission is isotropic, target rotation and the reaction kinematics of
    !! the detector response are not modelled.
    !!
    !! Channel c returns scale_c sum_m w_m eps(x_m) over its points: with
    !! volume weights [m^3] the total rate [1/s] (neutron_volume_points), with
    !! chord weights dl / (4 pi) [m/sr] the line-integrated emission
    !! [m^-2 s^-1 sr^-1] (neutron_chord_points). A fast species other than
    !! deuterium contributes no D-D neutrons.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t, Q_S, Q_NE, &
        Q_TI, Q_ZEFF, N_QUANTITIES
    use tiago_sight_geometry, only: chord_t, chord_rule, gauss_legendre01
    use tiago_fast_ion_distribution, only: fast_ion_distribution_t
    use tiago_fusion_rates, only: thermal_rule_t, thermal_rule, maxwellian_rate, &
        reaction_sigma, dd_thermal_reactivity, REACTION_DD_N, DEUTERON_MASS
    implicit none
    private

    integer, parameter, public :: NEUTRON_THERMONUCLEAR = 1, NEUTRON_BEAM_TARGET = 2, &
        NEUTRON_BEAM_BEAM = 3

    type, public, extends(profile_diagnostic_t) :: neutron_t
        type(fast_ion_distribution_t) :: dist
        type(thermal_rule_t) :: rule
        logical :: fast_deuterium = .false.
        real(dp) :: z_imp = 7.0_dp !! impurity charge (nitrogen)
        real(dp), allocatable :: speed(:) !! v_i [m/s]
        real(dp), allocatable :: fe(:, :) !! (ne, ns) fast density per energy cell [m^-3]
        real(dp), allocatable :: bb_node(:), bb_next(:) !! B(k, k), B(k, k + 1) [m^-3 s^-1]
        real(dp), allocatable :: x0(:, :), w(:) !! sample points (R, phi, Z), weights
        integer, allocatable :: first(:) !! points of channel c: first(c) .. first(c + 1) - 1
        real(dp), allocatable :: scale(:)
    contains
        procedure :: setup => neutron_setup
        procedure :: init => neutron_init
        procedure :: emissivity => neutron_emissivity
        procedure :: components => neutron_components
        procedure :: n_signals => neutron_n_signals
        procedure :: n_nuisance => neutron_n_nuisance
        procedure :: label => neutron_label
        procedure :: sample_points => neutron_points
        procedure :: predict => neutron_predict
        procedure :: jvp => neutron_jvp
        procedure :: vjp => neutron_vjp
    end type neutron_t

    public :: neutron_chord_points, neutron_volume_points

contains

    subroutine neutron_setup(self, dist, z_imp, ngyro, nthermal)
        !! Emissivity model without sample points. ngyro: Gauss points of the
        !! beam-beam gyro-phase average on [0, pi] (default 24); nthermal:
        !! nodes of the thermal-target average (default 96).
        class(neutron_t), intent(inout) :: self
        type(fast_ion_distribution_t), intent(in) :: dist
        real(dp), intent(in), optional :: z_imp
        integer, intent(in), optional :: ngyro, nthermal
        real(dp), allocatable :: kernel(:, :), fk(:), fl(:), xg(:), wg(:), p(:), q(:), v(:), dw(:)
        real(dp) :: cospsi, u2, sg, dsg, acc
        integer :: ne, np, ns, nc, a, b, k, g, ng

        self%dist = dist
        if (present(z_imp)) self%z_imp = z_imp
        self%rule = thermal_rule(96)
        if (present(nthermal)) self%rule = thermal_rule(nthermal)
        ng = 24
        if (present(ngyro)) ng = ngyro
        ne = size(dist%energy)
        np = size(dist%pitch)
        ns = size(dist%s)
        self%fast_deuterium = dist%is_deuterium()
        self%speed = dist%speed(dist%energy)
        allocate(self%fe(ne, ns), self%bb_node(ns), self%bb_next(ns))
        do k = 1, ns
            self%fe(:, k) = matmul(dist%f(:, :, k), dist%dpitch) * dist%denergy
        end do
        self%bb_node = 0.0_dp
        self%bb_next = 0.0_dp
        if (.not. self%fast_deuterium) return

        ! Cells a = (i, j) in array element order; kernel(a, b) is the gyro
        ! average of sigma u times both cell volumes dE dp.
        nc = ne * np
        allocate(kernel(nc, nc), p(nc), q(nc), v(nc), dw(nc))
        do a = 1, nc
            p(a) = dist%pitch((a - 1) / ne + 1)
            q(a) = sqrt(max(0.0_dp, 1.0_dp - p(a)**2))
            v(a) = self%speed(mod(a - 1, ne) + 1)
            dw(a) = dist%denergy(mod(a - 1, ne) + 1) * dist%dpitch((a - 1) / ne + 1)
        end do
        call gauss_legendre01(ng, xg, wg)
        !$omp parallel do private(a, b, g, cospsi, u2, sg, dsg, acc) schedule(dynamic)
        do b = 1, nc
            do a = 1, b
                acc = 0.0_dp
                do g = 1, ng
                    cospsi = cos(acos(-1.0_dp) * xg(g))
                    u2 = v(a)**2 + v(b)**2 - 2.0_dp * v(a) * v(b) * (p(a) * p(b) + q(a) * q(b) * cospsi)
                    u2 = sqrt(max(u2, 0.0_dp))
                    call reaction_sigma(REACTION_DD_N, u2, sg, dsg)
                    acc = acc + wg(g) * sg * u2
                end do
                kernel(a, b) = acc * dw(a) * dw(b)
                kernel(b, a) = kernel(a, b)
            end do
        end do
        !$omp end parallel do
        do k = 1, ns
            fk = reshape(dist%f(:, :, k), [nc])
            self%bb_node(k) = 0.5_dp * dot_product(fk, matmul(kernel, fk))
            if (k < ns) then
                fl = reshape(dist%f(:, :, k + 1), [nc])
                self%bb_next(k) = 0.5_dp * dot_product(fk, matmul(kernel, fl))
            end if
        end do
    end subroutine neutron_setup

    subroutine neutron_init(self, dist, points, weights, first, z_imp)
        !! points (3, n): (R, phi, Z) [m, rad, m]; weights (n); channel c owns
        !! points first(c) .. first(c + 1) - 1.
        class(neutron_t), intent(inout) :: self
        type(fast_ion_distribution_t), intent(in) :: dist
        real(dp), intent(in) :: points(:, :), weights(:)
        integer, intent(in) :: first(:)
        real(dp), intent(in), optional :: z_imp
        call self%setup(dist, z_imp)
        self%x0 = points
        self%w = weights
        self%first = first
        if (allocated(self%scale)) deallocate(self%scale)
        allocate(self%scale(size(first) - 1))
        self%scale = 1.0_dp
    end subroutine neutron_init

    pure subroutine neutron_emissivity(self, q, eps, deps)
        !! eps(3): thermonuclear, beam-target, beam-beam [m^-3 s^-1] at the
        !! sampled quantities q(N_QUANTITIES); deps(3, N_QUANTITIES): d eps / d q.
        class(neutron_t), intent(in) :: self
        real(dp), intent(in) :: q(:)
        real(dp), intent(out) :: eps(3), deps(3, N_QUANTITIES)
        real(dp) :: nd, dnd_dne, dnd_dz, sv, dsv, t, dt_ds, r, dr, fe, dfe, bt, dbt_ds, dbt_dti
        integer :: k, i

        eps = 0.0_dp
        deps = 0.0_dp
        dnd_dne = (self%z_imp - q(Q_ZEFF)) / (self%z_imp - 1.0_dp)
        dnd_dz = -q(Q_NE) / (self%z_imp - 1.0_dp)
        nd = q(Q_NE) * dnd_dne

        call dd_thermal_reactivity(q(Q_TI), sv, dsv)
        eps(1) = 0.5_dp * nd**2 * sv
        deps(1, Q_NE) = nd * sv * dnd_dne
        deps(1, Q_ZEFF) = nd * sv * dnd_dz
        deps(1, Q_TI) = 0.5_dp * nd**2 * dsv
        if (.not. self%fast_deuterium) return

        call self%dist%locate(q(Q_S), k, t, dt_ds)
        if (k == 0) return
        bt = 0.0_dp
        dbt_ds = 0.0_dp
        dbt_dti = 0.0_dp
        do i = 1, size(self%speed)
            if (k < size(self%dist%s)) then
                fe = (1.0_dp - t) * self%fe(i, k) + t * self%fe(i, k + 1)
                dfe = (self%fe(i, k + 1) - self%fe(i, k)) * dt_ds
            else
                fe = self%fe(i, k)
                dfe = 0.0_dp
            end if
            if (fe == 0.0_dp .and. dfe == 0.0_dp) cycle
            call maxwellian_rate(self%rule, REACTION_DD_N, self%speed(i), q(Q_TI), &
                DEUTERON_MASS, r, dr)
            bt = bt + fe * r
            dbt_ds = dbt_ds + dfe * r
            dbt_dti = dbt_dti + fe * dr
        end do
        eps(2) = nd * bt
        deps(2, Q_NE) = dnd_dne * bt
        deps(2, Q_ZEFF) = dnd_dz * bt
        deps(2, Q_TI) = nd * dbt_dti
        deps(2, Q_S) = nd * dbt_ds

        if (k < size(self%dist%s)) then
            eps(3) = (1.0_dp - t)**2 * self%bb_node(k) + 2.0_dp * t * (1.0_dp - t) * self%bb_next(k) &
                + t**2 * self%bb_node(k + 1)
            deps(3, Q_S) = dt_ds * (-2.0_dp * (1.0_dp - t) * self%bb_node(k) &
                + 2.0_dp * (1.0_dp - 2.0_dp * t) * self%bb_next(k) + 2.0_dp * t * self%bb_node(k + 1))
        else
            eps(3) = self%bb_node(k)
        end if
    end subroutine neutron_emissivity

    subroutine neutron_components(self, sample, y)
        !! y(3, n_signals): thermonuclear, beam-target and beam-beam parts.
        class(neutron_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(out) :: y(:, :)
        real(dp) :: eps(3), deps(3, N_QUANTITIES)
        integer :: c, m
        y = 0.0_dp
        do c = 1, size(self%first) - 1
            do m = self%first(c), self%first(c + 1) - 1
                call self%emissivity(sample%q(:, m), eps, deps)
                y(:, c) = y(:, c) + self%scale(c) * self%w(m) * eps
            end do
        end do
    end subroutine neutron_components

    integer function neutron_n_signals(self)
        class(neutron_t), intent(in) :: self
        neutron_n_signals = size(self%first) - 1
    end function neutron_n_signals

    integer function neutron_n_nuisance(self)
        class(neutron_t), intent(in) :: self
        associate (unused => self)
        end associate
        neutron_n_nuisance = 0
    end function neutron_n_nuisance

    function neutron_label(self, i) result(text)
        class(neutron_t), intent(in) :: self
        integer, intent(in) :: i
        character(len=:), allocatable :: text
        character(len=16) :: buffer
        associate (unused => self)
        end associate
        write (buffer, '(I0)') i
        text = 'neutron_'//trim(buffer)
    end function neutron_label

    subroutine neutron_points(self, eta, x, dx_deta)
        class(neutron_t), intent(in) :: self
        real(dp), intent(in) :: eta(:)
        real(dp), allocatable, intent(out) :: x(:, :), dx_deta(:, :, :)
        associate (unused => eta)
        end associate
        x = self%x0
        allocate(dx_deta(3, size(x, 2), 0))
    end subroutine neutron_points

    subroutine neutron_predict(self, sample, eta, y)
        class(neutron_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: parts(3, size(y))
        associate (unused => eta)
        end associate
        call self%components(sample, parts)
        y = sum(parts, dim=1)
    end subroutine neutron_predict

    subroutine neutron_jvp(self, sample, eta, q_dot, eta_dot, y_dot)
        class(neutron_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        real(dp) :: eps(3), deps(3, N_QUANTITIES)
        integer :: c, m
        associate (unused => eta, unused2 => eta_dot)
        end associate
        y_dot = 0.0_dp
        do c = 1, size(self%first) - 1
            do m = self%first(c), self%first(c + 1) - 1
                call self%emissivity(sample%q(:, m), eps, deps)
                y_dot(c) = y_dot(c) + self%scale(c) * self%w(m) * sum(matmul(deps, q_dot(:, m)))
            end do
        end do
    end subroutine neutron_jvp

    subroutine neutron_vjp(self, sample, eta, y_bar, q_bar, eta_bar)
        class(neutron_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        real(dp) :: eps(3), deps(3, N_QUANTITIES)
        integer :: c, m
        associate (unused => eta)
        end associate
        q_bar = 0.0_dp
        eta_bar = 0.0_dp
        do c = 1, size(self%first) - 1
            do m = self%first(c), self%first(c + 1) - 1
                call self%emissivity(sample%q(:, m), eps, deps)
                q_bar(:, m) = q_bar(:, m) + self%scale(c) * self%w(m) * y_bar(c) * sum(deps, dim=1)
            end do
        end do
    end subroutine neutron_vjp

    subroutine neutron_chord_points(chords, points, weights, first)
        !! Sight-line rule of each chord with weights dl / (4 pi) [m/sr].
        type(chord_t), intent(in) :: chords(:)
        real(dp), allocatable, intent(out) :: points(:, :), weights(:)
        integer, allocatable, intent(out) :: first(:)
        real(dp), allocatable :: x(:, :), w(:)
        integer :: c, n

        allocate(points(3, 0), weights(0), first(size(chords) + 1))
        first(1) = 1
        do c = 1, size(chords)
            call chord_rule(chords(c), x, w)
            n = size(w)
            points = reshape([points, x], [3, size(weights) + n])
            weights = [weights, w / (4.0_dp * acos(-1.0_dp))]
            first(c + 1) = first(c) + n
        end do
    end subroutine neutron_chord_points

    subroutine neutron_volume_points(r_min, r_max, z_min, z_max, nr, nz, points, weights, first)
        !! One channel: midpoints of an nr x nz grid in (R, Z) at phi = 0 with
        !! volume weights 2 pi R dR dZ [m^3] (axisymmetric total rate).
        real(dp), intent(in) :: r_min, r_max, z_min, z_max
        integer, intent(in) :: nr, nz
        real(dp), allocatable, intent(out) :: points(:, :), weights(:)
        integer, allocatable, intent(out) :: first(:)
        real(dp) :: dr, dz, r
        integer :: i, j, m

        dr = (r_max - r_min) / nr
        dz = (z_max - z_min) / nz
        allocate(points(3, nr * nz), weights(nr * nz))
        m = 0
        do j = 1, nz
            do i = 1, nr
                m = m + 1
                r = r_min + (i - 0.5_dp) * dr
                points(:, m) = [r, 0.0_dp, z_min + (j - 0.5_dp) * dz]
                weights(m) = 2.0_dp * acos(-1.0_dp) * r * dr * dz
            end do
        end do
        first = [1, nr * nz + 1]
    end subroutine neutron_volume_points

end module tiago_neutron_rate
