module fida_test_providers
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_beam_neutrals, only: beam_plasma_t
    implicit none
    private

    type, public, extends(beam_plasma_t) :: slab_plasma_t
        !! Uniform plasma in 1.3 <= R <= 1.9 m, |Z| <= 0.4 m.
    contains
        procedure :: eval => slab_eval
    end type slab_plasma_t

    public :: tilted_field

contains

    subroutine slab_eval(self, x, ne, te, ti, zeff)
        class(slab_plasma_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: ne, te, ti, zeff
        real(dp) :: r
        associate (unused => self)
        end associate
        r = hypot(x(1), x(2))
        ne = 0.0_dp
        te = 0.0_dp
        ti = 0.0_dp
        zeff = 1.0_dp
        if (r < 1.3_dp .or. r > 1.9_dp .or. abs(x(3)) > 0.4_dp) return
        ne = 4.0e19_dp
        te = 2.0e3_dp
        ti = 1.5e3_dp
        zeff = 1.6_dp
    end subroutine slab_eval

    pure function tilted_field(x) result(b)
        !! Uniform unit field at 50 degrees to the sight lines along +x.
        real(dp), intent(in) :: x(3)
        real(dp) :: b(3)
        real(dp), parameter :: angle = 50.0_dp * acos(-1.0_dp) / 180.0_dp
        associate (unused => x)
        end associate
        b = [cos(angle), sin(angle), 0.0_dp]
    end function tilted_field
end module fida_test_providers

program test_fida
    !! FIDA spectra (tiago_fida). Oracles: the closed-form gyro-phase
    !! distribution of the projected velocity u = v (p cos a + q sin a cos g)
    !! (Salewski et al., Nucl. Fusion 51, 083014 (2011)), whose probability in
    !! [u1, u2] is (arccos x1 - arccos x2) / pi with x = (u - v p cos a) /
    !! (v q sin a); the flat Doppler box lambda_0 (1 +- v / c) of an isotropic
    !! monoenergetic population; the weight-function identity S = sum W n for
    !! a distribution uniform along the sight line; the contract's derivative
    !! check with beam and halo neutrals.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_profile_diagnostic, only: plasma_sample_t, Q_S, Q_B, Q_NE, Q_TE, Q_TI, Q_ZEFF, &
        N_QUANTITIES, check_diagnostic_derivatives
    use tiago_sight_geometry, only: chord_t
    use tiago_fast_ion_distribution, only: fast_ion_distribution_t
    use tiago_beam_neutrals, only: neutral_beam_t
    use tiago_fusion_rates, only: SPEED_OF_LIGHT
    use tiago_fida, only: fida_t, D_ALPHA_NM
    use fida_test_providers, only: slab_plasma_t, tilted_field
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: failures

    failures = 0
    call check_gyro_distribution()
    call check_isotropic_box()
    call check_weights_and_derivatives()
    if (failures > 0) then
        write (error_unit, '(I0,A)') failures, ' FIDA checks failed'
        error stop 1
    end if
    print '(A)', 'FIDA checks passed'

contains

    subroutine expect(ok, what, value)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        real(dp), intent(in) :: value
        if (.not. ok) then
            failures = failures + 1
            write (error_unit, '(A,A,ES12.4)') 'FAIL ', what, value
        else
            print '(A,A,ES12.4)', 'ok   ', what, value
        end if
    end subroutine expect

    function one_chord() result(chords)
        type(chord_t) :: chords(1)
        chords(1)%origin = [1.5_dp, 0.0_dp, 0.0_dp]
        chords(1)%direction = [1.0_dp, 0.0_dp, 0.0_dp]
        chords(1)%length = 0.2_dp
        chords(1)%npanel = 2
        chords(1)%ngauss = 2
    end function one_chord

    function cell_distribution(energy, pitch, np) result(dist)
        !! One energy cell; np pitch cells on [-1, 1] if pitch is outside it,
        !! otherwise one cell at pitch. Constant in s.
        real(dp), intent(in) :: energy, pitch
        integer, intent(in) :: np
        type(fast_ion_distribution_t) :: dist
        integer :: j
        dist%s = [0.0_dp, 2.0_dp]
        dist%energy = [energy]
        dist%denergy = [1.0e3_dp]
        if (abs(pitch) <= 1.0_dp) then
            dist%pitch = [pitch]
            dist%dpitch = [0.01_dp]
        else
            dist%dpitch = [(2.0_dp / np, j = 1, np)]
            dist%pitch = [(-1.0_dp + (j - 0.5_dp) * 2.0_dp / np, j = 1, np)]
        end if
        allocate(dist%f(1, size(dist%pitch), 2))
        dist%f = 1.0e12_dp
    end function cell_distribution

    function uniform_sample(x) result(sample)
        real(dp), intent(in) :: x(:, :)
        type(plasma_sample_t) :: sample
        integer :: m
        sample%n = size(x, 2)
        allocate(sample%inside(sample%n), sample%q(N_QUANTITIES, sample%n), &
            sample%dq_dx(3, N_QUANTITIES, sample%n))
        sample%x = x
        sample%inside = .true.
        sample%dq_dx = 0.0_dp
        do m = 1, sample%n
            sample%q(:, m) = [0.3_dp + 0.05_dp * m, 2.5_dp, 4.0e19_dp, 2.0e3_dp, &
                1.5e3_dp + 100.0_dp * m, 1.6_dp]
        end do
    end function uniform_sample

    subroutine check_gyro_distribution()
        type(fida_t) :: fida
        type(plasma_sample_t) :: sample
        real(dp), allocatable :: x(:, :), dx(:, :, :), y(:), edges(:), expected(:)
        real(dp), parameter :: energy = 40.0e3_dp, pitch = 0.3_dp
        real(dp) :: v, a, u1, u2, width, err
        integer :: status, k, nb

        nb = 60
        v = sqrt(2.0_dp * energy * 1.602176634e-19_dp / (2.01410177812_dp * 1.66053906660e-27_dp))
        width = D_ALPHA_NM * v / SPEED_OF_LIGHT
        edges = [(D_ALPHA_NM - 1.05_dp * width + 2.1_dp * width * k / nb, k = 0, nb)]
        fida%ngyro = 20000
        call fida%init(one_chord(), edges, cell_distribution(energy, pitch, 1), tilted_field, status)
        call expect(status == 0, 'FIDA init status', real(status, dp))
        call fida%sample_points([real(dp) ::], x, dx)
        call halo_only(fida, size(x, 2))
        sample = uniform_sample(x)
        sample%q(Q_TI, :) = 1.5e3_dp
        allocate(y(fida%n_signals()), expected(nb))
        call fida%predict(sample, [real(dp) ::], y)
        a = 50.0_dp * pi / 180.0_dp
        do k = 1, nb
            u1 = SPEED_OF_LIGHT * (edges(k) / D_ALPHA_NM - 1.0_dp)
            u2 = SPEED_OF_LIGHT * (edges(k + 1) / D_ALPHA_NM - 1.0_dp)
            expected(k) = (acos(clip((u1 - v * pitch * cos(a)) / (v * sqrt(1 - pitch**2) * sin(a)))) &
                - acos(clip((u2 - v * pitch * cos(a)) / (v * sqrt(1 - pitch**2) * sin(a))))) / pi &
                / (edges(k + 1) - edges(k))
        end do
        y = y / sum(y * (edges(2:) - edges(:nb)))
        err = maxval(abs(y - expected)) / maxval(expected)
        call expect(err < 2.0e-3_dp, 'gyro-phase spectrum against the closed form, max rel. error', err)
    end subroutine check_gyro_distribution

    subroutine halo_only(fida, n)
        !! Uniform halo of 1e16 m^-3 and no beam neutrals.
        type(fida_t), intent(inout) :: fida
        integer, intent(in) :: n
        real(dp) :: nb(0, n), vb(3, 0, n), nh(n)
        nh = 1.0e16_dp
        call fida%set_neutrals(nb, vb, nh)
    end subroutine halo_only

    pure real(dp) function clip(x)
        real(dp), intent(in) :: x
        clip = max(-1.0_dp, min(1.0_dp, x))
    end function clip

    subroutine check_isotropic_box()
        type(fida_t) :: fida
        type(plasma_sample_t) :: sample
        real(dp), allocatable :: x(:, :), dx(:, :, :), y(:), edges(:)
        real(dp), parameter :: energy = 60.0e3_dp
        real(dp) :: v, width, lo, hi, mean, spread, outside
        integer :: status, k, nb, nin

        nb = 40
        v = sqrt(2.0_dp * energy * 1.602176634e-19_dp / (2.01410177812_dp * 1.66053906660e-27_dp))
        width = D_ALPHA_NM * v / SPEED_OF_LIGHT
        edges = [(D_ALPHA_NM - 1.25_dp * width + 2.5_dp * width * k / nb, k = 0, nb)]
        fida%ngyro = 400
        call fida%init(one_chord(), edges, cell_distribution(energy, 2.0_dp, 400), tilted_field, &
            status)
        call fida%sample_points([real(dp) ::], x, dx)
        call halo_only(fida, size(x, 2))
        sample = uniform_sample(x)
        allocate(y(nb))
        call fida%predict(sample, [real(dp) ::], y)
        mean = 0.0_dp
        spread = 0.0_dp
        outside = 0.0_dp
        nin = 0
        lo = D_ALPHA_NM - width
        hi = D_ALPHA_NM + width
        do k = 1, nb
            if (edges(k) >= lo .and. edges(k + 1) <= hi) then
                nin = nin + 1
                mean = mean + y(k)
            else if (edges(k + 1) <= lo .or. edges(k) >= hi) then
                outside = max(outside, y(k))
            end if
        end do
        mean = mean / nin
        do k = 1, nb
            if (edges(k) >= lo .and. edges(k + 1) <= hi) spread = max(spread, abs(y(k) / mean - 1))
        end do
        call expect(spread < 1.0e-2_dp, 'isotropic population: flat Doppler box, max deviation', spread)
        call expect(outside == 0.0_dp, 'isotropic population: no emission outside the box', outside)
    end subroutine check_isotropic_box

    subroutine check_weights_and_derivatives()
        type(fida_t) :: fida
        type(neutral_beam_t) :: beam(1)
        type(slab_plasma_t) :: plasma
        type(fast_ion_distribution_t) :: dist, uniform
        type(plasma_sample_t) :: sample
        real(dp), allocatable :: x(:, :), dx(:, :, :), y(:), edges(:), w(:, :, :), eta(:)
        real(dp) :: s_from_w, jvp_error, adjoint_error, err
        integer :: status, i, j, k, nb, ne, np, ns

        beam(1)%origin = [1.6_dp, -2.0_dp, 0.0_dp]
        beam(1)%direction = [0.0_dp, 1.0_dp, 0.0_dp]
        beam(1)%width = [0.12_dp, 0.12_dp]
        beam(1)%energy = 60.0e3_dp
        beam(1)%power = 2.5e6_dp
        beam(1)%power_fraction = [0.6_dp, 0.3_dp, 0.1_dp]
        ne = 6
        np = 5
        ns = 4
        allocate(dist%s(ns), dist%energy(ne), dist%denergy(ne), dist%pitch(np), dist%dpitch(np), &
            dist%f(ne, np, ns))
        dist%s = [0.0_dp, 0.3_dp, 0.7_dp, 1.0_dp]
        dist%denergy = 10.0e3_dp
        dist%energy = [(15.0e3_dp + 10.0e3_dp * (i - 1), i = 1, ne)]
        dist%dpitch = 0.4_dp
        dist%pitch = [(-0.8_dp + 0.4_dp * (j - 1), j = 1, np)]
        do k = 1, ns
            do j = 1, np
                do i = 1, ne
                    dist%f(i, j, k) = 1.0e14_dp * (1 - dist%s(k)) * (1 + 0.1_dp * i) &
                        * (1.2_dp + dist%pitch(j))
                end do
            end do
        end do
        nb = 30
        edges = [(649.0_dp + 14.0_dp * k / nb, k = 0, nb)]

        uniform = dist
        uniform%s = [0.0_dp, 2.0_dp]
        uniform%f = dist%f(:, :, 1:2)
        uniform%f(:, :, 2) = uniform%f(:, :, 1)
        call fida%init(one_chord(), edges, uniform, tilted_field, status)
        call fida%sample_points([real(dp) ::], x, dx)
        call fida%set_beams(beam, plasma, status)
        call expect(status == 0, 'FIDA beam neutrals status', real(status, dp))
        sample = uniform_sample(x)
        allocate(y(nb), w(nb, ne, np), eta(0))
        call fida%predict(sample, eta, y)
        call fida%weights(sample, 1, w)
        err = 0.0_dp
        do k = 1, nb
            s_from_w = 0.0_dp
            do j = 1, np
                do i = 1, ne
                    s_from_w = s_from_w + w(k, i, j) * uniform%f(i, j, 1) * dist%denergy(i) * dist%dpitch(j)
                end do
            end do
            err = max(err, abs(s_from_w - y(k)) / maxval(abs(y)))
        end do
        call expect(maxval(y) > 0.0_dp .and. err < 1.0e-12_dp, &
            'weight functions reproduce the spectrum, max rel. error', err)

        call fida%init(one_chord(), edges, dist, tilted_field, status)
        call fida%set_beams(beam, plasma, status)
        call check_diagnostic_derivatives(fida, sample, eta, 1.0e-6_dp, jvp_error, adjoint_error)
        call expect(jvp_error < 1.0e-6_dp, 'FIDA jvp against central differences', jvp_error)
        call expect(adjoint_error < 1.0e-12_dp, 'FIDA vjp transposition', adjoint_error)
    end subroutine check_weights_and_derivatives

end program test_fida
