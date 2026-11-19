program test_fast_ion_neutrons
    !! D-D neutron emission (tiago_neutron_rate, tiago_fusion_rates).
    !! Oracles: the Bosch-Hale cross-section and reactivity tables (Nucl.
    !! Fusion 32, 611 (1992), Tables V and VIII), which come from the R-matrix
    !! cross-sections, not from the fits evaluated here; a fast population that
    !! is itself a Maxwellian at T_i, for which the beam-target and beam-beam
    !! rates must equal n_f n_D <sigma v>(T_i) and n_f^2 <sigma v>(T_i) / 2 with
    !! the tabulated reactivity; the cold-target limit sigma(E / 2) v; the
    !! volume of the R-Z rule; the contract's derivative check.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_profile_diagnostic, only: plasma_sample_t, Q_S, Q_B, Q_NE, Q_TE, Q_TI, Q_ZEFF, &
        N_QUANTITIES, check_diagnostic_derivatives
    use tiago_fast_ion_distribution, only: fast_ion_distribution_t, read_fast_ion_distribution
    use tiago_fusion_rates, only: thermal_rule_t, thermal_rule, maxwellian_rate, reaction_sigma, &
        dd_thermal_reactivity, REACTION_DD_N, DEUTERON_MASS, ELEMENTARY_CHARGE
    use tiago_dd_cross_section, only: dd_cross_section
    use tiago_neutron_rate, only: neutron_t, neutron_volume_points
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: failures

    failures = 0
    call check_tables()
    call check_cold_target()
    call check_maxwellian_population()
    call check_volume_rule()
    call check_derivatives()
    call check_file_roundtrip()
    if (failures > 0) then
        write (error_unit, '(I0,A)') failures, ' fast-ion neutron checks failed'
        error stop 1
    end if
    print '(A)', 'fast-ion neutron checks passed'

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

    subroutine check_tables()
        ! Table V, D(d,n)3He [mb] at E_cm [keV]; Table VIII [cm^3/s] at T [keV].
        real(dp), parameter :: e(4) = [10.0_dp, 30.0_dp, 100.0_dp, 400.0_dp]
        real(dp), parameter :: sigma(4) = [2.779e-1_dp, 6.857_dp, 3.701e1_dp, 8.702e1_dp]
        real(dp), parameter :: t(6) = [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, 10.0_dp, 20.0_dp]
        real(dp), parameter :: sv(6) = [9.933e-23_dp, 3.110e-21_dp, 1.602e-20_dp, 9.128e-20_dp, &
            6.023e-19_dp, 2.603e-18_dp]
        real(dp) :: s, ds, r, dr, worst
        integer :: i
        worst = 0.0_dp
        do i = 1, size(e)
            call dd_cross_section(e(i), s, ds)
            worst = max(worst, abs(s / sigma(i) - 1.0_dp))
        end do
        call expect(worst < 1.0e-3_dp, 'cross-section against Table V, max rel. error', worst)
        worst = 0.0_dp
        do i = 1, size(t)
            call dd_thermal_reactivity(1.0e3_dp * t(i), r, dr)
            worst = max(worst, abs(1.0e6_dp * r / sv(i) - 1.0_dp))
        end do
        ! Bosch and Hale quote a maximum fit deviation of 0.3 % for this reaction.
        call expect(worst < 3.5e-3_dp, 'reactivity fit against Table VIII, max rel. error', worst)
    end subroutine check_tables

    subroutine check_cold_target()
        type(thermal_rule_t) :: rule
        real(dp) :: v, r, dr, sg, dsg
        rule = thermal_rule(96)
        v = sqrt(2.0_dp * 80.0e3_dp * ELEMENTARY_CHARGE / DEUTERON_MASS)
        call maxwellian_rate(rule, REACTION_DD_N, v, 0.01_dp, DEUTERON_MASS, r, dr)
        call reaction_sigma(REACTION_DD_N, v, sg, dsg)
        call expect(abs(r / (sg * v) - 1.0_dp) < 1.0e-6_dp, &
            'beam-target at T_i = 0.01 eV against sigma(E/2) v, rel. error', abs(r / (sg * v) - 1.0_dp))
    end subroutine check_cold_target

    function maxwellian_distribution(t_ev, n, ne, np, emax) result(dist)
        !! Isotropic Maxwellian, constant in s on two nodes.
        real(dp), intent(in) :: t_ev, n, emax
        integer, intent(in) :: ne, np
        type(fast_ion_distribution_t) :: dist
        integer :: i, j
        dist%s = [0.0_dp, 1.0_dp]
        allocate(dist%energy(ne), dist%denergy(ne), dist%pitch(np), dist%dpitch(np), &
            dist%f(ne, np, 2))
        dist%denergy = emax / ne
        dist%energy = [((i - 0.5_dp) * emax / ne, i = 1, ne)]
        dist%dpitch = 2.0_dp / np
        dist%pitch = [(-1.0_dp + (j - 0.5_dp) * 2.0_dp / np, j = 1, np)]
        do j = 1, np
            dist%f(:, j, 1) = n / sqrt(pi) * sqrt(dist%energy) * t_ev**(-1.5_dp) &
                * exp(-dist%energy / t_ev)
        end do
        dist%f(:, :, 2) = dist%f(:, :, 1)
    end function maxwellian_distribution

    subroutine check_maxwellian_population()
        type(fast_ion_distribution_t) :: dist
        type(neutron_t) :: model
        real(dp) :: q(N_QUANTITIES), eps(3), deps(3, N_QUANTITIES), nf, nd, sv, dsv, err
        real(dp), parameter :: t_ev = 5.0e3_dp
        integer :: k

        dist = maxwellian_distribution(t_ev, 1.0e19_dp, 320, 12, 30.0_dp * t_ev)
        nf = 0.0_dp
        do k = 1, size(dist%pitch)
            nf = nf + sum(dist%f(:, k, 1) * dist%denergy) * dist%dpitch(k)
        end do
        call model%setup(dist, z_imp=7.0_dp)
        q = 0.0_dp
        q(Q_S) = 0.5_dp
        q(Q_NE) = 4.0e19_dp
        q(Q_TI) = t_ev
        q(Q_ZEFF) = 1.6_dp
        nd = q(Q_NE) * (7.0_dp - q(Q_ZEFF)) / 6.0_dp
        call model%emissivity(q, eps, deps)
        call dd_thermal_reactivity(t_ev, sv, dsv)
        call expect(abs(eps(1) / (0.5_dp * nd**2 * sv) - 1.0_dp) < 1.0e-14_dp, &
            'thermonuclear emissivity n_D^2 <sigma v> / 2', eps(1))
        err = abs(eps(2) / (nf * nd * sv) - 1.0_dp)
        call expect(err < 5.0e-3_dp, 'beam-target of a Maxwellian population against Bosch-Hale', err)
        err = abs(eps(3) / (0.5_dp * nf**2 * sv) - 1.0_dp)
        call expect(err < 5.0e-3_dp, 'beam-beam of a Maxwellian population against Bosch-Hale', err)
    end subroutine check_maxwellian_population

    subroutine check_volume_rule()
        real(dp), allocatable :: points(:, :), weights(:)
        integer, allocatable :: first(:)
        real(dp) :: exact
        call neutron_volume_points(1.0_dp, 2.0_dp, -0.5_dp, 0.5_dp, 7, 3, points, weights, first)
        exact = pi * (2.0_dp**2 - 1.0_dp**2) * 1.0_dp
        call expect(abs(sum(weights) / exact - 1.0_dp) < 1.0e-13_dp, 'volume of the R-Z rule', &
            sum(weights))
    end subroutine check_volume_rule

    function slowing_down(ne, np, ns) result(dist)
        !! Beam-like test distribution: classical slowing-down in E, Gaussian
        !! in pitch about 0.6, falling to zero at s = 1.
        integer, intent(in) :: ne, np, ns
        type(fast_ion_distribution_t) :: dist
        integer :: i, j, k
        allocate(dist%s(ns), dist%energy(ne), dist%denergy(ne), dist%pitch(np), dist%dpitch(np), &
            dist%f(ne, np, ns))
        dist%s = [((k - 1.0_dp) / (ns - 1), k = 1, ns)]
        dist%denergy = 80.0e3_dp / ne
        dist%energy = [(5.0e3_dp + (i - 0.5_dp) * 80.0e3_dp / ne, i = 1, ne)]
        dist%dpitch = 2.0_dp / np
        dist%pitch = [(-1.0_dp + (j - 0.5_dp) * 2.0_dp / np, j = 1, np)]
        do k = 1, ns
            do j = 1, np
                do i = 1, ne
                    dist%f(i, j, k) = 1.0e14_dp * (1.0_dp - dist%s(k)) * (1.0_dp + dist%s(k)) &
                        * sqrt(dist%energy(i)) / (dist%energy(i)**1.5_dp + 2.0e4_dp**1.5_dp) &
                        * exp(-(dist%pitch(j) - 0.6_dp)**2 / 0.2_dp)
                end do
            end do
        end do
    end function slowing_down

    subroutine check_derivatives()
        type(neutron_t) :: model
        type(plasma_sample_t) :: sample
        real(dp), allocatable :: points(:, :), weights(:), eta(:)
        integer, allocatable :: first(:)
        real(dp) :: jvp_error, adjoint_error
        integer :: m

        call neutron_volume_points(1.2_dp, 2.0_dp, -0.4_dp, 0.4_dp, 4, 3, points, weights, first)
        call model%init(slowing_down(10, 6, 5), points, weights, [1, 7, 13])
        sample%n = size(weights)
        allocate(sample%x(3, sample%n), sample%inside(sample%n), &
            sample%q(N_QUANTITIES, sample%n), sample%dq_dx(3, N_QUANTITIES, sample%n), eta(0))
        sample%x = points
        sample%inside = .true.
        sample%dq_dx = 0.0_dp
        do m = 1, sample%n
            sample%q(Q_S, m) = 0.07_dp + 0.071_dp * m
            sample%q(Q_B, m) = 2.5_dp
            sample%q(Q_NE, m) = 5.0e19_dp * (1.0_dp - 0.5_dp * sample%q(Q_S, m))
            sample%q(Q_TE, m) = 3.0e3_dp
            sample%q(Q_TI, m) = 2.5e3_dp * (1.0_dp - 0.6_dp * sample%q(Q_S, m))
            sample%q(Q_ZEFF, m) = 1.5_dp + 0.1_dp * sample%q(Q_S, m)
        end do
        call check_diagnostic_derivatives(model, sample, eta, 1.0e-6_dp, jvp_error, adjoint_error)
        call expect(jvp_error < 1.0e-6_dp, 'neutron jvp against central differences', jvp_error)
        call expect(adjoint_error < 1.0e-12_dp, 'neutron vjp transposition', adjoint_error)
    end subroutine check_derivatives

    subroutine check_file_roundtrip()
        type(fast_ion_distribution_t) :: dist, back
        integer :: unit, status
        character(len=*), parameter :: path = 'fast_ion_distribution_test.txt'
        dist = slowing_down(4, 3, 3)
        open (newunit=unit, file=path, status='replace')
        write (unit, '(A)') '# test distribution'
        write (unit, *) dist%mass_amu
        write (unit, *) size(dist%s), size(dist%energy), size(dist%pitch)
        write (unit, *) dist%s
        write (unit, *) dist%energy, dist%denergy
        write (unit, *) dist%pitch, dist%dpitch
        write (unit, *) dist%f
        close (unit)
        call read_fast_ion_distribution(path, back, status)
        call expect(status == 0, 'distribution file read status', real(status, dp))
        if (status == 0) call expect(maxval(abs(back%f - dist%f)) <= 1.0e-12_dp * maxval(dist%f), &
            'distribution file round trip', maxval(abs(back%f - dist%f)))
        open (newunit=unit, file=path, status='old')
        close (unit, status='delete')
    end subroutine check_file_roundtrip

end program test_fast_ion_neutrons
