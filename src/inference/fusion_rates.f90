module tiago_fusion_rates
    !! Reaction cross-sections in SI units and their average over isotropic
    !! Maxwellian targets at rest. For a test particle of speed v in targets of
    !! mass m_t and temperature T, with thermal speed w = sqrt(2 T / m_t),
    !!
    !!   <sigma u>(v, T) = 1 / (sqrt(pi) v) int g(v + w t) exp(-t^2) dt,
    !!   g(u) = sgn(u) u^2 sigma(|u|),
    !!
    !! which is the angle-integrated form of int sigma(|v - v_t|) |v - v_t| f_M d^3v_t.
    !! The integral is a fixed Gauss-Legendre rule on |t| <= 7 (exp(-49)
    !! neglected). Reactions: D(d,n)3He (Bosch and Hale 1992) and H+ + H(1s)
    !! charge exchange (Janev and Smith 1993), kernels in src/generated.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_sight_geometry, only: gauss_legendre01
    use tiago_dd_cross_section, only: dd_cross_section
    use tiago_dd_reactivity, only: dd_reactivity
    use tiago_cx_hydrogen, only: cx_hydrogen
    implicit none
    private

    real(dp), parameter, public :: ATOMIC_MASS = 1.66053906660e-27_dp !! [kg]
    real(dp), parameter, public :: DEUTERON_MASS = 2.01410177812_dp * ATOMIC_MASS !! atom [kg]
    real(dp), parameter, public :: ELEMENTARY_CHARGE = 1.602176634e-19_dp !! [C], J per eV
    real(dp), parameter, public :: SPEED_OF_LIGHT = 299792458.0_dp !! [m/s]
    integer, parameter, public :: REACTION_DD_N = 1 !! D(d,n)3He
    integer, parameter, public :: REACTION_CX_H = 2 !! H+ + H(1s) charge exchange

    real(dp), parameter :: pi = acos(-1.0_dp), t_max = 7.0_dp
    real(dp), parameter :: kev = 1.0e3_dp * ELEMENTARY_CHARGE

    type, public :: thermal_rule_t
        !! Nodes t and weights exp(-t^2) dt / sqrt(pi) of the thermal average.
        real(dp), allocatable :: t(:), w(:)
    end type thermal_rule_t

    public :: thermal_rule, reaction_sigma, maxwellian_rate, dd_thermal_reactivity

contains

    function thermal_rule(n) result(rule)
        integer, intent(in) :: n
        type(thermal_rule_t) :: rule
        real(dp), allocatable :: x(:), w(:)
        call gauss_legendre01(n, x, w)
        rule%t = t_max * (2.0_dp * x - 1.0_dp)
        rule%w = 2.0_dp * t_max * w * exp(-rule%t**2) / sqrt(pi)
    end function thermal_rule

    pure subroutine reaction_sigma(kind, u, sigma, dsigma_du)
        !! Cross-section [m^2] at the relative speed u >= 0 [m/s] and its u derivative.
        !! D-D: centre-of-mass energy m_D u^2 / 4, zero below 0.05 keV, where
        !! the Gamow factor is below 1e-55. Charge exchange: collision energy
        !! per amu, m_u u^2 / 2, floored at 1e-8 keV/amu where u^2 sigma vanishes.
        integer, intent(in) :: kind
        real(dp), intent(in) :: u
        real(dp), intent(out) :: sigma, dsigma_du
        real(dp) :: e, de_du, s, ds

        sigma = 0.0_dp
        dsigma_du = 0.0_dp
        select case (kind)
        case (REACTION_DD_N)
            e = 0.25_dp * DEUTERON_MASS * u**2 / kev
            if (e < 0.05_dp) return
            de_du = 0.5_dp * DEUTERON_MASS * u / kev
            call dd_cross_section(e, s, ds)
            sigma = 1.0e-31_dp * s
            dsigma_du = 1.0e-31_dp * ds * de_du
        case (REACTION_CX_H)
            e = 0.5_dp * ATOMIC_MASS * u**2 / kev
            if (e < 1.0e-8_dp) return
            de_du = ATOMIC_MASS * u / kev
            call cx_hydrogen(e, s, ds)
            sigma = 1.0e-4_dp * s
            dsigma_du = 1.0e-4_dp * ds * de_du
        end select
    end subroutine reaction_sigma

    pure subroutine maxwellian_rate(rule, kind, v, temperature, mass, rate, drate_dt)
        !! <sigma u> [m^3/s] of a test particle of speed v > 0 [m/s] in targets
        !! of mass [kg] at temperature [eV], and its temperature derivative.
        !! A non-positive temperature gives the cold-target rate sigma(v) v.
        type(thermal_rule_t), intent(in) :: rule
        integer, intent(in) :: kind
        real(dp), intent(in) :: v, temperature, mass
        real(dp), intent(out) :: rate, drate_dt
        real(dp) :: w, u, sg, dsg, g, dg
        integer :: i

        drate_dt = 0.0_dp
        if (.not. temperature > 0.0_dp) then
            call reaction_sigma(kind, v, sg, dsg)
            rate = sg * v
            return
        end if
        w = sqrt(2.0_dp * temperature * ELEMENTARY_CHARGE / mass)
        rate = 0.0_dp
        do i = 1, size(rule%t)
            u = v + w * rule%t(i)
            call reaction_sigma(kind, abs(u), sg, dsg)
            g = sign(1.0_dp, u) * u**2 * sg
            dg = 2.0_dp * abs(u) * sg + u**2 * dsg
            rate = rate + rule%w(i) * g
            drate_dt = drate_dt + rule%w(i) * rule%t(i) * dg
        end do
        rate = rate / v
        ! d/dT = (d/dw) (w / 2T)
        drate_dt = drate_dt / v * w / (2.0_dp * temperature)
    end subroutine maxwellian_rate

    pure subroutine dd_thermal_reactivity(temperature, rate, drate_dt)
        !! Bosch-Hale D(d,n)3He reactivity [m^3/s] at the ion temperature [eV]
        !! and its derivative [m^3/s/eV]; zero below 1 eV.
        real(dp), intent(in) :: temperature
        real(dp), intent(out) :: rate, drate_dt
        rate = 0.0_dp
        drate_dt = 0.0_dp
        if (temperature < 1.0_dp) return
        call dd_reactivity(1.0e-3_dp * temperature, rate, drate_dt)
        rate = 1.0e-6_dp * rate
        drate_dt = 1.0e-9_dp * drate_dt
    end subroutine dd_thermal_reactivity

end module tiago_fusion_rates
