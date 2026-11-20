module tiago_sauter_bootstrap
    !! Bootstrap current of Sauter, Angioni and Lin-Liu, Phys. Plasmas 6, 2834
    !! (1999), with the erratum Phys. Plasmas 9, 5140 (2002) (docs/bootstrap.md).
    !! Densities [m^-3], temperatures [eV], R [m], psi [Wb/rad] in TIAGO's
    !! convention B_p = grad psi x grad phi, in which Sauter's psi is -psi:
    !!     <j.B>_bs = F p_e [L31 p/p_e dln p/dpsi + L32 dln T_e/dpsi
    !!                       + L34 alpha (1 - R_pe)/R_pe dln T_i/dpsi],   R_pe = p_e / p.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    real(dp), parameter :: e_charge = 1.602176634e-19_dp

    public :: sauter_l31, sauter_l32, sauter_l34, sauter_alpha
    public :: sauter_nu_e, sauter_nu_i, sauter_jdotb

contains

    elemental real(dp) function sauter_l31(ft, nue, z) result(l)
        !! Eqs. 14a, b.
        real(dp), intent(in) :: ft, nue, z
        real(dp) :: x
        x = ft / (1 + (1 - 0.1_dp * ft) * sqrt(nue) + 0.5_dp * (1 - ft) * nue / z)
        l = (1 + 1.4_dp / (z + 1)) * x - 1.9_dp / (z + 1) * x**2 &
            + 0.3_dp / (z + 1) * x**3 + 0.2_dp / (z + 1) * x**4
    end function sauter_l31

    elemental real(dp) function sauter_l32(ft, nue, z) result(l)
        !! Eqs. 15a-e.
        real(dp), intent(in) :: ft, nue, z
        real(dp) :: xe, xi
        xe = ft / (1 + 0.26_dp * (1 - ft) * sqrt(nue) &
            + 0.18_dp * (1 - 0.37_dp * ft) * nue / sqrt(z))
        xi = ft / (1 + (1 + 0.6_dp * ft) * sqrt(nue) &
            + 0.85_dp * (1 - 0.37_dp * ft) * nue * (1 + z))
        l = (0.05_dp + 0.62_dp * z) / (z * (1 + 0.44_dp * z)) * (xe - xe**4) &
            + 1 / (1 + 0.22_dp * z) * (xe**2 - xe**4 - 1.2_dp * (xe**3 - xe**4)) &
            + 1.2_dp / (1 + 0.5_dp * z) * xe**4 &
            - (0.56_dp + 1.93_dp * z) / (z * (1 + 0.44_dp * z)) * (xi - xi**4) &
            + 4.95_dp / (1 + 2.48_dp * z) * (xi**2 - xi**4 - 0.55_dp * (xi**3 - xi**4)) &
            - 1.2_dp / (1 + 0.5_dp * z) * xi**4
    end function sauter_l32

    elemental real(dp) function sauter_l34(ft, nue, z) result(l)
        !! Eqs. 16a, b.
        real(dp), intent(in) :: ft, nue, z
        real(dp) :: x
        x = ft / (1 + (1 - 0.1_dp * ft) * sqrt(nue) + 0.5_dp * (1 - 0.5_dp * ft) * nue / z)
        l = (1 + 1.4_dp / (z + 1)) * x - 1.9_dp / (z + 1) * x**2 &
            + 0.3_dp / (z + 1) * x**3 + 0.2_dp / (z + 1) * x**4
    end function sauter_l34

    elemental real(dp) function sauter_alpha(ft, nui) result(a)
        !! Eqs. 17a, b.
        real(dp), intent(in) :: ft, nui
        real(dp) :: a0
        a0 = -1.17_dp * (1 - ft) / (1 - 0.22_dp * ft - 0.19_dp * ft**2)
        a = ((a0 + 0.25_dp * (1 - ft**2) * sqrt(nui)) / (1 + 0.5_dp * sqrt(nui)) &
            + 0.315_dp * nui**2 * ft**6) / (1 + 0.15_dp * nui**2 * ft**6)
    end function sauter_alpha

    elemental real(dp) function sauter_nu_e(q, r, eps, ne, te, z) result(nu)
        !! Eq. 18b with ln Lambda_e of Eq. 18d.
        real(dp), intent(in) :: q, r, eps, ne, te, z
        nu = 6.921e-18_dp * q * r * ne * z * (31.3_dp - log(sqrt(ne) / te)) / (te**2 * eps**1.5_dp)
    end function sauter_nu_e

    elemental real(dp) function sauter_nu_i(q, r, eps, ni, ti, zi) result(nu)
        !! Eq. 18c with the ion charge (erratum) and ln Lambda_ii of Eq. 18e.
        real(dp), intent(in) :: q, r, eps, ni, ti, zi
        nu = 4.90e-18_dp * q * r * ni * zi**4 * (30.0_dp - log(zi**3 * sqrt(ni) / ti**1.5_dp)) &
            / (ti**2 * eps**1.5_dp)
    end function sauter_nu_i

    elemental real(dp) function sauter_jdotb(f, ft, q, r, eps, zeff, zi, ne, ni, te, ti, &
            dne, dni, dte, dti) result(jb)
        !! <j.B>_bs [A T/m^2]; dne ... dti are d/dpsi of n_e, n_i, T_e, T_i
        !! [per Wb/rad], f = F [T m].
        real(dp), intent(in) :: f, ft, q, r, eps, zeff, zi, ne, ni, te, ti, dne, dni, dte, dti
        real(dp) :: pe, pi_, nue, nui, l31, l32, l34
        pe = e_charge * ne * te
        pi_ = e_charge * ni * ti
        nue = sauter_nu_e(q, r, eps, ne, te, zeff)
        nui = sauter_nu_i(q, r, eps, ni, ti, zi)
        l31 = sauter_l31(ft, nue, zeff)
        l32 = sauter_l32(ft, nue, zeff)
        l34 = sauter_l34(ft, nue, zeff)
        jb = f * (l31 * (pe * dne / ne + pi_ * dni / ni) + (l31 + l32) * pe * dte / te &
            + (l31 + sauter_alpha(ft, nui) * l34) * pi_ * dti / ti)
    end function sauter_jdotb

end module tiago_sauter_bootstrap
