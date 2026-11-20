module tiago_geqdsk_state
    !! Provider side of the profile-diagnostic contract for a fixed
    !! axisymmetric equilibrium from a G-EQDSK file: the rows Q_S and Q_B of a
    !! plasma_sample_t and their spatial derivatives at points (R, phi, Z).
    !! s = Phi(psi_N) / Phi(1) is the normalized toroidal flux of
    !! tiago_geqdsk_plasma, continued linearly in psi_N outside the plasma, and
    !! |B| = sqrt(|grad psi|^2 + F^2) / R. Points are located directly from
    !! psi_N(R, Z), so no Newton inversion is needed. A point is inside when it
    !! lies inside the last closed surface (tiago_geqdsk_plasma%inside).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_geqdsk_plasma, only: geqdsk_plasma_t
    use tiago_profile_diagnostic, only: plasma_sample_t, Q_S, Q_B, N_QUANTITIES
    implicit none
    private

    real(dp), parameter :: fd_step = 1.0e-5_dp !! [m], central differences of |B|

    public :: sample_geqdsk

contains

    subroutine sample_geqdsk(eq, x, sample)
        type(geqdsk_plasma_t), intent(in) :: eq
        real(dp), intent(in) :: x(:, :) !! (3, n): R, phi, Z
        type(plasma_sample_t), intent(out) :: sample
        real(dp) :: r, z, psin, a, b, s, ds, br, bz
        integer :: k, n
        n = size(x, 2)
        sample%n = n
        sample%x = x
        allocate(sample%inside(n), sample%q(N_QUANTITIES, n), sample%dq_dx(3, N_QUANTITIES, n))
        sample%q = 0
        sample%dq_dx = 0
        do k = 1, n
            r = x(1, k)
            z = x(3, k)
            call eq%label(r, z, psin, a, b)
            call eq%toroidal_flux(psin, s, ds)
            sample%q(Q_S, k) = s
            sample%dq_dx(:, Q_S, k) = [ds * a, 0.0_dp, ds * b]
            sample%q(Q_B, k) = eq%modb(r, z)
            br = (eq%modb(r + fd_step, z) - eq%modb(r - fd_step, z)) / (2 * fd_step)
            bz = (eq%modb(r, z + fd_step) - eq%modb(r, z - fd_step)) / (2 * fd_step)
            sample%dq_dx(:, Q_B, k) = [br, 0.0_dp, bz]
            sample%inside(k) = eq%inside(r, z)
        end do
    end subroutine sample_geqdsk

end module tiago_geqdsk_state
