program plasma_surface_fixture
    !! Analytic circular torus (R0 = 3 m, a = 0.5 m) with total boundary field
    !! B = Bp e_theta + B0 R0/R e_phi, as a --plasma-surface file, plus B-probe
    !! chains whose summed signals are loop integrals of B with closed forms:
    !!   poloidal circle of radius 0.9 m around the plasma: 2 pi a Bp (Ampere),
    !!   toroidal circle at R0 inside the boundary: -2 pi R0 B0 (minus coil field),
    !!   toroidal circle at 4 m outside: 0 (no net poloidal current linked).
    !! Usage: plasma_surface_fixture <dir> [tilted]; "tilted" adds a normal field.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 3.0_dp, a = 0.5_dp
    real(dp), parameter :: bp = 0.1_dp, b0 = 2.0_dp
    integer, parameter :: ntheta = 64, nphi = 256, nloop = 400
    character(len=512) :: dir, mode
    real(dp) :: th, ph, r, xt(3), xp(3), x(3), b(3), t, dl
    integer :: unit, j, k

    call get_command_argument(1, dir)
    call get_command_argument(2, mode)
    open(newunit=unit, file=trim(dir)//'/surface.txt', status='replace', action='write')
    write(unit, '(2I6)') ntheta, nphi
    do k = 0, nphi - 1
        ph = 2.0_dp * pi * real(k, dp) / real(nphi, dp)
        do j = 0, ntheta - 1
            th = 2.0_dp * pi * real(j, dp) / real(ntheta, dp)
            r = r0 + a * cos(th)
            x = [r * cos(ph), r * sin(ph), a * sin(th)]
            xt = [-a * sin(th) * cos(ph), -a * sin(th) * sin(ph), a * cos(th)]
            xp = [-r * sin(ph), r * cos(ph), 0.0_dp]
            b = bp * xt / a + b0 * r0 / r * xp / r
            if (trim(mode) == 'tilted') b = b + 0.01_dp * [cos(th) * cos(ph), cos(th) * sin(ph), sin(th)]
            if (trim(mode) == 'nonfinite') b(1) = ieee_value(0.0_dp, ieee_quiet_nan)
            if (trim(mode) == 'degenerate') xt = 0.0_dp
            if (j == ntheta/2 .and. k == nphi/2 .and. trim(mode) == 'null') then
                write(unit, '(A,11ES25.16E3)') ',', x(2:3), xt, xp, b
            else if (j == ntheta/2 .and. k == nphi/2 .and. trim(mode) == 'extra') then
                write(unit, '(13ES25.16E3)') x, xt, xp, b, 1.0_dp
            else
                write(unit, '(12ES25.16E3)') x, xt, xp, b
            end if
        end do
    end do
    close(unit)

    dl = 2.0_dp * pi * 0.9_dp / real(nloop, dp)
    open(newunit=unit, file=trim(dir)//'/poloidal.diagno', status='replace', action='write')
    write(unit, '(I0)') nloop
    do k = 1, nloop
        t = 2.0_dp * pi * (real(k, dp) - 0.5_dp) / real(nloop, dp)
        ! normal (-sin t, 0, cos t): polar angle acos(cos t), azimuth 180 or 0 deg
        write(unit, '(6ES25.16E3)') r0 + 0.9_dp * cos(t), 0.0_dp, 0.9_dp * sin(t), &
            merge(180.0_dp, 0.0_dp, sin(t) > 0.0_dp), acos(cos(t)) * 180.0_dp / pi, dl
    end do
    close(unit)
    call toroidal('toroidal_in.diagno', r0)
    call toroidal('toroidal_out.diagno', 4.0_dp)
    open(newunit=unit, file=trim(dir)//'/expected.txt', status='replace', action='write')
    write(unit, '(3ES25.16E3)') 2.0_dp * pi * a * bp, -2.0_dp * pi * r0 * b0, 0.0_dp
    close(unit)
contains
    subroutine toroidal(name, radius)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: radius
        integer :: u, i
        real(dp) :: p

        open(newunit=u, file=trim(dir)//'/'//name, status='replace', action='write')
        write(u, '(I0)') nloop
        do i = 1, nloop
            p = 2.0_dp * pi * (real(i, dp) - 0.5_dp) / real(nloop, dp)
            write(u, '(6ES25.16E3)') radius * cos(p), radius * sin(p), 0.0_dp, &
                p * 180.0_dp / pi + 90.0_dp, 90.0_dp, 2.0_dp * pi * radius / real(nloop, dp)
        end do
        close(u)
    end subroutine toroidal
end program plasma_surface_fixture
