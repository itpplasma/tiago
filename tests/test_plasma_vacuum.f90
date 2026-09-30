program test_plasma_vacuum
    !! Plasma sheet-current model on a vacuum VMEC equilibrium (no plasma currents).
    !! Usage: test_plasma_vacuum <wout.nc>
    !! The boundary sheet current mu0 K = n x B reproduces the field of the plasma
    !! currents outside the boundary and minus the coil field inside. Without plasma
    !! currents therefore
    !!  1. its field vanishes outside the plasma, and
    !!  2. a poloidal loop around the plasma sees -phiedge, so the loop signal with
    !!     DIAGNO's idia = 1 correction (+ phiedge), the diamagnetic flux, vanishes.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_get_var, nf90_nowrite
    use tiago_plasma_support, only: plasma_support_t
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: nphi = 64, ntheta = 64  !! grid per field period
    integer, parameter :: nside = 64              !! A samples per side of the poloidal loop
    type(plasma_support_t) :: plasma
    character(len=512) :: wout
    real(dp) :: b0, flux, d, volume, lo(2), hi(2), corner(2, 5), p(2), dl(2)
    real(dp), allocatable :: pts(:, :), normal(:, :), field(:, :), loop(:, :), dls(:, :), a(:, :)
    real(dp), allocatable :: r(:), z(:)
    logical, allocatable :: plane(:)
    integer :: i, k, n, failures

    failures = 0
    call get_command_argument(1, wout)
    call read_b0(trim(wout), b0)
    call plasma%init_from_vmec(trim(wout), nphi, ntheta)
    d = 3.0_dp * plasma%spacing

    ! 1. Field of the sheet 3 grid spacings outside the boundary, along the normal.
    n = size(plasma%xs, 2)
    allocate(normal(3, n))
    do k = 1, n
        normal(:, k) = cross(plasma%xs(:, neighbour(k, 1, 0)) - plasma%xs(:, neighbour(k, -1, 0)), &
            plasma%xs(:, neighbour(k, 0, 1)) - plasma%xs(:, neighbour(k, 0, -1)))
        normal(:, k) = normal(:, k) / norm2(normal(:, k))
    end do
    volume = sum(plasma%xs * normal)
    if (volume < 0.0_dp) normal = -normal
    pts = plasma%xs + d * normal
    allocate(field, mold=pts)
    call plasma%sample_bfield(pts, field)
    call check_small('max |B_sheet| outside / |B0|', maxval(norm2(field, dim=1)), abs(b0), 2.0e-5_dp)

    ! 2. Rectangle around the cross-section at phi = 0, d outside its bounding box.
    plane = plasma%phi == 0.0_dp
    r = pack(hypot(plasma%xs(1, :), plasma%xs(2, :)), plane)
    z = pack(plasma%xs(3, :), plane)
    lo = [minval(r), minval(z)] - d
    hi = [maxval(r), maxval(z)] + d
    corner = reshape([lo(1), lo(2), hi(1), lo(2), hi(1), hi(2), lo(1), hi(2), lo(1), lo(2)], [2, 5])
    allocate(loop(3, 4 * nside), dls(3, 4 * nside), a(3, 4 * nside))
    k = 0
    do i = 1, 4
        dl = (corner(:, i + 1) - corner(:, i)) / real(nside, dp)
        do n = 1, nside
            k = k + 1
            p = corner(:, i) + (real(n, dp) - 0.5_dp) * dl
            loop(:, k) = [p(1), 0.0_dp, p(2)]
            dls(:, k) = [dl(1), 0.0_dp, dl(2)]
        end do
    end do
    call plasma%sample_vector_potential(loop, a)
    flux = sum(a * dls)
    call check_small('|loop flux + phiedge| / phiedge (diamagnetic flux)', &
        abs(flux + plasma%diamagnetic_flux), abs(plasma%diamagnetic_flux), 1.0e-4_dp)

    call plasma%finalize()
    if (failures > 0) error stop 1
    print '(A)', 'test_plasma_vacuum passed'

contains

    integer function neighbour(k, dt, dp_) result(j)
        !! Index of the grid point dt steps in theta and dp_ in phi from point k.
        integer, intent(in) :: k, dt, dp_
        integer :: it, ip, nphi_total

        nphi_total = size(plasma%xs, 2) / ntheta
        it = modulo(k - 1, ntheta)
        ip = (k - 1) / ntheta
        j = modulo(ip + dp_, nphi_total) * ntheta + modulo(it + dt, ntheta) + 1
    end function neighbour

    pure function cross(u, v) result(w)
        real(dp), intent(in) :: u(3), v(3)
        real(dp) :: w(3)
        w = [u(2) * v(3) - u(3) * v(2), u(3) * v(1) - u(1) * v(3), u(1) * v(2) - u(2) * v(1)]
    end function cross

    subroutine read_b0(path, b0)
        character(len=*), intent(in) :: path
        real(dp), intent(out) :: b0
        integer :: ncid, varid, status

        status = nf90_open(path, nf90_nowrite, ncid)
        if (status /= 0) error stop 'cannot open wout file'
        status = nf90_inq_varid(ncid, 'b0', varid)
        status = nf90_get_var(ncid, varid, b0)
        status = nf90_close(ncid)
    end subroutine read_b0

    subroutine check_small(name, value, scale, rtol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: value, scale, rtol

        if (value > rtol * scale) then
            write(error_unit, '(A,ES12.4,A,ES9.2)') 'FAIL '//name//': ', value / scale, ' > ', rtol
            failures = failures + 1
        else
            print '(A,ES12.4)', 'ok   '//name//': ', value / scale
        end if
    end subroutine check_small
end program test_plasma_vacuum
