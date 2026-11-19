module tiago_geqdsk_plasma
    !! Plasma current and flux labels of an axisymmetric equilibrium from a
    !! G-EQDSK file, in the convention of docs/axisymmetric-free-boundary.md:
    !! (R, phi, Z) right-handed, currents positive along +phi, psi = R A_phi
    !! per radian. The file's current is taken along +phi. Its psi is converted
    !! as psi = sigma psi_file / (2 pi)^e, with e in {0, 1} chosen so that the
    !! integral of j_phi matches the file's current and sigma from the sign.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use geqdsk_tools, only: geqdsk_t, geqdsk_read
    use batch_interpolate, only: BatchSplineData1D, BatchSplineData2D, &
        construct_batch_splines_1d, construct_batch_splines_2d, &
        evaluate_batch_splines_1d, evaluate_batch_splines_2d_der
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: mu0 = 4.0e-7_dp * pi
    ! libneo stores G-EQDSK data in CGS.
    real(dp), parameter :: cm = 1.0e-2_dp, maxwell = 1.0e-8_dp, barye = 0.1_dp, &
        gauss = 1.0e-4_dp
    integer, parameter :: nflux = 2001

    type, public :: geqdsk_plasma_t
        integer :: nw = 0, nh = 0
        real(dp), allocatable :: r(:), z(:), psi(:, :) !! grid [m], psi [Wb/rad]
        real(dp), allocatable :: rb(:), zb(:) !! last closed surface [m]
        real(dp) :: psi_axis = 0, psi_b = 0 !! [Wb/rad]
        real(dp) :: r_axis = 0, z_axis = 0, current = 0 !! [m], file current [A]
        real(dp) :: sigma = 1 !! sign of the psi conversion
        integer :: exponent = 0 !! e of the 2 pi factor
        real(dp), allocatable :: phi(:) !! toroidal flux on psi_N [Wb]
        type(BatchSplineData2D) :: psi_spline
        type(BatchSplineData1D) :: profile_spline !! p', F F', |q| on psi_N
    contains
        procedure :: label => plasma_label
        procedure :: rho_tor => plasma_rho_tor
        procedure :: inside => plasma_inside
        procedure :: filaments => plasma_filaments
    end type geqdsk_plasma_t

    public :: read_geqdsk_plasma

contains

    subroutine read_geqdsk_plasma(path, eq, ierr, message)
        character(*), intent(in) :: path
        type(geqdsk_plasma_t), intent(out) :: eq
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        type(geqdsk_t) :: g
        real(dp), allocatable :: prof(:, :), rf(:), zf(:), wf(:), x(:)
        real(dp) :: raw, ratio(0:1), dx
        integer :: i
        message = ''
        call geqdsk_read(g, path, ierr)
        if (ierr /= 0) then
            message = 'cannot read G-EQDSK file '//trim(path)
            return
        end if
        eq%nw = g%nw
        eq%nh = g%nh
        eq%r = g%r_eqd * cm
        eq%z = g%z_eqd * cm
        eq%psi = g%psirz * maxwell
        eq%psi_axis = g%simag * maxwell
        eq%psi_b = g%sibry * maxwell
        eq%r_axis = g%rmaxis * cm
        eq%z_axis = g%zmaxis * cm
        eq%current = g%current / (0.1_dp * 2.99792458e10_dp)
        eq%rb = g%rbbbs * cm
        eq%zb = g%zbbbs * cm
        allocate(prof(g%nw, 3))
        prof(:, 1) = g%pprime * barye / maxwell
        prof(:, 2) = g%ffprim * gauss
        prof(:, 3) = abs(g%qpsi)
        call construct_batch_splines_1d(0.0_dp, 1.0_dp, prof, 3, .false., eq%profile_spline)
        call construct_batch_splines_2d([eq%r(1), eq%z(1)], [eq%r(eq%nw), eq%z(eq%nh)], &
            reshape(eq%psi, [eq%nw, eq%nh, 1]), [3, 3], [.false., .false.], eq%psi_spline)
        ! File units first, then the convention from the plasma current.
        call eq%filaments(rf, zf, wf)
        raw = sum(wf)
        if (raw == 0 .or. eq%current == 0) then
            ierr = 2
            message = 'no plasma current in '//trim(path)
            return
        end if
        ratio(0) = abs(raw / eq%current)
        ratio(1) = abs(2 * pi * raw / eq%current)
        eq%exponent = minloc(abs(log(ratio)), 1) - 1
        if (abs(log(ratio(eq%exponent))) > 0.1_dp) then
            ierr = 3
            message = 'integral of j_phi does not match the current of '//trim(path)
            return
        end if
        eq%sigma = sign(1.0_dp, raw * eq%current)
        eq%psi = eq%sigma * eq%psi / (2 * pi)**eq%exponent
        eq%psi_axis = eq%sigma * eq%psi_axis / (2 * pi)**eq%exponent
        eq%psi_b = eq%sigma * eq%psi_b / (2 * pi)**eq%exponent
        prof(:, 1:2) = eq%sigma * (2 * pi)**eq%exponent * prof(:, 1:2)
        call construct_batch_splines_1d(0.0_dp, 1.0_dp, prof, 3, .false., eq%profile_spline)
        call construct_batch_splines_2d([eq%r(1), eq%z(1)], [eq%r(eq%nw), eq%z(eq%nh)], &
            reshape(eq%psi, [eq%nw, eq%nh, 1]), [3, 3], [.false., .false.], eq%psi_spline)
        ! Phi(psi_N) = 2 pi |psi_b - psi_axis| int_0^psi_N |q| dpsi_N, trapezoidal.
        allocate(eq%phi(nflux), x(nflux))
        dx = 1.0_dp / (nflux - 1)
        do i = 1, nflux
            x(i) = q_at(eq, (i - 1) * dx)
        end do
        eq%phi(1) = 0
        do i = 2, nflux
            eq%phi(i) = eq%phi(i - 1) + (x(i - 1) + x(i)) * dx / 2
        end do
        eq%phi = 2 * pi * abs(eq%psi_b - eq%psi_axis) * eq%phi
    end subroutine read_geqdsk_plasma

    real(dp) function q_at(eq, x)
        type(geqdsk_plasma_t), intent(in) :: eq
        real(dp), intent(in) :: x
        real(dp) :: y(3)
        call evaluate_batch_splines_1d(eq%profile_spline, x, y)
        q_at = y(3)
    end function q_at

    subroutine plasma_label(self, r, z, psin, dpsin_dr, dpsin_dz)
        !! psi_N = (psi - psi_axis) / (psi_b - psi_axis) and its gradient [1/m].
        class(geqdsk_plasma_t), intent(in) :: self
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: psin, dpsin_dr, dpsin_dz
        real(dp) :: y(1), dy(2, 1), scale
        call evaluate_batch_splines_2d_der(self%psi_spline, [r, z], y, dy)
        scale = 1 / (self%psi_b - self%psi_axis)
        psin = (y(1) - self%psi_axis) * scale
        dpsin_dr = dy(1, 1) * scale
        dpsin_dz = dy(2, 1) * scale
    end subroutine plasma_label

    real(dp) function plasma_rho_tor(self, psin) result(rho)
        !! sqrt(Phi(psi_N) / Phi(1)), linear in the tabulated Phi.
        class(geqdsk_plasma_t), intent(in) :: self
        real(dp), intent(in) :: psin
        real(dp) :: t, f
        integer :: i
        t = min(max(psin, 0.0_dp), 1.0_dp) * (nflux - 1)
        i = min(int(t) + 1, nflux - 1)
        f = self%phi(i) + (t - (i - 1)) * (self%phi(i + 1) - self%phi(i))
        rho = sqrt(f / self%phi(nflux))
    end function plasma_rho_tor

    logical function plasma_inside(self, r, z) result(in)
        !! Inside the last closed surface polygon and below psi_N = 1.
        class(geqdsk_plasma_t), intent(in) :: self
        real(dp), intent(in) :: r, z
        real(dp) :: psin, a, b
        integer :: i, j
        in = .false.
        j = size(self%rb)
        do i = 1, size(self%rb)
            if ((self%zb(i) > z) .neqv. (self%zb(j) > z)) then
                if (r < self%rb(i) + (z - self%zb(i)) * (self%rb(j) - self%rb(i)) &
                    / (self%zb(j) - self%zb(i))) in = .not. in
            end if
            j = i
        end do
        if (.not. in) return
        call self%label(r, z, psin, a, b)
        in = psin >= 0 .and. psin < 1
    end function plasma_inside

    subroutine plasma_filaments(self, rf, zf, wf)
        !! Grid nodes inside the plasma with the current j_phi dR dZ [A] of
        !! j_phi = R p'(psi) + F F'(psi) / (mu0 R).
        class(geqdsk_plasma_t), intent(in) :: self
        real(dp), allocatable, intent(out) :: rf(:), zf(:), wf(:)
        real(dp) :: area, psin, a, b, y(3), rbuf(self%nw * self%nh), &
            zbuf(self%nw * self%nh), wbuf(self%nw * self%nh)
        integer :: i, j, n
        area = (self%r(2) - self%r(1)) * (self%z(2) - self%z(1))
        n = 0
        do j = 1, self%nh
            do i = 1, self%nw
                if (.not. self%inside(self%r(i), self%z(j))) cycle
                call self%label(self%r(i), self%z(j), psin, a, b)
                call evaluate_batch_splines_1d(self%profile_spline, psin, y)
                n = n + 1
                rbuf(n) = self%r(i)
                zbuf(n) = self%z(j)
                wbuf(n) = (self%r(i) * y(1) + y(2) / (mu0 * self%r(i))) * area
            end do
        end do
        rf = rbuf(:n)
        zf = zbuf(:n)
        wf = wbuf(:n)
    end subroutine plasma_filaments

end module tiago_geqdsk_plasma
