module tiago_cxrs
    !! Charge-exchange recombination spectroscopy (CXRS) on straight sight lines
    !! crossing a neutral heating beam. Channel c returns
    !!
    !!   b_c  = s_c int_c em dl                 (brightness, scale s_c),
    !!   Ti_c = int_c em Ti dl / int_c em dl     (emission-weighted Ti),
    !!
    !! signals (b_1 .. b_n, Ti_1 .. Ti_n). The emission at a sight point x is
    !!
    !!   em = g(x) exp(-a(x)) ne (Zeff - 1) / (Z (Z - 1))   (FortSym cxrs_emission),
    !!   g  = exp(-r_perp^2 / w^2),  a = sigma int_0^xi ne ds along the beam axis,
    !!
    !! with r_perp the distance of x from the beam axis, xi its projection on the
    !! axis from the beam origin, w the 1/e beam radius and sigma a constant
    !! effective stopping cross-section. The attenuation is that of the beam axis
    !! up to the projection of x. The charge-exchange rate coefficient, beam
    !! energy fractions and etendue are constant per channel and absorbed in s_c;
    !! halo neutrals and the dependence of sigma on ne, Te and beam energy are not
    !! modelled. Ti_c is the second spectral moment of a sum of Gaussian lines;
    !! a single-Gaussian fit of a broadened spectrum approximates it.
    !!
    !! Sample points are the sight-line rule followed, for each sight point, by
    !! its beam-axis rule on [0, xi]. Each point block belongs to one channel, so
    !! the channel loops are free of write conflicts.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t, Q_NE, &
        Q_TI, Q_ZEFF
    use tiago_sight_geometry, only: chord_t, gauss_legendre01, to_rpz
    use tiago_cxrs_emission, only: cxrs_emission
    implicit none
    private

    type, public :: cxrs_beam_t
        real(dp) :: origin(3) = 0.0_dp !! Cartesian [m], outside the plasma
        real(dp) :: direction(3) = [1.0_dp, 0.0_dp, 0.0_dp] !! normalized on use
        real(dp) :: width = 0.1_dp !! 1/e radius w [m]
        real(dp) :: stopping = 0.0_dp !! effective stopping cross-section [m^2]
        integer :: npanel = 40 !! panels of the beam-axis rule
        integer :: ngauss = 4 !! Gauss points per panel
    end type cxrs_beam_t

    type, public, extends(profile_diagnostic_t) :: cxrs_t
        type(chord_t), allocatable :: chords(:)
        type(cxrs_beam_t), allocatable :: beams(:) !! beam seen by each channel
        real(dp) :: z_imp = 6.0_dp !! impurity charge (carbon)
        real(dp), allocatable :: scale(:) !! brightness scale per channel
        real(dp), allocatable :: x0(:, :) !! sample points (R, phi, Z)
        real(dp), allocatable :: w(:), g(:) !! sight weights [m], beam profile
        real(dp), allocatable :: u(:) !! sigma times beam-axis weights [m^3]
        integer, allocatable :: first(:) !! sight points of channel c
        integer, allocatable :: beam_first(:) !! beam points of sight point k
    contains
        procedure :: init => cxrs_init
        procedure :: n_signals => cxrs_n_signals
        procedure :: n_nuisance => cxrs_n_nuisance
        procedure :: label => cxrs_label
        procedure :: sample_points => cxrs_points
        procedure :: predict => cxrs_predict
        procedure :: jvp => cxrs_jvp
        procedure :: vjp => cxrs_vjp
    end type cxrs_t

contains

    subroutine cxrs_init(self, chords, beams, status)
        !! status 1: chords and beams differ in size; 2: a beam has a
        !! non-positive width, a negative stopping cross-section, no direction
        !! or an empty rule.
        class(cxrs_t), intent(inout) :: self
        type(chord_t), intent(in) :: chords(:)
        type(cxrs_beam_t), intent(in) :: beams(:)
        integer, intent(out) :: status
        real(dp), allocatable :: xs(:), ws(:), xb(:), wb(:), sight(:, :), axis(:, :)
        real(dp) :: d(3), db(3), x(3), rel(3), xi, h
        integer :: c, i, j, k, m, ns, nb, nsight, nbeam

        status = 1
        if (size(chords) /= size(beams)) return
        status = 2
        do c = 1, size(beams)
            if (.not. beams(c)%width > 0.0_dp) return
            if (.not. beams(c)%stopping >= 0.0_dp) return
            if (.not. norm2(beams(c)%direction) > 0.0_dp) return
            if (beams(c)%npanel < 1) return
            if (beams(c)%ngauss < 1) return
        end do
        status = 0
        self%chords = chords
        self%beams = beams
        if (allocated(self%scale)) then
            if (size(self%scale) /= size(chords)) deallocate(self%scale)
        end if
        if (.not. allocated(self%scale)) then
            allocate(self%scale(size(chords)))
            self%scale = 1.0_dp
        end if

        nsight = sum(chords%npanel * chords%ngauss)
        nbeam = 0
        do c = 1, size(chords)
            nbeam = nbeam + chords(c)%npanel * chords(c)%ngauss &
                * beams(c)%npanel * beams(c)%ngauss
        end do
        if (allocated(self%first)) deallocate(self%first)
        if (allocated(self%beam_first)) deallocate(self%beam_first)
        if (allocated(self%x0)) deallocate(self%x0, self%w, self%g, self%u)
        allocate(self%first(size(chords) + 1), self%beam_first(nsight + 1), &
            sight(3, nsight), axis(3, nbeam), self%w(nsight), self%g(nsight), &
            self%u(nbeam))

        ns = 0
        nb = 0
        self%first(1) = 1
        do c = 1, size(chords)
            call gauss_legendre01(chords(c)%ngauss, xs, ws)
            call gauss_legendre01(beams(c)%ngauss, xb, wb)
            d = chords(c)%direction / norm2(chords(c)%direction)
            db = beams(c)%direction / norm2(beams(c)%direction)
            h = chords(c)%length / chords(c)%npanel
            do i = 1, chords(c)%npanel
                do j = 1, chords(c)%ngauss
                    ns = ns + 1
                    x = chords(c)%origin + (i - 1 + xs(j)) * h * d
                    sight(:, ns) = to_rpz(x)
                    self%w(ns) = h * ws(j)
                    rel = x - beams(c)%origin
                    xi = dot_product(rel, db)
                    self%g(ns) = exp(-sum((rel - xi * db)**2) / beams(c)%width**2)
                    ! Points behind the beam origin see the unattenuated beam.
                    xi = max(xi, 0.0_dp) / beams(c)%npanel
                    self%beam_first(ns) = nb + 1
                    do k = 1, beams(c)%npanel
                        do m = 1, beams(c)%ngauss
                            nb = nb + 1
                            axis(:, nb) = to_rpz(beams(c)%origin + (k - 1 + xb(m)) * xi * db)
                            self%u(nb) = beams(c)%stopping * xi * wb(m)
                        end do
                    end do
                end do
            end do
            self%first(c + 1) = ns + 1
        end do
        self%beam_first(nsight + 1) = nb + 1
        self%beam_first = self%beam_first + nsight
        allocate(self%x0(3, nsight + nbeam))
        self%x0(:, :nsight) = sight
        self%x0(:, nsight + 1:) = axis
    end subroutine cxrs_init

    integer function cxrs_n_signals(self)
        class(cxrs_t), intent(in) :: self
        cxrs_n_signals = 2 * size(self%chords)
    end function cxrs_n_signals

    integer function cxrs_n_nuisance(self)
        class(cxrs_t), intent(in) :: self
        cxrs_n_nuisance = 0
    end function cxrs_n_nuisance

    function cxrs_label(self, i) result(text)
        class(cxrs_t), intent(in) :: self
        integer, intent(in) :: i
        character(len=:), allocatable :: text
        integer :: n
        n = size(self%chords)
        if (i <= n) then
            text = 'brightness CXRS '//self%chords(i)%label
        else
            text = 'Ti_CXRS '//self%chords(i - n)%label
        end if
    end function cxrs_label

    subroutine cxrs_points(self, eta, x, dx_deta)
        class(cxrs_t), intent(in) :: self
        real(dp), intent(in) :: eta(:)
        real(dp), allocatable, intent(out) :: x(:, :), dx_deta(:, :, :)
        x = self%x0
        allocate(dx_deta(3, size(x, 2), 0))
    end subroutine cxrs_points

    pure subroutine emission(self, sample, k, em, de)
        !! em and d em / d (ne, Zeff, a) at sight point k; zero outside.
        class(cxrs_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        integer, intent(in) :: k
        real(dp), intent(out) :: em, de(3)
        real(dp) :: a
        integer :: j
        em = 0.0_dp
        de = 0.0_dp
        if (.not. sample%inside(k)) return
        a = 0.0_dp
        do j = self%beam_first(k), self%beam_first(k + 1) - 1
            if (sample%inside(j)) a = a + self%u(j - self%beam_first(1) + 1) &
                * sample%q(Q_NE, j)
        end do
        call cxrs_emission(sample%q(Q_NE, k), sample%q(Q_ZEFF, k), a, self%g(k), &
            self%z_imp, em, de(1), de(2), de(3))
    end subroutine emission

    pure real(dp) function attenuation_dot(self, sample, k, q_dot)
        class(cxrs_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        integer, intent(in) :: k
        real(dp), intent(in) :: q_dot(:, :)
        integer :: j
        attenuation_dot = 0.0_dp
        do j = self%beam_first(k), self%beam_first(k + 1) - 1
            if (sample%inside(j)) attenuation_dot = attenuation_dot &
                + self%u(j - self%beam_first(1) + 1) * q_dot(Q_NE, j)
        end do
    end function attenuation_dot

    subroutine cxrs_sums(self, sample, b, et)
        class(cxrs_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(out) :: b(:), et(:)
        real(dp) :: em, de(3)
        integer :: c, k
        !$omp parallel do default(shared) private(c, k, em, de) &
        !$omp if(size(self%x0, 2)>4096)
        do c = 1, size(self%chords)
            b(c) = 0.0_dp
            et(c) = 0.0_dp
            do k = self%first(c), self%first(c + 1) - 1
                call emission(self, sample, k, em, de)
                b(c) = b(c) + self%w(k) * em
                et(c) = et(c) + self%w(k) * em * sample%q(Q_TI, k)
            end do
        end do
        !$omp end parallel do
    end subroutine cxrs_sums

    subroutine cxrs_predict(self, sample, eta, y)
        class(cxrs_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: b(size(self%chords)), et(size(self%chords))
        integer :: n
        n = size(self%chords)
        call cxrs_sums(self, sample, b, et)
        y(:n) = self%scale * b
        y(n + 1:2 * n) = 0.0_dp
        where (b > 0.0_dp) y(n + 1:2 * n) = et / b
    end subroutine cxrs_predict

    subroutine cxrs_jvp(self, sample, eta, q_dot, eta_dot, y_dot)
        class(cxrs_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        real(dp) :: b(size(self%chords)), et(size(self%chords))
        real(dp) :: db(size(self%chords)), det(size(self%chords)), em, de(3), dem
        integer :: c, k, n

        n = size(self%chords)
        call cxrs_sums(self, sample, b, et)
        !$omp parallel do default(shared) private(c, k, em, de, dem) &
        !$omp if(size(self%x0, 2)>4096)
        do c = 1, n
            db(c) = 0.0_dp
            det(c) = 0.0_dp
            do k = self%first(c), self%first(c + 1) - 1
                call emission(self, sample, k, em, de)
                if (.not. sample%inside(k)) cycle
                dem = de(1) * q_dot(Q_NE, k) + de(2) * q_dot(Q_ZEFF, k) &
                    + de(3) * attenuation_dot(self, sample, k, q_dot)
                db(c) = db(c) + self%w(k) * dem
                det(c) = det(c) + self%w(k) * (dem * sample%q(Q_TI, k) &
                    + em * q_dot(Q_TI, k))
            end do
        end do
        !$omp end parallel do
        y_dot(:n) = self%scale * db
        y_dot(n + 1:2 * n) = 0.0_dp
        where (b > 0.0_dp) y_dot(n + 1:2 * n) = (det - et / b * db) / b
    end subroutine cxrs_jvp

    subroutine cxrs_vjp(self, sample, eta, y_bar, q_bar, eta_bar)
        class(cxrs_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        real(dp) :: b(size(self%chords)), et(size(self%chords))
        real(dp) :: bbar(size(self%chords)), etbar(size(self%chords)), em, de(3), embar
        integer :: c, j, k, n, off

        n = size(self%chords)
        call cxrs_sums(self, sample, b, et)
        bbar = self%scale * y_bar(:n)
        etbar = 0.0_dp
        where (b > 0.0_dp)
            bbar = bbar - y_bar(n + 1:2 * n) * et / b**2
            etbar = y_bar(n + 1:2 * n) / b
        end where
        q_bar = 0.0_dp
        eta_bar = 0.0_dp
        off = self%beam_first(1) - 1
        !$omp parallel do default(shared) private(c, j, k, em, de, embar) &
        !$omp if(size(self%x0, 2)>4096)
        do c = 1, n
            do k = self%first(c), self%first(c + 1) - 1
                call emission(self, sample, k, em, de)
                if (.not. sample%inside(k)) cycle
                embar = self%w(k) * (bbar(c) + etbar(c) * sample%q(Q_TI, k))
                q_bar(Q_NE, k) = q_bar(Q_NE, k) + embar * de(1)
                q_bar(Q_ZEFF, k) = q_bar(Q_ZEFF, k) + embar * de(2)
                q_bar(Q_TI, k) = q_bar(Q_TI, k) + self%w(k) * etbar(c) * em
                do j = self%beam_first(k), self%beam_first(k + 1) - 1
                    if (sample%inside(j)) q_bar(Q_NE, j) = q_bar(Q_NE, j) &
                        + embar * de(3) * self%u(j - off)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine cxrs_vjp
end module tiago_cxrs
