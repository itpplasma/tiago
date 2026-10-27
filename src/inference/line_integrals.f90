module tiago_line_integrals
    !! Line-integrated profile diagnostics on straight chords, with composite
    !! Gauss-Legendre rules (tiago_sight_geometry):
    !!
    !!   interferometer_t   y_i = int ne dl;
    !!   xics_t             b_i = s_i int eps dl and T_i = int eps T dl / int eps dl,
    !!                      signals (b_1 .. b_n, T_1 .. T_n), T = Ti (default) or Te.
    !!
    !! The XICS emissivity is the coronal excitation scaling of one impurity
    !! line, eps = ne n_Z exp(-E_line / Te) / sqrt(Te) with
    !! n_Z = ne (Zeff - 1) / (Z (Z - 1)) (FortSym kernel xics_emissivity); the
    !! brightness scale s_i absorbs atomic data and etendue. A more detailed
    !! emissivity model replaces that kernel with the same partial derivatives.
    !!
    !! Outside the plasma the kinetic quantities vanish, so the chords need no
    !! boundary crossings; the rule does not resolve the kink of the integrand
    !! at the boundary, which limits its accuracy to O(h^2) there.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t, Q_NE, &
        Q_TE, Q_TI, Q_ZEFF
    use tiago_sight_geometry, only: chord_t, chord_rule
    use tiago_xics_emissivity, only: xics_emissivity
    implicit none
    private

    type, public, extends(profile_diagnostic_t) :: chord_set_t
        !! Common chord sampling of interferometer_t and xics_t.
        type(chord_t), allocatable :: chords(:)
        real(dp), allocatable :: x0(:, :), w(:)
        integer, allocatable :: owner(:), first(:)
    contains
        procedure :: init => chords_init
        procedure :: n_nuisance => chords_n_nuisance
        procedure :: sample_points => chords_points
        procedure :: n_signals => interferometer_n_signals
        procedure :: label => interferometer_label
        procedure :: predict => interferometer_predict
        procedure :: jvp => interferometer_jvp
        procedure :: vjp => interferometer_vjp
    end type chord_set_t

    type, public, extends(chord_set_t) :: interferometer_t
    end type interferometer_t

    type, public, extends(chord_set_t) :: xics_t
        integer :: temperature = Q_TI
        real(dp) :: z_imp = 18.0_dp !! impurity charge (argon)
        real(dp) :: e_line = 3.1e3_dp !! line excitation energy [eV]
        real(dp), allocatable :: scale(:) !! brightness scale per chord
    contains
        procedure :: n_signals => xics_n_signals
        procedure :: label => xics_label
        procedure :: predict => xics_predict
        procedure :: jvp => xics_jvp
        procedure :: vjp => xics_vjp
    end type xics_t

contains

    subroutine chords_init(self, chords)
        class(chord_set_t), intent(inout) :: self
        type(chord_t), intent(in) :: chords(:)
        real(dp), allocatable :: p(:, :), w(:)
        integer :: c

        self%chords = chords
        if (allocated(self%first)) deallocate(self%first)
        allocate(self%first(size(chords) + 1))
        self%first(1) = 1
        do c = 1, size(chords)
            self%first(c + 1) = self%first(c) + chords(c)%npanel*chords(c)%ngauss
        end do
        if (allocated(self%x0)) deallocate(self%x0, self%w, self%owner)
        allocate(self%x0(3, 0), self%w(0), self%owner(0))
        do c = 1, size(chords)
            call chord_rule(chords(c), p, w)
            self%x0 = reshape([self%x0, p], [3, size(self%w) + size(w)])
            self%w = [self%w, w]
            self%owner = [self%owner, spread(c, 1, size(w))]
        end do
        select type (self)
            type is (xics_t)
            if (allocated(self%scale)) then
                if (size(self%scale) /= size(chords)) deallocate(self%scale)
            end if
            if (.not. allocated(self%scale)) then
                allocate(self%scale(size(chords)))
                self%scale = 1.0_dp
            end if
        end select
    end subroutine chords_init

    integer function chords_n_nuisance(self)
        class(chord_set_t), intent(in) :: self
        chords_n_nuisance = 0
    end function chords_n_nuisance

    subroutine chords_points(self, eta, x, dx_deta)
        class(chord_set_t), intent(in) :: self
        real(dp), intent(in) :: eta(:)
        real(dp), allocatable, intent(out) :: x(:, :), dx_deta(:, :, :)
        x = self%x0
        allocate(dx_deta(3, size(x, 2), 0))
    end subroutine chords_points

    integer function interferometer_n_signals(self)
        class(chord_set_t), intent(in) :: self
        interferometer_n_signals = size(self%chords)
    end function interferometer_n_signals

    function interferometer_label(self, i) result(text)
        class(chord_set_t), intent(in) :: self
        integer, intent(in) :: i
        character(len=:), allocatable :: text
        text = 'int ne dl '//self%chords(i)%label
    end function interferometer_label

    subroutine interferometer_predict(self, sample, eta, y)
        class(chord_set_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        integer :: c, lo, hi
        !$omp parallel do default(shared) private(c, lo, hi) &
        !$omp if(size(self%w)>4096)
        do c = 1, size(self%chords)
            lo = self%first(c)
            hi = self%first(c + 1) - 1
            y(c) = dot_product(self%w(lo:hi), sample%q(Q_NE, lo:hi))
        end do
        !$omp end parallel do
    end subroutine interferometer_predict

    subroutine interferometer_jvp(self, sample, eta, q_dot, eta_dot, y_dot)
        class(chord_set_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        integer :: c, lo, hi
        !$omp parallel do default(shared) private(c, lo, hi) &
        !$omp if(size(self%w)>4096)
        do c = 1, size(self%chords)
            lo = self%first(c)
            hi = self%first(c + 1) - 1
            y_dot(c) = dot_product(self%w(lo:hi), q_dot(Q_NE, lo:hi))
        end do
        !$omp end parallel do
    end subroutine interferometer_jvp

    subroutine interferometer_vjp(self, sample, eta, y_bar, q_bar, eta_bar)
        class(chord_set_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        integer :: k
        eta_bar = 0.0_dp
        !$omp parallel do default(shared) private(k) if(size(self%w)>4096)
        do k = 1, size(self%w)
            q_bar(:, k) = 0.0_dp
            q_bar(Q_NE, k) = self%w(k) * y_bar(self%owner(k))
        end do
        !$omp end parallel do
    end subroutine interferometer_vjp

    integer function xics_n_signals(self)
        class(xics_t), intent(in) :: self
        xics_n_signals = 2 * size(self%chords)
    end function xics_n_signals

    function xics_label(self, i) result(text)
        class(xics_t), intent(in) :: self
        integer, intent(in) :: i
        character(len=:), allocatable :: text
        integer :: n
        n = size(self%chords)
        if (i <= n) then
            text = 'brightness '//self%chords(i)%label
        else
            text = merge('Ti_XICS ', 'Te_XICS ', self%temperature == Q_TI)// &
                self%chords(i - n)%label
        end if
    end function xics_label

    pure subroutine emissivity(self, sample, k, eps, de)
        !! eps and d eps / d (ne, Te, Zeff) at sample k; zero outside.
        class(xics_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        integer, intent(in) :: k
        real(dp), intent(out) :: eps, de(3)
        eps = 0.0_dp
        de = 0.0_dp
        if (.not. sample%inside(k)) return
        if (sample%q(Q_TE, k) <= 0.0_dp) return
        call xics_emissivity(sample%q(Q_NE, k), sample%q(Q_TE, k), sample%q(Q_ZEFF, k), &
            self%z_imp, self%e_line, eps, de(1), de(2), de(3))
    end subroutine emissivity

    subroutine xics_sums(self, sample, b, et)
        class(xics_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(out) :: b(:), et(:)
        real(dp) :: eps, de(3)
        integer :: k, c
        b = 0.0_dp
        et = 0.0_dp
        do k = 1, size(self%w)
            c = self%owner(k)
            call emissivity(self, sample, k, eps, de)
            b(c) = b(c) + self%w(k) * eps
            et(c) = et(c) + self%w(k) * eps * sample%q(self%temperature, k)
        end do
    end subroutine xics_sums

    subroutine xics_predict(self, sample, eta, y)
        class(xics_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: b(size(self%chords)), et(size(self%chords))
        integer :: n
        n = size(self%chords)
        call xics_sums(self, sample, b, et)
        y(:n) = self%scale * b
        y(n + 1:2 * n) = 0.0_dp
        where (b > 0.0_dp) y(n + 1:2 * n) = et / b
    end subroutine xics_predict

    subroutine xics_jvp(self, sample, eta, q_dot, eta_dot, y_dot)
        class(xics_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        real(dp) :: b(size(self%chords)), et(size(self%chords))
        real(dp) :: db(size(self%chords)), det(size(self%chords)), eps, de(3), deps
        integer :: k, c, n

        n = size(self%chords)
        call xics_sums(self, sample, b, et)
        db = 0.0_dp
        det = 0.0_dp
        do k = 1, size(self%w)
            c = self%owner(k)
            call emissivity(self, sample, k, eps, de)
            deps = de(1) * q_dot(Q_NE, k) + de(2) * q_dot(Q_TE, k) + de(3) * q_dot(Q_ZEFF, k)
            db(c) = db(c) + self%w(k) * deps
            det(c) = det(c) + self%w(k) * (deps * sample%q(self%temperature, k) &
                + eps * q_dot(self%temperature, k))
        end do
        y_dot(:n) = self%scale * db
        y_dot(n + 1:2 * n) = 0.0_dp
        where (b > 0.0_dp) y_dot(n + 1:2 * n) = (det - et / b * db) / b
    end subroutine xics_jvp

    subroutine xics_vjp(self, sample, eta, y_bar, q_bar, eta_bar)
        class(xics_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        real(dp) :: b(size(self%chords)), et(size(self%chords))
        real(dp) :: bbar(size(self%chords)), etbar(size(self%chords)), eps, de(3), epsbar
        integer :: k, c, n

        n = size(self%chords)
        call xics_sums(self, sample, b, et)
        bbar = self%scale * y_bar(:n)
        etbar = 0.0_dp
        where (b > 0.0_dp)
            bbar = bbar - y_bar(n + 1:2 * n) * et / b**2
            etbar = y_bar(n + 1:2 * n) / b
        end where
        q_bar = 0.0_dp
        eta_bar = 0.0_dp
        do k = 1, size(self%w)
            c = self%owner(k)
            call emissivity(self, sample, k, eps, de)
            epsbar = self%w(k) * (bbar(c) + etbar(c) * sample%q(self%temperature, k))
            q_bar(Q_NE, k) = q_bar(Q_NE, k) + epsbar * de(1)
            q_bar(Q_TE, k) = q_bar(Q_TE, k) + epsbar * de(2)
            q_bar(Q_ZEFF, k) = q_bar(Q_ZEFF, k) + epsbar * de(3)
            q_bar(self%temperature, k) = q_bar(self%temperature, k) &
                + self%w(k) * etbar(c) * eps
        end do
    end subroutine xics_vjp
end module tiago_line_integrals
