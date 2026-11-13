module tiago_spectre_source
    !! SPECTRE adapter: libneo owns file reading, geometry, inversion and field
    !! evaluation. Tiago's existing virtual casing owns exterior plasma B/A.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    use tiago_plasma_field_source, only: plasma_field_source_t, field_hint_t
    use tiago_plasma_support, only: plasma_support_t
    use field_spectre, only: spectre_field_t, create_spectre_field
    use field_base, only: magnetic_field_t
    use libneo_coordinates, only: SPECTRE_SUCCESS, SPECTRE_OUTSIDE_DOMAIN, &
        spectre_coordinate_system_t
    implicit none
    private
    public :: spectre_source_t

    real(dp), parameter :: twopi = 2.0_dp*acos(-1.0_dp)

    type, extends(plasma_field_source_t) :: spectre_source_t
        private
        type(spectre_field_t) :: field
        type(plasma_support_t) :: casing
        integer :: plasma_volumes = 0
    contains
        procedure :: init => spectre_source_init
        procedure :: sample_bfield => spectre_source_bfield
        procedure :: sample_vector_potential => spectre_source_potential
        procedure :: sample_interior_bfield => spectre_source_interior
        procedure :: sample_interior_bfield_status => spectre_source_interior_status
        procedure :: sample_interior_bfield_with_hint => spectre_source_with_hint
    end type spectre_source_t

contains

    subroutine spectre_source_init(self, filename, nphi, ntheta, linked_flux)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        class(spectre_source_t), intent(inout) :: self
        character(len=*), intent(in) :: filename
        integer(i32), intent(in) :: nphi, ntheta
        real(dp), intent(in), optional :: linked_flux
        real(dp), allocatable :: xs(:, :), xt(:, :), xp(:, :), b(:, :)
        real(dp) :: u(3), e(3, 3), a(3), h(3), bmod, sqgb(3)
        real(dp) :: dtheta, dphi, normal_fraction
        real(dp) :: g(3, 3), ginv(3, 3), sqrtg
        integer :: ierr, ntor, j, k, p

        self%enabled = .false.
        self%linked_flux_known = present(linked_flux)
        self%diamagnetic_flux = 0.0_dp
        if (present(linked_flux)) then
            if (.not. ieee_is_finite(linked_flux)) error stop 'nonfinite linked flux'
            self%diamagnetic_flux = linked_flux
        end if
        call create_spectre_field(self%field, filename, ierr)
        if (ierr /= 0) error stop 'cannot initialize libneo SPECTRE field'
        self%plasma_volumes = self%field%data%Nvol
        if (nphi < 4 .or. ntheta < 4) error stop 'plasma grid needs at least 4x4 points'
        ntor = nphi*self%field%data%Nfp
        dtheta = twopi/real(ntheta, dp)
        dphi = twopi/real(ntor, dp)
        allocate (xs(3, ntor*ntheta), xt(3, ntor*ntheta))
        allocate (xp(3, ntor*ntheta), b(3, ntor*ntheta))
        ! The casing is the last PLASMA interface, not the computational wall.
        ! Use the vacuum trace when present, including any physical boundary
        ! sheet current. Fixed-boundary files only supply the plasma trace and
        ! therefore require the usual no-boundary-sheet assumption.
        u(1) = real(self%plasma_volumes, dp) &
               - 64.0_dp*epsilon(1.0_dp)*real(self%plasma_volumes, dp)
        if (self%field%data%Mvol > self%plasma_volumes) &
            u(1) = real(self%plasma_volumes, dp) &
                   + 64.0_dp*epsilon(1.0_dp)*real(self%plasma_volumes, dp)
        normal_fraction = 0.0_dp
        p = 0
        do k = 0, ntor - 1
            u(3) = real(k, dp)*dphi
            do j = 0, ntheta - 1
                u(2) = real(j, dp)*dtheta
                p = p + 1
                call self%field%coords%evaluate_cart(u, xs(:, p))
                call self%field%coords%covariant_basis(u, e)
                xt(:, p) = e(:, 2)
                xp(:, p) = e(:, 3)
                call self%field%evaluate(u, a, h, bmod, sqgb)
                call self%field%coords%metric_tensor(u, g, ginv, sqrtg)
                call self%field%coords%cov_to_cart(u, h*bmod, b(:, p))
                ! sqrt(g) B^rho is the normal flux density. Require a flux
                ! surface: the current-only virtual casing omits magnetic charge.
                normal_fraction = max(normal_fraction, &
                                      abs(sqgb(1))/ &
                                      max(sqrtg*sqrt(ginv(1, 1))*bmod, tiny(bmod)))
            end do
        end do
        ! Compute B.n directly from the chart metric instead of a physical
        ! extrapolation outside SPECTRE's represented domain.
        if (normal_fraction > 1.0e-6_dp) &
            error stop 'SPECTRE outer plasma boundary is not a flux surface'
        call self%casing%init_from_surface(xs, xt, xp, b, dtheta, dphi)
        ! libneo currently exposes normalized interface fluxes, but not the
        ! signed global toroidal flux. A_theta has a separate gauge in each
        ! volume: its outer circulation alone is NOT the global linked flux.
        ! Until libneo supplies it, require explicit metadata for idia=1 loops.
        self%enabled = .true.
    end subroutine spectre_source_init

    subroutine spectre_source_bfield(self, points, bfield)
        class(spectre_source_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: bfield(:, :)

        call require_exterior(self, points)
        call self%casing%sample_bfield(points, bfield)
    end subroutine spectre_source_bfield

    subroutine spectre_source_potential(self, points, avec)
        class(spectre_source_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: avec(:, :)

        call require_exterior(self, points)
        call self%casing%sample_vector_potential(points, avec)
    end subroutine spectre_source_potential

    subroutine spectre_source_interior(self, points, values)
        class(spectre_source_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer :: status(size(points, 2))

        if (.not. self%enabled) error stop 'SPECTRE source not initialized'
        call self%sample_interior_bfield_status(points, values, status)
        if (any(status /= SPECTRE_SUCCESS)) error stop 'cannot locate interior SPECTRE point'
    end subroutine spectre_source_interior

    subroutine spectre_source_interior_status(self, points, values, status)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        class(spectre_source_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: status(:)
        real(dp) :: u(3)
        integer :: j, ierr

        values = 0.0_dp
        status = 1
        if (.not. self%enabled) return
        !$omp parallel do private(j, u, ierr) schedule(static)
        do j = 1, size(points, 2)
            if (.not. all(ieee_is_finite(points(:, j)))) cycle
            call locate(self, points(:, j), u, ierr)
            status(j) = ierr
            if (ierr /= SPECTRE_SUCCESS) cycle
            if (u(1) >= real(self%plasma_volumes, dp)) then
                status(j) = SPECTRE_OUTSIDE_DOMAIN
                cycle
            end if
            call chart_bfield(self%field, u, values(:, j))
            if (.not. all(ieee_is_finite(values(:, j)))) then
                values(:, j) = 0.0_dp
                status(j) = 1
            end if
        end do
        !$omp end parallel do
    end subroutine spectre_source_interior_status

    subroutine spectre_source_with_hint(self, point, value, status, hint)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        class(spectre_source_t), intent(in) :: self
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: value(3)
        integer, intent(out) :: status
        type(field_hint_t), intent(inout) :: hint
        real(dp) :: u(3), cyl(3)

        value = 0.0_dp
        status = 1
        if (.not. self%enabled .or. .not. all(ieee_is_finite(point))) then
            hint = field_hint_t()
            return
        end if
        if (hint%volume >= 1 .and. hint%volume <= self%field%data%Mvol &
            .and. all(ieee_is_finite(hint%u))) then
            cyl = [sqrt(point(1)**2 + point(2)**2), &
                   atan2(point(2), point(1)), point(3)]
            u = hint%u
            select type (cs => self%field%coords)
            type is (spectre_coordinate_system_t)
                call cs%from_cyl_warm(cyl, u, hint%volume, status)
            end select
            if (status == SPECTRE_SUCCESS) then
                ! Warm inversion snaps near interfaces; preserve the cold sheet side.
                if (abs(u(1) - real(nint(u(1)), dp)) <= 0.5e-9_dp) status = 1
            end if
            if (status == SPECTRE_SUCCESS) &
                call interior_at_labels(self, u, value, status)
        end if
        ! A previous volume is only a seed, never a domain restriction.
        if (status /= SPECTRE_SUCCESS) then
            call locate(self, point, u, status)
            if (status == SPECTRE_SUCCESS) &
                call interior_at_labels(self, u, value, status)
        end if
        if (status == SPECTRE_SUCCESS) then
            hint%u = u
            hint%volume = min(int(u(1)) + 1, self%field%data%Mvol)
        else
            hint = field_hint_t()
        end if
    end subroutine spectre_source_with_hint

    subroutine interior_at_labels(self, u, value, status)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        class(spectre_source_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: value(3)
        integer, intent(out) :: status

        value = 0.0_dp
        status = SPECTRE_OUTSIDE_DOMAIN
        if (u(1) >= real(self%plasma_volumes, dp)) return
        call chart_bfield(self%field, u, value)
        status = SPECTRE_SUCCESS
        if (.not. all(ieee_is_finite(value))) then
            value = 0.0_dp
            status = 1
        end if
    end subroutine interior_at_labels

    subroutine chart_bfield(field, u, b)
        class(magnetic_field_t), intent(in) :: field
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: b(3)
        real(dp) :: a(3), h(3), bmod

        call field%evaluate(u, a, h, bmod)
        call field%coords%cov_to_cart(u, h*bmod, b)
    end subroutine chart_bfield

    subroutine locate(self, point, u, ierr)
        class(spectre_source_t), intent(in) :: self
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: ierr
        real(dp) :: cyl(3)

        cyl = [sqrt(point(1)**2 + point(2)**2), &
               atan2(point(2), point(1)), point(3)]
        call self%field%coords%from_cyl(cyl, u, ierr)
    end subroutine locate

    subroutine require_exterior(self, points)
        !! Accept certified chart exterior or a located vacuum-volume point.
        !! Numerical failure and unresolved membership do not establish exterior.
        class(spectre_source_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp) :: u(3)
        integer :: j, ierr

        if (.not. self%enabled) error stop 'SPECTRE source not initialized'
        do j = 1, size(points, 2)
            call locate(self, points(:, j), u, ierr)
            if (ierr == SPECTRE_OUTSIDE_DOMAIN) cycle
            if (ierr /= SPECTRE_SUCCESS) &
                error stop 'cannot determine SPECTRE diagnostic domain'
            if (u(1) <= real(self%plasma_volumes, dp) + 1.0e-9_dp) &
                error stop 'virtual casing requires exterior diagnostic points'
        end do
    end subroutine require_exterior
end module tiago_spectre_source
