module tiago_axisym_magnetics
    !! Magnetic sensors of an axisymmetric device from toroidal current filaments,
    !! see docs/axisymmetric-free-boundary.md. (R, phi, Z) is right-handed,
    !! currents are positive along +phi and psi = R A_phi is the flux per radian.
    !! A filament set is (r, z, w): positions [m] and currents or, for coil
    !! groups, signed turns per unit group current.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_special_elliptic, only: ellip_k_from_complement, ellip_e_from_complement
    use fortnum_quadrature, only: gauss_legendre_ab
    use tiago_axisym_green_modulus, only: axisym_green_modulus
    use tiago_axisym_green, only: axisym_green
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: mu0 = 4.0e-7_dp * pi

    !> Pickup coil: centre (r, z), unit axis (nr, nz), length along the unit
    !> direction (tr, tz), effective area N A [m^2], Gauss points along the length.
    type, public :: axisym_pickup_t
        character(len=:), allocatable :: label
        real(dp) :: r = 0, z = 0, nr = 0, nz = 1, tr = 1, tz = 0
        real(dp) :: length = 0, area = 1
        integer :: points = 1
    end type axisym_pickup_t

    type, public :: axisym_flux_loop_t
        character(len=:), allocatable :: label
        real(dp) :: r = 0, z = 0
    end type axisym_flux_loop_t

    !> Rogowski coil along the closed polygon (r(i), z(i)) in the poloidal plane.
    type, public :: axisym_rogowski_t
        character(len=:), allocatable :: label
        real(dp), allocatable :: r(:), z(:)
    end type axisym_rogowski_t

    !> Signals in the order flux loops [Wb], pickup coils [Wb], Rogowski coils [A].
    type, public :: axisym_sensors_t
        type(axisym_flux_loop_t), allocatable :: flux(:)
        type(axisym_pickup_t), allocatable :: pickup(:)
        type(axisym_rogowski_t), allocatable :: rogowski(:)
    contains
        procedure :: count => sensors_count
        procedure :: signals => sensors_signals
        procedure :: response => sensors_response
    end type axisym_sensors_t

    public :: axisym_flux_field, axisym_rectangle

contains

    pure subroutine axisym_flux_field(r, z, rf, zf, wf, psi, br, bz)
        !! psi [Wb/rad] and B_R, B_Z [T] at (r, z) of the filaments (rf, zf)
        !! carrying wf [A].
        real(dp), intent(in) :: r, z, rf(:), zf(:), wf(:)
        real(dp), intent(out) :: psi, br, bz
        real(dp) :: m, q, g, gr, gz
        integer :: k
        psi = 0
        br = 0
        bz = 0
        do k = 1, size(rf)
            if (wf(k) == 0) cycle
            call axisym_green_modulus(r, z, rf(k), zf(k), m, q)
            if (q <= 0) error stop 'axisym_flux_field: field point on a filament'
            call axisym_green(r, z, rf(k), zf(k), ellip_k_from_complement(q), &
                ellip_e_from_complement(q), g, gr, gz)
            psi = psi + wf(k) * g
            br = br + wf(k) * gr
            bz = bz + wf(k) * gz
        end do
        psi = mu0 * psi
        br = mu0 * br
        bz = mu0 * bz
    end subroutine axisym_flux_field

    subroutine axisym_rectangle(rc, zc, dr, dz, turns, n, rf, zf, wf)
        !! n x n Gauss-Legendre filaments of a rectangular coil of uniform
        !! current density, centre (rc, zc), widths dr, dz [m].
        real(dp), intent(in) :: rc, zc, dr, dz, turns
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: rf(:), zf(:), wf(:)
        real(dp) :: xr(n), wr(n), xz(n), wz(n)
        integer :: i, j
        call gauss_legendre_ab(n, rc - dr / 2, rc + dr / 2, xr, wr)
        call gauss_legendre_ab(n, zc - dz / 2, zc + dz / 2, xz, wz)
        allocate(rf(n * n), zf(n * n), wf(n * n))
        do j = 1, n
            do i = 1, n
                rf(i + (j - 1) * n) = xr(i)
                zf(i + (j - 1) * n) = xz(j)
                wf(i + (j - 1) * n) = turns * wr(i) * wz(j) / (dr * dz)
            end do
        end do
    end subroutine axisym_rectangle

    pure integer function sensors_count(self) result(n)
        class(axisym_sensors_t), intent(in) :: self
        n = 0
        if (allocated(self%flux)) n = n + size(self%flux)
        if (allocated(self%pickup)) n = n + size(self%pickup)
        if (allocated(self%rogowski)) n = n + size(self%rogowski)
    end function sensors_count

    function sensors_signals(self, rf, zf, wf) result(s)
        !! Signals of the filaments (rf, zf) carrying wf [A].
        class(axisym_sensors_t), intent(in) :: self
        real(dp), intent(in) :: rf(:), zf(:), wf(:)
        real(dp), allocatable :: s(:)
        real(dp) :: psi, br, bz
        real(dp), allocatable :: x(:), w(:)
        integer :: i, j, k
        allocate(s(self%count()))
        k = 0
        if (allocated(self%flux)) then
            do i = 1, size(self%flux)
                k = k + 1
                call axisym_flux_field(self%flux(i)%r, self%flux(i)%z, rf, zf, wf, psi, br, bz)
                s(k) = 2 * pi * psi
            end do
        end if
        if (allocated(self%pickup)) then
            do i = 1, size(self%pickup)
                k = k + 1
                associate (p => self%pickup(i))
                    allocate(x(max(p%points, 1)), w(max(p%points, 1)))
                    if (p%length > 0 .and. p%points > 1) then
                        call gauss_legendre_ab(p%points, -p%length / 2, p%length / 2, x, w)
                        w = w / p%length
                    else
                        x = 0
                        w = 1
                    end if
                    s(k) = 0
                    do j = 1, size(x)
                        call axisym_flux_field(p%r + x(j) * p%tr, p%z + x(j) * p%tz, &
                            rf, zf, wf, psi, br, bz)
                        s(k) = s(k) + w(j) * (br * p%nr + bz * p%nz)
                    end do
                    s(k) = p%area * s(k)
                    deallocate(x, w)
                end associate
            end do
        end if
        if (allocated(self%rogowski)) then
            do i = 1, size(self%rogowski)
                k = k + 1
                s(k) = 0
                do j = 1, size(rf)
                    if (inside(self%rogowski(i)%r, self%rogowski(i)%z, rf(j), zf(j))) &
                        s(k) = s(k) + wf(j)
                end do
            end do
        end if
    end function sensors_signals

    function sensors_response(self, rf, zf, wf, group, ngroups) result(a)
        !! Signals per unit current of each coil group: a(:, g) are the signals
        !! of the filaments with group(c) = g and signed turns wf(c).
        class(axisym_sensors_t), intent(in) :: self
        real(dp), intent(in) :: rf(:), zf(:), wf(:)
        integer, intent(in) :: group(:), ngroups
        real(dp), allocatable :: a(:, :)
        integer :: g
        allocate(a(self%count(), ngroups))
        do g = 1, ngroups
            a(:, g) = self%signals(rf, zf, merge(wf, 0.0_dp, group == g))
        end do
    end function sensors_response

    pure logical function inside(pr, pz, r, z)
        !! Even-odd rule for the closed polygon (pr, pz).
        real(dp), intent(in) :: pr(:), pz(:), r, z
        integer :: i, j
        inside = .false.
        j = size(pr)
        do i = 1, size(pr)
            if ((pz(i) > z) .neqv. (pz(j) > z)) then
                if (r < pr(i) + (z - pz(i)) * (pr(j) - pr(i)) / (pz(j) - pz(i))) &
                    inside = .not. inside
            end if
            j = i
        end do
    end function inside

end module tiago_axisym_magnetics
