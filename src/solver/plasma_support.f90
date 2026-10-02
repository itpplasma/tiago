module tiago_plasma_support
    !! Plasma contribution to magnetic diagnostics from a VMEC equilibrium.
    !!
    !! The VMEC boundary is a flux surface (B.n = 0), so by the virtual-casing
    !! principle the field of the plasma currents outside the boundary equals the
    !! field of the surface current  mu0 K = n x B  on the boundary (n outward,
    !! B the total field there). Both the field and the vector potential follow
    !! from the Biot-Savart law of that sheet current,
    !!     A(x) = 1/(4 pi) sum_k (n x B)_k dS_k / |x - x_k|
    !!     B(x) = 1/(4 pi) sum_k (n x B)_k dS_k x (x - x_k) / |x - x_k|^3,
    !! evaluated with the trapezoidal rule on a uniform full-torus (theta, phi)
    !! grid, which converges spectrally for points away from the boundary.
    !! This is the same representation as DIAGNO's vecpot/bfield_virtual_casing.
    !!
    !! Only valid outside the plasma. Inside, the sheet-current field equals minus
    !! the coil field; a flux loop that links the plasma poloidally therefore
    !! misses the toroidal flux inside the boundary, which DIAGNO (and Tiago) add
    !! as phiedge for loops flagged idia = 1.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32, error_unit
    use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_get_var, &
        nf90_inquire_dimension, nf90_inq_dimid, nf90_noerr, nf90_nowrite
    implicit none
    private
    public :: read_vmec_boundary

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: inv_four_pi = 0.25_dp / pi
    !> Warn when a sensor point is closer to the boundary than this many grid spacings.
    real(dp), parameter :: min_distance_in_spacings = 2.0_dp

    type, public :: vmec_boundary_t
        !! Last-surface Fourier data of a VMEC equilibrium (wout conventions).
        integer :: nfp = 1, signgs = 1
        logical :: lasym = .false.
        logical :: covariant = .false., conservative = .false.
        real(dp) :: phiedge = 0.0_dp
        real(dp), allocatable :: xm(:), xn(:), xm_nyq(:), xn_nyq(:)
        real(dp), allocatable :: rmnc(:), zmns(:), rmns(:), zmnc(:)
        real(dp), allocatable :: bumnc(:), bvmnc(:), bumns(:), bvmns(:)
    end type vmec_boundary_t

    type, public :: plasma_support_t
        real(dp), allocatable :: xs(:, :)      !! (3, n) boundary points [m]
        real(dp), allocatable :: sheet(:, :)   !! (3, n) (n x B) dS / (4 pi) [T m^2]
        real(dp) :: spacing = 0.0_dp           !! largest grid spacing [m]
        real(dp) :: diamagnetic_flux = 0.0_dp  !! phiedge * signgs [Wb]
        logical :: enabled = .false.
        logical :: covariant = .false., conservative = .false.
        real(dp) :: curl_norm = 0.0_dp, current_ripple = 0.0_dp
        real(dp) :: projection_norm = 0.0_dp
        ! Boundary-field Jacobian: sheet = sum_mn b^u_mn c_mn jt + b^v_mn c_mn jp
        real(dp), allocatable :: jt(:, :), jp(:, :)  !! (3, n) dS x x_theta / (4 pi), dS x x_phi / (4 pi)
        real(dp), allocatable :: theta(:), phi(:)    !! (n) grid angles
        real(dp), allocatable :: xm_nyq(:), xn_nyq(:)
        real(dp), allocatable :: coefficients(:)     !! boundary values in mode-column order
        logical :: lasym = .false.
        ! Shape Jacobian: geometry modes and the per-point data the sheet depends on.
        real(dp), allocatable :: xm(:), xn(:)            !! geometry modes (wout xm, xn)
        real(dp), allocatable :: xt(:, :), xp(:, :)      !! (3, n) d x / d theta, d x / d phi
        real(dp), allocatable :: bu(:), bv(:)            !! (n) B^u, B^v on the boundary
        real(dp) :: sheet_scale = 0.0_dp                 !! +-dtheta dphi / (4 pi), outward
    contains
        procedure :: init_from_vmec => plasma_init_from_vmec
        procedure :: init_from_boundary => plasma_init_from_boundary
        procedure :: shape_response => plasma_shape_response
        procedure :: n_shape_columns => plasma_n_shape_columns
        procedure :: shape_column_name => plasma_shape_column_name
        procedure :: finalize => plasma_finalize
        procedure :: sample_bfield => plasma_sample_bfield
        procedure :: sample_vector_potential => plasma_sample_vector_potential
        procedure :: has_data => plasma_has_data
        procedure :: potential_weights => plasma_potential_weights
        procedure :: field_weights => plasma_field_weights
        procedure :: mode_response => plasma_mode_response
        procedure :: n_mode_columns => plasma_n_mode_columns
        procedure :: mode_column_name => plasma_mode_column_name
    end type plasma_support_t

contains

    subroutine plasma_init_from_vmec(self, wout_file, nphi, ntheta, &
            covariant, conservative)
        !! nphi: toroidal grid points per field period; ntheta: poloidal points.
        class(plasma_support_t), intent(inout) :: self
        character(len=*), intent(in) :: wout_file
        integer(i32), intent(in) :: nphi
        integer(i32), intent(in) :: ntheta
        logical, intent(in), optional :: covariant, conservative
        type(vmec_boundary_t) :: vb

        call read_vmec_boundary(wout_file, vb, covariant, conservative)
        call self%init_from_boundary(vb, nphi, ntheta)
    end subroutine plasma_init_from_vmec

    subroutine plasma_init_from_boundary(self, vb, nphi, ntheta)
        !! Geometry and field come from the boundary Fourier coefficients.
        !! Covariant mode uses (n x B) dS = (B_u x_phi - B_v x_theta) dtheta dphi.
        !! Conservative mode also projects the covariant field to a closed form.
        class(plasma_support_t), intent(inout) :: self
        type(vmec_boundary_t), intent(in) :: vb
        integer(i32), intent(in) :: nphi
        integer(i32), intent(in) :: ntheta

        type(vmec_boundary_t) :: field_boundary
        integer :: nphi_total, iphi, itheta, k
        real(dp) :: theta, phi, dtheta, dphi, volume
        real(dp) :: x(3), x_t(3), x_p(3), ds(3), b(3)

        call self%finalize()
        field_boundary = vb
        field_boundary%covariant = vb%covariant .or. vb%conservative
        self%projection_norm = 0.0_dp
        if (vb%conservative) then
            call project_covariant_boundary(field_boundary, self%projection_norm)
        end if
        self%curl_norm = 0.0_dp
        self%current_ripple = 0.0_dp
        if (field_boundary%covariant) then
            call covariant_diagnostics(field_boundary, self%curl_norm, &
                self%current_ripple)
        end if
        if (nphi < 4 .or. ntheta < 4) error stop 'plasma grid needs at least 4x4 points'
        self%diamagnetic_flux = vb%phiedge * real(vb%signgs, dp)

        nphi_total = vb%nfp * nphi
        dphi = 2.0_dp * pi / real(nphi_total, dp)
        dtheta = 2.0_dp * pi / real(ntheta, dp)
        allocate(self%xs(3, nphi_total * ntheta), self%sheet(3, nphi_total * ntheta))
        allocate(self%jt(3, nphi_total * ntheta), self%jp(3, nphi_total * ntheta))
        allocate(self%theta(nphi_total * ntheta), self%phi(nphi_total * ntheta))
        allocate(self%xt(3, nphi_total * ntheta), self%xp(3, nphi_total * ntheta))
        allocate(self%bu(nphi_total * ntheta), self%bv(nphi_total * ntheta))
        self%xm = vb%xm
        self%xn = vb%xn
        self%xm_nyq = vb%xm_nyq
        self%xn_nyq = vb%xn_nyq
        self%lasym = vb%lasym
        self%covariant = field_boundary%covariant
        self%conservative = vb%conservative
        if (vb%lasym) then
            self%coefficients = [vb%bumnc, vb%bvmnc, vb%bumns, vb%bvmns]
        else
            self%coefficients = [vb%bumnc, vb%bvmnc]
        end if

        volume = 0.0_dp
        k = 0
        do iphi = 1, nphi_total
            phi = dphi * real(iphi - 1, dp)
            do itheta = 1, ntheta
                theta = dtheta * real(itheta - 1, dp)
                k = k + 1
                call boundary_point(field_boundary, theta, phi, x, x_t, x_p, b, &
                    self%bu(k), self%bv(k))
                ds = cross(x_t, x_p) * (dtheta * dphi)
                self%xs(:, k) = x
                if (self%covariant) then
                    self%jt(:, k) = x_p * (dtheta * dphi) * inv_four_pi
                    self%jp(:, k) = -x_t * (dtheta * dphi) * inv_four_pi
                    self%sheet(:, k) = self%bu(k) * self%jt(:, k) &
                        + self%bv(k) * self%jp(:, k)
                else
                    self%sheet(:, k) = cross(ds, b) * inv_four_pi
                    self%jt(:, k) = cross(ds, x_t) * inv_four_pi
                    self%jp(:, k) = cross(ds, x_p) * inv_four_pi
                end if
                self%theta(k) = theta
                self%phi(k) = phi
                self%xt(:, k) = x_t
                self%xp(:, k) = x_p
                volume = volume + dot_product(x, ds) / 3.0_dp
                self%spacing = max(self%spacing, norm2(x_t) * dtheta, norm2(x_p) * dphi)
            end do
        end do
        ! (theta, phi) orientation decides whether x_t x x_p points outward.
        self%sheet_scale = dtheta * dphi * inv_four_pi
        if (volume < 0.0_dp) then
            self%sheet = -self%sheet
            self%jt = -self%jt
            self%jp = -self%jp
            self%sheet_scale = -self%sheet_scale
        end if
        self%enabled = .true.
    end subroutine plasma_init_from_boundary

    subroutine project_covariant_boundary(vb, correction_norm)
        !! Orthogonal projection in the unweighted Euclidean Fourier-coefficient
        !! norm: enforce n B_u + m B_v = 0; keep the two constant circulations.
        type(vmec_boundary_t), intent(inout) :: vb
        real(dp), intent(out) :: correction_norm
        real(dp) :: m, n, denom, delta
        integer :: k

        correction_norm = 0.0_dp
        do k = 1, size(vb%xm_nyq)
            m = vb%xm_nyq(k)
            n = vb%xn_nyq(k)
            denom = m * m + n * n
            if (denom == 0.0_dp) cycle
            delta = (n * vb%bumnc(k) + m * vb%bvmnc(k)) / denom
            vb%bumnc(k) = vb%bumnc(k) - n * delta
            vb%bvmnc(k) = vb%bvmnc(k) - m * delta
            correction_norm = correction_norm + denom * delta * delta
            if (vb%lasym) then
                delta = (n * vb%bumns(k) + m * vb%bvmns(k)) / denom
                vb%bumns(k) = vb%bumns(k) - n * delta
                vb%bvmns(k) = vb%bvmns(k) - m * delta
                correction_norm = correction_norm + denom * delta * delta
            end if
        end do
        correction_norm = sqrt(correction_norm)
    end subroutine project_covariant_boundary

    subroutine covariant_diagnostics(vb, curl_norm, current_ripple)
        !! L2 norm of Fourier curl coefficients [T m], and an upper bound on
        !! poloidal-section toroidal-current variation [A]. xn includes nfp.
        type(vmec_boundary_t), intent(in) :: vb
        real(dp), intent(out) :: curl_norm, current_ripple
        real(dp) :: m, n, residual
        integer :: k

        curl_norm = 0.0_dp
        current_ripple = 0.0_dp
        do k = 1, size(vb%xm_nyq)
            m = vb%xm_nyq(k)
            n = vb%xn_nyq(k)
            residual = n * vb%bumnc(k) + m * vb%bvmnc(k)
            curl_norm = curl_norm + residual * residual
            if (m == 0.0_dp .and. n /= 0.0_dp) then
                current_ripple = current_ripple + abs(vb%bumnc(k))
            end if
            if (vb%lasym) then
                residual = n * vb%bumns(k) + m * vb%bvmns(k)
                curl_norm = curl_norm + residual * residual
                if (m == 0.0_dp .and. n /= 0.0_dp) then
                    current_ripple = current_ripple + abs(vb%bumns(k))
                end if
            end if
        end do
        curl_norm = sqrt(curl_norm)
        current_ripple = current_ripple / 2.0e-7_dp
    end subroutine covariant_diagnostics

    subroutine boundary_point(vb, theta, phi, x, x_t, x_p, b, bu, bv)
        !! Position, tangents d/dtheta, d/dphi [m], total field [T] and its
        !! selected contravariant or covariant components at s = 1.
        type(vmec_boundary_t), intent(in) :: vb
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: x(3), x_t(3), x_p(3), b(3), bu, bv
        real(dp) :: r, z, r_t, r_p, z_t, z_p, c, s, guu, guv, gvv, det
        real(dp) :: arg(size(vb%xm)), arg_nyq(size(vb%xm_nyq))

        arg = vb%xm * theta - vb%xn * phi
        arg_nyq = vb%xm_nyq * theta - vb%xn_nyq * phi
        r = sum(vb%rmnc * cos(arg))
        z = sum(vb%zmns * sin(arg))
        r_t = -sum(vb%rmnc * vb%xm * sin(arg))
        r_p = sum(vb%rmnc * vb%xn * sin(arg))
        z_t = sum(vb%zmns * vb%xm * cos(arg))
        z_p = -sum(vb%zmns * vb%xn * cos(arg))
        bu = sum(vb%bumnc * cos(arg_nyq))
        bv = sum(vb%bvmnc * cos(arg_nyq))
        if (vb%lasym) then
            r = r + sum(vb%rmns * sin(arg))
            z = z + sum(vb%zmnc * cos(arg))
            r_t = r_t + sum(vb%rmns * vb%xm * cos(arg))
            r_p = r_p - sum(vb%rmns * vb%xn * cos(arg))
            z_t = z_t - sum(vb%zmnc * vb%xm * sin(arg))
            z_p = z_p + sum(vb%zmnc * vb%xn * sin(arg))
            bu = bu + sum(vb%bumns * sin(arg_nyq))
            bv = bv + sum(vb%bvmns * sin(arg_nyq))
        end if
        c = cos(phi)
        s = sin(phi)
        x = [r * c, r * s, z]
        x_t = [r_t * c, r_t * s, z_t]
        x_p = [r_p * c - r * s, r_p * s + r * c, z_p]
        if (vb%covariant) then
            guu = dot_product(x_t, x_t)
            guv = dot_product(x_t, x_p)
            gvv = dot_product(x_p, x_p)
            det = guu * gvv - guv * guv
            b = ((gvv * bu - guv * bv) * x_t &
                + (guu * bv - guv * bu) * x_p) / det
        else
            b = bu * x_t + bv * x_p
        end if
    end subroutine boundary_point

    subroutine read_vmec_boundary(path, vb, covariant, conservative)
        !! Last-surface Fourier data from a VMEC wout file. B^u, B^v live on the
        !! half mesh and are extrapolated to s = 1 as 1.5 b(ns) - 0.5 b(ns-1).
        character(len=*), intent(in) :: path
        type(vmec_boundary_t), intent(out) :: vb
        logical, intent(in), optional :: covariant, conservative
        character(len=8) :: bu_name, bv_name, bus_name, bvs_name
        integer :: ncid, ns, mnmax, mnmax_nyq, lasym_int
        real(dp), allocatable :: buf(:, :)

        if (present(covariant)) vb%covariant = covariant
        if (present(conservative)) vb%conservative = conservative
        vb%covariant = vb%covariant .or. vb%conservative
        bu_name = 'bsupumnc'
        bv_name = 'bsupvmnc'
        bus_name = 'bsupumns'
        bvs_name = 'bsupvmns'
        if (vb%covariant) then
            bu_name = 'bsubumnc'
            bv_name = 'bsubvmnc'
            bus_name = 'bsubumns'
            bvs_name = 'bsubvmns'
        end if
        call nc(nf90_open(trim(path), nf90_nowrite, ncid), 'open '//trim(path))
        ns = dim_len(ncid, 'radius')
        mnmax = dim_len(ncid, 'mn_mode')
        mnmax_nyq = dim_len(ncid, 'mn_mode_nyq')
        call get_int(ncid, 'nfp', vb%nfp)
        call get_int(ncid, 'signgs', vb%signgs)
        call get_int(ncid, 'lasym__logical__', lasym_int)
        vb%lasym = lasym_int /= 0
        allocate(vb%xm(mnmax), vb%xn(mnmax), vb%xm_nyq(mnmax_nyq), vb%xn_nyq(mnmax_nyq))
        call get_1d(ncid, 'xm', vb%xm)
        call get_1d(ncid, 'xn', vb%xn)
        call get_1d(ncid, 'xm_nyq', vb%xm_nyq)
        call get_1d(ncid, 'xn_nyq', vb%xn_nyq)
        allocate(buf(mnmax, ns))
        call get_2d(ncid, 'rmnc', buf); vb%rmnc = buf(:, ns)
        call get_2d(ncid, 'zmns', buf); vb%zmns = buf(:, ns)
        if (vb%lasym) then
            call get_2d(ncid, 'rmns', buf); vb%rmns = buf(:, ns)
            call get_2d(ncid, 'zmnc', buf); vb%zmnc = buf(:, ns)
        end if
        deallocate(buf)
        allocate(buf(mnmax_nyq, ns))
        call get_2d(ncid, bu_name, buf)
        vb%bumnc = 1.5_dp * buf(:, ns) - 0.5_dp * buf(:, ns - 1)
        call get_2d(ncid, bv_name, buf)
        vb%bvmnc = 1.5_dp * buf(:, ns) - 0.5_dp * buf(:, ns - 1)
        if (vb%lasym) then
            call get_2d(ncid, bus_name, buf)
            vb%bumns = 1.5_dp * buf(:, ns) - 0.5_dp * buf(:, ns - 1)
            call get_2d(ncid, bvs_name, buf)
            vb%bvmns = 1.5_dp * buf(:, ns) - 0.5_dp * buf(:, ns - 1)
        end if
        deallocate(buf)
        allocate(buf(ns, 1))
        call get_1d(ncid, 'phi', buf(:, 1))
        vb%phiedge = buf(ns, 1)
        call nc(nf90_close(ncid), 'close')
    contains
        integer function dim_len(id, name)
            integer, intent(in) :: id
            character(len=*), intent(in) :: name
            integer :: dimid
            call nc(nf90_inq_dimid(id, name, dimid), name)
            call nc(nf90_inquire_dimension(id, dimid, len=dim_len), name)
        end function dim_len
        subroutine get_int(id, name, value)
            integer, intent(in) :: id
            character(len=*), intent(in) :: name
            integer, intent(out) :: value
            integer :: varid
            call nc(nf90_inq_varid(id, name, varid), name)
            call nc(nf90_get_var(id, varid, value), name)
        end subroutine get_int
        subroutine get_1d(id, name, values)
            integer, intent(in) :: id
            character(len=*), intent(in) :: name
            real(dp), intent(out) :: values(:)
            integer :: varid
            call nc(nf90_inq_varid(id, name, varid), name)
            call nc(nf90_get_var(id, varid, values), name)
        end subroutine get_1d
        subroutine get_2d(id, name, values)
            integer, intent(in) :: id
            character(len=*), intent(in) :: name
            real(dp), intent(out) :: values(:, :)
            integer :: varid
            call nc(nf90_inq_varid(id, name, varid), name)
            call nc(nf90_get_var(id, varid, values), name)
        end subroutine get_2d
        subroutine nc(status, what)
            integer, intent(in) :: status
            character(len=*), intent(in) :: what
            if (status /= nf90_noerr) then
                write(error_unit, '(A)') 'netCDF error reading VMEC '//what//' from '//trim(path)
                error stop 1
            end if
        end subroutine nc
    end subroutine read_vmec_boundary

    subroutine plasma_finalize(self)
        class(plasma_support_t), intent(inout) :: self
        if (allocated(self%xs)) deallocate(self%xs)
        if (allocated(self%sheet)) deallocate(self%sheet)
        if (allocated(self%jt)) deallocate(self%jt, self%jp, self%theta, self%phi)
        if (allocated(self%xt)) deallocate(self%xt, self%xp, self%bu, self%bv)
        self%spacing = 0.0_dp
        self%diamagnetic_flux = 0.0_dp
        self%enabled = .false.
    end subroutine plasma_finalize

    subroutine plasma_sample_bfield(self, points, bfield)
        !! points(3, n) [m] -> plasma field bfield(3, n) [T]
        class(plasma_support_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: bfield(:, :)
        integer :: i, k
        real(dp) :: d(3), r2, r, rmin, acc(3)

        bfield = 0.0_dp
        if (.not. self%enabled) return
        rmin = huge(1.0_dp)
!$omp parallel do schedule(static) private(i, k, d, r2, r, acc) reduction(min:rmin)
        do i = 1, size(points, 2)
            acc = 0.0_dp
            do k = 1, size(self%xs, 2)
                d = points(:, i) - self%xs(:, k)
                r2 = dot_product(d, d)
                r = sqrt(r2)
                rmin = min(rmin, r)
                acc = acc + cross(self%sheet(:, k), d) / (r2 * r)
            end do
            bfield(:, i) = acc
        end do
!$omp end parallel do
        call warn_if_close(self, rmin)
    end subroutine plasma_sample_bfield

    subroutine plasma_sample_vector_potential(self, points, avec)
        !! points(3, n) [m] -> plasma vector potential avec(3, n) [T m]
        class(plasma_support_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: avec(:, :)
        integer :: i, k
        real(dp) :: d(3), r, rmin, acc(3)

        avec = 0.0_dp
        if (.not. self%enabled) return
        rmin = huge(1.0_dp)
!$omp parallel do schedule(static) private(i, k, d, r, acc) reduction(min:rmin)
        do i = 1, size(points, 2)
            acc = 0.0_dp
            do k = 1, size(self%xs, 2)
                d = points(:, i) - self%xs(:, k)
                r = sqrt(dot_product(d, d))
                rmin = min(rmin, r)
                acc = acc + self%sheet(:, k) / r
            end do
            avec(:, i) = acc
        end do
!$omp end parallel do
        call warn_if_close(self, rmin)
    end subroutine plasma_sample_vector_potential

    subroutine plasma_potential_weights(self, points, dls, owner, nsig, w, gx)
        !! A-based signals s = sum_j A(x_j) . dl_j = sum_k sheet_k . w(:, k, s),
        !! w(:, k, s) = sum_{j in s} dl_j / |x_j - x_k|.
        !! gx(:, k, s) = d s / d x_k at fixed sheet (for the shape Jacobian).
        class(plasma_support_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :), dls(:, :)
        integer, intent(in) :: owner(:), nsig
        real(dp), allocatable, intent(out) :: w(:, :, :)
        real(dp), allocatable, intent(out), optional :: gx(:, :, :)
        integer :: j, k
        real(dp) :: d(3), r

        allocate(w(3, size(self%xs, 2), nsig))
        w = 0.0_dp
        if (present(gx)) then
            allocate(gx(3, size(self%xs, 2), nsig))
            gx = 0.0_dp
        end if
!$omp parallel do schedule(static) private(j, k, d, r)
        do k = 1, size(self%xs, 2)
            do j = 1, size(owner)
                d = points(:, j) - self%xs(:, k)
                r = norm2(d)
                w(:, k, owner(j)) = w(:, k, owner(j)) + dls(:, j) / r
                if (present(gx)) gx(:, k, owner(j)) = gx(:, k, owner(j)) + &
                    dot_product(self%sheet(:, k), dls(:, j)) * d / (r * r * r)
            end do
        end do
!$omp end parallel do
    end subroutine plasma_potential_weights

    subroutine plasma_field_weights(self, points, dls, owner, nsig, w, gx)
        !! B-based signals s = sum_j B(x_j) . dl_j = sum_k sheet_k . w(:, k, s),
        !! w(:, k, s) = sum_{j in s} (x_j - x_k) x dl_j / |x_j - x_k|^3.
        !! gx(:, k, s) = d s / d x_k at fixed sheet (for the shape Jacobian).
        class(plasma_support_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :), dls(:, :)
        integer, intent(in) :: owner(:), nsig
        real(dp), allocatable, intent(out) :: w(:, :, :)
        real(dp), allocatable, intent(out), optional :: gx(:, :, :)
        integer :: j, k
        real(dp) :: d(3), r, v(3)

        allocate(w(3, size(self%xs, 2), nsig))
        w = 0.0_dp
        if (present(gx)) then
            allocate(gx(3, size(self%xs, 2), nsig))
            gx = 0.0_dp
        end if
!$omp parallel do schedule(static) private(j, k, d, r, v)
        do k = 1, size(self%xs, 2)
            do j = 1, size(owner)
                d = points(:, j) - self%xs(:, k)
                r = norm2(d)
                w(:, k, owner(j)) = w(:, k, owner(j)) + cross(d, dls(:, j)) / (r * r * r)
                if (present(gx)) then
                    ! s = d . (dl x sheet) / r^3 with d = x_j - x_k
                    v = cross(dls(:, j), self%sheet(:, k))
                    gx(:, k, owner(j)) = gx(:, k, owner(j)) - v / r**3 + &
                        3.0_dp * dot_product(d, v) * d / r**5
                end if
            end do
        end do
!$omp end parallel do
    end subroutine plasma_field_weights

    integer function plasma_n_shape_columns(self)
        class(plasma_support_t), intent(in) :: self
        plasma_n_shape_columns = merge(4, 2, self%lasym) * size(self%xm)
    end function plasma_n_shape_columns

    subroutine plasma_shape_column_name(self, column, coefficient, m, n)
        !! Column c -> VMEC boundary geometry coefficient name and mode numbers.
        class(plasma_support_t), intent(in) :: self
        integer, intent(in) :: column
        character(len=:), allocatable, intent(out) :: coefficient
        integer, intent(out) :: m, n
        character(len=4), parameter :: names(4) = [character(len=4) :: 'rmnc', 'zmns', 'rmns', 'zmnc']
        integer :: nm

        nm = size(self%xm)
        coefficient = trim(names((column - 1) / nm + 1))
        m = nint(self%xm(mod(column - 1, nm) + 1))
        n = nint(self%xn(mod(column - 1, nm) + 1))
    end subroutine plasma_shape_column_name

    subroutine plasma_shape_response(self, w, gx, resp)
        !! resp(s, c): derivative of the plasma part of signal s with respect to
        !! the boundary geometry coefficient c (rmnc, zmns[, rmns, zmnc] at s = 1),
        !! at fixed input field coefficients of the selected representation.
        !! Conservative mode uses their fixed spectral projection. w and gx are
        !! the weights and kernel-position cotangents of the signals.
        !! Reverse mode through sheet = c (alpha x_p - beta x_t),
        !! alpha = x_t . b, beta = x_p . b, b = B^u x_t + B^v x_p.
        class(plasma_support_t), intent(in) :: self
        real(dp), intent(in) :: w(:, :, :), gx(:, :, :)
        real(dp), allocatable, intent(out) :: resp(:, :)
        real(dp), allocatable :: rb(:, :), rtb(:, :), rpb(:, :), zb(:, :), ztb(:, :), zpb(:, :)
        real(dp) :: xtb(3), xpb(3), xb(3), alpha, beta, abar, bbar, c, cp, sp
        real(dp) :: xt(3), xp(3), bu, bv, arg, co, si, m, n
        integer :: k, s, mode, nm, npts, nsig

        npts = size(self%xs, 2)
        nsig = size(w, 3)
        nm = size(self%xm)
        c = self%sheet_scale
        allocate(rb(nsig, npts), rtb(nsig, npts), rpb(nsig, npts))
        allocate(zb(nsig, npts), ztb(nsig, npts), zpb(nsig, npts))
!$omp parallel do schedule(static) private(k, s, xt, xp, bu, bv, alpha, beta, xtb, xpb, xb, &
!$omp abar, bbar, cp, sp)
        do k = 1, npts
            xt = self%xt(:, k)
            xp = self%xp(:, k)
            bu = self%bu(k)
            bv = self%bv(k)
            if (self%covariant) then
                alpha = bu
                beta = bv
            else
                alpha = bu * dot_product(xt, xt) + bv * dot_product(xt, xp)
                beta = bu * dot_product(xt, xp) + bv * dot_product(xp, xp)
            end if
            cp = cos(self%phi(k))
            sp = sin(self%phi(k))
            do s = 1, nsig
                xpb = c * alpha * w(:, k, s)
                xtb = -c * beta * w(:, k, s)
                if (.not. self%covariant) then
                    abar = c * dot_product(w(:, k, s), xp)
                    bbar = -c * dot_product(w(:, k, s), xt)
                    xtb = xtb + abar * (2.0_dp * bu * xt + bv * xp) &
                        + bbar * bu * xp
                    xpb = xpb + abar * bv * xt &
                        + bbar * (bu * xt + 2.0_dp * bv * xp)
                end if
                xb = gx(:, k, s)
                ! x = (R c, R s, Z), x_t = (R_t c, R_t s, Z_t),
                ! x_p = (R_p c - R s, R_p s + R c, Z_p)
                rb(s, k) = xb(1) * cp + xb(2) * sp - xpb(1) * sp + xpb(2) * cp
                rtb(s, k) = xtb(1) * cp + xtb(2) * sp
                rpb(s, k) = xpb(1) * cp + xpb(2) * sp
                zb(s, k) = xb(3)
                ztb(s, k) = xtb(3)
                zpb(s, k) = xpb(3)
            end do
        end do
!$omp end parallel do
        allocate(resp(nsig, self%n_shape_columns()))
        resp = 0.0_dp
!$omp parallel do schedule(static) private(mode, k, m, n, arg, co, si)
        do mode = 1, nm
            m = self%xm(mode)
            n = self%xn(mode)
            do k = 1, npts
                arg = m * self%theta(k) - n * self%phi(k)
                co = cos(arg)
                si = sin(arg)
                ! R = sum rmnc cos, Z = sum zmns sin (and rmns sin, zmnc cos)
                resp(:, mode) = resp(:, mode) + rb(:, k) * co - m * rtb(:, k) * si &
                    + n * rpb(:, k) * si
                resp(:, nm + mode) = resp(:, nm + mode) + zb(:, k) * si + m * ztb(:, k) * co &
                    - n * zpb(:, k) * co
                if (self%lasym) then
                    resp(:, 2 * nm + mode) = resp(:, 2 * nm + mode) + rb(:, k) * si &
                        + m * rtb(:, k) * co - n * rpb(:, k) * co
                    resp(:, 3 * nm + mode) = resp(:, 3 * nm + mode) + zb(:, k) * co &
                        - m * ztb(:, k) * si + n * zpb(:, k) * si
                end if
            end do
        end do
!$omp end parallel do
    end subroutine plasma_shape_response

    integer function plasma_n_mode_columns(self)
        class(plasma_support_t), intent(in) :: self
        plasma_n_mode_columns = merge(4, 2, self%lasym) * size(self%xm_nyq)
    end function plasma_n_mode_columns

    subroutine plasma_mode_column_name(self, column, coefficient, m, n)
        !! Column c -> VMEC boundary coefficient name and mode numbers (xm, xn).
        class(plasma_support_t), intent(in) :: self
        integer, intent(in) :: column
        character(len=:), allocatable, intent(out) :: coefficient
        integer, intent(out) :: m, n
        character(len=8) :: names(4)
        integer :: nm

        names = [character(len=8) :: 'bsupumnc', 'bsupvmnc', &
            'bsupumns', 'bsupvmns']
        if (self%covariant) names = [character(len=8) :: &
            'bsubumnc', 'bsubvmnc', 'bsubumns', 'bsubvmns']
        nm = size(self%xm_nyq)
        coefficient = trim(names((column - 1) / nm + 1))
        m = nint(self%xm_nyq(mod(column - 1, nm) + 1))
        n = nint(self%xn_nyq(mod(column - 1, nm) + 1))
    end subroutine plasma_mode_column_name

    subroutine plasma_mode_response(self, w, resp)
        !! resp(s, c): derivative of signal s with respect to the boundary value
        !! (s = 1, half-mesh extrapolated) of coefficient column c. The plasma
        !! part of every signal is exactly sum_c resp(s, c) * coefficient_c.
        class(plasma_support_t), intent(in) :: self
        real(dp), intent(in) :: w(:, :, :)
        real(dp), allocatable, intent(out) :: resp(:, :)
        real(dp), allocatable :: cu(:), cv(:)
        real(dp) :: arg, m, n, denom
        integer :: k, mode, s, nm

        nm = size(self%xm_nyq)
        allocate(resp(size(w, 3), self%n_mode_columns()), cu(size(w, 3)), cv(size(w, 3)))
        resp = 0.0_dp
!$omp parallel do schedule(static) private(mode, k, s, arg, cu, cv, m, n, denom)
        do mode = 1, nm
            do k = 1, size(self%xs, 2)
                arg = self%xm_nyq(mode) * self%theta(k) - self%xn_nyq(mode) * self%phi(k)
                do s = 1, size(w, 3)
                    cu(s) = dot_product(w(:, k, s), self%jt(:, k))
                    cv(s) = dot_product(w(:, k, s), self%jp(:, k))
                end do
                resp(:, mode) = resp(:, mode) + cu * cos(arg)
                resp(:, nm + mode) = resp(:, nm + mode) + cv * cos(arg)
                if (self%lasym) then
                    resp(:, 2 * nm + mode) = resp(:, 2 * nm + mode) + cu * sin(arg)
                    resp(:, 3 * nm + mode) = resp(:, 3 * nm + mode) + cv * sin(arg)
                end if
            end do
            if (self%conservative) then
                m = self%xm_nyq(mode)
                n = self%xn_nyq(mode)
                denom = m * m + n * n
                if (denom > 0.0_dp) then
                    cu = (n * resp(:, mode) + m * resp(:, nm + mode)) / denom
                    resp(:, mode) = resp(:, mode) - n * cu
                    resp(:, nm + mode) = resp(:, nm + mode) - m * cu
                    if (self%lasym) then
                        cu = (n * resp(:, 2 * nm + mode) &
                            + m * resp(:, 3 * nm + mode)) / denom
                        resp(:, 2 * nm + mode) = resp(:, 2 * nm + mode) - n * cu
                        resp(:, 3 * nm + mode) = resp(:, 3 * nm + mode) - m * cu
                    end if
                end if
            end if
        end do
!$omp end parallel do
    end subroutine plasma_mode_response

    subroutine warn_if_close(self, rmin)
        type(plasma_support_t), intent(in) :: self
        real(dp), intent(in) :: rmin

        if (rmin < min_distance_in_spacings * self%spacing) then
            write(error_unit, '(A,ES9.2,A,ES9.2,A)') 'WARNING: sensor point ', rmin, &
                ' m from the plasma boundary (grid spacing ', self%spacing, &
                ' m); increase --plasma-nphi/--plasma-ntheta'
        end if
    end subroutine warn_if_close

    logical function plasma_has_data(self)
        class(plasma_support_t), intent(in) :: self
        plasma_has_data = self%enabled
    end function plasma_has_data

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
    end function cross
end module tiago_plasma_support
