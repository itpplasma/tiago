module tiago_surface_biot_savart
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    use interpolate, only: SplineData2D, construct_splines_2d, &
        destroy_splines_2d, evaluate_splines_2d
    implicit none
    private

    integer, parameter :: quad_order = 4

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: two_pi = 2.0_dp * pi
    real(dp), parameter :: mu0 = 4.0_dp * pi * 1.0e-7_dp
    real(dp), parameter :: inv_mu0 = 1.0_dp / mu0
    real(dp), parameter :: gauss_nodes(quad_order) = &
        (/ -0.8611363115940526_dp, -0.3399810435848563_dp, &
           0.3399810435848563_dp,  0.8611363115940526_dp /)
    real(dp), parameter :: gauss_weights(quad_order) = &
        (/ 0.3478548451374539_dp, 0.6521451548625461_dp, &
           0.6521451548625461_dp, 0.3478548451374539_dp /)

    type :: surface_biot_savart_t
        type(SplineData2D) :: x_spline(3)
        type(SplineData2D) :: k_spline(3)
        type(SplineData2D) :: normal_spline(3)
        type(SplineData2D) :: bn_spline
        type(SplineData2D) :: jacobian_spline
        real(dp) :: phi_min = 0.0_dp
        real(dp) :: phi_max = 0.0_dp
        real(dp) :: theta_min = 0.0_dp
        real(dp) :: theta_max = 0.0_dp
        real(dp) :: dphi = 0.0_dp
        real(dp) :: dtheta = 0.0_dp
        integer(i32) :: nphi = 0_i32
        integer(i32) :: ntheta = 0_i32
        real(dp) :: norm_const = 0.0_dp
        logical :: initialized = .false.
    contains
        procedure :: init => biot_savart_init
        procedure :: finalize => biot_savart_finalize
        procedure :: evaluate_b => biot_savart_eval_b
        procedure :: evaluate_a => biot_savart_eval_a
        procedure :: sample_b_field => biot_savart_sample_b
        procedure :: sample_vector_potential => biot_savart_sample_a
    end type surface_biot_savart_t

    public :: surface_biot_savart_t

contains

    subroutine biot_savart_init(self, x_surf, b_total, nfp)
        class(surface_biot_savart_t), intent(inout) :: self
        real(dp), intent(in) :: x_surf(:, :, :)
        real(dp), intent(in) :: b_total(:, :, :)
        integer, intent(in) :: nfp
        real(dp) :: nfp_real

        real(dp), allocatable :: normals(:, :, :)
        real(dp), allocatable :: k(:, :, :)
        real(dp), allocatable :: bn(:, :)
        real(dp), allocatable :: jac(:, :)
        real(dp), allocatable :: x_extended(:, :, :)
        real(dp), allocatable :: k_extended(:, :, :)
        real(dp), allocatable :: normals_extended(:, :, :)
        real(dp), allocatable :: bn_extended(:, :)
        real(dp), allocatable :: jac_extended(:, :)
        integer :: phi_count
        integer :: theta_count
        integer :: phi_count_extended
        integer :: i
        integer :: j
        integer :: ifp
        integer :: ii
        integer :: idx_start
        integer :: idx_end
        real(dp) :: phi_step
        real(dp) :: theta_step
        real(dp) :: phi_period
        real(dp) :: theta_period
        integer :: ip
        integer :: im
        integer :: jp
        integer :: jm
        real(dp) :: dx_dphi(3)
        real(dp) :: dx_dtheta(3)
        real(dp) :: normal_vec(3)
        real(dp) :: norm_mag
        integer :: spline_order(2)
        logical :: periodic(2)
        real(dp) :: limits_min(2)
        real(dp) :: limits_max(2)
        real(dp) :: factor
        real(dp) :: cop
        real(dp) :: sip

        if (size(x_surf, 1) < 4 .or. size(x_surf, 2) < 4) then
            error stop 'surface_biot_savart requires at least 4x4 surface grid'
        end if
        if (size(x_surf, 3) /= 3 .or. size(b_total, 3) /= 3) then
            error stop 'surface_biot_savart expects 3D Cartesian arrays'
        end if
        if (size(x_surf, 1) /= size(b_total, 1) .or. &
            size(x_surf, 2) /= size(b_total, 2)) then
            error stop 'surface_biot_savart grid mismatch between x_surf and b_total'
        end if

        call self%finalize()

        phi_count = size(x_surf, 1)
        theta_count = size(x_surf, 2)

        allocate(normals(phi_count, theta_count, 3))
        allocate(k(phi_count, theta_count, 3))
        allocate(bn(phi_count, theta_count))
        allocate(jac(phi_count, theta_count))

        phi_period = two_pi / real(2 * nfp, dp)
        theta_period = two_pi
        phi_step = phi_period / real(phi_count, dp)
        theta_step = theta_period / real(theta_count, dp)

        do i = 1, phi_count
            ip = i + 1
            if (ip > phi_count) ip = 1
            im = i - 1
            if (im < 1) im = phi_count

            do j = 1, theta_count
                jp = j + 1
                if (jp > theta_count) jp = 1
                jm = j - 1
                if (jm < 1) jm = theta_count

                dx_dphi(:) = (x_surf(ip, j, :) - x_surf(im, j, :)) / (2.0_dp * phi_step)
                dx_dtheta(:) = (x_surf(i, jp, :) - x_surf(i, jm, :)) / (2.0_dp * theta_step)

                normal_vec = cross_product(dx_dphi, dx_dtheta)
                norm_mag = sqrt(sum(normal_vec**2))
                if (norm_mag < 1.0e-12_dp) error stop 'degenerate surface normal encountered'
                normals(i, j, :) = normal_vec / norm_mag
                jac(i, j) = norm_mag

                k(i, j, :) = cross_product(b_total(i, j, :), normal_vec)
                bn(i, j) = sum(normal_vec * b_total(i, j, :))
            end do
        end do

        phi_count_extended = phi_count * nfp
        allocate(x_extended(phi_count_extended, theta_count, 3))
        allocate(k_extended(phi_count_extended, theta_count, 3))
        allocate(normals_extended(phi_count_extended, theta_count, 3))
        allocate(bn_extended(phi_count_extended, theta_count))
        allocate(jac_extended(phi_count_extended, theta_count))

        x_extended(1:phi_count, :, :) = x_surf
        k_extended(1:phi_count, :, :) = k
        normals_extended(1:phi_count, :, :) = normals
        bn_extended(1:phi_count, :) = bn
        jac_extended(1:phi_count, :) = jac

        factor = two_pi / real(nfp, dp)
        do ifp = 2, nfp
            cop = cos(real(ifp - 1, dp) * factor)
            sip = sin(real(ifp - 1, dp) * factor)
            idx_start = (ifp - 1) * phi_count + 1
            idx_end = ifp * phi_count

            do j = 1, theta_count
                do i = 1, phi_count
                    ii = idx_start + i - 1
                    x_extended(ii, j, 1) = x_surf(i, j, 1) * cop - x_surf(i, j, 2) * sip
                    x_extended(ii, j, 2) = x_surf(i, j, 1) * sip + x_surf(i, j, 2) * cop
                    x_extended(ii, j, 3) = x_surf(i, j, 3)

                    k_extended(ii, j, 1) = k(i, j, 1) * cop - k(i, j, 2) * sip
                    k_extended(ii, j, 2) = k(i, j, 1) * sip + k(i, j, 2) * cop
                    k_extended(ii, j, 3) = k(i, j, 3)

                    normals_extended(ii, j, 1) = normals(i, j, 1) * cop - normals(i, j, 2) * sip
                    normals_extended(ii, j, 2) = normals(i, j, 1) * sip + normals(i, j, 2) * cop
                    normals_extended(ii, j, 3) = normals(i, j, 3)

                    bn_extended(ii, j) = bn(i, j)
                    jac_extended(ii, j) = jac(i, j)
                end do
            end do
        end do

        limits_min = [0.0_dp, 0.0_dp]
        limits_max = [two_pi, theta_period]
        periodic = [.true., .true.]
        spline_order = [3, 3]

        do i = 1, 3
            call construct_splines_2d(limits_min, limits_max, x_extended(:, :, i), &
                spline_order, periodic, self%x_spline(i))
            call construct_splines_2d(limits_min, limits_max, k_extended(:, :, i), &
                spline_order, periodic, self%k_spline(i))
            call construct_splines_2d(limits_min, limits_max, normals_extended(:, :, i), &
                spline_order, periodic, self%normal_spline(i))
        end do
        call construct_splines_2d(limits_min, limits_max, bn_extended, &
            spline_order, periodic, self%bn_spline)
        call construct_splines_2d(limits_min, limits_max, jac_extended, &
            spline_order, periodic, self%jacobian_spline)

        self%phi_min = 0.0_dp
        self%phi_max = two_pi
        self%theta_min = 0.0_dp
        self%theta_max = theta_period
        self%nphi = phi_count_extended
        self%ntheta = theta_count
        self%dphi = two_pi / real(phi_count_extended, dp)
        self%dtheta = theta_step
        nfp_real = real(nfp, dp)
        self%norm_const = 1.0_dp / (2.0_dp * pi * pi * real(phi_count * theta_count, dp))
        self%initialized = .true.

        deallocate(normals, k, bn, jac)
        deallocate(x_extended, k_extended, normals_extended, bn_extended, jac_extended)
    end subroutine biot_savart_init

    subroutine biot_savart_finalize(self)
        class(surface_biot_savart_t), intent(inout) :: self
        integer :: i
        do i = 1, 3
            call destroy_splines_2d(self%x_spline(i))
            call destroy_splines_2d(self%k_spline(i))
            call destroy_splines_2d(self%normal_spline(i))
        end do
        call destroy_splines_2d(self%bn_spline)
        call destroy_splines_2d(self%jacobian_spline)
        self%phi_min = 0.0_dp
        self%phi_max = 0.0_dp
        self%theta_min = 0.0_dp
        self%theta_max = 0.0_dp
        self%dphi = 0.0_dp
        self%dtheta = 0.0_dp
        self%nphi = 0_i32
        self%ntheta = 0_i32
        self%norm_const = 0.0_dp
        self%initialized = .false.
    end subroutine biot_savart_finalize

    subroutine biot_savart_sample_b(self, points, bfield)
        class(surface_biot_savart_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: bfield(:, :)
        integer :: idx
        if (.not. self%initialized) then
            bfield = 0.0_dp
            return
        end if
        if (size(points, 2) /= 3) error stop 'biot_savart_sample_b expects 3D points'
        if (size(bfield, 1) /= size(points, 1) .or. size(bfield, 2) /= 3) then
            error stop 'biot_savart_sample_b output array mismatch'
        end if

        do idx = 1, size(points, 1)
            call biot_savart_eval_b(self, points(idx, :), bfield(idx, :))
        end do
    end subroutine biot_savart_sample_b

    subroutine biot_savart_sample_a(self, points, avec)
        class(surface_biot_savart_t), intent(in) :: self
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: avec(:, :)
        integer :: idx
        if (.not. self%initialized) then
            avec = 0.0_dp
            return
        end if
        if (size(points, 2) /= 3) error stop 'biot_savart_sample_a expects 3D points'
        if (size(avec, 1) /= size(points, 1) .or. size(avec, 2) /= 3) then
            error stop 'biot_savart_sample_a output array mismatch'
        end if

        do idx = 1, size(points, 1)
            call biot_savart_eval_a(self, points(idx, :), avec(idx, :))
        end do
    end subroutine biot_savart_sample_a

    subroutine biot_savart_eval_b(self, point, bvec)
        class(surface_biot_savart_t), intent(in) :: self
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: bvec(3)

        call integrate_kernel(self, point, bvec, .true.)
    end subroutine biot_savart_eval_b

    subroutine biot_savart_eval_a(self, point, avec)
        class(surface_biot_savart_t), intent(in) :: self
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: avec(3)

        call integrate_kernel(self, point, avec, .false.)
    end subroutine biot_savart_eval_a

    subroutine integrate_kernel(self, point, output, compute_b)
        class(surface_biot_savart_t), intent(in) :: self
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: output(3)
        logical, intent(in) :: compute_b

        integer :: i
        integer :: j
        integer :: gp
        integer :: gt
        real(dp) :: phi_a
        real(dp) :: phi_b
        real(dp) :: theta_a
        real(dp) :: theta_b
        real(dp) :: phi_value
        real(dp) :: theta_value
        real(dp) :: weight_phi
        real(dp) :: weight_theta
        real(dp) :: surface_point(3)
        real(dp) :: kvec(3)
        real(dp) :: normal_vec(3)
        real(dp) :: jacobian
        real(dp) :: bn_value
        real(dp) :: r_vec(3)
        real(dp) :: dist2
        real(dp) :: dist
        real(dp) :: inv_r
        real(dp) :: inv_r3
        real(dp) :: contrib(3)
        real(dp) :: weight

        output = 0.0_dp
        if (.not. self%initialized) return

        do i = 1, self%nphi
            phi_a = self%phi_min + self%dphi * real(i - 1, dp)
            phi_b = phi_a + self%dphi
            do j = 1, self%ntheta
                theta_a = self%theta_min + self%dtheta * real(j - 1, dp)
                theta_b = theta_a + self%dtheta

                do gp = 1, quad_order
                    phi_value = 0.5_dp * ((phi_b - phi_a) * gauss_nodes(gp) + phi_b + phi_a)
                    weight_phi = 0.5_dp * (phi_b - phi_a) * gauss_weights(gp)
                    do gt = 1, quad_order
                        theta_value = 0.5_dp * ((theta_b - theta_a) * gauss_nodes(gt) + theta_b + theta_a)
                        weight_theta = 0.5_dp * (theta_b - theta_a) * gauss_weights(gt)

                        call evaluate_surface_state(self, phi_value, theta_value, &
                            surface_point, kvec, normal_vec, bn_value, jacobian)

                        r_vec = point - surface_point
                        dist2 = sum(r_vec**2)
                        if (dist2 < 1.0e-20_dp) cycle
                        dist = sqrt(dist2)
                        inv_r = 1.0_dp / dist
                        weight = self%norm_const * weight_phi * weight_theta * jacobian

                        if (compute_b) then
                            inv_r3 = inv_r / dist2
                            contrib = cross_product(kvec, r_vec) + bn_value * r_vec
                            output = output + weight * contrib * inv_r3
                        else
                            contrib = kvec + bn_value * normal_vec
                            output = output + weight * contrib * inv_r
                        end if
                    end do
                end do
            end do
        end do
    end subroutine integrate_kernel

    subroutine evaluate_surface_state(self, phi, theta, xs, kvec, normal_vec, bn_value, jac_value)
        class(surface_biot_savart_t), intent(in) :: self
        real(dp), intent(in) :: phi
        real(dp), intent(in) :: theta
        real(dp), intent(out) :: xs(3)
        real(dp), intent(out) :: kvec(3)
        real(dp), intent(out) :: normal_vec(3)
        real(dp), intent(out) :: bn_value
        real(dp), intent(out) :: jac_value

        real(dp) :: coords(2)
        integer :: comp

        coords = [phi, theta]
        do comp = 1, 3
            call evaluate_splines_2d(self%x_spline(comp), coords, xs(comp))
            call evaluate_splines_2d(self%k_spline(comp), coords, kvec(comp))
            call evaluate_splines_2d(self%normal_spline(comp), coords, normal_vec(comp))
        end do
        call evaluate_splines_2d(self%bn_spline, coords, bn_value)
        call evaluate_splines_2d(self%jacobian_spline, coords, jac_value)
    end subroutine evaluate_surface_state

    pure function cross_product(a, b) result(c)
        real(dp), intent(in) :: a(3)
        real(dp), intent(in) :: b(3)
        real(dp) :: c(3)
        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end function cross_product

end module tiago_surface_biot_savart
