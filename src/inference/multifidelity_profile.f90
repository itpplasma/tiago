module tiago_multifidelity_profile
    !! Physical-temperature AR(1) coefficients in standard-normal coordinates.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_ar1, only: rbf_covariance, covariance_factor
    use fortnum_generated_mf_autoregressive, only: mf_autoregressive
    use fortnum_gaussian_constraints, only: constraint_nullspace
    use fortnum_status, only: fortnum_status_t
    use tiago_calibration, only: calibration_prior_t
    implicit none
    private
    type, public :: ar1_profile_t
        logical :: active = .false., discrepancy = .true., fixed_scale = .false.
        integer :: offset = 0, n = 0, n_shape = 0, n_params = 0
        type(calibration_prior_t) :: calibration
        real(dp), allocatable :: nodes(:), mean_low(:), shape_basis(:, :)
        real(dp) :: log_mean(4) = 0, log_sd(4) = 1
        character(:), allocatable :: hyper_source, constraint_source
    contains
        procedure :: init => ar1_profile_init
        procedure :: coefficients => ar1_profile_coefficients
        procedure :: scale => ar1_profile_scale
        procedure :: scale_index => ar1_profile_scale_index
        procedure :: hyper_index => ar1_profile_hyper_index
    end type
contains
    subroutine ar1_profile_init(self, nodes, mean_low, calibration, log_mean, log_sd, &
            hyper_source, constraint_source, status, &
            protected, discrepancy, fixed_scale)
        class(ar1_profile_t), intent(out) :: self
        real(dp), intent(in) :: nodes(:), mean_low(:), log_mean(4), log_sd(4)
        type(calibration_prior_t), intent(in) :: calibration
        character(*), intent(in) :: hyper_source, constraint_source
        integer, intent(out) :: status
        real(dp), intent(in), optional :: protected(:, :)
        logical, intent(in), optional :: discrepancy, fixed_scale
        real(dp), allocatable :: columns(:, :)
        integer :: n, nc
        status = 1
        n = size(nodes)
        if (n < 3 .or. size(mean_low) /= n) return
        if (.not. all(ieee_is_finite(nodes)) .or. .not. all(ieee_is_finite(mean_low))) return
        if (any(nodes(2:) <= nodes(:n - 1)) .or. any(mean_low <= 0)) return
        if (.not. all(ieee_is_finite(log_mean)) .or. .not. all(ieee_is_finite(log_sd))) return
        if (any(log_sd <= 0)) return
        if (len_trim(hyper_source) == 0 .or. len_trim(constraint_source) == 0) return
        if (calibration%n_params() /= 1) return
        if (.not. allocated(calibration%factor)) return
        if (.not. allocated(calibration%record%record_id)) return
        if (len_trim(calibration%record%record_id) == 0) return
        if (calibration%mean(1) <= 0) return
        if (present(discrepancy)) self%discrepancy = discrepancy
        if (present(fixed_scale)) self%fixed_scale = fixed_scale
        self%n = n
        self%nodes = nodes; self%mean_low = mean_low
        self%calibration = calibration
        self%log_mean = log_mean; self%log_sd = log_sd
        self%hyper_source = hyper_source; self%constraint_source = constraint_source
        nc = 2
        if (present(protected)) then
            if (size(protected, 1) /= n) return
            nc = nc + size(protected, 2)
        end if
        allocate (columns(n, nc))
        columns(:, 1) = 1
        columns(:, 2) = mean_low
        if (present(protected)) columns(:, 3:) = protected
        if (self%discrepancy) then
            call constraint_nullspace(columns, self%shape_basis, status)
            if (status /= 0) return
            self%n_shape = size(self%shape_basis, 2)
            if (self%n_shape == 0) then
                status = 1
                return
            end if
        else
            allocate (self%shape_basis(n, 0))
        end if
        self%n_params = n + self%n_shape + 2
        if (self%discrepancy) self%n_params = self%n_params + 2
        if (.not. self%fixed_scale) self%n_params = self%n_params + 1
        self%active = .true.
        status = 0
    end subroutine

    integer function ar1_profile_hyper_index(self) result(i)
        class(ar1_profile_t), intent(in) :: self
        i = self%n + self%n_shape
    end function

    integer function ar1_profile_scale_index(self) result(i)
        class(ar1_profile_t), intent(in) :: self
        i = 0
        if (.not. self%fixed_scale) i = self%n_params
    end function

    real(dp) function ar1_profile_scale(self, z) result(c)
        class(ar1_profile_t), intent(in) :: self
        real(dp), intent(in) :: z(:)
        c = self%calibration%mean(1)
        if (.not. self%fixed_scale) c = c + self%calibration%factor(1, 1)*z(self%n_params)
    end function

    subroutine ar1_profile_coefficients(self, z, low, high, status, low_jac, high_jac)
        class(ar1_profile_t), intent(in) :: self
        real(dp), intent(in) :: z(:)
        real(dp), intent(out) :: low(:), high(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(out), optional :: low_jac(:, :), high_jac(:, :)
        real(dp) :: hyper(4), c, dr, dl, dd
        real(dp), allocatable :: k(:, :), dk(:, :, :), l(:, :), ld(:, :, :)
        real(dp), allocatable :: delta(:), lj(:, :), dj(:, :), projected(:, :), dprojected(:, :, :)
        integer :: n, ns, h, i, p, nh
        status%code = 1; status%msg = 'AR1 profile: invalid state'
        if (.not. self%active) return
        n = self%n; ns = self%n_shape; h = self%hyper_index()
        if (size(z) /= self%n_params) return
        if (size(low) /= n .or. size(high) /= n) return
        if (present(low_jac)) then
            if (any(shape(low_jac) /= [n, self%n_params])) return
        end if
        if (present(high_jac)) then
            if (any(shape(high_jac) /= [n, self%n_params])) return
        end if
        if (.not. all(ieee_is_finite(z))) return
        nh = 2
        if (self%discrepancy) nh = 4
        hyper = 1
        hyper(:nh) = exp(self%log_mean(:nh) + self%log_sd(:nh)*z(h + 1:h + nh))
        if (.not. all(ieee_is_finite(hyper)) .or. any(hyper <= tiny(1._dp))) return
        c = self%scale(z)
        allocate (k(n, n), dk(n, n, 2), l(n, n), ld(n, n, 2), delta(n), &
            lj(n, self%n_params), dj(n, self%n_params))
        call rbf_covariance(self%nodes, self%nodes, hyper(1), hyper(2), k, dk, status)
        if (status%code /= 0) return
        call covariance_factor(k, dk, l, ld, status)
        if (status%code /= 0) return
        low = self%mean_low + matmul(l, z(:n))
        lj = 0; dj = 0; delta = 0
        lj(:, :n) = l
        do p = 1, 2
            lj(:, h + p) = matmul(ld(:, :, p), z(:n))*self%log_sd(p)
        end do
        if (self%discrepancy) then
            call rbf_covariance(self%nodes, self%nodes, hyper(3), hyper(4), k, dk, status)
            if (status%code /= 0) return
            projected = matmul(transpose(self%shape_basis), matmul(k, self%shape_basis))
            allocate (dprojected(ns, ns, 2))
            do p = 1, 2
                dprojected(:, :, p) = &
                    matmul(transpose(self%shape_basis), matmul(dk(:, :, p), self%shape_basis))
            end do
            deallocate (l, ld)
            allocate (l(ns, ns), ld(ns, ns, 2))
            call covariance_factor(projected, dprojected, l, ld, status)
            if (status%code /= 0) return
            delta = matmul(self%shape_basis, matmul(l, z(n + 1:n + ns)))
            dj(:, n + 1:n + ns) = matmul(self%shape_basis, l)
            do p = 1, 2
                dj(:, h + 2 + p) = matmul(self%shape_basis, &
                    matmul(ld(:, :, p), z(n + 1:n + ns)))*self%log_sd(2 + p)
            end do
        end if
        do i = 1, n
            call mf_autoregressive(c, low(i), delta(i), high(i), dr, dl, dd)
            dj(i, :) = dl*lj(i, :) + dd*dj(i, :)
            if (.not. self%fixed_scale) dj(i, self%n_params) = &
                dj(i, self%n_params) + dr*self%calibration%factor(1, 1)
        end do
        if (present(low_jac)) low_jac = lj
        if (present(high_jac)) high_jac = dj
        if (.not. all(ieee_is_finite(high)) .or. .not. all(ieee_is_finite(low))) then
            status%code = 1; status%msg = 'AR1 profile: nonfinite coefficients'
        end if
    end subroutine
end module tiago_multifidelity_profile
