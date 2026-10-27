module tiago_calibration
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_signal_source, only: signal_source_t, signal_request_t, &
                                   signal_window_t, signal_provenance_t
    implicit none
    private
    type, public :: calibration_prior_t
        real(dp), allocatable :: mean(:), covariance(:, :), factor(:, :)
        type(signal_provenance_t) :: record
        real(dp) :: valid_from = -huge(1.0_dp), valid_to = huge(1.0_dp)
    contains
        procedure :: init => prior_init
        procedure :: load => prior_load
        procedure :: n_params => prior_count
        procedure :: physical => prior_physical
        procedure :: jvp => prior_jvp
        procedure :: vjp => prior_vjp
    end type
contains
    subroutine prior_init(self, mean, covariance, status)
        class(calibration_prior_t), intent(inout) :: self
        real(dp), intent(in) :: mean(:), covariance(:, :)
        integer, intent(out) :: status
        integer :: n, i
        n = size(mean)
        self%record = signal_provenance_t()
        self%valid_from = -huge(1.0_dp)
        self%valid_to = huge(1.0_dp)
        status = 1
        if (allocated(self%mean)) deallocate (self%mean)
        if (allocated(self%factor)) deallocate (self%factor)
        if (allocated(self%covariance)) deallocate (self%covariance)
        if (n == 0) return
        if (any(shape(covariance) /= [n, n])) return
        if (.not. all(ieee_is_finite(mean))) return
        if (.not. all(ieee_is_finite(covariance))) return
        if (any(covariance /= transpose(covariance))) return
        self%factor = covariance
        call dpotrf('L', n, self%factor, n, status)
        if (status /= 0) then
            deallocate (self%factor)
            return
        end if
        do i = 1, n
            self%factor(i, i + 1:) = 0.0_dp
        end do
        self%mean = mean
        self%covariance = covariance
    end subroutine

    subroutine prior_load(self, source, request, status)
        class(calibration_prior_t), intent(inout) :: self
        class(signal_source_t), intent(inout) :: source
        type(signal_request_t), intent(in) :: request
        integer, intent(out) :: status
        type(signal_window_t) :: data
        real(dp), allocatable :: covariance(:, :)
        status = 1
        if (.not. request%valid()) return
        call source%read(request, data, status)
        if (status /= 0) return
        status = 1
        if (.not. data%valid()) return
        if (.not. data%has_provenance()) return
        if (request%first < data%valid_from .or. request%last > data%valid_to) return
        if (.not. allocated(data%provenance%record_id)) return
        if (len(data%provenance%record_id) == 0) return
        if (size(data%time) == 0) return
        if (any(data%time /= data%time(1))) return
        covariance = data%total_covariance()
        call self%init(data%samples, covariance, status)
        if (status /= 0) return
        self%record = data%provenance
        self%valid_from = data%valid_from
        self%valid_to = data%valid_to
    end subroutine

    integer function prior_count(self)
        class(calibration_prior_t), intent(in) :: self
        prior_count = 0
        if (allocated(self%mean)) prior_count = size(self%mean)
    end function

    function prior_physical(self, z) result(c)
        class(calibration_prior_t), intent(in) :: self
        real(dp), intent(in) :: z(:)
        real(dp) :: c(size(z))
        c = self%mean + matmul(self%factor, z)
    end function

    function prior_jvp(self, z_dot) result(c_dot)
        class(calibration_prior_t), intent(in) :: self
        real(dp), intent(in) :: z_dot(:)
        real(dp) :: c_dot(size(z_dot))
        c_dot = matmul(self%factor, z_dot)
    end function

    function prior_vjp(self, c_bar) result(z_bar)
        class(calibration_prior_t), intent(in) :: self
        real(dp), intent(in) :: c_bar(:)
        real(dp) :: z_bar(size(c_bar))
        z_bar = matmul(transpose(self%factor), c_bar)
    end function
end module tiago_calibration
