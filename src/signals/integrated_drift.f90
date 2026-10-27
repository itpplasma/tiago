module tiago_integrated_drift
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_signal_source, only: signal_window_t
    use tiago_signal_reduction, only: reduction_step_t, linear_reduction_t
    use tiago_integrator_signal, only: integrator_signal
    implicit none
    private
    public :: integrated_drift_t

    type, extends(reduction_step_t) :: integrated_drift_t
        private
        type(linear_reduction_t) :: linear
        real(dp), allocatable :: times(:), fit_covariance(:, :)
        integer, allocatable :: keep(:)
    contains
        procedure :: init
        procedure :: baseline_covariance
        procedure :: apply
        procedure :: jvp
        procedure :: vjp
    end type
contains
    subroutine init(self, input, baseline, keep, origin, id, status)
        class(integrated_drift_t), intent(inout) :: self
        type(signal_window_t), intent(in) :: input
        integer, intent(in) :: baseline(:), keep(:)
        real(dp), intent(in) :: origin
        character(*), intent(in) :: id
        integer, intent(out) :: status
        real(dp), allocatable :: c(:, :), cb(:, :), h(:, :), k(:, :), &
                                 normal(:, :), fit(:, :), f(:, :)
        real(dp) :: unused, d_raw, d_gain, one, dt, inverse(2, 2)
        integer :: n, nb, m, i, j, info
        status = 1
        if (allocated(self%times)) deallocate (self%times)
        if (.not. input%valid()) return
        if (.not. input%has_provenance()) return
        if (.not. ieee_is_finite(origin)) return
        n = size(input%samples)
        nb = size(baseline)
        m = size(keep)
        if (nb < 2 .or. m == 0 .or. len_trim(id) == 0) return
        if (any(baseline < 1) .or. any(baseline > n)) return
        if (any(keep < 1) .or. any(keep > n)) return
        do i = 1, nb
            if (any(baseline(:i - 1) == baseline(i))) return
        end do
        do i = 1, m
            if (any(keep(:i - 1) >= keep(i))) return
            if (any(baseline == keep(i))) return
        end do
        c = input%total_covariance()
        cb = c(baseline, baseline)
        allocate (h(nb, 2), f(m, 2))
        do i = 1, nb
            call integrator_signal(0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                        input%time(baseline(i)), origin, unused, d_raw, d_gain, one, dt)
            h(i, :) = [one, dt]
        end do
        do i = 1, m
            call integrator_signal(0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                            input%time(keep(i)), origin, unused, d_raw, d_gain, one, dt)
            f(i, :) = [one, dt]
        end do
        k = h
        call dposv('L', nb, 2, cb, nb, k, nb, info)
        if (info /= 0) return
        normal = matmul(transpose(h), k)
        fit = transpose(k)
        call dposv('L', 2, nb, normal, 2, fit, 2, info)
        if (info /= 0) return
        inverse = 0.0_dp
        inverse(1, 1) = 1.0_dp
        inverse(2, 2) = 1.0_dp
        call dpotrs('L', 2, 2, normal, 2, inverse, 2, info)
        if (info /= 0) return
        self%fit_covariance = inverse
        if (allocated(self%linear%weights)) deallocate (self%linear%weights)
        if (allocated(self%linear%noise_covariance)) then
            deallocate (self%linear%noise_covariance)
        end if
        allocate (self%linear%weights(m, n), self%linear%noise_covariance(m, m))
        self%linear%weights = 0.0_dp
        self%linear%noise_covariance = 0.0_dp
        do i = 1, m
            self%linear%weights(i, keep(i)) = 1.0_dp
        end do
        !$omp parallel do default(shared) private(j) if(nb*m>1024)
        do j = 1, nb
            self%linear%weights(:, baseline(j)) = -matmul(f, fit(:, j))
        end do
        !$omp end parallel do
        self%id = id
        self%linear%id = id
        self%linear%units = input%units
        self%times = input%time
        self%keep = keep
        status = 0
    end subroutine

    function baseline_covariance(self) result(covariance)
        class(integrated_drift_t), intent(in) :: self
        real(dp) :: covariance(2, 2)
        covariance = self%fit_covariance
    end function

    subroutine apply(self, input, output, status)
        class(integrated_drift_t), intent(in) :: self
        type(signal_window_t), intent(in) :: input
        type(signal_window_t), intent(out) :: output
        integer, intent(out) :: status
        status = 1
        if (.not. allocated(self%times)) return
        if (.not. input%valid()) return
        if (size(input%time) /= size(self%times)) return
        if (any(input%time /= self%times)) return
        if (input%units /= self%linear%units) return
        call self%linear%apply(input, output, status)
        if (status /= 0) return
        output%time = input%time(self%keep)
        call output%seal()
    end subroutine

    subroutine jvp(self, input, samples_dot, output_dot)
        class(integrated_drift_t), intent(in) :: self
        type(signal_window_t), intent(in) :: input
        real(dp), intent(in) :: samples_dot(:)
        real(dp), intent(out) :: output_dot(:)
        call self%linear%jvp(input, samples_dot, output_dot)
    end subroutine

    subroutine vjp(self, input, output_bar, samples_bar)
        class(integrated_drift_t), intent(in) :: self
        type(signal_window_t), intent(in) :: input
        real(dp), intent(in) :: output_bar(:)
        real(dp), intent(out) :: samples_bar(:)
        call self%linear%vjp(input, output_bar, samples_bar)
    end subroutine
end module
