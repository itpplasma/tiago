module tiago_iq_phase_reduction
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_signal_source, only: signal_window_t
    use tiago_signal_reduction, only: reduction_step_t
    use tiago_dispersion_atan, only: dispersion_atan
    implicit none
    private
    real(dp), parameter :: twopi = 2.0_dp*acos(-1.0_dp)

    type, extends(reduction_step_t), public :: iq_phase_reduction_t
        real(dp) :: amplitude(2) = 1.0_dp, offset(2) = 0.0_dp
        real(dp) :: anchor = 0.0_dp
        real(dp) :: noise_variance = 0.0_dp
        real(dp) :: min_radius = 1.0e-12_dp
    contains
        procedure :: apply => phase_apply
        procedure :: jvp => phase_jvp
        procedure :: vjp => phase_vjp
    end type
contains
    subroutine rows(self, input, phase, jac)
        class(iq_phase_reduction_t), intent(in) :: self
        type(signal_window_t), intent(in) :: input
        real(dp), intent(out) :: phase(:), jac(:, :)
        real(dp) :: value(1), partial(1, 6), previous
        integer :: k
        !$omp parallel do default(shared) private(k, value, partial) &
        !$omp if(size(phase)>1024)
        do k = 1, size(phase)
            call dispersion_atan(input%samples(2*k - 1), input%samples(2*k), &
                                 self%amplitude(1), self%amplitude(2), self%offset(1), &
                                 self%offset(2), value, partial)
            phase(k) = value(1)
            jac(:, k) = partial(1, :2)
        end do
        !$omp end parallel do
        previous = self%anchor
        do k = 1, size(phase)
            phase(k) = phase(k) + twopi*anint((previous - phase(k))/twopi)
            previous = phase(k)
        end do
    end subroutine

    subroutine phase_apply(self, input, output, status)
        class(iq_phase_reduction_t), intent(in) :: self
        type(signal_window_t), intent(in) :: input
        type(signal_window_t), intent(out) :: output
        integer, intent(out) :: status
        real(dp), allocatable :: jac(:, :)
        real(dp) :: a, b, block(2, 2)
        character(len=256) :: config
        integer :: n, k, j
        status = 1
        if (.not. input%valid()) return
        if (.not. input%has_provenance()) return
        if (.not. allocated(self%id)) return
        if (len_trim(self%id) == 0) return
        n = size(input%samples)/2
        if (n == 0 .or. mod(size(input%samples), 2) /= 0) return
        if (.not. all(ieee_is_finite(self%amplitude))) return
        if (.not. all(ieee_is_finite(self%offset))) return
        if (.not. ieee_is_finite(self%anchor)) return
        if (.not. ieee_is_finite(self%noise_variance)) return
        if (.not. ieee_is_finite(self%min_radius)) return
        if (any(self%amplitude <= 0.0_dp)) return
        if (self%noise_variance < 0.0_dp .or. self%min_radius <= 0.0_dp) return
        do k = 1, n
            if (input%time(2*k - 1) /= input%time(2*k)) return
            a = (input%samples(2*k - 1) - self%offset(1))/self%amplitude(1)
            b = (input%samples(2*k) - self%offset(2))/self%amplitude(2)
            if (hypot(a, b) <= self%min_radius) return
        end do
        allocate (output%samples(n), output%time(n), output%variance(n), jac(2, n))
        call rows(self, input, output%samples, jac)
        output%time = input%time(2::2)
        !$omp parallel do default(shared) private(k) if(n>1024)
        do k = 1, n
            output%variance(k) = dot_product(jac(:, k)**2, &
                                      input%variance(2*k - 1:2*k)) + self%noise_variance
        end do
        !$omp end parallel do
        if (allocated(input%covariance)) then
            allocate (output%covariance(n, n))
            !$omp parallel do default(shared) private(k, j, block) if(n>128)
            do j = 1, n
                do k = 1, j
                    block = input%covariance(2*k - 1:2*k, 2*j - 1:2*j)
                    output%covariance(k, j) = &
                        dot_product(jac(:, k), matmul(block, jac(:, j)))
                    output%covariance(j, k) = output%covariance(k, j)
                end do
            end do
            !$omp end parallel do
        end if
        output%units = 'rad'
        output%valid_from = input%valid_from
        output%valid_to = input%valid_to
        output%provenance = input%provenance
        output%provenance%request_key = input%provenance%request_key// &
                                        '|reduction:'//self%id
        write (config, '(A,8ES24.16E3)') 'iq_phase: ', self%amplitude, self%offset, &
            self%anchor, self%noise_variance, self%min_radius
        output%provenance%metadata = input%provenance%metadata// &
                        new_line('a')//'parent_checksum='//input%provenance%checksum// &
                   new_line('a')//'reduction='//self%id//new_line('a')//trim(config)// &
                         new_line('a')//'noise=conditional first-order atan2 covariance'
        call output%seal()
        status = 0
    end subroutine

    subroutine phase_jvp(self, input, samples_dot, output_dot)
        class(iq_phase_reduction_t), intent(in) :: self
        type(signal_window_t), intent(in) :: input
        real(dp), intent(in) :: samples_dot(:)
        real(dp), intent(out) :: output_dot(:)
        real(dp) :: phase(size(output_dot)), jac(2, size(output_dot))
        integer :: k
        call rows(self, input, phase, jac)
        !$omp parallel do default(shared) private(k) if(size(output_dot)>1024)
        do k = 1, size(output_dot)
            output_dot(k) = dot_product(jac(:, k), samples_dot(2*k - 1:2*k))
        end do
        !$omp end parallel do
    end subroutine

    subroutine phase_vjp(self, input, output_bar, samples_bar)
        class(iq_phase_reduction_t), intent(in) :: self
        type(signal_window_t), intent(in) :: input
        real(dp), intent(in) :: output_bar(:)
        real(dp), intent(out) :: samples_bar(:)
        real(dp) :: phase(size(output_bar)), jac(2, size(output_bar))
        integer :: k
        call rows(self, input, phase, jac)
        !$omp parallel do default(shared) private(k) if(size(output_bar)>1024)
        do k = 1, size(output_bar)
            samples_bar(2*k - 1:2*k) = jac(:, k)*output_bar(k)
        end do
        !$omp end parallel do
    end subroutine
end module
