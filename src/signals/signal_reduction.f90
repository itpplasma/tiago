module tiago_signal_reduction
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_signal_source, only: signal_source_t, signal_request_t, signal_window_t
    implicit none
    private

    type, abstract, public :: reduction_step_t
        character(:), allocatable :: id
    contains
        procedure(apply_i), deferred :: apply
        procedure(jvp_i), deferred :: jvp
        procedure(vjp_i), deferred :: vjp
    end type
    abstract interface
        subroutine apply_i(self, input, output, status)
            import :: reduction_step_t, signal_window_t
            class(reduction_step_t), intent(in) :: self
            type(signal_window_t), intent(in) :: input
            type(signal_window_t), intent(out) :: output
            integer, intent(out) :: status
        end subroutine
        subroutine jvp_i(self, input, samples_dot, output_dot)
            import :: reduction_step_t, signal_window_t, dp
            class(reduction_step_t), intent(in) :: self
            type(signal_window_t), intent(in) :: input
            real(dp), intent(in) :: samples_dot(:)
            real(dp), intent(out) :: output_dot(:)
        end subroutine
        subroutine vjp_i(self, input, output_bar, samples_bar)
            import :: reduction_step_t, signal_window_t, dp
            class(reduction_step_t), intent(in) :: self
            type(signal_window_t), intent(in) :: input
            real(dp), intent(in) :: output_bar(:)
            real(dp), intent(out) :: samples_bar(:)
        end subroutine
    end interface

    type, extends(reduction_step_t), public :: linear_reduction_t
        real(dp), allocatable :: weights(:, :)
        real(dp), allocatable :: noise_covariance(:, :)
        character(:), allocatable :: units
    contains
        procedure :: apply => linear_apply
        procedure :: jvp => linear_jvp
        procedure :: vjp => linear_vjp
    end type

    type, extends(signal_source_t), public :: reduced_source_t
        class(signal_source_t), pointer :: upstream => null()
        class(reduction_step_t), pointer :: step => null()
    contains
        procedure :: read => reduced_read
    end type
contains
    subroutine reduced_read(self, request, window, status)
        class(reduced_source_t), intent(inout) :: self
        type(signal_request_t), intent(in) :: request
        type(signal_window_t), intent(out) :: window
        integer, intent(out) :: status
        type(signal_window_t) :: raw
        status = 1
        if (.not. request%valid()) return
        if (.not. associated(self%upstream)) return
        if (.not. associated(self%step)) return
        call self%upstream%read(request, raw, status)
        if (status /= 0) return
        call self%step%apply(raw, window, status)
    end subroutine

    subroutine linear_apply(self, input, output, status)
        class(linear_reduction_t), intent(in) :: self
        type(signal_window_t), intent(in) :: input
        type(signal_window_t), intent(out) :: output
        integer, intent(out) :: status
        real(dp), allocatable :: scaled(:, :), eig_matrix(:, :), eig(:), work(:)
        integer :: n, m, i, info
        status = 1
        if (.not. input%valid()) return
        if (.not. input%has_provenance()) return
        if (.not. allocated(self%id)) return
        if (.not. allocated(self%weights)) return
        if (.not. allocated(self%noise_covariance)) return
        if (.not. allocated(self%units)) return
        n = size(input%samples)
        m = size(self%weights, 1)
        if (n == 0 .or. m == 0) return
        if (size(self%weights, 2) /= n) return
        if (any(shape(self%noise_covariance) /= [m, m])) return
        if (.not. all(ieee_is_finite(self%weights))) return
        if (.not. all(ieee_is_finite(self%noise_covariance))) return
        if (any(self%noise_covariance /= transpose(self%noise_covariance))) return
        eig_matrix = self%noise_covariance
        allocate (eig(m), work(max(1, 3*m - 1)))
        call dsyev('N', 'L', m, eig_matrix, m, eig, work, size(work), info)
        if (info /= 0) return
        if (any(eig < 0.0_dp)) return
        output%samples = matmul(self%weights, input%samples)
        allocate (output%time(m), output%variance(m), scaled(m, n))
        output%time = input%time(n)
        output%variance = 0.0_dp
        !$omp parallel do schedule(static) if(n>1024)
        do i = 1, n
            scaled(:, i) = self%weights(:, i)*input%variance(i)
        end do
        !$omp end parallel do
        output%covariance = matmul(scaled, transpose(self%weights)) + &
                            self%noise_covariance
        if (allocated(input%covariance)) then
            output%covariance = output%covariance + matmul(self%weights, &
                                      matmul(input%covariance, transpose(self%weights)))
        end if
        ! Roundoff from matrix products need not preserve exact symmetry.
        do i = 1, m
            output%covariance(i, :i - 1) = output%covariance(:i - 1, i)
        end do
        output%units = self%units
        output%valid_from = input%valid_from
        output%valid_to = input%valid_to
        output%provenance = input%provenance
        output%provenance%request_key = input%provenance%request_key// &
                                        '|reduction:'//self%id
        output%provenance%metadata = input%provenance%metadata// &
                            new_line('a')//'parent_checksum='//input%provenance%checksum
        output%provenance%metadata = output%provenance%metadata// &
                                     new_line('a')//'reduction='//self%id
        call output%seal()
        status = 0
    end subroutine

    subroutine linear_jvp(self, input, samples_dot, output_dot)
        class(linear_reduction_t), intent(in) :: self
        type(signal_window_t), intent(in) :: input
        real(dp), intent(in) :: samples_dot(:)
        real(dp), intent(out) :: output_dot(:)
        output_dot = matmul(self%weights, samples_dot)
    end subroutine

    subroutine linear_vjp(self, input, output_bar, samples_bar)
        class(linear_reduction_t), intent(in) :: self
        type(signal_window_t), intent(in) :: input
        real(dp), intent(in) :: output_bar(:)
        real(dp), intent(out) :: samples_bar(:)
        samples_bar = matmul(transpose(self%weights), output_bar)
    end subroutine
end module tiago_signal_reduction
