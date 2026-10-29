module tiago_ece_state_medium
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_ece_ray, only: cold_medium_t
    use tiago_real_space_state, only: real_space_state_t
    use tiago_kinetic_profiles, only: kinetic_profiles_t
    use tiago_state_sampler, only: sample_medium
    implicit none
    private
    type, extends(cold_medium_t), public :: state_medium_t
        class(real_space_state_t), allocatable :: state
        type(kinetic_profiles_t) :: profiles
        real(dp), allocatable :: z(:)
    contains
        procedure :: n_params => state_medium_count
        procedure :: sample => state_medium_sample
    end type
contains
    pure integer function state_medium_count(self) result(n)
        class(state_medium_t), intent(in) :: self
        n = self%state%n_params() + self%profiles%n_latent
    end function
    subroutine state_medium_sample(self, x, q, dq, d2q, qp, dqp)
        class(state_medium_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: q(4), dq(3, 4), d2q(3, 3, 4)
        real(dp), intent(out), optional :: qp(:, :), dqp(:, :, :)
        integer :: status
        call sample_medium(self%state, self%profiles, self%z, x, q, dq, d2q, qp, dqp, &
                           status)
    end subroutine
end module
