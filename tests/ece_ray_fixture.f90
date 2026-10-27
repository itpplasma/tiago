module ece_ray_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_ece_ray, only: cold_medium_t
    use tiago_ece_test_medium, only: ece_test_medium
    implicit none
    type, extends(cold_medium_t) :: polynomial_medium_t
        real(dp) :: p(8) = 0.0_dp
    contains
        procedure :: n_params => fixture_n_params
        procedure :: sample
    end type
contains
    pure integer function fixture_n_params(self)
        class(polynomial_medium_t), intent(in) :: self
        fixture_n_params = size(self%p)
    end function
    subroutine sample(self, x, q, dq, d2q, qp, dqp)
        class(polynomial_medium_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: q(4), dq(3, 4), d2q(3, 3, 4)
        real(dp), intent(out) :: qp(:, :), dqp(:, :, :)
        call ece_test_medium(x(1), x(2), x(3), &
                             self%p(1), self%p(2), self%p(3), self%p(4), &
                        self%p(5), self%p(6), self%p(7), self%p(8), q, dq, d2q, qp, dqp)
    end subroutine
end module
