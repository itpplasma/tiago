module tiago_robust_signals
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_student_residual, only: student_residual
    use tiago_student_residual_small, only: student_residual_small
    implicit none
    private
    public :: robust_rows, robust_normalization
contains
    subroutine robust_rows(residual, nu, slope)
        real(dp), intent(inout)::residual(:)
        real(dp), intent(in)::nu
        real(dp), intent(out), optional::slope(:)
        real(dp)::value, derivative, r
        integer::i
        if (present(slope)) slope = 1
        if (nu <= 0) return
        !$omp parallel do default(shared) private(i,r,value,derivative) if(size(residual)>1024)
        do i = 1, size(residual)
            r = residual(i)
            if (abs(r) < 1.e-4_dp*sqrt(nu)) then
                call student_residual_small(r, nu, value, derivative)
            else
                call student_residual(abs(r), nu, value, derivative)
                value = sign(value, r)
            end if
            residual(i) = value
            if (present(slope)) slope(i) = derivative
        end do
        !$omp end parallel do
    end subroutine
    real(dp) function robust_normalization(nu) result(value)
        real(dp), intent(in)::nu
        value = 0
        if (nu <= 0) return
        value = log_gamma((nu + 1)/2) - log_gamma(nu/2) - .5_dp*log(nu/2)
    end function
end module
