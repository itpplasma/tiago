module tiago_signal_correlation
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR
    implicit none
    private
    type, public :: signal_correlation_t
        private
        real(dp), allocatable :: location(:), scale(:, :)
        real(dp) :: df = 0
    contains
        procedure :: fit
        procedure :: predict
    end type
contains
    subroutine fit(self, samples, mean0, psi0, kappa0, nu0, status)
        class(signal_correlation_t), intent(inout) :: self
        real(dp), intent(in) :: samples(:, :), mean0(:), psi0(:, :), kappa0, nu0
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: average(:), centered(:, :), delta(:), psi(:, :)
        real(dp) :: kappa, count
        integer :: d, n
        if (allocated(self%location)) deallocate(self%location, self%scale)
        self%df = 0
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid correlation prior or samples')
        d = size(mean0)
        n = size(samples, 2)
        if (d < 1 .or. size(samples, 1) /= d) return
        if (size(psi0, 1) /= d .or. size(psi0, 2) /= d) return
        if (.not. all(ieee_is_finite(samples))) return
        if (.not. all(ieee_is_finite(mean0))) return
        if (.not. all(ieee_is_finite(psi0))) return
        if (.not. ieee_is_finite(kappa0) .or. .not. ieee_is_finite(nu0)) return
        if (kappa0 <= 0 .or. nu0 <= d-1) return
        if (maxval(abs(psi0-transpose(psi0))) > &
            32*epsilon(1.0_dp)*maxval(abs(psi0))) return
        if (.not. positive_definite(psi0)) return
        count = real(n, dp)
        kappa = kappa0+count
        self%df = nu0+count-d+1
        average = mean0
        if (n > 0) average = sum(samples, dim=2)/count
        centered = samples-spread(average, 2, n)
        delta = average-mean0
        psi = psi0+matmul(centered, transpose(centered)) &
            + kappa0*count/kappa*spread(delta, 2, d)*spread(delta, 1, d)
        self%scale = psi*((kappa+1)/kappa/self%df)
        self%location = mean0+count/kappa*delta
        if (.not. all(ieee_is_finite(self%scale))) then
            deallocate(self%location, self%scale)
            return
        end if
        call status_set(status, FORTNUM_OK, '')
    end subroutine

    subroutine predict(self, observed, values, missing, location, scale, df, log_observed, status)
        class(signal_correlation_t), intent(in) :: self
        integer, intent(in) :: observed(:)
        real(dp), intent(in) :: values(:)
        integer, allocatable, intent(out) :: missing(:)
        real(dp), allocatable, intent(out) :: location(:), scale(:, :)
        real(dp), intent(out) :: df, log_observed
        type(fortnum_status_t), intent(out) :: status
        logical, allocatable :: used(:)
        real(dp), allocatable :: lower(:, :), residual(:), rhs(:, :), cross(:, :)
        real(dp) :: delta, logdet
        integer :: d, no, nm, i, info
        external :: dpotrf, dpotrs
        df = 0
        log_observed = -huge(1.0_dp)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid observed signal subset')
        if (.not. allocated(self%location)) return
        d = size(self%location)
        no = size(observed)
        if (no > d .or. size(values) /= no) return
        if (any(observed < 1) .or. any(observed > d)) return
        if (.not. all(ieee_is_finite(values))) return
        allocate(used(d))
        used = .false.
        do i = 1, no
            if (used(observed(i))) return
            used(observed(i)) = .true.
        end do
        missing = pack([(i, i=1,d)], .not. used)
        nm = size(missing)
        location = self%location(missing)
        scale = self%scale(missing, missing)
        df = self%df+no
        log_observed = 0
        if (no > 0) then
            lower = self%scale(observed, observed)
            residual = values-self%location(observed)
            allocate(rhs(no, 1+nm))
            rhs(:, 1) = residual
            rhs(:, 2:) = self%scale(observed, missing)
            call dpotrf('L', no, lower, no, info)
            if (info /= 0) return
            logdet = 2*sum([(log(lower(i, i)), i=1,no)])
            call dpotrs('L', no, 1+nm, lower, no, rhs, no, info)
            if (info /= 0) return
            delta = dot_product(residual, rhs(:, 1))
            cross = self%scale(missing, observed)
            location = location+matmul(cross, rhs(:, 1))
            scale = (scale-matmul(cross, rhs(:, 2:)))*((self%df+delta)/df)
            log_observed = log_gamma(df/2)-log_gamma(self%df/2) &
                - (no*log(self%df*acos(-1.0_dp))+logdet)/2 &
                - df/2*log(1+delta/self%df)
        end if
        if (.not. all(ieee_is_finite(location))) return
        if (.not. all(ieee_is_finite(scale))) return
        if (.not. ieee_is_finite(log_observed)) return
        if (nm > 0) then
            if (.not. positive_definite(scale)) return
        end if
        call status_set(status, FORTNUM_OK, '')
    end subroutine

    logical function positive_definite(matrix) result(ok)
        real(dp), intent(in) :: matrix(:, :)
        real(dp), allocatable :: lower(:, :)
        integer :: info
        external :: dpotrf
        lower = matrix
        call dpotrf('L', size(matrix, 1), lower, size(matrix, 1), info)
        ok = info == 0
    end function
end module
