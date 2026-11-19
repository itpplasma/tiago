module profile_reflectometry_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_profile_reflectometry, only: sight_density_t
    implicit none
    !! ne = p1 s + p2 s^2 + p3 s^3 + p4 (tanh((s - s0)/w) + tanh(s0/w))
    type, extends(sight_density_t) :: polynomial_density_t
        real(dp) :: p(4) = 0.0_dp
        real(dp) :: centre = 0.05_dp, width = 0.02_dp
    contains
        procedure :: n_params
        procedure :: sample
    end type
contains
    pure integer function n_params(self)
        class(polynomial_density_t), intent(in) :: self
        n_params = 4
    end function
    subroutine sample(self, s, ne, dne_ds, ne_p, dne_ds_p)
        class(polynomial_density_t), intent(in) :: self
        real(dp), intent(in) :: s
        real(dp), intent(out) :: ne, dne_ds
        real(dp), intent(out), optional :: ne_p(:), dne_ds_p(:)
        real(dp) :: step, slope
        step = tanh((s - self%centre)/self%width) + tanh(self%centre/self%width)
        slope = (1 - tanh((s - self%centre)/self%width)**2)/self%width
        ne = self%p(1)*s + self%p(2)*s**2 + self%p(3)*s**3 + self%p(4)*step
        dne_ds = self%p(1) + 2*self%p(2)*s + 3*self%p(3)*s**2 + self%p(4)*slope
        if (present(ne_p)) ne_p = [s, s**2, s**3, step]
        if (present(dne_ds_p)) dne_ds_p = [1.0_dp, 2*s, 3*s**2, slope]
    end subroutine
end module profile_reflectometry_fixture

program test_profile_reflectometry
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use profile_reflectometry_fixture, only: polynomial_density_t
    use tiago_profile_reflectometry, only: omode_reflection_t, omode_reflection, &
        invert_omode_group_delay
    use fortnum_status, only: fortnum_status_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), light = 299792458.0_dp
    real(dp), parameter :: electron_mass = 9.1093837015e-31_dp
    real(dp), parameter :: charge = 1.602176634e-19_dp, eps0 = 8.8541878128e-12_dp
    real(dp), parameter :: vacuum = 0.3_dp
    type(polynomial_density_t) :: profile, trial
    type(omode_reflection_t) :: r, rp, rm
    type(fortnum_status_t) :: status
    real(dp) :: f, critical, sc, error, worst, step, fd(3)
    real(dp) :: frequency(180), delay(180), cutoff(180), exact(180)
    integer :: i, j

    ! Linear layer ne = nc(f) s/L (Ginzburg): s_c = L, int N ds = 2L/3 and
    ! int ds/N = 2L, so phase = (4 pi f/c)(d + 2 s_c/3) - pi/2 and
    ! group delay = (2/c)(d + 2 s_c).
    profile%p = 0.0_dp
    profile%p(1) = 4.0e20_dp
    error = 0.0_dp
    do i = 1, 5
        f = 2.0e10_dp + 1.0e10_dp*i
        critical = eps0*electron_mass*(2*pi*f)**2/charge**2
        sc = critical/profile%p(1)
        call omode_reflection(profile, f, vacuum, 1.0_dp, r, status)
        call check(status%code == 0, 'linear layer '//status%msg)
        error = max(error, abs(r%cutoff - sc)/sc)
        error = max(error, abs(r%phase - (4*pi*f/light*(vacuum + 2*sc/3) - pi/2)) &
            /abs(r%phase))
        error = max(error, abs(r%group_delay - 2/light*(vacuum + 2*sc))/r%group_delay)
    end do
    call check(error < 1.0e-9_dp, 'linear layer phase and group delay')
    print '(A,ES12.4)', 'linear layer relative error: ', error

    ! Group delay is the frequency derivative of the phase over 2 pi.
    profile%p = [1.0e20_dp, 2.0e21_dp, -1.0e21_dp, 1.5e19_dp]
    f = 5.0e10_dp
    call omode_reflection(profile, f, vacuum, 1.0_dp, r, status)
    call check(status%code == 0, 'pedestal profile '//status%msg)
    step = 1.0e-5_dp*f
    call omode_reflection(profile, f + step, vacuum, 1.0_dp, rp, status)
    call omode_reflection(profile, f - step, vacuum, 1.0_dp, rm, status)
    error = abs((rp%phase - rm%phase)/(4*pi*step) - r%group_delay)/r%group_delay
    call check(error < 1.0e-7_dp, 'group delay against phase slope')
    print '(A,ES12.4)', 'group delay vs phase slope: ', error

    ! Profile derivatives against central differences.
    worst = 0.0_dp
    do j = 1, 4
        step = 3.0e-3_dp*abs(profile%p(j))
        trial = profile
        trial%p(j) = profile%p(j) + step
        call omode_reflection(trial, f, vacuum, 1.0_dp, rp, status)
        trial%p(j) = profile%p(j) - step
        call omode_reflection(trial, f, vacuum, 1.0_dp, rm, status)
        fd = [(rp%cutoff - rm%cutoff), (rp%phase - rm%phase), &
            (rp%group_delay - rm%group_delay)]/(2*step)
        worst = max(worst, abs(fd(1) - r%cutoff_p(j))/abs(r%cutoff_p(j)))
        worst = max(worst, abs(fd(2) - r%phase_p(j))/abs(r%phase_p(j)))
        worst = max(worst, abs(fd(3) - r%group_delay_p(j))/abs(r%group_delay_p(j)))
    end do
    call check(worst < 2.0e-5_dp, 'profile derivatives')
    print '(A,ES12.4)', 'profile derivative relative error: ', worst

    ! Inversion of synthetic group delays recovers the cutoff positions.
    do i = 1, size(frequency)
        frequency(i) = 2.0e9_dp + 0.25e9_dp*(i - 1)
        call omode_reflection(profile, frequency(i), 0.0_dp, 1.0_dp, r, status)
        call check(status%code == 0, 'inversion data '//status%msg)
        delay(i) = r%group_delay
        exact(i) = r%cutoff
    end do
    call invert_omode_group_delay(frequency, delay, exact(1), cutoff, status)
    call check(status%code == 0, 'inversion '//status%msg)
    error = maxval(abs(cutoff - exact))
    call check(error < 1.0e-4_dp, 'inverted cutoff positions')
    print '(A,ES12.4)', 'inversion max position error [m]: ', error

    call omode_reflection(profile, 4.0e11_dp, vacuum, 1.0_dp, r, status)
    call check(status%code /= 0, 'frequency above the peak density has no cutoff')
contains
    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not. condition) then
            print '(A)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine
end program test_profile_reflectometry
