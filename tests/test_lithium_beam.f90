program test_lithium_beam
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_lithium_beam, only: beam_atomic_data_t, beam_emission_t, &
        lithium_beam_emission
    use fortnum_status, only: fortnum_status_t
    implicit none
    type(beam_atomic_data_t) :: two, three
    type(beam_emission_t) :: r, rp, rm
    type(fortnum_status_t) :: status
    real(dp), allocatable :: z(:), ne(:), te(:), trial(:)
    real(dp) :: errors(2), worst, step, fd, a, b, c, d, lam1, lam2, length, total
    real(dp) :: n2_exact
    integer :: n, j, k, level

    ! Two levels in a uniform plasma: dN/dz = M N / v has the closed form
    ! of a 2x2 matrix exponential; Crank-Nicolson converges at second order.
    call two_level(two)
    length = 0.2_dp
    a = -(2.0e-14_dp + 1.0e-14_dp)*3.0e19_dp
    b = 5.0e-15_dp*3.0e19_dp + 3.7e7_dp
    c = 2.0e-14_dp*3.0e19_dp
    d = -(5.0e-15_dp + 4.0e-14_dp)*3.0e19_dp - 3.7e7_dp
    lam1 = 0.5_dp*(a + d) + sqrt(0.25_dp*(a - d)**2 + b*c)
    lam2 = 0.5_dp*(a + d) - sqrt(0.25_dp*(a - d)**2 + b*c)
    ! N2(t) for N(0) = (1, 0): c (exp(lam1 t) - exp(lam2 t))/(lam1 - lam2).
    n2_exact = c*(exp(lam1*length/two%speed) - exp(lam2*length/two%speed)) &
        /(lam1 - lam2)
    do level = 1, 2
        n = 1 + 100*2**(level - 1)
        call uniform(n, length, 3.0e19_dp, 100.0_dp)
        call lithium_beam_emission(two, z, ne, te, 0.0_dp, r, status)
        call check(status%code == 0, 'two-level beam '//status%msg)
        errors(level) = abs(r%population(2, n) - n2_exact)/n2_exact
    end do
    call check(errors(2) < 1.0e-5_dp, 'two-level population')
    call check(abs(errors(1)/errors(2) - 4) < 0.1_dp, 'second-order convergence')
    print '(A,2ES12.4)', 'two-level relative error at h, h/2: ', errors

    ! Without ionisation the total population is a linear invariant.
    call three_level(three)
    three%electron_loss = 0.0_dp
    three%ion_loss = 0.0_dp
    call profile(80)
    call lithium_beam_emission(three, z, ne, te, 1.0_dp, r, status)
    call check(status%code == 0, 'conserving beam '//status%msg)
    total = maxval(abs(sum(r%population, dim=1) - 1))
    call check(total < 1.0e-13_dp, 'population conservation without ionisation')

    ! Derivatives in ne and te at every node against central differences.
    call three_level(three)
    call lithium_beam_emission(three, z, ne, te, 1.0_dp, r, status)
    call check(status%code == 0, 'attenuated beam '//status%msg)
    call check(sum(r%population(:, size(z))) < 0.99_dp, 'beam attenuation')
    worst = 0.0_dp
    do j = 1, size(z), 7
        trial = ne
        step = 1.0e-4_dp*ne(j)
        trial(j) = ne(j) + step
        call lithium_beam_emission(three, z, trial, te, 1.0_dp, rp, status)
        trial(j) = ne(j) - step
        call lithium_beam_emission(three, z, trial, te, 1.0_dp, rm, status)
        do k = 1, size(z)
            fd = (rp%emission(k) - rm%emission(k))/(2*step)
            worst = max(worst, abs(fd - r%emission_ne(k, j)) &
                /maxval(abs(r%emission_ne(:, j))))
        end do
        trial = te
        step = 1.0e-4_dp*te(j)
        trial(j) = te(j) + step
        call lithium_beam_emission(three, z, ne, trial, 1.0_dp, rp, status)
        trial(j) = te(j) - step
        call lithium_beam_emission(three, z, ne, trial, 1.0_dp, rm, status)
        do k = 1, size(z)
            fd = (rp%emission(k) - rm%emission(k))/(2*step)
            worst = max(worst, abs(fd - r%emission_te(k, j)) &
                /maxval(abs(r%emission_te(:, j))))
        end do
    end do
    call check(worst < 1.0e-6_dp, 'emission derivatives')
    print '(A,ES12.4)', 'emission derivative relative error: ', worst

    te(10) = 1.0e5_dp
    call lithium_beam_emission(three, z, ne, te, 1.0_dp, r, status)
    call check(status%code /= 0, 'temperature outside the table is rejected')
contains
    subroutine uniform(n, l, density, temperature)
        integer, intent(in) :: n
        real(dp), intent(in) :: l, density, temperature
        integer :: i
        z = [(l*(i - 1)/(n - 1), i=1, n)]
        ne = [(density, i=1, n)]
        te = [(temperature, i=1, n)]
    end subroutine

    subroutine profile(n)
        integer, intent(in) :: n
        integer :: i
        z = [(0.15_dp*(i - 1)/(n - 1), i=1, n)]
        ne = 1.0e18_dp + 4.0e19_dp*(1 + tanh((z - 0.08_dp)/0.02_dp))
        te = 20.0_dp + 400.0_dp*(1 + tanh((z - 0.09_dp)/0.02_dp))
    end subroutine

    subroutine two_level(atomic)
        type(beam_atomic_data_t), intent(out) :: atomic
        atomic%temperature = [1.0_dp, 1.0e4_dp]
        allocate (atomic%electron_rate(2, 2, 2), atomic%electron_loss(2, 2))
        atomic%electron_rate = 0.0_dp
        atomic%electron_rate(2, 1, :) = 2.0e-14_dp
        atomic%electron_rate(1, 2, :) = 5.0e-15_dp
        atomic%electron_loss(1, :) = 1.0e-14_dp
        atomic%electron_loss(2, :) = 4.0e-14_dp
        allocate (atomic%ion_rate(2, 2), atomic%ion_loss(2), atomic%einstein(2, 2))
        atomic%ion_rate = 0.0_dp
        atomic%ion_loss = 0.0_dp
        atomic%einstein = 0.0_dp
        atomic%einstein(1, 2) = 3.7e7_dp
        atomic%speed = 1.2e6_dp
    end subroutine

    subroutine three_level(atomic)
        !! Synthetic rates of the size of lithium 2s, 2p, 3s at 50 keV.
        type(beam_atomic_data_t), intent(out) :: atomic
        integer :: i
        atomic%temperature = [5.0_dp, 50.0_dp, 500.0_dp, 5000.0_dp]
        allocate (atomic%electron_rate(3, 3, 4), atomic%electron_loss(3, 4))
        atomic%electron_rate = 0.0_dp
        do i = 1, 4
            atomic%electron_rate(2, 1, i) = 1.0e-14_dp*(0.3_dp + 0.2_dp*i)
            atomic%electron_rate(3, 1, i) = 1.0e-15_dp*i
            atomic%electron_rate(1, 2, i) = 3.0e-15_dp*(0.5_dp + 0.1_dp*i)
            atomic%electron_rate(3, 2, i) = 2.0e-14_dp*(0.2_dp + 0.2_dp*i)
            atomic%electron_rate(2, 3, i) = 3.0e-14_dp*(0.4_dp + 0.1_dp*i)
            atomic%electron_loss(:, i) = [5.0e-15_dp, 2.0e-14_dp, 5.0e-14_dp]*i
        end do
        allocate (atomic%ion_rate(3, 3), atomic%ion_loss(3), atomic%einstein(3, 3))
        atomic%ion_rate = 0.0_dp
        atomic%ion_rate(2, 1) = 4.0e-15_dp
        atomic%ion_rate(3, 2) = 6.0e-15_dp
        atomic%ion_loss = [1.0e-15_dp, 5.0e-15_dp, 1.0e-14_dp]
        atomic%einstein = 0.0_dp
        atomic%einstein(1, 2) = 3.7e7_dp
        atomic%einstein(2, 3) = 3.5e7_dp
        atomic%speed = 1.2e6_dp
    end subroutine

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not. condition) then
            print '(A)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine
end program test_lithium_beam
