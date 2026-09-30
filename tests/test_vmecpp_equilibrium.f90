program test_vmecpp_equilibrium
    !! VMEC++ through Tiago's C adapter, without Python.
    !! Usage: test_vmecpp_equilibrium <VMEC input (INDATA or JSON)> <workdir> [adjoint-only]
    !! (adjoint-only: checks 1-2 and the Jacobian timing, no finite differences,
    !! for timing large equilibria)
    !!  1. Tiago's edge field y equals VMEC++'s wout (1.5 b(ns) - 0.5 b(ns-1)).
    !!  2. The edge VJP matches directional finite differences of the edge map.
    !!  3. The full Jacobian cot (dy/dx) (edge VJP + VMEC++ adjoint + profile
    !!     chain + PHIEDGE scaling) matches central differences of re-solved
    !!     equilibria for every parameter kind.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_equilibrium, only: equilibrium_t
    implicit none

    character(len=16), parameter :: names(7) = [character(len=16) :: 'PRES_SCALE', 'CURTOR', &
        'AM(1)', 'AC(1)', 'RBC(1,1)', 'ZBS(0,1)', 'PHIEDGE']
    type(equilibrium_t) :: eq
    character(len=512) :: input, workdir, mode
    real(dp), allocatable :: x0(:), x(:), y0(:), yp(:), ym(:), cot(:, :), jac(:, :), fd(:)
    real(dp), allocatable :: ref(:), b(:, :), coef_bar(:, :, :, :), dcoef(:, :, :, :)
    real(dp), allocatable :: iota_bar(:), current_bar(:), ybar(:), cp(:, :, :, :), cm(:, :, :, :)
    real(dp) :: h, err, t0, t1, eps
    real(dp), allocatable :: xp(:), yhot(:)
    integer :: failures, k, ns, nb, ng, mnq, mnmax
    logical :: ok

    failures = 0
    call get_command_argument(1, input)
    call get_command_argument(2, workdir)
    call get_command_argument(3, mode)
    call execute_command_line('mkdir -p "'//trim(workdir)//'"')
    call eq%init(trim(input), trim(workdir), names)
    if (trim(mode) /= 'adjoint-only') then  ! else FTOL and NITER of the input
        call eq%vmec%set_input('ftol', 1.0e-14_dp)
        call eq%vmec%set_input('niter', 20000.0_dp)
    end if
    x0 = eq%values()
    call cpu_time(t0)
    call eq%solve(x0, ok)
    call cpu_time(t1)
    if (.not. ok) error stop 'VMEC++ did not converge'
    print '(A,F7.2,A)', 'solve: ', t1 - t0, ' s'

    ! 1. edge field vs VMEC++'s own wout
    ns = eq%ns
    mnq = eq%vmec%get_int('mnmax_nyq')
    mnmax = eq%vmec%get_int('mnmax')
    nb = eq%edge%mnmax_nyq
    ng = eq%edge%mnmax
    call check_equal('Nyquist mode count', real(nb, dp), real(mnq, dp), 0.0_dp)
    call check_equal('mode count', real(ng, dp), real(mnmax, dp), 0.0_dp)
    call check_equal('xn_nyq table', sum(abs(eq%edge%xn_nyq - eq%vmec%wout('xn_nyq', mnq))), 0.0_dp, 0.0_dp)
    b = reshape(eq%vmec%wout('bsupumnc', ns * mnq), [mnq, ns])
    ref = 1.5_dp * b(:, ns) - 0.5_dp * b(:, ns - 1)
    call check_vector('edge B^u vs VMEC++ wout', eq%y(1:nb), ref, 1.0e-9_dp)
    b = reshape(eq%vmec%wout('bsupvmnc', ns * mnq), [mnq, ns])
    ref = 1.5_dp * b(:, ns) - 0.5_dp * b(:, ns - 1)
    call check_vector('edge B^v vs VMEC++ wout', eq%y(nb + 1:2 * nb), ref, 1.0e-9_dp)
    b = reshape(eq%vmec%wout('rmnc', ns * mnmax), [mnmax, ns])
    call check_vector('edge rmnc vs VMEC++ wout', eq%y(2 * nb + 1:2 * nb + ng), b(:, ns), 1.0e-12_dp)
    b = reshape(eq%vmec%wout('zmns', ns * mnmax), [mnmax, ns])
    call check_vector('edge zmns vs VMEC++ wout', eq%y(2 * nb + ng + 1:), b(:, ns), 1.0e-12_dp)

    ! 2. edge VJP vs finite differences along a random geometry direction
    allocate(ybar(size(eq%y)))
    allocate(coef_bar, dcoef, mold=eq%coef)
    allocate(iota_bar(ns - 1), current_bar(ns - 1), yp(size(eq%y)), ym(size(eq%y)))
    block
        integer :: n_seed
        call random_seed(size=n_seed)
        call random_seed(put=[(17 * k + 5, k = 1, n_seed)])
    end block
    call random_number(ybar)
    ybar = ybar - 0.5_dp
    call random_number(dcoef)
    dcoef = (dcoef - 0.5_dp) * max(abs(eq%coef), 1.0e-3_dp * maxval(abs(eq%coef)))
    dcoef(:, :, :ns - 3, :) = 0.0_dp
    call eq%edge%vjp(eq%coef, eq%phip_full, eq%phip_half, eq%iota_half, eq%current_half, &
        ybar, coef_bar, iota_bar, current_bar)
    eps = 1.0e-6_dp
    cp = eq%coef + eps * dcoef
    cm = eq%coef - eps * dcoef
    call eq%edge%forward(cp, eq%phip_full, eq%phip_half, eq%iota_half, eq%current_half, yp)
    call eq%edge%forward(cm, eq%phip_full, eq%phip_half, eq%iota_half, eq%current_half, ym)
    ! relative to the size of the terms: the random sum cancels strongly
    call check_close('edge VJP vs FD (geometry)', sum(coef_bar * dcoef), &
        dot_product(ybar, yp - ym) / (2.0_dp * eps), 1.0e-7_dp * sum(abs(coef_bar * dcoef)))
    block
        real(dp), allocatable :: cur(:)
        cur = eq%current_half
        cur(ns - 1) = cur(ns - 1) * (1.0_dp + eps)
        call eq%edge%forward(eq%coef, eq%phip_full, eq%phip_half, eq%iota_half, cur, yp)
        cur(ns - 1) = eq%current_half(ns - 1) * (1.0_dp - eps)
        call eq%edge%forward(eq%coef, eq%phip_full, eq%phip_half, eq%iota_half, cur, ym)
        if (eq%ncurr == 1) call check_equal('edge VJP vs FD (current)', &
            current_bar(ns - 1) * eq%current_half(ns - 1), dot_product(ybar, yp - ym) / (2.0_dp * eps), &
            1.0e-7_dp)
    end block

    ! 3. a hot-restarted solve (from the last equilibrium) agrees with a cold one
    if (trim(mode) /= 'adjoint-only') then
        xp = x0 * (1.0_dp + 1.0e-3_dp)
        call eq%solve(xp, ok)
        yhot = eq%y
        call eq%vmec%set_input('hot_restart', 0.0_dp)
        call eq%solve(xp, ok)
        call check_vector('hot restart vs cold solve', yhot, eq%y, 1.0e-6_dp)
        call eq%solve(x0, ok)
    end if

    ! 4. full Jacobian vs central differences of re-solved equilibria. Cold
    !    solves: at FTOL = 1e-14 a solve is accurate to ~1e-7, smoothly in x
    !    along the same (cold) path, but not across hot restarts from
    !    different start states, which would swamp differences of 1e-5 steps.
    y0 = eq%y
    allocate(cot(4, size(y0)))
    call random_number(cot)
    cot = (cot - 0.5_dp) / spread(max(abs(y0), 1.0e-3_dp * maxval(abs(y0))), 1, 4)
    call cpu_time(t0)
    call eq%jacobian(x0, cot, jac)
    call cpu_time(t1)
    print '(A,F7.2,A)', 'Jacobian (4 cotangents, incl. factorization): ', t1 - t0, ' s'
    if (trim(mode) == 'adjoint-only') then
        call cpu_time(t0)
        call eq%jacobian(x0, cot, jac)
        call cpu_time(t1)
        print '(A,F7.2,A)', 'Jacobian (4 cotangents, factorization cached): ', t1 - t0, ' s'
        call eq%vmec%destroy()
        if (failures > 0) error stop 1
        stop
    end if
    do k = 1, size(names)
        h = 1.0e-5_dp * max(abs(x0(k)), 1.0e-2_dp)
        x = x0
        x(k) = x0(k) + h
        call eq%solve(x, ok)
        yp = eq%y
        x(k) = x0(k) - h
        call eq%solve(x, ok)
        ym = eq%y
        fd = matmul(cot, yp - ym) / (2.0_dp * h)
        err = maxval(abs(jac(:, k) - fd)) / maxval(abs(fd))
        if (err > 2.0e-3_dp) then
            write(error_unit, '(A,A12,A,ES10.3)') 'FAIL Jacobian ', trim(names(k)), &
                ': max |adjoint - FD| / max |FD| = ', err
            failures = failures + 1
        else
            print '(A,A12,A,ES10.3)', 'ok   Jacobian ', trim(names(k)), &
                ': max |adjoint - FD| / max |FD| = ', err
        end if
    end do
    call eq%vmec%destroy()
    if (failures > 0) error stop 1
    print '(A)', 'test_vmecpp_equilibrium passed'

contains

    subroutine check_vector(name, actual, expected, rtol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:), expected(:), rtol
        real(dp) :: e
        e = maxval(abs(actual - expected)) / max(maxval(abs(expected)), tiny(1.0_dp))
        if (e > rtol) then
            write(error_unit, '(A,ES10.3)') 'FAIL '//name//': max rel diff ', e
            failures = failures + 1
        else
            print '(A,ES10.3)', 'ok   '//name//': max rel diff ', e
        end if
    end subroutine check_vector

    subroutine check_close(name, actual, expected, atol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual, expected, atol
        if (abs(actual - expected) > atol) then
            write(error_unit, '(A,3ES24.16)') 'FAIL '//name//': ', actual, expected, atol
            failures = failures + 1
        else
            print '(A,2ES24.16)', 'ok   '//name//': ', actual, expected
        end if
    end subroutine check_close

    subroutine check_equal(name, actual, expected, rtol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual, expected, rtol
        if (abs(actual - expected) > rtol * max(abs(expected), tiny(1.0_dp))) then
            write(error_unit, '(A,2ES24.16)') 'FAIL '//name//': ', actual, expected
            failures = failures + 1
        else
            print '(A,2ES24.16)', 'ok   '//name//': ', actual, expected
        end if
    end subroutine check_equal
end program test_vmecpp_equilibrium
