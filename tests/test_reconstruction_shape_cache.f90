program test_reconstruction_shape_cache
    !! A current change at fixed boundary must refresh the shape cotangent.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_reconstruction, only: rec, KIND_BPROBE, jacobian_at, model_signals
    implicit none

    character(len=512) :: input, workdir
    character(len=16), parameter :: names(2) = [character(len=16) :: &
        'CURTOR', 'RBC(1,1)']
    real(dp), allocatable :: x(:), shifted(:), jac(:, :), initial(:, :)
    real(dp), allocatable :: flux(:), seg(:), plus(:), minus(:), fd(:)
    real(dp) :: step, error, norm
    integer :: k
    logical :: ok

    call get_command_argument(1, input)
    call get_command_argument(2, workdir)
    call execute_command_line('mkdir -p "'//trim(workdir)//'"')
    call rec%plasma%init('')
    call rec%eq%init(trim(input), trim(workdir), names, conservative=.true.)
    call rec%eq%vmec%set_input('ftol', 1.0e-14_dp)
    call rec%eq%vmec%set_input('niter', 20000.0_dp)
    rec%n = 2
    rec%names = names
    rec%is_coil = [.false., .false.]
    rec%eq_index = [1, 2]
    allocate(rec%extcur(0), rec%probes(3))
    rec%n_meas = 3
    rec%n_diag_probes = 3
    rec%meas_kind = [KIND_BPROBE, KIND_BPROBE, KIND_BPROBE]
    rec%meas_index = [1, 2, 3]
    rec%meas_sigma = [1.0_dp, 1.0_dp, 1.0_dp]
    rec%probes(1)%position = [2.5_dp, 0.0_dp, 0.2_dp]
    rec%probes(2)%position = [2.2_dp, 0.3_dp, -0.3_dp]
    rec%probes(3)%position = [0.3_dp, 2.3_dp, 0.3_dp]
    do k = 1, 3
        rec%probes(k)%normal = [0.3_dp, 0.5_dp, 0.8_dp]
        rec%probes(k)%normal = rec%probes(k)%normal / norm2(rec%probes(k)%normal)
    end do

    x = rec%eq%values()
    call jacobian_at(x, initial)
    x(1) = 0.5_dp * x(1)
    call jacobian_at(x, jac)
    call rec%eq%vmec%set_input('hot_restart', 0.0_dp)
    step = 1.0e-4_dp * max(abs(x(2)), 0.1_dp)
    shifted = x
    shifted(2) = x(2) + step
    call model_signals(shifted, flux, seg, plus, ok)
    if (.not. ok) error stop 'positive FD equilibrium failed'
    shifted(2) = x(2) - step
    call model_signals(shifted, flux, seg, minus, ok)
    if (.not. ok) error stop 'negative FD equilibrium failed'
    fd = (plus - minus) / (2.0_dp * step)
    if (.not. all(ieee_is_finite(fd))) error stop 'nonfinite FD derivative'
    if (.not. all(ieee_is_finite(jac))) error stop 'nonfinite adjoint derivative'
    norm = maxval(abs(fd))
    if (norm < 1.0e-8_dp) error stop 'shape derivative oracle is too small'
    error = maxval(abs(jac(:, 2) - fd)) / norm
    print '(A,ES12.4)', &
        'shape derivative after fixed-boundary current step, relative FD error: ', error
    if (error > 1.0e-3_dp) error stop 'stale or incorrect shape cotangent'
end program test_reconstruction_shape_cache
