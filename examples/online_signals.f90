program online_signals
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_rest_adapter, only: rest_adapter_t
    use tiago_online_source, only: online_source_t
    use tiago_signal_source, only: cached_source_t, signal_request_t, &
        signal_window_t, SIGNAL_OK
    use tiago_signal_reduction, only: linear_reduction_t
    use tiago_point_probe, only: point_probe_t
    use tiago_profile_diagnostic, only: plasma_sample_t, N_QUANTITIES, Q_TE
    implicit none
    type(rest_adapter_t), target :: adapter
    type(online_source_t), target :: online
    type(linear_reduction_t), target :: mean_step
    type(cached_source_t) :: cache
    type(signal_request_t) :: request
    type(signal_window_t) :: observation
    type(point_probe_t) :: probe
    type(plasma_sample_t) :: sample
    real(dp) :: prediction(1)
    integer :: status
    character(2048) :: config

    call get_command_argument(1, config)
    call adapter%configure(trim(config), status)
    if (status /= SIGNAL_OK) stop 1
    online%adapter => adapter
    online%policy = adapter%policy
    mean_step%id = 'synthetic-mean-v1'
    mean_step%units = 'eV'
    allocate (mean_step%weights(1, 3), mean_step%noise_covariance(1, 1))
    mean_step%weights = 1.0_dp/3.0_dp
    mean_step%noise_covariance = 0.0_dp
    online%reduction => mean_step
    call cache%init(online, 8)
    request%channel = 'ok'
    request%version = 'v1'
    request%first = 0.0_dp
    request%last = 1.0_dp
    call cache%read(request, observation, status)
    if (status /= SIGNAL_OK) stop 2
    call cache%read(request, observation, status)
    if (status /= SIGNAL_OK) stop 3

    ! An independently specified synthetic provider predicts 5 eV here.
    allocate (sample%q(N_QUANTITIES, 1))
    sample%q = 0.0_dp
    sample%q(Q_TE, 1) = 5.0_dp
    probe%quantity = [Q_TE]
    probe%x0 = reshape([1.0_dp, 0.0_dp, 0.0_dp], [3, 1])
    call probe%predict(sample, [real(dp) ::], prediction)
    print '(A,F8.4)', 'Observed window mean [eV]: ', observation%samples(1)
    print '(A,F8.4)', 'Point-probe prediction [eV]: ', prediction(1)
    if (abs(prediction(1) - observation%samples(1)) > 1e-13_dp) stop 4
end program
