program test_harmonic_lockin
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use tiago_signal_source, only: signal_window_t
    use tiago_iq_phase_reduction, only: iq_phase_reduction_t
    use tiago_harmonic_lockin, only: harmonic_lockin_t, &
        harmonic_lockin_config_t, harmonic_lockin_output_t
    implicit none
    integer, parameter :: n = 64, periods = 4, total = n*periods
    real(dp), parameter :: pi = acos(-1.0_dp)
    type(harmonic_lockin_config_t) :: c, bad
    type(harmonic_lockin_t) :: whole, chunks
    type(harmonic_lockin_output_t) :: out, piece
    integer(int64) :: time(total), corrupt(1)
    real(dp) :: detector(total), reference(total), expected(2, periods)
    real(dp) :: theta, phase, phi, amplitude, rho(3), got(2, periods)
    integer :: k, j, trial, status, start, finish, count, part, iteration
    rho = [1.2_dp, 2.5_dp, 4.5_dp]
    do trial = 1, 3
        c%cadence_ns = 1000
        c%samples_per_period = n
        c%frequency_hz = 1e9_dp/real(n*c%cadence_ns, dp)
        c%epoch_ns = 1700000000000000000_int64
        c%valid_from_ns = c%epoch_ns
        c%valid_to_ns = c%epoch_ns + 1000000_int64*n*c%cadence_ns
        c%delay_ns = int(3*trial, int64)*c%cadence_ns
        c%reference_phase = 0.37_dp*trial
        c%reference_polarity = (-1)**trial
        c%reference_gain = 1.7_dp
        c%reference_offset = -0.4_dp
        c%reference_tolerance = 1e-11_dp
        c%detector_gain = -2.3_dp
        c%detector_offset = 0.7_dp
        c%modulation_depth = rho(trial)
        c%iq_gain = [1.2_dp, -0.8_dp]
        c%iq_offset = [0.12_dp, -0.23_dp]
        c%minimum_normalizer = 1e-8_dp
        c%alias_tolerance = 1e-11_dp
        c%minimum_radius = 1e-10_dp
        do k = 1, total
            j = (k - 1)/n + 1
            theta = 2*pi*real(modulo(k - 1, n), dp)/n + c%reference_phase
            phase = theta - 2*pi*real(c%delay_ns, dp)/real(n*c%cadence_ns, dp)
            phi = -2.8_dp + 1.7_dp*j
            amplitude = 0.4_dp + 0.6_dp*j
            ! Oracle is the public detector model, independent of mixer kernels.
            detector(k) = c%detector_offset + c%detector_gain*(&
                20.0_dp*j + amplitude*cos(rho(trial)*sin(phase) + phi))
            reference(k) = c%reference_offset + c%reference_gain &
                *c%reference_polarity*sin(theta)
            time(k) = c%epoch_ns + int(k - 1, int64)*c%cadence_ns
            expected(:, j) = ([amplitude*cos(phi), amplitude*sin(phi)] &
                - c%iq_offset)/c%iq_gain
        end do
        call whole%initialize(c, status)
        call require(status == 0, 'initialize')
        call whole%consume(time, detector, reference, out, status)
        call require(status == 0, 'whole period processing')
        call require(size(out%time_ns) == periods, 'output count')
        call require(all(out%valid), 'output validity')
        call require(maxval(abs(out%iq - expected)) < 2e-11_dp, 'analytic IQ')
        call require(all(out%time_ns == time(n::n)), 'int64 endpoint timestamps')
        do part = 1, 3
            call chunks%initialize(c, status)
            call require(status == 0, 'chunk initialize')
            start = 1
            count = 0
            do while (start <= total)
                select case (part)
                case (1)
                    finish = start
                case (2)
                    finish = min(total, start + 16)
                case default
                    finish = min(total, start + modulo(11*start, 79))
                end select
                call chunks%consume(time(start:finish), detector(start:finish), &
                    reference(start:finish), piece, status)
                call require(status == 0, 'chunk processing')
                j = size(piece%time_ns)
                got(:, count + 1:count + j) = piece%iq
                count = count + j
                start = finish + 1
                call require(chunks%buffered_samples() < n, 'bounded state')
            end do
            call require(count == periods, 'chunk output count')
            call require(maxval(abs(got - out%iq)) < 1e-14_dp, 'partition invariance')
        end do
        call chunks%initialize(c, status)
        call chunks%consume(time(:7), detector(:7), reference(:7), piece, status)
        call require(status == 0, 'partial initial period')
        corrupt(1) = time(8) + c%cadence_ns
        call chunks%consume(corrupt, detector(8:8), reference(8:8), piece, status)
        call require(status /= 0, 'gap rejection')
        corrupt(1) = time(7)
        call chunks%consume(corrupt, detector(8:8), reference(8:8), piece, status)
        call require(status /= 0, 'overlap rejection')
        call chunks%consume(time(8:8), detector(8:8), reference(8:8)+1, piece, status)
        call require(status /= 0, 'reference rejection')
        call require(chunks%buffered_samples() == 7, 'rejection preserves state')
        call chunks%consume(time(8:), detector(8:), reference(8:), piece, status)
        call require(status == 0, 'retry rejected chunk')
        call require(maxval(abs(piece%iq - expected)) < 2e-11_dp, 'retry oracle')
        bad = c
        bad%iq_gain = 1
        bad%iq_offset = 0
        call whole%initialize(bad, status)
        call whole%consume(time, detector, reference, out, status)
        call require(status == 0, 'consumer input')
        call consumer(out, c%epoch_ns)
    end do
    bad = c
    bad%modulation_depth = 3.8317059702075123_dp
    call whole%initialize(bad, status)
    call require(status /= 0, 'singular Bessel rejection')
    bad = c
    bad%frequency_hz = 1.01_dp*c%frequency_hz
    call whole%initialize(bad, status)
    call require(status /= 0, 'noncoherent frequency rejection')
    bad = c
    bad%samples_per_period = 8
    bad%frequency_hz = 1e9_dp/real(8*c%cadence_ns, dp)
    call whole%initialize(bad, status)
    call require(status /= 0, 'harmonic alias rejection')
    ! Long streams reuse a single period of caller storage and bounded state.
    c%iq_gain = 1
    c%iq_offset = 0
    call whole%initialize(c, status)
    call require(status == 0, 'long stream initialize')
    do iteration = 1, 1000
        time(:n) = c%epoch_ns + int((iteration - 1)*n, int64)*c%cadence_ns &
            + [(int(k - 1, int64)*c%cadence_ns, k=1,n)]
        call whole%consume(time(:n), detector(:n), reference(:n), piece, status)
        call require(status == 0, 'long stream processing')
        call require(size(piece%time_ns) == 1, 'long stream output count')
        call require(whole%buffered_samples() == 0, 'long stream bounded state')
    end do
    call whole%consume(time(:0), detector(:0), reference(:0), piece, status)
    call require(status == 0, 'empty chunk')
    call require(size(piece%time_ns) == 0, 'empty output')
    call whole%initialize(c, status)
    call whole%consume(time(:1), [ieee_value(0.0_dp, ieee_quiet_nan)], &
        reference(:1), piece, status)
    call require(status /= 0, 'nonfinite detector rejection')
    call require(whole%buffered_samples() == 0, 'nonfinite rejection preserves state')
    time(:n) = c%epoch_ns + [(int(k - 1, int64)*c%cadence_ns, k=1,n)]
    call whole%consume(time(:n), spread(c%detector_offset, 1, n), &
        reference(:n), piece, status)
    call require(status == 0, 'zero signal processing')
    call require(.not. piece%valid(1), 'zero amplitude invalidity')
    print *, 'harmonic lock-in synthetic oracle passed'
contains
    subroutine consumer(output, epoch)
        type(harmonic_lockin_output_t), intent(in) :: output
        integer(int64), intent(in) :: epoch
        type(signal_window_t) :: input, reduced
        type(iq_phase_reduction_t) :: reduction
        integer :: k, status, nout
        nout = size(output%time_ns)
        allocate (input%samples(2*nout), input%time(2*nout), input%variance(2*nout))
        input%samples = reshape(output%iq, [2*nout])
        do k = 1, nout
            input%time(2*k - 1:2*k) = real(output%time_ns(k) - epoch, dp)*1e-9_dp
        end do
        input%variance = 1e-8_dp
        input%units = 'synthetic detector units'
        input%provenance%source_id = 'synthetic cosine detector'
        input%provenance%version = '1'
        input%provenance%request_key = 'harmonic-consumer'
        input%provenance%record_id = 'synthetic'
        input%provenance%metadata = 'calibrated harmonic IQ'
        call input%seal()
        reduction%id = 'harmonic-phase'
        reduction%anchor = -1.1_dp
        call reduction%apply(input, reduced, status)
        call require(status == 0, 'real phase consumer')
        call require(maxval(abs(reduced%samples &
            - [(-2.8_dp + 1.7_dp*k, k=1,nout)])) < 2e-11_dp, &
            'detector to unwrapped phase oracle')
    end subroutine
    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (condition) return
        print *, 'FAIL: ', message
        error stop 1
    end subroutine
end program
