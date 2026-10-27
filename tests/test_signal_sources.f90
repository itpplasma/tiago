module synthetic_source
    use tiago_signal_source, only: signal_source_t, signal_request_t, signal_window_t
    implicit none
    type, extends(signal_source_t) :: counting_source_t
        type(signal_window_t) :: data
        integer :: calls = 0
        logical :: fail = .false.
    contains
        procedure :: read => counting_read
    end type
contains
    subroutine counting_read(self, request, window, status)
        class(counting_source_t), intent(inout) :: self
        type(signal_request_t), intent(in) :: request
        type(signal_window_t), intent(out) :: window
        integer, intent(out) :: status
        self%calls = self%calls + 1
        status = 1
        if (self%fail) return
        window = self%data
        window%provenance%request_key = request%key()
        status = 0
    end subroutine
end module

module synthetic_adapter
    use, intrinsic :: iso_fortran_env, only: int8, dp => real64
    use tiago_signal_source, only: signal_adapter_t, signal_request_t, &
                backend_request_t, signal_window_t, SIGNAL_OK, SIGNAL_ERROR, SIGNAL_MORE
    implicit none
    type, extends(signal_adapter_t) :: fixture_adapter_t
        character(32) :: version = ''
        integer(int8), allocatable :: buffered(:)
    contains
        procedure :: configure => fixture_configure
        procedure :: map_request => fixture_map
        procedure :: decode => fixture_decode
    end type
contains
    subroutine fixture_configure(self, path, status)
        class(fixture_adapter_t), intent(inout) :: self
        character(*), intent(in) :: path
        integer, intent(out) :: status
        integer :: unit
        open (newunit=unit, file=path, status='old', iostat=status)
        if (status /= 0) return
        read (unit, '(A)', iostat=status) self%version
        close (unit)
        if (allocated(self%buffered)) deallocate (self%buffered)
        allocate (self%buffered(0))
    end subroutine
    subroutine fixture_map(self, request, backend, status)
        class(fixture_adapter_t), intent(in) :: self
        type(signal_request_t), intent(in) :: request
        type(backend_request_t), intent(out) :: backend
        integer, intent(out) :: status
        status = SIGNAL_ERROR
        if (.not. request%valid()) return
        backend%endpoint = 'synthetic:'//request%channel
        backend%method = 'GET'
        backend%body = trim(self%version)
        backend%request_key = request%key()
        status = SIGNAL_OK
    end subroutine
    subroutine fixture_decode(self, request, bytes, final, window, status)
        class(fixture_adapter_t), intent(inout) :: self
        type(signal_request_t), intent(in) :: request
        integer(int8), intent(in) :: bytes(:)
        logical, intent(in) :: final
        type(signal_window_t), intent(out) :: window
        integer, intent(out) :: status
        integer :: i
        self%buffered = [self%buffered, bytes]
        status = SIGNAL_MORE
        if (.not. final) return
        window%samples = real(self%buffered, dp)
        window%time = [(real(i - 1, dp), i=1, size(self%buffered))]
        allocate (window%variance(size(self%buffered)))
        window%variance = 1.0_dp
        window%units = 'V'
        window%provenance%source_id = 'synthetic-adapter'
        window%provenance%version = trim(self%version)
        window%provenance%record_id = 'synthetic-chunks'
        window%provenance%metadata = 'one signed byte per synthetic sample'
        window%provenance%request_key = request%key()
        call window%seal()
        deallocate (self%buffered)
        allocate (self%buffered(0))
        status = SIGNAL_OK
    end subroutine
end module

program test_signal_sources
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit, int8
    use netcdf, only: nf90_create, nf90_close, nf90_clobber, nf90_double, &
                  nf90_def_dim, nf90_def_var, nf90_put_att, nf90_enddef, nf90_put_var, &
                      nf90_global
    use tiago_signal_source, only: signal_request_t, signal_window_t, &
                                   cached_source_t, backend_request_t, SIGNAL_MORE
    use tiago_netcdf_source, only: netcdf_source_t
    use tiago_signal_reduction, only: linear_reduction_t, reduced_source_t
    use tiago_calibration, only: calibration_prior_t
    use synthetic_source, only: counting_source_t
    use synthetic_adapter, only: fixture_adapter_t
    implicit none
    type(signal_request_t) :: request, other
    type(fixture_adapter_t) :: adapter
    type(backend_request_t) :: backend
    type(signal_window_t) :: raw, got, reduced, plus, minus, yp, ym
    type(netcdf_source_t), target :: file
    type(counting_source_t), target :: source
    type(cached_source_t), target :: cache
    type(cached_source_t) :: product_cache
    type(linear_reduction_t), target :: step, second
    type(reduced_source_t), target :: pipeline
    type(reduced_source_t) :: chain
    type(calibration_prior_t) :: prior
    real(dp) :: dot(4), tangent(2), bar(2), adjoint(4), expected(2, 2)
    real(dp) :: fd_error, adj_error, h, c(2), cdot(2), zbar(2)
    real(dp) :: covariance_error, prior_error
    integer :: status, fails, i, unit

    fails = 0
    request%channel = 'synthetic-trace'
    request%version = 'v1'
    request%first = 0.0_dp
    request%last = 3.0_dp
    open (newunit=unit, file='synthetic_adapter.cfg', status='replace')
    write (unit, '(A)') 'v1'
    close (unit)
    call adapter%configure('synthetic_adapter.cfg', status)
    call check(status == 0, 'adapter configured by file')
    call adapter%map_request(request, backend, status)
    call check(status == 0, 'adapter maps logical request')
    call check(backend%request_key == request%key(), 'backend request key')
    call adapter%decode(request, [1_int8, 2_int8], .false., got, status)
    call check(status == SIGNAL_MORE, 'incremental adapter needs more bytes')
    call adapter%decode(request, [3_int8], .true., got, status)
    call check(status == 0, 'incremental adapter completes window')
    call check(all(got%samples == [1.0_dp, 2.0_dp, 3.0_dp]), 'decoded chunks')
    call check(got%has_provenance(), 'adapter provenance')
    open (newunit=unit, file='synthetic_adapter.cfg')
    close (unit, status='delete')
    raw%time = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    raw%samples = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    raw%variance = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    raw%units = 'V'
    raw%provenance%source_id = 'synthetic-source'
    raw%provenance%version = 'v1'
    raw%provenance%record_id = 'synthetic-record'
    raw%provenance%metadata = 'analytic synthetic ramp'
    raw%provenance%request_key = request%key()
    call raw%seal()
    call write_file('synthetic_signal.nc', raw)
    file%path = 'synthetic_signal.nc'
    call file%read(request, got, status)
    call check(status == 0, 'netCDF round trip')
    if (status /= 0) error stop 'cannot continue without round trip'
    call check(all(got%samples == raw%samples), 'samples unchanged')
    call check(all(got%time == raw%time), 'time unchanged')
    call check(got%provenance%checksum == raw%provenance%checksum, 'checksum')
    call check(got%provenance%source_id == raw%provenance%source_id, 'source id')
    call check(got%units == 'V', 'units')
    other = request
    other%first = 0.5_dp
    other%last = 2.0_dp
    call file%read(other, got, status)
    call check(status == 0, 'partial window')
    call check(all(got%samples == [2.0_dp, 3.0_dp]), 'closed window endpoints')
    call check(got%provenance%request_key == other%key(), 'request provenance')
    other%first = 8.0_dp
    other%last = 9.0_dp
    call file%read(other, got, status)
    call check(status == 0, 'empty window succeeds')
    call check(size(got%samples) == 0, 'empty allocation')
    other = request
    other%version = 'v2'
    call file%read(other, got, status)
    call check(status /= 0, 'pinned version mismatch rejected')
    other = request
    other%channel = 'unknown'
    call file%read(other, got, status)
    call check(status /= 0, 'unknown channel rejected')
    other = request
    other%first = 4.0_dp
    call file%read(other, got, status)
    call check(status /= 0, 'reversed window rejected')
    deallocate (other%version)
    other%first = 0.0_dp
    call file%read(other, got, status)
    call check(status == 0, 'unpinned resolves file version')
    raw%samples(1) = 99.0_dp
    call write_file('synthetic_bad.nc', raw)
    file%path = 'synthetic_bad.nc'
    call file%read(request, got, status)
    call check(status /= 0, 'bad checksum rejected')
    raw%samples(1) = 1.0_dp
    call write_file('synthetic_bad.nc', raw, .true.)
    call file%read(request, got, status)
    call check(status /= 0, 'nonscalar validity attribute rejected')

    source%data = raw
    call cache%init(source, 2)
    call cache%read(request, got, status)
    call check(status == 0, 'cache miss fetches')
    source%fail = .true.
    call cache%read(request, got, status)
    call check(status == 0 .and. source%calls == 1, 'read cache before fetch')
    got%samples(1) = -1.0_dp
    call cache%read(request, got, status)
    call check(got%samples(1) == 1.0_dp, 'cache returns independent copy')
    other = request
    other%last = 2.0_dp
    call cache%read(other, got, status)
    call check(status /= 0 .and. source%calls == 2, 'window is part of key')
    other = request
    other%version = 'v2'
    call cache%read(other, got, status)
    call check(status /= 0 .and. source%calls == 3, 'version is part of key')
    source%fail = .false.
    call cache%read(other, got, status)
    call check(status /= 0, 'cache rejects wrong upstream version')
    call cache%clear()
    call cache%read(request, got, status)
    call check(status == 0 .and. source%calls == 5, 'clear invalidates cache')

    ! Linear propagation of a Gaussian covariance: W C W^T + N.
    ! JCGM 100:2008, section 5.2 (correlated input quantities).
    step%id = 'synthetic-integral-v1'
    step%units = 'V s'
    step%weights = reshape([1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                            0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], [2, 4])
    step%noise_covariance = reshape([0.5_dp, 0.1_dp, 0.1_dp, 0.5_dp], [2, 2])
    allocate (raw%covariance(4, 4))
    raw%covariance = 0.0_dp
    raw%covariance(1, 3) = 0.25_dp
    raw%covariance(3, 1) = 0.25_dp
    call raw%seal()
    call step%apply(raw, reduced, status)
    call check(status == 0, 'linear reduction')
    call check(all(reduced%samples == [3.0_dp, 7.0_dp]), 'analytic pulse sum')
    expected = reshape([3.5_dp, 0.35_dp, 0.35_dp, 7.5_dp], [2, 2])
    covariance_error = maxval(abs(reduced%covariance - expected))
    call check(covariance_error < 1.0e-14_dp, &
               'correlated covariance and step noise')
    call check(index(reduced%provenance%metadata, raw%provenance%checksum) > 0, &
               'reduction preserves parent checksum')
    dot = [0.3_dp, -0.2_dp, 0.7_dp, -0.4_dp]
    h = 1.0e-5_dp
    plus = raw
    minus = raw
    plus%samples = raw%samples + h*dot
    minus%samples = raw%samples - h*dot
    call step%apply(plus, yp, status)
    call step%apply(minus, ym, status)
    call step%jvp(raw, dot, tangent)
    fd_error = maxval(abs((yp%samples - ym%samples)/(2*h) - tangent))
    bar = [0.7_dp, -0.3_dp]
    call step%vjp(raw, bar, adjoint)
    adj_error = abs(dot_product(bar, tangent) - dot_product(adjoint, dot))
    call check(fd_error < 1.0e-9_dp, 'reduction central differences')
    call check(adj_error < 1.0e-14_dp, 'reduction adjoint identity')
    step%noise_covariance(1, 1) = -1.0_dp
    call step%apply(raw, got, status)
    call check(status /= 0, 'indefinite reduction noise rejected')
    step%noise_covariance(1, 1) = 0.5_dp
    source%data = raw
    call cache%clear()
    pipeline%upstream => cache
    pipeline%step => step
    second%id = 'synthetic-second-v1'
    second%units = 'V s'
    second%weights = reshape([1.0_dp, -1.0_dp], [1, 2])
    second%noise_covariance = reshape([0.2_dp], [1, 1])
    chain%upstream => pipeline
    chain%step => second
    call chain%read(request, got, status)
    call check(status == 0, 'chained source reductions')
    call check(abs(got%samples(1) + 4.0_dp) < 1.0e-14_dp, 'chain observable')
    call check(abs(got%covariance(1, 1) - 10.5_dp) < 1.0e-14_dp, 'chain covariance')
    call product_cache%init(chain, 1)
    call product_cache%read(request, got, status)
    source%fail = .true.
    call product_cache%read(request, got, status)
    call check(status == 0, 'cache stores reduced products')

    raw%time = [1.0_dp, 1.0_dp]
    raw%samples = [2.0_dp, 0.3_dp]
    raw%variance = [0.04_dp, 0.09_dp]
    raw%covariance = reshape([0.0_dp, 0.01_dp, 0.01_dp, 0.0_dp], [2, 2])
    raw%units = 'mixed; gain=1; offset=V'
    raw%valid_from = 0.0_dp
    raw%valid_to = 3.0_dp
    call raw%seal()
    call write_file('synthetic_calibration.nc', raw)
    file%path = 'synthetic_calibration.nc'
    call prior%load(file, request, status)
    call check(status == 0, 'calibration loaded through source')
    call check(prior%record%record_id == 'synthetic-record', 'calibration record id')
    prior_error = maxval(abs(matmul(prior%factor, transpose(prior%factor)) - &
                             prior%covariance))
    call check(prior_error < 1.0e-14_dp, 'prior reproduces correlated covariance')
    c = prior%physical([0.0_dp, 0.0_dp])
    call check(all(c == raw%samples), 'zero latent is prior mean')
    cdot = prior%jvp(bar)
    zbar = prior%vjp(c)
    call check(abs(dot_product(c, cdot) - dot_product(zbar, bar)) < 1.0e-14_dp, &
               'calibration latent transpose')
    other = request
    other%last = 4.0_dp
    call prior%load(file, other, status)
    call check(status /= 0, 'calibration validity window enforced')
    call prior%init([1.0_dp, 2.0_dp], &
                    reshape([1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp], [2, 2]), status)
    call check(status /= 0, 'indefinite prior rejected')
    do i = 1, 3
        select case (i)
        case (1)
            open (newunit=unit, file='synthetic_signal.nc')
        case (2)
            open (newunit=unit, file='synthetic_bad.nc')
        case (3)
            open (newunit=unit, file='synthetic_calibration.nc')
        end select
        close (unit, status='delete')
    end do
    print '(A,ES12.4)', 'correlated covariance max absolute error: ', covariance_error
    print '(A,ES12.4)', 'calibration covariance max absolute error: ', prior_error
    print '(A,ES12.4)', 'reduction FD max absolute error: ', fd_error
    print '(A,ES12.4)', 'reduction adjoint absolute error: ', adj_error
    if (fails > 0) error stop 1
    print '(A)', 'signal sources, covariance reductions and calibration: PASS'
contains
    subroutine check(ok, text)
        logical, intent(in) :: ok
        character(*), intent(in) :: text
        if (ok) return
        write (error_unit, '(A)') 'FAIL: '//text
        fails = fails + 1
    end subroutine

    subroutine nc(code)
        integer, intent(in) :: code
        if (code /= 0) then
            write (error_unit, *) 'netCDF fixture error', code
            error stop 1
        end if
    end subroutine

    subroutine write_file(path, data, malformed)
        character(*), intent(in) :: path
        type(signal_window_t), intent(in) :: data
        logical, intent(in), optional :: malformed
        integer :: id, dim, tv, sv, vv, cv
        call nc(nf90_create(path, nf90_clobber, id))
        call nc(nf90_def_dim(id, 'sample', size(data%time), dim))
        call nc(nf90_def_var(id, 'time', nf90_double, [dim], tv))
        call nc(nf90_def_var(id, 'samples', nf90_double, [dim], sv))
        call nc(nf90_def_var(id, 'variance', nf90_double, [dim], vv))
        if (allocated(data%covariance)) &
            call nc(nf90_def_var(id, 'covariance', nf90_double, [dim, dim], cv))
        call nc(nf90_put_att(id, nf90_global, 'format', 'tiago-signal-1'))
        call nc(nf90_put_att(id, nf90_global, 'channel', request%channel))
        call nc(nf90_put_att(id, nf90_global, 'source_id', data%provenance%source_id))
        call nc(nf90_put_att(id, nf90_global, 'version', data%provenance%version))
        call nc(nf90_put_att(id, nf90_global, 'record_id', data%provenance%record_id))
        call nc(nf90_put_att(id, nf90_global, 'metadata', data%provenance%metadata))
        call nc(nf90_put_att(id, nf90_global, 'checksum', data%provenance%checksum))
        call nc(nf90_put_att(id, nf90_global, 'units', data%units))
        if (present(malformed)) then
            call nc(nf90_put_att(id, nf90_global, 'valid_from', &
                                 [data%valid_from, data%valid_from]))
        else
            call nc(nf90_put_att(id, nf90_global, 'valid_from', data%valid_from))
        end if
        call nc(nf90_put_att(id, nf90_global, 'valid_to', data%valid_to))
        call nc(nf90_enddef(id))
        call nc(nf90_put_var(id, tv, data%time))
        call nc(nf90_put_var(id, sv, data%samples))
        call nc(nf90_put_var(id, vv, data%variance))
        if (allocated(data%covariance)) call nc(nf90_put_var(id, cv, data%covariance))
        call nc(nf90_close(id))
    end subroutine
end program
