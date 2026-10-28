module joint_signal_io
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf
    use tiago_signal_source, only: signal_request_t, signal_window_t
    use tiago_netcdf_source, only: netcdf_source_t
    use tiago_rest_adapter, only: rest_adapter_t
    use tiago_online_source, only: online_source_t
    use joint_recorded_fixture, only: require
    implicit none
contains
    subroutine signal_roundtrip(directory, name, units, values, variance)
        character(*), intent(in)::directory, name, units
        real(dp), intent(inout)::values(:)
        real(dp), intent(in)::variance(:)
        type(signal_window_t)::raw, offline, online
        type(signal_request_t)::request
        type(netcdf_source_t)::file
        type(rest_adapter_t), target::adapter
        type(online_source_t)::source
        integer::i, code, id, dim, tv, sv, vv, u
        request%channel = name; request%version = 'v1'
        request%first = 0; request%last = real(size(values) - 1, dp)
        raw%time = [(real(i - 1, dp), i=1, size(values))]
        raw%samples = values; raw%variance = variance; raw%units = units
        raw%valid_from = request%first; raw%valid_to = request%last
        raw%provenance%source_id = 'synthetic-recorded'
        raw%provenance%version = 'v1'; raw%provenance%record_id = 'synthetic'
        raw%provenance%metadata = 'public analytic synthetic signal packet'
        raw%provenance%request_key = request%key()
        call raw%seal()
        file%path = directory//'/'//name//'.nc'
        call nc(nf90_create(file%path, nf90_clobber, id))
        call nc(nf90_def_dim(id, 'sample', size(values), dim))
        call nc(nf90_def_var(id, 'time', nf90_double, [dim], tv))
        call nc(nf90_def_var(id, 'samples', nf90_double, [dim], sv))
        call nc(nf90_def_var(id, 'variance', nf90_double, [dim], vv))
        call nc(nf90_put_att(id, nf90_global, 'format', 'tiago-signal-1'))
        call nc(nf90_put_att(id, nf90_global, 'channel', name))
        call nc(nf90_put_att(id, nf90_global, 'source_id', raw%provenance%source_id))
        call nc(nf90_put_att(id, nf90_global, 'version', 'v1'))
        call nc(nf90_put_att(id, nf90_global, 'record_id', raw%provenance%record_id))
        call nc(nf90_put_att(id, nf90_global, 'metadata', raw%provenance%metadata))
        call nc(nf90_put_att(id, nf90_global, 'checksum', raw%provenance%checksum))
        call nc(nf90_put_att(id, nf90_global, 'units', units))
        call nc(nf90_put_att(id, nf90_global, 'valid_from', raw%valid_from))
        call nc(nf90_put_att(id, nf90_global, 'valid_to', raw%valid_to))
        call nc(nf90_enddef(id))
        call nc(nf90_put_var(id, tv, raw%time))
        call nc(nf90_put_var(id, sv, values))
        call nc(nf90_put_var(id, vv, variance))
        call nc(nf90_close(id))
        open (newunit=u, file=directory//'/'//name//'.json', status='replace')
        write (u, '(A)', advance='no') '{"version":"v1","data":{"time":'
        call vector(u, raw%time)
        write (u, '(A)', advance='no') ',"values":'
        call vector(u, values)
        write (u, '(A)', advance='no') ',"variance":'
        call vector(u, variance)
        write (u, '(A)') '}}'
        close (u)
        call file%read(request, offline, code)
        call require(code == 0, 'joint offline source '//name)
        call adapter%configure(directory//'/adapter.json', code)
        call require(code == 0, 'joint adapter configuration')
        source%adapter => adapter; source%policy = adapter%policy
        source%policy%max_window = 2
        call source%read(request, online, code)
        call require(code == 0, 'joint online source '//name)
        call require(size(online%samples) == size(values), 'joint online sample count')
        call require(all(offline%samples == online%samples), 'offline/HTTP signal parity')
        call require(all(offline%variance == online%variance), 'offline/HTTP noise parity')
        call require(all(offline%time == online%time), 'offline/HTTP clock parity')
        values = online%samples
        print '(A,ES12.4)', 'offline/HTTP maximum signal difference '//name//'=', &
            maxval(abs(offline%samples - online%samples))
    end subroutine
    subroutine vector(unit, x)
        integer, intent(in)::unit
        real(dp), intent(in)::x(:)
        integer::i
        write (unit, '(A)', advance='no') '['
        do i = 1, size(x)
            if (i > 1) write (unit, '(A)', advance='no') ','
            write (unit, '(ES25.17E3)', advance='no') x(i)
        end do
        write (unit, '(A)', advance='no') ']'
    end subroutine
    subroutine nc(code)
        integer, intent(in)::code
        if (code /= nf90_noerr) print *, nf90_strerror(code)
        call require(code == nf90_noerr, 'joint netCDF fixture')
    end subroutine
end module
