module tiago_netcdf_source
    use, intrinsic :: iso_c_binding, only: &
        c_int, c_size_t, c_char, c_double, c_null_char
    use tiago_signal_source, only: signal_source_t, signal_request_t, &
        signal_window_t, version_matches
    implicit none
    private
    type, extends(signal_source_t), public :: netcdf_source_t
        character(:), allocatable :: path
        character(:), allocatable :: group
    contains
        procedure :: read => file_read
    end type
    interface
        integer(c_int) function nc_open(path, mode, id) bind(C)
            import :: c_int, c_size_t, c_char, c_double
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int), value :: mode
            integer(c_int), intent(out) :: id
        end function
        integer(c_int) function nc_close(id) bind(C)
            import :: c_int, c_size_t, c_char, c_double
            integer(c_int), value :: id
        end function
        integer(c_int) function nc_inq_grp_full_ncid(grpid, grpname, ncgrpid) &
                bind(C)
            import :: c_int, c_size_t, c_char, c_double
            integer(c_int), value :: grpid
            character(kind=c_char), intent(in) :: grpname(*)
            integer(c_int), intent(out) :: ncgrpid
        end function
        integer(c_int) function nc_inq_dimid(id, name, dim) bind(C)
            import :: c_int, c_size_t, c_char, c_double
            integer(c_int), value :: id
            character(kind=c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: dim
        end function
        integer(c_int) function nc_inq_dimlen(id, dim, n) bind(C)
            import :: c_int, c_size_t, c_char, c_double
            integer(c_int), value :: id, dim
            integer(c_size_t), intent(out) :: n
        end function
        integer(c_int) function nc_inq_varid(id, name, var) bind(C)
            import :: c_int, c_size_t, c_char, c_double
            integer(c_int), value :: id
            character(kind=c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: var
        end function
        integer(c_int) function nc_inq_varndims(id, var, n) bind(C)
            import :: c_int, c_size_t, c_char, c_double
            integer(c_int), value :: id, var
            integer(c_int), intent(out) :: n
        end function
        integer(c_int) function nc_inq_vardimid(id, var, dims) bind(C)
            import :: c_int, c_size_t, c_char, c_double
            integer(c_int), value :: id, var
            integer(c_int), intent(out) :: dims(*)
        end function
        integer(c_int) function nc_get_vara_double( &
                id, var, first, count, values) bind(C)
            import :: c_int, c_size_t, c_char, c_double
            integer(c_int), value :: id, var
            integer(c_size_t), intent(in) :: first(*), count(*)
            real(c_double), intent(out) :: values(*)
        end function
        integer(c_int) function nc_inq_attlen(id, var, name, n) bind(C)
            import :: c_int, c_size_t, c_char, c_double
            integer(c_int), value :: id, var
            character(kind=c_char), intent(in) :: name(*)
            integer(c_size_t), intent(out) :: n
        end function
        integer(c_int) function nc_get_att_text(id, var, name, value) bind(C)
            import :: c_int, c_size_t, c_char, c_double
            integer(c_int), value :: id, var
            character(kind=c_char), intent(in) :: name(*)
            character(kind=c_char), intent(out) :: value(*)
        end function
        integer(c_int) function nc_get_att_double(id, var, name, value) bind(C)
            import :: c_int, c_size_t, c_char, c_double
            integer(c_int), value :: id, var
            character(kind=c_char), intent(in) :: name(*)
            real(c_double), intent(out) :: value(*)
        end function
    end interface
contains
    subroutine file_read(self, request, window, status)
        class(netcdf_source_t), intent(inout) :: self
        type(signal_request_t), intent(in) :: request
        type(signal_window_t), intent(out) :: window
        integer, intent(out) :: status
        integer(c_int) :: root, id, dim, tv, sv, vv, cv, rc, dims(2), rank
        integer(c_size_t) :: total, lo, hi, mid, first, after, n
        real(c_double) :: value(1)
        character(:), allocatable :: text, original_checksum

        status = 1
        if (.not. request%valid()) return
        if (.not. allocated(self%path)) return
        rc = nc_open(self%path//c_null_char, 0_c_int, root)
        if (rc /= 0) return
        id = root
        if (allocated(self%group)) then
            if (len(self%group) > 0) then
                rc = nc_inq_grp_full_ncid(root, self%group//c_null_char, id)
                if (rc /= 0) goto 900
            end if
        end if
        call attribute('format', text)
        if (rc /= 0) goto 900
        if (text /= 'tiago-signal-1') goto 900
        call attribute('channel', text)
        if (rc /= 0) goto 900
        if (text /= request%channel) goto 900
        call attribute('source_id', window%provenance%source_id)
        if (rc /= 0) goto 900
        call attribute('version', window%provenance%version)
        if (rc /= 0) goto 900
        if (.not. version_matches(request, window%provenance%version)) goto 900
        call attribute('record_id', window%provenance%record_id)
        if (rc /= 0) goto 900
        call attribute('metadata', window%provenance%metadata)
        if (rc /= 0) goto 900
        call attribute('checksum', original_checksum)
        if (rc /= 0) goto 900
        call attribute('units', window%units)
        if (rc /= 0) goto 900
        call scalar_attribute('valid_from', value)
        if (rc /= 0) goto 900
        window%valid_from = value(1)
        call scalar_attribute('valid_to', value)
        if (rc /= 0) goto 900
        window%valid_to = value(1)
        rc = nc_inq_dimid(id, 'sample'//c_null_char, dim)
        if (rc /= 0) goto 900
        rc = nc_inq_dimlen(id, dim, total)
        if (rc /= 0) goto 900
        call variable('time', tv)
        if (rc /= 0) goto 900
        call variable('samples', sv)
        if (rc /= 0) goto 900
        call variable('variance', vv)
        if (rc /= 0) goto 900
        ! Binary searches keep the time-base working set bounded as well.
        lo = 0
        hi = total
        do while (lo < hi)
            mid = lo + (hi - lo)/2
            rc = nc_get_vara_double(id, tv, [mid], [1_c_size_t], value)
            if (rc /= 0) goto 900
            if (value(1) < request%first) then
                lo = mid + 1
            else
                hi = mid
            end if
        end do
        first = lo
        hi = total
        do while (lo < hi)
            mid = lo + (hi - lo)/2
            rc = nc_get_vara_double(id, tv, [mid], [1_c_size_t], value)
            if (rc /= 0) goto 900
            if (value(1) <= request%last) then
                lo = mid + 1
            else
                hi = mid
            end if
        end do
        after = lo
        n = after - first
        allocate (window%time(n), window%samples(n), window%variance(n))
        if (n > 0) then
            rc = nc_get_vara_double(id, tv, [first], [n], window%time)
            if (rc /= 0) goto 900
            rc = nc_get_vara_double(id, sv, [first], [n], window%samples)
            if (rc /= 0) goto 900
            rc = nc_get_vara_double(id, vv, [first], [n], window%variance)
            if (rc /= 0) goto 900
        end if
        rc = nc_inq_varid(id, 'covariance'//c_null_char, cv)
        if (rc == 0) then
            rc = nc_inq_varndims(id, cv, rank)
            if (rc /= 0) goto 900
            if (rank /= 2) goto 900
            rc = nc_inq_vardimid(id, cv, dims)
            if (rc /= 0) goto 900
            if (any(dims /= dim)) goto 900
            allocate (window%covariance(n, n))
            if (n > 0) then
                rc = nc_get_vara_double(id, cv, [first, first], [n, n], &
                    window%covariance)
                if (rc /= 0) goto 900
            end if
        else if (rc /= -49) then ! NC_ENOTVAR
            goto 900
        end if
        if (.not. window%valid()) goto 900
        call window%seal()
        if (n == total) then
            if (window%provenance%checksum /= original_checksum) goto 900
        end if
        window%provenance%metadata = window%provenance%metadata// &
            new_line('a')//'record_checksum='//original_checksum
        if (any(window%time < request%first)) goto 900
        if (any(window%time > request%last)) goto 900
        window%provenance%request_key = request%key()
        if (.not. window%has_provenance()) goto 900
        status = 0
        900     rc = nc_close(root)
        if (rc /= 0) status = 1
        if (status /= 0) window = signal_window_t()
    contains
        subroutine attribute(name, value)
            character(*), intent(in) :: name
            character(:), allocatable, intent(out) :: value
            integer(c_size_t) :: length
            rc = nc_inq_attlen(id, -1_c_int, name//c_null_char, length)
            if (rc /= 0) return
            allocate (character(length) :: value)
            rc = nc_get_att_text(id, -1_c_int, name//c_null_char, value)
        end subroutine
        subroutine scalar_attribute(name, value)
            character(*), intent(in) :: name
            real(c_double), intent(out) :: value(1)
            integer(c_size_t) :: length
            rc = nc_inq_attlen(id, -1_c_int, name//c_null_char, length)
            if (rc /= 0) return
            if (length /= 1) then
                rc = -1
                return
            end if
            rc = nc_get_att_double(id, -1_c_int, name//c_null_char, value)
        end subroutine
        subroutine variable(name, var)
            character(*), intent(in) :: name
            integer(c_int), intent(out) :: var
            rc = nc_inq_varid(id, name//c_null_char, var)
            if (rc /= 0) return
            rc = nc_inq_varndims(id, var, rank)
            if (rc /= 0) return
            if (rank /= 1) then
                rc = -1
                return
            end if
            rc = nc_inq_vardimid(id, var, dims)
            if (rc /= 0) return
            if (dims(1) /= dim) rc = -1
        end subroutine
    end subroutine
end module tiago_netcdf_source
