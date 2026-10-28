module tiago_simple_wall_source
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: iso_c_binding, only: c_int, c_size_t, c_char, c_double, c_null_char
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_wall_synthetic, only: surface_deposition_t, loss_event_block_t, &
        injection_history_t, accumulate_losses, WALL_OK, WALL_ERROR, &
        PARTICLE_COUNT, PARTICLE_RATE
    use tiago_wall_legacy_energy, only: wall_legacy_energy
    implicit none
    private
    public :: accumulate_simple_netcdf
    interface
        integer(c_int) function nc_open(path, mode, id) bind(C)
            import :: c_int, c_char
            character(c_char), intent(in) :: path(*)
            integer(c_int), value :: mode
            integer(c_int), intent(out) :: id
        end function
        integer(c_int) function nc_close(id) bind(C)
            import :: c_int
            integer(c_int), value :: id
        end function
        integer(c_int) function nc_inq_varid(id, name, var) bind(C)
            import :: c_int, c_char
            integer(c_int), value :: id
            character(c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: var
        end function
        integer(c_int) function nc_inq_varndims(id, var, rank) bind(C)
            import :: c_int
            integer(c_int), value :: id, var
            integer(c_int), intent(out) :: rank
        end function
        integer(c_int) function nc_inq_vardimid(id, var, dims) bind(C)
            import :: c_int
            integer(c_int), value :: id, var
            integer(c_int), intent(out) :: dims(*)
        end function
        integer(c_int) function nc_inq_dimlen(id, dim, n) bind(C)
            import :: c_int, c_size_t
            integer(c_int), value :: id, dim
            integer(c_size_t), intent(out) :: n
        end function
        integer(c_int) function nc_get_vara_double(id, var, first, count, values) bind(C)
            import :: c_int, c_size_t, c_double
            integer(c_int), value :: id, var
            integer(c_size_t), intent(in) :: first(*), count(*)
            real(c_double), intent(out) :: values(*)
        end function
        integer(c_int) function nc_get_vara_int(id, var, first, count, values) bind(C)
            import :: c_int, c_size_t
            integer(c_int), value :: id, var
            integer(c_size_t), intent(in) :: first(*), count(*)
            integer(c_int), intent(out) :: values(*)
        end function
        integer(c_int) function nc_inq_attlen(id, var, name, n) bind(C)
            import :: c_int, c_char, c_size_t
            integer(c_int), value :: id, var
            character(c_char), intent(in) :: name(*)
            integer(c_size_t), intent(out) :: n
        end function
        integer(c_int) function nc_get_att_text(id, var, name, value) bind(C)
            import :: c_int, c_char
            integer(c_int), value :: id, var
            character(c_char), intent(in) :: name(*)
            character(c_char), intent(out) :: value(*)
        end function
        integer(c_int) function nc_get_att_double(id, var, name, value) bind(C)
            import :: c_int, c_char, c_double
            integer(c_int), value :: id, var
            character(c_char), intent(in) :: name(*)
            real(c_double), intent(out) :: value(*)
        end function
        integer(c_int) function nc_get_att_int(id, var, name, value) bind(C)
            import :: c_int, c_char
            integer(c_int), value :: id, var
            character(c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: value(*)
        end function
    end interface
contains
    subroutine accumulate_simple_netcdf(path, surface, source_strength, normalization, &
            source_id, source_version, frame, deposited_fraction, status, history, &
            external_weight, block_size)
        character(*), intent(in) :: path, source_id, source_version, frame
        type(surface_deposition_t), intent(inout) :: surface
        real(dp), intent(in) :: source_strength, deposited_fraction(:)
        integer, intent(in) :: normalization
        integer, intent(out) :: status
        type(injection_history_t), optional, intent(in) :: history
        real(dp), optional, intent(in) :: external_weight(:)
        integer, optional, intent(in) :: block_size
        type(surface_deposition_t) :: candidate
        type(loss_event_block_t) :: events
        integer(c_int) :: id, rc, cv, ev, wv, tv, fv, hv, xv, zv, dims(2), rank
        integer(c_int) :: particle_dim, ivalue(1)
        integer(c_size_t) :: total, first(2), count(2), start, n, phase_size, length
        integer :: block, i
        real(dp) :: weight_sum, energy0(1), deriv
        real(dp), allocatable :: work(:)
        character(:), allocatable :: text
        logical :: legacy, weighted
        status = WALL_ERROR
        if (.not. allocated(surface%area_m2)) return
        if (size(deposited_fraction) /= size(surface%area_m2)) return
        if (.not. all(ieee_is_finite(deposited_fraction))) return
        if (any(deposited_fraction < 0) .or. any(deposited_fraction > 1)) return
        if (.not. ieee_is_finite(source_strength)) return
        if (source_strength <= 0) return
        if (normalization /= PARTICLE_COUNT .and. normalization /= PARTICLE_RATE) return
        if (len_trim(source_id) == 0 .or. len_trim(source_version) == 0) return
        if (len_trim(frame) == 0) return
        block = 65536
            if (present(block_size)) block = block_size
            if (block < 1) return
            rc = nc_open(path//c_null_char, 0_c_int, id)
            if (rc /= 0) return
            call text_attribute(-1_c_int, 'results_schema_version', text)
            if (rc == 0) then
                if (text /= '1.1') goto 900
                call text_attribute(-1_c_int, 'cartesian_frame', text)
                if (rc /= 0) goto 900
                if (text /= frame) goto 900
            end if
            rc = nc_inq_varid(id, 'times_lost'//c_null_char, tv)
            if (rc /= 0) goto 900
            rc = nc_inq_varndims(id, tv, rank)
            if (rc /= 0) goto 900
            if (rank /= 1) goto 900
            rc = nc_inq_vardimid(id, tv, dims)
            if (rc /= 0) goto 900
            particle_dim = dims(1)
            rc = nc_inq_dimlen(id, particle_dim, total)
            if (rc /= 0) goto 900
            if (total < 1) goto 900
            call require_units(tv, 's')
            if (rc /= 0) goto 900
            call variable('wall_hit_triangle', fv)
            if (rc /= 0) goto 900
            rc = nc_inq_attlen(id, fv, 'start_index'//c_null_char, length)
            if (rc /= 0) goto 900
            if (length /= 1) goto 900
            rc = nc_get_att_int(id, fv, 'start_index'//c_null_char, ivalue)
            if (rc /= 0) goto 900
            if (ivalue(1) /= 0) goto 900
            rc = nc_inq_attlen(id, -1_c_int, 'stl_facet_count'//c_null_char, length)
            if (rc /= 0) goto 900
            if (length /= 1) goto 900
            rc = nc_get_att_int(id, -1_c_int, 'stl_facet_count'//c_null_char, ivalue)
            if (rc /= 0) goto 900
            if (ivalue(1) /= size(surface%area_m2)) goto 900
            call variable('wall_hit', hv)
            if (rc /= 0) goto 900
            call variable('orbit_exit_code', xv)
            if (rc /= 0) goto 900
            rc = nc_inq_varid(id, 'kinetic_energy_eV'//c_null_char, ev)
            legacy = rc /= 0
            if (.not. legacy) then
                call variable('kinetic_energy_eV', ev)
                if (rc /= 0) goto 900
                call require_units(ev, 'eV')
                if (rc /= 0) goto 900
            else
                rc = nc_inq_attlen(id, -1_c_int, 'particle_energy_eV'//c_null_char, length)
                if (rc /= 0) goto 900
                if (length /= 1) goto 900
                rc = nc_get_att_double(id, -1_c_int, 'particle_energy_eV'//c_null_char, energy0)
                if (rc /= 0) goto 900
                if (.not. ieee_is_finite(energy0(1))) goto 900
                if (energy0(1) < 0) goto 900
                rc = nc_inq_varid(id, 'zend'//c_null_char, zv)
                if (rc /= 0) goto 900
                rc = nc_inq_varndims(id, zv, rank)
                if (rc /= 0) goto 900
                if (rank /= 2) goto 900
                rc = nc_inq_vardimid(id, zv, dims)
                if (rc /= 0) goto 900
                if (dims(1) /= particle_dim) goto 900
                rc = nc_inq_dimlen(id, dims(2), phase_size)
                if (rc /= 0) goto 900
                if (phase_size < 4) goto 900
            end if
            weighted = present(external_weight)
            if (weighted) then
                if (size(external_weight, kind=c_size_t) /= total) goto 900
            else
                call variable('marker_weight', wv)
                if (rc /= 0) goto 900
                call require_units(wv, '1')
                if (rc /= 0) goto 900
                call text_attribute(wv, 'source', text)
                if (rc /= 0) goto 900
                if (text /= 'supplied_relative_weights') goto 900
            end if
            allocate (work(min(int(block, c_size_t), total)))
            weight_sum = 0
            start = 0
            do while (start < total)
                status = WALL_ERROR
                n = min(int(block, c_size_t), total - start)
                first = [start, 0_c_size_t]
                count = [n, 1_c_size_t]
                if (weighted) then
                    work(:n) = external_weight(start + 1:start + n)
                else
                    rc = nc_get_vara_double(id, wv, first, count, work)
                    if (rc /= 0) goto 900
                end if
                if (.not. all(ieee_is_finite(work(:n)))) goto 900
                if (any(work(:n) < 0)) goto 900
                weight_sum = weight_sum + sum(work(:n))
                start = start + n
            end do
            if (.not. ieee_is_finite(weight_sum)) goto 900
            if (weight_sum <= 0) goto 900
            candidate = surface
            start = 0
            do while (start < total)
                status = WALL_ERROR
                n = min(int(block, c_size_t), total - start)
                first = [start, 0_c_size_t]
                count = [n, 1_c_size_t]
                call event_block(int(n))
                rc = nc_get_vara_double(id, tv, first, count, events%time_s)
                if (rc /= 0) goto 900
                rc = nc_get_vara_int(id, fv, first, count, events%facet)
                if (rc /= 0) goto 900
                rc = nc_get_vara_int(id, hv, first, count, events%wall_hit)
                if (rc /= 0) goto 900
                rc = nc_get_vara_int(id, xv, first, count, events%exit_code)
                if (rc /= 0) goto 900
                do i = 1, int(n)
                    if (events%wall_hit(i) /= 1 .or. events%exit_code(i) /= 2) cycle
                    if (events%facet(i) < 0 .or. &
                        events%facet(i) >= size(deposited_fraction)) goto 900
                    events%facet(i) = events%facet(i) + 1
                    events%deposited_fraction(i) = deposited_fraction(events%facet(i))
                end do
                if (weighted) then
                    events%weight = external_weight(start + 1:start + n)
                else
                    rc = nc_get_vara_double(id, wv, first, count, events%weight)
                    if (rc /= 0) goto 900
                end if
                events%weight = events%weight/weight_sum*source_strength
                if (legacy) then
                    first(2) = 3
                    rc = nc_get_vara_double(id, zv, first, count, work)
                    if (rc /= 0) goto 900
                    do i = 1, int(n)
                        call wall_legacy_energy(energy0(1), work(i), &
                            events%kinetic_energy_eV(i), deriv)
                    end do
                else
                    rc = nc_get_vara_double(id, ev, first, count, events%kinetic_energy_eV)
                    if (rc /= 0) goto 900
                end if
                call accumulate_losses(candidate, events, status, history)
                if (status /= WALL_OK) goto 900
                start = start + n
            end do
            surface = candidate
            status = WALL_OK
            900     continue
            cv = nc_close(id)
            if (cv /= 0) status = WALL_ERROR
        contains
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
                if (dims(1) /= particle_dim) rc = -1
            end subroutine
            subroutine text_attribute(var, name, value)
                integer(c_int), intent(in) :: var
                character(*), intent(in) :: name
                character(:), allocatable, intent(out) :: value
                integer(c_size_t) :: length
                rc = nc_inq_attlen(id, var, name//c_null_char, length)
                if (rc /= 0) return
                if (length > 1048576) then
                    rc = -1
                    return
                end if
                allocate (character(int(length)) :: value)
                rc = nc_get_att_text(id, var, name//c_null_char, value)
            end subroutine
            subroutine require_units(var, units)
                integer(c_int), intent(in) :: var
                character(*), intent(in) :: units
                character(:), allocatable :: value
                call text_attribute(var, 'units', value)
                if (rc /= 0) return
                if (value /= units) rc = -1
            end subroutine
            subroutine event_block(n)
                integer, intent(in) :: n
                type(loss_event_block_t) :: empty
                events = empty
                allocate (events%time_s(n), events%kinetic_energy_eV(n), events%weight(n))
                allocate (events%deposited_fraction(n), source=1.0_dp)
                allocate (events%facet(n), events%wall_hit(n), events%exit_code(n))
                events%normalization = normalization
                events%source_id = source_id
                events%version = source_version
                events%frame = frame
                events%position_kind = 'guiding-centre'
            end subroutine
        end subroutine
    end module
