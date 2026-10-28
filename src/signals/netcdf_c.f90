module tiago_netcdf_c
    !! netCDF-4 groups through the netCDF C library. The Fortran `netcdf`
    !! module of the build (fortio) has no groups; see tiago_netcdf_source.
    !! Array variables are passed in C order, so a (n, 3) file variable is a
    !! Fortran (3, n) array.
    use, intrinsic :: iso_c_binding, only: c_int, c_size_t, c_char, c_double, c_null_char
    implicit none
    private

    integer(c_int), parameter, public :: NC_NOERR = 0, NC_NOWRITE = 0, NC_CLOBBER = 0, &
        NC_NETCDF4 = int(z'1000', c_int), NC_GLOBAL = -1, NC_CHAR = 2, NC_INT = 4, NC_DOUBLE = 6

    interface
        integer(c_int) function nc_open(path, mode, id) bind(C)
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int), value :: mode
            integer(c_int), intent(out) :: id
        end function
        integer(c_int) function nc_create(path, mode, id) bind(C)
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int), value :: mode
            integer(c_int), intent(out) :: id
        end function
        integer(c_int) function nc_close(id) bind(C)
            import :: c_int
            integer(c_int), value :: id
        end function
        integer(c_int) function nc_inq_ncid(id, name, group) bind(C)
            import :: c_int, c_char
            integer(c_int), value :: id
            character(kind=c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: group
        end function
        integer(c_int) function nc_inq_dimid(id, name, dim) bind(C)
            import :: c_int, c_char
            integer(c_int), value :: id
            character(kind=c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: dim
        end function
        integer(c_int) function nc_def_grp(id, name, group) bind(C)
            import :: c_int, c_char
            integer(c_int), value :: id
            character(kind=c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: group
        end function
        integer(c_int) function nc_def_dim(id, name, n, dim) bind(C)
            import :: c_int, c_char, c_size_t
            integer(c_int), value :: id
            character(kind=c_char), intent(in) :: name(*)
            integer(c_size_t), value :: n
            integer(c_int), intent(out) :: dim
        end function
        integer(c_int) function nc_def_var(id, name, xtype, ndims, dims, var) bind(C)
            import :: c_int, c_char
            integer(c_int), value :: id, xtype, ndims
            character(kind=c_char), intent(in) :: name(*)
            integer(c_int), intent(in) :: dims(*)
            integer(c_int), intent(out) :: var
        end function
        integer(c_int) function nc_put_att_text(id, var, name, n, text) bind(C)
            import :: c_int, c_char, c_size_t
            integer(c_int), value :: id, var
            character(kind=c_char), intent(in) :: name(*), text(*)
            integer(c_size_t), value :: n
        end function
        integer(c_int) function nc_put_att_double(id, var, name, xtype, n, values) bind(C)
            import :: c_int, c_char, c_size_t, c_double
            integer(c_int), value :: id, var, xtype
            character(kind=c_char), intent(in) :: name(*)
            integer(c_size_t), value :: n
            real(c_double), intent(in) :: values(*)
        end function
        integer(c_int) function nc_put_var_double(id, var, values) bind(C)
            import :: c_int, c_double
            integer(c_int), value :: id, var
            real(c_double), intent(in) :: values(*)
        end function
        integer(c_int) function nc_put_var_text(id, var, text) bind(C)
            import :: c_int, c_char
            integer(c_int), value :: id, var
            character(kind=c_char), intent(in) :: text(*)
        end function
        integer(c_int) function nc_inq_varid(id, name, var) bind(C)
            import :: c_int, c_char
            integer(c_int), value :: id
            character(kind=c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: var
        end function
        integer(c_int) function nc_inq_varndims(id, var, n) bind(C)
            import :: c_int
            integer(c_int), value :: id, var
            integer(c_int), intent(out) :: n
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
        integer(c_int) function nc_get_var_double(id, var, values) bind(C)
            import :: c_int, c_double
            integer(c_int), value :: id, var
            real(c_double), intent(out) :: values(*)
        end function
        integer(c_int) function nc_get_var_text(id, var, text) bind(C)
            import :: c_int, c_char
            integer(c_int), value :: id, var
            character(kind=c_char), intent(out) :: text(*)
        end function
        integer(c_int) function nc_get_att_double(id, var, name, values) bind(C)
            import :: c_int, c_char, c_double
            integer(c_int), value :: id, var
            character(kind=c_char), intent(in) :: name(*)
            real(c_double), intent(out) :: values(*)
        end function
    end interface

    public :: nc_open, nc_create, nc_close, nc_inq_ncid, nc_inq_dimid, nc_def_grp, nc_def_dim, nc_def_var, &
        nc_put_att_text, nc_put_att_double, nc_put_var_double, nc_put_var_text, nc_inq_varid, &
        nc_get_var_double, nc_get_att_double
    public :: cstr, nc_shape, nc_read_text, put_text_att

contains

    pure function cstr(text) result(c)
        character(len=*), intent(in) :: text
        character(kind=c_char, len=len(text) + 1) :: c
        c = text//c_null_char
    end function cstr

    integer(c_int) function nc_shape(id, var, shape_) result(code)
        !! Dimension lengths of a variable in file (C) order.
        integer(c_int), intent(in) :: id, var
        integer, allocatable, intent(out) :: shape_(:)
        integer(c_int) :: ndims, dims(16)
        integer(c_size_t) :: n
        integer :: k
        code = nc_inq_varndims(id, var, ndims)
        if (code /= NC_NOERR) return
        if (ndims > 16) then
            code = -1
            return
        end if
        code = nc_inq_vardimid(id, var, dims)
        if (code /= NC_NOERR) return
        allocate(shape_(ndims))
        do k = 1, ndims
            code = nc_inq_dimlen(id, dims(k), n)
            if (code /= NC_NOERR) return
            shape_(k) = int(n)
        end do
    end function nc_shape

    integer(c_int) function nc_read_text(id, var, width, n, text) result(code)
        !! A (n, width) character variable as n strings, NUL padding blanked.
        integer(c_int), intent(in) :: id, var
        integer, intent(in) :: width, n
        character(len=width), intent(out) :: text(n)
        character(kind=c_char) :: buf(width * n)
        integer :: i, j
        code = nc_get_var_text(id, var, buf)
        if (code /= NC_NOERR) return
        do i = 1, n
            do j = 1, width
                text(i)(j:j) = buf((i - 1) * width + j)
                if (iachar(text(i)(j:j)) == 0) text(i)(j:j) = ' '
            end do
        end do
    end function nc_read_text

    integer(c_int) function put_text_att(id, var, name, text) result(code)
        integer(c_int), intent(in) :: id, var
        character(len=*), intent(in) :: name, text
        code = nc_put_att_text(id, var, cstr(name), int(len(text), c_size_t), text)
    end function put_text_att
end module tiago_netcdf_c
