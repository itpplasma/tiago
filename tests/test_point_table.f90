program test_point_table
    !! Point-table provider (tiago_point_table) through the synthetic driver.
    !! Usage: test_point_table <workdir> <wout.nc>
    !! 1. Island fixture: concentric circles (R0 = 3, a = 1, rho = r) with an
    !!    island disk of radius 0.1 at (R, Z) = (3.5, 0) whose interior carries
    !!    the plateau label rho = 0.5. A Thomson point inside the island gives
    !!    the plateau ne exactly; the vertical chord through the island equals
    !!    the hand integral with the plateau (Simpson quadrature here).
    !! 2. The LI383 VMEC state sampled on a grid: Thomson points on grid nodes
    !!    give the s, ne and Te of the native VMEC path to round-off; a chord
    !!    inside the boundary agrees within the trilinear interpolation error of
    !!    the 1 cm grid.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use, intrinsic :: iso_c_binding, only: c_int, c_size_t
    use tiago_netcdf_c, only: NC_NOERR, NC_NOWRITE, NC_CLOBBER, NC_NETCDF4, NC_GLOBAL, NC_CHAR, &
        NC_DOUBLE, nc_open, nc_create, nc_close, nc_inq_ncid, nc_def_grp, nc_def_dim, nc_def_var, &
        nc_put_var_double, nc_put_att_double, nc_inq_varid, nc_get_var_double, cstr
    use tiago_real_space_state, only: located_t, locate_ok
    use tiago_vmec_state, only: vmec_state_t
    use tiago_synthetic_forward, only: synth_config_t, synth_state_t, synth_load, synth_write
    use tiago_profile_geometry_io, only: profile_table_t
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp), n0 = 5.0e19_dp, t0 = 3.0e3_dp
    integer :: failures
    character(len=1024) :: workdir, wout

    failures = 0
    call get_command_argument(1, workdir)
    call get_command_argument(2, wout)
    call write_profiles(trim(workdir)//'/table_profiles.nc')
    call island_case()
    call vmec_case()
    if (failures > 0) then
        write(error_unit, '(I0,A)') failures, ' check(s) failed'
        error stop 1
    end if
    print '(A)', 'test_point_table passed'

contains

    pure real(dp) function ne_of(s)
        real(dp), intent(in) :: s
        ne_of = n0 * exp(-0.4_dp * s - 1.1_dp * s**2 + 0.3_dp * s**3)
    end function ne_of

    pure real(dp) function te_of(s)
        real(dp), intent(in) :: s
        te_of = t0 * exp(-1.0_dp * s - 0.9_dp * s**2 + 0.2_dp * s**3)
    end function te_of

    subroutine island_case()
        integer, parameter :: nr = 401, nz = 401
        real(dp) :: r(nr), z(nz), rho(nr, nz, 1), b(nr, nz, 1), st(nr, nz, 1), isl(nr, nz, 1)
        real(dp) :: chord(3, 2), ts(3, 2), u0, hand, h, zz
        real(dp), allocatable :: v(:)
        integer :: i, k, m

        r = [(2.0_dp + 0.005_dp * (i - 1), i = 1, nr)]
        z = [(-1.0_dp + 0.005_dp * (i - 1), i = 1, nz)]
        do k = 1, nz
            do i = 1, nr
                rho(i, k, 1) = hypot(r(i) - 3.0_dp, z(k))
                b(i, k, 1) = 2.5_dp * 3.0_dp / r(i)
                st(i, k, 1) = merge(0.0_dp, 1.0_dp, rho(i, k, 1) <= 1.0_dp)
                isl(i, k, 1) = merge(1.0_dp, 0.0_dp, hypot(r(i) - 3.5_dp, z(k)) < 0.1_dp)
            end do
        end do
        call write_table(trim(workdir)//'/island_table.nc', 1, r, z, [0.0_dp], rho, &
            0.0_dp * b, b, 0.0_dp * b, st, isl, [0.5_dp])

        u0 = sqrt(1.0_dp - 0.25_dp)
        chord(:, 1) = [3.5_dp, 0.0_dp, -u0]
        chord(:, 2) = [3.5_dp, 0.0_dp, u0]
        ts(:, 1) = [3.5_dp, 0.0_dp, 0.0225_dp] ! inside the island, between nodes
        ts(:, 2) = [3.5_dp * cos(2.0_dp), 3.5_dp * sin(2.0_dp), 0.4_dp] ! outside it
        call write_geometry(trim(workdir)//'/island_geometry.nc', chord, ts)
        call run('island_table.nc', 'island_geometry.nc', '', 'island_synthetic.nc')

        call read_var(trim(workdir)//'/island_synthetic.nc', 'thomson', 'ne', v)
        call close('island plateau ne', v(1), ne_of(0.25_dp), 1.0e-13_dp)
        call close('ne outside the island', v(2), ne_of(0.25_dp + 0.16_dp), 1.0e-4_dp)
        ! Hand integral: ne(s = 0.25 + z^2) outside |z| < 0.1, plateau ne(0.25) inside.
        m = 200000
        h = 2.0_dp * u0 / m
        hand = 0.0_dp
        do i = 0, m
            zz = -u0 + h * i
            hand = hand + merge(1.0_dp, merge(4.0_dp, 2.0_dp, mod(i, 2) == 1), i == 0 .or. i == m) &
                * merge(ne_of(0.25_dp), ne_of(0.25_dp + zz**2), abs(zz) < 0.1_dp)
        end do
        hand = hand * h / 3.0_dp
        call read_var(trim(workdir)//'/island_synthetic.nc', 'interferometer', 'line_density', v)
        call close('chord through the island vs hand integral', v(1), hand, 2.0e-4_dp)
        call check(abs(v(1) - nested_integral(u0)) > 10.0_dp * abs(v(1) - hand), &
            'the plateau changes the chord integral by more than the interpolation error')
    end subroutine island_case

    real(dp) function nested_integral(u0)
        real(dp), intent(in) :: u0
        integer, parameter :: m = 200000
        real(dp) :: h
        integer :: i
        h = 2.0_dp * u0 / m
        nested_integral = 0.0_dp
        do i = 0, m
            nested_integral = nested_integral + merge(1.0_dp, merge(4.0_dp, 2.0_dp, mod(i, 2) == 1), &
                i == 0 .or. i == m) * ne_of(0.25_dp + (-u0 + h * i)**2)
        end do
        nested_integral = nested_integral * h / 3.0_dp
    end function nested_integral

    subroutine vmec_case()
        integer, parameter :: nr = 91, nz = 141, np = 12
        type(vmec_state_t) :: vm
        type(located_t) :: loc
        real(dp) :: r(nr), z(nz), phi(np), rho(nr, nz, np), br(nr, nz, np), bp(nr, nz, np), &
            bz(nr, nz, np), st(nr, nz, np), x(3), chord(3, 2), ts(3, 6)
        real(dp), allocatable :: a(:), b(:)
        integer :: i, k, j, q, ierr
        character(len=:), allocatable :: message

        call vm%load(trim(wout), ierr, message, vector=.true.)
        call check(ierr == 0, 'load '//trim(wout))
        r = [(1.0_dp + 0.01_dp * (i - 1), i = 1, nr)]
        z = [(-0.7_dp + 0.01_dp * (i - 1), i = 1, nz)]
        phi = [(2.0_dp * pi / 3.0_dp * (i - 1) / np, i = 1, np)]
        !$omp parallel do private(i, k, x, loc) collapse(2)
        do j = 1, np
            do k = 1, nz
                do i = 1, nr
                    x = [r(i) * cos(phi(j)), r(i) * sin(phi(j)), z(k)]
                    call vm%locate(x, loc)
                    st(i, k, j) = 2.0_dp
                    rho(i, k, j) = 0.0_dp
                    br(i, k, j) = 0.0_dp
                    bp(i, k, j) = 0.0_dp
                    bz(i, k, j) = 0.0_dp
                    if (loc%status /= locate_ok) cycle
                    st(i, k, j) = merge(0.0_dp, 1.0_dp, loc%u(1) <= 1.0_dp)
                    rho(i, k, j) = loc%u(1)
                    ! |B| of the state's scalar field, as the native path samples it.
                    br(i, k, j) = loc%b
                end do
            end do
        end do
        !$omp end parallel do
        call write_table(trim(workdir)//'/vmec_table.nc', 3, r, z, phi, rho, br, bp, bz, st)
        ! Thomson points on grid nodes (inside), chord in the phi = 0 plane.
        ! Six interior nodes, 0.2 < rho < 0.9, in different planes.
        q = 0
        outer: do j = 1, np, 2
            do k = 1, nz, 7
                do i = 1, nr, 5
                    if (st(i, k, j) /= 0.0_dp .or. rho(i, k, j) < 0.2_dp .or. rho(i, k, j) > 0.9_dp) cycle
                    q = q + 1
                    ts(:, q) = [r(i) * cos(phi(j)), r(i) * sin(phi(j)), z(k)]
                    if (q == 6) exit outer
                    exit
                end do
            end do
        end do outer
        call check(q == 6, 'six interior grid nodes')
        ! Vertical chord inside the boundary at phi = 0.
        chord(:, 1) = [1.45_dp, 0.0_dp, 0.1_dp]
        chord(:, 2) = [1.45_dp, 0.0_dp, 0.35_dp]
        call write_geometry(trim(workdir)//'/vmec_table_geometry.nc', chord, ts)
        call run('vmec_table.nc', 'vmec_table_geometry.nc', '', 'vmec_table_synthetic.nc')
        call run('', 'vmec_table_geometry.nc', trim(wout), 'vmec_native_synthetic.nc')
        call compare('s', 6)
        call compare('ne', 6)
        call compare('te', 6)
        call read_var(trim(workdir)//'/vmec_table_synthetic.nc', 'interferometer', 'line_density', a)
        call read_var(trim(workdir)//'/vmec_native_synthetic.nc', 'interferometer', 'line_density', b)
        call close('chord on the 1 cm grid vs native VMEC', a(1), b(1), 1.0e-3_dp)
    end subroutine vmec_case

    subroutine compare(name, n)
        character(len=*), intent(in) :: name
        integer, intent(in) :: n
        real(dp), allocatable :: a(:), b(:)
        integer :: q
        call read_var(trim(workdir)//'/vmec_table_synthetic.nc', 'thomson', name, a)
        call read_var(trim(workdir)//'/vmec_native_synthetic.nc', 'thomson', name, b)
        do q = 1, n
            call close('grid node '//name//' vs native VMEC', a(q), b(q), 1.0e-12_dp)
        end do
    end subroutine compare

    subroutine run(table, geometry, wout_path, output)
        character(len=*), intent(in) :: table, geometry, wout_path, output
        type(synth_config_t) :: cfg
        type(synth_state_t) :: st
        integer :: ierr
        character(len=:), allocatable :: message
        cfg%wout = wout_path
        cfg%point_table = ''
        if (len(table) > 0) cfg%point_table = trim(workdir)//'/'//table
        cfg%geometry = trim(workdir)//'/'//geometry
        cfg%profiles = trim(workdir)//'/table_profiles.nc'
        cfg%output = trim(workdir)//'/'//output
        cfg%sensors = ''
        cfg%coils = ''
        cfg%extcur = ''
        cfg%sheet = 'contravariant'
        call synth_load(cfg, st, ierr, message)
        call check(ierr == 0, 'synth_load '//table//': '//message)
        if (ierr /= 0) return
        call synth_write(cfg, st, ierr, message)
        call check(ierr == 0, 'synth_write '//table//': '//message)
    end subroutine run

    subroutine write_table(path, nfp, r, z, phi, rho, br, bp, bz, st, isl, plateau)
        character(len=*), intent(in) :: path
        integer, intent(in) :: nfp
        real(dp), intent(in) :: r(:), z(:), phi(:), rho(:, :, :), br(:, :, :), bp(:, :, :), &
            bz(:, :, :), st(:, :, :)
        real(dp), intent(in), optional :: isl(:, :, :), plateau(:)
        integer(c_int) :: ncid, dr, dz, dph, di, s, d3(3)
        s = nc_create(cstr(path), ior(NC_CLOBBER, NC_NETCDF4), ncid)
        s = nc_put_att_double(ncid, NC_GLOBAL, cstr('nfp'), NC_DOUBLE, 1_c_size_t, [real(nfp, dp)])
        s = nc_def_dim(ncid, cstr('r'), int(size(r), c_size_t), dr)
        s = nc_def_dim(ncid, cstr('z'), int(size(z), c_size_t), dz)
        s = nc_def_dim(ncid, cstr('phi'), int(size(phi), c_size_t), dph)
        d3 = [dph, dz, dr]
        call put_var(ncid, 'R', [dr], r)
        call put_var(ncid, 'Z', [dz], z)
        call put_var(ncid, 'phi', [dph], phi)
        call put_var(ncid, 'rho', d3, reshape(rho, [size(rho)]))
        call put_var(ncid, 'br', d3, reshape(br, [size(br)]))
        call put_var(ncid, 'bphi', d3, reshape(bp, [size(bp)]))
        call put_var(ncid, 'bz', d3, reshape(bz, [size(bz)]))
        call put_var(ncid, 'status', d3, reshape(st, [size(st)]))
        if (present(isl)) then
            call put_var(ncid, 'island', d3, reshape(isl, [size(isl)]))
            s = nc_def_dim(ncid, cstr('island_count'), int(size(plateau), c_size_t), di)
            call put_var(ncid, 'island_rho', [di], plateau)
        end if
        s = nc_close(ncid)
    end subroutine write_table

    subroutine put_var(ncid, name, dims, values)
        integer(c_int), intent(in) :: ncid, dims(:)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: values(:)
        integer(c_int) :: var, s
        s = nc_def_var(ncid, cstr(name), NC_DOUBLE, int(size(dims), c_int), dims, var)
        if (s == NC_NOERR) s = nc_put_var_double(ncid, var, values)
        call check(s == NC_NOERR, 'write '//name)
    end subroutine put_var

    subroutine write_geometry(path, chord, centres)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: chord(3, 2), centres(:, :)
        integer(c_int) :: ncid, gid, dn, dx, var, s
        real(dp) :: dir(3, size(centres, 2))
        s = nc_create(cstr(path), ior(NC_CLOBBER, NC_NETCDF4), ncid)
        s = nc_def_grp(ncid, cstr('interferometer'), gid)
        s = nc_def_dim(gid, cstr('chord'), 1_c_size_t, dn)
        s = nc_def_dim(gid, cstr('xyz'), 3_c_size_t, dx)
        s = nc_def_var(gid, cstr('start'), NC_DOUBLE, 2_c_int, [dn, dx], var)
        s = nc_put_var_double(gid, var, chord(:, 1))
        s = nc_def_var(gid, cstr('end'), NC_DOUBLE, 2_c_int, [dn, dx], var)
        s = nc_put_var_double(gid, var, chord(:, 2))
        s = nc_def_grp(ncid, cstr('thomson'), gid)
        s = nc_def_dim(gid, cstr('volume'), int(size(centres, 2), c_size_t), dn)
        s = nc_def_dim(gid, cstr('xyz'), 3_c_size_t, dx)
        s = nc_def_var(gid, cstr('center'), NC_DOUBLE, 2_c_int, [dn, dx], var)
        s = nc_put_var_double(gid, var, centres)
        dir = 0.0_dp
        dir(3, :) = 1.0_dp
        s = nc_def_var(gid, cstr('direction'), NC_DOUBLE, 2_c_int, [dn, dx], var)
        s = nc_put_var_double(gid, var, dir)
        s = nc_def_var(gid, cstr('length'), NC_DOUBLE, 1_c_int, [dn], var)
        s = nc_put_var_double(gid, var, spread(0.0_dp, 1, size(centres, 2)))
        call check(s == NC_NOERR, 'write '//path)
        s = nc_close(ncid)
    end subroutine write_geometry

    subroutine write_profiles(path)
        character(len=*), intent(in) :: path
        integer(c_int) :: ncid, dn, var, s
        real(dp) :: sv(31)
        integer :: k
        sv = [(real(k, dp) / 30, k = 0, 30)]
        s = nc_create(cstr(path), ior(NC_CLOBBER, NC_NETCDF4), ncid)
        s = nc_def_dim(ncid, cstr('s'), 31_c_size_t, dn)
        s = nc_def_var(ncid, cstr('s'), NC_DOUBLE, 1_c_int, [dn], var)
        s = nc_put_var_double(ncid, var, sv)
        s = nc_def_var(ncid, cstr('ne'), NC_DOUBLE, 1_c_int, [dn], var)
        s = nc_put_var_double(ncid, var, [(ne_of(sv(k)), k = 1, 31)])
        s = nc_def_var(ncid, cstr('te'), NC_DOUBLE, 1_c_int, [dn], var)
        s = nc_put_var_double(ncid, var, [(te_of(sv(k)), k = 1, 31)])
        s = nc_def_var(ncid, cstr('ti'), NC_DOUBLE, 1_c_int, [dn], var)
        s = nc_put_var_double(ncid, var, [(te_of(sv(k)), k = 1, 31)])
        call check(s == NC_NOERR, 'write '//path)
        s = nc_close(ncid)
    end subroutine write_profiles

    subroutine read_var(path, group, name, v)
        character(len=*), intent(in) :: path, group, name
        real(dp), allocatable, intent(out) :: v(:)
        integer(c_int) :: ncid, gid, var, s
        allocate(v(64))
        v = huge(1.0_dp)
        s = nc_open(cstr(path), NC_NOWRITE, ncid)
        if (s == NC_NOERR) s = nc_inq_ncid(ncid, cstr(group), gid)
        if (s == NC_NOERR) s = nc_inq_varid(gid, cstr(name), var)
        if (s == NC_NOERR) s = nc_get_var_double(gid, var, v)
        call check(s == NC_NOERR, 'read '//path//' '//group//'/'//name)
        s = nc_close(ncid)
    end subroutine read_var

    subroutine check(ok, name)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: name
        if (.not. ok) then
            write(error_unit, '(A)') 'FAIL '//name
            failures = failures + 1
        end if
    end subroutine check

    subroutine close(name, a, b, rtol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: a, b, rtol
        if (abs(a - b) > rtol * abs(b) .or. a /= a) then
            write(error_unit, '(A,2ES24.16)') 'FAIL '//name//': ', a, b
            failures = failures + 1
        else
            print '(A,2ES24.16)', 'ok   '//name//': ', a, b
        end if
    end subroutine close
end program test_point_table
