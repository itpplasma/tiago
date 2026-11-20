program tiago_axisym_fit
    !! Fit an axisymmetric free-boundary equilibrium to magnetic sensors, coil
    !! currents and profile rows through an equilibrium provider process
    !! (docs/axisymmetric-free-boundary.md). Usage: tiago_axisym_fit input.nml
    !!
    !! &axisym_fit
    !!   machine, provider, workdir      machine file, provider command, scratch directory
    !!   measurements                     CSV label,value,sigma: every sensor label of
    !!                                    the machine, and rows pressure@<psi_N> [Pa]
    !!                                    or jphi@<psi_N> [A/m^2]
    !!   nx, ny, rmin, rmax, zmin, zmax, tolerance, fvac, alpha_m, alpha_n
    !!   coil_current(g), plasma_current, pressure_axis
    !!   profile ('paxis' or 'knots'), knot_psin(k), pprime_knot(k), ffprime_knot(k)
    !!   parameters(i), start_values(i), parameter_scale(i), prior_values(i), prior_sigma(i)
    !!   max_iterations, fd_relative, draws, seed
    !!   edge_psin(j)                     labels of the reported outboard j_phi
    !!   output                           prefix of the result files
    !!   predict_geqdsk                   only predict the sensors of this G-EQDSK
    !!   initial                          plasma flux (.npy) every provider solve starts from
    !! /
    !!
    !! Writes <output>_parameters.csv (MAP, local standard deviation, start),
    !! <output>_signals.csv (measured, predicted, sigma), <output>_map.geqdsk and,
    !! with draws > 0, <output>_draw_NNNN.geqdsk of Laplace draws re-solved by the
    !! provider and <output>_edge.csv with j_phi at edge_psin for the MAP and
    !! every draw.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, error_unit
    use fortnum_status, only: fortnum_status_t
    use fortnum_rng, only: rng_t, rng_seed, rng_normal
    use tiago_axisym_free_boundary, only: axisym_free_boundary_t
    use tiago_axisym_machine, only: read_axisym_machine
    use tiago_geqdsk_plasma, only: geqdsk_plasma_t, read_geqdsk_plasma
    use tiago_inference_engines, only: map_fit, laplace, laplace_t
    implicit none

    integer, parameter :: mp = 64, mk = 32
    character(len=1024) :: machine, provider, workdir, measurements, output, nml, predict_geqdsk, initial
    character(len=16) :: profile
    character(len=32) :: parameters(mp)
    integer :: nx, ny, max_iterations, draws, seed, nknots
    real(dp) :: rmin, rmax, zmin, zmax, tolerance, fvac, alpha_m, alpha_n
    real(dp) :: coil_current(mp), plasma_current, pressure_axis, fd_relative
    real(dp) :: knot_psin(mk), pprime_knot(mk), ffprime_knot(mk), edge_psin(16)
    real(dp) :: start_values(mp), parameter_scale(mp), prior_values(mp), prior_sigma(mp)
    namelist /axisym_fit/ machine, provider, workdir, measurements, output, nx, ny, &
        rmin, rmax, zmin, zmax, tolerance, fvac, alpha_m, alpha_n, coil_current, &
        plasma_current, pressure_axis, profile, knot_psin, pprime_knot, ffprime_knot, &
        parameters, start_values, parameter_scale, prior_values, prior_sigma, &
        max_iterations, fd_relative, draws, seed, edge_psin, predict_geqdsk, initial

    type(axisym_free_boundary_t) :: model
    type(laplace_t) :: lap
    type(geqdsk_plasma_t) :: eq
    type(rng_t) :: rng
    type(fortnum_status_t) :: status
    character(len=:), allocatable :: message
    character(len=32) :: suffix
    real(dp), allocatable :: theta(:), start(:), s(:), z(:), td(:)
    integer :: unit, ierr, n, ng, ne, k, d, steps
    logical :: ok, conv

    machine = ''
    provider = ''
    workdir = 'axisym_work'
    measurements = ''
    output = 'axisym_fit'
    nx = 33
    ny = 33
    rmin = 0.1_dp
    rmax = 2.0_dp
    zmin = -1.0_dp
    zmax = 1.0_dp
    tolerance = 1.0e-6_dp
    fvac = 1.0_dp
    alpha_m = 1.0_dp
    alpha_n = 2.0_dp
    coil_current = 0
    plasma_current = 0
    pressure_axis = 0
    profile = 'paxis'
    knot_psin = -1
    pprime_knot = 0
    ffprime_knot = 0
    parameters = ''
    start_values = huge(1.0_dp)
    parameter_scale = 1
    prior_values = 0
    prior_sigma = 0
    max_iterations = 20
    fd_relative = 1.0e-3_dp
    draws = 0
    seed = 1
    edge_psin = -1
    predict_geqdsk = ''
    initial = ''
    call get_command_argument(1, nml)
    if (len_trim(nml) == 0) error stop 'usage: tiago_axisym_fit input.nml'
    open(newunit=unit, file=trim(nml), status='old', action='read')
    read(unit, nml=axisym_fit)
    close(unit)

    model%machine_path = trim(machine)
    call read_axisym_machine(model%machine_path, model%machine, ierr, message)
    if (ierr /= 0) error stop message
    model%command = trim(provider)
    model%workdir = trim(workdir)
    call execute_command_line('mkdir -p '//model%workdir)
    model%nx = nx
    model%ny = ny
    model%rmin = rmin
    model%rmax = rmax
    model%zmin = zmin
    model%zmax = zmax
    model%tolerance = tolerance
    model%fvac = fvac
    model%alpha_m = alpha_m
    model%alpha_n = alpha_n
    ng = model%machine%ngroups
    model%coil_current = coil_current(:ng)
    model%plasma_current = plasma_current
    model%pressure_axis = pressure_axis
    model%profile = trim(profile)
    nknots = count(knot_psin >= 0)
    model%knot_psin = knot_psin(:nknots)
    model%pprime_knot = pprime_knot(:nknots)
    model%ffprime_knot = ffprime_knot(:nknots)
    model%fd_relative = fd_relative
    model%initial = trim(initial)
    n = count(len_trim(parameters) > 0)
    model%names = parameters(:n)
    model%scale = parameter_scale(:n)
    model%prior_mean = prior_values(:n)
    model%prior_sigma = prior_sigma(:n)
    call read_measurements(trim(measurements))
    if (len_trim(predict_geqdsk) > 0) then
        call predict(trim(predict_geqdsk))
        stop
    end if

    allocate(theta(n))
    do k = 1, n
        if (start_values(k) == huge(1.0_dp)) then
            theta(k) = default_value(model%names(k)) / model%scale(k)
        else
            theta(k) = start_values(k) / model%scale(k)
        end if
    end do
    start = theta
    call map_fit(model, theta, max_iterations, steps, message, conv)
    write(*, '(A,I0,A,I0,A,L1)') 'fit: ', steps, ' steps, ', model%solves, &
        ' provider solves, converged ', conv
    if (.not. conv) write(error_unit, '(A)') 'fit: '//message
    call laplace(model, theta, lap)
    if (lap%info /= 0) error stop 'Laplace covariance failed'

    open(newunit=unit, file=trim(output)//'_parameters.csv', status='replace')
    write(unit, '(A)') 'name,value,sd,start'
    do k = 1, n
        write(unit, '(A,3(",",ES24.16))') trim(model%names(k)), theta(k) * model%scale(k), &
            sqrt(lap%cov(k, k)) * model%scale(k), start(k) * model%scale(k)
    end do
    close(unit)

    call model%signals(theta, s, ok, eq)
    if (.not. ok) error stop 'no equilibrium at the MAP'
    call execute_command_line('cp '//model%workdir//'/equilibrium.geqdsk '//trim(output) &
        //'_map.geqdsk')
    open(newunit=unit, file=trim(output)//'_signals.csv', status='replace')
    write(unit, '(A)') 'label,measured,predicted,sigma'
    do k = 1, size(s)
        write(unit, '(A,3(",",ES24.16))') trim(row_label(k)), model%data(k), s(k), model%sigma(k)
    end do
    close(unit)

    ne = count(edge_psin >= 0)
    if (ne > 0) then
        open(newunit=unit, file=trim(output)//'_edge.csv', status='replace')
        write(unit, '(A)') 'sample,psin,r,jphi'
        call write_edge(unit, 0, eq)
        call rng_seed(rng, int(seed, int64), status)
        allocate(z(n))
        do d = 1, draws
            do k = 1, n
                call rng_normal(rng, z(k))
            end do
            td = theta + matmul(lap%chol, z)
            call model%signals(td, s, ok, eq)
            if (.not. ok) then
                write(error_unit, '(A,I0,A)') 'draw ', d, ' did not converge, skipped'
                cycle
            end if
            write(suffix, '(A,I4.4,A)') '_draw_', d, '.geqdsk'
            call execute_command_line('cp '//model%workdir//'/equilibrium.geqdsk '//trim(output) &
                //trim(suffix))
            call write_edge(unit, d, eq)
        end do
        close(unit)
    end if

contains

    subroutine predict(path)
        !! Signals of a given G-EQDSK and the coil currents of the input,
        !! written to <output>_signals.csv, without a fit.
        character(*), intent(in) :: path
        type(geqdsk_plasma_t) :: given
        real(dp), allocatable :: rp(:), zp(:), wp(:), sig(:)
        real(dp) :: j, r
        integer :: unit, k, m
        call read_geqdsk_plasma(path, given, ierr, message)
        if (ierr /= 0) error stop message
        call given%filaments(rp, zp, wp)
        sig = model%machine%sensors%signals([model%machine%r, rp], [model%machine%z, zp], &
            [model%machine%turns * model%coil_current(model%machine%group), wp])
        m = size(sig)
        if (allocated(model%row_kind)) then
            do k = 1, size(model%row_kind)
                if (model%row_kind(k) == 1) then
                    sig = [sig, given%pressure(model%row_psin(k))]
                else
                    call given%midplane_jphi(model%row_psin(k), j, r)
                    sig = [sig, j]
                end if
            end do
        end if
        write(*, '(A,ES14.6,A)') 'plasma current of the G-EQDSK: ', sum(wp), ' A'
        open(newunit=unit, file=trim(output)//'_signals.csv', status='replace')
        write(unit, '(A)') 'label,measured,predicted,sigma'
        do k = 1, size(sig)
            write(unit, '(A,3(",",ES24.16))') trim(row_label(k)), model%data(k), sig(k), model%sigma(k)
        end do
        close(unit)
    end subroutine predict

    subroutine write_edge(unit, sample, eq)
        integer, intent(in) :: unit, sample
        type(geqdsk_plasma_t), intent(in) :: eq
        real(dp) :: j, r
        integer :: i
        do i = 1, ne
            call eq%midplane_jphi(edge_psin(i), j, r)
            write(unit, '(I0,3(",",ES24.16))') sample, edge_psin(i), r, j
        end do
    end subroutine write_edge

    real(dp) function default_value(name) result(v)
        character(*), intent(in) :: name
        integer :: g
        if (name == 'PLASMA_CURRENT') then
            v = plasma_current
        else if (name == 'PRESSURE_AXIS') then
            v = pressure_axis
        else if (index(name, 'CURRENT(') == 1) then
            read(name(9:index(name, ')') - 1), *) g
            v = coil_current(g)
        else if (index(name, 'PPRIME(') == 1) then
            read(name(8:index(name, ')') - 1), *) g
            v = pprime_knot(g)
        else if (index(name, 'FFPRIME(') == 1) then
            read(name(9:index(name, ')') - 1), *) g
            v = ffprime_knot(g)
        else
            error stop 'unknown parameter '//trim(name)
        end if
    end function default_value

    function row_label(k) result(label)
        integer, intent(in) :: k
        character(len=:), allocatable :: label
        integer :: nf, np, nr
        character(len=32) :: buf
        nf = size(model%machine%sensors%flux)
        np = size(model%machine%sensors%pickup)
        nr = size(model%machine%sensors%rogowski)
        if (k <= nf) then
            label = model%machine%sensors%flux(k)%label
        else if (k <= nf + np) then
            label = model%machine%sensors%pickup(k - nf)%label
        else if (k <= nf + np + nr) then
            label = model%machine%sensors%rogowski(k - nf - np)%label
        else
            write(buf, '(F0.4)') model%row_psin(k - nf - np - nr)
            if (model%row_kind(k - nf - np - nr) == 1) then
                label = 'pressure@'//trim(buf)
            else
                label = 'jphi@'//trim(buf)
            end if
        end if
    end function row_label

    subroutine read_measurements(path)
        !! Sensor rows in machine order, then the profile rows in file order.
        character(*), intent(in) :: path
        character(len=256) :: line, label
        real(dp) :: value, sigma
        real(dp), allocatable :: rv(:), rs(:), rp(:)
        integer, allocatable :: rk(:)
        logical, allocatable :: seen(:)
        integer :: unit, ios, i, nsens, at
        nsens = model%machine%sensors%count()
        allocate(model%data(nsens), model%sigma(nsens), seen(nsens), rv(0), rs(0), rp(0), rk(0))
        seen = .false.
        open(newunit=unit, file=path, status='old', action='read')
        read(unit, '(A)') line
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            do i = 1, len_trim(line)
                if (line(i:i) == ',') line(i:i) = ' '
            end do
            read(line, *) label, value, sigma
            at = index(label, '@')
            if (at > 0) then
                rk = [rk, merge(1, 2, label(:at - 1) == 'pressure')]
                if (label(:at - 1) /= 'pressure' .and. label(:at - 1) /= 'jphi') &
                    error stop 'unknown profile row '//trim(label)
                read(label(at + 1:), *) value
                rp = [rp, value]
                read(line, *) label, value, sigma
                rv = [rv, value]
                rs = [rs, sigma]
                cycle
            end if
            i = sensor_index(trim(label))
            if (i == 0) error stop 'measurement of an unknown sensor: '//trim(label)
            model%data(i) = value
            model%sigma(i) = sigma
            seen(i) = .true.
        end do
        close(unit)
        if (.not. all(seen)) error stop 'every sensor of the machine needs a measurement'
        if (size(rk) > 0) then
            model%row_kind = rk
            model%row_psin = rp
            model%data = [model%data, rv]
            model%sigma = [model%sigma, rs]
        end if
    end subroutine read_measurements

    integer function sensor_index(label) result(i)
        character(*), intent(in) :: label
        integer :: k, nf, np
        nf = size(model%machine%sensors%flux)
        np = size(model%machine%sensors%pickup)
        do k = 1, nf
            i = k
            if (model%machine%sensors%flux(k)%label == label) return
        end do
        do k = 1, np
            i = nf + k
            if (model%machine%sensors%pickup(k)%label == label) return
        end do
        do k = 1, size(model%machine%sensors%rogowski)
            i = nf + np + k
            if (model%machine%sensors%rogowski(k)%label == label) return
        end do
        i = 0
    end function sensor_index

end program tiago_axisym_fit
