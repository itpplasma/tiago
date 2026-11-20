module tiago_axisym_free_boundary
    !! Magnetic reconstruction on an axisymmetric free-boundary equilibrium from
    !! an external provider process (docs/axisymmetric-free-boundary.md).
    !! Each evaluation writes a request file, runs
    !!     <command> request.txt equilibrium.geqdsk
    !! in the work directory, reads the G-EQDSK and predicts the sensors of the
    !! machine from the coil currents and the plasma current. Parameters are
    !! CURRENT(g) [A], PLASMA_CURRENT [A], PRESSURE_AXIS [Pa] for the
    !! ConstrainPaxisIp profiles, or PPRIME(k) and FFPRIME(k), the values at knot
    !! k of p' and F F' in psi_N, which the provider scales together to the
    !! plasma current. They are fitted as value / scale. Besides the sensors,
    !! data rows may hold the pressure [Pa] or the outboard-midplane j_phi
    !! [A/m^2] at given psi_N. Derivatives are central differences with a cold provider
    !! solve per side, one-sided where one side does not converge.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_rng, only: rng_t, rng_normal
    use tiago_local_posterior, only: local_posterior_t
    use tiago_axisym_machine, only: axisym_machine_t
    use tiago_geqdsk_plasma, only: geqdsk_plasma_t, read_geqdsk_plasma
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    type, public, extends(local_posterior_t) :: axisym_free_boundary_t
        type(axisym_machine_t) :: machine
        character(len=:), allocatable :: machine_path, command, workdir
        !> Request values: grid nx, ny, rmin, rmax, zmin, zmax and the solver
        !> tolerance, then the physics fvac [T m], alpha_m, alpha_n.
        integer :: nx = 33, ny = 33
        real(dp) :: rmin = 0.1_dp, rmax = 2.0_dp, zmin = -1.0_dp, zmax = 1.0_dp
        real(dp) :: tolerance = 1.0e-6_dp, fvac = 1.0_dp, alpha_m = 1.0_dp, alpha_n = 2.0_dp
        !> Values of the parameters that are not fitted.
        real(dp), allocatable :: coil_current(:)
        real(dp) :: plasma_current = 0, pressure_axis = 0
        !> 'paxis' or 'knots', with the knots in psi_N and the shape values of
        !> p' and F F' there.
        character(len=:), allocatable :: profile
        real(dp), allocatable :: knot_psin(:), pprime_knot(:), ffprime_knot(:)
        !> Profile data rows after the sensors: row_kind 1 pressure, 2 j_phi
        !> at the outboard midplane, at the labels row_psin.
        integer, allocatable :: row_kind(:)
        real(dp), allocatable :: row_psin(:)
        !> Plasma flux the provider starts from (.npy). Without one, the first
        !> converged solve provides it, so that every later solve starts from
        !> the same state.
        character(len=:), allocatable :: initial
        character(len=32), allocatable :: names(:)
        real(dp), allocatable :: scale(:), prior_mean(:), prior_sigma(:)
        real(dp), allocatable :: data(:), sigma(:)
        real(dp) :: fd_relative = 1.0e-3_dp
        integer :: solves = 0
    contains
        procedure :: residual => fb_residual
        procedure :: jacobian => fb_jacobian
        procedure :: n_residuals => fb_n_residuals
        procedure :: log_normalization => fb_log_normalization
        procedure :: proper_prior => fb_proper_prior
        procedure :: draw_prior => fb_draw_prior
        procedure :: signals => fb_signals
    end type axisym_free_boundary_t

contains

    subroutine fb_signals(self, theta, s, ok, eq)
        !! Sensor signals of the equilibrium at theta.
        class(axisym_free_boundary_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), allocatable, intent(out) :: s(:)
        logical, intent(out) :: ok
        type(geqdsk_plasma_t), intent(out), optional :: eq
        type(geqdsk_plasma_t) :: solved
        real(dp), allocatable :: current(:), rp(:), zp(:), wp(:), pk(:), fk(:), rows(:)
        real(dp) :: ip, paxis, rmid
        character(len=:), allocatable :: message, request, output
        integer :: k, g, status, ierr, attempt
        current = self%coil_current
        ip = self%plasma_current
        paxis = self%pressure_axis
        if (.not. allocated(self%profile)) self%profile = 'paxis'
        if (allocated(self%pprime_knot)) then
            pk = self%pprime_knot
            fk = self%ffprime_knot
        else
            allocate(pk(0), fk(0))
        end if
        do k = 1, size(self%names)
            associate (name => self%names(k))
                if (name == 'PLASMA_CURRENT') then
                    ip = theta(k) * self%scale(k)
                else if (name == 'PRESSURE_AXIS') then
                    paxis = theta(k) * self%scale(k)
                else if (index(name, 'CURRENT(') == 1) then
                    read(name(9:index(name, ')') - 1), *) g
                    current(g) = theta(k) * self%scale(k)
                else if (index(name, 'PPRIME(') == 1) then
                    read(name(8:index(name, ')') - 1), *) g
                    pk(g) = theta(k) * self%scale(k)
                else if (index(name, 'FFPRIME(') == 1) then
                    read(name(9:index(name, ')') - 1), *) g
                    fk(g) = theta(k) * self%scale(k)
                else
                    error stop 'axisym free boundary: unknown parameter '//trim(name)
                end if
            end associate
        end do
        request = self%workdir//'/request.txt'
        output = self%workdir//'/equilibrium.geqdsk'
        call write_request(request, .true.)
        ! FreeGSNKE's Newton-Krylov directions are random, so a failed solve
        ! is repeated once, and then once more from the provider's own initial
        ! flux, before the point counts as failed.
        do attempt = 1, 3
            if (attempt == 3) call write_request(request, .false.)
            call execute_command_line(self%command//' '//request//' '//output, exitstat=status)
            self%solves = self%solves + 1
            ok = status == 0
            if (ok) ok = converged(output//'.status')
            if (ok) then
                call read_geqdsk_plasma(output, solved, ierr, message)
                ok = ierr == 0
                if (.not. ok) then
                    write(error_unit, '(A)') 'axisym free boundary: '//message
                    call execute_command_line('cp '//output//' '//output//'.failed')
                end if
            end if
            if (ok) exit
        end do
        if (.not. ok) return
        if (.not. allocated(self%initial)) self%initial = ''
        if (len(self%initial) == 0) then
            self%initial = self%workdir//'/initial.psi.npy'
            call execute_command_line('cp '//output//'.psi.npy '//self%initial)
        end if
        call solved%filaments(rp, zp, wp)
        s = self%machine%sensors%signals([self%machine%r, rp], [self%machine%z, zp], &
            [self%machine%turns * current(self%machine%group), wp])
        if (allocated(self%row_kind)) then
            allocate(rows(size(self%row_kind)))
            do k = 1, size(rows)
                if (self%row_kind(k) == 1) then
                    rows(k) = solved%pressure(self%row_psin(k))
                else
                    call solved%midplane_jphi(self%row_psin(k), rows(k), rmid)
                end if
            end do
            s = [s, rows]
        end if
        if (present(eq)) eq = solved
    contains
        subroutine write_request(path, use_initial)
            character(*), intent(in) :: path
            logical, intent(in) :: use_initial
            integer :: unit, i
            open(newunit=unit, file=path, status='replace', action='write')
            write(unit, '(A)') 'machine = '//self%machine_path
            write(unit, '(A,I0)') 'nx = ', self%nx
            write(unit, '(A,I0)') 'ny = ', self%ny
            write(unit, '(A,ES24.16)') 'rmin = ', self%rmin
            write(unit, '(A,ES24.16)') 'rmax = ', self%rmax
            write(unit, '(A,ES24.16)') 'zmin = ', self%zmin
            write(unit, '(A,ES24.16)') 'zmax = ', self%zmax
            write(unit, '(A,ES24.16)') 'tolerance = ', self%tolerance
            write(unit, '(A,ES24.16)') 'fvac = ', self%fvac
            write(unit, '(A,ES24.16)') 'alpha_m = ', self%alpha_m
            write(unit, '(A,ES24.16)') 'alpha_n = ', self%alpha_n
            write(unit, '(A,ES24.16)') 'plasma_current = ', ip
            write(unit, '(A,ES24.16)') 'pressure_axis = ', paxis
            if (allocated(self%initial) .and. use_initial) then
                if (len(self%initial) > 0) write(unit, '(A)') 'initial = '//self%initial
            end if
            if (self%profile == 'knots') then
                write(unit, '(A)') 'profile = knots'
                write(unit, '(A,*(ES24.16,:,","))') 'psin_knots = ', self%knot_psin
                write(unit, '(A,*(ES24.16,:,","))') 'pprime_knots = ', pk
                write(unit, '(A,*(ES24.16,:,","))') 'ffprime_knots = ', fk
            end if
            do i = 1, size(current)
                write(unit, '(A,I0,A,ES24.16)') 'current(', i, ') = ', current(i)
            end do
            close(unit)
        end subroutine write_request
    end subroutine fb_signals

    logical function converged(path)
        character(*), intent(in) :: path
        character(len=32) :: key
        integer :: unit, ios, flag
        converged = .false.
        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        read(unit, *, iostat=ios) key, flag
        close(unit)
        converged = ios == 0 .and. trim(key) == 'converged' .and. flag == 1
    end function converged

    subroutine fb_residual(self, theta, r, ok)
        class(axisym_free_boundary_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), intent(out) :: r(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: s(:)
        integer :: n, k, j
        call self%signals(theta, s, ok)
        if (.not. ok) return
        n = size(self%data)
        r(:n) = (s - self%data) / self%sigma
        j = n
        do k = 1, size(theta)
            if (self%prior_sigma(k) <= 0) cycle
            j = j + 1
            r(j) = (theta(k) * self%scale(k) - self%prior_mean(k)) / self%prior_sigma(k)
        end do
    end subroutine fb_residual

    subroutine fb_jacobian(self, theta, jac)
        class(axisym_free_boundary_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), allocatable, intent(out) :: jac(:, :)
        real(dp), allocatable :: rp(:), rm(:), t(:), r0(:)
        real(dp) :: h
        logical :: okp, okm, ok0
        integer :: k, attempt
        allocate(jac(self%n_residuals(), size(theta)), rp(self%n_residuals()), &
            rm(self%n_residuals()))
        do k = 1, size(theta)
            h = self%fd_relative * max(abs(theta(k)), 1.0_dp)
            ! Where a side does not converge, shrink the step twice by 4, then
            ! fall back to one side.
            do attempt = 1, 3
                t = theta
                t(k) = theta(k) + h
                call self%residual(t, rp, okp)
                t(k) = theta(k) - h
                call self%residual(t, rm, okm)
                if (okp .and. okm) exit
                h = h / 4
            end do
            if (okp .and. okm) then
                jac(:, k) = (rp - rm) / (2 * h)
                cycle
            end if
            if (.not. allocated(r0)) then
                allocate(r0(self%n_residuals()))
                call self%residual(theta, r0, ok0)
                if (.not. ok0) error stop 'axisym free boundary: no equilibrium at the point'
            end if
            if (okp) then
                jac(:, k) = (rp - r0) / h
            else if (okm) then
                jac(:, k) = (r0 - rm) / h
            else
                write(error_unit, '(A)') 'axisym free boundary: no derivative for '// &
                    trim(self%names(k))
                error stop 'axisym free boundary: derivative solve failed'
            end if
        end do
    end subroutine fb_jacobian

    integer function fb_n_residuals(self) result(n)
        class(axisym_free_boundary_t), intent(in) :: self
        n = size(self%data) + count(self%prior_sigma > 0)
    end function fb_n_residuals

    real(dp) function fb_log_normalization(self) result(v)
        class(axisym_free_boundary_t), intent(in) :: self
        v = -sum(log(sqrt(2 * pi) * self%sigma)) &
            - sum(log(sqrt(2 * pi) * self%prior_sigma), mask=self%prior_sigma > 0)
    end function fb_log_normalization

    logical function fb_proper_prior(self)
        class(axisym_free_boundary_t), intent(in) :: self
        fb_proper_prior = all(self%prior_sigma > 0)
    end function fb_proper_prior

    subroutine fb_draw_prior(self, rng, theta, ok)
        class(axisym_free_boundary_t), intent(in) :: self
        type(rng_t), intent(inout) :: rng
        real(dp), intent(out) :: theta(:)
        logical, intent(out) :: ok
        real(dp) :: z
        integer :: k
        ok = self%proper_prior()
        if (.not. ok) return
        do k = 1, size(theta)
            call rng_normal(rng, z)
            theta(k) = (self%prior_mean(k) + self%prior_sigma(k) * z) / self%scale(k)
        end do
    end subroutine fb_draw_prior

end module tiago_axisym_free_boundary
