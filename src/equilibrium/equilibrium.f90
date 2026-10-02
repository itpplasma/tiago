module tiago_equilibrium
    !! Fixed-boundary VMEC++ equilibria as functions of named input parameters,
    !! the edge vector y Tiago's plasma model reads, and the exact Jacobian
    !! C (dy/dx) for a cotangent matrix C = dS/dy:
    !!
    !!   edge VJP (tiago_vmec_edge)          C -> geometry and direct profile cotangents
    !!   VMEC++ implicit adjoint (adapter)   geometry -> boundary (rbc, zbs) and implicit
    !!                                       half-grid (mu0 p, iota, current) cotangents
    !!   power-series profiles (here)        -> PRES_SCALE, AM(k), AI(k), AC(k), CURTOR
    !!   MHD scaling (here)                  -> PHIEDGE
    !!
    !! PHIEDGE is not an input VMEC++ differentiates. Ideal-MHD equilibria are
    !! invariant under B -> l B, p -> l^2 p, I -> l I at fixed geometry and iota,
    !! so d y / d ln(phiedge) = y_B - 2 d y / d ln(PRES_SCALE) - d y / d ln(CURTOR)
    !! (the last term for ncurr = 1), with y_B the field part of y.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_vmecpp, only: vmecpp_t
    use tiago_vmec_edge, only: vmec_edge_t
    use tiago_plasma_support, only: vmec_boundary_t
    use tiago_vmecpp_paths, only: indata2json_path
    implicit none
    private

    real(dp), parameter :: mu0 = 4.0e-7_dp * acos(-1.0_dp)
    integer, parameter :: K_SCALAR = 1, K_ARRAY = 2, K_BOUNDARY = 3

    type, public :: equilibrium_t
        type(vmecpp_t) :: vmec
        type(vmec_edge_t) :: edge
        character(len=64), allocatable :: names(:)
        character(len=16), allocatable :: key(:)     !! adapter input name
        integer, allocatable :: kind(:), i1(:), i2(:)
        integer :: ns = 0, mpol = 0, ntor = 0, nfp = 1, ncurr = 0
        logical :: solved = .false., conservative = .false.
        integer :: max_block = 8000
        ! the last solve
        real(dp), allocatable :: coef(:, :, :, :), y(:)
        real(dp), allocatable :: phip_full(:), phip_half(:), iota_half(:), current_half(:)
        real(dp), allocatable :: mass_half(:)
        real(dp) :: lamscale = 0.0_dp
        integer :: signgs = 1
        ! bookkeeping
        integer :: solves = 0, adjoint_solves = 0
        real(dp) :: t_solve = 0.0_dp, t_adjoint = 0.0_dp
    contains
        procedure :: init => equilibrium_init
        procedure :: values => equilibrium_values
        procedure :: value_of => equilibrium_value_of
        procedure :: solve => equilibrium_solve
        procedure :: boundary => equilibrium_boundary
        procedure :: jacobian => equilibrium_jacobian
        procedure :: write_wout => equilibrium_write_wout
        procedure :: n_y => equilibrium_n_y
    end type equilibrium_t

contains

    subroutine equilibrium_init(self, input_path, workdir, names, covariant, &
            conservative)
        !! input_path: VMEC++ JSON or classic &INDATA (converted with indata2json).
        class(equilibrium_t), intent(inout) :: self
        character(len=*), intent(in) :: input_path, workdir
        character(len=*), intent(in) :: names(:)
        logical, intent(in), optional :: covariant, conservative
        character(len=:), allocatable :: json
        integer :: k
        logical :: use_covariant

        self%conservative = .false.
        if (present(conservative)) self%conservative = conservative
        use_covariant = self%conservative
        if (present(covariant)) use_covariant = use_covariant .or. covariant
        json = json_input(input_path, workdir)
        call self%vmec%create(json)
        if (self%vmec%get_int('lfreeb') /= 0) error stop 'VMEC input must be fixed-boundary'
        self%ns = self%vmec%get_int('ns')
        self%mpol = self%vmec%get_int('mpol')
        self%ntor = self%vmec%get_int('ntor')
        self%nfp = self%vmec%get_int('nfp')
        self%ncurr = self%vmec%get_int('ncurr')
        call self%edge%init(self%ns, self%mpol, self%ntor, self%nfp, self%ncurr, &
            self%vmec%get_int('ntheta'), self%vmec%get_int('nzeta'), use_covariant)
        allocate(self%names(size(names)), self%key(size(names)), self%kind(size(names)), &
            self%i1(size(names)), self%i2(size(names)))
        do k = 1, size(names)
            self%names(k) = names(k)
            call parse(names(k), self%key(k), self%kind(k), self%i1(k), self%i2(k))
        end do
    contains
        subroutine parse(name, key, kind, i1, i2)
            !! PRES_SCALE, CURTOR, PHIEDGE; AM(k), AI(k), AC(k) (0-based k);
            !! RBC(n,m), ZBS(n,m) as in the VMEC namelist.
            character(len=*), intent(in) :: name
            character(len=16), intent(out) :: key
            integer, intent(out) :: kind, i1, i2
            integer :: open, comma, ios, a, b
            character(len=64) :: base

            open = index(name, '(')
            base = lower(adjustl(name))
            if (open == 0) then
                key = trim(base)
                kind = K_SCALAR
                i1 = 1
                i2 = 1
                if (all(key /= [character(len=16) :: 'pres_scale', 'curtor', 'phiedge'])) then
                    error stop 'unsupported equilibrium parameter: '//trim(name)
                end if
                return
            end if
            key = lower(adjustl(name(:open - 1)))
            comma = index(name, ',')
            if (comma == 0) then
                read(name(open + 1:index(name, ')') - 1), *, iostat=ios) a
                if (ios /= 0 .or. a < 0) error stop 'bad parameter: '//trim(name)
                if (all(key /= [character(len=16) :: 'am', 'ai', 'ac'])) then
                    error stop 'unsupported equilibrium parameter: '//trim(name)
                end if
                kind = K_ARRAY
                i1 = a + 1
                i2 = 1
            else
                read(name(open + 1:comma - 1), *, iostat=ios) a      ! n
                if (ios == 0) read(name(comma + 1:index(name, ')') - 1), *, iostat=ios) b  ! m
                if (ios /= 0) error stop 'bad parameter: '//trim(name)
                if (all(key /= [character(len=16) :: 'rbc', 'zbs'])) then
                    error stop 'unsupported equilibrium parameter: '//trim(name)
                end if
                if (b < 0 .or. b >= self%mpol .or. abs(a) > self%ntor) then
                    error stop 'boundary mode outside mpol/ntor: '//trim(name)
                end if
                kind = K_BOUNDARY
                i1 = b + 1                 ! m + 1
                i2 = a + self%ntor + 1     ! n + ntor + 1
            end if
        end subroutine parse
    end subroutine equilibrium_init

    function json_input(path, workdir) result(json)
        character(len=*), intent(in) :: path, workdir
        character(len=:), allocatable :: json
        character(len=1) :: first
        integer :: unit, ios, status

        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'cannot open VMEC input '//path
        first = ' '
        do
            read(unit, '(A1)', advance='no', iostat=ios) first
            if (ios /= 0 .or. (first /= ' ' .and. first /= achar(9))) exit
        end do
        close(unit)
        if (first == '{') then
            json = path
            return
        end if
        call execute_command_line('cp "'//path//'" "'//workdir//'/input.tiago" && cd "'// &
            workdir//'" && "'//indata2json_path//'" input.tiago > indata2json.log 2>&1', &
            exitstat=status)
        if (status /= 0) error stop 'indata2json failed on '//path//' (see '//workdir//'/indata2json.log)'
        json = workdir//'/tiago.json'
    end function json_input

    function equilibrium_values(self) result(x)
        !! Current input values of the parameters.
        class(equilibrium_t), intent(in) :: self
        real(dp) :: x(size(self%names))
        integer :: k
        do k = 1, size(x)
            x(k) = self%vmec%get_input(trim(self%key(k)), self%i1(k), self%i2(k))
        end do
    end function equilibrium_values

    real(dp) function equilibrium_value_of(self, key, i) result(value)
        class(equilibrium_t), intent(in) :: self
        character(len=*), intent(in) :: key
        integer, intent(in), optional :: i
        value = self%vmec%get_input(key, i)
    end function equilibrium_value_of

    subroutine equilibrium_solve(self, x, ok, message)
        !! Solve at parameter values x; ok = .false. if VMEC++ found no equilibrium.
        class(equilibrium_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        logical, intent(out) :: ok
        character(len=:), allocatable, intent(out), optional :: message
        real(dp), allocatable :: tflux(:), pflux(:)
        integer :: k, start, finish, rate

        do k = 1, size(x)
            call self%vmec%set_input(trim(self%key(k)), x(k), self%i1(k), self%i2(k))
        end do
        call system_clock(start, rate)
        ok = self%vmec%solve(message)
        call system_clock(finish)
        self%t_solve = self%t_solve + real(finish - start, dp) / rate
        self%solves = self%solves + 1
        self%solved = ok
        if (.not. ok) return
        call self%vmec%geometry(self%coef, tflux, pflux)
        call self%vmec%radial(self%phip_full, self%phip_half, self%iota_half, self%current_half, &
            self%mass_half, self%lamscale)
        self%signgs = self%vmec%get_int('signgs')
        if (allocated(self%y)) deallocate(self%y)
        allocate(self%y(self%edge%n_y()))
        call self%edge%forward(self%coef, self%phip_full, self%phip_half, self%iota_half, &
            self%current_half, self%y)
    end subroutine equilibrium_solve

    integer function equilibrium_n_y(self)
        class(equilibrium_t), intent(in) :: self
        equilibrium_n_y = self%edge%n_y()
    end function equilibrium_n_y

    subroutine equilibrium_boundary(self, vb)
        !! The last solve's s = 1 data for Tiago's plasma model.
        class(equilibrium_t), intent(in) :: self
        type(vmec_boundary_t), intent(out) :: vb
        integer :: nb, ng

        nb = self%edge%mnmax_nyq
        ng = self%edge%mnmax
        vb%nfp = self%nfp
        vb%signgs = self%signgs
        vb%lasym = .false.
        vb%covariant = self%edge%covariant
        vb%conservative = self%conservative
        vb%phiedge = self%vmec%get_input('phiedge')
        vb%xm = self%edge%xm
        vb%xn = self%edge%xn
        vb%xm_nyq = self%edge%xm_nyq
        vb%xn_nyq = self%edge%xn_nyq
        vb%bumnc = self%y(1:nb)
        vb%bvmnc = self%y(nb + 1:2 * nb)
        vb%rmnc = self%y(2 * nb + 1:2 * nb + ng)
        vb%zmns = self%y(2 * nb + ng + 1:)
    end subroutine equilibrium_boundary

    subroutine equilibrium_jacobian(self, x, cot, jac)
        !! jac = cot (dy/dx) at the last solve (which must be at x); cot(k, n_y).
        class(equilibrium_t), intent(inout) :: self
        real(dp), intent(in) :: x(:), cot(:, :)
        real(dp), allocatable, intent(out) :: jac(:, :)
        real(dp), allocatable :: coef_bar(:, :, :, :), gbar(:, :), bbar(:, :), pbar(:, :)
        real(dp), allocatable :: iota_d(:, :), current_d(:, :), ps_col(:), cur_col(:), yb(:)
        real(dp), allocatable :: g_mass_ps(:), g_mass_am(:, :), g_iota_ai(:, :), g_cur(:, :)
        real(dp) :: pres_scale, curtor, phiedge
        integer :: r, k, nh, nrow, start, finish, rate, flip

        if (.not. self%solved) error stop 'equilibrium_jacobian: no equilibrium solved'
        nrow = size(cot, 1)
        nh = self%ns - 1
        allocate(coef_bar, mold=self%coef)
        allocate(gbar(size(self%coef), nrow), iota_d(nh, nrow), current_d(nh, nrow))
        call system_clock(start, rate)
        do r = 1, nrow
            call self%edge%vjp(self%coef, self%phip_full, self%phip_half, self%iota_half, &
                self%current_half, cot(r, :), coef_bar, iota_d(:, r), current_d(:, r))
            gbar(:, r) = reshape(coef_bar, [size(coef_bar)])
        end do
        call self%vmec%adjoint(gbar, bbar, pbar, self%max_block)
        call system_clock(finish)
        self%t_adjoint = self%t_adjoint + real(finish - start, dp) / rate
        self%adjoint_solves = self%adjoint_solves + nrow
        ! direct iota cotangent: VMEC++ reports the implicit one for the input
        ! iota (before a theta flip), the direct one is for the solver's iota
        flip = merge(-1, 1, self%vmec%get_int('have_to_flip_theta') == 1)
        pbar(nh + 1:2 * nh, :) = pbar(nh + 1:2 * nh, :) + flip * iota_d
        pbar(2 * nh + 1:3 * nh, :) = pbar(2 * nh + 1:3 * nh, :) + current_d

        call profile_gradients(self, g_mass_ps, g_mass_am, g_iota_ai, g_cur)
        pres_scale = self%vmec%get_input('pres_scale')
        curtor = self%vmec%get_input('curtor')
        allocate(jac(nrow, size(self%names)))
        jac = 0.0_dp
        ps_col = matmul(g_mass_ps, pbar(1:nh, :))
        cur_col = matmul(g_cur(:, 0), pbar(2 * nh + 1:3 * nh, :))
        do k = 1, size(self%names)
            select case (trim(self%key(k)))
            case ('rbc')
                jac(:, k) = bbar(bindex(0, self%i1(k), self%i2(k)), :)
            case ('zbs')
                jac(:, k) = bbar(bindex(1, self%i1(k), self%i2(k)), :)
            case ('pres_scale')
                jac(:, k) = ps_col
            case ('am')
                jac(:, k) = matmul(g_mass_am(:, self%i1(k) - 1), pbar(1:nh, :))
            case ('ai')
                if (self%ncurr == 1) error stop 'AI parameters need ncurr = 0'
                jac(:, k) = matmul(g_iota_ai(:, self%i1(k) - 1), pbar(nh + 1:2 * nh, :))
            case ('ac')
                if (self%ncurr /= 1) error stop 'AC parameters need ncurr = 1'
                jac(:, k) = matmul(g_cur(:, self%i1(k)), pbar(2 * nh + 1:3 * nh, :))
            case ('curtor')
                if (self%ncurr /= 1) error stop 'CURTOR needs ncurr = 1'
                jac(:, k) = cur_col
            case ('phiedge')
                phiedge = x(k)
                yb = self%y
                yb(2 * self%edge%mnmax_nyq + 1:) = 0.0_dp
                jac(:, k) = matmul(cot, yb) - 2.0_dp * pres_scale * ps_col
                if (self%ncurr == 1) jac(:, k) = jac(:, k) - curtor * cur_col
                jac(:, k) = jac(:, k) / phiedge
            end select
        end do
    contains
        integer function bindex(block, i1, i2)
            !! (block, m, n + ntor) in the adapter's layout; i1 = m + 1, i2 = n + ntor + 1
            integer, intent(in) :: block, i1, i2
            bindex = (block * self%mpol + (i1 - 1)) * (2 * self%ntor + 1) + i2
        end function bindex
    end subroutine equilibrium_jacobian

    subroutine profile_gradients(self, g_mass_ps, g_mass_am, g_iota_ai, g_cur)
        !! Derivatives of the half-grid profiles (VMEC++ RadialProfiles, power
        !! series, gamma = 0) with respect to the input coefficients:
        !!   mu0 p_h = mu0 PRES_SCALE sum_k AM(k) x_h^k,  x_h = min(|torflux(min(s_h, spres_ped)) bloat|, 1)
        !!   iota_h  = sum_k AI(k) torflux(s_h)^k
        !!   I_h     = itor x_h sum_k AC(k) x_h^k / (k + 1),  itor ~ CURTOR / (edge value)
        !! g_cur(:, 0) is d I_h / d CURTOR, g_cur(:, k + 1) is d I_h / d AC(k).
        class(equilibrium_t), intent(in) :: self
        real(dp), allocatable, intent(out) :: g_mass_ps(:), g_mass_am(:, :), g_iota_ai(:, :)
        real(dp), allocatable, intent(out) :: g_cur(:, :)
        real(dp), allocatable :: s(:), xm(:), xc(:), xi(:), am(:), ac(:)
        real(dp) :: bloat, spres_ped, pres_scale, curtor, e, integral_edge, itor
        integer :: nh, k, n_am, n_ac, n_ai, h

        nh = self%ns - 1
        s = [((h - 0.5_dp) / nh, h = 1, nh)]
        bloat = self%vmec%get_input('bloat')
        spres_ped = self%vmec%get_input('spres_ped')
        pres_scale = self%vmec%get_input('pres_scale')
        curtor = self%vmec%get_input('curtor')
        xm = min(abs(torflux(min(s, spres_ped)) * bloat), 1.0_dp)
        xi = torflux(s)
        xc = min(abs(torflux(s) * bloat), 1.0_dp)
        n_am = max(self%vmec%input_length('am'), 1) + 10
        n_ai = max(self%vmec%input_length('ai'), 1) + 10
        n_ac = self%vmec%input_length('ac')
        am = [(self%vmec%get_input('am', k), k = 1, n_am)]
        allocate(g_mass_am(nh, 0:n_am - 1), g_iota_ai(nh, 0:n_ai - 1), g_cur(nh, 0:n_ac + 10))
        g_mass_ps = mu0 * [(sum([(am(k + 1) * xm(h)**k, k = 0, n_am - 1)]), h = 1, nh)]
        do k = 0, n_am - 1
            g_mass_am(:, k) = mu0 * pres_scale * xm**k
        end do
        do k = 0, n_ai - 1
            g_iota_ai(:, k) = xi**k
        end do
        g_cur = 0.0_dp
        if (self%ncurr == 1 .and. n_ac > 0) then
            ac = [(self%vmec%get_input('ac', k), k = 1, n_ac + 10)]
            e = min(abs(bloat), 1.0_dp)
            integral_edge = e * sum([(ac(k + 1) * e**k / (k + 1), k = 0, size(ac) - 1)])
            if (abs(integral_edge) > abs(epsilon(1.0_dp) * curtor)) then
                ! I_h = itor * integral(x_h), itor = signgs mu0 curtor / (2 pi integral_edge)
                itor = self%current_half(nh) / integral(xc(nh))
                if (curtor /= 0.0_dp) g_cur(:, 0) = self%current_half / curtor
                do k = 0, size(ac) - 1
                    g_cur(:, k + 1) = itor * xc**(k + 1) / (k + 1) &
                        - self%current_half * e**(k + 1) / (k + 1) / integral_edge
                end do
            end if
        end if
    contains
        impure elemental real(dp) function torflux(x)
            real(dp), intent(in) :: x
            real(dp) :: t
            integer :: j, n
            n = self%vmec%input_length('aphi')
            if (n == 0) then
                t = x
            else
                t = 0.0_dp
                do j = 1, n
                    t = t + self%vmec%get_input('aphi', j) * x**j
                end do
            end if
            torflux = min(t, 1.0_dp)
        end function torflux
        real(dp) function integral(x)
            real(dp), intent(in) :: x
            integral = x * sum([(ac(k + 1) * x**k / (k + 1), k = 0, size(ac) - 1)])
        end function integral
    end subroutine profile_gradients

    subroutine equilibrium_write_wout(self, path)
        !! The wout fields Tiago and its plots use (netCDF, VMEC names).
        use netcdf, only: nf90_create, nf90_def_dim, nf90_def_var, nf90_enddef, nf90_put_var, &
            nf90_close, nf90_clobber, nf90_double, nf90_int
        class(equilibrium_t), intent(in) :: self
        character(len=*), intent(in) :: path
        character(len=8), parameter :: ints(4) = [character(len=8) :: 'nfp', 'signgs', 'ns', &
            'lasym']
        character(len=8), parameter :: modes(4) = [character(len=8) :: 'xm', 'xn', 'xm_nyq', &
            'xn_nyq']
        character(len=8), parameter :: radial(6) = [character(len=8) :: 'phi', 'presf', &
            'iotaf', 'jcurv', 'buco', 'chi']
        character(len=8), parameter :: fields(4) = [character(len=8) :: 'rmnc', 'zmns', &
            'bsupumnc', 'bsupvmnc']
        character(len=8), parameter :: scalars(4) = [character(len=8) :: 'ctor', 'b0', &
            'Rmajor_p', 'Aminor_p']
        integer :: ncid, dr, dm, dn, status, ns, mnmax, mnq, k, length
        integer :: id_int(4), id_mode(4), id_radial(6), id_field(4), id_scalar(4)

        ns = self%ns
        mnmax = self%vmec%get_int('mnmax')
        mnq = self%vmec%get_int('mnmax_nyq')
        status = nf90_create(path, nf90_clobber, ncid)
        if (status /= 0) error stop 'equilibrium_write_wout: cannot create '//path
        status = nf90_def_dim(ncid, 'radius', ns, dr)
        status = nf90_def_dim(ncid, 'mn_mode', mnmax, dm)
        status = nf90_def_dim(ncid, 'mn_mode_nyq', mnq, dn)
        do k = 1, 4
            status = nf90_def_var(ncid, trim(wout_name(ints(k))), nf90_int, id_int(k))
            status = nf90_def_var(ncid, trim(modes(k)), nf90_double, [merge(dm, dn, k <= 2)], &
                id_mode(k))
            status = nf90_def_var(ncid, trim(fields(k)), nf90_double, [merge(dm, dn, k <= 2), dr], &
                id_field(k))
            status = nf90_def_var(ncid, trim(scalars(k)), nf90_double, id_scalar(k))
        end do
        do k = 1, size(radial)
            status = nf90_def_var(ncid, trim(radial(k)), nf90_double, [dr], id_radial(k))
        end do
        status = nf90_enddef(ncid)
        status = nf90_put_var(ncid, id_int(1), self%nfp)
        status = nf90_put_var(ncid, id_int(2), self%signgs)
        status = nf90_put_var(ncid, id_int(3), ns)
        status = nf90_put_var(ncid, id_int(4), 0)
        do k = 1, 4
            length = merge(mnmax, mnq, k <= 2)
            status = nf90_put_var(ncid, id_mode(k), self%vmec%wout(trim(modes(k)), length))
            status = nf90_put_var(ncid, id_field(k), &
                reshape(self%vmec%wout(trim(fields(k)), ns * length), [length, ns]))
            status = nf90_put_var(ncid, id_scalar(k), scalar(self%vmec%wout(trim(scalars(k)), 1)))
        end do
        do k = 1, size(radial)
            status = nf90_put_var(ncid, id_radial(k), self%vmec%wout(trim(radial(k)), ns))
        end do
        status = nf90_close(ncid)
        if (status /= 0) error stop 'equilibrium_write_wout: cannot write '//path
    contains
        pure function wout_name(name) result(out)
            character(len=*), intent(in) :: name
            character(len=:), allocatable :: out
            out = trim(name)
            if (out == 'lasym') out = 'lasym__logical__'
        end function wout_name
        pure real(dp) function scalar(values)
            real(dp), intent(in) :: values(:)
            scalar = values(1)
        end function scalar
    end subroutine equilibrium_write_wout

    pure function lower(text) result(out)
        character(len=*), intent(in) :: text
        character(len=len(text)) :: out
        integer :: i
        out = text
        do i = 1, len(text)
            if (text(i:i) >= 'A' .and. text(i:i) <= 'Z') out(i:i) = achar(iachar(text(i:i)) + 32)
        end do
    end function lower
end module tiago_equilibrium
