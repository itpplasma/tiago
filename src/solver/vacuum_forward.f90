module tiago_vacuum_forward
    !! Forward model for flux loops, segmented Rogowskis and magnetic probes.
    !!
    !! Every diagnostic is reduced to quadrature points x_j with weighted line
    !! elements dl_j; all points are evaluated in one batch by the coil kernels
    !! (and the plasma sheet-current model if enabled):
    !!     flux loop   sum_j A(x_j) . dl_j          (closed polygon, DIAGNO iflflg)
    !!     Rogowski    sum_j eff_area B(x_j) . dl_j (open path)
    !!     probe       eff_area B(x) . n
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    use neo_biotsavart_field, only: biotsavart_field_t
    use tiago_coil_loader, only: load_coils_into_field
    use tiago_coil_kernels, only: coil_set_t, build_coil_set
    use tiago_diagnostic_types, only: flux_loop_t, segmented_rogowski_t, bprobe_t
    use tiago_plasma_support, only: plasma_support_t, vmec_boundary_t
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    type :: quadrature_rule_t
        !> Points per segment: midpoint rule (DIAGNO int_type='midpoint') or
        !> Gauss-Legendre, which converges much faster for smooth integrands.
        integer(i32) :: samples_per_segment = 6_i32
        logical :: gauss = .false.
    end type quadrature_rule_t

    type :: vacuum_solver_t
        type(coil_set_t) :: coils
        logical :: is_ready = .false.
        integer(i32) :: nfp = 1_i32
        type(plasma_support_t) :: plasma
        real(dp), allocatable :: x(:), y(:), z(:)  !! coil nodes [m] (response matrices)
        integer, allocatable :: group(:)           !! coil group of every node
        real(dp), allocatable :: unit_current(:)   !! current per unit EXTCUR [A]
        real(dp), allocatable :: extcur(:)         !! EXTCUR of every coil group
    contains
        procedure :: init => vacuum_solver_init
        procedure :: finalize => vacuum_solver_finalize
        procedure :: set_nfp => vacuum_solver_set_nfp
        procedure :: enable_plasma_from_vmec => vacuum_solver_enable_plasma_from_vmec
        procedure :: enable_plasma_from_boundary => vacuum_solver_enable_plasma_from_boundary
        procedure :: flux_loops => vacuum_solver_flux_loops
        procedure :: segrog => vacuum_solver_segrog
        procedure :: flux_and_segrog => vacuum_solver_flux_and_segrog
        procedure :: bprobes => vacuum_solver_bprobes
        procedure :: response => vacuum_solver_response
        procedure :: plasma_response => vacuum_solver_plasma_response
        procedure :: n_groups => vacuum_solver_n_groups
    end type vacuum_solver_t

    public :: vacuum_solver_t
    public :: quadrature_rule_t

contains

    subroutine vacuum_solver_init(self, coil_file, coil_extcur)
        !! An empty coil_file gives a coil-free (plasma-only) solver.
        class(vacuum_solver_t), intent(inout) :: self
        character(len=*), intent(in) :: coil_file
        character(len=*), intent(in), optional :: coil_extcur
        type(biotsavart_field_t) :: field
        integer :: g

        if (len_trim(coil_file) == 0) then
            allocate(self%x(0), self%y(0), self%z(0), self%group(0), self%unit_current(0))
            allocate(self%extcur(0))
            self%coils = build_coil_set(self%x, self%y, self%z, self%unit_current)
            self%is_ready = .true.
            return
        end if

        if (present(coil_extcur)) then
            call load_coils_into_field(field, trim(coil_file), trim(coil_extcur), &
                self%group, self%unit_current)
        else
            call load_coils_into_field(field, trim(coil_file), groups=self%group, &
                unit_current=self%unit_current)
        end if
        self%x = field%coils%x
        self%y = field%coils%y
        self%z = field%coils%z
        self%coils = build_coil_set(self%x, self%y, self%z, field%coils%current)
        allocate(self%extcur(self%n_groups()))
        self%extcur = 0.0_dp
        do g = size(self%group), 1, -1
            if (self%group(g) < 1 .or. self%unit_current(g) == 0.0_dp) cycle
            self%extcur(self%group(g)) = field%coils%current(g) / self%unit_current(g)
        end do
        self%is_ready = .true.
    end subroutine vacuum_solver_init

    subroutine vacuum_solver_set_nfp(self, value)
        class(vacuum_solver_t), intent(inout) :: self
        integer(i32), intent(in) :: value

        if (value > 0_i32) self%nfp = value
    end subroutine vacuum_solver_set_nfp

    subroutine vacuum_solver_finalize(self)
        class(vacuum_solver_t), intent(inout) :: self
        self%is_ready = .false.
        call self%plasma%finalize()
    end subroutine vacuum_solver_finalize

    subroutine vacuum_solver_enable_plasma_from_vmec(self, wout_file, nphi, ntheta, &
            covariant, conservative)
        class(vacuum_solver_t), intent(inout) :: self
        character(len=*), intent(in) :: wout_file
        integer(i32), intent(in) :: nphi
        integer(i32), intent(in) :: ntheta

        logical, intent(in), optional :: covariant, conservative

        call self%plasma%init_from_vmec(trim(wout_file), nphi, ntheta, &
            covariant, conservative)
    end subroutine vacuum_solver_enable_plasma_from_vmec

    subroutine vacuum_solver_enable_plasma_from_boundary(self, vb, nphi, ntheta)
        !! Plasma model from in-memory boundary data (e.g. a VMEC++ solve).
        class(vacuum_solver_t), intent(inout) :: self
        type(vmec_boundary_t), intent(in) :: vb
        integer(i32), intent(in) :: nphi
        integer(i32), intent(in) :: ntheta

        call self%plasma%init_from_boundary(vb, nphi, ntheta)
    end subroutine vacuum_solver_enable_plasma_from_boundary

    integer function vacuum_solver_n_groups(self)
        class(vacuum_solver_t), intent(in) :: self
        vacuum_solver_n_groups = 0
        if (allocated(self%group)) then
            if (size(self%group) > 0) vacuum_solver_n_groups = maxval(self%group)
        end if
    end function vacuum_solver_n_groups

    subroutine vacuum_solver_flux_loops(self, loops, fluxes, rule)
        class(vacuum_solver_t), intent(in) :: self
        type(flux_loop_t), allocatable, intent(in) :: loops(:)
        real(dp), allocatable, intent(out) :: fluxes(:)
        type(quadrature_rule_t), intent(in), optional :: rule

        call assert_ready(self)
        if (.not. allocated(loops)) error stop 'flux loop array not set'
        fluxes = loop_fluxes(self, self%coils, loops, rule_or_default(rule), .true.)
    end subroutine vacuum_solver_flux_loops

    subroutine vacuum_solver_segrog(self, diagnostics, voltages, rule)
        class(vacuum_solver_t), intent(in) :: self
        type(segmented_rogowski_t), allocatable, intent(in) :: diagnostics(:)
        real(dp), allocatable, intent(out) :: voltages(:)
        type(quadrature_rule_t), intent(in), optional :: rule

        call assert_ready(self)
        if (.not. allocated(diagnostics)) error stop 'segmented Rogowski array not set'
        voltages = segrog_signals(self, self%coils, diagnostics, rule_or_default(rule), .true.)
    end subroutine vacuum_solver_segrog

    subroutine vacuum_solver_flux_and_segrog(self, loops, fluxes, diagnostics, &
            voltages, rule)
        class(vacuum_solver_t), intent(in) :: self
        type(flux_loop_t), allocatable, intent(in) :: loops(:)
        real(dp), allocatable, intent(out) :: fluxes(:)
        type(segmented_rogowski_t), allocatable, intent(in) :: diagnostics(:)
        real(dp), allocatable, intent(out) :: voltages(:)
        type(quadrature_rule_t), intent(in), optional :: rule

        call self%flux_loops(loops, fluxes, rule)
        call self%segrog(diagnostics, voltages, rule)
    end subroutine vacuum_solver_flux_and_segrog

    subroutine vacuum_solver_bprobes(self, probes, signals)
        !! eff_area * B . normal at every probe (coils + plasma), before turns.
        class(vacuum_solver_t), intent(in) :: self
        type(bprobe_t), intent(in) :: probes(:)
        real(dp), allocatable, intent(out) :: signals(:)

        call assert_ready(self)
        signals = probe_signals(self, self%coils, probes, .true.)
    end subroutine vacuum_solver_bprobes

    subroutine vacuum_solver_response(self, loops, segs, probes, flux_resp, seg_resp, &
            probe_resp, rule)
        !! Coil signals per unit EXTCUR of each coil group, column g = group g
        !! (DIAGNO -mutual). Vacuum only: the plasma part is not linear in EXTCUR.
        !! Before turns/idia post-processing, which the caller applies per column.
        class(vacuum_solver_t), intent(in) :: self
        type(flux_loop_t), allocatable, intent(in) :: loops(:)
        type(segmented_rogowski_t), allocatable, intent(in) :: segs(:)
        type(bprobe_t), allocatable, intent(in) :: probes(:)
        real(dp), allocatable, intent(out) :: flux_resp(:, :), seg_resp(:, :), probe_resp(:, :)
        type(quadrature_rule_t), intent(in), optional :: rule

        type(coil_set_t) :: group_coils
        type(quadrature_rule_t) :: q
        logical, allocatable :: in_group(:)
        integer :: g, ng

        call assert_ready(self)
        q = rule_or_default(rule)
        ng = self%n_groups()
        allocate(flux_resp(0, ng), seg_resp(0, ng), probe_resp(0, ng))
        if (allocated(loops)) then
            deallocate(flux_resp)
            allocate(flux_resp(size(loops), ng))
        end if
        if (allocated(segs)) then
            deallocate(seg_resp)
            allocate(seg_resp(size(segs), ng))
        end if
        if (allocated(probes)) then
            deallocate(probe_resp)
            allocate(probe_resp(size(probes), ng))
        end if
        do g = 1, ng
            ! A group's coils end with a zero-current node, so dropping the other
            ! groups' nodes cannot join coils; the cost is one evaluation in total.
            in_group = self%group == g
            group_coils = build_coil_set(pack(self%x, in_group), pack(self%y, in_group), &
                pack(self%z, in_group), pack(self%unit_current, in_group))
            if (allocated(loops)) flux_resp(:, g) = loop_fluxes(self, group_coils, loops, q, .false.)
            if (allocated(segs)) seg_resp(:, g) = segrog_signals(self, group_coils, segs, q, .false.)
            if (allocated(probes)) probe_resp(:, g) = probe_signals(self, group_coils, probes, .false.)
        end do
    end subroutine vacuum_solver_response

    subroutine vacuum_solver_plasma_response(self, loops, segs, probes, flux_resp, seg_resp, &
            probe_resp, rule, flux_shape, seg_shape, probe_shape)
        !! Derivatives of the plasma part of every signal with respect to the VMEC
        !! boundary field coefficients (columns: plasma%mode_column_name). Exact,
        !! since the plasma part is linear in them at fixed boundary shape.
        !! Optionally also with respect to the boundary geometry coefficients
        !! (*_shape, columns: plasma%shape_column_name) at fixed field coefficients.
        class(vacuum_solver_t), intent(in) :: self
        type(flux_loop_t), allocatable, intent(in) :: loops(:)
        type(segmented_rogowski_t), allocatable, intent(in) :: segs(:)
        type(bprobe_t), allocatable, intent(in) :: probes(:)
        real(dp), allocatable, intent(out) :: flux_resp(:, :), seg_resp(:, :), probe_resp(:, :)
        type(quadrature_rule_t), intent(in), optional :: rule
        real(dp), allocatable, intent(out), optional :: flux_shape(:, :), seg_shape(:, :), &
            probe_shape(:, :)
        real(dp), allocatable :: points(:, :), dls(:, :), w(:, :, :), gx(:, :, :)
        integer, allocatable :: owner(:)
        integer :: j, nc, ng
        logical :: shape

        if (.not. self%plasma%has_data()) error stop 'plasma response needs --plasma-wout'
        shape = present(flux_shape)
        nc = self%plasma%n_mode_columns()
        ng = self%plasma%n_shape_columns()
        allocate(flux_resp(0, nc), seg_resp(0, nc), probe_resp(0, nc))
        if (shape) allocate(flux_shape(0, ng), seg_shape(0, ng), probe_shape(0, ng))
        if (allocated(loops)) then
            call sample_loops(loops, self%nfp, rule_or_default(rule), points, dls, owner)
            if (shape) then
                call self%plasma%potential_weights(points, dls, owner, size(loops), w, gx)
                call self%plasma%shape_response(w, gx, flux_shape)
            else
                call self%plasma%potential_weights(points, dls, owner, size(loops), w)
            end if
            call self%plasma%mode_response(w, flux_resp)
            do j = 1, size(loops)
                if (.not. loops(j)%one_period) cycle
                flux_resp(j, :) = flux_resp(j, :) * real(self%nfp, dp)
                if (shape) flux_shape(j, :) = flux_shape(j, :) * real(self%nfp, dp)
            end do
        end if
        if (allocated(segs)) then
            call sample_segrogs(segs, rule_or_default(rule), points, dls, owner)
            if (shape) then
                call self%plasma%field_weights(points, dls, owner, size(segs), w, gx)
                call self%plasma%shape_response(w, gx, seg_shape)
            else
                call self%plasma%field_weights(points, dls, owner, size(segs), w)
            end if
            call self%plasma%mode_response(w, seg_resp)
        end if
        if (allocated(probes)) then
            if (allocated(points)) deallocate(points, dls, owner)
            allocate(points(3, size(probes)), dls(3, size(probes)), owner(size(probes)))
            do j = 1, size(probes)
                points(:, j) = probes(j)%position
                dls(:, j) = probes(j)%eff_area * probes(j)%normal
                owner(j) = j
            end do
            if (shape) then
                call self%plasma%field_weights(points, dls, owner, size(probes), w, gx)
                call self%plasma%shape_response(w, gx, probe_shape)
            else
                call self%plasma%field_weights(points, dls, owner, size(probes), w)
            end if
            call self%plasma%mode_response(w, probe_resp)
        end if
    end subroutine vacuum_solver_plasma_response

    function loop_fluxes(self, coils, loops, rule, with_plasma) result(fluxes)
        class(vacuum_solver_t), intent(in) :: self
        type(coil_set_t), intent(in) :: coils
        type(flux_loop_t), intent(in) :: loops(:)
        type(quadrature_rule_t), intent(in) :: rule
        logical, intent(in) :: with_plasma
        real(dp) :: fluxes(size(loops))
        real(dp), allocatable :: points(:, :), dls(:, :), a(:, :), ap(:, :)
        integer, allocatable :: owner(:)
        integer :: j

        call sample_loops(loops, self%nfp, rule, points, dls, owner)
        allocate(a(3, size(owner)))
        call coils%vector_potential(points, a)
        if (with_plasma .and. self%plasma%has_data()) then
            allocate(ap(3, size(owner)))
            call self%plasma%sample_vector_potential(points, ap)
            a = a + ap
        end if
        fluxes = 0.0_dp
        do j = 1, size(owner)
            fluxes(owner(j)) = fluxes(owner(j)) + dot_product(a(:, j), dls(:, j))
        end do
        do j = 1, size(loops)
            if (loops(j)%one_period) fluxes(j) = fluxes(j) * real(self%nfp, dp)
        end do
    end function loop_fluxes

    function segrog_signals(self, coils, segs, rule, with_plasma) result(signals)
        class(vacuum_solver_t), intent(in) :: self
        type(coil_set_t), intent(in) :: coils
        type(segmented_rogowski_t), intent(in) :: segs(:)
        type(quadrature_rule_t), intent(in) :: rule
        logical, intent(in) :: with_plasma
        real(dp) :: signals(size(segs))
        real(dp), allocatable :: points(:, :), dls(:, :), b(:, :), bp(:, :)
        integer, allocatable :: owner(:)
        integer :: j

        call sample_segrogs(segs, rule, points, dls, owner)
        allocate(b(3, size(owner)))
        call coils%field(points, b)
        if (with_plasma .and. self%plasma%has_data()) then
            allocate(bp(3, size(owner)))
            call self%plasma%sample_bfield(points, bp)
            b = b + bp
        end if
        signals = 0.0_dp
        do j = 1, size(owner)
            signals(owner(j)) = signals(owner(j)) + dot_product(b(:, j), dls(:, j))
        end do
    end function segrog_signals

    function probe_signals(self, coils, probes, with_plasma) result(signals)
        class(vacuum_solver_t), intent(in) :: self
        type(coil_set_t), intent(in) :: coils
        type(bprobe_t), intent(in) :: probes(:)
        logical, intent(in) :: with_plasma
        real(dp) :: signals(size(probes))
        real(dp) :: points(3, size(probes)), b(3, size(probes)), bp(3, size(probes))
        integer :: j

        do j = 1, size(probes)
            points(:, j) = probes(j)%position
        end do
        call coils%field(points, b)
        if (with_plasma .and. self%plasma%has_data()) then
            call self%plasma%sample_bfield(points, bp)
            b = b + bp
        end if
        do j = 1, size(probes)
            signals(j) = probes(j)%eff_area * dot_product(b(:, j), probes(j)%normal)
        end do
    end function probe_signals

    subroutine sample_loops(loops, nfp, rule, points, dls, owner)
        !! Quadrature points of closed flux loops. DIAGNO semantics: the last
        !! point connects back to the first, or for iflflg=1 to the first point
        !! rotated by one field period.
        type(flux_loop_t), intent(in) :: loops(:)
        integer(i32), intent(in) :: nfp
        type(quadrature_rule_t), intent(in) :: rule
        real(dp), allocatable, intent(out) :: points(:, :), dls(:, :)
        integer, allocatable, intent(out) :: owner(:)
        real(dp), allocatable :: t(:), w(:)
        real(dp) :: a(3), b(3)
        integer :: i, seg, k, n, m

        call segment_rule(rule, t, w)
        m = size(t)
        n = 0
        do i = 1, size(loops)
            n = n + size(loops(i)%points) * m
        end do
        allocate(points(3, n), dls(3, n), owner(n))
        n = 0
        do i = 1, size(loops)
            do seg = 1, size(loops(i)%points)
                a = xyz(loops(i), seg)
                if (seg < size(loops(i)%points)) then
                    b = xyz(loops(i), seg + 1)
                else if (loops(i)%one_period) then
                    b = rotate_z(xyz(loops(i), 1), 2.0_dp * pi / real(nfp, dp))
                else
                    b = xyz(loops(i), 1)
                end if
                do k = 1, m
                    n = n + 1
                    points(:, n) = a + t(k) * (b - a)
                    dls(:, n) = w(k) * (b - a)
                    owner(n) = i
                end do
            end do
        end do
    end subroutine sample_loops

    subroutine sample_segrogs(segs, rule, points, dls, owner)
        !! Quadrature points of open Rogowski paths; line elements carry eff_area.
        type(segmented_rogowski_t), intent(in) :: segs(:)
        type(quadrature_rule_t), intent(in) :: rule
        real(dp), allocatable, intent(out) :: points(:, :), dls(:, :)
        integer, allocatable, intent(out) :: owner(:)
        real(dp), allocatable :: t(:), w(:)
        real(dp) :: a(3), b(3)
        integer :: i, seg, k, n, m

        call segment_rule(rule, t, w)
        m = size(t)
        n = 0
        do i = 1, size(segs)
            n = n + (size(segs(i)%path) - 1) * m
        end do
        allocate(points(3, n), dls(3, n), owner(n))
        n = 0
        do i = 1, size(segs)
            do seg = 1, size(segs(i)%path) - 1
                a = [segs(i)%path(seg)%x, segs(i)%path(seg)%y, segs(i)%path(seg)%z]
                b = [segs(i)%path(seg + 1)%x, segs(i)%path(seg + 1)%y, segs(i)%path(seg + 1)%z]
                do k = 1, m
                    n = n + 1
                    points(:, n) = a + t(k) * (b - a)
                    dls(:, n) = w(k) * segs(i)%segment_area(seg) * (b - a)
                    owner(n) = i
                end do
            end do
        end do
    end subroutine sample_segrogs

    subroutine segment_rule(rule, t, w)
        !! Nodes t in (0, 1) and weights w (sum 1) of the per-segment rule.
        type(quadrature_rule_t), intent(in) :: rule
        real(dp), allocatable, intent(out) :: t(:), w(:)
        integer :: n, k

        n = max(1_i32, rule%samples_per_segment)
        allocate(t(n), w(n))
        if (rule%gauss) then
            call gauss_legendre(n, t, w)
        else
            t = [((real(k, dp) - 0.5_dp) / real(n, dp), k = 1, n)]
            w = 1.0_dp / real(n, dp)
        end if
    end subroutine segment_rule

    subroutine gauss_legendre(n, t, w)
        !! n-point Gauss-Legendre rule mapped to [0, 1] (Newton on P_n).
        integer, intent(in) :: n
        real(dp), intent(out) :: t(n), w(n)
        real(dp) :: x, p0, p1, p2, dp_dx
        integer :: i, k, iter

        do i = 1, n
            x = cos(pi * (real(i, dp) - 0.25_dp) / (real(n, dp) + 0.5_dp))
            do iter = 1, 100
                p0 = 1.0_dp
                p1 = x
                do k = 2, n
                    p2 = ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
                    p0 = p1
                    p1 = p2
                end do
                dp_dx = n * (x * p1 - p0) / (x * x - 1.0_dp)
                x = x - p1 / dp_dx
                if (abs(p1 / dp_dx) < 1.0e-15_dp) exit
            end do
            t(i) = 0.5_dp * (1.0_dp - x)
            w(i) = 1.0_dp / ((1.0_dp - x * x) * dp_dx * dp_dx)
        end do
    end subroutine gauss_legendre

    function xyz(loop, i) result(p)
        type(flux_loop_t), intent(in) :: loop
        integer, intent(in) :: i
        real(dp) :: p(3)
        p = [loop%points(i)%x, loop%points(i)%y, loop%points(i)%z]
    end function xyz

    pure function rotate_z(p, angle) result(q)
        real(dp), intent(in) :: p(3), angle
        real(dp) :: q(3)
        q = [p(1) * cos(angle) - p(2) * sin(angle), p(1) * sin(angle) + p(2) * cos(angle), p(3)]
    end function rotate_z

    type(quadrature_rule_t) function rule_or_default(rule)
        type(quadrature_rule_t), intent(in), optional :: rule
        if (present(rule)) rule_or_default = rule
    end function rule_or_default

    subroutine assert_ready(self)
        class(vacuum_solver_t), intent(in) :: self
        if (.not. self%is_ready) error stop 'vacuum solver not initialised'
    end subroutine assert_ready
end module tiago_vacuum_forward
