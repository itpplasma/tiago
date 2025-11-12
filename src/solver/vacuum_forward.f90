module tiago_vacuum_forward
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    use neo_biotsavart_field, only: biotsavart_field_t
    use tiago_coil_loader, only: load_coils_into_field
    use tiago_diagnostic_types, only: flux_loop_t, segmented_rogowski_t, &
        loop_point_t
    implicit none
    private

    real(dp), parameter :: closure_tolerance = 1.0e-10_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: two_pi = 2.0_dp * pi
    real(dp), parameter :: meters_to_cm = 100.0_dp
    real(dp), parameter :: amps_to_statamp = 2.9979245368431e9_dp
    real(dp), parameter :: gauss_to_tesla = 1.0e-4_dp
    real(dp), parameter :: maxwell_to_weber = 1.0e-8_dp

    type :: quadrature_rule_t
        integer(i32) :: samples_per_segment = 4_i32
    end type quadrature_rule_t

    type :: quadrature_override_t
        character(len=:), allocatable :: label
        type(quadrature_rule_t) :: rule
    end type quadrature_override_t

    type :: vacuum_solver_t
        type(biotsavart_field_t) :: field
        logical :: is_ready = .false.
        integer(i32) :: nfp = 1_i32
    contains
        procedure :: init => vacuum_solver_init
        procedure :: finalize => vacuum_solver_finalize
        procedure :: flux_loops => vacuum_solver_flux_loops
        procedure :: segrog => vacuum_solver_segrog
        procedure :: set_nfp => vacuum_solver_set_nfp
        procedure :: flux_and_segrog => vacuum_solver_flux_and_segrog
    end type vacuum_solver_t

    public :: vacuum_solver_t
    public :: quadrature_rule_t
    public :: quadrature_override_t

contains
    subroutine vacuum_solver_init(self, coil_file, coil_extcur)
        class(vacuum_solver_t), intent(inout) :: self
        character(len=*), intent(in) :: coil_file
        character(len=*), intent(in), optional :: coil_extcur

        if (present(coil_extcur)) then
            if (len_trim(coil_extcur) > 0) then
                call load_coils_into_field(self%field, trim(coil_file), trim(coil_extcur))
            else
                call load_coils_into_field(self%field, trim(coil_file))
            end if
        else
            call load_coils_into_field(self%field, trim(coil_file))
        end if

        call scale_coils_to_cgs(self%field)
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
    end subroutine vacuum_solver_finalize

    subroutine vacuum_solver_flux_loops(self, loops, fluxes, default_rule, &
            overrides)
        class(vacuum_solver_t), intent(in) :: self
        type(flux_loop_t), allocatable, intent(in) :: loops(:)
        real(dp), allocatable, intent(out) :: fluxes(:)
        type(quadrature_rule_t), intent(in), optional :: default_rule
        type(quadrature_override_t), intent(in), optional :: overrides(:)

        integer :: i
        type(quadrature_rule_t) :: rule

        call assert_ready(self)
        if (.not. allocated(loops)) call abort_with('flux loop array not set')
        if (allocated(fluxes)) deallocate(fluxes)
        allocate(fluxes(size(loops)))

!$omp parallel do default(shared) private(i, rule) collapse(1)
        do i = 1, size(loops)
            rule = select_rule(loops(i)%label, default_rule, overrides)
            fluxes(i) = evaluate_loop_flux(self%field, loops(i), rule, &
                self%nfp)
        end do
!$omp end parallel do
    end subroutine vacuum_solver_flux_loops

    subroutine vacuum_solver_segrog(self, diagnostics, voltages, default_rule, &
            overrides)
        class(vacuum_solver_t), intent(in) :: self
        type(segmented_rogowski_t), allocatable, intent(in) :: diagnostics(:)
        real(dp), allocatable, intent(out) :: voltages(:)
        type(quadrature_rule_t), intent(in), optional :: default_rule
        type(quadrature_override_t), intent(in), optional :: overrides(:)

        integer :: i
        type(quadrature_rule_t) :: rule

        call assert_ready(self)
        if (.not. allocated(diagnostics)) then
            call abort_with('segmented Rogowski array not set')
        end if
        if (allocated(voltages)) deallocate(voltages)
        allocate(voltages(size(diagnostics)))

!$omp parallel do default(shared) private(i, rule)
        do i = 1, size(diagnostics)
            rule = select_rule(diagnostics(i)%label, default_rule, overrides)
            voltages(i) = evaluate_segrog_signal(self%field, diagnostics(i), &
                rule)
        end do
!$omp end parallel do
    end subroutine vacuum_solver_segrog

    subroutine vacuum_solver_flux_and_segrog(self, loops, fluxes, diagnostics, &
            voltages, default_rule, overrides)
        class(vacuum_solver_t), intent(in) :: self
        type(flux_loop_t), allocatable, intent(in) :: loops(:)
        real(dp), allocatable, intent(out) :: fluxes(:)
        type(segmented_rogowski_t), allocatable, intent(in) :: diagnostics(:)
        real(dp), allocatable, intent(out) :: voltages(:)
        type(quadrature_rule_t), intent(in), optional :: default_rule
        type(quadrature_override_t), intent(in), optional :: overrides(:)

        integer :: i
        type(quadrature_rule_t) :: rule

        call assert_ready(self)
        if (.not. allocated(loops)) call abort_with('flux loop array not set')
        if (.not. allocated(diagnostics)) then
            call abort_with('segmented Rogowski array not set')
        end if
        if (allocated(fluxes)) deallocate(fluxes)
        if (allocated(voltages)) deallocate(voltages)
        allocate(fluxes(size(loops)))
        allocate(voltages(size(diagnostics)))

!$omp parallel default(shared) private(i, rule)
!$omp do nowait
        do i = 1, size(loops)
            rule = select_rule(loops(i)%label, default_rule, overrides)
            fluxes(i) = evaluate_loop_flux(self%field, loops(i), rule, &
                self%nfp)
        end do
!$omp end do
!$omp do
        do i = 1, size(diagnostics)
            rule = select_rule(diagnostics(i)%label, default_rule, overrides)
            voltages(i) = evaluate_segrog_signal(self%field, diagnostics(i), &
                rule)
        end do
!$omp end do
!$omp end parallel
    end subroutine vacuum_solver_flux_and_segrog

    subroutine assert_ready(self)
        class(vacuum_solver_t), intent(in) :: self
        if (.not. self%is_ready) call abort_with('vacuum solver not initialised')
    end subroutine assert_ready

    function evaluate_loop_flux(field, loop, rule, nfp) result(flux)
        type(biotsavart_field_t), intent(in) :: field
        type(flux_loop_t), intent(in) :: loop
        type(quadrature_rule_t), intent(in) :: rule
        integer(i32), intent(in) :: nfp
        real(dp) :: flux

        integer :: seg
        integer :: samples
        integer :: period
        real(dp) :: angle
        real(dp) :: dl(3)
        real(dp) :: start_point(3)
        real(dp) :: end_point(3)
        real(dp) :: weight
        real(dp) :: rotated_start(3)
        real(dp) :: rotated_end(3)
        real(dp) :: rotated_dl(3)
        real(dp) :: seg_flux

        flux = 0.0_dp
        samples = max(1_i32, rule%samples_per_segment)

        do seg = 1, segment_count(loop)
            call segment_endpoints(loop, seg, start_point, end_point)
            dl = end_point - start_point
            weight = 1.0_dp / real(samples, dp)
            seg_flux = integrate_segment(field, start_point, dl, samples, weight)
            flux = flux + seg_flux
            if (loop%repeat_count > 0 .and. nfp > 1) then
                do period = 1, nfp - 1
                    angle = real(period, dp) * two_pi / real(nfp, dp)
                    call rotate_point(start_point, angle, rotated_start)
                    call rotate_point(end_point, angle, rotated_end)
                    rotated_dl = rotated_end - rotated_start
                    flux = flux + integrate_segment(field, rotated_start, &
                        rotated_dl, samples, weight)
                end do
            end if
        end do

        flux = flux * maxwell_to_weber

        if (loop%subtract_toroidal_flux) then
            flux = flux - estimate_toroidal_flux(field, loop)
            if (loop%repeat_count > 0 .and. nfp > 1) then
                do period = 1, nfp - 1
                    angle = real(period, dp) * two_pi / real(nfp, dp)
                    call rotate_point(loop_centroid(loop), angle, rotated_start)
                    flux = flux - toroidal_flux_at_point(field, rotated_start, &
                        polygon_area_xy(loop))
                end do
            end if
        end if
    end function evaluate_loop_flux

    real(dp) function integrate_segment(field, start_point, dl, samples, weight)
        type(biotsavart_field_t), intent(in) :: field
        real(dp), intent(in) :: start_point(3)
        real(dp), intent(in) :: dl(3)
        integer(i32), intent(in) :: samples
        real(dp), intent(in) :: weight

        integer :: s
        real(dp) :: a_field(3)
        real(dp) :: step
        real(dp) :: sample_point_cm(3)
        real(dp) :: dl_cm(3)
        real(dp) :: start_point_cm(3)

        integrate_segment = 0.0_dp
        start_point_cm = start_point * meters_to_cm
        dl_cm = dl * meters_to_cm
        ! Avoid recomputing conversion: sample_point_cm = (start_point + step*dl) * meters_to_cm
        !                             = start_point_cm + step * dl_cm
        do s = 1, samples
            step = (real(s, dp) - 0.5_dp) / real(samples, dp)
            sample_point_cm = start_point_cm + step * dl_cm
            call field%compute_afield(sample_point_cm, a_field)
            integrate_segment = integrate_segment + &
                weight * dot_product(a_field, dl_cm)
        end do
    end function integrate_segment

    real(dp) function estimate_toroidal_flux(field, loop)
        type(biotsavart_field_t), intent(in) :: field
        type(flux_loop_t), intent(in) :: loop

        real(dp) :: centroid(3)
        real(dp) :: area

        centroid = loop_centroid(loop)
        area = polygon_area_xy(loop)
        estimate_toroidal_flux = toroidal_flux_at_point(field, centroid, area)
    end function estimate_toroidal_flux

    real(dp) function toroidal_flux_at_point(field, position, area)
        type(biotsavart_field_t), intent(in) :: field
        real(dp), intent(in) :: position(3)
        real(dp), intent(in) :: area

        real(dp) :: position_cm(3)
        real(dp) :: b_gauss(3)

        position_cm = position * meters_to_cm
        call field%compute_bfield(position_cm, b_gauss)
        toroidal_flux_at_point = b_gauss(3) * gauss_to_tesla * area
    end function toroidal_flux_at_point

    function loop_centroid(loop) result(center)
        type(flux_loop_t), intent(in) :: loop
        real(dp) :: center(3)
        integer :: i

        center = 0.0_dp
        do i = 1, size(loop%points)
            center(1) = center(1) + loop%points(i)%x
            center(2) = center(2) + loop%points(i)%y
            center(3) = center(3) + loop%points(i)%z
        end do
        center = center / real(size(loop%points), dp)
    end function loop_centroid

    real(dp) function polygon_area_xy(loop)
        type(flux_loop_t), intent(in) :: loop
        integer :: npts
        integer :: i
        real(dp) :: x1
        real(dp) :: y1
        real(dp) :: x2
        real(dp) :: y2

        npts = size(loop%points)
        polygon_area_xy = 0.0_dp
        do i = 1, npts
            x1 = loop%points(i)%x
            y1 = loop%points(i)%y
            x2 = loop%points(next_index(loop, i))%x
            y2 = loop%points(next_index(loop, i))%y
            polygon_area_xy = polygon_area_xy + (x1 * y2 - x2 * y1)
        end do
        polygon_area_xy = 0.5_dp * polygon_area_xy
    end function polygon_area_xy

    function evaluate_segrog_signal(field, diagnostic, rule) result(voltage)
        type(biotsavart_field_t), intent(in) :: field
        type(segmented_rogowski_t), intent(in) :: diagnostic
        type(quadrature_rule_t), intent(in) :: rule
        real(dp) :: voltage

        integer :: seg
        integer :: samples
        real(dp) :: start_point(3)
        real(dp) :: end_point(3)
        real(dp) :: dl(3)
        real(dp) :: tangent(3)
        real(dp) :: norm_dl
        real(dp) :: weight
        integer(i32) :: effective_segments

        samples = max(1_i32, rule%samples_per_segment)
        weight = 1.0_dp / real(samples, dp)
        voltage = 0.0_dp
        effective_segments = max(1_i32, diagnostic%segments)

        do seg = 1, size(diagnostic%path) - 1
            call extract_path_segment(diagnostic, seg, start_point, end_point)
            dl = end_point - start_point
            norm_dl = max(closure_tolerance, sqrt(sum(dl**2)))
            tangent = dl / norm_dl
            voltage = voltage + integrate_segrog_segment(field, start_point, dl, &
                norm_dl, tangent, samples, weight)
        end do

        voltage = voltage * diagnostic%effective_area / &
            real(effective_segments, dp)
    end function evaluate_segrog_signal

    real(dp) function integrate_segrog_segment(field, start_point, dl, &
            norm_dl, tangent, samples, weight)
        type(biotsavart_field_t), intent(in) :: field
        real(dp), intent(in) :: start_point(3)
        real(dp), intent(in) :: dl(3)
        real(dp), intent(in) :: norm_dl
        real(dp), intent(in) :: tangent(3)
        integer(i32), intent(in) :: samples
        real(dp), intent(in) :: weight

        integer :: s
        real(dp) :: step
        real(dp) :: sample_point(3)
        real(dp) :: b_field(3)
        real(dp) :: sample_point_cm(3)
        real(dp) :: b_gauss(3)

        integrate_segrog_segment = 0.0_dp
        do s = 1, samples
            step = (real(s, dp) - 0.5_dp) / real(samples, dp)
            sample_point = start_point + step * dl
            sample_point_cm = sample_point * meters_to_cm
            call field%compute_bfield(sample_point_cm, b_gauss)
            b_field = b_gauss * gauss_to_tesla
            integrate_segrog_segment = integrate_segrog_segment + weight * &
                dot_product(b_field, tangent) * norm_dl
        end do
    end function integrate_segrog_segment

    type(quadrature_rule_t) function select_rule(label, default_rule, overrides)
        character(len=*), intent(in) :: label
        type(quadrature_rule_t), intent(in), optional :: default_rule
        type(quadrature_override_t), intent(in), optional :: overrides(:)

        integer :: i

        if (present(overrides)) then
            do i = 1, size(overrides)
                if (labels_equal(label, overrides(i)%label)) then
                    select_rule = overrides(i)%rule
                    return
                end if
            end do
        end if

        if (present(default_rule)) then
            select_rule = default_rule
        else
            select_rule%samples_per_segment = 4_i32
        end if
    end function select_rule

    logical function labels_equal(a, b)
        character(len=*), intent(in) :: a
        character(len=*), intent(in) :: b

        labels_equal = trim(a) == trim(b)
    end function labels_equal

    integer function segment_count(loop)
        type(flux_loop_t), intent(in) :: loop
        integer :: npts

        npts = size(loop%points)
        if (loop%is_open) then
            segment_count = max(0, npts - 1)
        else if (points_match(loop%points(1), loop%points(npts))) then
            segment_count = max(1, npts - 1)
        else
            segment_count = npts
        end if
    end function segment_count

    subroutine segment_endpoints(loop, index, start_point, end_point)
        type(flux_loop_t), intent(in) :: loop
        integer, intent(in) :: index
        real(dp), intent(out) :: start_point(3)
        real(dp), intent(out) :: end_point(3)

        type(loop_point_t) :: p1
        type(loop_point_t) :: p2
        integer :: start_idx
        integer :: end_idx

        start_idx = index
        end_idx = next_index(loop, index)
        p1 = loop%points(start_idx)
        p2 = loop%points(end_idx)
        start_point = [p1%x, p1%y, p1%z]
        end_point = [p2%x, p2%y, p2%z]
    end subroutine segment_endpoints

    subroutine extract_path_segment(diagnostic, index, start_point, end_point)
        type(segmented_rogowski_t), intent(in) :: diagnostic
        integer, intent(in) :: index
        real(dp), intent(out) :: start_point(3)
        real(dp), intent(out) :: end_point(3)

        start_point = [diagnostic%path(index)%x, diagnostic%path(index)%y, &
            diagnostic%path(index)%z]
        end_point = [diagnostic%path(index + 1)%x, &
            diagnostic%path(index + 1)%y, diagnostic%path(index + 1)%z]
    end subroutine extract_path_segment

    integer function next_index(loop, current)
        type(flux_loop_t), intent(in) :: loop
        integer, intent(in) :: current
        integer :: npts

        npts = size(loop%points)
        if (current == npts) then
            next_index = 1
        else
            next_index = current + 1
        end if
    end function next_index

    logical function points_match(a, b)
        type(loop_point_t), intent(in) :: a
        type(loop_point_t), intent(in) :: b
        real(dp) :: dx
        real(dp) :: dy
        real(dp) :: dz

        dx = a%x - b%x
        dy = a%y - b%y
        dz = a%z - b%z
        points_match = sqrt(dx * dx + dy * dy + dz * dz) < closure_tolerance
    end function points_match

    subroutine abort_with(message)
        character(len=*), intent(in) :: message
        error stop trim(message)
    end subroutine abort_with

    subroutine rotate_point(point, angle, rotated)
        real(dp), intent(in) :: point(3)
        real(dp), intent(in) :: angle
        real(dp), intent(out) :: rotated(3)

        real(dp) :: cang
        real(dp) :: sang

        cang = cos(angle)
        sang = sin(angle)
        rotated(1) = point(1) * cang - point(2) * sang
        rotated(2) = point(1) * sang + point(2) * cang
        rotated(3) = point(3)
    end subroutine rotate_point

    subroutine scale_coils_to_cgs(field)
        type(biotsavart_field_t), intent(inout) :: field

        if (.not. allocated(field%coils%x)) return
        field%coils%x = field%coils%x * meters_to_cm
        field%coils%y = field%coils%y * meters_to_cm
        field%coils%z = field%coils%z * meters_to_cm
        field%coils%current = field%coils%current * amps_to_statamp
    end subroutine scale_coils_to_cgs
end module tiago_vacuum_forward
