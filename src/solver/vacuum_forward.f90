module tiago_vacuum_forward
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    use neo_biotsavart_field, only: biotsavart_field_t
    use tiago_coil_loader, only: load_coils_into_field
    use neo_biotsavart, only: coils_init
    use tiago_diagnostic_types, only: flux_loop_t, segmented_rogowski_t, &
        loop_point_t
    use tiago_plasma_support, only: plasma_support_t
    implicit none
    private

    real(dp), parameter :: min_segment_length = 1.0e-10_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: two_pi = 2.0_dp * pi
    real(dp), parameter :: meters_to_cm = 100.0_dp
    real(dp), parameter :: gauss_to_tesla = 1.0e-4_dp
    real(dp), parameter :: amps_to_statamp = 2.9979245368431e9_dp
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
        type(plasma_support_t) :: plasma
    contains
        procedure :: init => vacuum_solver_init
        procedure :: finalize => vacuum_solver_finalize
        procedure :: flux_loops => vacuum_solver_flux_loops
        procedure :: segrog => vacuum_solver_segrog
        procedure :: set_nfp => vacuum_solver_set_nfp
        procedure :: flux_and_segrog => vacuum_solver_flux_and_segrog
        procedure :: enable_plasma_from_vmec => vacuum_solver_enable_plasma_from_vmec
        procedure :: disable_plasma => vacuum_solver_disable_plasma
        procedure :: has_plasma => vacuum_solver_has_plasma
    end type vacuum_solver_t

    public :: vacuum_solver_t
    public :: quadrature_rule_t
    public :: quadrature_override_t

contains
    subroutine vacuum_solver_init(self, coil_file, coil_extcur)
        class(vacuum_solver_t), intent(inout) :: self
        character(len=*), intent(in) :: coil_file
        character(len=*), intent(in), optional :: coil_extcur
        real(dp) :: no_points(0)

        if (len_trim(coil_file) == 0) then
            ! No coils: plasma-only evaluation.
            call coils_init(no_points, no_points, no_points, no_points, self%field%coils)
            self%is_ready = .true.
            return
        end if

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
        call self%plasma%finalize()
    end subroutine vacuum_solver_finalize

    subroutine vacuum_solver_enable_plasma_from_vmec(self, wout_file, nphi, ntheta)
        class(vacuum_solver_t), intent(inout) :: self
        character(len=*), intent(in) :: wout_file
        integer(i32), intent(in) :: nphi
        integer(i32), intent(in) :: ntheta

        call self%plasma%init_from_vmec(trim(wout_file), nphi, ntheta)
    end subroutine vacuum_solver_enable_plasma_from_vmec

    subroutine vacuum_solver_disable_plasma(self)
        class(vacuum_solver_t), intent(inout) :: self
        call self%plasma%finalize()
    end subroutine vacuum_solver_disable_plasma

    logical function vacuum_solver_has_plasma(self)
        class(vacuum_solver_t), intent(in) :: self
        vacuum_solver_has_plasma = self%plasma%has_data()
    end function vacuum_solver_has_plasma

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
        if (self%plasma%has_data()) then
            call add_plasma_flux(self, loops, fluxes, default_rule, overrides)
        end if
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
            voltages(i) = evaluate_segrog_signal(self%field, diagnostics(i), rule)
        end do
!$omp end parallel do
        if (self%plasma%has_data()) then
            call add_plasma_segrog(self, diagnostics, voltages, default_rule, overrides)
        end if
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

        ! Single flattened loop: compute all flux loops, then all segrog diagnostics
        ! This distributes work evenly across threads without intermediate barriers
!$omp parallel do default(shared) private(i, rule)
        do i = 1, size(loops) + size(diagnostics)
            if (i <= size(loops)) then
                rule = select_rule(loops(i)%label, default_rule, overrides)
                fluxes(i) = evaluate_loop_flux(self%field, loops(i), rule, &
                    self%nfp)
            else
                rule = select_rule(diagnostics(i - size(loops))%label, default_rule, &
                    overrides)
                voltages(i - size(loops)) = evaluate_segrog_signal(self%field, &
                    diagnostics(i - size(loops)), rule)
            end if
        end do
!$omp end parallel do
        if (self%plasma%has_data()) then
            call add_plasma_flux(self, loops, fluxes, default_rule, overrides)
            call add_plasma_segrog(self, diagnostics, voltages, default_rule, overrides)
        end if
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
        real(dp) :: dl(3)
        real(dp) :: start_point(3)
        real(dp) :: end_point(3)
        real(dp) :: weight

        flux = 0.0_dp
        samples = max(1_i32, rule%samples_per_segment)
        weight = 1.0_dp / real(samples, dp)

        ! DIAGNO semantics: a flux loop is always closed. For iflflg=1 the
        ! closing point is the first point rotated by one field period and the
        ! result is multiplied by nfp.
        do seg = 1, size(loop%points)
            call loop_segment(loop, seg, nfp, start_point, end_point)
            dl = end_point - start_point
            flux = flux + integrate_segment(field, start_point, dl, samples, weight)
        end do

        flux = flux * maxwell_to_weber
        if (loop%one_period) flux = flux * real(nfp, dp)
    end function evaluate_loop_flux

    real(dp) function integrate_segment(field, start_point, dl, samples, weight)
        type(biotsavart_field_t), intent(in) :: field
        real(dp), intent(in) :: start_point(3)
        real(dp), intent(in) :: dl(3)
        integer(i32), intent(in) :: samples
        real(dp), intent(in) :: weight

        integer :: s, coil_i
        real(dp) :: a_field(3)
        real(dp) :: step
        real(dp) :: sample_point_cm(3)
        real(dp) :: dl_cm(3)
        real(dp) :: start_point_cm(3)
        real(dp) :: dx_i(3), dx_f(3), dl_coil(3)
        real(dp) :: R_i, R_f, L, eps, log_term
        real(dp) :: clight_param
        integer :: n_coil_segments

        integrate_segment = 0.0_dp
        start_point_cm = start_point * meters_to_cm
        dl_cm = dl * meters_to_cm

        ! Cache coil parameters
        n_coil_segments = size(field%coils%x) - 1
        clight_param = 2.99792458d10  ! Speed of light in CGS

        ! Coil loop OUTERMOST - keep coil data in L1D cache
        ! This achieves ~95% L1D hit rate vs ~10% with sample loop outermost
        do coil_i = 1, n_coil_segments
            ! Get coil segment vector (stays in L1D for all samples)
            dl_coil(1) = field%coils%x(coil_i + 1) - field%coils%x(coil_i)
            dl_coil(2) = field%coils%y(coil_i + 1) - field%coils%y(coil_i)
            dl_coil(3) = field%coils%z(coil_i + 1) - field%coils%z(coil_i)
            L = sqrt(dl_coil(1)**2 + dl_coil(2)**2 + dl_coil(3)**2)
            ! Duplicated points and inter-coil jumpers carry no field (and L=0 gives 0/0).
            if (L == 0.0_dp .or. field%coils%current(coil_i) == 0.0_dp) cycle

            ! Sample loop inner - compute contribution from coil_i to all samples
            do s = 1, samples
                step = (real(s, dp) - 0.5_dp) / real(samples, dp)
                sample_point_cm = start_point_cm + step * dl_cm

                ! Hanson-Hirshman formula for vector potential contribution
                ! from coil segment coil_i to sample_point_cm
                dx_i(1) = sample_point_cm(1) - field%coils%x(coil_i)
                dx_i(2) = sample_point_cm(2) - field%coils%y(coil_i)
                dx_i(3) = sample_point_cm(3) - field%coils%z(coil_i)
                R_i = sqrt(dx_i(1)**2 + dx_i(2)**2 + dx_i(3)**2)

                dx_f(1) = sample_point_cm(1) - field%coils%x(coil_i + 1)
                dx_f(2) = sample_point_cm(2) - field%coils%y(coil_i + 1)
                dx_f(3) = sample_point_cm(3) - field%coils%z(coil_i + 1)
                R_f = sqrt(dx_f(1)**2 + dx_f(2)**2 + dx_f(3)**2)

                eps = L / (R_i + R_f)
                log_term = log((1.0_dp + eps) / (1.0_dp - eps))

                ! Accumulate vector potential from this coil segment
                a_field(1) = (field%coils%current(coil_i) / clight_param) * &
                    (dl_coil(1) / L) * log_term
                a_field(2) = (field%coils%current(coil_i) / clight_param) * &
                    (dl_coil(2) / L) * log_term
                a_field(3) = (field%coils%current(coil_i) / clight_param) * &
                    (dl_coil(3) / L) * log_term

                integrate_segment = integrate_segment + &
                    weight * (a_field(1) * dl_cm(1) + &
                              a_field(2) * dl_cm(2) + &
                              a_field(3) * dl_cm(3))
            end do
        end do
    end function integrate_segment

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

        samples = max(1_i32, rule%samples_per_segment)
        weight = 1.0_dp / real(samples, dp)
        voltage = 0.0_dp

        do seg = 1, size(diagnostic%path) - 1
            call extract_path_segment(diagnostic, seg, start_point, end_point)
            dl = end_point - start_point
            norm_dl = max(min_segment_length, sqrt(sum(dl**2)))
            tangent = dl / norm_dl
            voltage = voltage + diagnostic%segment_area(seg) * &
                integrate_segrog_segment(field, start_point, dl, norm_dl, tangent, &
                samples, weight)
        end do
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

        integer :: s, coil_i
        real(dp) :: step
        real(dp) :: sample_point(3)
        real(dp) :: b_field(3)
        real(dp) :: sample_point_cm(3)
        real(dp) :: dx_i(3), dx_f(3), dl_coil(3), dl_coil_hat(3), dx_i_hat(3)
        real(dp) :: R_i, R_f, L, eps, cross_prod(3)
        real(dp) :: clight_param
        integer :: n_coil_segments

        integrate_segrog_segment = 0.0_dp

        ! Cache coil parameters
        n_coil_segments = size(field%coils%x) - 1
        clight_param = 2.99792458d10  ! Speed of light in CGS

        ! Coil loop OUTERMOST - keep coil data in L1D cache
        ! This achieves ~95% L1D hit rate vs ~10% with sample loop outermost
        do coil_i = 1, n_coil_segments
            ! Get coil segment vector (stays in L1D for all samples)
            dl_coil(1) = field%coils%x(coil_i + 1) - field%coils%x(coil_i)
            dl_coil(2) = field%coils%y(coil_i + 1) - field%coils%y(coil_i)
            dl_coil(3) = field%coils%z(coil_i + 1) - field%coils%z(coil_i)
            L = sqrt(dl_coil(1)**2 + dl_coil(2)**2 + dl_coil(3)**2)
            ! Duplicated points and inter-coil jumpers carry no field (and L=0 gives 0/0).
            if (L == 0.0_dp .or. field%coils%current(coil_i) == 0.0_dp) cycle
            dl_coil_hat(1) = dl_coil(1) / L
            dl_coil_hat(2) = dl_coil(2) / L
            dl_coil_hat(3) = dl_coil(3) / L

            ! Sample loop inner - compute contribution from coil_i to all samples
            do s = 1, samples
                step = (real(s, dp) - 0.5_dp) / real(samples, dp)
                sample_point = start_point + step * dl
                sample_point_cm = sample_point * meters_to_cm

                ! Hanson-Hirshman formula for magnetic field contribution
                ! from coil segment coil_i to sample_point_cm
                dx_i(1) = sample_point_cm(1) - field%coils%x(coil_i)
                dx_i(2) = sample_point_cm(2) - field%coils%y(coil_i)
                dx_i(3) = sample_point_cm(3) - field%coils%z(coil_i)
                R_i = sqrt(dx_i(1)**2 + dx_i(2)**2 + dx_i(3)**2)

                dx_f(1) = sample_point_cm(1) - field%coils%x(coil_i + 1)
                dx_f(2) = sample_point_cm(2) - field%coils%y(coil_i + 1)
                dx_f(3) = sample_point_cm(3) - field%coils%z(coil_i + 1)
                R_f = sqrt(dx_f(1)**2 + dx_f(2)**2 + dx_f(3)**2)

                ! Normalized vector from segment start to sample point
                dx_i_hat(1) = dx_i(1) / R_i
                dx_i_hat(2) = dx_i(2) / R_i
                dx_i_hat(3) = dx_i(3) / R_i

                ! Cross product: dl_hat × dx_i_hat
                cross_prod(1) = dl_coil_hat(2) * dx_i_hat(3) - &
                                dl_coil_hat(3) * dx_i_hat(2)
                cross_prod(2) = dl_coil_hat(3) * dx_i_hat(1) - &
                                dl_coil_hat(1) * dx_i_hat(3)
                cross_prod(3) = dl_coil_hat(1) * dx_i_hat(2) - &
                                dl_coil_hat(2) * dx_i_hat(1)

                eps = L / (R_i + R_f)

                ! Accumulate magnetic field from this coil segment
                b_field(1) = (field%coils%current(coil_i) / clight_param) * &
                    cross_prod(1) * (1.0_dp / R_f) * &
                    (2.0_dp * eps / (1.0_dp - eps**2))
                b_field(2) = (field%coils%current(coil_i) / clight_param) * &
                    cross_prod(2) * (1.0_dp / R_f) * &
                    (2.0_dp * eps / (1.0_dp - eps**2))
                b_field(3) = (field%coils%current(coil_i) / clight_param) * &
                    cross_prod(3) * (1.0_dp / R_f) * &
                    (2.0_dp * eps / (1.0_dp - eps**2))

                ! Convert from Gauss to Tesla and accumulate
                b_field = b_field * gauss_to_tesla

                integrate_segrog_segment = integrate_segrog_segment + weight * &
                    (b_field(1) * tangent(1) + &
                     b_field(2) * tangent(2) + &
                     b_field(3) * tangent(3)) * norm_dl
            end do
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

    subroutine loop_segment(loop, seg, nfp, start_point, end_point)
        !! Segment seg of a closed flux loop; for iflflg=1 loops the closing
        !! point is the first point rotated by one field period.
        type(flux_loop_t), intent(in) :: loop
        integer, intent(in) :: seg
        integer(i32), intent(in) :: nfp
        real(dp), intent(out) :: start_point(3), end_point(3)

        call segment_endpoints(loop, seg, start_point, end_point)
        if (seg == size(loop%points) .and. loop%one_period) then
            call rotate_point(start_point_of(loop), two_pi / real(nfp, dp), end_point)
        end if
    end subroutine loop_segment

    subroutine add_plasma_flux(self, loops, fluxes, default_rule, overrides)
        !! Plasma part of every flux loop, sum over samples of A_plasma . dl,
        !! evaluated in one batch (the plasma kernel is itself parallel).
        class(vacuum_solver_t), intent(in) :: self
        type(flux_loop_t), intent(in) :: loops(:)
        real(dp), intent(inout) :: fluxes(:)
        type(quadrature_rule_t), intent(in), optional :: default_rule
        type(quadrature_override_t), intent(in), optional :: overrides(:)

        real(dp), allocatable :: points(:, :), dls(:, :), avec(:, :)
        integer, allocatable :: owner(:)
        real(dp) :: a(3), b(3), plasma_flux(size(loops))
        integer :: i, seg, s, n, samples

        n = 0
        do i = 1, size(loops)
            n = n + size(loops(i)%points) * &
                sample_count(select_rule(loops(i)%label, default_rule, overrides))
        end do
        allocate(points(n, 3), dls(n, 3), owner(n), avec(n, 3))
        n = 0
        do i = 1, size(loops)
            samples = sample_count(select_rule(loops(i)%label, default_rule, overrides))
            do seg = 1, size(loops(i)%points)
                call loop_segment(loops(i), seg, self%nfp, a, b)
                do s = 1, samples
                    n = n + 1
                    points(n, :) = a + (real(s, dp) - 0.5_dp) / real(samples, dp) * (b - a)
                    dls(n, :) = (b - a) / real(samples, dp)
                    owner(n) = i
                end do
            end do
        end do

        call self%plasma%sample_vector_potential(points, avec)
        plasma_flux = 0.0_dp
        do s = 1, n
            plasma_flux(owner(s)) = plasma_flux(owner(s)) + dot_product(avec(s, :), dls(s, :))
        end do
        do i = 1, size(loops)
            if (loops(i)%one_period) plasma_flux(i) = plasma_flux(i) * real(self%nfp, dp)
        end do
        fluxes = fluxes + plasma_flux
    end subroutine add_plasma_flux

    subroutine add_plasma_segrog(self, diagnostics, voltages, default_rule, overrides)
        !! Plasma part of every segmented Rogowski, sum of eff_area * B_plasma . dl.
        class(vacuum_solver_t), intent(in) :: self
        type(segmented_rogowski_t), intent(in) :: diagnostics(:)
        real(dp), intent(inout) :: voltages(:)
        type(quadrature_rule_t), intent(in), optional :: default_rule
        type(quadrature_override_t), intent(in), optional :: overrides(:)

        real(dp), allocatable :: points(:, :), dls(:, :), bvec(:, :)
        integer, allocatable :: owner(:)
        real(dp) :: a(3), b(3)
        integer :: i, seg, s, n, samples

        n = 0
        do i = 1, size(diagnostics)
            n = n + (size(diagnostics(i)%path) - 1) * &
                sample_count(select_rule(diagnostics(i)%label, default_rule, overrides))
        end do
        allocate(points(n, 3), dls(n, 3), owner(n), bvec(n, 3))
        n = 0
        do i = 1, size(diagnostics)
            samples = sample_count(select_rule(diagnostics(i)%label, default_rule, overrides))
            do seg = 1, size(diagnostics(i)%path) - 1
                call extract_path_segment(diagnostics(i), seg, a, b)
                do s = 1, samples
                    n = n + 1
                    points(n, :) = a + (real(s, dp) - 0.5_dp) / real(samples, dp) * (b - a)
                    dls(n, :) = (b - a) / real(samples, dp) * diagnostics(i)%segment_area(seg)
                    owner(n) = i
                end do
            end do
        end do

        call self%plasma%sample_bfield(points, bvec)
        do s = 1, n
            voltages(owner(s)) = voltages(owner(s)) + dot_product(bvec(s, :), dls(s, :))
        end do
    end subroutine add_plasma_segrog

    integer function sample_count(rule)
        type(quadrature_rule_t), intent(in) :: rule
        sample_count = max(1_i32, rule%samples_per_segment)
    end function sample_count

    function start_point_of(loop) result(point)
        type(flux_loop_t), intent(in) :: loop
        real(dp) :: point(3)

        point = [loop%points(1)%x, loop%points(1)%y, loop%points(1)%z]
    end function start_point_of

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
