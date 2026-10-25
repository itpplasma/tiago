module tiago_sensor_signals
    !! Signals and exact Jacobians of a sensor set (tiago_sensor_set).
    !!
    !! Raw loop signals follow DIAGNO: the loop integral of A, plus the linked
    !! plasma toroidal flux for diamagnetic loops (idia = 1), times the turn
    !! scale. Saddle and compensation loops are plain flux loops. Segmented
    !! Rogowski coils and B-probes use the existing evaluators. Output signals
    !! are y = C raw, so every Jacobian is C times the raw Jacobian.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_sensor_set, only: sensor_set_t
    use tiago_diagnostic_types, only: flux_loop_t, segmented_rogowski_t, bprobe_t
    use tiago_vacuum_forward, only: vacuum_solver_t, quadrature_rule_t
    implicit none
    private

    public :: sensor_set_signals
    public :: sensor_set_coil_jacobian
    public :: sensor_set_plasma_jacobian
    public :: sensor_set_linked_flux_row

contains

    subroutine sensor_set_signals(solver, set, y, rule, linked_flux)
        !! linked_flux overrides the plasma toroidal flux added to diamagnetic
        !! loops (default: the solver's plasma flux, zero in vacuum).
        type(vacuum_solver_t), intent(in) :: solver
        type(sensor_set_t), intent(in) :: set
        real(dp), allocatable, intent(out) :: y(:)
        type(quadrature_rule_t), intent(in), optional :: rule
        real(dp), intent(in), optional :: linked_flux
        real(dp), allocatable :: raw(:), part(:)
        type(flux_loop_t), allocatable :: loops(:)
        type(segmented_rogowski_t), allocatable :: segs(:)
        real(dp) :: phi
        integer :: n, j

        allocate(raw(set%n_raw()))
        n = 0
        if (size(set%loops) > 0) then
            loops = set%loops
            call solver%flux_loops(loops, part, rule)
            if (present(linked_flux)) then
                phi = linked_flux
            else
                phi = solver%plasma_flux(required=.false.)
            end if
            do j = 1, size(loops)
                if (loops(j)%idia == 1) part(j) = part(j) + phi
                part(j) = part(j) * loops(j)%turn_scale
            end do
            raw(n + 1:n + size(part)) = part
            n = n + size(part)
        end if
        if (size(set%segs) > 0) then
            segs = set%segs
            call solver%segrog(segs, part, rule)
            raw(n + 1:n + size(part)) = part * segs%turn_scale
            n = n + size(part)
        end if
        if (size(set%probes) > 0) then
            call solver%bprobes(set%probes, part)
            raw(n + 1:n + size(part)) = part * set%probes%turn_scale
        end if
        allocate(y(set%n_signals()))
        call set%combine(raw, y)
    end subroutine sensor_set_signals

    subroutine sensor_set_coil_jacobian(solver, set, jac, rule)
        !! d y / d EXTCUR(g), one column per coil group (exact: y is linear).
        type(vacuum_solver_t), intent(in) :: solver
        type(sensor_set_t), intent(in) :: set
        real(dp), allocatable, intent(out) :: jac(:, :)
        type(quadrature_rule_t), intent(in), optional :: rule
        type(flux_loop_t), allocatable :: loops(:)
        type(segmented_rogowski_t), allocatable :: segs(:)
        type(bprobe_t), allocatable :: probes(:)
        real(dp), allocatable :: fr(:, :), sr(:, :), pr(:, :)

        call allocated_copies(set, loops, segs, probes)
        call solver%response(loops, segs, probes, fr, sr, pr, rule)
        call set%combine_matrix(stack(set, fr, sr, pr), jac)
    end subroutine sensor_set_coil_jacobian

    subroutine sensor_set_plasma_jacobian(solver, set, mode_jac, shape_jac, rule)
        !! d y / d (VMEC boundary field coefficients), and optionally d y / d
        !! (boundary geometry coefficients) at fixed field coefficients; columns as
        !! solver%plasma%mode_column_name / shape_column_name. The linked flux of
        !! diamagnetic loops (phiedge) is a separate parameter, see
        !! sensor_set_linked_flux_row.
        type(vacuum_solver_t), intent(in) :: solver
        type(sensor_set_t), intent(in) :: set
        real(dp), allocatable, intent(out) :: mode_jac(:, :)
        real(dp), allocatable, intent(out), optional :: shape_jac(:, :)
        type(quadrature_rule_t), intent(in), optional :: rule
        type(flux_loop_t), allocatable :: loops(:)
        type(segmented_rogowski_t), allocatable :: segs(:)
        type(bprobe_t), allocatable :: probes(:)
        real(dp), allocatable :: fr(:, :), sr(:, :), pr(:, :), fs(:, :), ss(:, :), ps(:, :)

        call allocated_copies(set, loops, segs, probes)
        if (present(shape_jac)) then
            call solver%plasma_response(loops, segs, probes, fr, sr, pr, rule, fs, ss, ps)
            call set%combine_matrix(stack(set, fs, ss, ps), shape_jac)
        else
            call solver%plasma_response(loops, segs, probes, fr, sr, pr, rule)
        end if
        call set%combine_matrix(stack(set, fr, sr, pr), mode_jac)
    end subroutine sensor_set_plasma_jacobian

    subroutine sensor_set_linked_flux_row(set, column)
        !! d y / d (linked plasma toroidal flux).
        type(sensor_set_t), intent(in) :: set
        real(dp), allocatable, intent(out) :: column(:)
        real(dp), allocatable :: raw(:)
        integer :: j

        allocate(raw(set%n_raw()), column(set%n_signals()))
        raw = 0.0_dp
        do j = 1, size(set%loops)
            if (set%loops(j)%idia == 1) raw(j) = set%loops(j)%turn_scale
        end do
        call set%combine(raw, column)
    end subroutine sensor_set_linked_flux_row

    subroutine allocated_copies(set, loops, segs, probes)
        !! The solver treats unallocated arrays as absent families.
        type(sensor_set_t), intent(in) :: set
        type(flux_loop_t), allocatable, intent(out) :: loops(:)
        type(segmented_rogowski_t), allocatable, intent(out) :: segs(:)
        type(bprobe_t), allocatable, intent(out) :: probes(:)

        if (size(set%loops) > 0) loops = set%loops
        if (size(set%segs) > 0) segs = set%segs
        if (size(set%probes) > 0) probes = set%probes
    end subroutine allocated_copies

    function stack(set, fr, sr, pr) result(raw)
        !! Raw Jacobian rows in sensor-set order, with the turn scales applied.
        type(sensor_set_t), intent(in) :: set
        real(dp), intent(in) :: fr(:, :), sr(:, :), pr(:, :)
        real(dp), allocatable :: raw(:, :)
        integer :: nl, ns, j, nc

        nl = size(set%loops)
        ns = size(set%segs)
        nc = max(size(fr, 2), size(sr, 2), size(pr, 2))
        allocate(raw(set%n_raw(), nc))
        raw = 0.0_dp
        do j = 1, nl
            raw(j, :) = fr(j, :) * set%loops(j)%turn_scale
        end do
        do j = 1, ns
            raw(nl + j, :) = sr(j, :) * set%segs(j)%turn_scale
        end do
        do j = 1, size(set%probes)
            raw(nl + ns + j, :) = pr(j, :) * set%probes(j)%turn_scale
        end do
    end function stack
end module tiago_sensor_signals
