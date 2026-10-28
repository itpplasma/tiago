module tiago_wall_synthetic
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_signal_source, only: signal_window_t
    use tiago_wall_energy, only: wall_energy
    use tiago_wall_flux, only: wall_flux
    use tiago_wall_thermal, only: wall_thermal
    use tiago_wall_radiance, only: wall_radiance
    use tiago_recorded_affine, only: recorded_affine
    implicit none
    private
    integer, parameter, public :: WALL_OK = 0, WALL_ERROR = 1
    integer, parameter, public :: PARTICLE_COUNT = 1, PARTICLE_RATE = 2
    integer, parameter, public :: IR_VALID = 0, IR_UNMAPPED = 1, IR_SATURATED = 2
    character(*), parameter, public :: SEMI_INFINITE_REGIME = &
        'semi-infinite-constant-properties-no-cooling-no-lateral'

    type, public :: loss_event_block_t
        real(dp), allocatable :: time_s(:), kinetic_energy_eV(:), weight(:)
        real(dp), allocatable :: deposited_fraction(:)
        integer, allocatable :: facet(:), wall_hit(:), exit_code(:)
        integer :: normalization = 0
        character(:), allocatable :: source_id, version, frame, position_kind
    end type
    type, public :: injection_history_t
        real(dp), allocatable :: edge_s(:), amplitude(:)
    end type
    type, public :: surface_deposition_t
        real(dp), allocatable :: edge_s(:), area_m2(:), energy_J(:, :)
        character(:), allocatable :: mesh_id, version, source_id, frame, position_kind
        integer :: normalization = 0
    end type
    type, public :: semi_infinite_material_t
        real(dp), allocatable :: conductivity_W_mK(:), capacity_J_m3K(:), initial_K(:)
        character(:), allocatable :: regime, version
    end type
    type, public :: ir_camera_t
        real(dp), allocatable :: exposure_edge_s(:, :)
        real(dp), allocatable :: quadrature_time_s(:), quadrature_weight_s(:)
        integer, allocatable :: quadrature_exposure(:)
        real(dp), allocatable :: wavelength_m(:), wavelength_weight_m(:), transmission(:)
        real(dp), allocatable :: emissivity(:, :), background_W_m3sr(:, :)
        integer, allocatable :: link_pixel(:), link_cell(:)
        real(dp), allocatable :: etendue_m2sr(:)
        logical, allocatable :: link_visible(:)
        real(dp), allocatable :: gain(:), offset(:), saturation(:), noise_variance(:, :)
        integer :: npixel = 0
        logical :: integrated = .true.
        character(:), allocatable :: units, geometry_version, calibration_version
        character(:), allocatable :: source_version, thermal_version
    end type
    type, public :: wall_ir_result_t
        type(signal_window_t), allocatable :: channel(:)
        logical, allocatable :: valid(:, :)
        integer, allocatable :: flags(:, :)
        real(dp), allocatable :: exposure_edge_s(:, :)
    end type
    public :: init_deposition, accumulate_losses, deposition_flux
    public :: surface_temperature, temperature_energy_jvp, temperature_energy_vjp
    public :: synthesize_ir, ir_jvp, ir_vjp, ir_valid_window
contains
    logical function present_text(s)
        character(:), allocatable, intent(in) :: s
        present_text = .false.
        if (.not. allocated(s)) return
        present_text = len_trim(s) > 0
    end function

    logical function increasing(x)
        real(dp), intent(in) :: x(:)
        increasing = .false.
        if (size(x) < 2) return
        if (.not. all(ieee_is_finite(x))) return
        increasing = all(x(2:) > x(:size(x) - 1))
    end function

    subroutine init_deposition(surface, edge_s, area_m2, mesh_id, version, status)
        type(surface_deposition_t), intent(out) :: surface
        real(dp), intent(in) :: edge_s(:), area_m2(:)
        character(*), intent(in) :: mesh_id, version
        integer, intent(out) :: status
        status = WALL_ERROR
        if (.not. increasing(edge_s)) return
        if (size(area_m2) < 1) return
        if (.not. all(ieee_is_finite(area_m2))) return
        if (any(area_m2 <= 0)) return
        if (len_trim(mesh_id) == 0 .or. len_trim(version) == 0) return
        surface%edge_s = edge_s
        surface%area_m2 = area_m2
        allocate (surface%energy_J(size(area_m2), size(edge_s) - 1), source=0.0_dp)
        surface%mesh_id = mesh_id
        surface%version = version
        status = WALL_OK
    end subroutine

    logical function valid_surface(s)
        type(surface_deposition_t), intent(in) :: s
        valid_surface = .false.
        if (.not. allocated(s%edge_s)) return
        if (.not. allocated(s%area_m2)) return
        if (.not. allocated(s%energy_J)) return
        if (.not. present_text(s%mesh_id)) return
        if (.not. present_text(s%version)) return
        if (.not. increasing(s%edge_s)) return
        if (size(s%area_m2) < 1) return
        if (any(shape(s%energy_J) /= [size(s%area_m2), size(s%edge_s) - 1])) return
        if (.not. all(ieee_is_finite(s%area_m2))) return
        if (.not. all(ieee_is_finite(s%energy_J))) return
        if (any(s%area_m2 <= 0)) return
        if (any(s%energy_J < 0)) return
        valid_surface = .true.
    end function

    subroutine accumulate_losses(surface, events, status, history)
        type(surface_deposition_t), intent(inout) :: surface
        type(loss_event_block_t), intent(in) :: events
        integer, intent(out) :: status
        type(injection_history_t), optional, intent(in) :: history
        real(dp) :: energy, deriv, overlap, count, lo, hi
        real(dp), allocatable :: delta(:, :)
        integer :: i, j, k, n, cell, nb
        status = WALL_ERROR
        if (.not. valid_surface(surface)) return
        if (.not. allocated(events%time_s)) return
        if (.not. allocated(events%kinetic_energy_eV)) return
        if (.not. allocated(events%weight)) return
        if (.not. allocated(events%deposited_fraction)) return
        if (.not. allocated(events%facet)) return
        if (.not. allocated(events%wall_hit)) return
        if (.not. allocated(events%exit_code)) return
        if (.not. present_text(events%source_id)) return
        if (.not. present_text(events%version)) return
        if (.not. present_text(events%frame)) return
        if (.not. present_text(events%position_kind)) return
        if (events%position_kind /= 'guiding-centre' .and. &
            events%position_kind /= 'physical-impact') return
        if (allocated(surface%source_id)) then
            if (.not. present_text(surface%frame)) return
            if (.not. present_text(surface%position_kind)) return
            if (surface%normalization /= events%normalization) return
            if (surface%frame /= events%frame) return
            if (surface%position_kind /= events%position_kind) return
            if (surface%source_id /= events%source_id//'@'//events%version) return
        end if
        n = size(events%time_s)
        if (size(events%kinetic_energy_eV) /= n) return
        if (size(events%weight) /= n) return
        if (size(events%deposited_fraction) /= n) return
        if (size(events%facet) /= n .or. size(events%wall_hit) /= n) return
        if (size(events%exit_code) /= n) return
        if (events%normalization /= PARTICLE_COUNT .and. &
            events%normalization /= PARTICLE_RATE) return
        if (events%normalization == PARTICLE_COUNT) then
            if (present(history)) return
        else
            if (.not. present(history)) return
            if (.not. allocated(history%edge_s)) return
            if (.not. allocated(history%amplitude)) return
            if (.not. increasing(history%edge_s)) return
            if (size(history%amplitude) /= size(history%edge_s) - 1) return
            if (.not. all(ieee_is_finite(history%amplitude))) return
            if (any(history%amplitude < 0)) return
        end if
        ! Validate the entire block before mutating the reduction.
        do i = 1, n
            if (events%wall_hit(i) /= 1 .or. events%exit_code(i) /= 2) cycle
            if (events%facet(i) < 1 .or. events%facet(i) > size(surface%area_m2)) return
            if (.not. ieee_is_finite(events%time_s(i))) return
            if (.not. ieee_is_finite(events%kinetic_energy_eV(i))) return
            if (.not. ieee_is_finite(events%weight(i))) return
            if (.not. ieee_is_finite(events%deposited_fraction(i))) return
            if (events%kinetic_energy_eV(i) < 0 .or. events%weight(i) < 0) return
            if (events%deposited_fraction(i) < 0 .or. &
                events%deposited_fraction(i) > 1) return
            if (events%normalization == PARTICLE_RATE) then
                if (events%time_s(i) < 0) return
                do j = 1, size(history%amplitude)
                    if (history%amplitude(j) == 0) cycle
                    if (history%edge_s(j) + events%time_s(i) < surface%edge_s(1)) return
                    if (history%edge_s(j + 1) + events%time_s(i) &
                        > surface%edge_s(size(surface%edge_s))) return
                end do
            else
                if (locate_bin(surface%edge_s, events%time_s(i)) == 0) return
            end if
        end do
        nb = size(surface%edge_s) - 1
        allocate (delta(size(surface%area_m2), nb), source=0.0_dp)
        do i = 1, n
            if (events%wall_hit(i) /= 1 .or. events%exit_code(i) /= 2) cycle
            cell = events%facet(i)
            if (events%normalization == PARTICLE_COUNT) then
                ! Bins are left-closed/right-open, including the final edge in the last.
                k = locate_bin(surface%edge_s, events%time_s(i))
                if (k == 0) cycle
                call wall_energy(events%kinetic_energy_eV(i), events%weight(i), &
                    events%deposited_fraction(i), energy, deriv)
                delta(cell, k) = delta(cell, k) + energy
            else
                do k = 1, nb
                    count = 0
                    do j = 1, size(history%amplitude)
                        lo = max(surface%edge_s(k), history%edge_s(j) + events%time_s(i))
                        hi = min(surface%edge_s(k + 1), &
                            history%edge_s(j + 1) + events%time_s(i))
                        overlap = max(0.0_dp, hi - lo)
                        count = count + events%weight(i)*history%amplitude(j)*overlap
                    end do
                    call wall_energy(events%kinetic_energy_eV(i), count, &
                        events%deposited_fraction(i), energy, deriv)
                    delta(cell, k) = delta(cell, k) + energy
                end do
            end if
        end do
        if (.not. all(ieee_is_finite(delta))) return
        if (.not. all(ieee_is_finite(surface%energy_J + delta))) return
        surface%energy_J = surface%energy_J + delta
        surface%source_id = events%source_id//'@'//events%version
        surface%frame = events%frame
        surface%position_kind = events%position_kind
        surface%normalization = events%normalization
        status = WALL_OK
    end subroutine

    integer function locate_bin(edge, t) result(k)
        real(dp), intent(in) :: edge(:), t
        integer :: lo, hi, mid
        k = 0
        if (t < edge(1) .or. t > edge(size(edge))) return
        lo = 1
        hi = size(edge) - 1
        do while (lo < hi)
            mid = (lo + hi)/2
            if (t < edge(mid + 1)) then
                hi = mid
            else
                lo = mid + 1
            end if
        end do
        k = lo
    end function

    subroutine deposition_flux(surface, flux, status)
        type(surface_deposition_t), intent(in) :: surface
        real(dp), allocatable, intent(out) :: flux(:, :)
        integer, intent(out) :: status
        integer :: j, i
        real(dp) :: deriv
        status = WALL_ERROR
        if (.not. valid_surface(surface)) return
        allocate (flux(size(surface%area_m2), size(surface%edge_s) - 1))
        !$omp parallel do collapse(2) private(deriv)
        do j = 1, size(flux, 2)
            do i = 1, size(flux, 1)
                call wall_flux(surface%energy_J(i, j), surface%area_m2(i), &
                    surface%edge_s(j + 1) - surface%edge_s(j), flux(i, j), deriv)
            end do
        end do
        !$omp end parallel do
        if (.not. all(ieee_is_finite(flux))) return
        status = WALL_OK
    end subroutine

    logical function valid_material(material, n)
        type(semi_infinite_material_t), intent(in) :: material
        integer, intent(in) :: n
        valid_material = .false.
        if (.not. present_text(material%regime)) return
        if (.not. present_text(material%version)) return
        if (material%regime /= SEMI_INFINITE_REGIME) return
        if (.not. allocated(material%conductivity_W_mK)) return
        if (.not. allocated(material%capacity_J_m3K)) return
        if (.not. allocated(material%initial_K)) return
        if (size(material%conductivity_W_mK) /= n) return
        if (size(material%capacity_J_m3K) /= n .or. size(material%initial_K) /= n) return
        if (.not. all(ieee_is_finite(material%conductivity_W_mK))) return
        if (.not. all(ieee_is_finite(material%capacity_J_m3K))) return
        if (.not. all(ieee_is_finite(material%initial_K))) return
        if (any(material%conductivity_W_mK <= 0)) return
        if (any(material%capacity_J_m3K <= 0) .or. any(material%initial_K <= 0)) return
        valid_material = .true.
    end function

    subroutine surface_temperature(surface, material, time_s, temperature, status)
        type(surface_deposition_t), intent(in) :: surface
        type(semi_infinite_material_t), intent(in) :: material
        real(dp), intent(in) :: time_s(:)
        real(dp), allocatable, intent(out) :: temperature(:, :)
        integer, intent(out) :: status
        integer :: i, j
        status = WALL_ERROR
        if (.not. valid_surface(surface)) return
        if (.not. valid_material(material, size(surface%area_m2))) return
        call temperature_energy_jvp(surface, material, time_s, surface%energy_J, &
            temperature, status)
        if (status /= WALL_OK) return
        do j = 1, size(time_s)
            do i = 1, size(surface%area_m2)
                temperature(i, j) = temperature(i, j) + material%initial_K(i)
            end do
        end do
        if (.not. all(ieee_is_finite(temperature))) status = WALL_ERROR
    end subroutine

    subroutine temperature_energy_jvp(surface, material, time_s, denergy, dtemp, status)
        type(surface_deposition_t), intent(in) :: surface
        type(semi_infinite_material_t), intent(in) :: material
        real(dp), intent(in) :: time_s(:), denergy(:, :)
        real(dp), allocatable, intent(out) :: dtemp(:, :)
        integer, intent(out) :: status
        integer :: i, j, k, n
        real(dp) :: q, dq, v, dv, lag0, lag1
        status = WALL_ERROR
        if (.not. valid_surface(surface)) return
        n = size(surface%area_m2)
        if (.not. valid_material(material, n)) return
        if (any(shape(denergy) /= shape(surface%energy_J))) return
        if (.not. all(ieee_is_finite(denergy))) return
        if (.not. all(ieee_is_finite(time_s))) return
        if (any(time_s < surface%edge_s(1))) return
        allocate (dtemp(n, size(time_s)), source=0.0_dp)
        !$omp parallel do collapse(2) private(k,q,dq,v,dv,lag0,lag1)
        do j = 1, size(time_s)
            do i = 1, n
                do k = 1, size(surface%energy_J, 2)
                    lag0 = max(0.0_dp, time_s(j) - surface%edge_s(k))
                    lag1 = max(0.0_dp, time_s(j) - surface%edge_s(k + 1))
                    call wall_flux(denergy(i, k), surface%area_m2(i), &
                        surface%edge_s(k + 1) - surface%edge_s(k), q, dq)
                    call wall_thermal(q, material%conductivity_W_mK(i), &
                        material%capacity_J_m3K(i), lag0, lag1, v, dv)
                    dtemp(i, j) = dtemp(i, j) + v
                end do
            end do
        end do
        !$omp end parallel do
        if (.not. all(ieee_is_finite(dtemp))) return
        status = WALL_OK
    end subroutine

    subroutine temperature_energy_vjp(surface, material, time_s, bar_temp, &
            bar_energy, status)
        type(surface_deposition_t), intent(in) :: surface
        type(semi_infinite_material_t), intent(in) :: material
        real(dp), intent(in) :: time_s(:), bar_temp(:, :)
        real(dp), allocatable, intent(out) :: bar_energy(:, :)
        integer, intent(out) :: status
        integer :: i, j, k, n
        real(dp) :: q, dq, v, dv, lag0, lag1
        status = WALL_ERROR
        if (.not. valid_surface(surface)) return
        n = size(surface%area_m2)
        if (.not. valid_material(material, n)) return
        if (any(shape(bar_temp) /= [n, size(time_s)])) return
        if (.not. all(ieee_is_finite(bar_temp))) return
        if (.not. all(ieee_is_finite(time_s))) return
        if (any(time_s < surface%edge_s(1))) return
        allocate (bar_energy(n, size(surface%energy_J, 2)), source=0.0_dp)
        !$omp parallel do collapse(2) private(j,q,dq,v,dv,lag0,lag1)
        do k = 1, size(surface%energy_J, 2)
            do i = 1, n
                call wall_flux(1.0_dp, surface%area_m2(i), &
                    surface%edge_s(k + 1) - surface%edge_s(k), q, dq)
                do j = 1, size(time_s)
                    lag0 = max(0.0_dp, time_s(j) - surface%edge_s(k))
                    lag1 = max(0.0_dp, time_s(j) - surface%edge_s(k + 1))
                    call wall_thermal(q, material%conductivity_W_mK(i), &
                        material%capacity_J_m3K(i), lag0, lag1, v, dv)
                    bar_energy(i, k) = bar_energy(i, k) + v*bar_temp(i, j)
                end do
            end do
        end do
        !$omp end parallel do
        if (.not. all(ieee_is_finite(bar_energy))) return
        status = WALL_OK
    end subroutine

    logical function valid_camera(cam, temp)
        type(ir_camera_t), intent(in) :: cam
        real(dp), intent(in) :: temp(:, :)
        integer :: n, e, l, j, ncell, nq, nw, p
        real(dp) :: duration, total
        valid_camera = .false.
        if (cam%npixel < 1) return
        if (.not. present_text(cam%units)) return
        if (cam%units /= 'counts' .and. cam%units /= 'V') return
        if (.not. present_text(cam%geometry_version)) return
        if (.not. present_text(cam%calibration_version)) return
        if (.not. present_text(cam%source_version)) return
        if (.not. present_text(cam%thermal_version)) return
        if (.not. allocated(cam%exposure_edge_s)) return
        if (size(cam%exposure_edge_s, 1) /= 2) return
        n = size(cam%exposure_edge_s, 2)
        if (n < 1) return
        if (.not. all(ieee_is_finite(cam%exposure_edge_s))) return
        if (any(cam%exposure_edge_s(2, :) <= cam%exposure_edge_s(1, :))) return
        do e = 2, n
            if (cam%exposure_edge_s(1, e) < cam%exposure_edge_s(2, e - 1)) return
        end do
        if (.not. allocated(cam%quadrature_time_s)) return
        if (.not. allocated(cam%quadrature_weight_s)) return
        if (.not. allocated(cam%quadrature_exposure)) return
        nq = size(cam%quadrature_time_s)
        if (nq < n) return
        if (size(cam%quadrature_weight_s) /= nq) return
        if (size(cam%quadrature_exposure) /= nq) return
        if (.not. all(ieee_is_finite(cam%quadrature_time_s))) return
        if (.not. all(ieee_is_finite(cam%quadrature_weight_s))) return
        if (any(cam%quadrature_weight_s <= 0)) return
        if (any(cam%quadrature_exposure < 1) .or. &
            any(cam%quadrature_exposure > n)) return
        do e = 1, n
            total = 0
            do j = 1, nq
                if (cam%quadrature_exposure(j) /= e) cycle
                if (cam%quadrature_time_s(j) < cam%exposure_edge_s(1, e)) return
                if (cam%quadrature_time_s(j) > cam%exposure_edge_s(2, e)) return
                total = total + cam%quadrature_weight_s(j)
            end do
            duration = cam%exposure_edge_s(2, e) - cam%exposure_edge_s(1, e)
            if (abs(total - duration) > 1e-12_dp*duration) return
        end do
        ncell = size(temp, 1)
        if (ncell < 1 .or. size(temp, 2) /= nq) return
        if (.not. all(ieee_is_finite(temp))) return
        if (any(temp <= 0)) return
        if (.not. allocated(cam%wavelength_m)) return
        if (.not. allocated(cam%wavelength_weight_m)) return
        if (.not. allocated(cam%transmission)) return
        nw = size(cam%wavelength_m)
        if (nw < 1) return
        if (size(cam%wavelength_weight_m) /= nw) return
        if (size(cam%transmission) /= nw) return
        if (.not. all(ieee_is_finite(cam%wavelength_m))) return
        if (.not. all(ieee_is_finite(cam%wavelength_weight_m))) return
        if (.not. all(ieee_is_finite(cam%transmission))) return
        if (any(cam%wavelength_m <= tiny(1.0_dp)**0.2_dp)) return
        if (any(cam%wavelength_weight_m <= 0)) return
        if (any(cam%transmission < 0) .or. any(cam%transmission > 1)) return
        ! The emitted exp(-x)/(1-exp(-x)) form requires resolved subtraction.
        if (maxval(temp)*maxval(cam%wavelength_m) > 100.0_dp) return
        if (.not. allocated(cam%emissivity)) return
        if (.not. allocated(cam%background_W_m3sr)) return
        if (any(shape(cam%emissivity) /= [nw, ncell])) return
        if (any(shape(cam%background_W_m3sr) /= [nw, ncell])) return
        if (.not. all(ieee_is_finite(cam%emissivity))) return
        if (.not. all(ieee_is_finite(cam%background_W_m3sr))) return
        if (any(cam%emissivity < 0) .or. any(cam%emissivity > 1)) return
        if (any(cam%background_W_m3sr < 0)) return
        if (.not. allocated(cam%link_pixel)) return
        if (.not. allocated(cam%link_cell)) return
        if (.not. allocated(cam%etendue_m2sr)) return
        if (.not. allocated(cam%link_visible)) return
        l = size(cam%link_pixel)
        if (size(cam%link_cell) /= l) return
        if (size(cam%etendue_m2sr) /= l .or. size(cam%link_visible) /= l) return
        if (any(cam%link_pixel < 1) .or. any(cam%link_pixel > cam%npixel)) return
        if (any(cam%link_cell < 1) .or. any(cam%link_cell > ncell)) return
        if (.not. all(ieee_is_finite(cam%etendue_m2sr))) return
        if (any(cam%etendue_m2sr <= 0)) return
        do j = 2, l
            if (cam%link_pixel(j) < cam%link_pixel(j - 1)) return
        end do
        if (.not. allocated(cam%gain)) return
        if (.not. allocated(cam%offset)) return
        if (.not. allocated(cam%saturation)) return
        if (.not. allocated(cam%noise_variance)) return
        p = cam%npixel
        if (size(cam%gain) /= p .or. size(cam%offset) /= p) return
        if (size(cam%saturation) /= p) return
        if (any(shape(cam%noise_variance) /= [p, n])) return
        if (.not. all(ieee_is_finite(cam%gain))) return
        if (.not. all(ieee_is_finite(cam%offset))) return
        if (.not. all(ieee_is_finite(cam%saturation))) return
        if (.not. all(ieee_is_finite(cam%noise_variance))) return
        if (any(cam%gain <= 0)) return
        if (any(cam%saturation <= cam%offset)) return
        if (any(cam%noise_variance < 0)) return
        valid_camera = .true.
    end function

    subroutine pixel_offsets(cam, first)
        type(ir_camera_t), intent(in) :: cam
        integer, allocatable, intent(out) :: first(:)
        integer :: p, l
        allocate (first(cam%npixel + 1))
        l = 1
        do p = 1, cam%npixel
            first(p) = l
            do while (l <= size(cam%link_pixel))
                if (cam%link_pixel(l) /= p) exit
                l = l + 1
            end do
        end do
        first(cam%npixel + 1) = l
    end subroutine

    subroutine camera_integral(cam, temp, raw, visible)
        type(ir_camera_t), intent(in) :: cam
        real(dp), intent(in) :: temp(:, :)
        real(dp), allocatable, intent(out) :: raw(:, :)
        logical, allocatable, intent(out) :: visible(:)
        integer, allocatable :: first(:)
        integer :: p, l, j, w, c, e, n
        real(dp) :: radiance, dtemp, factor
        n = size(cam%exposure_edge_s, 2)
        allocate (raw(cam%npixel, n), source=0.0_dp)
        allocate (visible(cam%npixel), source=.false.)
        call pixel_offsets(cam, first)
        !$omp parallel do private(l,j,w,c,e,radiance,dtemp,factor)
        do p = 1, cam%npixel
            do l = first(p), first(p + 1) - 1
                if (.not. cam%link_visible(l)) cycle
                visible(p) = .true.
                c = cam%link_cell(l)
                do j = 1, size(cam%quadrature_time_s)
                    e = cam%quadrature_exposure(j)
                    factor = exposure_factor(cam, j)*cam%etendue_m2sr(l)
                    do w = 1, size(cam%wavelength_m)
                        call wall_radiance(temp(c, j), cam%wavelength_m(w), &
                            cam%emissivity(w, c), cam%background_W_m3sr(w, c), &
                            radiance, dtemp)
                        raw(p, e) = raw(p, e) + radiance*factor &
                            *cam%transmission(w)*cam%wavelength_weight_m(w)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine

    pure real(dp) function exposure_factor(cam, j) result(factor)
        type(ir_camera_t), intent(in) :: cam
        integer, intent(in) :: j
        integer :: e
        factor = cam%quadrature_weight_s(j)
        if (.not. cam%integrated) then
            e = cam%quadrature_exposure(j)
            factor = factor/(cam%exposure_edge_s(2, e) - cam%exposure_edge_s(1, e))
        end if
    end function

    subroutine synthesize_ir(camera, temperature, result, status)
        type(ir_camera_t), intent(in) :: camera
        real(dp), intent(in) :: temperature(:, :)
        type(wall_ir_result_t), intent(out) :: result
        integer, intent(out) :: status
        real(dp), allocatable :: raw(:, :)
        logical, allocatable :: visible(:)
        real(dp) :: dr, dg, doff, value
        integer :: p, e, n
        character(32) :: pixel
        status = WALL_ERROR
        if (.not. valid_camera(camera, temperature)) return
        call camera_integral(camera, temperature, raw, visible)
        if (.not. all(ieee_is_finite(raw))) return
        n = size(camera%exposure_edge_s, 2)
        allocate (result%channel(camera%npixel))
        allocate (result%valid(camera%npixel, n), source=.true.)
        allocate (result%flags(camera%npixel, n), source=IR_VALID)
        result%exposure_edge_s = camera%exposure_edge_s
        !$omp parallel do collapse(2) private(value,dr,dg,doff)
        do e = 1, n
            do p = 1, camera%npixel
                call recorded_affine(raw(p, e), camera%gain(p), camera%offset(p), &
                    value, dr, dg, doff)
                raw(p, e) = value
            end do
        end do
        !$omp end parallel do
        if (.not. all(ieee_is_finite(raw))) return
        do p = 1, camera%npixel
            allocate (result%channel(p)%samples(n))
            result%channel(p)%time = sum(camera%exposure_edge_s, dim=1)/2
            result%channel(p)%variance = camera%noise_variance(p, :)
            result%channel(p)%units = camera%units
            do e = 1, n
                result%channel(p)%samples(e) = raw(p, e)
                if (.not. visible(p)) then
                    result%flags(p, e) = IR_UNMAPPED
                else if (raw(p, e) >= camera%saturation(p)) then
                    result%flags(p, e) = IR_SATURATED
                end if
                result%valid(p, e) = result%flags(p, e) == IR_VALID
            end do
            write (pixel, '(I0)') p
            result%channel(p)%provenance%source_id = 'wall-ir:'//camera%source_version
            result%channel(p)%provenance%version = camera%calibration_version
            result%channel(p)%provenance%request_key = 'pixel='//trim(pixel)
            result%channel(p)%provenance%record_id = camera%geometry_version
            result%channel(p)%provenance%metadata = 'wall-ir-1; thermal=' &
                //camera%thermal_version//'; geometry='//camera%geometry_version &
                //'; validity=wall_ir_result_t%valid; conditional-fixed-geometry'
            call result%channel(p)%seal()
        end do
        status = WALL_OK
    end subroutine

    subroutine ir_jvp(camera, temperature, dtemperature, dgain, doffset, dsignal, status)
        type(ir_camera_t), intent(in) :: camera
        real(dp), intent(in) :: temperature(:, :), dtemperature(:, :), dgain(:), doffset(:)
        real(dp), allocatable, intent(out) :: dsignal(:, :)
        integer, intent(out) :: status
        type(wall_ir_result_t) :: result
        real(dp), allocatable :: raw(:, :)
        logical, allocatable :: visible(:)
        integer, allocatable :: first(:)
        integer :: p, l, j, w, c, e, n
        real(dp) :: radiance, dt, factor
        status = WALL_ERROR
        if (.not. valid_camera(camera, temperature)) return
        if (any(shape(dtemperature) /= shape(temperature))) return
        if (size(dgain) /= camera%npixel .or. size(doffset) /= camera%npixel) return
        if (.not. all(ieee_is_finite(dtemperature))) return
        if (.not. all(ieee_is_finite(dgain))) return
        if (.not. all(ieee_is_finite(doffset))) return
        call synthesize_ir(camera, temperature, result, status)
        if (status /= WALL_OK) return
        status = WALL_ERROR
        if (.not. allocated(result%valid)) return
        if (.not. all(result%valid)) return
        call camera_integral(camera, temperature, raw, visible)
        n = size(raw, 2)
        allocate (dsignal(camera%npixel, n), source=0.0_dp)
        call pixel_offsets(camera, first)
        !$omp parallel do private(l,j,w,c,e,radiance,dt,factor)
        do p = 1, camera%npixel
            dsignal(p, :) = dgain(p)*raw(p, :) + doffset(p)
            do l = first(p), first(p + 1) - 1
                if (.not. camera%link_visible(l)) cycle
                c = camera%link_cell(l)
                do j = 1, size(camera%quadrature_time_s)
                    e = camera%quadrature_exposure(j)
                    factor = exposure_factor(camera, j)*camera%etendue_m2sr(l)
                    do w = 1, size(camera%wavelength_m)
                        call wall_radiance(temperature(c, j), camera%wavelength_m(w), &
                            camera%emissivity(w, c), camera%background_W_m3sr(w, c), &
                            radiance, dt)
                        dsignal(p, e) = dsignal(p, e) + camera%gain(p)*factor &
                            *camera%transmission(w)*camera%wavelength_weight_m(w) &
                            *dt*dtemperature(c, j)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        if (.not. all(ieee_is_finite(dsignal))) return
        status = WALL_OK
    end subroutine

    subroutine ir_vjp(camera, temperature, bar_signal, bar_temperature, &
            bar_gain, bar_offset, status)
        type(ir_camera_t), intent(in) :: camera
        real(dp), intent(in) :: temperature(:, :), bar_signal(:, :)
        real(dp), allocatable, intent(out) :: bar_temperature(:, :), bar_gain(:), bar_offset(:)
        integer, intent(out) :: status
        type(wall_ir_result_t) :: result
        real(dp), allocatable :: raw(:, :)
        logical, allocatable :: visible(:)
        integer :: p, l, j, w, c, e
        real(dp) :: radiance, dt, factor
        status = WALL_ERROR
        if (.not. valid_camera(camera, temperature)) return
        if (any(shape(bar_signal) /= [camera%npixel, size(camera%exposure_edge_s, 2)])) return
        if (.not. all(ieee_is_finite(bar_signal))) return
        call synthesize_ir(camera, temperature, result, status)
        if (status /= WALL_OK) return
        status = WALL_ERROR
        if (.not. allocated(result%valid)) return
        if (.not. all(result%valid)) return
        call camera_integral(camera, temperature, raw, visible)
        allocate (bar_temperature(size(temperature, 1), size(temperature, 2)), source=0.0_dp)
        allocate (bar_gain(camera%npixel), bar_offset(camera%npixel))
        do p = 1, camera%npixel
            bar_gain(p) = sum(bar_signal(p, :)*raw(p, :))
            bar_offset(p) = sum(bar_signal(p, :))
        end do
        ! Each time column is independent; shared optical cells require accumulation.
        !$omp parallel do private(l,p,c,e,w,radiance,dt,factor)
        do j = 1, size(camera%quadrature_time_s)
            e = camera%quadrature_exposure(j)
            do l = 1, size(camera%link_pixel)
                if (.not. camera%link_visible(l)) cycle
                p = camera%link_pixel(l)
                c = camera%link_cell(l)
                factor = exposure_factor(camera, j)*camera%etendue_m2sr(l) &
                    *camera%gain(p)*bar_signal(p, e)
                do w = 1, size(camera%wavelength_m)
                    call wall_radiance(temperature(c, j), camera%wavelength_m(w), &
                        camera%emissivity(w, c), camera%background_W_m3sr(w, c), &
                        radiance, dt)
                    bar_temperature(c, j) = bar_temperature(c, j) + factor &
                        *camera%transmission(w)*camera%wavelength_weight_m(w)*dt
                end do
            end do
        end do
        !$omp end parallel do
        if (.not. all(ieee_is_finite(bar_temperature))) return
        if (.not. all(ieee_is_finite(bar_gain))) return
        if (.not. all(ieee_is_finite(bar_offset))) return
        status = WALL_OK
    end subroutine

    subroutine ir_valid_window(result, pixel, window, status)
        type(wall_ir_result_t), intent(in) :: result
        integer, intent(in) :: pixel
        type(signal_window_t), intent(out) :: window
        integer, intent(out) :: status
        integer, allocatable :: indices(:)
        integer :: j, n
        status = WALL_ERROR
        if (.not. allocated(result%channel)) return
        if (.not. allocated(result%valid)) return
        if (pixel < 1 .or. pixel > size(result%channel)) return
        if (size(result%valid, 1) /= size(result%channel)) return
        if (.not. result%channel(pixel)%valid()) return
        n = size(result%channel(pixel)%samples)
        if (size(result%valid, 2) /= n) return
        indices = pack([(j, j=1, n)], result%valid(pixel, :))
        if (size(indices) == 0) return
        window = result%channel(pixel)
        window%time = window%time(indices)
        window%samples = window%samples(indices)
        window%variance = window%variance(indices)
        if (allocated(window%covariance)) then
            window%covariance = window%covariance(indices, indices)
        end if
        window%provenance%metadata = window%provenance%metadata//'; invalid-exposures-removed'
        call window%seal()
        status = WALL_OK
    end subroutine
end module
