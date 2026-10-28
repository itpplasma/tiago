program test_wall_synthetic
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use tiago_wall_synthetic
    use tiago_wall_radiance, only: wall_radiance
    use tiago_signal_reduction, only: linear_reduction_t
    use tiago_signal_source, only: signal_window_t
    implicit none
    type(surface_deposition_t) :: s, rate, part, perturbed
    type(loss_event_block_t) :: events
    type(injection_history_t) :: history
    type(semi_infinite_material_t) :: material
    type(ir_camera_t) :: camera, cam2
    type(wall_ir_result_t) :: result, other
    type(signal_window_t) :: valid_window, reduced
    type(linear_reduction_t) :: reduction
    real(dp), allocatable :: flux(:, :), temp(:, :), tplus(:, :), tminus(:, :)
    real(dp), allocatable :: dt(:, :), denergy(:, :), be(:, :), bs(:, :), bt(:, :)
    real(dp), allocatable :: dg(:), doff(:), bg(:), bo(:), ds(:, :)
    real(dp) :: energy, expected, v, deriv, h, left, right
    integer :: status, i
    energy = 1.602176634e-19_dp
    call init_deposition(s, [0.0_dp, 1.0_dp, 3.0_dp], [2.0_dp, 4.0_dp], &
        'synthetic-two-facets', '1', status)
    call check(status == WALL_OK, 'init')
    events%time_s = [0.0_dp, 1.0_dp, 3.0_dp, 0.4_dp, 0.7_dp]
    events%kinetic_energy_eV = [2.0_dp, 3.0_dp, 5.0_dp, 1e20_dp, 1e20_dp]
    events%weight = [10.0_dp, 20.0_dp, 30.0_dp, 1.0_dp, 1.0_dp]
    events%deposited_fraction = [1.0_dp, 0.5_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    events%facet = [1, 2, 1, -1, -1]
    events%wall_hit = [1, 1, 1, 0, 1]
    events%exit_code = [2, 2, 2, 2, 3]
    events%normalization = PARTICLE_COUNT
    events%source_id = 'synthetic-ensemble'
    events%version = '1'
    events%frame = 'synthetic-Cartesian-metres-right-handed'
    events%position_kind = 'guiding-centre'
    call accumulate_losses(s, events, status)
    call check(status == WALL_OK, 'finite count accumulation with status filtering')
    call close (sum(s%energy_J)/energy, 200.0_dp, 1e-13_dp, 'conservative energy')
    call close (s%energy_J(1, 1)/energy, 20.0_dp, 1e-13_dp, 'left edge')
    call close (s%energy_J(1, 2)/energy, 150.0_dp, 1e-13_dp, 'last edge')
    call close (s%energy_J(2, 2)/energy, 30.0_dp, 1e-13_dp, 'interior edge')
    call deposition_flux(s, flux, status)
    call check(status == WALL_OK, 'flux')
    call close (flux(1, 2)/energy, 37.5_dp, 1e-13_dp, 'unequal-bin area-time units')
    part = s
    events%kinetic_energy_eV(2) = -1
    call accumulate_losses(s, events, status)
    call check(status == WALL_ERROR, 'reject bad event before mutation')
    call check(all(s%energy_J == part%energy_J), 'transactional invalid block')
    events%kinetic_energy_eV(2) = 3
    call init_deposition(rate, [0.0_dp, 1.0_dp, 3.0_dp], [2.0_dp], &
        'synthetic-one-facet', '1', status)
    deallocate (events%time_s, events%kinetic_energy_eV, events%weight)
    deallocate (events%deposited_fraction, events%facet, events%wall_hit, events%exit_code)
    events%time_s = [0.5_dp]
    events%kinetic_energy_eV = [10.0_dp]
    events%weight = [4.0_dp]
    events%deposited_fraction = [1.0_dp]
    events%facet = [1]
    events%wall_hit = [1]
    events%exit_code = [2]
    events%normalization = PARTICLE_RATE
    history%edge_s = [0.0_dp, 1.0_dp, 2.0_dp]
    history%amplitude = [1.0_dp, 2.0_dp]
    call accumulate_losses(rate, events, status)
    call check(status == WALL_ERROR, 'rate requires injection history')
    call accumulate_losses(rate, events, status, history)
    call check(status == WALL_OK, 'source convolution')
    call close (rate%energy_J(1, 1)/energy, 20.0_dp, 1e-13_dp, 'rate overlap first bin')
    call close (rate%energy_J(1, 2)/energy, 100.0_dp, 1e-13_dp, 'rate not divided twice')
    call init_deposition(part, [0.0_dp, 1.0_dp, 3.0_dp], [2.0_dp], &
        'synthetic-one-facet', '1', status)
    events%normalization = PARTICLE_COUNT
    call accumulate_losses(part, events, status, history)
    call check(status == WALL_ERROR, 'counts reject source convolution')
    call accumulate_losses(part, events, status)
    call close (sum(part%energy_J)/energy, 40.0_dp, 1e-13_dp, 'counts distinct from rates')

    s%energy_J = reshape([4000.0_dp, 8000.0_dp, 8000.0_dp, 16000.0_dp], [2, 2])
    material%conductivity_W_mK = [10.0_dp, 20.0_dp]
    material%capacity_J_m3K = [1e6_dp, 2e6_dp]
    material%initial_K = [300.0_dp, 350.0_dp]
    material%regime = SEMI_INFINITE_REGIME
    material%version = 'synthetic-material-1'
    call surface_temperature(s, material, [0.0_dp, 1.0_dp, 3.0_dp, 4.0_dp], temp, status)
    call check(status == WALL_OK, 'thermal')
    call close (temp(1, 1), 300.0_dp, 1e-14_dp, 'initial state')
    expected = 300 + 2*2000.0_dp*sqrt(3.0_dp)/sqrt(acos(-1.0_dp)*1e7_dp)
    call close (temp(1, 3), expected, 1e-13_dp, 'analytic constant flux sqrt law')
    expected = 300 + 2*2000.0_dp*(2 - 1)/sqrt(acos(-1.0_dp)*1e7_dp)
    call close (temp(1, 4), expected, 1e-13_dp, 'analytic cooling after flux pulse')
    material%regime = 'actively-cooled'
    call surface_temperature(s, material, [1.0_dp], temp, status)
    call check(status == WALL_ERROR, 'reject unstated thermal regime')
    material%regime = SEMI_INFINITE_REGIME
    denergy = reshape([0.3_dp, -0.2_dp, 0.8_dp, 0.5_dp], [2, 2])
    call temperature_energy_jvp(s, material, [0.2_dp, 1.7_dp, 4.0_dp], denergy, dt, status)
    call check(status == WALL_OK, 'thermal JVP')
    bs = reshape([0.2_dp, 0.7_dp, -0.3_dp, 0.2_dp, 0.4_dp, -0.1_dp], [2, 3])
    call temperature_energy_vjp(s, material, [0.2_dp, 1.7_dp, 4.0_dp], bs, be, status)
    call close (sum(bs*dt), sum(be*denergy), 1e-13_dp, 'thermal transpose')
    perturbed = s
    h = 1e-2_dp
    perturbed%energy_J = s%energy_J + h*denergy
    call surface_temperature(perturbed, material, [0.2_dp, 1.7_dp, 4.0_dp], tplus, status)
    perturbed%energy_J = s%energy_J - h*denergy
    call surface_temperature(perturbed, material, [0.2_dp, 1.7_dp, 4.0_dp], tminus, status)
    call check(maxval(abs((tplus - tminus)/(2*h) - dt)) < 3e-12_dp, 'thermal FD')

    call wall_radiance(500.0_dp, 1e-5_dp, 1.0_dp, 0.0_dp, v, deriv)
    expected = 2*6.62607015e-34_dp*299792458.0_dp**2 &
        /(1e-25_dp*(exp(6.62607015e-34_dp*299792458.0_dp &
        /(1e-5_dp*1.380649e-23_dp*500)) - 1))
    call close (v, expected, 2e-14_dp, 'independent Planck oracle')
    call camera_fixture(camera)
    call surface_temperature(s, material, camera%quadrature_time_s, temp, status)
    call synthesize_ir(camera, temp, result, status)
    call check(status == WALL_OK, 'complete raw IR path')
    call check(all(result%valid), 'valid raw windows')
    call check(result%channel(1)%valid(), 'existing signal window contract')
    call check(result%channel(1)%has_provenance(), 'signal provenance')
    call check(result%channel(1)%units == 'counts', 'raw counts units')
    result%channel(1)%covariance = reshape([0.1_dp, 0.02_dp, 0.02_dp, 0.2_dp], [2, 2])
    call result%channel(1)%seal()
    call ir_valid_window(result, 1, valid_window, status)
    call check(status == WALL_OK, 'valid window adapter')
    reduction%id = 'synthetic-exposure-average'
    reduction%units = 'counts'
    reduction%weights = reshape([0.5_dp, 0.5_dp], [1, 2])
    reduction%noise_covariance = reshape([0.0_dp], [1, 1])
    call reduction%apply(valid_window, reduced, status)
    call check(status == WALL_OK, 'existing linear reduction')
    call close (reduced%covariance(1, 1), 0.21_dp, 1e-13_dp, &
        'existing correlated covariance')
    result%valid(1, 2) = .false.
    call ir_valid_window(result, 1, valid_window, status)
    call check(size(valid_window%samples) == 1, 'validity filtering')
    call close (valid_window%covariance(1, 1), 0.1_dp, 1e-14_dp, 'masked covariance')
    cam2 = camera
    cam2%emissivity = 0
    cam2%background_W_m3sr = 2.0_dp
    call synthesize_ir(cam2, temp, other, status)
    ! Constant radiance: 2 * bandwidth 3e-6 * etendue 0.5 * duration * gain + offset.
    call close (other%channel(1)%samples(1), 7 + 3e-6_dp*3*2, 1e-13_dp, 'exposure integral')
    call close (other%channel(1)%samples(2), 7 + 3e-6_dp*1*2, 1e-13_dp, 'unequal exposure')
    cam2%integrated = .false.
    cam2%units = 'V'
    call synthesize_ir(cam2, temp, other, status)
    call close (other%channel(1)%samples(1), 7 + 3e-6_dp*2, 1e-13_dp, 'voltage mean power')
    call close (other%channel(1)%samples(2), other%channel(1)%samples(1), &
        1e-14_dp, 'average does not depend on duration')
    cam2 = camera
    cam2%link_visible(2) = .false.
    call synthesize_ir(cam2, temp, other, status)
    call check(all(other%flags(2, :) == IR_UNMAPPED), 'occlusion invalid')
    cam2 = camera
    cam2%saturation(1) = 7.0000000001_dp
    call synthesize_ir(cam2, temp, other, status)
    call check(all(other%flags(1, :) == IR_SATURATED), 'saturation invalid')
    cam2 = camera
    deallocate (cam2%calibration_version)
    call synthesize_ir(cam2, temp, other, status)
    call check(status == WALL_ERROR, 'missing calibration rejected')
    cam2 = camera
    cam2%quadrature_weight_s(1) = 1
    call synthesize_ir(cam2, temp, other, status)
    call check(status == WALL_ERROR, 'bad exposure quadrature rejected')
    cam2 = camera
    cam2%gain(1) = ieee_value(1.0_dp, ieee_quiet_nan)
    call synthesize_ir(cam2, temp, other, status)
    call check(status == WALL_ERROR, 'NaN calibration rejected')
    dt = reshape([0.2_dp, -0.1_dp, 0.4_dp, 0.6_dp], [2, 2])
    dg = [0.03_dp, -0.01_dp]
    doff = [0.2_dp, -0.3_dp]
    call ir_jvp(camera, temp, dt, dg, doff, ds, status)
    call check(status == WALL_OK, 'IR JVP')
    bs = reshape([0.2_dp, -0.1_dp, 0.4_dp, 0.6_dp], [2, 2])
    call ir_vjp(camera, temp, bs, bt, bg, bo, status)
    call check(status == WALL_OK, 'IR VJP')
    left = sum(bs*ds)
    right = sum(bt*dt) + sum(bg*dg) + sum(bo*doff)
    call close (left, right, 1e-12_dp, 'IR transpose')
    h = 1e-3_dp
    cam2 = camera
    cam2%gain = camera%gain + h*dg
    cam2%offset = camera%offset + h*doff
    call synthesize_ir(cam2, temp + h*dt, result, status)
    cam2%gain = camera%gain - h*dg
    cam2%offset = camera%offset - h*doff
    call synthesize_ir(cam2, temp - h*dt, other, status)
    do i = 1, 2
        call check(maxval(abs((result%channel(i)%samples - other%channel(i)%samples) &
            /(2*h) - ds(i, :))) < 1e-9_dp, 'IR directional FD')
    end do
    cam2%saturation(1) = 7.0000000001_dp
    call ir_jvp(cam2, temp, dt, dg, doff, ds, status)
    call check(status == WALL_ERROR, 'reject derivative across saturated state')
    print *, 'wall synthetic independent energy/thermal/Planck/exposure/adjoint checks PASS'
contains
    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print *, 'FAIL: ', label
            error stop 1
        end if
    end subroutine
    subroutine close (actual, oracle, tolerance, label)
        real(dp), intent(in) :: actual, oracle, tolerance
        character(*), intent(in) :: label
        call check(abs(actual - oracle) <= tolerance*max(abs(oracle), 1e-20_dp), label)
    end subroutine
    subroutine camera_fixture(cam)
        type(ir_camera_t), intent(out) :: cam
        cam%npixel = 2
        cam%exposure_edge_s = reshape([0.0_dp, 3.0_dp, 3.0_dp, 4.0_dp], [2, 2])
        cam%quadrature_time_s = [1.5_dp, 3.5_dp]
        cam%quadrature_weight_s = [3.0_dp, 1.0_dp]
        cam%quadrature_exposure = [1, 2]
        cam%wavelength_m = [8e-6_dp, 1e-5_dp]
        cam%wavelength_weight_m = [1e-6_dp, 2e-6_dp]
        cam%transmission = [1.0_dp, 1.0_dp]
        allocate (cam%emissivity(2, 2), source=0.8_dp)
        allocate (cam%background_W_m3sr(2, 2), source=1.0_dp)
        cam%link_pixel = [1, 2]
        cam%link_cell = [1, 2]
        cam%etendue_m2sr = [0.5_dp, 0.25_dp]
        cam%link_visible = [.true., .true.]
        cam%gain = [2.0_dp, 3.0_dp]
        cam%offset = [7.0_dp, -2.0_dp]
        cam%saturation = [1e10_dp, 1e10_dp]
        allocate (cam%noise_variance(2, 2), source=0.25_dp)
        cam%units = 'counts'
        cam%geometry_version = 'synthetic-geometry-1'
        cam%calibration_version = 'synthetic-calibration-1'
        cam%source_version = 'synthetic-ensemble-1'
        cam%thermal_version = 'synthetic-material-1'
    end subroutine
end program
