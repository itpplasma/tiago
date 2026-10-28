program benchmark_wall
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use tiago_wall_synthetic
    implicit none
    type(surface_deposition_t) :: surface
    type(semi_infinite_material_t) :: material
    type(ir_camera_t) :: camera
    type(wall_ir_result_t) :: result
    real(dp), allocatable :: temp(:, :), edge(:), area(:)
    integer, parameter :: nc = 2048, nb = 64, nq = 128, nw = 8
    integer :: i, j, status, rep
    integer(int64) :: tick0, tick1, rate, phase0, phase1, thermal_ticks, ir_ticks
    allocate (edge(nb + 1), area(nc), source=1.0_dp)
    edge = [(real(i - 1, dp)/nb, i=1, nb + 1)]
    call init_deposition(surface, edge, area, 'synthetic-benchmark', '1', status)
    surface%energy_J = 100
    allocate (material%conductivity_W_mK(nc), source=100.0_dp)
    allocate (material%capacity_J_m3K(nc), source=1e6_dp)
    allocate (material%initial_K(nc), source=500.0_dp)
    material%regime = SEMI_INFINITE_REGIME
    material%version = 'synthetic'
    camera%npixel = nc
    allocate (camera%exposure_edge_s(2, nq))
    camera%exposure_edge_s(1, :) = [(real(i - 1, dp)/nq, i=1, nq)]
    camera%exposure_edge_s(2, :) = [(real(i, dp)/nq, i=1, nq)]
    camera%quadrature_time_s = [(real(i, dp)/nq, i=1, nq)]
    allocate (camera%quadrature_weight_s(nq), source=1.0_dp/nq)
    camera%quadrature_exposure = [(i, i=1, nq)]
    camera%wavelength_m = [(7e-6_dp + real(i, dp)*5e-7_dp, i=1, nw)]
    allocate (camera%wavelength_weight_m(nw), source=5e-7_dp)
    allocate (camera%transmission(nw), source=0.9_dp)
    allocate (camera%emissivity(nw, nc), source=0.8_dp)
    allocate (camera%background_W_m3sr(nw, nc), source=1.0_dp)
    camera%link_pixel = [(i, i=1, nc)]
    camera%link_cell = [(i, i=1, nc)]
    allocate (camera%etendue_m2sr(nc), source=1e-6_dp)
    allocate (camera%link_visible(nc), source=.true.)
    allocate (camera%gain(nc), source=1e6_dp)
    allocate (camera%offset(nc), source=1.0_dp)
    allocate (camera%saturation(nc), source=1e20_dp)
    allocate (camera%noise_variance(nc, nq), source=1.0_dp)
    camera%units = 'counts'
    camera%geometry_version = 'synthetic'
    camera%calibration_version = 'synthetic'
    camera%source_version = 'synthetic'
    camera%thermal_version = 'synthetic'
    call system_clock(tick0, rate)
    thermal_ticks = 0
    ir_ticks = 0
    do rep = 1, 4
        call system_clock(phase0)
        call surface_temperature(surface, material, camera%quadrature_time_s, &
            temp, status)
        if (status /= WALL_OK) error stop 'benchmark thermal'
        call system_clock(phase1)
        thermal_ticks = thermal_ticks + phase1 - phase0
        call system_clock(phase0)
        call synthesize_ir(camera, temp, result, status)
        if (status /= WALL_OK) error stop 'benchmark IR'
        call system_clock(phase1)
        ir_ticks = ir_ticks + phase1 - phase0
    end do
    call system_clock(tick1)
    print '(A,F12.6,A,F12.6)', 'thermal_seconds=', real(thermal_ticks, dp)/rate, &
        ' IR_seconds=', real(ir_ticks, dp)/rate
    print '(A,F12.6,A,ES24.16)', 'seconds=', real(tick1 - tick0, dp)/real(rate, dp), &
        ' checksum=', sum([(sum(result%channel(j)%samples), j=1, nc)])
end program
