program test_simple_wall_source
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_wall_synthetic
    use tiago_simple_wall_source, only: accumulate_simple_netcdf
    implicit none
    type(surface_deposition_t) :: surface, legacy, unknown
    type(injection_history_t) :: history
    character(4096) :: dir
    integer :: status
    real(dp), parameter :: ev = 1.602176634e-19_dp
    call get_command_argument(1, dir)
    call init_deposition(surface, [0.0_dp, 1.0_dp, 3.0_dp], [2.0_dp, 4.0_dp], &
        'synthetic-two-triangles', '1', status)
    legacy = surface
    unknown = surface
    call accumulate_simple_netcdf(trim(dir)//'/events.nc', surface, 100.0_dp, &
        PARTICLE_COUNT, 'synthetic-source', '1', 'synthetic-Cartesian', &
        [1.0_dp, 0.5_dp], status, block_size=1)
    if (status /= WALL_OK) error stop 'read explicit energy/weights'
    if (abs(surface%energy_J(1, 1)/ev - 250) > 1e-10_dp) error stop 'all-source denominator'
    if (abs(surface%energy_J(2, 2)/ev - 960) > 1e-10_dp) &
        error stop 'energy/absorption/facet'
    call accumulate_simple_netcdf(trim(dir)//'/legacy.nc', legacy, 100.0_dp, &
        PARTICLE_COUNT, 'synthetic-source', '1', 'synthetic-Cartesian', &
        [1.0_dp, 0.5_dp], status, &
        external_weight=[1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp], block_size=3)
    if (status /= WALL_OK) error stop 'legacy E0*p2 with external weights'
    if (maxval(abs(surface%energy_J - legacy%energy_J)) > 1e-28_dp) &
        error stop 'legacy equivalence'
    call accumulate_simple_netcdf(trim(dir)//'/unknown.nc', unknown, 100.0_dp, &
        PARTICLE_COUNT, 'synthetic-source', '1', 'synthetic-Cartesian', &
        [1.0_dp, 0.5_dp], status)
    if (status /= WALL_ERROR) error stop 'unknown source sampling must reject'
    if (any(unknown%energy_J /= 0)) error stop 'failed file mutated deposition'
    call accumulate_simple_netcdf(trim(dir)//'/unknown.nc', unknown, 100.0_dp, &
        PARTICLE_COUNT, 'synthetic-source', '1', 'synthetic-Cartesian', &
        [1.0_dp, 0.5_dp], status, external_weight=[1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp])
    if (status /= WALL_OK) error stop 'external weights override unknown sampling'
    history%edge_s = [0.0_dp, 0.5_dp]
    history%amplitude = [1.0_dp]
    call init_deposition(unknown, [0.0_dp, 1.0_dp, 3.0_dp], [2.0_dp, 4.0_dp], &
        'synthetic-two-triangles', '1', status)
    call accumulate_simple_netcdf(trim(dir)//'/events.nc', unknown, 100.0_dp, &
        PARTICLE_RATE, 'synthetic-source', '1', 'synthetic-Cartesian', &
        [1.0_dp, 0.5_dp], status, history=history, block_size=2)
    if (status /= WALL_OK) error stop 'rate NetCDF integration'
    if (abs(sum(unknown%energy_J)/ev - 605) > 1e-10_dp) &
        error stop 'particle rate normalization'
    call init_deposition(unknown, [0.0_dp, 1.0_dp, 3.0_dp], [2.0_dp, 4.0_dp], &
        'synthetic-two-triangles', '1', status)
    call accumulate_simple_netcdf(trim(dir)//'/bad-late.nc', unknown, 100.0_dp, &
        PARTICLE_COUNT, 'synthetic-source', '1', 'synthetic-Cartesian', &
        [1.0_dp, 0.5_dp], status, block_size=1)
    if (status /= WALL_ERROR) error stop 'invalid later block must reject'
    if (any(unknown%energy_J /= 0)) error stop 'later block mutated caller surface'
    call accumulate_simple_netcdf(trim(dir)//'/events.nc', unknown, 100.0_dp, &
        PARTICLE_COUNT, 'synthetic-source', '1', 'wrong-frame', &
        [1.0_dp, 0.5_dp], status)
    if (status /= WALL_ERROR) error stop 'wrong Cartesian frame'
    print *, 'SIMPLE-shaped NetCDF block/legacy/normalization checks PASS'
end program
