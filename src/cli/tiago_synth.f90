program tiago_synth
    !! Synthetic diagnostics from one namelist; see docs/synthetic-forward.md.
    !! Usage: tiago_synth <input.nml>
    use, intrinsic :: iso_fortran_env, only: error_unit
    use tiago_synthetic_forward, only: synth_config_t, synth_run
    implicit none

    character(len=4096) :: wout = '', output = 'synthetic.nc', profiles = '', geometry = '', &
        sensors = '', coils = '', extcur = ''
    character(len=32) :: plasma_sheet = 'contravariant'
    integer :: npanel = 200, thomson_points = 4, ece_points = 4000, profile_knots = 12, &
        plasma_nphi = 64, plasma_ntheta = 64
    namelist /synth/ wout, output, profiles, geometry, sensors, coils, extcur, plasma_sheet, &
        npanel, thomson_points, ece_points, profile_knots, plasma_nphi, plasma_ntheta
    character(len=4096) :: path
    character(len=:), allocatable :: message
    type(synth_config_t) :: cfg
    integer :: unit, ios, ierr

    if (command_argument_count() /= 1) call die('usage: tiago_synth <input.nml>')
    call get_command_argument(1, path)
    open(newunit=unit, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) call die('cannot open '//trim(path))
    read(unit, nml=synth, iostat=ios)
    close(unit)
    if (ios /= 0) call die(trim(path)//': cannot read namelist &synth')
    if (len_trim(wout) == 0) call die(trim(path)//': wout is required')
    select case (trim(plasma_sheet))
    case ('contravariant', 'covariant', 'conservative')
    case default
        call die(trim(path)//': plasma_sheet must be contravariant, covariant or conservative')
    end select
    if (npanel < 1 .or. thomson_points < 1 .or. ece_points < 2 .or. profile_knots < 1) &
        call die(trim(path)//': npanel, thomson_points, profile_knots >= 1 and ece_points >= 2')
    cfg%wout = trim(wout)
    cfg%output = trim(output)
    cfg%profiles = trim(profiles)
    cfg%geometry = trim(geometry)
    cfg%sensors = trim(sensors)
    cfg%coils = trim(coils)
    cfg%extcur = trim(extcur)
    cfg%sheet = trim(plasma_sheet)
    cfg%npanel = npanel
    cfg%thomson_points = thomson_points
    cfg%ece_points = ece_points
    cfg%profile_knots = profile_knots
    cfg%plasma_nphi = plasma_nphi
    cfg%plasma_ntheta = plasma_ntheta
    call synth_run(cfg, ierr, message)
    if (ierr /= 0) call die(message)
    print '(A)', 'wrote '//cfg%output

contains

    subroutine die(text)
        character(len=*), intent(in) :: text
        write(error_unit, '(A)') 'tiago_synth: '//text
        error stop 1
    end subroutine die
end program tiago_synth
