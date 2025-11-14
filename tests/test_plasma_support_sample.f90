program test_plasma_support_sample
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_plasma_support, only: plasma_support_t
    implicit none

    type(plasma_support_t) :: plasma
    real(dp), allocatable :: points(:, :)
    real(dp), allocatable :: bfield(:, :)
    real(dp), allocatable :: reference(:, :)
    real(dp), allocatable :: bfield_surface(:, :)
    real(dp), allocatable :: avec(:, :)
    real(dp) :: max_err
    real(dp) :: t0
    real(dp) :: t1
    real(dp) :: t2
    character(len=512) :: wout_path
    character(len=32) :: flag
    integer :: nargs
    logical :: dump_only
    character(len=*), parameter :: reference_file = &
                                   'tests/data/diagno_bfield_reference.dat'
    integer, parameter :: sample_nphi = 4
    integer, parameter :: sample_ntheta = 4

    dump_only = .false.
    flag = ''
    nargs = command_argument_count()
    if (nargs < 1) then
        write (*, '(A)') 'Usage: test_plasma_support_sample <wout_file> [--dump]'
        stop 1
    end if
    call get_command_argument(1, wout_path)
    if (nargs >= 2) then
        call get_command_argument(2, flag)
        dump_only = trim(flag) == '--dump'
    end if

    if (.not. file_exists(wout_path)) then
        write (*, '(A)') 'specified VMEC file not found'
        stop 1
    end if
    write (*, '(A)') 'using wout file: '//trim(wout_path)

    call cpu_time(t0)
    call plasma%init_from_vmec(trim(wout_path), sample_nphi, sample_ntheta)
    call cpu_time(t1)
    write (*, '(A,1pe12.5)') 'plasma init seconds=', t1 - t0
    if (.not. plasma%has_data()) error stop 'plasma context failed to initialize'

    call load_reference(points, reference)
    allocate (bfield(size(points, 1), 3))
    allocate (bfield_surface(size(points, 1), 3))
    allocate (avec(size(points, 1), 3))

    write (*, '(A)') 'sampling Biot-Savart field for diagnostics...'
    call plasma%sample_surface_bfield(points, bfield)
    call cpu_time(t2)
    write (*, '(A,1pe12.5)') 'sample seconds=', t2 - t1

    if (dump_only) then
        call dump_samples(points, bfield)
        call plasma%finalize()
        stop
    end if

    bfield_surface = bfield
    write (*, '(A)') 'sampling vector potential...'
    call plasma%sample_vector_potential(points, avec)

    max_err = maxval(abs(bfield - reference))
    if (max_err > 5.0e-6_dp) then
        write (*, '(A,1pe12.5)') 'ERROR: plasma sample mismatch; max error=', max_err
        error stop 1
    end if

    call verify_vector_potential(plasma, points, bfield)

    call plasma%finalize()

contains

    logical function file_exists(path)
        character(len=*), intent(in) :: path
        inquire (file=trim(path), exist=file_exists)
    end function file_exists

    subroutine load_reference(points, ref_values)
        real(dp), allocatable, intent(out) :: points(:, :)
        real(dp), allocatable, intent(out) :: ref_values(:, :)
        character(len=512) :: candidate
        integer :: unit
        integer :: npts
        integer :: idx

        candidate = reference_file
        if (.not. file_exists(candidate)) then
            candidate = '../'//reference_file
            if (.not. file_exists(candidate)) then
                write (*, '(A)') 'unable to locate DIAGNO reference dataset'
                error stop 1
            end if
        end if

        open (newunit=unit, file=trim(candidate), status='old', action='read')
        read (unit, *) npts
        if (npts <= 0) then
            write (*, '(A)') 'invalid DIAGNO reference dataset'
            close (unit)
            error stop 1
        end if
        allocate (points(npts, 3))
        allocate (ref_values(npts, 3))
        do idx = 1, npts
            read (unit, *) points(idx, 1), points(idx, 2), points(idx, 3), &
                ref_values(idx, 1), ref_values(idx, 2), ref_values(idx, 3)
        end do
        close (unit)
    end subroutine load_reference

    subroutine dump_samples(points, bfield)
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(in) :: bfield(:, :)
        integer :: idx

        do idx = 1, size(points, 1)
            write (*, '(A,I0)') 'sample ', idx
            write (*, '(3(1pe20.12))') points(idx, :)
            write (*, '(3(1pe20.12))') bfield(idx, :)
        end do
    end subroutine dump_samples

    subroutine verify_vector_potential(plasma, points, target_b)
        type(plasma_support_t), intent(in) :: plasma
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(in) :: target_b(:, :)
        real(dp), allocatable :: approx(:, :)
        real(dp) :: delta
        real(dp) :: max_err_local
        integer :: idx

        delta = 2.5e-3_dp
        allocate(approx(size(points, 1), 3))
        do idx = 1, size(points, 1)
            call curl_from_vector_potential(plasma, points(idx, :), delta, approx(idx, :))
        end do
        max_err_local = maxval(abs(approx - target_b))
        if (max_err_local > 5.0e-2_dp) then
            write (*, '(A,1pe12.5)') 'ERROR: vector potential curl mismatch; max error=', max_err_local
            error stop 1
        end if
        deallocate(approx)
    end subroutine verify_vector_potential

    subroutine curl_from_vector_potential(plasma, point, delta, curl_value)
        type(plasma_support_t), intent(in) :: plasma
        real(dp), intent(in) :: point(3)
        real(dp), intent(in) :: delta
        real(dp), intent(out) :: curl_value(3)

        real(dp) :: eval_points(6, 3)
        real(dp) :: a_eval(6, 3)
        real(dp) :: inv_twodelta

        inv_twodelta = 1.0_dp / (2.0_dp * delta)
        eval_points(1, :) = point
        eval_points(1, 3) = eval_points(1, 3) + delta
        eval_points(2, :) = point
        eval_points(2, 3) = eval_points(2, 3) - delta
        eval_points(3, :) = point
        eval_points(3, 1) = eval_points(3, 1) + delta
        eval_points(4, :) = point
        eval_points(4, 1) = eval_points(4, 1) - delta
        eval_points(5, :) = point
        eval_points(5, 2) = eval_points(5, 2) + delta
        eval_points(6, :) = point
        eval_points(6, 2) = eval_points(6, 2) - delta

        call plasma%sample_vector_potential(eval_points, a_eval)

        curl_value(1) = (a_eval(5, 3) - a_eval(6, 3)) * inv_twodelta - &
            (a_eval(1, 2) - a_eval(2, 2)) * inv_twodelta
        curl_value(2) = (a_eval(1, 1) - a_eval(2, 1)) * inv_twodelta - &
            (a_eval(3, 3) - a_eval(4, 3)) * inv_twodelta
        curl_value(3) = (a_eval(3, 2) - a_eval(4, 2)) * inv_twodelta - &
            (a_eval(5, 1) - a_eval(6, 1)) * inv_twodelta
    end subroutine curl_from_vector_potential

end program test_plasma_support_sample
