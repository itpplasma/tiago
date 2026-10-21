module tiago_mse
    !! Point MSE pi-polarization; definitions and derivation: docs/mse.md.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, &
        ieee_value, ieee_quiet_nan
    use tiago_diagnostic_types, only: mse_channel_t
    use tiago_mse_rows, only: mse_rows
    use tiago_mse_local, only: mse_local
    implicit none
    private
    public :: read_mse_file, mse_geometry, mse_angles
    real(dp), parameter :: pi = acos(-1.0_dp)
contains
    subroutine read_mse_file(path, channels, ierr, message)
        character(len=*), intent(in) :: path
        type(mse_channel_t), allocatable, intent(out) :: channels(:)
        integer(i32), intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        type(mse_channel_t), allocatable :: parsed(:)
        integer :: unit, ios, count, i, j
        character(len=4096) :: line
        character(len=128) :: label, extra
        real(dp) :: row(16), point(3), coefficients(6)

        ierr = 1_i32
        message = 'unable to open MSE file: '//trim(path)
        open (newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        call data_line(unit, line, ios)
        count = 0
        if (ios == 0) read (line, *, iostat=ios) count
        message = 'invalid MSE channel count: '//trim(path)
        if (ios /= 0 .or. count <= 0) then
            close (unit)
            return
        end if
        allocate (parsed(count))
        do i = 1, count
            call data_line(unit, line, ios)
            row = ieee_value(0.0_dp, ieee_quiet_nan)
            label = ''
            if (ios == 0) read (line, *, iostat=ios) label, row
            message = 'malformed MSE row '//int_string(i)//' in '//trim(path)
            if (ios /= 0) exit
            if (.not. all(ieee_is_finite(row))) exit
            read (line, *, iostat=ios) label, row, extra
            if (ios >= 0) exit ! Reject extra columns.
            if (len_trim(label) == 0 .or. len_trim(label) >= len(label)) exit
            if (verify(trim(label), &
                'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-') &
                /= 0) exit
            do j = 1, i - 1
                if (parsed(j)%label == trim(label)) exit
            end do
            if (j < i) then
                message = 'duplicate MSE label: '//trim(label)
                exit
            end if
            parsed(i)%label = trim(label)
            parsed(i)%beam_point = row(1:3)
            parsed(i)%beam_direction = row(4:6)
            parsed(i)%sight_point = row(7:9)
            parsed(i)%sight_direction = row(10:12)
            parsed(i)%reference = row(13:15)
            parsed(i)%max_gap = row(16)
            call mse_geometry(parsed(i), point, coefficients, ierr, message)
            if (ierr /= 0) then
                message = trim(label)//': '//message
                exit
            end if
        end do
        if (i > count) then
            call data_line(unit, line, ios)
            if (ios < 0) then
                call move_alloc(parsed, channels)
                ierr = 0_i32
                message = ''
            else
                ierr = 1_i32
                message = 'extra MSE rows in '//trim(path)
            end if
        else
            ierr = 1_i32
        end if
        close (unit)
    end subroutine read_mse_file

    subroutine mse_geometry(channel, point, coefficients, ierr, message, gap, rows)
        type(mse_channel_t), intent(in) :: channel
        real(dp), intent(out) :: point(3), coefficients(6)
        integer(i32), intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        real(dp), intent(out), optional :: gap, rows(2, 3)
        real(dp) :: v(3), n(3), zero(3), h(3), w(3), a(3), b(3)
        real(dp) :: c, determinant, t, u, pb(3), ps(3), distance, radius
        real(dp) :: er(3), ep(3)
        logical :: ok

        ierr = 1_i32
        message = 'nonfinite MSE geometry or negative max_gap'
        point = 0.0_dp
        coefficients = 0.0_dp
        if (.not. all(ieee_is_finite(channel%beam_point))) return
        if (.not. all(ieee_is_finite(channel%sight_point))) return
        if (.not. ieee_is_finite(channel%max_gap)) return
        if (channel%max_gap < 0.0_dp) return
        message = 'zero or nonfinite MSE direction/reference'
        call unit_vector(channel%beam_direction, v, ok)
        if (.not. ok) return
        call unit_vector(channel%sight_direction, n, ok)
        if (.not. ok) return
        call unit_vector(channel%reference, zero, ok)
        if (.not. ok) return
        zero = zero - dot_product(zero, n)*n
        message = 'MSE reference parallel to sightline'
        if (norm2(zero) <= 1.0e-12_dp) return
        zero = zero/norm2(zero)
        h = cross(n, zero)
        c = dot_product(v, n)
        determinant = sum(cross(v, n)**2)
        message = 'parallel or nearly parallel MSE lines'
        if (determinant <= 1.0e-12_dp) return
        w = channel%sight_point - channel%beam_point
        t = (dot_product(v, w) - c*dot_product(n, w))/determinant
        u = (c*dot_product(v, w) - dot_product(n, w))/determinant
        pb = channel%beam_point + t*v
        ps = channel%sight_point + u*n
        point = 0.5_dp*pb + 0.5_dp*ps
        distance = norm2(pb - ps)
        if (present(gap)) gap = distance
        message = 'MSE closest-approach gap exceeds max_gap'
        if (.not. ieee_is_finite(distance)) return
        if (distance > channel%max_gap) return
        message = 'MSE intersection is nonfinite or on cylindrical axis'
        if (.not. all(ieee_is_finite(point))) return
        radius = norm2(point(1:2))
        if (radius <= 64*epsilon(radius)*max(1.0_dp, &
            norm2(channel%beam_point), norm2(channel%sight_point))) return
        er = [point(1)/radius, point(2)/radius, 0.0_dp]
        ep = [-er(2), er(1), 0.0_dp]
        call mse_rows(v(1), v(2), v(3), h(1), h(2), h(3), &
            zero(1), zero(2), zero(3), a(1), a(2), a(3), b(1), b(2), b(3))
        coefficients = [a(3), dot_product(a, er), dot_product(a, ep), &
            b(3), dot_product(b, er), dot_product(b, ep)]
        if (present(rows)) then
            rows(1, :) = a
            rows(2, :) = b
        end if
        ierr = 0_i32
        message = ''
    end subroutine mse_geometry

    subroutine mse_angles(channels, field, angles, ierr, message, field_jacobian)
        type(mse_channel_t), intent(in) :: channels(:)
        ! Cartesian B [T] at mse_geometry's midpoint for each channel.
        real(dp), intent(in) :: field(:, :)
        real(dp), allocatable, intent(out) :: angles(:)
        integer(i32), intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        real(dp), allocatable, intent(out), optional :: field_jacobian(:, :)
        real(dp) :: point(3), coefficients(6), rows(2, 3), scale, b(3), g(3)
        real(dp) :: num, den
        integer :: j

        ierr = 1_i32
        message = 'MSE field must have shape (3, number of channels)'
        if (size(field, 1) /= 3 .or. size(field, 2) /= size(channels)) return
        allocate (angles(size(channels)))
        if (present(field_jacobian)) allocate (field_jacobian(3, size(channels)))
        do j = 1, size(channels)
            call mse_geometry(channels(j), point, coefficients, ierr, message, &
                rows=rows)
            if (ierr /= 0) return
            ierr = 1_i32
            message = 'undefined MSE polarization at channel '//int_string(j)
            if (.not. all(ieee_is_finite(field(:, j)))) return
            scale = maxval(abs(field(:, j)))
            if (scale <= tiny(scale)) return
            b = field(:, j)/scale
            num = dot_product(rows(1, :), b)
            den = dot_product(rows(2, :), b)
            if (num*num + den*den <= 1.0e-24_dp) return
            call mse_local(b(1), b(2), b(3), rows(1, 1), rows(1, 2), rows(1, 3), &
                rows(2, 1), rows(2, 2), rows(2, 3), angles(j), g(1), g(2), g(3))
            ! Linear polarization is modulo pi; keep [-pi/2, pi/2).
            angles(j) = modulo(angles(j) + pi/2, pi) - pi/2
            if (present(field_jacobian)) field_jacobian(:, j) = g/scale
        end do
        ierr = 0_i32
        message = ''
    end subroutine mse_angles

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
            a(1)*b(2) - a(2)*b(1)]
    end function cross

    subroutine unit_vector(a, b, ok)
        real(dp), intent(in) :: a(3)
        real(dp), intent(out) :: b(3)
        logical, intent(out) :: ok
        real(dp) :: scale
        ok = .false.
        if (.not. all(ieee_is_finite(a))) return
        scale = maxval(abs(a))
        if (scale <= tiny(scale)) return
        b = a/scale
        b = b/norm2(b)
        ok = .true.
    end subroutine unit_vector

    subroutine data_line(unit, line, ios)
        integer, intent(in) :: unit
        character(len=*), intent(out) :: line
        integer, intent(out) :: ios
        integer :: k
        do
            read (unit, '(A)', iostat=ios) line
            if (ios /= 0) return
            k = scan(line, '#!')
            if (k > 0) line(k:) = ' '
            if (len_trim(line) > 0) return
        end do
    end subroutine data_line

    function int_string(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=32) :: buffer
        write (buffer, '(I0)') i
        s = trim(buffer)
    end function int_string
end module tiago_mse
