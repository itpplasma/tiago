module tiago_sensor_set
    !! One ingestion format for every magnetic sensor family.
    !!
    !! A sensor-set manifest (format "tiago-sensors 1", docs/sensor-sets.md)
    !! names geometry files in the existing DIAGNO formats and tags each group
    !! with a sensor kind. It may add linear combinations of raw signals, such as
    !! compensated diamagnetic loops. The geometry itself stays in the named
    !! files, so private machine geometry is read at run time from paths outside
    !! the repository (root may reference an environment variable).
    !!
    !! Raw signals are ordered loops, then segmented Rogowski coils, then probes.
    !! Every output signal is y = C raw with a sparse matrix C; plain sensors are
    !! unit rows. Since the map is linear, its JVP is C and its VJP is C^T.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    use tiago_diagnostic_types, only: flux_loop_t, segmented_rogowski_t, bprobe_t
    use tiago_flux_loops, only: read_flux_loop_file
    use tiago_segmented_rogowski, only: read_segmented_rogowski_file
    use tiago_bprobes, only: read_bprobe_file
    implicit none
    private

    integer, parameter, public :: kind_flux_loop = 1, kind_diamagnetic_loop = 2, &
        kind_compensation_loop = 3, kind_saddle_coil = 4, &
        kind_segmented_rogowski = 5, kind_bprobe = 6, kind_combination = 7
    character(len=*), parameter :: kind_names(7) = [character(len=20) :: &
        'flux_loop', 'diamagnetic_loop', 'compensation_loop', 'saddle_coil', &
        'segmented_rogowski', 'bprobe', 'combination']

    type, public :: label_t
        character(len=:), allocatable :: s
    end type label_t

    type, public :: sensor_set_t
        type(flux_loop_t), allocatable :: loops(:)
        integer, allocatable :: loop_kind(:)
        type(segmented_rogowski_t), allocatable :: segs(:)
        type(bprobe_t), allocatable :: probes(:)
        !> Output signals: label, kind, and CSR rows of C over raw signals.
        type(label_t), allocatable :: labels(:)
        integer, allocatable :: kinds(:)
        integer, allocatable :: row_start(:), col(:)
        real(dp), allocatable :: val(:)
    contains
        procedure :: n_raw => set_n_raw
        procedure :: n_signals => set_n_signals
        procedure :: combine => set_combine
        procedure :: combine_vjp => set_combine_vjp
        procedure :: combine_matrix => set_combine_matrix
        procedure :: index_of => set_index_of
    end type sensor_set_t

    public :: read_sensor_set
    public :: sensor_kind_name

contains

    function sensor_kind_name(kind) result(name)
        integer, intent(in) :: kind
        character(len=:), allocatable :: name
        name = trim(kind_names(kind))
    end function sensor_kind_name

    integer function set_n_raw(self)
        class(sensor_set_t), intent(in) :: self
        set_n_raw = 0
        if (allocated(self%loops)) set_n_raw = set_n_raw + size(self%loops)
        if (allocated(self%segs)) set_n_raw = set_n_raw + size(self%segs)
        if (allocated(self%probes)) set_n_raw = set_n_raw + size(self%probes)
    end function set_n_raw

    integer function set_n_signals(self)
        class(sensor_set_t), intent(in) :: self
        set_n_signals = 0
        if (allocated(self%labels)) set_n_signals = size(self%labels)
    end function set_n_signals

    subroutine set_combine(self, raw, y)
        !! y = C raw. Also the JVP for raw tangents.
        class(sensor_set_t), intent(in) :: self
        real(dp), intent(in) :: raw(:)
        real(dp), intent(out) :: y(:)
        integer :: i, k

        do i = 1, self%n_signals()
            y(i) = 0.0_dp
            do k = self%row_start(i), self%row_start(i + 1) - 1
                y(i) = y(i) + self%val(k) * raw(self%col(k))
            end do
        end do
    end subroutine set_combine

    subroutine set_combine_vjp(self, ybar, rawbar)
        !! rawbar = C^T ybar.
        class(sensor_set_t), intent(in) :: self
        real(dp), intent(in) :: ybar(:)
        real(dp), intent(out) :: rawbar(:)
        integer :: i, k

        rawbar = 0.0_dp
        do i = 1, self%n_signals()
            do k = self%row_start(i), self%row_start(i + 1) - 1
                rawbar(self%col(k)) = rawbar(self%col(k)) + self%val(k) * ybar(i)
            end do
        end do
    end subroutine set_combine_vjp

    subroutine set_combine_matrix(self, raw_jac, jac)
        !! jac = C raw_jac for a raw Jacobian with one row per raw signal.
        class(sensor_set_t), intent(in) :: self
        real(dp), intent(in) :: raw_jac(:, :)
        real(dp), allocatable, intent(out) :: jac(:, :)
        integer :: i, k

        allocate(jac(self%n_signals(), size(raw_jac, 2)))
        jac = 0.0_dp
        do i = 1, self%n_signals()
            do k = self%row_start(i), self%row_start(i + 1) - 1
                jac(i, :) = jac(i, :) + self%val(k) * raw_jac(self%col(k), :)
            end do
        end do
    end subroutine set_combine_matrix

    integer function set_index_of(self, label) result(idx)
        class(sensor_set_t), intent(in) :: self
        character(len=*), intent(in) :: label
        integer :: i

        idx = 0
        do i = 1, self%n_signals()
            if (self%labels(i)%s == label) then
                idx = i
                return
            end if
        end do
    end function set_index_of

    subroutine read_sensor_set(path, set, ierr, message)
        !! Manifest lines (blank lines and text after '#' are ignored):
        !!   tiago-sensors 1
        !!   root <dir>                       (default: manifest directory;
        !!                                     ${NAME} expands an environment variable)
        !!   group <kind> <file> [select=a,b,...] [turns=x] [area=x] [rphiz]
        !!   combine <label> <coef> <label> [<coef> <label> ...]
        character(len=*), intent(in) :: path
        type(sensor_set_t), intent(out) :: set
        integer(i32), intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message

        integer :: unit, ios, lineno, ntok, kind, i
        character(len=4096) :: line
        character(len=:), allocatable :: root, tokens(:)
        character(len=512), allocatable :: options(:)
        logical :: header_seen
        type(label_t), allocatable :: raw_labels(:)
        integer, allocatable :: raw_kinds(:)
        ! combinations, kept as label lists until all raw sensors are known
        type(label_t), allocatable :: combo_label(:), combo_terms(:)
        real(dp), allocatable :: combo_coef(:)
        integer, allocatable :: combo_owner(:)

        ierr = 0_i32
        message = ''
        allocate(set%loops(0), set%loop_kind(0), set%segs(0), set%probes(0))
        allocate(combo_label(0), combo_terms(0), combo_coef(0), combo_owner(0))
        root = dirname(path)
        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            call fail(1, 'unable to open sensor set: '//trim(path))
            return
        end if
        header_seen = .false.
        lineno = 0
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            lineno = lineno + 1
            i = index(line, '#')
            if (i > 0) line(i:) = ' '
            call split_tokens(line, tokens, ntok)
            if (ntok == 0) cycle
            if (.not. header_seen) then
                if (ntok /= 2 .or. tokens(1) /= 'tiago-sensors' .or. tokens(2) /= '1') then
                    call fail(2, trim(path)//': first line must be "tiago-sensors 1"')
                    exit
                end if
                header_seen = .true.
                cycle
            end if
            select case (tokens(1))
            case ('root')
                if (ntok /= 2) then
                    call fail(3, at(lineno)//'root expects one directory')
                    exit
                end if
                call expand_env(trim(tokens(2)), root, ierr, message)
                if (ierr /= 0) then
                    message = at(lineno)//message
                    exit
                end if
                if (.not. is_absolute(root)) root = dirname(path)//'/'//root
            case ('group')
                if (ntok < 3) then
                    call fail(3, at(lineno)//'group expects a kind and a file')
                    exit
                end if
                kind = kind_of(trim(tokens(2)))
                if (kind == 0 .or. kind == kind_combination) then
                    call fail(3, at(lineno)//'unknown sensor kind '//trim(tokens(2)))
                    exit
                end if
                ! Copy: gfortran mispasses sections of deferred-length arrays.
                allocate(options(ntok - 3))
                do i = 4, ntok
                    options(i - 3) = tokens(i)
                end do
                call add_group(kind, trim(tokens(3)), options)
                deallocate(options)
                if (ierr /= 0) exit
            case ('combine')
                if (ntok < 4 .or. mod(ntok, 2) /= 0) then
                    call fail(3, at(lineno)//'combine expects a label and coefficient-label pairs')
                    exit
                end if
                do i = 3, ntok, 2
                    combo_owner = [combo_owner, size(combo_label) + 1]
                    combo_terms = [combo_terms, label_t(trim(tokens(i + 1)))]
                    combo_coef = [combo_coef, real_value(tokens(i))]
                    if (ierr /= 0) exit
                end do
                if (ierr /= 0) exit
                combo_label = [combo_label, label_t(trim(tokens(2)))]
            case default
                call fail(3, at(lineno)//'unknown directive '//trim(tokens(1)))
                exit
            end select
        end do
        close(unit)
        if (ierr /= 0) return
        if (.not. header_seen) then
            call fail(2, trim(path)//': empty sensor set')
            return
        end if
        call build_outputs()

    contains

        subroutine add_group(kind, file, options)
            integer, intent(in) :: kind
            character(len=*), intent(in) :: file
            character(len=*), intent(in) :: options(:)
            type(flux_loop_t), allocatable :: loops(:)
            type(segmented_rogowski_t), allocatable :: segs(:)
            type(bprobe_t), allocatable :: probes(:)
            character(len=:), allocatable :: full, select_list
            logical, allocatable :: keep(:)
            real(dp) :: turns, area
            logical :: rphiz, has_area
            integer :: k, j

            full = file
            if (.not. is_absolute(full)) full = root//'/'//file
            select_list = ''
            turns = 1.0_dp
            area = 0.0_dp
            has_area = .false.
            rphiz = .false.
            do k = 1, size(options)
                j = index(options(k), '=')
                if (trim(options(k)) == 'rphiz') then
                    rphiz = .true.
                else if (j == 0) then
                    call fail(3, at(lineno)//'unknown group option '//trim(options(k)))
                else if (options(k)(:j - 1) == 'select') then
                    select_list = trim(options(k)(j + 1:))
                else if (options(k)(:j - 1) == 'turns') then
                    turns = real_value(options(k)(j + 1:))
                else if (options(k)(:j - 1) == 'area') then
                    area = real_value(options(k)(j + 1:))
                    has_area = .true.
                else
                    call fail(3, at(lineno)//'unknown group option '//trim(options(k)))
                end if
                if (ierr /= 0) return
            end do
            select case (kind)
            case (kind_segmented_rogowski)
                if (has_area) then
                    call read_segmented_rogowski_file(full, segs, ierr, message, area)
                else
                    call read_segmented_rogowski_file(full, segs, ierr, message)
                end if
                if (ierr /= 0) return
                call selection(select_list, [(label_t(segs(k)%label), k = 1, size(segs))], keep)
                if (ierr /= 0) return
                segs = pack(segs, keep)
                do k = 1, size(segs)
                    segs(k)%turn_scale = turns
                    if (all(segs(k)%segment_area == 0.0_dp)) then
                        call fail(4, at(lineno)//'segmented Rogowski '//segs(k)%label// &
                            ' has no segment areas; give area=')
                        return
                    end if
                end do
                set%segs = [set%segs, segs]
            case (kind_bprobe)
                call read_bprobe_file(full, probes, ierr, message, rphiz)
                if (ierr /= 0) return
                call selection(select_list, [(label_t(probes(k)%label), k = 1, size(probes))], keep)
                if (ierr /= 0) return
                probes = pack(probes, keep)
                probes%turn_scale = turns
                set%probes = [set%probes, probes]
            case default
                call read_flux_loop_file(full, loops, ierr, message)
                if (ierr /= 0) return
                call selection(select_list, [(label_t(loops(k)%label), k = 1, size(loops))], keep)
                if (ierr /= 0) return
                loops = pack(loops, keep)
                do k = 1, size(loops)
                    if (loops(k)%idia < 0) then
                        call fail(4, at(lineno)//'loop '//loops(k)%label// &
                            ' uses idia < 0; express differences with combine')
                        return
                    end if
                    ! The kind fixes whether the loop links the plasma poloidally.
                    if (kind == kind_diamagnetic_loop) then
                        loops(k)%idia = 1
                    else
                        if (loops(k)%idia == 1) then
                            call fail(4, at(lineno)//'loop '//loops(k)%label// &
                                ' has idia = 1; declare it as diamagnetic_loop')
                            return
                        end if
                    end if
                    loops(k)%turn_scale = turns
                end do
                set%loops = [set%loops, loops]
                set%loop_kind = [set%loop_kind, [(kind, k = 1, size(loops))]]
            end select
        end subroutine add_group

        subroutine selection(list, labels, keep)
            character(len=*), intent(in) :: list
            type(label_t), intent(in) :: labels(:)
            logical, allocatable, intent(out) :: keep(:)
            integer :: first, last, k
            logical :: found

            allocate(keep(size(labels)))
            keep = len_trim(list) == 0
            first = 1
            do while (first <= len_trim(list))
                last = index(list(first:), ',')
                last = merge(len_trim(list), first + last - 2, last == 0)
                found = .false.
                do k = 1, size(labels)
                    if (labels(k)%s == list(first:last)) then
                        keep(k) = .true.
                        found = .true.
                    end if
                end do
                if (.not. found) then
                    call fail(4, at(lineno)//'selected sensor '//list(first:last)//' not in file')
                    return
                end if
                first = last + 2
            end do
        end subroutine selection

        subroutine build_outputs()
            integer :: nraw, nout, k, j, c, m, n
            integer, allocatable :: count_terms(:)

            nraw = set%n_raw()
            allocate(raw_labels(nraw), raw_kinds(nraw))
            n = 0
            do k = 1, size(set%loops)
                n = n + 1
                raw_labels(n)%s = set%loops(k)%label
                raw_kinds(n) = set%loop_kind(k)
            end do
            do k = 1, size(set%segs)
                n = n + 1
                raw_labels(n)%s = set%segs(k)%label
                raw_kinds(n) = kind_segmented_rogowski
            end do
            do k = 1, size(set%probes)
                n = n + 1
                raw_labels(n)%s = set%probes(k)%label
                raw_kinds(n) = kind_bprobe
            end do
            do k = 1, nraw
                do j = 1, k - 1
                    if (raw_labels(j)%s == raw_labels(k)%s) then
                        call fail(5, trim(path)//': duplicate sensor label '//raw_labels(k)%s)
                        return
                    end if
                end do
            end do
            do k = 1, size(combo_label)
                if (find(raw_labels, combo_label(k)%s) /= 0) then
                    call fail(5, trim(path)//': combination label '//combo_label(k)%s// &
                        ' repeats a sensor label')
                    return
                end if
            end do
            nout = nraw + size(combo_label)
            allocate(count_terms(size(combo_label)))
            count_terms = [(count(combo_owner == k), k = 1, size(combo_label))]
            allocate(set%labels(nout), set%kinds(nout), set%row_start(nout + 1))
            allocate(set%col(nraw + size(combo_terms)), set%val(nraw + size(combo_terms)))
            set%labels(:nraw) = raw_labels
            set%kinds(:nraw) = raw_kinds
            set%row_start(1) = 1
            do k = 1, nraw
                set%col(k) = k
                set%val(k) = 1.0_dp
                set%row_start(k + 1) = k + 1
            end do
            m = nraw
            do k = 1, size(combo_label)
                set%labels(nraw + k) = combo_label(k)
                set%kinds(nraw + k) = kind_combination
                do j = 1, size(combo_terms)
                    if (combo_owner(j) /= k) cycle
                    c = find(raw_labels, combo_terms(j)%s)
                    if (c == 0) then
                        call fail(5, trim(path)//': combination '//combo_label(k)%s// &
                            ' refers to unknown sensor '//combo_terms(j)%s)
                        return
                    end if
                    m = m + 1
                    set%col(m) = c
                    set%val(m) = combo_coef(j)
                end do
                set%row_start(nraw + k + 1) = m + 1
            end do
        end subroutine build_outputs

        real(dp) function real_value(text)
            character(len=*), intent(in) :: text
            integer :: ios2
            read(text, *, iostat=ios2) real_value
            if (ios2 /= 0) then
                call fail(3, at(lineno)//'not a number: '//trim(text))
                real_value = 0.0_dp
            end if
        end function real_value

        subroutine fail(code, text)
            integer, intent(in) :: code
            character(len=*), intent(in) :: text
            ierr = int(code, i32)
            message = text
        end subroutine fail

        function at(n) result(s)
            integer, intent(in) :: n
            character(len=:), allocatable :: s
            character(len=16) :: buf
            write(buf, '(I0)') n
            s = trim(path)//':'//trim(buf)//': '
        end function at
    end subroutine read_sensor_set

    integer function kind_of(name)
        character(len=*), intent(in) :: name
        integer :: k
        kind_of = 0
        do k = 1, size(kind_names)
            if (trim(kind_names(k)) == name) kind_of = k
        end do
    end function kind_of

    integer function find(labels, label)
        type(label_t), intent(in) :: labels(:)
        character(len=*), intent(in) :: label
        integer :: k
        find = 0
        do k = 1, size(labels)
            if (labels(k)%s == label) then
                find = k
                return
            end if
        end do
    end function find

    subroutine split_tokens(line, tokens, n)
        character(len=*), intent(in) :: line
        character(len=:), allocatable, intent(out) :: tokens(:)
        integer, intent(out) :: n
        integer :: pos, first, last, width

        n = 0
        width = 1
        pos = 1
        do
            first = verify(line(pos:), ' '//achar(9))
            if (first == 0) exit
            first = pos + first - 1
            last = scan(line(first:), ' '//achar(9))
            last = merge(len(line), first + last - 2, last == 0)
            n = n + 1
            width = max(width, last - first + 1)
            pos = last + 1
            if (pos > len(line)) exit
        end do
        allocate(character(len=width) :: tokens(n))
        n = 0
        pos = 1
        do
            first = verify(line(pos:), ' '//achar(9))
            if (first == 0) exit
            first = pos + first - 1
            last = scan(line(first:), ' '//achar(9))
            last = merge(len(line), first + last - 2, last == 0)
            n = n + 1
            tokens(n) = line(first:last)
            pos = last + 1
            if (pos > len(line)) exit
        end do
    end subroutine split_tokens

    subroutine expand_env(text, out, ierr, message)
        !! Expands ${NAME}; an unset or empty variable is an error.
        character(len=*), intent(in) :: text
        character(len=:), allocatable, intent(out) :: out
        integer(i32), intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        integer :: pos, open_brace, close_brace, length, status
        character(len=4096) :: value

        ierr = 0_i32
        message = ''
        out = ''
        pos = 1
        do
            open_brace = index(text(pos:), '${')
            if (open_brace == 0) exit
            open_brace = pos + open_brace - 1
            close_brace = index(text(open_brace:), '}')
            if (close_brace == 0) then
                ierr = 3_i32
                message = 'unterminated ${ in '//text
                return
            end if
            close_brace = open_brace + close_brace - 1
            call get_environment_variable(text(open_brace + 2:close_brace - 1), value, &
                length, status)
            if (status /= 0 .or. length == 0) then
                ierr = 8_i32
                message = 'environment variable '//text(open_brace + 2:close_brace - 1)// &
                    ' is not set (it should point to the sensor geometry)'
                return
            end if
            out = out//text(pos:open_brace - 1)//value(:length)
            pos = close_brace + 1
        end do
        out = out//text(pos:)
    end subroutine expand_env

    pure logical function is_absolute(p)
        character(len=*), intent(in) :: p
        is_absolute = .false.
        if (len(p) > 0) is_absolute = p(1:1) == '/'
    end function is_absolute

    function dirname(p) result(d)
        character(len=*), intent(in) :: p
        character(len=:), allocatable :: d
        integer :: k
        k = index(trim(p), '/', back=.true.)
        if (k == 0) then
            d = '.'
        else if (k == 1) then
            d = '/'
        else
            d = p(:k - 1)
        end if
    end function dirname
end module tiago_sensor_set
