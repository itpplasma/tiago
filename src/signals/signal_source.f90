module tiago_signal_source
    use, intrinsic :: iso_fortran_env, only: dp => real64, int8, int32, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortio_checksum, only: lookup3_checksum
    implicit none
    private

    integer, parameter, public :: SIGNAL_OK = 0, SIGNAL_ERROR = 1, SIGNAL_MORE = 2

    type, public :: signal_request_t
        character(:), allocatable :: channel, version
        real(dp) :: first = 0.0_dp, last = 0.0_dp
    contains
        procedure :: key => request_key
        procedure :: valid => request_valid
    end type

    type, public :: signal_provenance_t
        character(:), allocatable :: source_id, version, request_key, checksum
        character(:), allocatable :: record_id, metadata
    end type

    type, public :: signal_window_t
        real(dp), allocatable :: time(:), samples(:), variance(:)
        real(dp), allocatable :: covariance(:, :) ! optional correlated contribution
        character(:), allocatable :: units
        type(signal_provenance_t) :: provenance
        real(dp) :: valid_from = -huge(1.0_dp), valid_to = huge(1.0_dp)
    contains
        procedure :: valid => window_valid
        procedure :: has_provenance => window_has_provenance
        procedure :: seal => window_seal
        procedure :: total_covariance => window_covariance
    end type

    type, abstract, public :: signal_source_t
    contains
        procedure(read_i), deferred :: read
    end type

    type, public :: backend_request_t
        character(:), allocatable :: endpoint, method, body, request_key
    end type

    type, abstract, public :: signal_adapter_t
    contains
        procedure(configure_i), deferred :: configure
        procedure(map_i), deferred :: map_request
        procedure(decode_i), deferred :: decode
    end type

    abstract interface
        subroutine read_i(self, request, window, status)
            import :: signal_source_t, signal_request_t, signal_window_t
            class(signal_source_t), intent(inout) :: self
            type(signal_request_t), intent(in) :: request
            type(signal_window_t), intent(out) :: window
            integer, intent(out) :: status
        end subroutine
        subroutine configure_i(self, path, status)
            import :: signal_adapter_t
            class(signal_adapter_t), intent(inout) :: self
            character(*), intent(in) :: path
            integer, intent(out) :: status
        end subroutine
        subroutine map_i(self, request, backend, status)
            import :: signal_adapter_t, signal_request_t, backend_request_t
            class(signal_adapter_t), intent(in) :: self
            type(signal_request_t), intent(in) :: request
            type(backend_request_t), intent(out) :: backend
            integer, intent(out) :: status
        end subroutine
        subroutine decode_i(self, request, bytes, final, window, status)
            import :: signal_adapter_t, signal_request_t, signal_window_t, int8
            class(signal_adapter_t), intent(inout) :: self
            type(signal_request_t), intent(in) :: request
            integer(int8), intent(in) :: bytes(:)
            logical, intent(in) :: final
            type(signal_window_t), intent(out) :: window
            integer, intent(out) :: status
        end subroutine
    end interface

    type :: cache_entry_t
        character(:), allocatable :: key
        type(signal_window_t) :: window
    end type

    type, extends(signal_source_t), public :: cached_source_t
        private
        class(signal_source_t), pointer :: upstream => null()
        type(cache_entry_t), allocatable :: entries(:)
        integer :: next = 1
    contains
        procedure :: init => cache_init
        procedure :: clear => cache_clear
        procedure :: read => cache_read
    end type

    public :: content_checksum, version_matches
contains
    logical function request_valid(self)
        class(signal_request_t), intent(in) :: self
        request_valid = .false.
        if (.not. allocated(self%channel)) return
        if (len(self%channel) == 0) return
        if (.not. ieee_is_finite(self%first)) return
        if (.not. ieee_is_finite(self%last)) return
        request_valid = self%first <= self%last
    end function

    function request_key(self) result(key)
        class(signal_request_t), intent(in) :: self
        character(:), allocatable :: key, version
        character(80) :: numbers
        character(32) :: length_text
        version = ''
        if (allocated(self%version)) version = self%version
        write (numbers, '(Z16.16,":",Z16.16)') &
            transfer(self%first, 0_int64), transfer(self%last, 0_int64)
        write (length_text, '(I0)') len(self%channel)
        key = trim(length_text)//':'//self%channel//':'//trim(numbers)//':'//version
    end function

    logical function version_matches(request, version)
        type(signal_request_t), intent(in) :: request
        character(*), intent(in) :: version
        version_matches = .true.
        if (allocated(request%version)) then
            if (len(request%version) > 0) version_matches = request%version == version
        end if
    end function

    logical function window_valid(self)
        class(signal_window_t), intent(in) :: self
        integer :: n
        window_valid = .false.
        if (.not. allocated(self%time)) return
        if (.not. allocated(self%samples)) return
        if (.not. allocated(self%variance)) return
        if (.not. allocated(self%units)) return
        n = size(self%time)
        if (size(self%samples) /= n .or. size(self%variance) /= n) return
        if (.not. all(ieee_is_finite(self%time))) return
        if (.not. all(ieee_is_finite(self%samples))) return
        if (.not. all(ieee_is_finite(self%variance))) return
        if (any(self%variance < 0.0_dp)) return
        if (n > 1) then
            if (any(self%time(2:) < self%time(:n - 1))) return
        end if
        if (allocated(self%covariance)) then
            if (any(shape(self%covariance) /= [n, n])) return
            if (.not. all(ieee_is_finite(self%covariance))) return
            if (any(self%covariance /= transpose(self%covariance))) return
        end if
        if (.not. ieee_is_finite(self%valid_from)) return
        if (.not. ieee_is_finite(self%valid_to)) return
        if (self%valid_to < self%valid_from) return
        window_valid = .true.
    end function

    logical function window_has_provenance(self)
        class(signal_window_t), intent(in) :: self
        window_has_provenance = .false.
        if (.not. allocated(self%provenance%source_id)) return
        if (.not. allocated(self%provenance%version)) return
        if (.not. allocated(self%provenance%request_key)) return
        if (.not. allocated(self%provenance%checksum)) return
        if (.not. allocated(self%provenance%metadata)) return
        if (.not. allocated(self%provenance%record_id)) return
        if (len(self%provenance%source_id) == 0) return
        if (len(self%provenance%version) == 0) return
        if (len(self%provenance%request_key) == 0) return
        if (len(self%provenance%checksum) == 0) return
        window_has_provenance = .true.
    end function

    function content_checksum(window) result(checksum)
        type(signal_window_t), intent(in) :: window
        character(:), allocatable :: checksum
        character(16) :: hex
        character(8) :: digest
        integer(int32) :: seed
        integer :: i, j
        seed = 0_int32
        call add_text(window%units)
        call add_real(window%valid_from)
        call add_real(window%valid_to)
        do i = 1, size(window%samples)
            call add_real(window%time(i))
            call add_real(window%samples(i))
            call add_real(window%variance(i))
        end do
        if (allocated(window%covariance)) then
            call add_text('covariance')
            do j = 1, size(window%samples)
                do i = 1, size(window%samples)
                    call add_real(window%covariance(i, j))
                end do
            end do
        end if
        write (digest, '(Z8.8)') seed
        checksum = 'lookup3-hex64:'//digest
    contains
        subroutine add_real(value)
            real(dp), intent(in) :: value
            write (hex, '(Z16.16)') transfer(value, 0_int64)
            call add_text(hex)
        end subroutine
        subroutine add_text(value)
            character(*), intent(in) :: value
            seed = lookup3_checksum(transfer(value, [0_int8], len(value)), seed)
        end subroutine
    end function

    function window_covariance(self) result(c)
        class(signal_window_t), intent(in) :: self
        real(dp), allocatable :: c(:, :)
        integer :: n, i
        n = size(self%samples)
        allocate (c(n, n))
        c = 0.0_dp
        if (allocated(self%covariance)) c = self%covariance
        do i = 1, n
            c(i, i) = c(i, i) + self%variance(i)
        end do
    end function

    subroutine window_seal(self)
        class(signal_window_t), intent(inout) :: self
        self%provenance%checksum = content_checksum(self)
    end subroutine

    subroutine cache_init(self, upstream, capacity)
        class(cached_source_t), intent(inout) :: self
        class(signal_source_t), target, intent(inout) :: upstream
        integer, intent(in) :: capacity
        if (capacity < 1) error stop 'signal cache: positive capacity required'
        self%upstream => upstream
        if (allocated(self%entries)) deallocate (self%entries)
        allocate (self%entries(capacity))
        self%next = 1
    end subroutine

    subroutine cache_clear(self)
        class(cached_source_t), intent(inout) :: self
        integer :: i
        if (.not. allocated(self%entries)) return
        do i = 1, size(self%entries)
            self%entries(i) = cache_entry_t()
        end do
        self%next = 1
    end subroutine

    subroutine cache_read(self, request, window, status)
        class(cached_source_t), intent(inout) :: self
        type(signal_request_t), intent(in) :: request
        type(signal_window_t), intent(out) :: window
        integer, intent(out) :: status
        character(:), allocatable :: key
        integer :: i
        status = 1
        if (.not. request%valid()) return
        if (.not. associated(self%upstream)) return
        key = request%key()
        do i = 1, size(self%entries)
            if (.not. allocated(self%entries(i)%key)) cycle
            if (self%entries(i)%key /= key) cycle
            window = self%entries(i)%window
            status = 0
            return
        end do
        call self%upstream%read(request, window, status)
        if (status /= 0) return
        status = 1
        if (.not. window%valid()) return
        if (.not. window%has_provenance()) return
        if (.not. allocated(window%provenance%version)) return
        if (.not. version_matches(request, window%provenance%version)) return
        if (.not. allocated(window%provenance%checksum)) return
        if (content_checksum(window) /= window%provenance%checksum) return
        self%entries(self%next)%key = key
        self%entries(self%next)%window = window
        self%next = mod(self%next, size(self%entries)) + 1
        status = 0
    end subroutine
end module tiago_signal_source
