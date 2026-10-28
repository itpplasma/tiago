module tiago_sha256
    !! SHA-256 (FIPS 180-4) of files, for input provenance in outputs.
    !! Words are held in 64-bit integers masked to 32 bits.
    use, intrinsic :: iso_fortran_env, only: i8 => int64, int8
    implicit none
    private

    integer(i8), parameter :: mask = int(z'FFFFFFFF', i8)
    integer(i8), parameter :: k(64) = [ &
        int(z'428a2f98', i8), int(z'71374491', i8), int(z'b5c0fbcf', i8), int(z'e9b5dba5', i8), &
        int(z'3956c25b', i8), int(z'59f111f1', i8), int(z'923f82a4', i8), int(z'ab1c5ed5', i8), &
        int(z'd807aa98', i8), int(z'12835b01', i8), int(z'243185be', i8), int(z'550c7dc3', i8), &
        int(z'72be5d74', i8), int(z'80deb1fe', i8), int(z'9bdc06a7', i8), int(z'c19bf174', i8), &
        int(z'e49b69c1', i8), int(z'efbe4786', i8), int(z'0fc19dc6', i8), int(z'240ca1cc', i8), &
        int(z'2de92c6f', i8), int(z'4a7484aa', i8), int(z'5cb0a9dc', i8), int(z'76f988da', i8), &
        int(z'983e5152', i8), int(z'a831c66d', i8), int(z'b00327c8', i8), int(z'bf597fc7', i8), &
        int(z'c6e00bf3', i8), int(z'd5a79147', i8), int(z'06ca6351', i8), int(z'14292967', i8), &
        int(z'27b70a85', i8), int(z'2e1b2138', i8), int(z'4d2c6dfc', i8), int(z'53380d13', i8), &
        int(z'650a7354', i8), int(z'766a0abb', i8), int(z'81c2c92e', i8), int(z'92722c85', i8), &
        int(z'a2bfe8a1', i8), int(z'a81a664b', i8), int(z'c24b8b70', i8), int(z'c76c51a3', i8), &
        int(z'd192e819', i8), int(z'd6990624', i8), int(z'f40e3585', i8), int(z'106aa070', i8), &
        int(z'19a4c116', i8), int(z'1e376c08', i8), int(z'2748774c', i8), int(z'34b0bcb5', i8), &
        int(z'391c0cb3', i8), int(z'4ed8aa4a', i8), int(z'5b9cca4f', i8), int(z'682e6ff3', i8), &
        int(z'748f82ee', i8), int(z'78a5636f', i8), int(z'84c87814', i8), int(z'8cc70208', i8), &
        int(z'90befffa', i8), int(z'a4506ceb', i8), int(z'bef9a3f7', i8), int(z'c67178f2', i8)]

    public :: sha256_file, sha256_bytes

contains

    function sha256_file(path) result(hex)
        !! Lower-case hex digest, or '' if the file cannot be read.
        character(len=*), intent(in) :: path
        character(len=:), allocatable :: hex
        integer(int8), allocatable :: data(:)
        integer :: unit, ios, n

        hex = ''
        open(newunit=unit, file=path, access='stream', form='unformatted', &
            status='old', action='read', iostat=ios)
        if (ios /= 0) return
        inquire(unit=unit, size=n)
        allocate(data(max(n, 0)))
        if (n > 0) read(unit, iostat=ios) data
        close(unit)
        if (ios /= 0) return
        hex = sha256_bytes(data)
    end function sha256_file

    function sha256_bytes(data) result(hex)
        integer(int8), intent(in) :: data(:)
        character(len=64) :: hex
        integer(i8) :: h(8), w(64), a(8), t1, t2, bitlen
        integer(int8), allocatable :: msg(:)
        integer :: n, padded, blk, j, i

        h = [int(z'6a09e667', i8), int(z'bb67ae85', i8), int(z'3c6ef372', i8), &
            int(z'a54ff53a', i8), int(z'510e527f', i8), int(z'9b05688c', i8), &
            int(z'1f83d9ab', i8), int(z'5be0cd19', i8)]
        n = size(data)
        padded = ((n + 8) / 64 + 1) * 64
        allocate(msg(padded))
        msg = 0_int8
        msg(:n) = data
        msg(n + 1) = int(z'80', int8)
        bitlen = 8_i8 * int(n, i8)
        do j = 0, 7
            msg(padded - j) = int(iand(shiftr(bitlen, 8 * j), 255_i8), int8)
        end do
        do blk = 0, padded / 64 - 1
            do j = 1, 16
                w(j) = 0_i8
                do i = 0, 3
                    w(j) = ior(shiftl(w(j), 8), &
                        iand(int(msg(64 * blk + 4 * (j - 1) + i + 1), i8), 255_i8))
                end do
            end do
            do j = 17, 64
                w(j) = iand(w(j - 16) + s0(w(j - 15)) + w(j - 7) + s1(w(j - 2)), mask)
            end do
            a = h
            do j = 1, 64
                t1 = iand(a(8) + big1(a(5)) + ior(iand(a(5), a(6)), &
                    iand(ieor(a(5), mask), a(7))) + k(j) + w(j), mask)
                t2 = iand(big0(a(1)) + ieor(ieor(iand(a(1), a(2)), iand(a(1), a(3))), &
                    iand(a(2), a(3))), mask)
                a(8) = a(7)
                a(7) = a(6)
                a(6) = a(5)
                a(5) = iand(a(4) + t1, mask)
                a(4) = a(3)
                a(3) = a(2)
                a(2) = a(1)
                a(1) = iand(t1 + t2, mask)
            end do
            h = iand(h + a, mask)
        end do
        do j = 1, 8
            write(hex(8 * j - 7:8 * j), '(Z8.8)') h(j)
        end do
        do j = 1, 64
            if (hex(j:j) >= 'A' .and. hex(j:j) <= 'F') &
                hex(j:j) = achar(iachar(hex(j:j)) + 32)
        end do
    end function sha256_bytes

    pure integer(i8) function rotr(x, n)
        integer(i8), intent(in) :: x
        integer, intent(in) :: n
        rotr = iand(ior(shiftr(x, n), shiftl(x, 32 - n)), mask)
    end function rotr

    pure integer(i8) function s0(x)
        integer(i8), intent(in) :: x
        s0 = ieor(ieor(rotr(x, 7), rotr(x, 18)), shiftr(x, 3))
    end function s0

    pure integer(i8) function s1(x)
        integer(i8), intent(in) :: x
        s1 = ieor(ieor(rotr(x, 17), rotr(x, 19)), shiftr(x, 10))
    end function s1

    pure integer(i8) function big0(x)
        integer(i8), intent(in) :: x
        big0 = ieor(ieor(rotr(x, 2), rotr(x, 13)), rotr(x, 22))
    end function big0

    pure integer(i8) function big1(x)
        integer(i8), intent(in) :: x
        big1 = ieor(ieor(rotr(x, 6), rotr(x, 11)), rotr(x, 25))
    end function big1
end module tiago_sha256
