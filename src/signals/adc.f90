module tiago_adc
    !! Ideal digitiser helpers for a generic two's-complement ADC.
    !! A recorded count c maps to volts by V = c*vres + voff per channel;
    !! quantisation rounds to the nearest code and clips to the signed range.
    !! The uniform quantisation noise variance vres^2/12 follows Widrow and
    !! Kollar, Quantization Noise, Cambridge University Press 2008, ch. 4,
    !! and holds when the signal spans many codes (see docs/digitizer.md).
    use, intrinsic :: iso_fortran_env, only: dp => real64, int32
    implicit none
    private

    integer(int32), parameter, public :: ADC_MIN_16BIT = -32768_int32
    integer(int32), parameter, public :: ADC_MAX_16BIT = 32767_int32

    public :: adc_decode, adc_quantize, adc_quantization_variance

    interface adc_decode
        module procedure adc_decode_scalar
        module procedure adc_decode_array
    end interface

    interface adc_quantize
        module procedure adc_quantize_scalar
        module procedure adc_quantize_array
    end interface
contains
    pure elemental function adc_decode_scalar(counts, vres, voff) result(volts)
        integer(int32), intent(in) :: counts
        real(dp), intent(in) :: vres, voff
        real(dp) :: volts
        volts = real(counts, dp)*vres + voff
    end function

    pure function adc_decode_array(counts, vres, voff) result(volts)
        integer(int32), intent(in) :: counts(:)
        real(dp), intent(in) :: vres(:), voff(:)
        real(dp) :: volts(size(counts))
        volts = real(counts, dp)*vres + voff
    end function

    pure elemental subroutine adc_quantize_scalar(volts, vres, voff, bits, &
            counts, saturated)
        real(dp), intent(in) :: volts, vres, voff
        integer, intent(in) :: bits
        integer(int32), intent(out) :: counts
        logical, intent(out), optional :: saturated
        integer(int32) :: lo, hi
        lo = shiftl(-1_int32, bits - 1)
        hi = shiftl(1_int32, bits - 1) - 1_int32
        counts = nint((volts - voff)/vres)
        if (present(saturated)) saturated = counts < lo .or. counts > hi
        counts = min(max(counts, lo), hi)
    end subroutine

    pure subroutine adc_quantize_array(volts, vres, voff, bits, counts, saturated)
        real(dp), intent(in) :: volts(:), vres(:), voff(:)
        integer, intent(in) :: bits
        integer(int32), intent(out) :: counts(:)
        logical, intent(out), optional :: saturated(:)
        integer(int32) :: lo, hi
        integer :: i
        lo = shiftl(-1_int32, bits - 1)
        hi = shiftl(1_int32, bits - 1) - 1_int32
        do i = 1, size(volts)
            counts(i) = nint((volts(i) - voff(i))/vres(i))
            if (present(saturated)) saturated(i) = counts(i) < lo .or. counts(i) > hi
            counts(i) = min(max(counts(i), lo), hi)
        end do
    end subroutine

    pure elemental function adc_quantization_variance(vres) result(var)
        real(dp), intent(in) :: vres
        real(dp) :: var
        var = vres*vres/12.0_dp
    end function
end module tiago_adc
