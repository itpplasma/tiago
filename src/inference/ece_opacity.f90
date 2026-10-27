module tiago_ece_opacity
    !! Perpendicular, weakly relativistic Maxwellian harmonic opacity.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_special_dnestrovskii, only: dnestrovskii
    use tiago_ece_opacity_coordinates, only: ece_opacity_coordinates
    use tiago_ece_opacity_coordinates2, only: ece_opacity_coordinates2
    use tiago_ece_opacity_o1, only: ece_opacity_o1
    use tiago_ece_opacity_o2, only: ece_opacity_o2
    use tiago_ece_opacity_o2_blue, only: ece_opacity_o2_blue
    use tiago_ece_opacity_x1_tenuous, only: ece_opacity_x1_tenuous
    use tiago_ece_opacity_x1_finite, only: ece_opacity_x1_finite
    use tiago_ece_opacity_x2, only: ece_opacity_x2
    use tiago_ece_opacity_emission, only: ece_opacity_emission
    use tiago_ece_ray, only: ECE_O_MODE, ECE_X_MODE
    implicit none
    private
    integer, parameter, public :: ECE_X1_TENUOUS = 1, ECE_X1_FINITE = 2
    integer, parameter, public :: ECE_OPACITY_OK = 0, ECE_OPACITY_INPUT = 1, &
                                  ECE_OPACITY_DOMAIN = 2
    public :: ece_opacity, ece_opacity_batch, ECE_O_MODE, ECE_X_MODE
contains
    pure subroutine ece_opacity(te, ne, bfield, frequency, mode, harmonic, &
                                alpha, emissivity, index2, jacobian, status, x1_regime)
        real(dp), intent(in) :: te, ne, bfield, frequency
        integer, intent(in) :: mode, harmonic
        real(dp), intent(out) :: alpha, emissivity, index2, jacobian(3, 4)
        integer, intent(out) :: status
        integer, intent(in), optional :: x1_regime
        complex(dp) :: f5, f7, d5, d7
        real(dp) :: coord(4), cr(4, 4), value(2), row(2, 8), &
                    chained(2, 4), emission(1), erow(1, 4), frow(4), ratio
        integer :: regime, k, code
        alpha = 0.0_dp
        emissivity = 0.0_dp
        index2 = 0.0_dp
        jacobian = 0.0_dp
        status = ECE_OPACITY_INPUT
        if (.not. all(ieee_is_finite([te, ne, bfield, frequency]))) return
        if (te <= 0.0_dp .or. ne < 0.0_dp .or. bfield <= 0.0_dp .or. frequency <= 0.0_dp) return
        if (mode /= ECE_O_MODE .and. mode /= ECE_X_MODE) return
        if (harmonic /= 1 .and. harmonic /= 2) return
        regime = ECE_X1_FINITE
        if (present(x1_regime)) regime = x1_regime
        if (regime /= ECE_X1_TENUOUS .and. regime /= ECE_X1_FINITE) return
        if (harmonic == 1) then
            call ece_opacity_coordinates(te, ne, bfield, frequency, coord, cr)
        else
            call ece_opacity_coordinates2(te, ne, bfield, frequency, coord, cr)
        end if
        status = ECE_OPACITY_DOMAIN
        if (.not. all(ieee_is_finite(coord))) return
        if (.not. all(ieee_is_finite(cr))) return
        if (coord(1) > 0.01_dp .or. coord(1) < 1.0e-12_dp) return
        if (abs(coord(1)*coord(4)) > 0.1_dp) return
        ! The asymptotic models exclude mode cutoffs and strong damping.
        if (mode == ECE_O_MODE) then
            if (coord(3) >= real(harmonic**2, dp)) return
        else if (harmonic == 2) then
            if (coord(3) >= 2.0_dp) return
        else if (regime == ECE_X1_FINITE) then
            if (coord(3) >= 2.0_dp .or. coord(3) < 20*coord(1)) return
        else
            if (coord(3) > 0.2_dp*coord(1)) return
        end if
        if (abs(coord(4)) > 1.0e80_dp) return
        call dnestrovskii(5, coord(4), f5, d5, code)
        if (code /= 0) return
        call dnestrovskii(7, coord(4), f7, d7, code)
        if (code /= 0) return
        if (mode == ECE_O_MODE) then
            if (harmonic == 1) then
                if (abs(coord(3)*aimag(f7)) > 0.1_dp) return
                call ece_opacity_o1(te, ne, bfield, frequency, real(f5, dp), aimag(f5), &
                                    real(f7, dp), aimag(f7), value, row)
            else if (coord(4) <= 0.0_dp) then
                call ece_opacity_o2(te, ne, bfield, frequency, real(f5, dp), aimag(f5), &
                                    real(f7, dp), aimag(f7), value, row)
            else
                call ece_opacity_o2_blue(te, ne, bfield, frequency, real(f5, dp), aimag(f5), &
                                         real(f7, dp), aimag(f7), value, row)
            end if
        else if (harmonic == 2) then
            if (abs(coord(3)*aimag(f7)) > 0.1_dp) return
            call ece_opacity_x2(te, ne, bfield, frequency, real(f5, dp), aimag(f5), &
                                real(f7, dp), aimag(f7), value, row)
        else if (regime == ECE_X1_TENUOUS) then
            call ece_opacity_x1_tenuous(te, ne, bfield, frequency, real(f5, dp), aimag(f5), &
                                        real(f7, dp), aimag(f7), value, row)
        else
            call ece_opacity_x1_finite(te, ne, bfield, frequency, real(f5, dp), aimag(f5), &
                                       real(f7, dp), aimag(f7), value, row)
        end if
        if (.not. all(ieee_is_finite(value))) return
        if (.not. all(ieee_is_finite(row))) return
        if (value(1) < 0.0_dp .or. value(2) <= 0.0_dp) return
        ratio = value(1)*299792458.0_dp/ &
                (4*acos(-1.0_dp)*frequency*sqrt(value(2)))
        if (ratio > 0.1_dp) return
        frow = [real(d5, dp), aimag(d5), real(d7, dp), aimag(d7)]
        chained = row(:, :4)
        do k = 1, 4
            chained(:, k) = chained(:, k) + matmul(row(:, 5:8), frow)*cr(4, k)
        end do
        call ece_opacity_emission(value(1), value(2), te, frequency, emission, erow)
        alpha = value(1)
        index2 = value(2)
        emissivity = emission(1)
        jacobian(1, :) = chained(1, :)
        jacobian(3, :) = chained(2, :)
        jacobian(2, :) = erow(1, 1)*chained(1, :) + erow(1, 2)*chained(2, :)
        jacobian(2, 1) = jacobian(2, 1) + erow(1, 3)
        jacobian(2, 4) = jacobian(2, 4) + erow(1, 4)
        if (.not. ieee_is_finite(emissivity)) return
        if (.not. all(ieee_is_finite(jacobian))) return
        status = ECE_OPACITY_OK
    end subroutine

    subroutine ece_opacity_batch(te, ne, bfield, frequency, mode, harmonic, alpha, &
                                 emissivity, index2, jacobian, status, x1_regime)
        real(dp), intent(in) :: te(:), ne(:), bfield(:), frequency(:)
        integer, intent(in) :: mode(:), harmonic(:)
        real(dp), intent(out) :: alpha(:), emissivity(:), index2(:), jacobian(:, :, :)
        integer, intent(out) :: status(:)
        integer, intent(in), optional :: x1_regime
        integer :: j, n
        n = size(te)
        status = ECE_OPACITY_INPUT
        if (size(ne) /= n .or. size(bfield) /= n .or. size(frequency) /= n) return
        if (size(mode) /= n .or. size(harmonic) /= n .or. size(status) /= n) return
        if (size(alpha) /= n .or. size(emissivity) /= n .or. size(index2) /= n) return
        if (any(shape(jacobian) /= [3, 4, n])) return
        !$omp parallel do schedule(static)
        do j = 1, n
            call ece_opacity(te(j), ne(j), bfield(j), frequency(j), mode(j), harmonic(j), &
                        alpha(j), emissivity(j), index2(j), jacobian(:, :, j), status(j), x1_regime)
        end do
        !$omp end parallel do
    end subroutine
end module
