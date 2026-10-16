module tiago_diagnostic_types
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    implicit none
    private

    type :: loop_point_t
        real(dp) :: x
        real(dp) :: y
        real(dp) :: z
    end type loop_point_t

    type :: flux_loop_t
        character(len=:), allocatable :: label
        !> DIAGNO idia: 1 = diamagnetic loop (adds phiedge; zero in vacuum),
        !> -k = subtract the flux of loop k.
        integer(i32) :: idia = 0
        !> DIAGNO iflflg=1: the loop spans one field period; it is closed to its
        !> first point rotated by 2*pi/nfp and the flux is multiplied by nfp.
        logical :: one_period = .false.
        real(dp) :: turn_scale = 1.0_dp
        type(loop_point_t), allocatable :: points(:)
    end type flux_loop_t

    type :: segmented_rogowski_t
        character(len=:), allocatable :: label
        real(dp) :: turn_scale = 1.0_dp
        type(loop_point_t), allocatable :: path(:)
        !> DIAGNO eff_area of each path segment (segment j uses point j's value)
        real(dp), allocatable :: segment_area(:)
    end type segmented_rogowski_t


    type :: bprobe_t
        !> DIAGNO magnetic probe: signal = eff_area * B . normal * turn_scale
        character(len=:), allocatable :: label
        real(dp) :: position(3) = 0.0_dp
        real(dp) :: normal(3) = 0.0_dp
        real(dp) :: eff_area = 1.0_dp
        real(dp) :: turn_scale = 1.0_dp
    end type bprobe_t

    type :: mse_channel_t
        character(len=:), allocatable :: label
        real(dp) :: beam_point(3) = 0.0_dp
        real(dp) :: beam_direction(3) = 0.0_dp
        real(dp) :: sight_point(3) = 0.0_dp
        ! Direction from observation volume towards the observer.
        real(dp) :: sight_direction(3) = 0.0_dp
        ! Image-plane zero axis: projection of this reference vector.
        real(dp) :: reference(3) = [0.0_dp, 0.0_dp, 1.0_dp]
        real(dp) :: max_gap = 1.0e-3_dp
    end type mse_channel_t

    public :: mse_channel_t
    public :: loop_point_t
    public :: flux_loop_t
    public :: segmented_rogowski_t
    public :: bprobe_t
end module tiago_diagnostic_types
