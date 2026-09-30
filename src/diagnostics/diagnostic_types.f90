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

    public :: loop_point_t
    public :: flux_loop_t
    public :: segmented_rogowski_t
    public :: bprobe_t
end module tiago_diagnostic_types
