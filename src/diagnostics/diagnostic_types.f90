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
        real(dp) :: effective_area = 0.0_dp
        integer(i32) :: segments = 0_i32
        real(dp) :: turn_scale = 1.0_dp
        type(loop_point_t), allocatable :: path(:)
    end type segmented_rogowski_t

    type :: metadata_entry_t
        character(len=:), allocatable :: name
        character(len=:), allocatable :: diag_type
        real(dp) :: effective_area = 0.0_dp
        integer(i32) :: segments = 0_i32
        real(dp) :: reference_orientation(3) = 0.0_dp
    end type metadata_entry_t

    public :: loop_point_t
    public :: flux_loop_t
    public :: segmented_rogowski_t
    public :: metadata_entry_t
end module tiago_diagnostic_types
