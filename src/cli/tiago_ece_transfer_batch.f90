program tiago_ece_transfer_batch
    !! Batch thermal ECE transfer on supplied cells, with the gradient of each
    !! brightness temperature with respect to cell te, ne, bfield and length.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_ece_opacity_transfer, only: thermal_ece_transfer, ece_opacity_transfer_t
    implicit none
    type :: ray_t
        real(dp), allocatable :: te(:), ne(:), bfield(:), length(:)
        real(dp), allocatable :: te_bar(:), ne_bar(:), b_bar(:), length_bar(:)
        real(dp) :: frequency, rho, wall_te, damping_limit, temperature, optical_depth
        integer :: mode, harmonic, status
    end type
    type(ray_t), allocatable :: rays(:)
    character(len=4096) :: input, output
    integer :: u, v, n_rays, k, i, n
    call get_command_argument(1, input)
    call get_command_argument(2, output)
    open (newunit=u, file=trim(input), status='old', action='read')
    read (u, *) n_rays
    if (n_rays < 1) error stop 'empty ECE transfer batch'
    allocate (rays(n_rays))
    do k = 1, n_rays
        associate (r => rays(k))
            read (u, *) n, r%frequency, r%mode, r%harmonic, r%rho, r%wall_te, r%damping_limit
            if (n < 1) error stop 'ECE transfer ray without cells'
            allocate (r%te(n), r%ne(n), r%bfield(n), r%length(n))
            do i = 1, n
                read (u, *) r%te(i), r%ne(i), r%bfield(i), r%length(i)
            end do
        end associate
    end do
    close (u)
    !$omp parallel do default(shared) private(k) schedule(dynamic)
    do k = 1, n_rays
        call solve(rays(k))
    end do
    !$omp end parallel do
    open (newunit=v, file=trim(output), status='replace', action='write')
    do k = 1, n_rays
        associate (r => rays(k))
            write (v, *) k, r%status, r%temperature, r%optical_depth
            if (r%status /= 0) cycle
            do i = 1, size(r%te)
                write (v, *) r%te_bar(i), r%ne_bar(i), r%b_bar(i), r%length_bar(i)
            end do
        end associate
    end do
    close (v)
contains
    subroutine solve(r)
        type(ray_t), intent(inout) :: r
        type(ece_opacity_transfer_t) :: trace
        real(dp) :: frequency_bar, rho_bar, wall_bar
        integer :: n
        n = size(r%te)
        r%temperature = 0.0_dp
        r%optical_depth = 0.0_dp
        call thermal_ece_transfer(r%te, r%ne, r%bfield, r%frequency, r%mode, r%harmonic, &
            r%length, r%rho, r%wall_te, trace, r%status, &
            damping_limit=r%damping_limit)
        if (r%status /= 0) return
        if (.not. trace%valid) then
            r%status = -1
            return
        end if
        r%temperature = trace%temperature
        r%optical_depth = trace%optical_depth
        allocate (r%te_bar(n), r%ne_bar(n), r%b_bar(n), r%length_bar(n))
        call trace%vjp(1.0_dp, r%te_bar, r%ne_bar, r%b_bar, frequency_bar, r%length_bar, &
            rho_bar, wall_bar)
    end subroutine
end program
