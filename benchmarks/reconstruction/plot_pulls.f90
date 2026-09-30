program tiago_plot_pulls
    !! Pull (fitted - truth) / sigma per noise seed and parameter from the pull
    !! study of benchmarks/reconstruction/run.sh.
    !! Usage: tiago_plot_pulls pulls.csv pulls.png
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortplot, only: figure_t
    implicit none

    character(len=512) :: csv, png
    character(len=256) :: line
    character(len=64) :: names(32), name
    real(dp) :: seed(4096), pull(4096), fitted, sigma, truth, value, lo, hi
    integer :: owner(4096), n, nnames, unit, ios, s, k
    real(dp), parameter :: colors(3, 6) = reshape([0.12_dp, 0.47_dp, 0.71_dp, &
        0.84_dp, 0.15_dp, 0.16_dp, 0.17_dp, 0.63_dp, 0.17_dp, 1.0_dp, 0.5_dp, 0.05_dp, &
        0.58_dp, 0.40_dp, 0.74_dp, 0.55_dp, 0.34_dp, 0.29_dp], [3, 6])
    character(len=1), parameter :: markers(6) = ['o', 's', '^', 'D', 'v', 'x']
    type(figure_t) :: fig

    call get_command_argument(1, csv)
    call get_command_argument(2, png)
    open(newunit=unit, file=trim(csv), status='old', action='read')
    read(unit, '(A)') line
    n = 0
    nnames = 0
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        ! seed,"name",fitted,sigma,truth,pull (the name may contain commas)
        read(line(:index(line, ',') - 1), *) s
        name = line(index(line, '"') + 1:index(line, '"', back=.true.) - 1)
        read(line(index(line, '"', back=.true.) + 2:), *) fitted, sigma, truth, value
        n = n + 1
        seed(n) = s
        pull(n) = value
        owner(n) = findloc(names(:nnames), name, dim=1)
        if (owner(n) == 0) then
            nnames = nnames + 1
            names(nnames) = name
            owner(n) = nnames
        end if
    end do
    close(unit)

    lo = 0.5_dp
    hi = maxval(seed(:n)) + 0.5_dp
    call fig%initialize(900, 500)
    do k = -2, 2, 4
        call fig%add_plot([lo, hi], [real(k, dp), real(k, dp)], color=[0.8_dp, 0.8_dp, 0.8_dp], &
            linestyle='--')
    end do
    do k = -1, 1, 2
        call fig%add_plot([lo, hi], [real(k, dp), real(k, dp)], color=[0.55_dp, 0.55_dp, 0.55_dp], &
            linestyle='--')
    end do
    call fig%add_plot([lo, hi], [0.0_dp, 0.0_dp], color=[0.0_dp, 0.0_dp, 0.0_dp])
    do k = 1, nnames
        call fig%scatter(pack(seed(:n), owner(:n) == k) + 0.12_dp * (k - 0.5_dp * (nnames + 1)), &
            pack(pull(:n), owner(:n) == k), label=trim(names(k)), &
            color=colors(:, mod(k - 1, 6) + 1), marker=markers(mod(k - 1, 6) + 1))
    end do
    call fig%set_xlabel('noise seed')
    call fig%set_ylabel('(fitted - truth) / sigma')
    call fig%set_title('Pulls of the synthetic LI383 reconstruction (dashed: 1 and 2 sigma)')
    call fig%legend()
    call fig%savefig(trim(png))
end program tiago_plot_pulls
