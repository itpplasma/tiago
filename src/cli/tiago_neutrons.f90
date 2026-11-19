program tiago_neutrons
    !! D-D neutron emission profiles and total rates from a fast-deuteron
    !! distribution and flux-surface profiles (tiago_neutron_rate).
    !!
    !! Usage: tiago_neutrons <distribution.txt> <profiles.txt> [Z_imp]
    !!
    !! profiles.txt: rows s, dV/ds [m^3], n_e [m^-3], T_e [eV], T_i [eV], Z_eff
    !! (lines starting with # skipped). Output on stdout: one row per profile row
    !! with the thermonuclear, beam-target and beam-beam emissivities
    !! [m^-3 s^-1], then the totals int eps dV/ds ds [1/s] (trapezoidal rule).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_profile_diagnostic, only: Q_S, Q_NE, Q_TE, Q_TI, Q_ZEFF, N_QUANTITIES
    use tiago_fast_ion_distribution, only: fast_ion_distribution_t, read_fast_ion_distribution
    use tiago_neutron_rate, only: neutron_t
    implicit none
    character(len=4096) :: dist_path, prof_path, arg
    type(fast_ion_distribution_t) :: dist
    type(neutron_t) :: model
    real(dp), allocatable :: rows(:, :), eps(:, :)
    real(dp) :: q(N_QUANTITIES), e3(3), de(3, N_QUANTITIES), z_imp, row(6), total(3)
    integer :: status, unit, ios, n, i
    character(len=1024) :: line

    if (command_argument_count() < 2) then
        write (error_unit, '(A)') 'usage: tiago_neutrons <distribution.txt> <profiles.txt> [Z_imp]'
        error stop 2
    end if
    call get_command_argument(1, dist_path)
    call get_command_argument(2, prof_path)
    z_imp = 7.0_dp
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read (arg, *) z_imp
    end if
    call read_fast_ion_distribution(trim(dist_path), dist, status)
    if (status /= 0) then
        write (error_unit, '(A,I0)') 'distribution unreadable or invalid, status ', status
        error stop 1
    end if

    allocate(rows(6, 0))
    open (newunit=unit, file=trim(prof_path), status='old', action='read')
    do
        read (unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line) == 0) cycle
        if (adjustl(line(1:1)) == '#') cycle
        read (line, *) row
        rows = reshape([rows, row], [6, size(rows, 2) + 1])
    end do
    close (unit)
    n = size(rows, 2)

    call model%setup(dist, z_imp)
    allocate(eps(3, n))
    write (*, '(A)') '# s eps_thermonuclear eps_beam_target eps_beam_beam [m^-3 s^-1]'
    do i = 1, n
        q = 0.0_dp
        q(Q_S) = rows(1, i)
        q(Q_NE) = rows(3, i)
        q(Q_TE) = rows(4, i)
        q(Q_TI) = rows(5, i)
        q(Q_ZEFF) = rows(6, i)
        call model%emissivity(q, e3, de)
        eps(:, i) = e3
        write (*, '(F10.6,3ES16.7)') rows(1, i), e3
    end do
    total = 0.0_dp
    do i = 1, n - 1
        total = total + 0.5_dp * (rows(1, i + 1) - rows(1, i)) &
            * (eps(:, i) * rows(2, i) + eps(:, i + 1) * rows(2, i + 1))
    end do
    write (*, '(A,ES16.7)') '# total_thermonuclear_per_s ', total(1)
    write (*, '(A,ES16.7)') '# total_beam_target_per_s ', total(2)
    write (*, '(A,ES16.7)') '# total_beam_beam_per_s ', total(3)
end program tiago_neutrons
