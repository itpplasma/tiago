program tiago_bootstrap
    !! Sauter bootstrap current on the flux surfaces of a G-EQDSK equilibrium
    !! (docs/bootstrap.md).
    !! Usage: tiago_bootstrap geqdsk kinetic.csv relative_sigma psi_N...
    !! kinetic.csv: header, then psi_N,n_e,T_e,T_i,n_i,Z_eff on a uniform psi_N
    !! grid from 0 to 1 [m^-3, eV], main ions of charge 1. Prints one row
    !! jdotb@<psi_N>,<j.B>_bs,sigma per label for tiago_axisym_fit, with sigma the
    !! relative_sigma times |<j.B>_bs|, and writes the surface quantities to
    !! stderr as psi_N,f_trap,epsilon,R_geo,q,nu_e*,nu_i*,<j.B>_bs,<j.B>_eq,F,<B^2>.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use batch_interpolate, only: BatchSplineData1D, construct_batch_splines_1d, &
        evaluate_batch_splines_1d_der
    use tiago_geqdsk_plasma, only: geqdsk_plasma_t, read_geqdsk_plasma
    use tiago_geqdsk_surfaces, only: flux_surface_t, flux_surface, surface_jdotb
    use tiago_sauter_bootstrap, only: sauter_jdotb, sauter_nu_e, sauter_nu_i
    implicit none

    type(geqdsk_plasma_t) :: eq
    type(flux_surface_t) :: surf
    type(BatchSplineData1D) :: kin
    character(len=1024) :: arg
    character(len=:), allocatable :: message
    character(len=32) :: buf
    real(dp), allocatable :: table(:, :)
    real(dp) :: rel, x, y(5), dy(5), dpsi, jb, row(6)
    integer :: ierr, unit, ios, n, i, k

    if (command_argument_count() < 4) &
        error stop 'usage: tiago_bootstrap geqdsk kinetic.csv relative_sigma psi_N...'
    call get_command_argument(1, arg)
    call read_geqdsk_plasma(trim(arg), eq, ierr, message)
    if (ierr /= 0) error stop message
    call get_command_argument(2, arg)
    open(newunit=unit, file=trim(arg), status='old', action='read')
    read(unit, *)
    n = 0
    do
        read(unit, *, iostat=ios) row
        if (ios /= 0) exit
        n = n + 1
    end do
    rewind(unit)
    read(unit, *)
    allocate(table(n, 6))
    do i = 1, n
        read(unit, *) table(i, :)
    end do
    close(unit)
    if (n < 4 .or. abs(table(1, 1)) > 1.0e-12_dp .or. abs(table(n, 1) - 1) > 1.0e-12_dp .or. &
        maxval(abs(table(2:, 1) - table(:n - 1, 1) - 1.0_dp / (n - 1))) > 1.0e-9_dp) &
        error stop 'kinetic profiles need a uniform psi_N grid from 0 to 1'
    call construct_batch_splines_1d(0.0_dp, 1.0_dp, table(:, 2:6), 3, .false., kin)
    call get_command_argument(3, arg)
    read(arg, *) rel
    dpsi = eq%psi_b - eq%psi_axis
    print '(A)', 'label,value,sigma'
    write(error_unit, '(A)') 'psin,f_trap,epsilon,r_geo,q,nu_e,nu_i,jdotb_bs,jdotb_eq,f,b2'
    do k = 4, command_argument_count()
        call get_command_argument(k, arg)
        read(arg, *) x
        call flux_surface(eq, x, surf)
        call evaluate_batch_splines_1d_der(kin, x, y, dy)
        ! y = n_e, T_e, T_i, n_i, Z_eff; d/dpsi = d/dpsi_N / (psi_b - psi_axis)
        jb = sauter_jdotb(surf%f, surf%f_trap, surf%q, surf%r_geo, surf%epsilon, y(5), 1.0_dp, &
            y(1), y(4), y(2), y(3), dy(1) / dpsi, dy(4) / dpsi, dy(2) / dpsi, dy(3) / dpsi)
        write(buf, '(F0.4)') x
        print '(A,2(",",ES24.16))', 'jdotb@'//trim(buf), jb, rel * abs(jb)
        write(error_unit, '(F0.4,10(",",ES14.6))') x, surf%f_trap, surf%epsilon, surf%r_geo, surf%q, &
            sauter_nu_e(surf%q, surf%r_geo, surf%epsilon, y(1), y(2), y(5)), &
            sauter_nu_i(surf%q, surf%r_geo, surf%epsilon, y(4), y(3), 1.0_dp), jb, &
            surface_jdotb(eq, surf), surf%f, surf%b2
    end do

end program tiago_bootstrap
