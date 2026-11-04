program pressure_reference
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_kinetic_profiles, only: kinetic_profiles_t, P_NE, P_TE, P_TI
    implicit none
    type(kinetic_profiles_t) :: profiles
    integer, parameter :: n = 101
    real(dp) :: s(n), thermal(n), total(n), reference(n), amplitude, exponent
    real(dp), allocatable :: z(:)
    character(512) :: output
    integer :: k, unit

    call profiles%init([0._dp, .5_dp, 1._dp])
    call profiles%set_profile(P_NE, [(1.e19_dp, k=1, profiles%ncoef)], .false.)
    call profiles%set_profile(P_TE, [(800._dp, k=1, profiles%ncoef)], .false.)
    call profiles%set_profile(P_TI, [(600._dp, k=1, profiles%ncoef)], .false.)
    s = [(real(k-1, dp)/(n-1), k=1, n)]
    allocate (z(profiles%n_latent))
    z = 0
    call profiles%pressure(z, s, thermal)
    call profiles%set_nonthermal([0._dp, 2._dp], [2._dp, 8._dp])
    z = [z, .3_dp, -.4_dp]
    call profiles%pressure(z, s, total)
    amplitude = 2/(1 + exp(-.3_dp))
    exponent = 2 + 6/(1 + exp(.4_dp))
    reference = amplitude*thermal(1)*(1-s)**exponent
    if (maxval(abs(total-thermal-reference)) > 1.e-10_dp) &
        error stop 'native pressure differs from the independent reference'
    output = 'pressure-reference.csv'
    if (command_argument_count() > 0) call get_command_argument(1, output)
    open (newunit=unit, file=trim(output), status='replace')
    write (unit, '(A)') 's,thermal_Pa,total_Pa,additional_Pa,formula_Pa'
    do k = 1, n
        write (unit, '(ES22.14,4(A,ES22.14))') s(k), ',', thermal(k), ',', total(k), &
            ',', total(k)-thermal(k), ',', reference(k)
    end do
    close (unit)
end program
