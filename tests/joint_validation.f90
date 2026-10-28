module joint_validation
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use tiago_joint_model, only: joint_model_t
    use tiago_reconstruction, only: rec, invalidate_caches
    use tiago_recorded_diagnostic, only: recorded_diagnostic_t
    use tiago_dispersion_atan, only: dispersion_atan
    use tiago_dispersion_iq, only: dispersion_iq
    use tiago_kinetic_profiles, only: N_PROFILES, P_NE, P_TE
    use fortnum_rng, only: rng_t, rng_seed, rng_normal
    use fortnum_status, only: fortnum_status_t
    implicit none
contains
    subroutine calibration_start(model, theta, residual)
        type(joint_model_t), intent(in) :: model
        real(dp), intent(inout) :: theta(:), residual(:)
        real(dp), allocatable :: prediction(:), sigma(:), a(:), b(:)
        real(dp) :: ratio, offset, phase(1), measured(1), row(1, 6), iq(2), iqrow(2, 5), shift
        integer :: d, first, last, k, j, nloop
        do d = 1, 3
            first = model%n_mag + model%n_fb + model%n_rec_prior + model%diags(d)%row_offset + 1
            last = first + size(model%diags(d)%data) - 1
            prediction = model%diags(d)%data + matmul(model%diags(d)%chol, residual(first:last))
            sigma = [(model%diags(d)%chol(j, j), j=1, size(prediction))]
            k = model%diags(d)%eta_offset
            select case (d)
            case (1, 3)
                offset = 0
                if (d == 3) offset = .02_dp
                a = (prediction - offset)/sigma; b = (model%diags(d)%data - offset)/sigma
                ratio = dot_product(a, b)/dot_product(a, a)
                if (d == 1) theta(k + 1) = (ratio - 1)/.2_dp
                if (d == 3) theta(k + 1) = (ratio - 1)/.1_dp
                residual(first:last) = (offset + ratio*(prediction - offset) - model%diags(d)%data)/sigma
            case (2)
                shift = 0
                do j = 1, size(prediction), 2
                    call dispersion_atan(prediction(j), prediction(j + 1), 1._dp, 1._dp, &
                                         0._dp, 0._dp, phase, row)
                    call dispersion_atan(model%diags(d)%data(j), model%diags(d)%data(j + 1), &
                                         1._dp, 1._dp, 0._dp, 0._dp, measured, row)
                    shift = shift + measured(1) - phase(1)
                end do
                shift = shift/(size(prediction)/2)
                theta(k + 6) = shift/.2_dp
                do j = 1, size(prediction), 2
                    call dispersion_atan(prediction(j), prediction(j + 1), 1._dp, 1._dp, &
                                         0._dp, 0._dp, phase, row)
                    call dispersion_iq(phase(1) + shift, 1._dp, 1._dp, 0._dp, 0._dp, iq, iqrow)
                    prediction(j:j + 1) = iq
                end do
                residual(first:last) = (prediction - model%diags(d)%data)/sigma
            end select
        end do
        nloop = model%n_mag/3
        k = model%mag_calibration_index
        theta(k + 3) = sum(rec%meas_value(2*nloop + 1:) - rec%meas_value(:nloop))/(nloop*2.e-5_dp)
        do j = 1, model%n_mag
            residual(j) = residual(j) + theta(k + 3)*1.e-5_dp*model%mag_time(j)/model%mag_chol(j, j)
        end do
        first = model%n_residuals() - (model%n - model%n_rec) + 1
        residual(first:) = theta(model%n_rec + 1:)
    end subroutine

    subroutine write_calibrations(model, theta, cov, truth, path)
        type(joint_model_t), intent(in) :: model
        real(dp), intent(in) :: theta(:), cov(:, :), truth(:)
        character(*), intent(in) :: path
        character(32), parameter :: names(3) = [character(32) :: &
                                                'Thomson absolute scale', 'Interferometer phase offset', 'ECE gain']
        character(8), parameter :: units(3) = [character(8) :: '1', 'rad', 'V/eV']
        integer, parameter :: index(3) = [1, 6, 1]
        real(dp), allocatable :: eta(:), physical(:), injected(:), g(:)
        integer :: d, k, ncal, first, last, u
        allocate (g(model%n))
        open (newunit=u, file=path, status='replace')
        write (u, '(A)') 'parameter,unit,truth,map,sd'
        do d = 1, 3
            select type (instrument => model%diags(d)%d)
            class is (recorded_diagnostic_t)
                ncal = instrument%calibration%n_params()
                call model%unpack(theta, eta=eta, d=d)
                physical = instrument%calibration%physical(eta(:ncal))
                call model%unpack(truth, eta=eta, d=d)
                injected = instrument%calibration%physical(eta(:ncal))
                first = model%diags(d)%eta_offset - model%n_rec - model%n_prof + 1
                last = first + ncal - 1
                k = index(d)
                g = 0
                g(model%n_rec + model%n_prof + 1:model%n_rec + model%n_prof + model%n_eta) = &
                    matmul(instrument%calibration%factor(k, :), model%eta_factor(first:last, :))
                write (u, '(A,",",A,3(",",ES25.17E3))') trim(names(d)), trim(units(d)), &
                    injected(k), physical(k), sqrt(dot_product(g, matmul(cov, g)))
            end select
        end do
        k = model%mag_calibration_index
        physical = model%mag_calibration%physical(theta(k + 1:k + 3))
        injected = model%mag_calibration%physical(truth(k + 1:k + 3))
        g = 0; g(k + 1:k + 3) = model%mag_calibration%factor(3, :)
        write (u, '(A,3(",",ES25.17E3))') 'Magnetic integrator drift,Wb/s', &
            injected(3), physical(3), sqrt(dot_product(g, matmul(cov, g)))
        close (u)
    end subroutine

    subroutine add_noise(model, seed)
        type(joint_model_t), intent(inout)::model
        integer, intent(in)::seed
        type(rng_t)::rng
        type(fortnum_status_t)::status
        real(dp), allocatable::e(:)
        real(dp)::normal
        integer::i, d
        call rng_seed(rng, int(seed, int64), status)
        do i = 1, model%n_mag
            call rng_normal(rng, normal)
            rec%meas_value(i) = rec%meas_value(i) + rec%meas_sigma(i)*normal
        end do
        do d = 1, size(model%diags)
            allocate (e(size(model%diags(d)%data)))
            do i = 1, size(e)
                call rng_normal(rng, e(i))
            end do
            model%diags(d)%data = model%diags(d)%data + matmul(model%diags(d)%chol, e)
            deallocate (e)
        end do
        call invalidate_caches()
    end subroutine

    subroutine write_profiles(model, theta, cov, truth, path, error)
        type(joint_model_t), intent(in)::model
        real(dp), intent(in)::theta(:), cov(:, :), truth(:)
        character(*), intent(in)::path
        real(dp), intent(out)::error
        integer, parameter::ng = 21
        real(dp)::s(ng), f(N_PROFILES, ng), ft(N_PROFILES, ng), &
                   seed(N_PROFILES, ng), bar(model%n_prof), sb(ng), variance, df(N_PROFILES, ng)
        integer::i, p, u, a, b
        a = model%n_rec + 1; b = model%n_rec + model%n_prof
        s = [(real(i - 1, dp)/(ng - 1), i=1, ng)]
        call model%prof%evaluate(theta(a:b), s, f, df)
        call model%prof%evaluate(truth(a:b), s, ft, df)
        open (newunit=u, file=path, status='replace')
        write (u, '(A)') 'profile,s,truth,map,sd'
        error = 0
        do p = P_NE, P_TE
            do i = 1, ng
                seed = 0; seed(p, i) = 1; bar = 0
                call model%prof%vjp(theta(a:b), s, seed, bar, sb)
                variance = dot_product(bar, matmul(cov(a:b, a:b), bar))
                write (u, '(I0,4(",",ES25.17E3))') p, s(i), ft(p, i), f(p, i), &
                    sqrt(max(0._dp, variance))
                error = error + ((f(p, i) - ft(p, i))/ft(p, i))**2
            end do
        end do
        close (u)
        error = sqrt(error/(2*ng))
    end subroutine

    subroutine read_start(path, theta)
        character(*), intent(in)::path
        real(dp), intent(out)::theta(:)
        integer::u, i, k
        real(dp)::truth, sd
        open (newunit=u, file=path, status='old')
        read (u, *)
        do i = 1, size(theta)
            read (u, *) k, truth, theta(i), sd
        end do
        close (u)
    end subroutine
end module
