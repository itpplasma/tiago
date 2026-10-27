module tiago_energy_consistency
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t, &
                                        Q_NE, Q_TE, Q_TI, Q_ZEFF
    use tiago_kinetic_pressure, only: kinetic_pressure
    use tiago_thermal_energy_density, only: thermal_energy_density
    use tiago_diamagnetic_energy, only: diamagnetic_energy
    implicit none
    private
    public :: energy_consistency_t
    real(dp), parameter :: charge = 1.602176634e-19_dp
    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 1.25663706212e-6_dp

    type, extends(profile_diagnostic_t) :: energy_consistency_t
        real(dp), allocatable :: points(:, :), volume_weights(:)
        real(dp) :: zimp = 6.0_dp
    contains
        procedure :: n_signals
        procedure :: n_nuisance
        procedure :: label
        procedure :: sample_points
        procedure :: predict
        procedure :: jvp
        procedure :: vjp
        procedure :: energies
    end type
contains
    integer function n_signals(self)
        class(energy_consistency_t), intent(in) :: self
        n_signals = 1
    end function
    integer function n_nuisance(self)
        class(energy_consistency_t), intent(in) :: self
        n_nuisance = 4
    end function
    function label(self, i) result(text)
        class(energy_consistency_t), intent(in) :: self
        integer, intent(in) :: i
        character(:), allocatable :: text
        text = 'Wdia-Wkin [J]'
    end function
    subroutine sample_points(self, eta, x, dx_deta)
        class(energy_consistency_t), intent(in) :: self
        real(dp), intent(in) :: eta(:)
        real(dp), allocatable, intent(out) :: x(:, :), dx_deta(:, :, :)
        x = self%points
        allocate (dx_deta(3, size(x, 2), 4))
        dx_deta = 0.0_dp
    end subroutine
    subroutine energies(self, sample, eta, wdia, wkin)
        class(energy_consistency_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: wdia, wkin
        real(dp) :: row(4), p, dp_dn, dp_dte, dp_dti, dp_dzf, dp_df, e, de
        integer :: i
        call diamagnetic_energy(eta(1), eta(2), eta(3), eta(4), mu0, pi, &
                                wdia, row(1), row(2), row(3), row(4))
        wkin = 0.0_dp
        !$omp parallel do default(shared) &
        !$omp private(i,p,dp_dn,dp_dte,dp_dti,dp_dzf,dp_df,e,de) &
        !$omp reduction(+:wkin) if(sample%n>1024)
        do i = 1, sample%n
            call kinetic_pressure(sample%q(Q_NE, i), sample%q(Q_TE, i), &
                            sample%q(Q_TI, i), sample%q(Q_ZEFF, i), 0.0_dp, self%zimp, &
                                  charge, p, dp_dn, dp_dte, dp_dti, dp_dzf, dp_df)
            call thermal_energy_density(p, e, de)
            wkin = wkin + self%volume_weights(i)*e
        end do
        !$omp end parallel do
    end subroutine
    subroutine predict(self, sample, eta, y)
        class(energy_consistency_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: wdia, wkin
        call self%energies(sample, eta, wdia, wkin)
        y(1) = wdia - wkin
    end subroutine
    subroutine jvp(self, sample, eta, q_dot, eta_dot, y_dot)
        class(energy_consistency_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        real(dp) :: w, row(4), p, p_row(4), df, e, de, dwkin
        integer :: i
        call diamagnetic_energy(eta(1), eta(2), eta(3), eta(4), mu0, pi, &
                                w, row(1), row(2), row(3), row(4))
        dwkin = 0.0_dp
        !$omp parallel do default(shared) private(i,p,p_row,df,e,de) &
        !$omp reduction(+:dwkin) if(sample%n>1024)
        do i = 1, sample%n
            call kinetic_pressure(sample%q(Q_NE, i), sample%q(Q_TE, i), &
                            sample%q(Q_TI, i), sample%q(Q_ZEFF, i), 0.0_dp, self%zimp, &
                                  charge, p, p_row(1), p_row(2), p_row(3), p_row(4), df)
            call thermal_energy_density(p, e, de)
            dwkin = dwkin + self%volume_weights(i)*de* &
                    dot_product(p_row, q_dot([Q_NE, Q_TE, Q_TI, Q_ZEFF], i))
        end do
        !$omp end parallel do
        y_dot(1) = dot_product(row, eta_dot) - dwkin
    end subroutine
    subroutine vjp(self, sample, eta, y_bar, q_bar, eta_bar)
        class(energy_consistency_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        real(dp) :: w, row(4), p, p_row(4), df, e, de
        integer :: i
        call diamagnetic_energy(eta(1), eta(2), eta(3), eta(4), mu0, pi, &
                                w, row(1), row(2), row(3), row(4))
        eta_bar = row*y_bar(1)
        q_bar = 0.0_dp
        !$omp parallel do default(shared) private(i,p,p_row,df,e,de) if(sample%n>1024)
        do i = 1, sample%n
            call kinetic_pressure(sample%q(Q_NE, i), sample%q(Q_TE, i), &
                            sample%q(Q_TI, i), sample%q(Q_ZEFF, i), 0.0_dp, self%zimp, &
                                  charge, p, p_row(1), p_row(2), p_row(3), p_row(4), df)
            call thermal_energy_density(p, e, de)
            q_bar([Q_NE, Q_TE, Q_TI, Q_ZEFF], i) = &
                -self%volume_weights(i)*de*p_row*y_bar(1)
        end do
        !$omp end parallel do
    end subroutine
end module
