module temporal_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_local_posterior, only: local_posterior_t
    use fortnum_rng, only: rng_t
    implicit none
    type, extends(local_posterior_t) :: temporal_t
        integer :: nt, n, m
        real(dp), allocatable :: design(:, :), y(:), sigma(:), rows(:, :), rhs(:)
        real(dp) :: normalization
    contains
        procedure :: activate
        procedure :: residual
        procedure :: jacobian
        procedure :: n_residuals
        procedure :: log_normalization
        procedure :: proper_prior
        procedure :: draw_prior
    end type
contains
    subroutine activate(self, drift_sd, drift_rho)
        class(temporal_t), intent(inout) :: self
        real(dp), intent(in) :: drift_sd, drift_rho
        real(dp) :: sd, rho, innovation
        integer :: i, j, k, offset
        if (allocated(self%rows)) deallocate (self%rows, self%rhs)
        allocate (self%rows(self%m + self%n, self%n), self%rhs(self%m + self%n))
        self%rows = 0
        self%rhs = 0
        do i = 1, self%m
            self%rows(i, :) = self%design(i, :)/self%sigma(i)
            self%rhs(i) = self%y(i)/self%sigma(i)
        end do
        self%normalization = -0.5_dp*(self%m + self%n)*log(2*acos(-1.0_dp)) &
            - sum(log(self%sigma))
        do j = 0, 3
            sd = drift_sd
            rho = drift_rho
            offset = self%nt + 3 + (j - 1)*self%nt
            if (j == 0) then
                sd = 0.28_dp
                rho = 0.88_dp
                offset = 0
            end if
            innovation = sd*sqrt(1 - rho**2)
            do i = 1, self%nt
                k = offset + i
                if (i == 1) then
                    self%rows(self%m + k, k) = 1/sd
                else
                    self%rows(self%m + k, k) = 1/innovation
                    self%rows(self%m + k, k - 1) = -rho/innovation
                end if
            end do
            self%normalization = self%normalization - log(sd) &
                - (self%nt - 1)*log(innovation)
        end do
        do i = self%nt + 1, self%nt + 3
            self%rows(self%m + i, i) = 1/0.35_dp
        end do
        self%normalization = self%normalization - 3*log(0.35_dp)
    end subroutine
    subroutine residual(self, theta, r, ok)
        class(temporal_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), intent(out) :: r(:)
        logical, intent(out) :: ok
        r = matmul(self%rows, theta) - self%rhs
        ok = .true.
    end subroutine
    subroutine jacobian(self, theta, jac)
        class(temporal_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), allocatable, intent(out) :: jac(:, :)
        jac = self%rows
    end subroutine
    integer function n_residuals(self) result(n)
        class(temporal_t), intent(in) :: self
        n = self%m + self%n
    end function
    real(dp) function log_normalization(self) result(c)
        class(temporal_t), intent(in) :: self
        c = self%normalization
    end function
    logical function proper_prior(self) result(proper)
        class(temporal_t), intent(in) :: self
        proper = .true.
    end function
    subroutine draw_prior(self, rng, theta, ok)
        class(temporal_t), intent(in) :: self
        type(rng_t), intent(inout) :: rng
        real(dp), intent(out) :: theta(:)
        logical, intent(out) :: ok
        ! This fixture supplies deterministic starts; prior drawing is not exercised.
        theta = 0
        ok = .false.
    end subroutine
end module

program test_temporal_inference
    use temporal_fixture
    use tiago_local_mixture, only: mode_start_t, mixture_options_t, &
        local_mixture_t, multi_start_map
    use fortnum_status, only: fortnum_status_t
    implicit none
    type(temporal_t) :: model
    type(mode_start_t) :: starts(8)
    type(mixture_options_t) :: options
    type(local_mixture_t) :: mixture
    type(fortnum_status_t) :: status
    character(len=1024) :: directory
    integer :: u, i
    call get_command_argument(1, directory)
    open (newunit=u, file=trim(directory)//'/input.txt', status='old')
    read (u, *) model%nt, model%m
    model%n = 4*model%nt + 3
    allocate (model%design(model%m, model%n), model%y(model%m), model%sigma(model%m))
    do i = 1, model%m
        read (u, *) model%y(i), model%sigma(i), model%design(i, :)
    end do
    close (u)
    do i = 1, 8
        allocate (starts(i)%theta(model%n))
        starts(i)%theta = 0.05_dp*(i - 1)
        write (starts(i)%family, '(A,I0)') 'covariance-', mod(i - 1, 4) + 1
        starts(i)%log_family_prior = log(0.25_dp)
    end do
    options%psis_draws = 20
    options%decrease_tolerance = 1.e-12_dp
    call multi_start_map(model, starts, model%n, trim(directory), options, mixture, &
        status, prepare)
    if (status%code /= 0) stop 1
    if (.not. mixture%evidence_valid) stop 2
    if (size(mixture%modes) /= 4) stop 3
    if (any(mixture%attempts%mode == 0)) stop 4
contains
    subroutine prepare(index, posterior, ok)
        integer, intent(in) :: index
        class(local_posterior_t), intent(inout) :: posterior
        logical, intent(out) :: ok
        real(dp), parameter :: sd(4) = [0.12_dp, 0.12_dp, 0.30_dp, 0.30_dp]
        real(dp), parameter :: rho(4) = [0.30_dp, 0.94_dp, 0.30_dp, 0.94_dp]
        integer :: family
        family = mod(index - 1, 4) + 1
        ok = .false.
        select type (posterior)
            type is (temporal_t)
            call posterior%activate(sd(family), rho(family))
            ok = .true.
        end select
    end subroutine
end program
