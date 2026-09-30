module tiago_reconstruction_plots
    !! PNG figures of a reconstruction (fortplot): measured vs. modelled signals,
    !! normalized residuals, convergence, profiles with 1-sigma bands, boundary
    !! cross-sections, and the Jacobian check.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortplot, only: figure_t
    use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_get_var, nf90_nowrite, &
        nf90_inq_dimid, nf90_inquire_dimension, nf90_noerr
    implicit none
    private

    real(dp), parameter :: c_fit(3) = [0.12_dp, 0.47_dp, 0.71_dp]
    real(dp), parameter :: c_initial(3) = [0.6_dp, 0.6_dp, 0.6_dp]
    real(dp), parameter :: c_truth(3) = [0.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: c_meas(3) = [0.85_dp, 0.37_dp, 0.0_dp]

    type, public :: equilibrium_profiles_t
        logical :: present = .false.
        integer :: nfp = 1
        real(dp), allocatable :: s(:), pres(:), iota(:), jcurv(:)
        real(dp), allocatable :: xm(:), xn(:), rmnc(:), zmns(:)
    end type equilibrium_profiles_t

    public :: read_profiles, plot_signals, plot_residuals, plot_convergence
    public :: plot_profiles, plot_boundary, plot_jacobian_check

contains

    subroutine read_profiles(path, prof)
        character(len=*), intent(in) :: path
        type(equilibrium_profiles_t), intent(out) :: prof
        integer :: ncid, ns, mn, i
        real(dp), allocatable :: buf(:, :), phi(:)

        if (nf90_open(path, nf90_nowrite, ncid) /= nf90_noerr) return
        ns = dim(ncid, 'radius')
        mn = dim(ncid, 'mn_mode')
        call scalar(ncid, 'nfp', prof%nfp)
        allocate(phi(ns), prof%pres(ns), prof%iota(ns), prof%jcurv(ns))
        call get1(ncid, 'phi', phi)
        call get1(ncid, 'presf', prof%pres)
        call get1(ncid, 'iotaf', prof%iota)
        call get1(ncid, 'jcurv', prof%jcurv)
        prof%s = phi / phi(ns)
        allocate(prof%xm(mn), prof%xn(mn), buf(mn, ns))
        call get1(ncid, 'xm', prof%xm)
        call get1(ncid, 'xn', prof%xn)
        if (nf90_inq_varid(ncid, 'rmnc', i) == nf90_noerr) then
            i = nf90_get_var(ncid, i, buf)
            prof%rmnc = buf(:, ns)
        end if
        if (nf90_inq_varid(ncid, 'zmns', i) == nf90_noerr) then
            i = nf90_get_var(ncid, i, buf)
            prof%zmns = buf(:, ns)
        end if
        i = nf90_close(ncid)
        prof%present = .true.
    contains
        integer function dim(ncid, name)
            integer, intent(in) :: ncid
            character(len=*), intent(in) :: name
            integer :: id, status
            dim = 0
            if (nf90_inq_dimid(ncid, name, id) /= nf90_noerr) return
            status = nf90_inquire_dimension(ncid, id, len=dim)
        end function dim
        subroutine scalar(ncid, name, value)
            integer, intent(in) :: ncid
            character(len=*), intent(in) :: name
            integer, intent(inout) :: value
            integer :: id, status
            if (nf90_inq_varid(ncid, name, id) == nf90_noerr) status = nf90_get_var(ncid, id, value)
        end subroutine scalar
        subroutine get1(ncid, name, values)
            integer, intent(in) :: ncid
            character(len=*), intent(in) :: name
            real(dp), intent(inout) :: values(:)
            integer :: id, status
            values = 0.0_dp
            if (nf90_inq_varid(ncid, name, id) == nf90_noerr) status = nf90_get_var(ncid, id, values)
        end subroutine get1
    end subroutine read_profiles

    subroutine plot_signals(path, title, ylabel, measured, sigma, initial, fitted)
        !! Measured +- sigma (band) and the initial and fitted model per signal.
        character(len=*), intent(in) :: path, title, ylabel
        real(dp), intent(in) :: measured(:), sigma(:), initial(:), fitted(:)
        type(figure_t) :: fig
        real(dp), allocatable :: index(:)
        integer :: i

        if (size(measured) == 0) return
        index = [(real(i, dp), i = 1, size(measured))]
        call fig%initialize(900, 500)
        call fig%add_fill_between(index, measured - 2 * sigma, measured + 2 * sigma, &
            color='orange', alpha=0.35_dp)
        call fig%add_plot(index, measured, label='measured (band: 2 sigma)', color=c_meas)
        call fig%scatter(index, initial, label='initial model', color=c_initial, marker='s')
        call fig%scatter(index, fitted, label='fitted model', color=c_fit, marker='o')
        call fig%set_xlabel('signal index')
        call fig%set_ylabel(ylabel)
        call fig%set_title(title)
        call fig%legend()
        call fig%savefig(path)
    end subroutine plot_signals

    subroutine plot_residuals(path, initial, fitted)
        !! Normalized residuals (S - m) / sigma before and after the fit.
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: initial(:), fitted(:)
        type(figure_t) :: fig
        real(dp), allocatable :: index(:)
        integer :: i

        if (size(fitted) == 0) return
        index = [(real(i, dp), i = 1, size(fitted))]
        call fig%initialize(900, 450)
        call fig%scatter(index, initial, label='initial', color=c_initial, marker='s')
        call fig%scatter(index, fitted, label='fitted', color=c_fit, marker='o')
        call fig%axhline(1.0_dp, color='gray', linestyle='--')
        call fig%axhline(-1.0_dp, color='gray', linestyle='--')
        call fig%set_xlabel('residual index (measurements, then consistency points)')
        call fig%set_ylabel('(model - measured) / sigma')
        call fig%set_title('Normalized residuals')
        call fig%legend()
        call fig%savefig(path)
    end subroutine plot_residuals

    subroutine plot_convergence(path, history)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: history(:)
        type(figure_t) :: fig
        real(dp), allocatable :: index(:), best(:)
        integer :: i

        if (size(history) == 0) return
        index = [(real(i, dp), i = 1, size(history))]
        allocate(best(size(history)))
        best(1) = history(1)
        do i = 2, size(history)
            best(i) = min(best(i - 1), history(i))
        end do
        call fig%initialize(800, 450)
        call fig%scatter(index, max(history, tiny(1.0_dp)), label='every evaluation', &
            color=c_initial, marker='o')
        call fig%add_plot(index, max(best, tiny(1.0_dp)), label='best so far', color=c_fit)
        call fig%set_yscale('log')
        call fig%set_xlabel('residual evaluation (one VMEC++ solve each)')
        call fig%set_ylabel('chi^2')
        call fig%set_title('Levenberg-Marquardt convergence')
        call fig%legend()
        call fig%savefig(path)
    end subroutine plot_convergence

    subroutine plot_profiles(path, initial, fitted, truth, pres_sigma)
        !! Pressure (with a 1-sigma band when pres_sigma is given), iota and the
        !! toroidal current density, initial / fitted / true.
        character(len=*), intent(in) :: path
        type(equilibrium_profiles_t), intent(in) :: initial, fitted, truth
        real(dp), intent(in), optional :: pres_sigma(:)
        type(figure_t) :: fig
        character(len=len(path) + 16) :: name

        call fig%initialize(800, 500)
        if (present(pres_sigma)) then
            call fig%add_fill_between(fitted%s, (fitted%pres - pres_sigma) / 1.0e3_dp, &
                (fitted%pres + pres_sigma) / 1.0e3_dp, color='steelblue', alpha=0.3_dp)
        end if
        if (truth%present) call fig%add_plot(truth%s, truth%pres / 1.0e3_dp, label='truth', &
            linestyle='--', color=c_truth)
        if (initial%present) call fig%add_plot(initial%s, initial%pres / 1.0e3_dp, &
            label='initial', color=c_initial)
        call fig%add_plot(fitted%s, fitted%pres / 1.0e3_dp, label='fitted (band: 1 sigma)', &
            color=c_fit)
        call fig%set_xlabel('s (normalized toroidal flux)')
        call fig%set_ylabel('pressure [kPa]')
        call fig%set_title('Pressure profile')
        call fig%legend()
        name = path(:len(path) - 4)//'_pressure.png'
        call fig%savefig(trim(name))

        call fig%initialize(800, 500)
        if (truth%present) call fig%add_plot(truth%s, truth%iota, label='truth', &
            linestyle='--', color=c_truth)
        if (initial%present) call fig%add_plot(initial%s, initial%iota, label='initial', &
            color=c_initial)
        call fig%add_plot(fitted%s, fitted%iota, label='fitted', color=c_fit)
        call fig%set_xlabel('s (normalized toroidal flux)')
        call fig%set_ylabel('rotational transform iota')
        call fig%set_title('Rotational transform')
        call fig%legend()
        name = path(:len(path) - 4)//'_iota.png'
        call fig%savefig(trim(name))

        call fig%initialize(800, 500)
        if (truth%present) call fig%add_plot(truth%s, truth%jcurv / 1.0e3_dp, label='truth', &
            linestyle='--', color=c_truth)
        if (initial%present) call fig%add_plot(initial%s, initial%jcurv / 1.0e3_dp, &
            label='initial', color=c_initial)
        call fig%add_plot(fitted%s, fitted%jcurv / 1.0e3_dp, label='fitted', color=c_fit)
        call fig%set_xlabel('s (normalized toroidal flux)')
        call fig%set_ylabel('toroidal current density <j.grad phi>-like jcurv [kA/m^2]')
        call fig%set_title('Toroidal current density')
        call fig%legend()
        name = path(:len(path) - 4)//'_current.png'
        call fig%savefig(trim(name))
    end subroutine plot_profiles

    subroutine plot_boundary(path, initial, fitted, truth)
        !! Boundary cross-sections at phi = 0 and phi = pi / nfp.
        character(len=*), intent(in) :: path
        type(equilibrium_profiles_t), intent(in) :: initial, fitted, truth
        type(figure_t) :: fig
        integer :: plane
        real(dp) :: phi
        character(len=16) :: tag

        call fig%initialize(700, 700)
        do plane = 0, 1
            phi = plane * acos(-1.0_dp) / real(fitted%nfp, dp)
            tag = merge(' phi=0     ', ' phi=pi/nfp', plane == 0)
            if (truth%present) call curve(truth, 'truth'//trim(tag), c_truth, '--')
            if (initial%present) call curve(initial, 'initial'//trim(tag), c_initial, ':')
            call curve(fitted, 'fitted'//trim(tag), c_fit, '-')
        end do
        call fig%set_xlabel('R [m]')
        call fig%set_ylabel('Z [m]')
        call fig%set_title('Plasma boundary')
        call fig%legend()
        call fig%savefig(path)
    contains
        subroutine curve(prof, label, color, style)
            type(equilibrium_profiles_t), intent(in) :: prof
            character(len=*), intent(in) :: label, style
            real(dp), intent(in) :: color(3)
            real(dp) :: r(201), z(201), theta
            integer :: i
            do i = 1, 201
                theta = 2.0_dp * acos(-1.0_dp) * real(i - 1, dp) / 200.0_dp
                r(i) = sum(prof%rmnc * cos(prof%xm * theta - prof%xn * phi))
                z(i) = sum(prof%zmns * sin(prof%xm * theta - prof%xn * phi))
            end do
            call fig%add_plot(r, z, label=label, color=color, linestyle=style)
        end subroutine curve
    end subroutine plot_boundary

    subroutine plot_jacobian_check(path, adjoint, finite_difference)
        !! |J| from the adjoint vs. central differences of re-solved equilibria.
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: adjoint(:), finite_difference(:)
        type(figure_t) :: fig
        real(dp), allocatable :: a(:), f(:)
        real(dp) :: lo, hi
        logical, allocatable :: keep(:)

        keep = abs(finite_difference) > 1.0e-12_dp * maxval(abs(finite_difference))
        a = log10(abs(pack(adjoint, keep)) + tiny(1.0_dp))
        f = log10(abs(pack(finite_difference, keep)) + tiny(1.0_dp))
        lo = min(minval(a), minval(f))
        hi = max(maxval(a), maxval(f))
        call fig%initialize(650, 600)
        call fig%add_plot([lo, hi], [lo, hi], label='equal', color=c_initial, linestyle='--')
        call fig%scatter(f, a, label='Jacobian entries', color=c_fit, marker='o')
        call fig%set_xlabel('log10 |dr/dx|, central differences of VMEC++ solves')
        call fig%set_ylabel('log10 |dr/dx|, VMEC++ adjoint + Tiago')
        call fig%set_title('Jacobian check')
        call fig%legend()
        call fig%savefig(path)
    end subroutine plot_jacobian_check
end module tiago_reconstruction_plots
