module tiago_soft_xray
    !! Line-integrated soft X-ray (SXR) emission seen by pinhole cameras
    !! (docs/soft-xray.md). Channel c views the straight chord from its pinhole
    !! along detector -> pinhole and records
    !!
    !!     y_c = G_c / (4 pi) int eps dl     [W],
    !!
    !! G_c the etendue [m^2 sr] (default 4 pi: y_c = int eps dl [W m^-2]) and
    !! eps the filtered emissivity [W m^-3], zero outside the last closed flux
    !! surface. Two emissivity models:
    !!
    !!   SXR_PROFILE_MODEL         eps = exp(u(s)), u a cubic B-spline in the flux
    !!                             label s with the kinetic-profile parametrisation
    !!                             and Gibbs/Matern prior (tiago_kinetic_profiles);
    !!                             its standard-normal latents are the nuisance
    !!                             parameters eta of this diagnostic. In a field
    !!                             with islands s is the provider's label (island
    !!                             plateaus of tiago_point_table).
    !!   SXR_BREMSSTRAHLUNG_MODEL  eps = C ne^2 Zeff Te^(-1/2) int t(E) exp(-E/Te) dE
    !!                             with C = 1.69e-38 W m^3 eV^-1/2 and the
    !!                             piecewise-linear transmission t(E) of the
    !!                             camera's filter and detector; no nuisances.
    !!
    !! Detector responsivity and amplifier gain belong to the recorded-signal
    !! calibration (calibrated_diagnostic_t with one log-gain group per camera).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t
    use fortnum_bspline, only: bspline_eval_deriv
    use tiago_profile_diagnostic, only: plasma_sample_t, Q_S, Q_NE, Q_TE, Q_ZEFF
    use tiago_sight_geometry, only: chord_t
    use tiago_line_integrals, only: chord_set_t
    use tiago_kinetic_profiles, only: kinetic_profiles_t, P_NE, N_PROFILES
    use tiago_sxr_bremsstrahlung, only: sxr_bremsstrahlung
    use tiago_sxr_filter_segment, only: sxr_filter_segment
    use tiago_sxr_filter_segment_small, only: sxr_filter_segment_small
    use tiago_sxr_filter_tail, only: sxr_filter_tail
    implicit none
    private

    integer, parameter, public :: SXR_PROFILE_MODEL = 1, SXR_BREMSSTRAHLUNG_MODEL = 2
    real(dp), parameter :: four_pi = 16.0_dp * atan(1.0_dp)
    real(dp), parameter :: small_width = 0.01_dp

    type, public :: sxr_filter_t
        !! Transmission t(E): 0 below energy(1), linear between the nodes,
        !! transmission(n) above energy(n). Default: t = 1 for E >= 0.
        real(dp), allocatable :: energy(:) !! [eV], strictly increasing, >= 0
        real(dp), allocatable :: transmission(:) !! in [0, 1]
    end type sxr_filter_t

    type, public :: pinhole_camera_t
        !! Detector array behind one pinhole, Cartesian [m].
        character(len=:), allocatable :: name
        real(dp) :: aperture(3) = 0.0_dp !! pinhole centre
        real(dp) :: aperture_normal(3) = 0.0_dp !! zero: normal incidence
        real(dp) :: aperture_area = 0.0_dp !! [m^2]; zero: no etendue
        real(dp), allocatable :: detector(:, :) !! (3, channels) detector centres
        real(dp) :: detector_normal(3) = 0.0_dp !! zero: normal incidence
        real(dp), allocatable :: detector_area(:) !! [m^2] per channel
        real(dp) :: length = 0.0_dp !! chord length from the pinhole [m]
        type(sxr_filter_t) :: filter
    end type pinhole_camera_t

    type, public, extends(chord_set_t) :: soft_xray_t
        integer :: model = SXR_PROFILE_MODEL
        type(kinetic_profiles_t) :: emissivity !! slot P_NE holds eps [W m^-3]
        real(dp), allocatable :: etendue(:) !! [m^2 sr] per chord
        integer, allocatable :: camera(:) !! filter index per chord
        type(sxr_filter_t), allocatable :: filters(:)
    contains
        procedure :: setup => sxr_setup
        procedure :: set_cameras => sxr_set_cameras
        procedure :: set_emissivity_basis => sxr_set_basis
        procedure :: set_emissivity_prior => sxr_set_prior
        procedure :: set_emissivity_table => sxr_set_table
        procedure :: n_nuisance => sxr_n_nuisance
        procedure :: label => sxr_label
        procedure :: predict => sxr_predict
        procedure :: jvp => sxr_jvp
        procedure :: vjp => sxr_vjp
    end type soft_xray_t

    public :: camera_chords, filter_response

contains

    subroutine camera_chords(camera, npanel, chords, etendue)
        !! Chords from the pinhole along detector -> pinhole, and the etendue
        !! A_d A_a cos(theta_d) cos(theta_a) / d^2 of each channel (pinhole
        !! and detector small compared with their distance d); 4 pi where an
        !! area is zero.
        type(pinhole_camera_t), intent(in) :: camera
        integer, intent(in) :: npanel
        type(chord_t), allocatable, intent(out) :: chords(:)
        real(dp), allocatable, intent(out) :: etendue(:)
        real(dp) :: u(3), d, cd, ca, ad
        character(len=16) :: number
        integer :: c, n

        n = size(camera%detector, 2)
        if (camera%length <= 0.0_dp) error stop 'soft X-ray: camera length must be positive'
        allocate(chords(n), etendue(n))
        do c = 1, n
            u = camera%aperture - camera%detector(:, c)
            d = norm2(u)
            if (d <= 0.0_dp) error stop 'soft X-ray: detector at the pinhole'
            u = u / d
            write(number, '(I0)') c
            chords(c)%label = camera%name//'_'//trim(number)
            chords(c)%origin = camera%aperture
            chords(c)%direction = u
            chords(c)%length = camera%length
            chords(c)%npanel = npanel
            chords(c)%ngauss = 4
            ad = 0.0_dp
            if (allocated(camera%detector_area)) ad = camera%detector_area(c)
            etendue(c) = four_pi
            if (ad <= 0.0_dp .or. camera%aperture_area <= 0.0_dp) cycle
            cd = incidence(camera%detector_normal, u)
            ca = incidence(camera%aperture_normal, u)
            etendue(c) = ad * camera%aperture_area * cd * ca / d**2
        end do
    end subroutine camera_chords

    pure real(dp) function incidence(normal, u)
        real(dp), intent(in) :: normal(3), u(3)
        incidence = 1.0_dp
        if (norm2(normal) > 0.0_dp) incidence = abs(dot_product(normal, u)) / norm2(normal)
    end function incidence

    subroutine sxr_set_cameras(self, cameras, npanel)
        !! Chords, etendues and filters of all cameras, in camera order.
        class(soft_xray_t), intent(inout) :: self
        type(pinhole_camera_t), intent(in) :: cameras(:)
        integer, intent(in) :: npanel
        type(chord_t), allocatable :: all(:), chords(:)
        real(dp), allocatable :: g(:), etendue(:)
        integer, allocatable :: index(:)
        type(sxr_filter_t), allocatable :: filters(:)
        integer :: i

        allocate(all(0), g(0), index(0), filters(size(cameras)))
        do i = 1, size(cameras)
            call camera_chords(cameras(i), npanel, chords, etendue)
            all = [all, chords]
            g = [g, etendue]
            index = [index, spread(i, 1, size(chords))]
            filters(i) = cameras(i)%filter
        end do
        call self%setup(all, index, g, filters)
    end subroutine sxr_set_cameras

    subroutine sxr_setup(self, chords, camera, etendue, filters)
        !! Chords with optional camera index (default 1), etendue [m^2 sr]
        !! (default 4 pi) and one filter per camera (default t = 1).
        class(soft_xray_t), intent(inout) :: self
        type(chord_t), intent(in) :: chords(:)
        integer, intent(in), optional :: camera(:)
        real(dp), intent(in), optional :: etendue(:)
        type(sxr_filter_t), intent(in), optional :: filters(:)
        integer :: i, n

        call self%chord_set_t%init(chords)
        n = size(chords)
        self%camera = [(1, i = 1, n)]
        if (present(camera)) self%camera = camera
        self%etendue = [(four_pi, i = 1, n)]
        if (present(etendue)) self%etendue = etendue
        if (size(self%camera) /= n .or. size(self%etendue) /= n) &
            error stop 'soft X-ray: one camera index and etendue per chord'
        if (any(self%etendue <= 0.0_dp)) error stop 'soft X-ray: etendue must be positive'
        if (any(self%camera < 1)) error stop 'soft X-ray: camera index below 1'
        if (allocated(self%filters)) deallocate(self%filters)
        allocate(self%filters(maxval([self%camera, 0])))
        if (present(filters)) then
            if (size(filters) < size(self%filters)) error stop 'soft X-ray: one filter per camera'
            self%filters = filters(:size(self%filters))
        end if
        do i = 1, size(self%filters)
            call check_filter(self%filters(i))
        end do
    end subroutine sxr_setup

    subroutine check_filter(filter)
        type(sxr_filter_t), intent(inout) :: filter
        integer :: n
        if (.not. allocated(filter%energy)) then
            filter%energy = [0.0_dp]
            filter%transmission = [1.0_dp]
        end if
        n = size(filter%energy)
        if (n < 1 .or. size(filter%transmission) /= n) &
            error stop 'soft X-ray: filter needs one transmission per energy'
        if (.not. all(ieee_is_finite(filter%energy)) .or. &
            .not. all(ieee_is_finite(filter%transmission))) &
            error stop 'soft X-ray: non-finite filter table'
        if (filter%energy(1) < 0.0_dp) error stop 'soft X-ray: negative filter energy'
        if (n > 1) then
            if (any(filter%energy(2:) <= filter%energy(:n - 1))) &
                error stop 'soft X-ray: filter energies must increase'
        end if
        if (any(filter%transmission < 0.0_dp) .or. any(filter%transmission > 1.0_dp)) &
            error stop 'soft X-ray: filter transmission outside [0, 1]'
    end subroutine check_filter

    pure subroutine filter_response(filter, te, f, df_dte)
        !! F(Te) = int t(E) exp(-E / Te) dE / Te and dF/dTe; F = 1 without filter.
        type(sxr_filter_t), intent(in) :: filter
        real(dp), intent(in) :: te
        real(dp), intent(out) :: f, df_dte
        real(dp) :: g, dg, h
        integer :: j, n

        n = size(filter%energy)
        call sxr_filter_tail(filter%energy(n), filter%transmission(n), te, f, df_dte)
        do j = 1, n - 1
            h = filter%energy(j + 1) - filter%energy(j)
            if (h < small_width * te) then
                call sxr_filter_segment_small(filter%energy(j), h, filter%transmission(j), &
                    filter%transmission(j + 1), te, g, dg)
            else
                call sxr_filter_segment(filter%energy(j), h, filter%transmission(j), &
                    filter%transmission(j + 1), te, g, dg)
            end if
            f = f + g
            df_dte = df_dte + dg
        end do
    end subroutine filter_response

    subroutine sxr_set_basis(self, breakpoints)
        !! Cubic B-spline basis of log eps on breakpoints in s; selects the
        !! emissivity-profile model. Unit: 1 W m^-3.
        class(soft_xray_t), intent(inout) :: self
        real(dp), intent(in) :: breakpoints(:)
        type(kinetic_profiles_t) :: fresh

        self%emissivity = fresh
        call self%emissivity%init(breakpoints)
        self%emissivity%p(P_NE)%unit = 1.0_dp
        self%model = SXR_PROFILE_MODEL
    end subroutine sxr_set_basis

    subroutine sxr_set_prior(self, mean_values, sigma_f, l_core, l_edge, s0, sw, kernel)
        !! Latent emissivity profile: prior mean eps [W m^-3] at the Greville
        !! abscissae (self%emissivity%greville) and the kernel of log eps.
        class(soft_xray_t), intent(inout) :: self
        real(dp), intent(in) :: mean_values(:), sigma_f, l_core, l_edge, s0, sw
        integer, intent(in), optional :: kernel
        call self%emissivity%set_profile(P_NE, mean_values, .true., sigma_f, l_core, l_edge, &
            s0, sw, kernel)
    end subroutine sxr_set_prior

    subroutine sxr_set_table(self, s, values, ierr)
        !! Fixed emissivity profile, least-squares fit of log eps to a table.
        class(soft_xray_t), intent(inout) :: self
        real(dp), intent(in) :: s(:), values(:)
        integer, intent(out) :: ierr
        call self%emissivity%set_table(P_NE, s, values, ierr)
    end subroutine sxr_set_table

    integer function sxr_n_nuisance(self)
        class(soft_xray_t), intent(in) :: self
        sxr_n_nuisance = 0
        if (self%model == SXR_PROFILE_MODEL) sxr_n_nuisance = self%emissivity%n_latent
    end function sxr_n_nuisance

    function sxr_label(self, i) result(text)
        class(soft_xray_t), intent(in) :: self
        integer, intent(in) :: i
        character(len=:), allocatable :: text
        text = 'sxr '//self%chords(i)%label
    end function sxr_label

    subroutine emission(self, sample, eta, eps, de)
        !! eps at every sample point and its partial derivatives
        !! de = d eps / d (s, ne, Te, Zeff); zero outside the plasma.
        class(soft_xray_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: eps(:), de(:, :)
        real(dp), allocatable :: f(:, :), df(:, :)
        real(dp) :: filt, dfilt
        integer :: k

        eps = 0.0_dp
        de = 0.0_dp
        if (self%model == SXR_PROFILE_MODEL) then
            if (self%emissivity%ncoef == 0) &
                error stop 'soft X-ray: call set_emissivity_basis before the profile model'
            allocate(f(N_PROFILES, size(eps)), df(N_PROFILES, size(eps)))
            call self%emissivity%evaluate(eta, sample%q(Q_S, :), f, df)
            where (sample%inside)
                eps = f(P_NE, :)
                de(1, :) = df(P_NE, :)
            end where
            return
        end if
        !$omp parallel do default(shared) private(k, filt, dfilt) if(size(eps)>4096)
        do k = 1, size(eps)
            if (.not. sample%inside(k)) cycle
            if (sample%q(Q_TE, k) <= 0.0_dp) cycle
            call filter_response(self%filters(self%camera(self%owner(k))), sample%q(Q_TE, k), &
                filt, dfilt)
            call sxr_bremsstrahlung(sample%q(Q_NE, k), sample%q(Q_TE, k), sample%q(Q_ZEFF, k), &
                filt, dfilt, eps(k), de(2, k), de(3, k), de(4, k))
        end do
        !$omp end parallel do
    end subroutine emission

    subroutine sxr_predict(self, sample, eta, y)
        class(soft_xray_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: eps(size(self%w)), de(4, size(self%w))
        call emission(self, sample, eta, eps, de)
        call chord_sums(self, self%w * eps, y)
    end subroutine sxr_predict

    subroutine chord_sums(self, v, y)
        !! y_c = G_c / (4 pi) sum over the points of chord c of v.
        class(soft_xray_t), intent(in) :: self
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: y(:)
        integer :: c
        do c = 1, size(self%chords)
            y(c) = self%etendue(c) / four_pi * sum(v(self%first(c):self%first(c + 1) - 1))
        end do
    end subroutine chord_sums

    subroutine sxr_jvp(self, sample, eta, q_dot, eta_dot, y_dot)
        class(soft_xray_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
        real(dp), intent(out) :: y_dot(:)
        real(dp) :: eps(size(self%w)), de(4, size(self%w)), v(size(self%w))
        real(dp), allocatable :: c_dot(:), b(:, :)
        type(fortnum_status_t) :: st
        integer :: k

        call emission(self, sample, eta, eps, de)
        v = de(1, :) * q_dot(Q_S, :) + de(2, :) * q_dot(Q_NE, :) + de(3, :) * q_dot(Q_TE, :) &
            + de(4, :) * q_dot(Q_ZEFF, :)
        if (self%n_nuisance() > 0) then
            c_dot = matmul(self%emissivity%p(P_NE)%chol, eta_dot)
            allocate(b(0:1, size(c_dot)))
            do k = 1, size(v)
                if (eps(k) == 0.0_dp) cycle
                call bspline_eval_deriv(self%emissivity%basis, &
                    min(max(sample%q(Q_S, k), 0.0_dp), 1.0_dp), 1, b, st)
                v(k) = v(k) + eps(k) * dot_product(b(0, :), c_dot)
            end do
        end if
        call chord_sums(self, self%w * v, y_dot)
    end subroutine sxr_jvp

    subroutine sxr_vjp(self, sample, eta, y_bar, q_bar, eta_bar)
        class(soft_xray_t), intent(in) :: self
        type(plasma_sample_t), intent(in) :: sample
        real(dp), intent(in) :: eta(:), y_bar(:)
        real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
        real(dp) :: eps(size(self%w)), de(4, size(self%w)), v(size(self%w))
        real(dp), allocatable :: f_bar(:, :), s_bar(:)
        integer :: k

        call emission(self, sample, eta, eps, de)
        do k = 1, size(v)
            v(k) = self%w(k) * self%etendue(self%owner(k)) / four_pi * y_bar(self%owner(k))
        end do
        q_bar = 0.0_dp
        q_bar(Q_S, :) = v * de(1, :)
        q_bar(Q_NE, :) = v * de(2, :)
        q_bar(Q_TE, :) = v * de(3, :)
        q_bar(Q_ZEFF, :) = v * de(4, :)
        eta_bar = 0.0_dp
        if (self%n_nuisance() == 0) return
        allocate(f_bar(N_PROFILES, size(v)), s_bar(size(v)))
        f_bar = 0.0_dp
        where (sample%inside) f_bar(P_NE, :) = v
        call self%emissivity%vjp(eta, sample%q(Q_S, :), f_bar, eta_bar, s_bar)
    end subroutine sxr_vjp
end module tiago_soft_xray
