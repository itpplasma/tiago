module tiago_dispersion_interferometer
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use tiago_profile_diagnostic, only: plasma_sample_t
   use tiago_recorded_diagnostic, only: recorded_diagnostic_t
   use tiago_calibration, only: calibration_prior_t
   use tiago_line_integrals, only: interferometer_t
   use tiago_sight_geometry, only: chord_t
   use tiago_dispersion_phase, only: dispersion_phase
   use tiago_dispersion_iq, only: dispersion_iq
   implicit none
   private
   integer, parameter, public :: DI_NCAL = 8
   integer, parameter, public :: DI_WAVELENGTH = 1, DI_AI = 2, DI_AQ = 3, &
                         DI_OI = 4, DI_OQ = 5, DI_PHASE0 = 6, DI_NEUTRAL = 7, DI_ENV = 8

   type, extends(recorded_diagnostic_t), public :: dispersion_interferometer_t
      type(interferometer_t) :: integral
      integer :: nchord = 0
      real(dp), allocatable :: neutral(:), environment(:)
   contains
      procedure :: init => dispersion_init
      procedure :: n_signals => dispersion_count
      procedure :: n_nuisance => dispersion_nuisance
      procedure :: label => dispersion_label
      procedure :: sample_points => dispersion_points
      procedure :: predict => dispersion_predict
      procedure :: jvp => dispersion_jvp
      procedure :: vjp => dispersion_vjp
      procedure :: quadrature_check
   end type
contains
   subroutine check_calibration(c)
      real(dp), intent(in) :: c(DI_NCAL)
      if (.not. all(ieee_is_finite(c))) error stop 'nonfinite DI calibration'
      if (c(DI_WAVELENGTH) <= 0.0_dp) error stop 'nonpositive DI wavelength'
      if (any(c(DI_AI:DI_AQ) <= 0.0_dp)) error stop 'nonpositive DI amplitude'
   end subroutine

   subroutine dispersion_init(self, chords, prior, neutral, environment, status)
      class(dispersion_interferometer_t), intent(inout) :: self
      type(chord_t), intent(in) :: chords(:)
      type(calibration_prior_t), intent(in) :: prior
      real(dp), intent(in) :: neutral(:), environment(:)
      integer, intent(out) :: status
      type(chord_t), allocatable :: rules(:)
      integer :: i, n
      status = 1
      n = size(chords)
      if (n < 1 .or. prior%n_params() /= DI_NCAL) return
      if (.not. allocated(prior%factor)) return
      if (size(neutral) /= n .or. size(environment) /= n) return
      if (.not. all(ieee_is_finite(neutral))) return
      if (.not. all(ieee_is_finite(environment))) return
      do i = 1, n
         if (.not. allocated(chords(i)%label)) return
         if (.not. all(ieee_is_finite(chords(i)%origin))) return
         if (.not. all(ieee_is_finite(chords(i)%direction))) return
         if (.not. ieee_is_finite(chords(i)%length)) return
         if (norm2(chords(i)%direction) <= tiny(1.0_dp)) return
         if (chords(i)%length <= 0.0_dp) return
         if (chords(i)%npanel < 1 .or. chords(i)%ngauss < 1) return
      end do
      if (prior%mean(DI_WAVELENGTH) <= 0.0_dp) return
      if (any(prior%mean(DI_AI:DI_AQ) <= 0.0_dp)) return
      ! Retain both rules so callers can enforce a convergence tolerance.
      rules = [chords, chords]
      do i = 1, n
         rules(i)%npanel = 2*chords(i)%npanel
      end do
      call self%integral%init(rules)
      self%nchord = n
      self%calibration = prior
      self%neutral = neutral
      self%environment = environment
      status = 0
   end subroutine

   integer function dispersion_count(self)
      class(dispersion_interferometer_t), intent(in) :: self
      dispersion_count = 2*self%nchord
   end function

   integer function dispersion_nuisance(self)
      class(dispersion_interferometer_t), intent(in) :: self
      dispersion_nuisance = DI_NCAL
   end function

   function dispersion_label(self, i) result(text)
      class(dispersion_interferometer_t), intent(in) :: self
      integer, intent(in) :: i
      character(:), allocatable :: text
      integer :: c
      c = (i + 1)/2
      if (mod(i, 2) == 1) then
         text = 'I:'//self%integral%chords(c)%label
      else
         text = 'Q:'//self%integral%chords(c)%label
      end if
   end function

   subroutine dispersion_points(self, eta, x, dx_deta)
      class(dispersion_interferometer_t), intent(in) :: self
      real(dp), intent(in) :: eta(:)
      real(dp), allocatable, intent(out) :: x(:, :), dx_deta(:, :, :)
      real(dp), allocatable :: unused(:, :, :)
      call self%integral%sample_points([real(dp) ::], x, unused)
      allocate (dx_deta(3, size(x, 2), DI_NCAL))
      dx_deta = 0.0_dp
   end subroutine

   subroutine response(self, column, c, k, value, dn, dc)
      class(dispersion_interferometer_t), intent(in) :: self
      real(dp), intent(in) :: column, c(DI_NCAL)
      integer, intent(in) :: k
      real(dp), intent(out) :: value(2), dn(2), dc(2, DI_NCAL)
      real(dp) :: phase(1), jp(1, 8), jq(2, 5)
      call dispersion_phase(column, c(DI_WAVELENGTH), self%neutral(k), &
                          c(DI_NEUTRAL), self%environment(k), c(DI_ENV), c(DI_PHASE0), &
                            0.0_dp, phase, jp)
      call dispersion_iq(phase(1), c(DI_AI), c(DI_AQ), c(DI_OI), c(DI_OQ), &
                         value, jq)
      dn = jq(:, 1)*jp(1, 1)
      dc(:, DI_WAVELENGTH) = jq(:, 1)*jp(1, 2)
      dc(:, DI_AI:DI_OQ) = jq(:, 2:5)
      dc(:, DI_PHASE0) = jq(:, 1)*jp(1, 7)
      dc(:, DI_NEUTRAL) = jq(:, 1)*jp(1, 4)
      dc(:, DI_ENV) = jq(:, 1)*jp(1, 6)
   end subroutine

   subroutine dispersion_predict(self, sample, eta, y)
      class(dispersion_interferometer_t), intent(in) :: self
      type(plasma_sample_t), intent(in) :: sample
      real(dp), intent(in) :: eta(:)
      real(dp), intent(out) :: y(:)
      real(dp) :: c(DI_NCAL), columns(2*self%nchord), dn(2), dc(2, DI_NCAL)
      integer :: k
      c = self%calibration%physical(eta)
      call check_calibration(c)
      call self%integral%predict(sample, [real(dp) ::], columns)
      !$omp parallel do default(shared) private(k, dn, dc) if(self%nchord>128)
      do k = 1, self%nchord
         call response(self, columns(k), c, k, y(2*k - 1:2*k), dn, dc)
      end do
      !$omp end parallel do
   end subroutine

   subroutine dispersion_jvp(self, sample, eta, q_dot, eta_dot, y_dot)
      class(dispersion_interferometer_t), intent(in) :: self
      type(plasma_sample_t), intent(in) :: sample
      real(dp), intent(in) :: eta(:), q_dot(:, :), eta_dot(:)
      real(dp), intent(out) :: y_dot(:)
      real(dp) :: c(DI_NCAL), c_dot(DI_NCAL), columns(2*self%nchord)
      real(dp) :: col_dot(2*self%nchord), value(2), dn(2), dc(2, DI_NCAL)
      integer :: k
      c = self%calibration%physical(eta)
      call check_calibration(c)
      c_dot = self%calibration%jvp(eta_dot)
      call self%integral%predict(sample, [real(dp) ::], columns)
      call self%integral%jvp(sample, [real(dp) ::], q_dot, [real(dp) ::], col_dot)
      !$omp parallel do default(shared) private(k, value, dn, dc) &
      !$omp if(self%nchord>128)
      do k = 1, self%nchord
         call response(self, columns(k), c, k, value, dn, dc)
         y_dot(2*k - 1:2*k) = dn*col_dot(k) + matmul(dc, c_dot)
      end do
      !$omp end parallel do
   end subroutine

   subroutine dispersion_vjp(self, sample, eta, y_bar, q_bar, eta_bar)
      class(dispersion_interferometer_t), intent(in) :: self
      type(plasma_sample_t), intent(in) :: sample
      real(dp), intent(in) :: eta(:), y_bar(:)
      real(dp), intent(out) :: q_bar(:, :), eta_bar(:)
      real(dp) :: c(DI_NCAL), c_bar(DI_NCAL), columns(2*self%nchord)
      real(dp) :: col_bar(2*self%nchord), value(2), dn(2), dc(2, DI_NCAL)
      real(dp) :: empty(0)
      integer :: k
      c = self%calibration%physical(eta)
      call check_calibration(c)
      call self%integral%predict(sample, [real(dp) ::], columns)
      c_bar = 0.0_dp
      col_bar = 0.0_dp
      !$omp parallel do default(shared) private(k, value, dn, dc) &
      !$omp reduction(+:c_bar) if(self%nchord>128)
      do k = 1, self%nchord
         call response(self, columns(k), c, k, value, dn, dc)
         col_bar(k) = dot_product(dn, y_bar(2*k - 1:2*k))
         c_bar = c_bar + matmul(transpose(dc), y_bar(2*k - 1:2*k))
      end do
      !$omp end parallel do
      call self%integral%vjp(sample, [real(dp) ::], col_bar, q_bar, empty)
      eta_bar = self%calibration%vjp(c_bar)
   end subroutine

   subroutine quadrature_check(self, sample, atol, rtol, columns, error, status)
      class(dispersion_interferometer_t), intent(in) :: self
      type(plasma_sample_t), intent(in) :: sample
      real(dp), intent(in) :: atol, rtol
      real(dp), intent(out) :: columns(:), error(:)
      integer, intent(out) :: status
      real(dp) :: both(2*self%nchord)
      status = 1
      if (.not. ieee_is_finite(atol) .or. .not. ieee_is_finite(rtol)) return
      if (atol < 0.0_dp .or. rtol < 0.0_dp) return
      call self%integral%predict(sample, [real(dp) ::], both)
      columns = both(:self%nchord)
      error = abs(columns - both(self%nchord + 1:))
      if (.not. all(ieee_is_finite(both))) return
      if (all(error <= atol + rtol*abs(columns))) status = 0
   end subroutine
end module
