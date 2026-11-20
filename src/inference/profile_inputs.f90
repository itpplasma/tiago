module tiago_profile_inputs
    !! Profile-diagnostic inputs of the command-line fits (docs/joint-fit.md):
    !! prior means of the kinetic profiles from a table, point probes grouped
    !! by instrument, and interferometer chords, as slots of
    !! tiago_profile_posterior.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_kinetic_profiles, only: kinetic_profiles_t, P_NE, P_TE, P_TI, P_ZEFF
    use tiago_profile_diagnostic, only: Q_NE, Q_TE, Q_TI
    use tiago_profile_posterior, only: profile_slot_t
    use tiago_point_probe, only: point_probe_t
    use tiago_line_integrals, only: interferometer_t
    use tiago_sight_geometry, only: chord_t
    implicit none
    private

    public :: set_profile_priors, read_probe_slots, read_chord_slot

contains

    subroutine set_profile_priors(prof, path, sigma_f, l_core, l_edge, s0, sw)
        !! Prior means of n_e, T_e, T_i (latent) and Z_eff (fixed) at the
        !! Greville abscissae, linear in s between the rows s,ne,Te,Ti,Zeff of path.
        type(kinetic_profiles_t), intent(inout) :: prof
        character(*), intent(in) :: path
        real(dp), intent(in) :: sigma_f(3), l_core, l_edge, s0, sw
        real(dp), allocatable :: t(:, :), m(:, :)
        integer :: u, n, j, c
        n = count_rows(path)
        if (n < 2) error stop 'prior_table: at least two rows'
        allocate(t(5, n), m(4, prof%ncoef))
        open(newunit=u, file=path, status='old', action='read')
        read(u, *)
        do j = 1, n
            read(u, *) t(:, j)
        end do
        close(u)
        do c = 1, 4
            do j = 1, prof%ncoef
                m(c, j) = interp(t(1, :), t(c + 1, :), prof%greville(j))
            end do
        end do
        call prof%set_profile(P_NE, m(1, :), .true., sigma_f(1), l_core, l_edge, s0, sw)
        call prof%set_profile(P_TE, m(2, :), .true., sigma_f(2), l_core, l_edge, s0, sw)
        call prof%set_profile(P_TI, m(3, :), .true., sigma_f(3), l_core, l_edge, s0, sw)
        call prof%set_profile(P_ZEFF, m(4, :), .false.)
    end subroutine set_profile_priors

    subroutine read_probe_slots(path, align_sigma, cal_sigma, model_rel, slots)
        !! Rows group,quantity,R,phi,Z,value,sigma with quantity 1 n_e [m^-3],
        !! 2 T_e [eV], 3 T_i [eV]; one point-probe slot per group 1..G, with a
        !! rigid (dR, dZ) alignment prior of sd align_sigma(g) [m] (none if 0),
        !! a log calibration prior cal_sigma(g) and relative model error model_rel(g).
        character(*), intent(in) :: path
        real(dp), intent(in) :: align_sigma(:), cal_sigma(:), model_rel(:)
        type(profile_slot_t), allocatable, intent(inout) :: slots(:)
        integer, parameter :: qmap(3) = [Q_NE, Q_TE, Q_TI]
        integer, allocatable :: group(:), quantity(:), idx(:)
        real(dp), allocatable :: x(:, :), val(:), sig(:)
        type(profile_slot_t), allocatable :: add(:)
        type(point_probe_t) :: p
        integer :: u, n, j, g
        n = count_rows(path)
        if (n < 1) error stop 'probes_file: no rows'
        allocate(group(n), quantity(n), x(3, n), val(n), sig(n))
        open(newunit=u, file=path, status='old', action='read')
        read(u, *)
        do j = 1, n
            read(u, *) group(j), quantity(j), x(:, j), val(j), sig(j)
        end do
        close(u)
        if (any(quantity < 1 .or. quantity > 3)) error stop 'probes_file: quantity must be 1, 2 or 3'
        if (any(sig <= 0)) error stop 'probes_file: sigma must be positive'
        if (minval(group) < 1 .or. maxval(group) > size(align_sigma)) &
            error stop 'probes_file: groups 1..16'
        allocate(add(maxval(group)))
        do g = 1, size(add)
            idx = pack([(j, j = 1, n)], group == g)
            if (size(idx) == 0) error stop 'probes_file: empty group'
            p%x0 = x(:, idx)
            p%quantity = qmap(quantity(idx))
            p%aligned = align_sigma(g) <= 0
            allocate(add(g)%d, source=p)
            if (align_sigma(g) > 0) then
                add(g)%eta_sigma = [align_sigma(g), align_sigma(g)]
            else
                allocate(add(g)%eta_sigma(0))
            end if
            add(g)%cal_sigma = cal_sigma(g)
            add(g)%model_rel = model_rel(g)
            add(g)%data = val(idx)
            add(g)%sigma = sig(idx)
        end do
        if (.not. allocated(slots)) allocate(slots(0))
        slots = [slots, add]
    end subroutine read_probe_slots

    subroutine read_chord_slot(path, cal_sigma, model_rel, slots)
        !! Rows x0,y0,z0,dx,dy,dz,length [m],value [m^-2],sigma [m^-2]: one
        !! interferometer slot with a common log calibration prior.
        character(*), intent(in) :: path
        real(dp), intent(in) :: cal_sigma, model_rel
        type(profile_slot_t), allocatable, intent(inout) :: slots(:)
        type(interferometer_t) :: interf
        type(chord_t), allocatable :: chords(:)
        type(profile_slot_t) :: add(1)
        real(dp) :: row(9)
        real(dp), allocatable :: v(:), s(:)
        character(len=32) :: tag
        integer :: u, j, n
        n = count_rows(path)
        if (n < 1) error stop 'chords_file: no rows'
        allocate(chords(n), v(n), s(n))
        open(newunit=u, file=path, status='old', action='read')
        read(u, *)
        do j = 1, n
            read(u, *) row
            write(tag, '("chord_",I0)') j
            chords(j)%label = trim(tag)
            chords(j)%origin = row(1:3)
            chords(j)%direction = row(4:6) / norm2(row(4:6))
            chords(j)%length = row(7)
            chords(j)%npanel = 80
            chords(j)%ngauss = 4
            v(j) = row(8)
            s(j) = row(9)
        end do
        close(u)
        call interf%init(chords)
        allocate(add(1)%d, source=interf)
        allocate(add(1)%eta_sigma(0))
        add(1)%cal_sigma = cal_sigma
        add(1)%model_rel = model_rel
        add(1)%data = v
        add(1)%sigma = s
        if (.not. allocated(slots)) allocate(slots(0))
        slots = [slots, add]
    end subroutine read_chord_slot

    integer function count_rows(path)
        !! Non-empty data rows after a header line.
        character(*), intent(in) :: path
        character(len=1024) :: line
        integer :: u, io
        count_rows = 0
        open(newunit=u, file=path, status='old', action='read', iostat=io)
        if (io /= 0) error stop 'cannot open '//path
        read(u, '(A)', iostat=io) line
        do
            read(u, '(A)', iostat=io) line
            if (io /= 0) exit
            if (len_trim(line) > 0) count_rows = count_rows + 1
        end do
        close(u)
    end function count_rows

    pure real(dp) function interp(x, y, x0)
        real(dp), intent(in) :: x(:), y(:), x0
        integer :: i
        if (x0 <= x(1)) then
            interp = y(1)
            return
        end if
        do i = 2, size(x)
            if (x0 <= x(i)) then
                interp = y(i - 1) + (y(i) - y(i - 1)) * (x0 - x(i - 1)) / (x(i) - x(i - 1))
                return
            end if
        end do
        interp = y(size(y))
    end function interp

end module tiago_profile_inputs
