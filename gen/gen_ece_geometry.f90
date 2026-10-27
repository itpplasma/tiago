program gen_ece_geometry
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
                            operator(*), operator(/), operator(**), sin, cos, exp
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs_many
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    implicit none
    type(arena_t), target :: arena
    type(kernel_spec_t) :: spec
    type(expr_t) :: roots(100), w(4), t(3), zero(3), angular, f, g, h(4)
    type(expr_t) :: xx, ff(2), gx(3), hx(3, 3), du(3), jp(3, 3), px, gp(3)
    integer :: i, j, k, l, n, mode, nr
    character(20), parameter :: modes(4) = [character(20) :: &
                                            'vmec_r_x', 'vmec_r_y', 'vmec_z', &
                                            'vmec_cos']

    call arena%init()
    do i = 1, 3
        t(i) = sym(arena, 't'//digit(i))
        zero(i) = num(arena, 0)
    end do
    do i = 1, 4
        w(i) = at('w', i)
    end do
    do mode = 1, 4
        call reset()
        spec%args = [str('w'), str('m'), str('n'), str('theta'), str('phi')]
        spec%arg_shapes = [str('(4)'), str(''), str(''), str(''), str('')]
        angular = cos(sym(arena, 'm')*sym(arena, 'theta') - &
                      sym(arena, 'n')*sym(arena, 'phi'))
        if (mode == 3) angular = sin(sym(arena, 'm')* &
                                     sym(arena, 'theta') - sym(arena, 'n')*sym(arena, &
                                                                               'phi'))
        if (mode == 1) angular = angular*cos(sym(arena, 'phi'))
        if (mode == 2) angular = angular*sin(sym(arena, 'phi'))
        n = 0
        call mode_root([0, 0, 0], 'value')
        do i = 1, 3
            call mode_root([i, 0, 0], index('d', i))
            do j = 1, 3
                call mode_root([i, j, 0], index('d2', i, j))
                do k = 1, 3
                    call mode_root([i, j, k], index('d3', i, j, k))
                end do
            end do
        end do
        spec%outputs = [str('value'), str('d'), str('d2'), str('d3')]
        spec%output_shapes = [str(''), str('(3)'), str('(3,3)'), str('(3,3,3)')]
        call emit(trim(modes(mode)))
    end do

    call reset()
    spec%args = [str('t')]
    f = sym(arena, 't')
    h = [2*f**3 - 3*f**2 + 1, f**3 - 2*f**2 + f, -2*f**3 + 3*f**2, f**3 - f**2]
    n = 0
    do i = 1, 4
        call root(diff(diff(h(i), f), f), index('d2', i))
        call root(diff(diff(diff(h(i), f), f), f), index('d3', i))
    end do
    spec%outputs = [str('d2'), str('d3')]
    spec%output_shapes = [str('(4)'), str('(4)')]
    call emit('hermite_higher')

    ! B = B^theta X_theta + B^phi X_phi, in the VMEC coordinate basis.
    call reset()
    spec%args = [str('xu'), str('xuu'), str('xuuu'), str('v'), str('vu'), str('vuu')]
    spec%arg_shapes = [str('(3)'), str('(3,3)'), str('(3,3,3)'), str('(2)'), &
                       str('(3,2)'), str('(3,3,2)')]
    xx = num(arena, 0)
    do i = 1, 3
        xx = xx + at('xu', i)*t(i)
        do j = 1, 3
            xx = xx + at('xuu', i, j)*t(i)*t(j)/2
            do k = 1, 3
                xx = xx + at('xuuu', i, j, k)*t(i)*t(j)*t(k)/6
            end do
        end do
    end do
    do l = 1, 2
        ff(l) = at('v', l)
        do i = 1, 3
            ff(l) = ff(l) + at('vu', i, l)*t(i)
            do j = 1, 3
                ff(l) = ff(l) + at('vuu', i, j, l)*t(i)*t(j)/2
            end do
        end do
    end do
    f = ff(1)*diff(xx, t(2)) + ff(2)*diff(xx, t(3))
    n = 0
    call root(subs_many(f, t, zero), 'value')
    do i = 1, 3
        call root(subs_many(diff(f, t(i)), t, zero), index('d', i))
        do j = 1, 3
            call root(subs_many(diff(diff(f, t(i)), t(j)), t, zero), index('d2', i, j))
        end do
    end do
    spec%outputs = [str('value'), str('d'), str('d2')]
    spec%output_shapes = [str(''), str('(3)'), str('(3,3)')]
    call emit('vmec_vector_component')

    ! Implicit inverse-map jet. A=(X_u)^-1; all parameter rows hold x fixed.
    call reset()
    spec%args = [str('a'), str('x2'), str('gu'), str('hu'), str('xp'), &
                 str('xup'), str('fp'), str('gup')]
    spec%arg_shapes = [str('(3,3)'), str('(3,3,3)'), str('(3)'), str('(3,3)'), &
                       str('(3)'), str('(3,3)'), str(''), str('(3)')]
    do i = 1, 3
        gx(i) = num(arena, 0)
        du(i) = num(arena, 0)
        do j = 1, 3
            gx(i) = gx(i) + at('gu', j)*at('a', j, i)
            du(i) = du(i) - at('a', i, j)*at('xp', j)
        end do
    end do
    do i = 1, 3
        do j = 1, 3
            hx(i, j) = num(arena, 0)
            do k = 1, 3
                do l = 1, 3
                    g = at('hu', k, l)
                    do nr = 1, 3
                        g = g - gx(nr)*at('x2', nr, k, l)
                    end do
                    hx(i, j) = hx(i, j) + at('a', k, i)*g*at('a', l, j)
                end do
            end do
            jp(i, j) = at('xup', i, j)
            do k = 1, 3
                jp(i, j) = jp(i, j) + at('x2', i, j, k)*du(k)
            end do
        end do
    end do
    px = sym(arena, 'fp')
    do i = 1, 3
        px = px + at('gu', i)*du(i)
    end do
    do i = 1, 3
        gp(i) = at('gup', i)
        do j = 1, 3
            gp(i) = gp(i) + at('hu', i, j)*du(j) - jp(j, i)*gx(j)
        end do
    end do
    n = 0
    call root(px, 'px')
    do i = 1, 3
        call root(gx(i), index('gx', i))
        g = num(arena, 0)
        do j = 1, 3
            g = g + at('a', j, i)*gp(j)
            call root(hx(i, j), index('hx', i, j))
        end do
        call root(g, index('gpx', i))
    end do
    spec%outputs = [str('gx'), str('hx'), str('px'), str('gpx')]
    spec%output_shapes = [str('(3)'), str('(3,3)'), str(''), str('(3)')]
    call emit('medium_pullback')

    call reset()
    spec%args = [str('rho'), str('u'), str('us'), str('uss'), str('unit'), &
                 str('up'), str('usp')]
    f = sym(arena, 'unit')*exp(sym(arena, 'u'))
    g = 2*sym(arena, 'rho')*sym(arena, 'us')
    n = 0
    call root(f, 'value')
    call root(f*g, 'dr')
    call root(f*(g*g + 2*sym(arena, 'us') + &
                 4*sym(arena, 'rho')**2*sym(arena, 'uss')), 'drr')
    call root(diff(f, sym(arena, 'u'))*sym(arena, 'up'), 'p')
    call root(f*g*sym(arena, 'up') + f*2*sym(arena, 'rho')*sym(arena, 'usp'), 'drp')
    spec%outputs = [str('value'), str('dr'), str('drr'), str('p'), str('drp')]
    call emit('density_label_jet')
    call reset()
    spec%args = [str('rcc'), str('rss'), str('zsc'), str('zcs')]
    n = 0
    call root((sym(arena, 'rcc') + sym(arena, 'rss'))/2, index('rc', 1))
    call root((sym(arena, 'rcc') - sym(arena, 'rss'))/2, index('rc', 2))
    call root((sym(arena, 'zsc') - sym(arena, 'zcs'))/2, index('zs', 1))
    call root((sym(arena, 'zsc') + sym(arena, 'zcs'))/2, index('zs', 2))
    do j = 1, 4
        do i = 1, 2
            call root(diff(roots(i), sym(arena, chars(spec%args(j)))), index('rr', i, &
                j))
            call root(diff(roots(i + 2), sym(arena, chars(spec%args(j)))), index('zr', &
                i, j))
        end do
    end do
    spec%outputs = [str('rc'), str('zs'), str('rr'), str('zr')]
    spec%output_shapes = [str('(2)'), str('(2)'), str('(2,4)'), str('(2,4)')]
    call emit('vmec_product_pair')

contains
    function digit(i) result(s)
        integer, intent(in)::i
        character(:), allocatable::s
        character(16)::b
        write (b, '(I0)') i
        s = trim(b)
    end function
    function index(name, i, j, k) result(s)
        character(*), intent(in)::name
        integer, intent(in)::i
        integer, intent(in), optional::j, k
        character(:), allocatable::s
        s = name//'('//digit(i)
        if (present(j)) s = s//','//digit(j)
        if (present(k)) s = s//','//digit(k)
        s = s//')'
    end function
    function at(name, i, j, k) result(e)
        character(*), intent(in)::name
        integer, intent(in)::i
        integer, intent(in), optional::j, k
        type(expr_t)::e
        e = sym(arena, index(name, i, j, k))
    end function
    subroutine reset()
        type(kernel_spec_t)::fresh
        spec = fresh
        allocate (spec%output_references(100))
    end subroutine
    subroutine root(e, name)
        type(expr_t), intent(in)::e
        character(*), intent(in)::name
        n = n + 1
        roots(n) = e
        spec%output_references(n) = str(name)
    end subroutine
    subroutine mode_root(ids, name)
        integer, intent(in)::ids(3)
        character(*), intent(in)::name
        type(expr_t)::e
        integer::j, radial
        e = angular
        radial = 0
        do j = 1, 3
            select case (ids(j))
            case (1)
                radial = radial + 1
            case (2)
                e = diff(e, sym(arena, 'theta'))
            case (3)
                e = diff(e, sym(arena, 'phi'))
            end select
        end do
        call root(w(radial + 1)*e, name)
    end subroutine
    subroutine emit(name)
        character(*), intent(in)::name
        character(:), allocatable::code
        integer::unit, last
        spec%output_references = spec%output_references(:n)
        spec%name = str(name)
        spec%module_name = str('tiago_'//name)
        spec%generator = str('gen/gen_ece_geometry.f90')
        spec%generator_revision = str('fortsym 21d157b394e775a999a8d9ed0be409eac1eb9574')
        spec%regenerate_command = str('build-codegen/gen_ece_geometry (see gen/README.md)')
        spec%pure_procedure = .true.
        code = chars(emit_kernel(roots(:n), spec))
        last = len_trim(code)
        do while (last > 0)
            if (code(last:last) /= new_line('a')) exit
            last = last - 1
        end do
        open (newunit=unit, file='src/generated/'//name//'.f90', status='replace')
        write (unit, '(A)') code(:last)
        close (unit)
    end subroutine
end program
