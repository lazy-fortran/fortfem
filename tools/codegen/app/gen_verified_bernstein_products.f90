program gen_verified_bernstein_products
    !! Interval leaves for the verified Bernstein energy audit on affine
    !! triangles. A degree-n polynomial is sum_a c_a B^n_a(l1, l2, l3) with
    !! B^n_a = n!/(a1! a2! a3!) l1^a1 l2^a2 l3^a3 and barycentric coordinates
    !! l1 = 1 - xi - eta, l2 = xi, l3 = eta. Multi-indices are ordered by
    !! descending a1, then descending a2 (the runtime uses the same order).
    !! Every product, elevation and derivative identity is proved natively as
    !! a homogeneous polynomial identity before emission. The weighted square
    !! integrals use the exact Bernstein moments
    !! int_T^ B^n_a B^n_b B^1_g = C(n,a) C(n,b) / (C(2n+1,a+b+g) (2n+2)(2n+3)).
    use, intrinsic :: iso_fortran_env, only: int64
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, rat, operator(+), operator(-), &
                            operator(*), operator(/), operator(**), sqrt
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_engine, only: engine_result_t, VERDICT_TRUE
    use fortsym_rigorous_emit, only: rigorous_kernel_spec_t, interval_runtime, &
                                     emit_rigorous_kernel
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: generated_path, fortsym_revision
    implicit none
    integer, parameter :: max_product = 5, max_gradient = 4, max_square = 6
    integer, parameter :: max_power = 4
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(rigorous_kernel_spec_t) :: spec
    type(expr_t) :: l(3), xi, eta
    character(:), allocatable :: code, message
    integer :: unit, n, m
    logical :: ok

    call arena%init()
    engine = make_native_engine(arena)
    l(1) = sym(arena, 'l1'); l(2) = sym(arena, 'l2'); l(3) = sym(arena, 'l3')
    xi = sym(arena, 'xi'); eta = sym(arena, 'eta')
    spec%generator = str('gen_verified_bernstein_products; '//fortsym_revision())
    spec%runtime = interval_runtime('fortnum_interval', 'interval_t')
    spec%temp_prefix = str('z')
    open (newunit=unit, status='replace', &
          file=generated_path('fortfem_verified_bernstein_products.f90'))
    write (unit, '(a)') 'module fortfem_verified_bernstein_products'
    write (unit, '(a)') 'use fortnum_interval, only: interval_t'
    write (unit, '(a)') 'implicit none'
    write (unit, '(a)') 'contains'
    call geometry_kernels()
    do n = 0, max_product
        call affine_product(n)
        call elevation(n)
    end do
    do n = 1, max_gradient
        call gradient(n)
    end do
    do n = 0, max_square
        call square_integral(n)
    end do
    do m = 1, max_power
        call chord(m)
    end do
    call lagrange_p2()
    call lagrange_p3()
    call source_remainder()
    call dispatchers()
    write (unit, '(a)') 'end module fortfem_verified_bernstein_products'
    close (unit)

contains

    integer function count_of(degree)
        integer, intent(in) :: degree
        count_of = (degree + 1)*(degree + 2)/2
    end function count_of

    function indices(degree) result(a)
        integer, intent(in) :: degree
        integer, allocatable :: a(:, :)
        integer :: i, j, k
        allocate (a(3, count_of(degree)))
        k = 0
        do i = degree, 0, -1
            do j = degree - i, 0, -1
                k = k + 1
                a(:, k) = [i, j, degree - i - j]
            end do
        end do
    end function indices

    integer function position(degree, a)
        integer, intent(in) :: degree, a(3)
        integer, allocatable :: all(:, :)
        integer :: k
        all = indices(degree)
        position = 0
        do k = 1, size(all, 2)
            if (all(1, k) == a(1) .and. all(2, k) == a(2)) position = k
        end do
        if (position == 0 .or. sum(a) /= degree) error stop 'Bernstein index'
    end function position

    integer(int64) function factorial(k)
        integer, intent(in) :: k
        integer :: i
        factorial = 1
        do i = 2, k
            factorial = factorial*i
        end do
    end function factorial

    integer(int64) function multinomial(degree, a)
        integer, intent(in) :: degree, a(3)
        multinomial = factorial(degree)/(factorial(a(1))*factorial(a(2)) &
                                         *factorial(a(3)))
    end function multinomial

    function bernstein(degree, a) result(b)
        integer, intent(in) :: degree, a(3)
        type(expr_t) :: b
        integer :: i
        b = num(arena, multinomial(degree, a))
        do i = 1, 3
            if (a(i) > 0) b = b*l(i)**a(i)
        end do
    end function bernstein

    function symbols(prefix, count) result(s)
        character(*), intent(in) :: prefix
        integer, intent(in) :: count
        type(expr_t), allocatable :: s(:)
        integer :: k
        allocate (s(count))
        do k = 1, count
            s(k) = sym(arena, prefix//itoa(k))
        end do
    end function symbols

    function polynomial(degree, c) result(p)
        integer, intent(in) :: degree
        type(expr_t), intent(in) :: c(:)
        type(expr_t) :: p
        integer, allocatable :: a(:, :)
        integer :: k
        a = indices(degree)
        p = num(arena, 0)
        do k = 1, size(a, 2)
            p = p + c(k)*bernstein(degree, a(:, k))
        end do
    end function polynomial

    subroutine prove(identity, what)
        type(expr_t), intent(in) :: identity
        character(*), intent(in) :: what
        type(engine_result_t) :: result
        result = engine%zero_test(identity)
        if (result%verdict /= VERDICT_TRUE) then
            print '(a)', 'identity not proved: '//what
            error stop 'Bernstein identity'
        end if
    end subroutine prove

    subroutine set_names(name, inputs, outputs)
        character(*), intent(in) :: name
        character(*), intent(in) :: inputs(:), outputs(:)
        integer :: k
        spec%name = str(name)
        if (allocated(spec%args)) deallocate (spec%args)
        if (allocated(spec%outputs)) deallocate (spec%outputs)
        allocate (spec%args(size(inputs)), spec%outputs(size(outputs)))
        do k = 1, size(inputs)
            spec%args(k) = str(trim(inputs(k)))
        end do
        do k = 1, size(outputs)
            spec%outputs(k) = str(trim(outputs(k)))
        end do
    end subroutine set_names

    function names(prefix, count) result(list)
        character(*), intent(in) :: prefix
        integer, intent(in) :: count
        character(16), allocatable :: list(:)
        integer :: k
        allocate (list(count))
        do k = 1, count
            list(k) = prefix//itoa(k)
        end do
    end function names

    subroutine emit(roots)
        type(expr_t), intent(in) :: roots(:)
        code = chars(emit_rigorous_kernel(roots, spec, ok=ok, message=message))
        if (.not. ok) then
            print *, message
            error stop 'Bernstein interval generation failed'
        end if
        write (unit, '(a)') code
    end subroutine emit

    subroutine geometry_kernels()
        type(expr_t) :: x(3), y(3), det, cx, cy, g, nx, ny, px, py
        type(expr_t) :: g1, g2, n1x, n1y, n2x, n2y, d
        character(16) :: inputs(6)
        x(1) = sym(arena, 'x1'); y(1) = sym(arena, 'y1')
        x(2) = sym(arena, 'x2'); y(2) = sym(arena, 'y2')
        x(3) = sym(arena, 'x3'); y(3) = sym(arena, 'y3')
        ! Physical gradient = G reference gradient, G = J^-T of the affine map
        ! (xi, eta) -> P1 + xi (P2 - P1) + eta (P3 - P1).
        det = (x(2) - x(1))*(y(3) - y(1)) - (x(3) - x(1))*(y(2) - y(1))
        inputs = [character(16) :: 'x1', 'y1', 'x2', 'y2', 'x3', 'y3']
        call set_names('verified_bernstein_geometry', inputs, &
                       [character(16) :: 'determinant', 'g11', 'g12', 'g21', 'g22'])
        call emit([det, (y(3) - y(1))/det, (y(1) - y(2))/det, &
                   (x(1) - x(3))/det, (x(2) - x(1))/det])
        ! Edge normal rotated from the directed edge a -> b.
        call set_names('verified_bernstein_edge_normal', &
                       [character(16) :: 'xa', 'ya', 'xb', 'yb'], &
                       [character(16) :: 'nx', 'ny'])
        call emit([sym(arena, 'yb') - sym(arena, 'ya'), &
                   sym(arena, 'xa') - sym(arena, 'xb')])
        ! Replace the normal component of c by the shared value g.
        cx = sym(arena, 'cx'); cy = sym(arena, 'cy'); g = sym(arena, 'g')
        nx = sym(arena, 'nx'); ny = sym(arena, 'ny')
        d = (g - nx*cx - ny*cy)/(nx*nx + ny*ny)
        px = cx + d*nx; py = cy + d*ny
        call prove(nx*px + ny*py - g, 'edge projection')
        call set_names('verified_bernstein_edge_projection', &
                       [character(16) :: 'cx', 'cy', 'g', 'nx', 'ny'], &
                       [character(16) :: 'px', 'py'])
        call emit([px, py])
        ! Corner vector with prescribed components along two edge normals.
        g1 = sym(arena, 'g1'); g2 = sym(arena, 'g2')
        n1x = sym(arena, 'n1x'); n1y = sym(arena, 'n1y')
        n2x = sym(arena, 'n2x'); n2y = sym(arena, 'n2y')
        d = n1x*n2y - n1y*n2x
        cx = (g1*n2y - g2*n1y)/d
        cy = (n1x*g2 - n2x*g1)/d
        call prove(n1x*cx + n1y*cy - g1, 'corner first normal')
        call prove(n2x*cx + n2y*cy - g2, 'corner second normal')
        call set_names('verified_bernstein_corner', &
                       [character(16) :: 'g1', 'n1x', 'n1y', 'g2', 'n2x', 'n2y'], &
                       [character(16) :: 'cx', 'cy'])
        call emit([cx, cy])
    end subroutine geometry_kernels

    subroutine affine_product(degree)
        !! d = c * w with c of the given degree and w = sum_i w_i l_i.
        integer, intent(in) :: degree
        type(expr_t), allocatable :: c(:), w(:), d(:)
        integer, allocatable :: a(:, :)
        integer :: k, i, b(3)
        c = symbols('c', count_of(degree)); w = symbols('w', 3)
        a = indices(degree + 1)
        allocate (d(size(a, 2)))
        do k = 1, size(a, 2)
            d(k) = num(arena, 0)
            do i = 1, 3
                if (a(i, k) == 0) cycle
                b = a(:, k); b(i) = b(i) - 1
                d(k) = d(k) + rat(arena, multinomial(degree, b), &
                                  multinomial(degree + 1, a(:, k))) &
                       *c(position(degree, b))*w(i)
            end do
        end do
        call prove(polynomial(degree, c)*(w(1)*l(1) + w(2)*l(2) + w(3)*l(3)) &
                   - polynomial(degree + 1, d), 'affine product')
        call set_names('verified_bernstein_affine_product_'//itoa(degree), &
                       [names('c', size(c)), names('w', 3)], names('d', size(d)))
        call emit(d)
    end subroutine affine_product

    subroutine elevation(degree)
        !! The same polynomial represented at one degree higher.
        integer, intent(in) :: degree
        type(expr_t), allocatable :: c(:), d(:)
        integer, allocatable :: a(:, :)
        integer :: k, i, b(3)
        c = symbols('c', count_of(degree))
        a = indices(degree + 1)
        allocate (d(size(a, 2)))
        do k = 1, size(a, 2)
            d(k) = num(arena, 0)
            do i = 1, 3
                if (a(i, k) == 0) cycle
                b = a(:, k); b(i) = b(i) - 1
                d(k) = d(k) + rat(arena, int(a(i, k), int64), int(degree + 1, int64)) &
                       *c(position(degree, b))
            end do
        end do
        call prove(polynomial(degree, c)*(l(1) + l(2) + l(3)) &
                   - polynomial(degree + 1, d), 'elevation')
        call set_names('verified_bernstein_elevate_'//itoa(degree), &
                       names('c', size(c)), names('d', size(d)))
        call emit(d)
    end subroutine elevation

    subroutine gradient(degree)
        !! Physical gradient coefficients (degree - 1) on an affine cell.
        integer, intent(in) :: degree
        type(expr_t), allocatable :: c(:), dxi(:), deta(:), roots(:)
        type(expr_t) :: g11, g12, g21, g22, p, q
        integer, allocatable :: a(:, :)
        integer :: k, nb
        c = symbols('c', count_of(degree))
        g11 = sym(arena, 'g11'); g12 = sym(arena, 'g12')
        g21 = sym(arena, 'g21'); g22 = sym(arena, 'g22')
        a = indices(degree - 1)
        nb = size(a, 2)
        allocate (dxi(nb), deta(nb), roots(2*nb))
        do k = 1, nb
            dxi(k) = num(arena, degree)*(c(position(degree, a(:, k) + [0, 1, 0])) &
                                         - c(position(degree, a(:, k) + [1, 0, 0])))
            deta(k) = num(arena, degree)*(c(position(degree, a(:, k) + [0, 0, 1])) &
                                          - c(position(degree, a(:, k) + [1, 0, 0])))
        end do
        p = reference(polynomial(degree, c))
        q = reference(polynomial(degree - 1, dxi))
        call prove(diff(p, xi) - q, 'xi derivative')
        q = reference(polynomial(degree - 1, deta))
        call prove(diff(p, eta) - q, 'eta derivative')
        do k = 1, nb
            roots(k) = g11*dxi(k) + g12*deta(k)
            roots(nb + k) = g21*dxi(k) + g22*deta(k)
        end do
        call set_names('verified_bernstein_gradient_'//itoa(degree), &
                       [names('c', size(c)), &
                        [character(16) :: 'g11', 'g12', 'g21', 'g22']], &
                       [names('dx', nb), names('dy', nb)])
        call emit(roots)
    end subroutine gradient

    function reference(p) result(r)
        type(expr_t), intent(in) :: p
        type(expr_t) :: r
        r = subs(p, l(1), num(arena, 1) - xi - eta)
        r = subs(r, l(2), xi)
        r = subs(r, l(3), eta)
    end function reference

    subroutine square_integral(degree)
        !! int_T^ t^2 w over the reference triangle (area 1/2), with t of the
        !! given degree and w = sum_g w_g l_g affine. With scaled coefficients
        !! tt_a = C(n,a) t_a, the scaled-monomial square p_m = sum_(a+b=m)
        !! tt_a tt_b and its product with w give the degree-(2n+1)
        !! coefficients; each scaled monomial l^d integrates to
        !! 1/(C(2n+1,d) (2n+2)(2n+3)) times its multinomial.
        integer, intent(in) :: degree
        type(expr_t), allocatable :: t(:), w(:), tt(:), p(:), q(:)
        type(expr_t) :: total, direct, inner
        integer, allocatable :: a(:, :), b(:, :), d(:, :)
        integer :: i, j, g, k, e(3)
        integer(int64) :: numerator, denominator
        t = symbols('t', count_of(degree)); w = symbols('w', 3)
        a = indices(degree); b = indices(2*degree); d = indices(2*degree + 1)
        allocate (tt(size(a, 2)), p(size(b, 2)), q(size(d, 2)))
        do i = 1, size(a, 2)
            tt(i) = num(arena, multinomial(degree, a(:, i)))*t(i)
        end do
        p = num(arena, 0)
        do i = 1, size(a, 2)
            do j = i, size(a, 2)
                k = position(2*degree, a(:, i) + a(:, j))
                if (j == i) then
                    p(k) = p(k) + tt(i)*tt(j)
                else
                    p(k) = p(k) + 2*tt(i)*tt(j)
                end if
            end do
        end do
        total = num(arena, 0)
        do k = 1, size(d, 2)
            q(k) = num(arena, 0)
            do g = 1, 3
                e = 0; e(g) = 1
                if (d(g, k) == 0) cycle
                q(k) = q(k) + w(g)*p(position(2*degree, d(:, k) - e))
            end do
            total = total + rat(arena, 1_int64, multinomial(2*degree + 1, d(:, k)) &
                                *(2*degree + 2)*(2*degree + 3))*q(k)
        end do
        ! Independent check against the termwise Bernstein moment formula.
        direct = num(arena, 0)
        do g = 1, 3
            e = 0; e(g) = 1
            inner = num(arena, 0)
            do i = 1, size(a, 2)
                do j = 1, size(a, 2)
                    numerator = multinomial(degree, a(:, i)) &
                                *multinomial(degree, a(:, j))
                    denominator = multinomial(2*degree + 1, a(:, i) + a(:, j) + e) &
                                  *(2*degree + 2)*(2*degree + 3)
                    inner = inner + rat(arena, numerator, denominator)*t(i)*t(j)
                end do
            end do
            direct = direct + w(g)*inner
        end do
        if (degree <= 3) call prove(total - direct, 'square integral')
        call set_names('verified_bernstein_square_integral_'//itoa(degree), &
                       [names('t', size(t)), names('w', 3)], &
                       [character(16) :: 'integral'])
        call emit([total])
    end subroutine square_integral

    subroutine chord(power)
        !! Chord of the convex w^-power over [lo, hi] evaluated at w (hi > lo),
        !! an affine upper bound of w^-power there; and the plain power.
        integer, intent(in) :: power
        type(expr_t) :: lo, hi, w
        lo = sym(arena, 'lo'); hi = sym(arena, 'hi'); w = sym(arena, 'w')
        call set_names('verified_bernstein_chord_'//itoa(power), &
                       [character(16) :: 'lo', 'hi', 'w'], [character(16) :: 'value'])
        call emit([((hi - w)/lo**power + (w - lo)/hi**power)/(hi - lo)])
        call set_names('verified_bernstein_reciprocal_'//itoa(power), &
                       [character(16) :: 'w'], [character(16) :: 'value'])
        call emit([num(arena, 1)/w**power])
    end subroutine chord

    subroutine lagrange_p2()
        !! Bernstein coefficients of the P2 Lagrange interpolant with vertex
        !! values u1..u3 and edge-midpoint values u4 (1-2), u5 (2-3), u6 (3-1).
        type(expr_t) :: u(6), b(6), p, half
        integer :: k
        do k = 1, 6
            u(k) = sym(arena, 'u'//itoa(k))
        end do
        half = rat(arena, 1_int64, 2_int64)
        b(position(2, [2, 0, 0])) = u(1)
        b(position(2, [0, 2, 0])) = u(2)
        b(position(2, [0, 0, 2])) = u(3)
        b(position(2, [1, 1, 0])) = 2*u(4) - half*(u(1) + u(2))
        b(position(2, [0, 1, 1])) = 2*u(5) - half*(u(2) + u(3))
        b(position(2, [1, 0, 1])) = 2*u(6) - half*(u(3) + u(1))
        p = polynomial(2, b)
        call prove(subs(subs(subs(p, l(1), half), l(2), half), l(3), num(arena, 0)) &
                   - u(4), 'P2 midpoint 1-2')
        call prove(subs(subs(subs(p, l(1), num(arena, 0)), l(2), half), l(3), half) &
                   - u(5), 'P2 midpoint 2-3')
        call prove(subs(subs(subs(p, l(1), half), l(2), num(arena, 0)), l(3), half) &
                   - u(6), 'P2 midpoint 3-1')
        call set_names('verified_bernstein_lagrange_p2', names('u', 6), names('b', 6))
        call emit(b)
    end subroutine lagrange_p2

    subroutine lagrange_p3()
        !! Vertices, two nodes per directed edge (1-2, 2-3, 3-1), centroid.
        type(expr_t) :: u(10), b(10), edge(6), p, at
        integer :: j, next, k, a(3), powers(3, 10)
        powers = reshape([3,0,0, 0,3,0, 0,0,3, 2,1,0, 1,2,0, &
                          0,2,1, 0,1,2, 1,0,2, 2,0,1, 1,1,1], [3,10])
        do k = 1, 10
            u(k) = sym(arena, 'u'//itoa(k))
        end do
        do j = 1, 3
            next = mod(j, 3) + 1
            k = 2*j + 2
            edge(2*j - 1) = (-5*u(j) + 18*u(k) - 9*u(k + 1) + 2*u(next))/6
            edge(2*j) = (2*u(j) - 9*u(k) + 18*u(k + 1) - 5*u(next))/6
            a = 0; a(j) = 3
            b(position(3, a)) = u(j)
            a(j) = 2; a(next) = 1
            b(position(3, a)) = edge(2*j - 1)
            a(j) = 1; a(next) = 2
            b(position(3, a)) = edge(2*j)
        end do
        b(position(3, [1,1,1])) = (27*u(10) - u(1) - u(2) - u(3) &
                                   - 3*sum_expr(edge))/6
        p = polynomial(3, b)
        do k = 1, 10
            at = p
            do j = 1, 3
                at = subs(at, l(j), rat(arena, int(powers(j,k), int64), 3_int64))
            end do
            call prove(at - u(k), 'P3 interpolation node '//itoa(k))
        end do
        call set_names('verified_bernstein_lagrange_p3', names('u', 10), names('b', 10))
        call emit(b)
    end subroutine lagrange_p3

    function sum_expr(values) result(total)
        type(expr_t), intent(in) :: values(:)
        type(expr_t) :: total
        integer :: k
        total = num(arena, 0)
        do k = 1, size(values)
            total = total + values(k)
        end do
    end function sum_expr

    subroutine source_remainder()
        !! Minkowski: ||r + delta/w^m||_T <= ||r||_T + rho sqrt(|T|)/wmin^m.
        type(expr_t) :: square, rho, determinant, wmin
        integer :: m
        square = sym(arena, 'square'); rho = sym(arena, 'rho')
        determinant = sym(arena, 'determinant'); wmin = sym(arena, 'wmin')
        do m = 0, 2
            if (m == 0) then
                call set_names('verified_bernstein_source_tail_0', &
                    [character(16) :: 'square', 'rho', 'determinant'], &
                    [character(16) :: 'enclosure'])
            else
                call set_names('verified_bernstein_source_tail_'//itoa(m), &
                    [character(16) :: 'square', 'rho', 'determinant', 'wmin'], &
                    [character(16) :: 'enclosure'])
            end if
            call emit([(sqrt(square) + rho*sqrt(determinant/2)/wmin**m)**2])
        end do
    end subroutine source_remainder

    subroutine dispatchers()
        !! Array front ends over the fixed-degree scalar leaves.
        integer :: n
        write (unit, '(a)') '! Generated by fortsym. Do not edit.'
        write (unit, '(a)') '! Generator: '//chars(spec%generator)
        write (unit, '(a)') 'pure subroutine '// &
            'verified_bernstein_affine_product(n, c, w, d)'
        write (unit, '(a)') '    integer, intent(in) :: n'
        write (unit, '(a)') '    type(interval_t), intent(in) :: c(:), w(3)'
        write (unit, '(a)') '    type(interval_t), intent(out) :: d(:)'
        write (unit, '(a)') '    select case (n)'
        do n = 0, max_product
            write (unit, '(a)') '    case ('//itoa(n)//')'
            write (unit, '(a)') '        call verified_bernstein_affine_product_'// &
                itoa(n)//'('//args('c', count_of(n))//', w(1), w(2), w(3), '// &
                args('d', count_of(n + 1))//')'
        end do
        write (unit, '(a)') '    end select'
        write (unit, '(a)') 'end subroutine verified_bernstein_affine_product'
        write (unit, '(a)') 'pure subroutine verified_bernstein_elevate(n, c, d)'
        write (unit, '(a)') '    integer, intent(in) :: n'
        write (unit, '(a)') '    type(interval_t), intent(in) :: c(:)'
        write (unit, '(a)') '    type(interval_t), intent(out) :: d(:)'
        write (unit, '(a)') '    select case (n)'
        do n = 0, max_product
            write (unit, '(a)') '    case ('//itoa(n)//')'
            write (unit, '(a)') '        call verified_bernstein_elevate_'//itoa(n)// &
                '('//args('c', count_of(n))//', '//args('d', count_of(n + 1))//')'
        end do
        write (unit, '(a)') '    end select'
        write (unit, '(a)') 'end subroutine verified_bernstein_elevate'
        write (unit, '(a)') 'pure subroutine '// &
            'verified_bernstein_gradient(n, c, g, dx, dy)'
        write (unit, '(a)') '    integer, intent(in) :: n'
        write (unit, '(a)') '    type(interval_t), intent(in) :: c(:), g(4)'
        write (unit, '(a)') '    type(interval_t), intent(out) :: dx(:), dy(:)'
        write (unit, '(a)') '    select case (n)'
        do n = 1, max_gradient
            write (unit, '(a)') '    case ('//itoa(n)//')'
            write (unit, '(a)') '        call verified_bernstein_gradient_'//itoa(n)// &
                '('//args('c', count_of(n))//', g(1), g(2), g(3), g(4), '// &
                args('dx', count_of(n - 1))//', '//args('dy', count_of(n - 1))//')'
        end do
        write (unit, '(a)') '    end select'
        write (unit, '(a)') 'end subroutine verified_bernstein_gradient'
        write (unit, '(a)') 'pure subroutine '// &
            'verified_bernstein_square_integral(n, t, w, value)'
        write (unit, '(a)') '    integer, intent(in) :: n'
        write (unit, '(a)') '    type(interval_t), intent(in) :: t(:), w(3)'
        write (unit, '(a)') '    type(interval_t), intent(out) :: value'
        write (unit, '(a)') '    select case (n)'
        do n = 0, max_square
            write (unit, '(a)') '    case ('//itoa(n)//')'
            write (unit, '(a)') '        call verified_bernstein_square_integral_'// &
                itoa(n)//'('//args('t', count_of(n))//', w(1), w(2), w(3), value)'
        end do
        write (unit, '(a)') '    end select'
        write (unit, '(a)') 'end subroutine verified_bernstein_square_integral'
        write (unit, '(a)') 'pure subroutine '// &
            'verified_bernstein_chord(m, lo, hi, w, value)'
        write (unit, '(a)') '    integer, intent(in) :: m'
        write (unit, '(a)') '    type(interval_t), intent(in) :: lo, hi, w'
        write (unit, '(a)') '    type(interval_t), intent(out) :: value'
        write (unit, '(a)') '    select case (m)'
        do n = 1, max_power
            write (unit, '(a)') '    case ('//itoa(n)//')'
            write (unit, '(a)') '        call verified_bernstein_chord_'//itoa(n)// &
                '(lo, hi, w, value)'
        end do
        write (unit, '(a)') '    end select'
        write (unit, '(a)') 'end subroutine verified_bernstein_chord'
        write (unit, '(a)') 'pure subroutine verified_bernstein_reciprocal(m, w, value)'
        write (unit, '(a)') '    integer, intent(in) :: m'
        write (unit, '(a)') '    type(interval_t), intent(in) :: w'
        write (unit, '(a)') '    type(interval_t), intent(out) :: value'
        write (unit, '(a)') '    select case (m)'
        do n = 1, max_power
            write (unit, '(a)') '    case ('//itoa(n)//')'
            write (unit, '(a)') '        call verified_bernstein_reciprocal_' &
                //itoa(n)// &
                '(w, value)'
        end do
        write (unit, '(a)') '    end select'
        write (unit, '(a)') 'end subroutine verified_bernstein_reciprocal'
    end subroutine dispatchers

    function args(prefix, count) result(text)
        character(*), intent(in) :: prefix
        integer, intent(in) :: count
        character(:), allocatable :: text
        integer :: k
        text = ''
        do k = 1, count
            if (k > 1) text = text//', '
            if (mod(k, 6) == 0) text = text//'&'//new_line('a')//'            '
            text = text//prefix//'('//itoa(k)//')'
        end do
    end function args

    function itoa(value) result(text)
        integer, intent(in) :: value
        character(:), allocatable :: text
        character(16) :: buffer
        write (buffer, '(i0)') value
        text = trim(buffer)
    end function itoa
end program gen_verified_bernstein_products
