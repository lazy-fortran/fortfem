program gen_triangle_feec_tables
    !! Exact reference data of the arbitrary-order triangle de Rham complex
    !! P_k -> N_k (first-kind Nedelec) -> RT_{k-1} -> P_{k-1}.
    !!
    !! The degrees of freedom are those of fortfem_triangle_nedelec_arbitrary_order:
    !! edge moments against shifted Legendre polynomials and interior moments
    !! against monomials. FortSym solves the moment system in exact rational
    !! arithmetic. The basis is emitted as Bernstein coefficients (all Bernstein
    !! polynomials are nonnegative, so evaluation is well conditioned) and every
    !! discrete operator table is emitted from exact values. Each emitted number
    !! is rounded to binary64 exactly once. Before emission the generator proves
    !! exactly: duality with the moments, the Bernstein representation, the
    !! gradient table and curl(grad) = 0, div(rot grad) = 0 at the table level.
    use fortsym_arena, only: arena_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
        operator(*), operator(**)
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_linalg, only: exact_linear_system_result_t, &
        solve_exact_linear_system
    use fortsym_parse, only: parse_expr
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    use fortfem_codegen_rational, only: q_t, q, q_frac, q_text, q_zero, &
        q_literal, q_factorial, q_binomial, operator(+), operator(-), &
        operator(*), operator(/), operator(==)
    implicit none

    integer, parameter :: max_order = 8

    type :: vpoly_t
        !! Vector polynomial sum c(a, b) x^a y^b, components 1 and 2.
        type(q_t), allocatable :: c(:, :, :)   ! (2, 0:n, 0:n)
    end type vpoly_t

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    integer :: table_unit, order

    call arena%init()
    engine = make_native_engine(arena)

    open(newunit=table_unit, &
        file=generated_path("fortfem_triangle_feec_tables.f90"), &
        status="replace", action="write")
    call write_table_header()
    call write_bernstein_section()
    call write_operator_section()
    write(table_unit, "(a)") ""
    write(table_unit, "(a)") "end module fortfem_generated_triangle_feec_tables"
    close(table_unit)

    call write_bernstein_jets()

contains

    ! ------------------------------------------------------------------
    ! Exact construction of one order
    ! ------------------------------------------------------------------

    subroutine build_order(k, bern, gref, cref, dref)
        !! bern(2, nb_k, ndof): Bernstein coefficients of the basis.
        !! gref(ndof, nl0): Nedelec coefficients of the P_k Lagrange gradients.
        !! cref(nld, ndof): curls of the basis at the P_{k-1} Lagrange nodes.
        !! dref(nld, ndof): divergences of RT_{k-1} at the same nodes.
        integer, intent(in) :: k
        type(q_t), allocatable, intent(out) :: bern(:, :, :), gref(:, :)
        type(q_t), allocatable, intent(out) :: cref(:, :), dref(:, :)

        type(vpoly_t), allocatable :: candidates(:), basis(:)
        type(q_t), allocatable :: moments(:, :), coefficients(:, :)
        type(q_t), allocatable :: nodes(:, :), lagrange(:, :), curl(:, :)
        type(q_t), allocatable :: grad_x(:, :), grad_y(:, :)
        type(vpoly_t) :: gradient, combination
        integer :: ndof, nb, nl0, nld, j, i, l, half, e3, target, source

        ndof = k*(k + 2)
        nb = (k + 1)*(k + 2)/2
        nl0 = nb
        nld = k*(k + 1)/2

        call nedelec_candidates(k, candidates)
        allocate (moments(ndof, ndof))
        do j = 1, ndof
            moments(:, j) = nedelec_moments(candidates(j), k)
        end do
        call exact_inverse(moments, coefficients)

        allocate (basis(ndof), bern(2, nb, ndof))
        do j = 1, ndof
            basis(j) = vpoly_new(k)
            do i = 1, ndof
                if (q_zero(coefficients(i, j))) cycle
                call vpoly_axpy(coefficients(i, j), candidates(i), basis(j))
            end do
            call check_identity(nedelec_moments(basis(j), k), j, &
                "basis is not dual to the moments")
            bern(1, :, j) = to_bernstein(basis(j)%c(1, :, :), k, k)
            bern(2, :, j) = to_bernstein(basis(j)%c(2, :, :), k, k)
            call check_bernstein(basis(j), bern(:, :, j), k)
        end do

        ! Lagrange gradients expanded in the Nedelec basis (exact dofs).
        call lagrange_monomials(k, lagrange)
        allocate (gref(ndof, nl0))
        do i = 1, nl0
            call scalar_gradient(lagrange(:, i), k, grad_x, grad_y)
            gradient = vpoly_new(k)
            call vpoly_set(gradient, grad_x, grad_y, k)
            gref(:, i) = nedelec_moments(gradient, k)
            combination = vpoly_new(k)
            do j = 1, ndof
                if (q_zero(gref(j, i))) cycle
                call vpoly_axpy(gref(j, i), basis(j), combination)
            end do
            call check_same(combination, gradient, k, &
                "gradient table does not reproduce the Lagrange gradient")
        end do

        ! Curls of the basis at the P_{k-1} nodes and RT_{k-1} divergences.
        call lagrange_nodes(k - 1, nodes)
        allocate (cref(nld, ndof), dref(nld, ndof), curl(0:k, 0:k))
        do j = 1, ndof
            curl = scalar_curl(basis(j), k)
            do l = 1, nld
                cref(l, j) = poly_value(curl, k, nodes(1, l), nodes(2, l))
            end do
        end do
        e3 = 3*k
        half = k*(k - 1)/2
        do j = 1, e3
            dref(:, j) = cref(:, j)
        end do
        do l = 1, half
            target = e3 + l
            source = e3 + half + l
            dref(:, target) = cref(:, source)
            dref(:, source) = negate(cref(:, target))
        end do
        call check_complex(k, gref, cref, dref)
    end subroutine build_order

    subroutine check_complex(k, gref, cref, dref)
        !! curl(grad) = 0 and div(rot grad) = 0 exactly at the table level.
        !! The rotation maps Nedelec dof j to RT dof pi(j) with sign sigma(j).
        integer, intent(in) :: k
        type(q_t), intent(in) :: gref(:, :), cref(:, :), dref(:, :)

        type(q_t) :: s, t
        integer :: l, i, j, half, e3, pj, sj

        e3 = 3*k
        half = k*(k - 1)/2
        do i = 1, size(gref, 2)
            do l = 1, size(cref, 1)
                s = q(0)
                t = q(0)
                do j = 1, size(gref, 1)
                    s = s + cref(l, j)*gref(j, i)
                    call rotation_target(j, e3, half, pj, sj)
                    t = t + sj*(dref(l, pj)*gref(j, i))
                end do
                if (.not. q_zero(s)) error stop "exact curl(grad) /= 0"
                if (.not. q_zero(t)) error stop "exact div(rot grad) /= 0"
            end do
        end do
    end subroutine check_complex

    pure subroutine rotation_target(j, e3, half, pj, sj)
        !! rot90 N_j = sj * RT_pj, the signed permutation of the module RT.
        integer, intent(in) :: j, e3, half
        integer, intent(out) :: pj, sj

        if (j <= e3) then
            pj = j
            sj = 1
        else if (j <= e3 + half) then
            pj = j + half
            sj = -1
        else
            pj = j - half
            sj = 1
        end if
    end subroutine rotation_target

    subroutine exact_inverse(matrix, inverse)
        type(q_t), intent(in) :: matrix(:, :)
        type(q_t), allocatable, intent(out) :: inverse(:, :)

        type(expr_t), allocatable :: a(:, :), identity(:, :)
        type(exact_linear_system_result_t) :: solution
        character(:), allocatable :: message
        logical :: ok
        integer :: i, j, n

        n = size(matrix, 1)
        allocate (a(n, n), identity(n, n), inverse(n, n))
        do j = 1, n
            do i = 1, n
                a(i, j) = parse_expr(arena, chars(matrix(i, j)%value), ok, &
                    message)
                if (.not. ok) error stop "cannot parse exact moment"
                identity(i, j) = num(arena, merge(1, 0, i == j))
            end do
        end do
        solution = solve_exact_linear_system(engine, a, identity)
        if (.not. solution%ok) then
            error stop "exact Nedelec moment system failed"
        end if
        do j = 1, n
            do i = 1, n
                inverse(i, j) = q_text(chars(solution%values(i, j)%exact_text()))
            end do
        end do
    end subroutine exact_inverse

    ! ------------------------------------------------------------------
    ! Polynomials over Q
    ! ------------------------------------------------------------------

    function vpoly_new(n) result(p)
        integer, intent(in) :: n
        type(vpoly_t) :: p

        allocate (p%c(2, 0:n, 0:n))
        p%c = q(0)
    end function vpoly_new

    subroutine vpoly_axpy(alpha, x, y)
        type(q_t), intent(in) :: alpha
        type(vpoly_t), intent(in) :: x
        type(vpoly_t), intent(inout) :: y

        integer :: a, b, c, n

        n = ubound(y%c, 2)
        do b = 0, n
            do a = 0, n - b
                do c = 1, 2
                    if (q_zero(x%c(c, a, b))) cycle
                    y%c(c, a, b) = y%c(c, a, b) + alpha*x%c(c, a, b)
                end do
            end do
        end do
    end subroutine vpoly_axpy

    subroutine vpoly_set(p, px, py, n)
        type(vpoly_t), intent(inout) :: p
        type(q_t), intent(in) :: px(0:, 0:), py(0:, 0:)
        integer, intent(in) :: n

        p%c(1, 0:n, 0:n) = px(0:n, 0:n)
        p%c(2, 0:n, 0:n) = py(0:n, 0:n)
    end subroutine vpoly_set

    subroutine nedelec_candidates(k, candidates)
        !! Same candidate order as triangle_nedelec_arbitrary_order.
        integer, intent(in) :: k
        type(vpoly_t), allocatable, intent(out) :: candidates(:)

        integer :: component, total, a, n

        allocate (candidates(k*(k + 2)))
        n = 0
        do component = 1, 2
            do total = 0, k - 1
                do a = 0, total
                    n = n + 1
                    candidates(n) = vpoly_new(k)
                    candidates(n)%c(component, a, total - a) = q(1)
                end do
            end do
        end do
        do a = 0, k - 1
            n = n + 1
            candidates(n) = vpoly_new(k)
            candidates(n)%c(1, a, k - a) = q(-1)
            candidates(n)%c(2, a + 1, k - 1 - a) = q(1)
        end do
    end subroutine nedelec_candidates

    function nedelec_moments(p, k) result(m)
        !! Edge moments against shifted Legendre polynomials on the oriented
        !! edges (1,0)->... as in reference_edge, then interior monomial
        !! moments, in the order of build_nedelec_moment_matrix.
        type(vpoly_t), intent(in) :: p
        integer, intent(in) :: k
        type(q_t) :: m(k*(k + 2))

        integer, parameter :: tangent(2, 3) = reshape([1, 0, -1, 1, 0, -1], &
            [2, 3])
        integer :: edge, j, component, total, a, b, row, n, alpha, beta

        n = ubound(p%c, 2)
        row = 0
        do edge = 1, 3
            do j = 0, k - 1
                row = row + 1
                m(row) = q(0)
                do beta = 0, n
                    do alpha = 0, n - beta
                        do component = 1, 2
                            if (tangent(component, edge) == 0) cycle
                            if (q_zero(p%c(component, alpha, beta))) cycle
                            m(row) = m(row) + tangent(component, edge)* &
                                (p%c(component, alpha, beta)* &
                                edge_integral(edge, alpha, beta, j))
                        end do
                    end do
                end do
            end do
        end do
        do component = 1, 2
            do total = 0, k - 2
                do a = 0, total
                    b = total - a
                    row = row + 1
                    m(row) = q(0)
                    do beta = 0, n
                        do alpha = 0, n - beta
                            if (q_zero(p%c(component, alpha, beta))) cycle
                            m(row) = m(row) + p%c(component, alpha, beta)* &
                                triangle_integral(alpha + a, beta + b)
                        end do
                    end do
                end do
            end do
        end do
    end function nedelec_moments

    function edge_integral(edge, a, b, j) result(v)
        !! int_0^1 L_j(s) x(s)^a y(s)^b ds with the edge parametrisations
        !! (s, 0), (1 - s, s) and (0, 1 - s).
        integer, intent(in) :: edge, a, b, j
        type(q_t) :: v

        integer :: s_power, r_power, mm

        select case (edge)
        case (1)
            if (b /= 0) then
                v = q(0)
                return
            end if
            s_power = a
            r_power = 0
        case (2)
            s_power = b
            r_power = a
        case default
            if (a /= 0) then
                v = q(0)
                return
            end if
            s_power = 0
            r_power = b
        end select
        v = q(0)
        do mm = 0, j
            v = v + (q((-1)**(j + mm))*q_binomial(j, mm)*q_binomial(j + mm, mm))* &
                beta_integral(mm + s_power, r_power)
        end do
    end function edge_integral

    function beta_integral(p, r) result(v)
        !! int_0^1 s^p (1 - s)^r ds
        integer, intent(in) :: p, r
        type(q_t) :: v

        v = q_factorial(p)*q_factorial(r)/q_factorial(p + r + 1)
    end function beta_integral

    function triangle_integral(a, b) result(v)
        !! int over the reference triangle of x^a y^b
        integer, intent(in) :: a, b
        type(q_t) :: v

        v = q_factorial(a)*q_factorial(b)/q_factorial(a + b + 2)
    end function triangle_integral

    subroutine lagrange_monomials(k, polys)
        !! Monomial coefficients polys(index(a, b), dof) of the equispaced
        !! Lagrange basis prod_c prod_f (k lambda_c - f)/(f + 1), in the order of
        !! triangle_lagrange_arbitrary_order (lambda_1 = 1 - x - y outer loop).
        integer, intent(in) :: k
        type(q_t), allocatable, intent(out) :: polys(:, :)

        type(q_t), allocatable :: p(:, :), factor(:, :)
        integer :: first, second, third, dof, power(3), c, f

        allocate (polys(0:(k + 1)*(k + 1) - 1, (k + 1)*(k + 2)/2))
        allocate (p(0:k, 0:k), factor(0:k, 0:k))
        dof = 0
        do first = 0, k
            do second = 0, k - first
                third = k - first - second
                power = [first, second, third]
                dof = dof + 1
                p = q(0)
                p(0, 0) = q(1)
                do c = 1, 3
                    do f = 0, power(c) - 1
                        factor = q(0)
                        select case (c)
                        case (1)
                            factor(0, 0) = q_frac(k - f, f + 1)
                            factor(1, 0) = q_frac(-k, f + 1)
                            factor(0, 1) = q_frac(-k, f + 1)
                        case (2)
                            factor(0, 0) = q_frac(-f, f + 1)
                            factor(1, 0) = q_frac(k, f + 1)
                        case (3)
                            factor(0, 0) = q_frac(-f, f + 1)
                            factor(0, 1) = q_frac(k, f + 1)
                        end select
                        p = poly_mul_linear(p, factor, k)
                    end do
                end do
                polys(:, dof) = reshape(p, [(k + 1)*(k + 1)])
            end do
        end do
    end subroutine lagrange_monomials

    function poly_mul_linear(p, l, n) result(r)
        type(q_t), intent(in) :: p(0:, 0:), l(0:, 0:)
        integer, intent(in) :: n
        type(q_t) :: r(0:n, 0:n)

        integer :: a, b

        r = q(0)
        do b = 0, n
            do a = 0, n - b
                if (q_zero(p(a, b))) cycle
                r(a, b) = r(a, b) + l(0, 0)*p(a, b)
                if (a + b < n) then
                    r(a + 1, b) = r(a + 1, b) + l(1, 0)*p(a, b)
                    r(a, b + 1) = r(a, b + 1) + l(0, 1)*p(a, b)
                end if
            end do
        end do
    end function poly_mul_linear

    subroutine scalar_gradient(flat, n, gx, gy)
        type(q_t), intent(in) :: flat(:)
        integer, intent(in) :: n
        type(q_t), allocatable, intent(out) :: gx(:, :), gy(:, :)

        type(q_t) :: p(0:n, 0:n)
        integer :: a, b

        p = reshape(flat, [n + 1, n + 1])
        allocate (gx(0:n, 0:n), gy(0:n, 0:n))
        gx = q(0)
        gy = q(0)
        do b = 0, n
            do a = 0, n - b
                if (q_zero(p(a, b))) cycle
                if (a > 0) gx(a - 1, b) = gx(a - 1, b) + a*p(a, b)
                if (b > 0) gy(a, b - 1) = gy(a, b - 1) + b*p(a, b)
            end do
        end do
    end subroutine scalar_gradient

    function scalar_curl(v, n) result(c)
        !! d v_y / dx - d v_x / dy
        type(vpoly_t), intent(in) :: v
        integer, intent(in) :: n
        type(q_t) :: c(0:n, 0:n)

        integer :: a, b

        c = q(0)
        do b = 0, n
            do a = 0, n - b
                if (a > 0 .and. .not. q_zero(v%c(2, a, b))) &
                    c(a - 1, b) = c(a - 1, b) + a*v%c(2, a, b)
                if (b > 0 .and. .not. q_zero(v%c(1, a, b))) &
                    c(a, b - 1) = c(a, b - 1) - b*v%c(1, a, b)
            end do
        end do
    end function scalar_curl

    function poly_value(p, n, x, y) result(v)
        type(q_t), intent(in) :: p(0:, 0:), x, y
        integer, intent(in) :: n
        type(q_t) :: v

        type(q_t) :: xp(0:n), yp(0:n)
        integer :: a, b

        xp(0) = q(1)
        yp(0) = q(1)
        do a = 1, n
            xp(a) = xp(a - 1)*x
            yp(a) = yp(a - 1)*y
        end do
        v = q(0)
        do b = 0, n
            do a = 0, n - b
                if (q_zero(p(a, b))) cycle
                v = v + p(a, b)*(xp(a)*yp(b))
            end do
        end do
    end function poly_value

    subroutine lagrange_nodes(degree, nodes)
        !! Exact nodes of triangle_lagrange_arbitrary_order.
        integer, intent(in) :: degree
        type(q_t), allocatable, intent(out) :: nodes(:, :)

        integer :: first, second, third, dof

        allocate (nodes(2, (degree + 1)*(degree + 2)/2))
        if (degree == 0) then
            nodes(:, 1) = [q_frac(1, 3), q_frac(1, 3)]
            return
        end if
        dof = 0
        do first = 0, degree
            do second = 0, degree - first
                third = degree - first - second
                dof = dof + 1
                nodes(1, dof) = q_frac(second, degree)
                nodes(2, dof) = q_frac(third, degree)
            end do
        end do
    end subroutine lagrange_nodes

    ! ------------------------------------------------------------------
    ! Bernstein form, B^m_beta = m!/(b1! b2! b3!) l1^b1 l2^b2 l3^b3 with
    ! l1 = 1 - x - y, l2 = x, l3 = y; index order: b2 outer, b3 inner.
    ! ------------------------------------------------------------------

    function to_bernstein(p, n, m) result(c)
        !! x^a y^b = sum_beta (m-a-b)! b2! b3! / (m! (b2-a)! (b3-b)!) B_beta.
        type(q_t), intent(in) :: p(0:, 0:)
        integer, intent(in) :: n, m
        type(q_t) :: c((m + 1)*(m + 2)/2)

        integer :: a, b, b2, b3, idx

        c = q(0)
        do b = 0, n
            do a = 0, n - b
                if (q_zero(p(a, b))) cycle
                if (a + b > m) error stop "Bernstein degree too small"
                idx = 0
                do b2 = 0, m
                    do b3 = 0, m - b2
                        idx = idx + 1
                        if (b2 < a .or. b3 < b) cycle
                        c(idx) = c(idx) + p(a, b)*(q_factorial(m - a - b)* &
                            q_factorial(b2)*q_factorial(b3)/(q_factorial(m)* &
                            q_factorial(b2 - a)*q_factorial(b3 - b)))
                    end do
                end do
            end do
        end do
    end function to_bernstein

    function bernstein_value(c, m, x, y) result(v)
        type(q_t), intent(in) :: c(:), x, y
        integer, intent(in) :: m
        type(q_t) :: v

        type(q_t) :: l(3)
        integer :: b2, b3, idx

        l = [q(1) - x - y, x, y]
        v = q(0)
        idx = 0
        do b2 = 0, m
            do b3 = 0, m - b2
                idx = idx + 1
                if (q_zero(c(idx))) cycle
                v = v + c(idx)*(q_factorial(m)/(q_factorial(m - b2 - b3)* &
                    q_factorial(b2)*q_factorial(b3))*(power(l(1), m - b2 - b3)* &
                    power(l(2), b2)*power(l(3), b3)))
            end do
        end do
    end function bernstein_value

    function power(x, e) result(v)
        type(q_t), intent(in) :: x
        integer, intent(in) :: e
        type(q_t) :: v

        integer :: i

        v = q(1)
        do i = 1, e
            v = v*x
        end do
    end function power

    subroutine check_bernstein(p, c, k)
        !! Both forms agree on the unisolvent P_k lattice, hence identically.
        type(vpoly_t), intent(in) :: p
        type(q_t), intent(in) :: c(:, :)
        integer, intent(in) :: k

        type(q_t), allocatable :: nodes(:, :)
        integer :: l, component

        call lagrange_nodes(k, nodes)
        do l = 1, size(nodes, 2)
            do component = 1, 2
                if (.not. (bernstein_value(c(component, :), k, nodes(1, l), &
                    nodes(2, l)) == poly_value(p%c(component, :, :), k, &
                    nodes(1, l), nodes(2, l)))) &
                    error stop "Bernstein form differs from the basis"
            end do
        end do
    end subroutine check_bernstein

    subroutine check_identity(m, j, message)
        type(q_t), intent(in) :: m(:)
        integer, intent(in) :: j
        character(*), intent(in) :: message

        integer :: i

        do i = 1, size(m)
            if (.not. (m(i) == q(merge(1, 0, i == j)))) error stop message
        end do
    end subroutine check_identity

    subroutine check_same(a, b, n, message)
        type(vpoly_t), intent(in) :: a, b
        integer, intent(in) :: n
        character(*), intent(in) :: message

        integer :: i, j, c

        do j = 0, n
            do i = 0, n - j
                do c = 1, 2
                    if (.not. (a%c(c, i, j) == b%c(c, i, j))) error stop message
                end do
            end do
        end do
    end subroutine check_same

    function negate(v) result(r)
        type(q_t), intent(in) :: v(:)
        type(q_t) :: r(size(v))

        integer :: i

        do i = 1, size(v)
            r(i) = q(0) - v(i)
        end do
    end function negate

    ! ------------------------------------------------------------------
    ! Emission
    ! ------------------------------------------------------------------

    subroutine write_table_header()
        call line("! Generated by fortsym. Do not edit.")
        call line("! Generator: gen_triangle_feec_tables")
        call line("! Generator revision: "//fortsym_revision())
        call line("! Regenerate with: cd tools/codegen && ./generate.sh")
        call line("! Every value is an exact rational rounded once to binary64.")
        call line("")
        call line("module fortfem_generated_triangle_feec_tables")
        call line("    use, intrinsic :: iso_fortran_env, only: dp => real64")
        call line("    implicit none")
        call line("    private")
        call line("")
        call line("    public :: triangle_feec_table_max_order")
        call line("    public :: triangle_nedelec_bernstein_table")
        call line("    public :: triangle_feec_operator_tables")
        call line("")
        call line("    integer, parameter :: triangle_feec_table_max_order = "// &
            itoa(max_order))
        call line("")
        call line("contains")
    end subroutine write_table_header

    subroutine write_bernstein_section()
        type(q_t), allocatable :: bern(:, :, :), gref(:, :)
        type(q_t), allocatable :: cref(:, :), dref(:, :)
        integer :: k, j, c

        call line("")
        call line("    subroutine triangle_nedelec_bernstein_table(order, "// &
            "coefficients, status)")
        call line("        !! coefficients(component, beta, dof) of the order-k")
        call line("        !! first-kind Nedelec basis in degree-k Bernstein form.")
        call line("        integer, intent(in) :: order")
        call line("        real(dp), allocatable, intent(out) :: "// &
            "coefficients(:, :, :)")
        call line("        integer, intent(out) :: status")
        call line("")
        call line("        status = 1")
        call line("        select case (order)")
        do k = 1, max_order
            call build_order(k, bern, gref, cref, dref)
            call store_operator_tables(k, gref, cref, dref)
            call line("        case ("//itoa(k)//")")
            call line("            allocate (coefficients(2, "// &
                itoa(size(bern, 2))//", "//itoa(size(bern, 3))//"))")
            do j = 1, size(bern, 3)
                do c = 1, 2
                    call array_assignment("            coefficients("// &
                        itoa(c)//", :, "//itoa(j)//")", bern(c, :, j))
                end do
            end do
        end do
        call line("        case default")
        call line("            return")
        call line("        end select")
        call line("        status = 0")
        call line("    end subroutine triangle_nedelec_bernstein_table")
    end subroutine write_bernstein_section

    subroutine store_operator_tables(k, gref, cref, dref)
        !! Buffer the operator tables per order; they are written after the
        !! Bernstein section so the exact construction runs once per order.
        integer, intent(in) :: k
        type(q_t), intent(in) :: gref(:, :), cref(:, :), dref(:, :)

        character(:), allocatable :: path
        integer :: u, i, j

        path = scratch_path(k)
        open(newunit=u, file=path, status="replace", action="write")
        do i = 1, size(gref, 2)
            call array_assignment_unit(u, "            gradient(:, "// &
                itoa(i)//")", gref(:, i))
        end do
        do j = 1, size(cref, 2)
            call array_assignment_unit(u, "            curl(:, "// &
                itoa(j)//")", cref(:, j))
        end do
        do j = 1, size(dref, 2)
            call array_assignment_unit(u, "            divergence(:, "// &
                itoa(j)//")", dref(:, j))
        end do
        close(u)
    end subroutine store_operator_tables

    function scratch_path(k) result(path)
        integer, intent(in) :: k
        character(:), allocatable :: path

        path = generated_path(".gen_triangle_feec_operators_"//itoa(k)//".tmp")
    end function scratch_path

    subroutine write_operator_section()
        character(4096) :: buffer
        integer :: k, u, ios, ndof, nl0, nld

        call line("")
        call line("    subroutine triangle_feec_operator_tables(order, "// &
            "gradient, curl, divergence, status)")
        call line("        !! gradient(ndof, nl0): Nedelec dofs of the P_k "// &
            "Lagrange gradients.")
        call line("        !! curl(nld, ndof): curl of the Nedelec basis at the "// &
            "P_{k-1} nodes.")
        call line("        !! divergence(nld, ndof): divergence of the RT_{k-1} "// &
            "basis at the P_{k-1} nodes.")
        call line("        integer, intent(in) :: order")
        call line("        real(dp), allocatable, intent(out) :: gradient(:, :), "// &
            "curl(:, :), divergence(:, :)")
        call line("        integer, intent(out) :: status")
        call line("")
        call line("        status = 1")
        call line("        select case (order)")
        do k = 1, max_order
            ndof = k*(k + 2)
            nl0 = (k + 1)*(k + 2)/2
            nld = k*(k + 1)/2
            call line("        case ("//itoa(k)//")")
            call line("            allocate (gradient("//itoa(ndof)//", "// &
                itoa(nl0)//"), curl("//itoa(nld)//", "//itoa(ndof)// &
                "), divergence("//itoa(nld)//", "//itoa(ndof)//"))")
            open(newunit=u, file=scratch_path(k), status="old", action="read")
            do
                read(u, "(a)", iostat=ios) buffer
                if (ios /= 0) exit
                call line(trim(buffer))
            end do
            close(u, status="delete")
        end do
        call line("        case default")
        call line("            return")
        call line("        end select")
        call line("        status = 0")
        call line("    end subroutine triangle_feec_operator_tables")
    end subroutine write_operator_section

    subroutine array_assignment(lhs, values)
        character(*), intent(in) :: lhs
        type(q_t), intent(in) :: values(:)

        call array_assignment_unit(table_unit, lhs, values)
    end subroutine array_assignment

    subroutine array_assignment_unit(u, lhs, values)
        !! Literals are packed into continuation lines of at most 88 columns.
        integer, intent(in) :: u
        character(*), intent(in) :: lhs
        type(q_t), intent(in) :: values(:)

        character(:), allocatable :: current, item
        integer :: i

        if (all_zero(values)) then
            write(u, "(a)") lhs//" = 0.0_dp"
            return
        end if
        write(u, "(a)") lhs//" = [ &"
        current = "                "
        do i = 1, size(values)
            item = q_literal(values(i))
            if (i < size(values)) then
                item = item//","
            else
                item = item//"]"
            end if
            if (len(current) + len(item) + 3 > 88) then
                write(u, "(a)") current//" &"
                current = "                "
            end if
            if (len(current) > 16) current = current//" "
            current = current//item
        end do
        write(u, "(a)") current
    end subroutine array_assignment_unit

    logical function all_zero(values)
        type(q_t), intent(in) :: values(:)

        integer :: i

        all_zero = .false.
        do i = 1, size(values)
            if (.not. q_zero(values(i))) return
        end do
        all_zero = .true.
    end function all_zero

    subroutine line(text)
        character(*), intent(in) :: text

        write(table_unit, "(a)") text
    end subroutine line

    function itoa(i) result(s)
        integer, intent(in) :: i
        character(:), allocatable :: s

        character(16) :: buffer

        write(buffer, "(i0)") i
        s = trim(buffer)
    end function itoa

    ! ------------------------------------------------------------------
    ! Bernstein jets: values, gradients and Hessians in (xi, eta)
    ! ------------------------------------------------------------------

    subroutine write_bernstein_jets()
        integer :: u, m

        open(newunit=u, file=generated_path("fortfem_triangle_bernstein_jets.f90"), &
            status="replace", action="write")
        do m = 0, max_order
            call bernstein_jet(u, m)
        end do
        close(u)
    end subroutine write_bernstein_jets

    subroutine bernstein_jet(u, m)
        !! B^m_beta and its xi/eta derivatives through second order, in terms
        !! of the barycentric inputs; d/dxi = d/dl2 - d/dl1, d/deta = d/dl3 - d/dl1.
        integer, intent(in) :: u, m

        type(kernel_spec_t) :: spec
        type(expr_t) :: l(3), b, bx, by
        type(expr_t), allocatable :: roots(:)
        type(str_t), allocatable :: refs(:)
        type(engine_result_t) :: result
        character(:), allocatable :: code
        integer :: nb, b2, b3, idx, i
        type(q_t) :: coefficient

        nb = (m + 1)*(m + 2)/2
        allocate (roots(6*nb), refs(6*nb))
        l(1) = sym(arena, "lambda_1")
        l(2) = sym(arena, "lambda_2")
        l(3) = sym(arena, "lambda_3")
        idx = 0
        do b2 = 0, m
            do b3 = 0, m - b2
                idx = idx + 1
                coefficient = q_factorial(m)/(q_factorial(m - b2 - b3)* &
                    q_factorial(b2)*q_factorial(b3))
                b = parse_coefficient(coefficient)*l(1)**(m - b2 - b3)* &
                    l(2)**b2*l(3)**b3
                bx = diff(b, l(2)) - diff(b, l(1))
                by = diff(b, l(3)) - diff(b, l(1))
                roots(idx) = b
                roots(nb + idx) = bx
                roots(2*nb + idx) = by
                roots(3*nb + idx) = diff(bx, l(2)) - diff(bx, l(1))
                roots(4*nb + idx) = diff(bx, l(3)) - diff(bx, l(1))
                roots(5*nb + idx) = diff(by, l(3)) - diff(by, l(1))
                refs(idx) = str("value("//itoa(idx)//")")
                refs(nb + idx) = str("gradient(1,"//itoa(idx)//")")
                refs(2*nb + idx) = str("gradient(2,"//itoa(idx)//")")
                refs(3*nb + idx) = str("hessian(1,"//itoa(idx)//")")
                refs(4*nb + idx) = str("hessian(2,"//itoa(idx)//")")
                refs(5*nb + idx) = str("hessian(3,"//itoa(idx)//")")
            end do
        end do
        do i = 1, size(roots)
            result = engine%simplify(roots(i))
            if (.not. result%ok) error stop "Bernstein jet simplification failed"
            roots(i) = result%value
        end do
        spec%name = str("generated_triangle_bernstein_jet_"//itoa(m))
        spec%module_name = str("fortfem_generated_triangle_bernstein_jet_"// &
            itoa(m))
        spec%mode = KERNEL_SUBROUTINE
        spec%pure_procedure = .true.
        spec%generator = str("gen_triangle_feec_tables")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        spec%args = [str("lambda_1"), str("lambda_2"), str("lambda_3")]
        spec%outputs = [str("value"), str("gradient"), str("hessian")]
        spec%output_shapes = [str("("//itoa(nb)//")"), &
            str("(2,"//itoa(nb)//")"), str("(3,"//itoa(nb)//")")]
        spec%output_references = refs
        code = chars(emit_kernel(roots, spec))
        write(u, "(a)") code(:len(code) - 1)
    end subroutine bernstein_jet

    function parse_coefficient(c) result(e)
        type(q_t), intent(in) :: c
        type(expr_t) :: e

        character(:), allocatable :: message
        logical :: ok

        e = parse_expr(arena, chars(c%value), ok, message)
        if (.not. ok) error stop "cannot parse Bernstein coefficient"
    end function parse_coefficient

end program gen_triangle_feec_tables
