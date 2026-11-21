module spline_cohomology_xpoint
    !! u - lap u = f with exact Neumann data, u = cos(0.7 R) sin(0.9 Z) + R Z,
    !! and a double-null geometry: Gamma_p made of four quadratic Bezier arcs
    !! with 90 degree corners at the X-points (R0, -h), (R0, h), straight legs
    !! with a slight bend, and an elliptic wall.
    use fortfem_kinds, only: dp
    use fortfem_spline_energy, only: pointwise_energy_t
    use fortfem_spline_xpoint, only: xpoint_curves_t
    implicit none
    private

    public :: problem_kernel_t, exact_solution, exact_gradient, source, &
        double_null, core_area

    real(dp), parameter, public :: r0 = 1.7_dp
    real(dp), parameter, public :: wall_a = 1.05_dp, wall_b = 1.6_dp

    type, extends(pointwise_energy_t) :: problem_kernel_t
        real(dp), allocatable :: f(:)
    contains
        procedure :: evaluate
    end type problem_kernel_t

contains

    subroutine evaluate(self, cell, first_point, npoint, jet, weight, need, w, &
            gradient, hessian, status)
        class(problem_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status

        real(dp) :: f
        integer :: q

        associate (unused => cell, unused_weight => weight)
        end associate
        do q = 1, npoint
            f = self%f(first_point + q - 1)
            w(q) = 0.5_dp*(jet(2, q)**2 + jet(3, q)**2 + jet(1, q)**2) - f*jet(1, q)
            if (need >= 1) gradient(:, q) = [jet(1, q) - f, jet(2, q), jet(3, q)]
            if (need >= 2) hessian(:, q) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
                0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
        end do
        status = 0
    end subroutine evaluate

    pure real(dp) function exact_solution(s) result(u)
        real(dp), intent(in) :: s(2)

        u = cos(0.7_dp*s(1))*sin(0.9_dp*s(2)) + s(1)*s(2)
    end function exact_solution

    pure function exact_gradient(s) result(g)
        real(dp), intent(in) :: s(2)
        real(dp) :: g(2)

        g = [-0.7_dp*sin(0.7_dp*s(1))*sin(0.9_dp*s(2)) + s(2), &
            0.9_dp*cos(0.7_dp*s(1))*cos(0.9_dp*s(2)) + s(1)]
    end function exact_gradient

    pure real(dp) function source(s) result(f)
        real(dp), intent(in) :: s(2)

        f = exact_solution(s) + 1.3_dp*cos(0.7_dp*s(1))*sin(0.9_dp*s(2))
    end function source

    pure function bezier_controls(a, h, shift) result(c)
        !! (2, 3, 4): control points of the four Gamma_p arcs.
        real(dp), intent(in) :: a, h, shift
        real(dp) :: c(2, 3, 4)
        real(dp) :: k(2, 4)

        k(:, 1) = [r0 + shift, -h]
        k(:, 2) = [r0 + a, 0.0_dp]
        k(:, 3) = [r0 - shift, h]
        k(:, 4) = [r0 - a, 0.0_dp]
        c(:, :, 1) = reshape([k(:, 1), k(:, 1) + [a - shift, a - shift], k(:, 2)], [2, 3])
        c(:, :, 2) = reshape([k(:, 2), k(:, 3) + [a + shift, -(a + shift)], k(:, 3)], &
            [2, 3])
        c(:, :, 3) = reshape([k(:, 3), k(:, 3) + [-(a - shift), -(a - shift)], &
            k(:, 4)], [2, 3])
        c(:, :, 4) = reshape([k(:, 4), k(:, 1) + [-(a + shift), a + shift], k(:, 1)], &
            [2, 3])
    end function bezier_controls

    function double_null(a, h, shift, legs) result(curves)
        !! legs(q): corner q is an X-point (only corners 1 and 3 may be).
        real(dp), intent(in) :: a, h, shift
        logical, intent(in) :: legs(4)
        type(xpoint_curves_t) :: curves

        real(dp) :: c(2, 3, 4), t, theta
        integer :: q, i
        integer, parameter :: n = 40

        c = bezier_controls(a, h, shift)
        do q = 1, 4
            allocate(curves%segment(q)%points(2, n))
            do i = 1, n
                t = real(i - 1, dp)/(n - 1)
                curves%segment(q)%points(:, i) = (1.0_dp - t)**2*c(:, 1, q) + &
                    2.0_dp*t*(1.0_dp - t)*c(:, 2, q) + t**2*c(:, 3, q)
            end do
        end do
        if (legs(1)) then
            curves%leg_after(1)%points = leg(c(:, 1, 1), [1.0_dp, -1.0_dp])
            curves%leg_before(1)%points = leg(c(:, 1, 1), [-1.0_dp, -1.0_dp])
        end if
        if (legs(3)) then
            curves%leg_after(3)%points = leg(c(:, 1, 3), [-1.0_dp, 1.0_dp])
            curves%leg_before(3)%points = leg(c(:, 1, 3), [1.0_dp, 1.0_dp])
        end if
        allocate(curves%wall(2, 300))
        do i = 1, 300
            theta = 2.0_dp*acos(-1.0_dp)*(i - 1)/300
            curves%wall(:, i) = [r0 + wall_a*cos(theta), wall_b*sin(theta)]
        end do
    end function double_null

    function leg(x, d) result(points)
        !! From x along d to the wall, bent by 3 % of its length.
        real(dp), intent(in) :: x(2), d(2)
        real(dp), allocatable :: points(:, :)

        real(dp) :: e(2), lo, hi, mid, s(2), nrm(2), t
        integer :: i, it

        e = d/norm2(d)
        lo = 0.0_dp
        hi = 5.0_dp
        do it = 1, 80
            mid = 0.5_dp*(lo + hi)
            s = x + mid*e
            if (((s(1) - r0)/wall_a)**2 + (s(2)/wall_b)**2 < 1.0_dp) then
                lo = mid
            else
                hi = mid
            end if
        end do
        nrm = [-e(2), e(1)]
        allocate(points(2, 25))
        do i = 1, 25
            t = real(i - 1, dp)/24
            points(:, i) = x + t*lo*e + 0.12_dp*lo*t*(1.0_dp - t)*nrm
        end do
    end function leg

    real(dp) function core_area(a, h, shift) result(area)
        !! int R dZ over Gamma_p (fine trapezoidal rule on the arcs).
        real(dp), intent(in) :: a, h, shift
        real(dp) :: c(2, 3, 4), t, x(2), dx(2)
        integer :: q, i
        integer, parameter :: n = 20000

        c = bezier_controls(a, h, shift)
        area = 0.0_dp
        do q = 1, 4
            do i = 0, n - 1
                t = (i + 0.5_dp)/n
                x = (1.0_dp - t)**2*c(:, 1, q) + 2.0_dp*t*(1.0_dp - t)*c(:, 2, q) + &
                    t**2*c(:, 3, q)
                dx = 2.0_dp*(1.0_dp - t)*(c(:, 2, q) - c(:, 1, q)) + &
                    2.0_dp*t*(c(:, 3, q) - c(:, 2, q))
                area = area + x(1)*dx(2)/n
            end do
        end do
    end function core_area

end module spline_cohomology_xpoint

module spline_cohomology_geometry
    !! Jacobians and coordinate fields of the current section space and chart.
    use fortfem_kinds, only: dp
    use fortfem_spline_space, only: spline_space_t
    use fortfem_spline_eval, only: spline_geometry_at
    use fortfem_curve_chart, only: curve_chart_t, curve_chart_map
    implicit none
    private

    type(spline_space_t), public :: geometry_space
    type(curve_chart_t), public :: geometry_chart
    logical, public :: three_d = .false.
    integer, public :: coordinate = 1

    public :: jacobian, coordinate_field

contains

    subroutine jacobian(block, xi, t, jac)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: jac(:, :)
        real(dp) :: jet(2, 6), x(3), jc(3, 3)
        integer :: status

        call spline_geometry_at(geometry_space, block, xi, jet, status)
        if (.not. three_d) then
            jac = jet(:, 2:3)
            return
        end if
        call curve_chart_map(geometry_chart, jet(:, 1), t, x, jc)
        jac(:, 1:2) = matmul(jc(:, 1:2), jet(:, 2:3))
        jac(:, 3) = jc(:, 3)
    end subroutine jacobian

    subroutine coordinate_field(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: jet(2, 6), x(3), jc(3, 3)
        integer :: status

        call spline_geometry_at(geometry_space, block, xi, jet, status)
        if (three_d) then
            call curve_chart_map(geometry_chart, jet(:, 1), t, x, jc)
            values(1) = x(coordinate)
        else
            values(1) = jet(coordinate, 1)
        end if
    end subroutine coordinate_field

end module spline_cohomology_geometry

program test_spline_cohomology
    !! Tree, cotree and first-cohomology representatives of the spline de
    !! Rham complex: the tree has (vertices - 1) edges, the number of
    !! representatives equals b1 (disc 0, ring 1, disc x circle 1,
    !! ring x circle 2), each is closed (curl h = 0), vanishes on the tree and
    !! lies outside the range of grad. The boundary surface of ring x circle
    !! has a surface curl that annihilates gradient traces and two closed,
    !! non-exact surface classes; mass matrices are symmetric and reproduce
    !! independently integrated norms.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, &
        layout_from_blocks, LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_tree, spline_derham_harmonic, spline_derham_surface, &
        spline_derham_mass, spline_derham_evaluate, spline_derham_boundary_dofs, &
        spline_derham_interpolation, spline_derham_project
    use fortfem_curve_chart, only: curve_chart_cylindrical
    use spline_cohomology_geometry, only: geometry_space, geometry_chart, three_d, &
        coordinate, jacobian, coordinate_field
    use fortfem_spline_xpoint, only: layout_xpoint
    use fortfem_spline_layout, only: LAYOUT_PLASMA_BOUNDARY
    use spline_cohomology_xpoint, only: double_null, r0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build
    implicit none

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine dsyev
    end interface

    type(block_layout_t) :: disc, ring, xpoint
    integer, allocatable :: gp_nodes(:), gp_edges(:)
    type(spline_derham_t) :: xd
    integer :: status

    call layout_ogrid([0.0_dp, 0.0_dp], 1.0_dp, 0.5_dp, disc, status, polygon=.true.)
    call ring_layout(ring)
    call check_condition(cohomology_ok(disc, -1, 0), "disc: tree and no loop")
    call check_condition(cohomology_ok(ring, -1, 1), "ring: one closed non-exact form")
    call check_condition(cohomology_ok(disc, 1, 1), "disc x circle: the toroidal class")
    call check_condition(cohomology_ok(ring, 1, 2), &
        "ring x circle: toroidal and poloidal classes")
    call layout_xpoint([r0, 0.0_dp], 0.4_dp, double_null(0.6_dp, 1.0_dp, 0.0_dp, &
        [.true., .false., .true., .false.]), xpoint, status)
    call check_condition(cohomology_ok(xpoint, -1, 0) .and. &
        cohomology_ok(xpoint, 1, 1), "double null: section and 3D cohomology")
    call build(xpoint, -1, xd)
    call spline_derham_boundary_dofs(xd, LAYOUT_PLASMA_BOUNDARY, 0, gp_nodes)
    call spline_derham_boundary_dofs(xd, LAYOUT_PLASMA_BOUNDARY, 1, gp_edges)
    call check_condition(size(gp_nodes) > 0 .and. size(gp_nodes) == size(gp_edges), &
        "double null: Gamma_p is a closed loop of nodes and tangential edges")
    call check_condition(surface_ok(ring), &
        "boundary torus: surface curl of gradients vanishes, two surface classes")
    call check_condition(mass_ok(ring, -1), "section mass matrices of all forms")
    call check_condition(mass_ok(ring, 1), "3D mass matrices of all forms")
    call check_condition(interpolation_ok(disc, -1), &
        "2D auxiliary-space interpolation maps constants to grad of coordinates")
    call layout_ogrid([2.0_dp, 0.0_dp], 0.5_dp, 0.5_dp, disc, status, polygon=.true.)
    call check_condition(interpolation_ok(disc, 1), &
        "3D auxiliary-space interpolation maps constants to grad of coordinates")
    call check_summary("spline_cohomology")

contains

    logical function cohomology_ok(layout, modes, b1) result(ok)
        type(block_layout_t), intent(in) :: layout
        integer, intent(in) :: modes, b1

        type(spline_space_t) :: space
        type(spline_space_options_t) :: options
        type(spline_derham_t) :: derham
        logical, allocatable :: tree(:)
        real(dp), allocatable :: h(:, :), g(:, :), c(:, :), gh(:, :)
        integer :: classes, status

        options%degree = 2
        allocate(options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = 3
        call spline_space_build(layout, options, space, status)
        if (modes >= 0) then
            call spline_derham_build(space, derham, status, modes)
        else
            call spline_derham_build(space, derham, status)
        end if
        call spline_derham_tree(derham, tree)
        call spline_derham_harmonic(derham, h, classes)
        g = dense(derham%grad, derham%dimension(0))
        c = dense(derham%curl, derham%dimension(1))
        ok = count(tree) == derham%dimension(0) - 1 .and. classes == b1
        if (.not. ok) print *, "tree", count(tree), derham%dimension(0), "classes", classes
        if (classes == 0) return
        ok = ok .and. maxval(abs(matmul(c, h))) < 1.0e-12_dp
        ok = ok .and. maxval(abs(h(pack([(status, status = 1, size(tree))], tree), &
            :))) == 0.0_dp
        gh = reshape([g, h], [size(g, 1), size(g, 2) + classes])
        ok = ok .and. rank(gh) == rank(g) + classes
    end function cohomology_ok

    subroutine build(layout, modes, derham)
        type(block_layout_t), intent(in) :: layout
        integer, intent(in) :: modes
        type(spline_derham_t), intent(out) :: derham

        type(spline_space_t) :: space
        type(spline_space_options_t) :: options
        integer :: status

        options%degree = 2
        allocate(options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = 3
        call spline_space_build(layout, options, space, status)
        if (modes >= 0) then
            call spline_derham_build(space, derham, status, modes)
        else
            call spline_derham_build(space, derham, status)
        end if
    end subroutine build

    logical function surface_ok(layout) result(ok)
        type(block_layout_t), intent(in) :: layout

        type(spline_derham_t) :: derham
        type(csr_matrix_t) :: sc
        integer, allocatable :: edges(:), faces(:), nodes(:)
        real(dp), allocatable :: h(:, :), phi(:), gs(:, :), s(:, :), hs(:, :), g(:, :)

        call build(layout, 1, derham)
        call spline_derham_surface(derham, LAYOUT_WALL, edges, faces, sc, h)
        call spline_derham_boundary_dofs(derham, LAYOUT_WALL, 0, nodes)
        g = dense(derham%grad, derham%dimension(0))
        gs = g(edges, nodes)
        s = dense(sc, size(edges))
        allocate(phi(size(nodes)))
        call random_number(phi)
        hs = h(edges, :)
        ok = maxval(abs(matmul(s, matmul(gs, phi)))) < 1.0e-12_dp
        ok = ok .and. maxval(abs(matmul(s, hs))) < 1.0e-12_dp
        ok = ok .and. rank(reshape([gs, hs], [size(gs, 1), size(gs, 2) + 2])) == &
            rank(gs) + 2
    end function surface_ok

    logical function mass_ok(layout, modes) result(ok)
        !! Identity metric: c^T M c equals the integral of |proxy|^2 by point
        !! evaluation on a different quadrature, and M is symmetric.
        type(block_layout_t), intent(in) :: layout
        integer, intent(in) :: modes

        type(spline_derham_t) :: derham
        type(spline_quadrature_t) :: quad, fine
        type(csr_matrix_t) :: m
        real(dp), allocatable :: c(:), md(:, :), v(:)
        real(dp) :: direct, tt
        integer :: form, status, q, j, ncomp, top, nt

        call build(layout, modes, derham)
        call spline_quadrature_build(derham%space, 3, quad, status)
        call spline_quadrature_build(derham%space, 4, fine, status)
        ok = .true.
        top = merge(3, 2, modes >= 0)
        nt = derham%plane_count
        do form = 0, top
            call spline_derham_mass(derham, quad, form, identity, m, status)
            ok = ok .and. status == 0
            allocate(c(derham%dimension(form)))
            call random_number(c)
            md = dense(m, size(c))
            ok = ok .and. maxval(abs(md - transpose(md))) < 1.0e-12_dp*maxval(abs(md))
            ncomp = merge(3, merge(2, 1, form == 1), modes >= 0 .and. &
                (form == 1 .or. form == 2))
            allocate(v(ncomp))
            direct = 0.0_dp
            do j = 1, merge(4*nt, 1, modes >= 0)
                tt = (j - 0.5_dp)/(4*nt)
                do q = 1, fine%point_count
                    call spline_derham_evaluate(derham, form, c, fine%block(q), &
                        fine%xi(:, q), tt, v)
                    direct = direct + fine%weight(q)*sum(v**2)/merge(4*nt, 1, modes >= 0)
                end do
            end do
            ok = ok .and. abs(dot_product(c, matmul(md, c)) - direct) < &
                1.0e-8_dp*direct
            if (.not. ok) print *, "form", form, dot_product(c, matmul(md, c)), direct
            deallocate(c, v)
        end do
    end function mass_ok

    logical function interpolation_ok(layout, modes) result(ok)
        type(block_layout_t), intent(in) :: layout
        integer, intent(in) :: modes

        type(spline_derham_t) :: derham
        type(csr_matrix_t) :: pi
        real(dp), allocatable :: ones(:), p0(:), ref(:), got(:)
        integer :: d, dimx, n0, status, i, k

        call build(layout, modes, derham)
        geometry_space = derham%space
        three_d = modes >= 0
        call curve_chart_cylindrical(geometry_chart, 1)
        call spline_derham_interpolation(derham, jacobian, pi, status, 1.0e-10_dp)
        ok = status == 0
        dimx = merge(3, 2, three_d)
        n0 = derham%dimension(0)
        allocate(p0(n0), got(pi%n))
        do d = 1, dimx
            coordinate = d
            call spline_derham_project(derham, 0, coordinate_field, p0, status, 12)
            ref = matvec(derham%grad, p0)
            ones = [(merge(1.0_dp, 0.0_dp, (i - 1)/n0 + 1 == d), i = 1, dimx*n0)]
            got = matvec(pi, ones)
            ok = ok .and. maxval(abs(got - ref)) < 1.0e-8_dp*max(1.0_dp, maxval(abs(ref)))
            if (.not. ok) print *, "d", d, maxval(abs(got - ref))
        end do
        print '(a,i0,a,f7.1)', "Pi columns ", dimx*n0, ", nonzeros per column ", &
            real(size(pi%val), dp)/(dimx*n0)
        associate (unused => k)
        end associate
    end function interpolation_ok

    function matvec(a, x) result(y)
        type(csr_matrix_t), intent(in) :: a
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: y(:)
        integer :: i, k

        allocate(y(a%n))
        y = 0.0_dp
        do i = 1, a%n
            do k = a%rowptr(i), a%rowptr(i + 1) - 1
                y(i) = y(i) + a%val(k)*x(a%colidx(k))
            end do
        end do
    end function matvec

    subroutine identity(block, xi, t, weight)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: weight(:, :)
        integer :: i

        weight = 0.0_dp
        do i = 1, size(weight, 1)
            weight(i, i) = 1.0_dp
        end do
        associate (unused => block, u2 => xi, u3 => t)
        end associate
    end subroutine identity

    subroutine ring_layout(layout)
        type(block_layout_t), intent(out) :: layout
        real(dp) :: v(2, 8)
        integer :: corners(4, 4), tags(4, 4), status

        v = reshape([-0.5_dp, -0.5_dp, 0.5_dp, -0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp, &
            0.5_dp, -1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, &
            1.0_dp], [2, 8])
        corners = reshape([5, 6, 2, 1, 6, 7, 3, 2, 7, 8, 4, 3, 8, 5, 1, 4], [4, 4])
        tags = LAYOUT_WALL
        call layout_from_blocks(v, corners, tags, layout, status)
    end subroutine ring_layout

    function dense(a, ncol) result(m)
        type(csr_matrix_t), intent(in) :: a
        integer, intent(in) :: ncol
        real(dp), allocatable :: m(:, :)
        integer :: i, k

        allocate(m(a%n, ncol))
        m = 0.0_dp
        do i = 1, a%n
            do k = a%rowptr(i), a%rowptr(i + 1) - 1
                m(i, a%colidx(k)) = m(i, a%colidx(k)) + a%val(k)
            end do
        end do
    end function dense

    integer function rank(m) result(r)
        real(dp), intent(in) :: m(:, :)
        real(dp), allocatable :: g(:, :), w(:), work(:)
        integer :: n, info

        n = size(m, 2)
        g = matmul(transpose(m), m)
        allocate(w(n), work(10*n))
        call dsyev('N', 'U', n, g, n, w, work, 10*n, info)
        r = count(w > 1.0e-9_dp*maxval(w))
    end function rank

end program test_spline_cohomology
