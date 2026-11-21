module spline_derham_fields
    !! Fields for the de Rham tests and their parametric proxies on the
    !! current section space (and chart, in 3D).
    use fortfem_kinds, only: dp
    use fortfem_spline_space, only: spline_space_t
    use fortfem_spline_eval, only: spline_geometry_at
    use fortfem_curve_chart, only: curve_chart_t, curve_chart_map
    implicit none
    private

    type(spline_space_t), public :: space
    type(curve_chart_t), public :: chart

    public :: f0, grad_f, a1, curl_a, f3d, grad_f3d, a3d, curl_a3d, b3d, div_b3d

contains

    subroutine section(block, xi, s, j)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2)
        real(dp), intent(out) :: s(2), j(2, 2)
        real(dp) :: jet(2, 6)
        integer :: status

        call spline_geometry_at(space, block, xi, jet, status)
        s = jet(:, 1)
        j = jet(:, 2:3)
    end subroutine section

    subroutine full(block, xi, t, x, j3)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: x(3), j3(3, 3)
        real(dp) :: s(2), j(2, 2), jc(3, 3)

        call section(block, xi, s, j)
        call curve_chart_map(chart, s, t, x, jc)
        j3(:, 1:2) = matmul(jc(:, 1:2), j)
        j3(:, 3) = jc(:, 3)
    end subroutine full

    subroutine f0(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: s(2), j(2, 2)

        call section(block, xi, s, j)
        values(1) = sin(1.3_dp*s(1))*cos(0.7_dp*s(2)) + s(1)*s(2)
        associate (unused => t)
        end associate
    end subroutine f0

    subroutine grad_f(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: s(2), j(2, 2), g(2)

        call section(block, xi, s, j)
        g = [1.3_dp*cos(1.3_dp*s(1))*cos(0.7_dp*s(2)) + s(2), &
            -0.7_dp*sin(1.3_dp*s(1))*sin(0.7_dp*s(2)) + s(1)]
        values(1:2) = matmul(transpose(j), g)
        associate (unused => t)
        end associate
    end subroutine grad_f

    subroutine a1(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: s(2), j(2, 2)

        call section(block, xi, s, j)
        values(1:2) = matmul(transpose(j), [-s(2)**2, s(1)*s(2) + sin(s(1))])
        associate (unused => t)
        end associate
    end subroutine a1

    subroutine curl_a(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: s(2), j(2, 2)

        call section(block, xi, s, j)
        values(1) = (3.0_dp*s(2) + cos(s(1)))*(j(1, 1)*j(2, 2) - j(1, 2)*j(2, 1))
        associate (unused => t)
        end associate
    end subroutine curl_a

    subroutine f3d(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: x(3), j3(3, 3)

        call full(block, xi, t, x, j3)
        values(1) = x(1) + x(2)*x(3) + x(1)**2
    end subroutine f3d

    subroutine grad_f3d(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: x(3), j3(3, 3)

        call full(block, xi, t, x, j3)
        values(1:3) = matmul(transpose(j3), [1.0_dp + 2.0_dp*x(1), x(3), x(2)])
    end subroutine grad_f3d

    subroutine a3d(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: x(3), j3(3, 3)

        call full(block, xi, t, x, j3)
        values(1:3) = matmul(transpose(j3), [x(2)*x(3), -x(1)**2, x(1)*x(2)])
    end subroutine a3d

    subroutine curl_a3d(block, xi, t, values)
        !! curl (x2 x3, -x1^2, x1 x2) = (x1, x2 - x2, -2 x1 - x3), Piola.
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: x(3), j3(3, 3)

        call full(block, xi, t, x, j3)
        values(1:3) = piola(j3, [x(1), 0.0_dp, -2.0_dp*x(1) - x(3)])
    end subroutine curl_a3d

    subroutine b3d(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: x(3), j3(3, 3)

        call full(block, xi, t, x, j3)
        values(1:3) = piola(j3, [x(1)*x(2), x(3), x(1)**2])
    end subroutine b3d

    subroutine div_b3d(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: x(3), j3(3, 3)

        call full(block, xi, t, x, j3)
        values(1) = x(2)*det3(j3)
    end subroutine div_b3d

    pure function piola(j, b) result(p)
        !! det(J) J^-1 b = adj(J) b; the columns of c are the rows of adj(J).
        real(dp), intent(in) :: j(3, 3), b(3)
        real(dp) :: p(3), c(3, 3)

        c(:, 1) = [j(2, 2)*j(3, 3) - j(3, 2)*j(2, 3), j(3, 2)*j(1, 3) - j(1, 2)*j(3, 3), &
            j(1, 2)*j(2, 3) - j(2, 2)*j(1, 3)]
        c(:, 2) = [j(3, 1)*j(2, 3) - j(2, 1)*j(3, 3), j(1, 1)*j(3, 3) - j(3, 1)*j(1, 3), &
            j(2, 1)*j(1, 3) - j(1, 1)*j(2, 3)]
        c(:, 3) = [j(2, 1)*j(3, 2) - j(3, 1)*j(2, 2), j(3, 1)*j(1, 2) - j(1, 1)*j(3, 2), &
            j(1, 1)*j(2, 2) - j(2, 1)*j(1, 2)]
        p = matmul(transpose(c), b)
    end function piola

    pure real(dp) function det3(j)
        real(dp), intent(in) :: j(3, 3)

        det3 = j(1, 1)*(j(2, 2)*j(3, 3) - j(3, 2)*j(2, 3)) - &
            j(1, 2)*(j(2, 1)*j(3, 3) - j(3, 1)*j(2, 3)) + &
            j(1, 3)*(j(2, 1)*j(3, 2) - j(3, 1)*j(2, 2))
    end function det3

end module spline_derham_fields

program test_spline_derham
    !! Gate 5 of doc/design/multipatch_spline_spaces.md: the spline de Rham
    !! complex of a multipatch section and its extension by a Fourier
    !! direction. d d = 0, Betti numbers from ranks (disc 1 0 0, ring 1 1 0,
    !! disc x circle 1 1 0 0, ring x circle 1 2 1 0), commuting projections
    !! d Pi = Pi d in 2D and 3D, projection convergence and boundary traces.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, &
        layout_from_blocks, LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_options_t, spline_space_build
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_project, spline_derham_evaluate, spline_derham_boundary_dofs, &
        spline_derham_cell_basis
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build
    use fortfem_curve_chart, only: curve_chart_cylindrical
    use spline_derham_fields, only: space, chart, f0, grad_f, a1, curl_a, f3d, &
        grad_f3d, a3d, curl_a3d, b3d, div_b3d
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

    type(block_layout_t) :: disc, ring
    type(spline_derham_t) :: derham
    real(dp), allocatable :: c0(:), c1(:), c2(:), c3(:), d1(:), d2(:), d3(:)
    integer, allocatable :: nodes(:), edges(:), faces(:)
    real(dp) :: errors(2, 2)
    integer :: status, level

    call layout_ogrid([0.0_dp, 0.0_dp], 1.0_dp, 0.5_dp, disc, status, polygon=.true.)
    call ring_layout(ring)

    call make(disc, 2, 3, -1)
    call check_condition(status == 0 .and. dd_zero(derham) .and. &
        all(betti(derham) == [1, 0, 0, 0]), "disc: d d = 0, Betti numbers 1 0 0")
    call make(ring, 2, 3, -1)
    call check_condition(status == 0 .and. dd_zero(derham) .and. &
        all(betti(derham) == [1, 1, 0, 0]), "ring: d d = 0, Betti numbers 1 1 0")
    call make(disc, 2, 2, 1)
    call check_condition(status == 0 .and. dd_zero(derham) .and. &
        all(betti(derham) == [1, 1, 0, 0]), "disc x circle: Betti numbers 1 1 0 0")
    call make(ring, 2, 2, 1)
    call check_condition(status == 0 .and. dd_zero(derham) .and. &
        all(betti(derham) == [1, 2, 1, 0]), "ring x circle: Betti numbers 1 2 1 0")

    do level = 1, 2
        call make(disc, 3, 4*level, -1)
        allocate(c0(derham%dimension(0)), c1(derham%dimension(1)), &
            c2(derham%dimension(2)), d1(derham%dimension(1)), d2(derham%dimension(2)))
        call spline_derham_project(derham, 0, f0, c0, status)
        call spline_derham_project(derham, 1, grad_f, d1, status)
        call spline_derham_project(derham, 1, a1, c1, status)
        call spline_derham_project(derham, 2, curl_a, d2, status)
        if (level == 1) then
            call check_condition(maxval(abs(apply(derham%grad, c0) - d1)) < 1.0e-10_dp &
                .and. maxval(abs(apply(derham%curl, c1) - d2)) < 1.0e-10_dp, &
                "2D commuting projections: grad Pi0 = Pi1 grad, curl Pi1 = Pi2 curl")
        end if
        errors(:, level) = projection_errors(c0, c1)
        if (level == 1) call check_condition(cell_basis_ok(c1), &
            "cell basis tables reproduce the evaluated 1-form proxies")
        deallocate(c0, c1, c2, d1, d2)
    end do
    print '(a,2es11.3,a,2es11.3)', "Pi0 f errors", errors(1, :), "  Pi1 A errors", &
        errors(2, :)
    call check_condition(errors(1, 2) < errors(1, 1)/6.0_dp .and. &
        errors(2, 2) < errors(2, 1)/6.0_dp, &
        "projections of 0- and 1-forms converge at order p in the maximum norm")

    call layout_ogrid([2.0_dp, 0.0_dp], 0.5_dp, 0.5_dp, disc, status, polygon=.true.)
    call curve_chart_cylindrical(chart, 1)
    call make(disc, 2, 2, 2)
    allocate(c0(derham%dimension(0)), c1(derham%dimension(1)), &
        c2(derham%dimension(2)), c3(derham%dimension(3)), d1(derham%dimension(1)), &
        d2(derham%dimension(2)), d3(derham%dimension(3)))
    call spline_derham_project(derham, 0, f3d, c0, status, 12)
    call spline_derham_project(derham, 1, grad_f3d, d1, status, 12)
    call spline_derham_project(derham, 1, a3d, c1, status, 12)
    call spline_derham_project(derham, 2, curl_a3d, d2, status, 12)
    call spline_derham_project(derham, 2, b3d, c2, status, 12)
    call spline_derham_project(derham, 3, div_b3d, d3, status, 12)
    print '(a,3es11.3)', "3D commuting defects", maxval(abs(apply(derham%grad, c0) - &
        d1)), maxval(abs(apply(derham%curl, c1) - d2)), &
        maxval(abs(apply(derham%div, c2) - d3))
    call check_condition(maxval(abs(apply(derham%grad, c0) - d1)) < 1.0e-9_dp .and. &
        maxval(abs(apply(derham%curl, c1) - d2)) < 1.0e-9_dp .and. &
        maxval(abs(apply(derham%div, c2) - d3)) < 1.0e-9_dp, &
        "3D commuting projections for grad, curl and div")

    call spline_derham_boundary_dofs(derham, LAYOUT_WALL, 0, nodes)
    call spline_derham_boundary_dofs(derham, LAYOUT_WALL, 1, edges)
    call spline_derham_boundary_dofs(derham, LAYOUT_WALL, 2, faces)
    call check_condition(size(edges) == 2*size(nodes) .and. &
        size(faces) == size(nodes), &
        "boundary traces: closed wall loop, tangential edges and normal fluxes")
    call check_summary("spline_derham")

contains

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

    subroutine make(layout, degree, cells, modes)
        type(block_layout_t), intent(in) :: layout
        integer, intent(in) :: degree, cells, modes

        type(spline_space_options_t) :: options

        options%degree = degree
        allocate(options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = cells
        call spline_space_build(layout, options, space, status)
        if (status /= 0) error stop "space"
        if (modes >= 0) then
            call spline_derham_build(space, derham, status, modes)
        else
            call spline_derham_build(space, derham, status)
        end if
    end subroutine make

    function apply(a, x) result(y)
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
    end function apply

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

    logical function dd_zero(dr) result(ok)
        type(spline_derham_t), intent(in) :: dr
        real(dp), allocatable :: g(:, :), c(:, :), d(:, :)

        g = dense(dr%grad, dr%dimension(0))
        c = dense(dr%curl, dr%dimension(1))
        ok = maxval(abs(matmul(c, g))) == 0.0_dp
        if (dr%toroidal) then
            d = dense(dr%div, dr%dimension(2))
            ok = ok .and. maxval(abs(matmul(d, c))) == 0.0_dp
        end if
    end function dd_zero

    function betti(dr) result(b)
        type(spline_derham_t), intent(in) :: dr
        integer :: b(0:3), rg, rc, rd

        rg = rank(dense(dr%grad, dr%dimension(0)))
        rc = rank(dense(dr%curl, dr%dimension(1)))
        rd = 0
        if (dr%toroidal) rd = rank(dense(dr%div, dr%dimension(2)))
        b = [dr%dimension(0) - rg, dr%dimension(1) - rg - rc, &
            dr%dimension(2) - rc - rd, dr%dimension(3) - rd]
    end function betti

    logical function cell_basis_ok(c1) result(ok)
        real(dp), intent(in) :: c1(:)
        type(spline_quadrature_t) :: quad
        real(dp), allocatable :: values(:, :, :)
        integer, allocatable :: dofs(:), signs(:)
        real(dp) :: v(2), w(2)
        integer :: cell, q, n, status, point
        real(dp) :: worst

        call spline_quadrature_build(space, 3, quad, status)
        worst = 0.0_dp
        do cell = 1, space%cell_count
            call spline_derham_cell_basis(derham, quad, cell, 1, values, dofs, signs)
            do q = 1, quad%points_per_cell
                point = quad%point_offset(cell) + q - 1
                w = 0.0_dp
                do n = 1, size(dofs)
                    w = w + signs(n)*c1(dofs(n))*values(:, n, q)
                end do
                call spline_derham_evaluate(derham, 1, c1, quad%block(point), &
                    quad%xi(:, point), 0.0_dp, v)
                worst = max(worst, maxval(abs(v - w)))
            end do
        end do
        ok = worst < 1.0e-13_dp
    end function cell_basis_ok

    function projection_errors(c0, c1) result(e)
        !! Largest errors of Pi0 f and of the proxy of Pi1 A at 7 x 7 points
        !! of every block.
        real(dp), intent(in) :: c0(:), c1(:)
        real(dp) :: e(2), v(2), w(2), xi(2)
        integer :: k, a, b

        e = 0.0_dp
        do k = 1, space%layout%block_count
            do b = 0, 6
                do a = 0, 6
                    xi = [a, b]/6.0_dp
                    call spline_derham_evaluate(derham, 0, c0, k, xi, 0.0_dp, v(1:1))
                    call f0(k, xi, 0.0_dp, w(1:1))
                    e(1) = max(e(1), abs(v(1) - w(1)))
                    call spline_derham_evaluate(derham, 1, c1, k, xi, 0.0_dp, v)
                    call a1(k, xi, 0.0_dp, w)
                    e(2) = max(e(2), maxval(abs(v - w)))
                end do
            end do
        end do
    end function projection_errors

end program test_spline_derham
