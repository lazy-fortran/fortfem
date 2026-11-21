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
        spline_derham_mass, spline_derham_evaluate, spline_derham_boundary_dofs
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

    type(block_layout_t) :: disc, ring
    integer :: status

    call layout_ogrid([0.0_dp, 0.0_dp], 1.0_dp, 0.5_dp, disc, status, polygon=.true.)
    call ring_layout(ring)
    call check_condition(cohomology_ok(disc, -1, 0), "disc: tree and no loop")
    call check_condition(cohomology_ok(ring, -1, 1), "ring: one closed non-exact form")
    call check_condition(cohomology_ok(disc, 1, 1), "disc x circle: the toroidal class")
    call check_condition(cohomology_ok(ring, 1, 2), &
        "ring x circle: toroidal and poloidal classes")
    call check_condition(surface_ok(ring), &
        "boundary torus: surface curl of gradients vanishes, two surface classes")
    call check_condition(mass_ok(ring, -1), "section mass matrices of all forms")
    call check_condition(mass_ok(ring, 1), "3D mass matrices of all forms")
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
