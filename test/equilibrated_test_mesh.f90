module equilibrated_test_mesh
    !! Structured unit-square meshes of isoparametric Lagrange P_m triangles
    !! for the equilibrated-flux tests. Lattice points (i, j)/(m n) are mapped
    !! by a smooth deformation that fixes the boundary, so the domain stays the
    !! unit square while interior cells and edges are curved. Geometry nodes
    !! are given both in the basis order of triangle_lagrange_basis_t (own
    !! evaluation, independent of the estimator) and in the estimator's
    !! documented order (vertices, edge nodes from local vertex e to e + 1,
    !! centroid).
    use fortfem_kinds, only: dp
    use fortfem_triangle_rt_arbitrary_order, only: triangle_rt_basis_t, &
        initialize_triangle_raviart_thomas, evaluate_triangle_raviart_thomas
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        triangle_lagrange_basis_t, initialize_triangle_lagrange_basis, &
        evaluate_triangle_lagrange_basis, triangle_lagrange_nodes
    implicit none
    private

    public :: test_mesh_t, build_test_mesh, cell_point, galerkin_solve
    public :: random_rt_flux, edge_flux, edge_jump

    real(dp), parameter :: pi = acos(-1.0_dp)

    type :: test_mesh_t
        integer :: m = 1, n = 0, nv = 0, nt = 0, nb = 0, ndof = 0
        real(dp), allocatable :: vertices(:, :)      !< (2, nv)
        integer, allocatable :: triangles(:, :)      !< (3, nt)
        real(dp), allocatable :: nodes(:, :, :)      !< estimator order (2, nloc, nt)
        real(dp), allocatable :: xb(:, :, :)         !< basis order (2, nb, nt)
        integer, allocatable :: gdof(:, :)           !< lattice dof (nb, nt)
        real(dp), allocatable :: dof_position(:, :)  !< (2, ndof), mapped
        type(triangle_lagrange_basis_t) :: basis
    end type test_mesh_t

contains

    pure function deform(x, y, amplitude) result(p)
        real(dp), intent(in) :: x, y, amplitude
        real(dp) :: p(2)

        p(1) = x + amplitude*sin(pi*x)*sin(2*pi*y)
        p(2) = y + amplitude*sin(2*pi*x)*sin(pi*y)
    end function deform

    subroutine build_test_mesh(n, m, amplitude, mesh, status)
        integer, intent(in) :: n, m
        real(dp), intent(in) :: amplitude
        type(test_mesh_t), intent(out) :: mesh
        integer, intent(out) :: status
        real(dp), allocatable :: ref(:, :)
        integer, allocatable :: est(:, :)
        integer :: i, j, v, c, t, k, nloc, lat(2), p(2, 3), a, b

        call initialize_triangle_lagrange_basis(m, mesh%basis, status)
        if (status /= 0) return
        call triangle_lagrange_nodes(mesh%basis, ref)
        mesh%m = m
        mesh%n = n
        mesh%nb = size(ref, 2)
        nloc = mesh%nb
        mesh%nv = (n + 1)**2
        mesh%nt = 2*n*n
        mesh%ndof = (m*n + 1)**2
        allocate (mesh%vertices(2, mesh%nv), mesh%triangles(3, mesh%nt))
        allocate (mesh%nodes(2, nloc, mesh%nt), mesh%xb(2, nloc, mesh%nt))
        allocate (mesh%gdof(nloc, mesh%nt), mesh%dof_position(2, mesh%ndof))
        do j = 0, m*n
            do i = 0, m*n
                mesh%dof_position(:, j*(m*n + 1) + i + 1) = deform( &
                    real(i, dp)/(m*n), real(j, dp)/(m*n), amplitude)
            end do
        end do
        do j = 0, n
            do i = 0, n
                mesh%vertices(:, j*(n + 1) + i + 1) = &
                    mesh%dof_position(:, (m*j)*(m*n + 1) + m*i + 1)
            end do
        end do
        c = 0
        do j = 0, n - 1
            do i = 0, n - 1
                v = j*(n + 1) + i + 1
                if (mod(i + j, 2) == 0) then
                    mesh%triangles(:, c + 1) = [v, v + 1, v + n + 2]
                    mesh%triangles(:, c + 2) = [v, v + n + 2, v + n + 1]
                else
                    mesh%triangles(:, c + 1) = [v, v + 1, v + n + 1]
                    mesh%triangles(:, c + 2) = [v + 1, v + n + 2, v + n + 1]
                end if
                c = c + 2
            end do
        end do
        call estimator_lattice_offsets(m, est)
        do t = 1, mesh%nt
            do k = 1, 3
                v = mesh%triangles(k, t)
                p(:, k) = [mod(v - 1, n + 1), (v - 1)/(n + 1)]
            end do
            ! Estimator order: integer reference coordinates (a, b) in units 1/m.
            do k = 1, nloc
                a = est(1, k)
                b = est(2, k)
                lat = m*p(:, 1) + a*(p(:, 2) - p(:, 1)) + b*(p(:, 3) - p(:, 1))
                mesh%nodes(:, k, t) = &
                    mesh%dof_position(:, lat(2)*(m*n + 1) + lat(1) + 1)
            end do
            ! Basis order: reference positions of the Lagrange basis nodes.
            do k = 1, nloc
                a = nint(m*ref(1, k))
                b = nint(m*ref(2, k))
                lat = m*p(:, 1) + a*(p(:, 2) - p(:, 1)) + b*(p(:, 3) - p(:, 1))
                mesh%gdof(k, t) = lat(2)*(m*n + 1) + lat(1) + 1
                mesh%xb(:, k, t) = mesh%dof_position(:, mesh%gdof(k, t))
            end do
        end do
        status = 0
    end subroutine build_test_mesh

    subroutine estimator_lattice_offsets(m, est)
        !! Integer reference coordinates (units 1/m) of the estimator's
        !! geometry-node order.
        integer, intent(in) :: m
        integer, allocatable, intent(out) :: est(:, :)

        select case (m)
        case (1)
            est = reshape([0, 0, 1, 0, 0, 1], [2, 3])
        case (2)
            est = reshape([0, 0, 2, 0, 0, 2, 1, 0, 1, 1, 0, 1], [2, 6])
        case default
            est = reshape([0, 0, 3, 0, 0, 3, 1, 0, 2, 0, 2, 1, 1, 2, 0, 2, 0, 1, &
                1, 1], [2, 10])
        end select
    end subroutine estimator_lattice_offsets

    subroutine cell_point(mesh, t, xi, eta, x, jac, det, phi, grad_ref)
        !! Position, Jacobian, determinant and Lagrange values / reference
        !! gradients of cell t at the reference point (xi, eta).
        type(test_mesh_t), intent(in) :: mesh
        integer, intent(in) :: t
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: x(2), jac(2, 2), det
        real(dp), intent(out) :: phi(:), grad_ref(:, :)
        integer :: status

        call evaluate_triangle_lagrange_basis(mesh%basis, xi, eta, phi, grad_ref, &
            status)
        if (status /= 0) error stop "lagrange evaluation failed"
        x = matmul(mesh%xb(:, :, t), phi)
        jac = matmul(mesh%xb(:, :, t), transpose(grad_ref))
        det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
    end subroutine cell_point

    subroutine random_rt_flux(mesh, k0, rt, dofs)
        !! Global conforming RT_k0 flux of the mesh from random edge moments (one
        !! array per global edge, flipped on cells traversing it backwards for
        !! even Legendre degrees, i.e. odd 1-based moments) and random interior
        !! moments. dofs(:, t) are the local coefficients of the Piola-mapped
        !! reference basis.
        type(test_mesh_t), intent(in) :: mesh
        integer, intent(in) :: k0
        type(triangle_rt_basis_t), intent(out) :: rt
        real(dp), allocatable, intent(out) :: dofs(:, :)
        real(dp), allocatable :: edge_val(:, :)
        integer, allocatable :: edge_id(:, :)
        integer :: nedge0, nint0, nt, nedges, t, e, a, b, j, i, stat, sign

        call initialize_triangle_raviart_thomas(k0, rt, stat)
        if (stat /= 0) error stop "rt"
        nedge0 = k0 + 1
        nint0 = k0*(k0 + 1)
        nt = mesh%nt
        allocate (edge_id(mesh%nv, mesh%nv), dofs((k0 + 1)*(k0 + 3), nt))
        edge_id = 0
        nedges = 0
        do t = 1, nt
            do e = 1, 3
                a = mesh%triangles(e, t)
                b = mesh%triangles(mod(e, 3) + 1, t)
                if (edge_id(min(a, b), max(a, b)) == 0) then
                    nedges = nedges + 1
                    edge_id(min(a, b), max(a, b)) = nedges
                end if
            end do
        end do
        allocate (edge_val(nedge0, nedges))
        do e = 1, nedges
            do j = 1, nedge0
                edge_val(j, e) = sin(2.3_dp*e + 1.1_dp*j)
            end do
        end do
        do t = 1, nt
            do e = 1, 3
                a = mesh%triangles(e, t)
                b = mesh%triangles(mod(e, 3) + 1, t)
                do j = 1, nedge0
                    sign = 1
                    if (a > b .and. mod(j, 2) == 1) sign = -1
                    dofs((e - 1)*nedge0 + j, t) = sign &
                        *edge_val(j, edge_id(min(a, b), max(a, b)))
                end do
            end do
            do i = 1, nint0
                dofs(3*nedge0 + i, t) = 0.5_dp*sin(0.7_dp*t + 1.3_dp*i)
            end do
        end do
    end subroutine random_rt_flux

    subroutine galerkin_solve(mesh, xi, eta, w, coef, src, fixed, extra, u, info)
        !! Dense Galerkin solve of int A grad u . grad v = int f v + extra(v) on
        !! the isoparametric P_m space; dofs flagged fixed are set to zero.
        type(test_mesh_t), intent(in) :: mesh
        real(dp), intent(in) :: xi(:), eta(:), w(:), coef(:, :), src(:, :)
        logical, intent(in) :: fixed(:)
        real(dp), intent(in) :: extra(:)
        real(dp), intent(out) :: u(:)
        integer, intent(out) :: info
        real(dp), allocatable :: k(:, :), rhs(:, :), phi(:), gref(:, :)
        real(dp) :: x(2), jac(2, 2), det, inv_t(2, 2)
        integer, allocatable :: ipiv(:)
        integer :: t, q, i, j, gi, gj
        interface
            subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                import :: dp
                integer, intent(in) :: n, nrhs, lda, ldb
                integer, intent(out) :: ipiv(*), info
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            end subroutine dgesv
        end interface
        real(dp), allocatable :: grad(:, :)

        allocate (k(mesh%ndof, mesh%ndof), rhs(mesh%ndof, 1), ipiv(mesh%ndof))
        allocate (phi(mesh%nb), gref(2, mesh%nb), grad(2, mesh%nb))
        k = 0.0_dp
        rhs(:, 1) = extra
        do t = 1, mesh%nt
            do q = 1, size(w)
                call cell_point(mesh, t, xi(q), eta(q), x, jac, det, phi, gref)
                inv_t(1, :) = [jac(2, 2), -jac(2, 1)]/det
                inv_t(2, :) = [-jac(1, 2), jac(1, 1)]/det
                grad = matmul(inv_t, gref)
                do i = 1, mesh%nb
                    gi = mesh%gdof(i, t)
                    rhs(gi, 1) = rhs(gi, 1) + w(q)*det*src(q, t)*phi(i)
                    do j = 1, mesh%nb
                        gj = mesh%gdof(j, t)
                        k(gi, gj) = k(gi, gj) + w(q)*det*coef(q, t) &
                            *dot_product(grad(:, i), grad(:, j))
                    end do
                end do
            end do
        end do
        do i = 1, mesh%ndof
            if (fixed(i)) then
                k(i, :) = 0.0_dp
                k(:, i) = 0.0_dp
                k(i, i) = 1.0_dp
                rhs(i, 1) = 0.0_dp
            end if
        end do
        call dgesv(mesh%ndof, 1, k, mesh%ndof, ipiv, rhs, mesh%ndof, info)
        u = rhs(:, 1)
    end subroutine galerkin_solve


    real(dp) function edge_jump(mesh, rt, dofs, t1, e1, t2, e2) result(jmp)
        !! Max over sample points of q_1(s) + q_2(1 - s), q the flux per unit
        !! edge parameter along the outward normal (t_y, -t_x).
        type(test_mesh_t), intent(in) :: mesh
        type(triangle_rt_basis_t), intent(in) :: rt
        real(dp), intent(in) :: dofs(:, :)
        integer, intent(in) :: t1, e1, t2, e2
        real(dp), parameter :: s(3) = [0.2_dp, 0.5_dp, 0.9_dp]
        integer :: l

        jmp = 0.0_dp
        do l = 1, 3
            jmp = max(jmp, abs(edge_flux(mesh, rt, dofs, t1, e1, s(l)) &
                + edge_flux(mesh, rt, dofs, t2, e2, 1.0_dp - s(l))))
        end do
    end function edge_jump

    real(dp) function edge_flux(mesh, rt, dofs, tc, ec, s) result(qq)
        type(test_mesh_t), intent(in) :: mesh
        type(triangle_rt_basis_t), intent(in) :: rt
        real(dp), intent(in) :: dofs(:, :), s
        integer, intent(in) :: tc, ec
        real(dp) :: p(2), tref(2), tphys(2), xx(2), jj(2, 2), dd, sg(2)
        real(dp) :: phi(mesh%nb), gref(2, mesh%nb), rt_val(2, size(dofs, 1))
        real(dp) :: rt_div(size(dofs, 1))
        integer :: stat

        select case (ec)
        case (1)
            p = [s, 0.0_dp]
            tref = [1.0_dp, 0.0_dp]
        case (2)
            p = [1.0_dp - s, s]
            tref = [-1.0_dp, 1.0_dp]
        case default
            p = [0.0_dp, 1.0_dp - s]
            tref = [0.0_dp, -1.0_dp]
        end select
        call cell_point(mesh, tc, p(1), p(2), xx, jj, dd, phi, gref)
        call evaluate_triangle_raviart_thomas(rt, p(1), p(2), rt_val, rt_div, stat)
        sg = matmul(jj, matmul(rt_val, dofs(:, tc)))/dd
        tphys = matmul(jj, tref)
        qq = sg(1)*tphys(2) - sg(2)*tphys(1)
    end function edge_flux

end module equilibrated_test_mesh
