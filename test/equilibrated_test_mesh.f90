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
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        triangle_lagrange_basis_t, initialize_triangle_lagrange_basis, &
        evaluate_triangle_lagrange_basis, triangle_lagrange_nodes
    implicit none
    private

    public :: test_mesh_t, build_test_mesh, cell_point, galerkin_solve

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

end module equilibrated_test_mesh
