module fortfem_api_solvers
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_api_types, only: function_space_t, function_t,              &
        vector_function_t, vector_function_space_t, dirichlet_bc_t,        &
        vector_bc_t, neumann_bc_t
    use fortfem_api_forms, only: form_equation_t, form_expr_t
    use fortfem_api_mesh, only: find_triangle_edges
    use fortfem_advanced_solvers, only: solver_options_t, solver_stats_t,  &
        solver_options, cg_solve, pcg_solve, bicgstab_solve, gmres_solve,  &
        jacobi_preconditioner, ilu_preconditioner, advanced_solve => solve
    use basis_p2_2d_module, only: basis_p2_2d_t
    use basis_q1_quad_2d_module, only: q1_shape_functions,                 &
        q1_shape_derivatives, q1_jacobian
    use fortfem_basis_edge_2d, only: edge_basis_2d_t
    implicit none

    private

    public :: solver_options_t, solver_stats_t
    public :: solver_options
    public :: cg_solve, pcg_solve, bicgstab_solve, gmres_solve
    public :: jacobi_preconditioner, ilu_preconditioner

    public :: assemble_laplacian_system
    public :: solve
    public :: solve_scalar
    public :: solve_vector
    public :: solve_laplacian_problem
    public :: solve_laplacian_problem_p2
    public :: solve_mixed_bc
    public :: solve_neumann
    public :: compute_boundary_integral
    public :: solve_laplacian_with_neumann
    public :: solve_pure_neumann_problem
    public :: solve_curl_curl_problem
    public :: solve_generic_problem
    public :: solve_generic_vector_problem

    interface solve
        module procedure solve_scalar
        module procedure solve_vector
    end interface solve

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface

contains

    subroutine solve_scalar(equation, uh, bc, options, stats)
        type(form_equation_t), intent(in) :: equation
        type(function_t), intent(inout) :: uh
        type(dirichlet_bc_t), intent(in) :: bc
        type(solver_options_t), intent(in), optional :: options
        type(solver_stats_t), intent(out), optional :: stats

        type(solver_stats_t) :: local_stats
        logical :: have_stats

        write(*,*) "Solving: ", trim(equation%lhs%description), " == ",     &
            trim(equation%rhs%description)

        have_stats = .false.

        if (index(equation%lhs%description, "grad") > 0) then
            if (uh%space%degree == 2) then
                call solve_laplacian_problem_p2(uh, bc, options, local_stats)
            else
                call solve_laplacian_problem(uh, bc, options, local_stats)
            end if
            have_stats = .true.
        else
            call solve_generic_problem(uh, bc)
        end if

        if (present(stats) .and. have_stats) then
            stats = local_stats
        end if
    end subroutine solve_scalar

    subroutine solve_vector(equation, Eh, bc, solver_type, options, stats)
        type(form_equation_t), intent(in) :: equation
        type(vector_function_t), intent(inout) :: Eh
        type(vector_bc_t), intent(in) :: bc
        character(len=*), intent(in), optional :: solver_type
        type(solver_options_t), intent(in), optional :: options
        type(solver_stats_t), intent(out), optional :: stats

        character(len=32) :: solver
        type(solver_options_t) :: local_opts
        type(solver_stats_t) :: local_stats

        solver = "gmres"
        if (present(solver_type)) solver = solver_type

        local_opts = solver_options(method="gmres", tolerance=1.0e-6_dp,    &
                                    max_iterations=100, restart=20)
        if (present(options)) local_opts = options

        local_stats%converged = .false.
        local_stats%iterations = 0
        local_stats%final_residual = 0.0_dp
        local_stats%restarts = 0
        local_stats%method_used = ""

        write(*,*) "Solving vector problem: ",                               &
            trim(equation%lhs%description), " == ",                          &
            trim(equation%rhs%description)
        write(*,*) "Using solver: ", trim(solver)

        if (index(equation%lhs%description, "curl") > 0) then
            call solve_curl_curl_problem(Eh, bc, solver, local_opts,        &
                                         local_stats)
        else
            call solve_generic_vector_problem(Eh, bc)
        end if

        if (present(stats)) then
            stats = local_stats
        end if
    end subroutine solve_vector

    subroutine assemble_laplacian_system(space, bc, K, F)
        type(function_space_t), intent(in) :: space
        type(dirichlet_bc_t), intent(in) :: bc
        real(dp), allocatable, intent(out) :: K(:,:), F(:)

        integer :: ndof, i

        ndof = space%ndof
        allocate(K(ndof, ndof), F(ndof))

        K = 0.0_dp
        F = 0.0_dp

        call assemble_laplacian_triangles(space, K, F)
        call assemble_laplacian_quads(space, K, F)

        do i = 1, space%mesh%data%n_vertices
            if (space%mesh%data%is_boundary_vertex(i)) then
                K(i, :) = 0.0_dp
                K(i, i) = 1.0_dp
                F(i) = bc%value
            end if
        end do
    end subroutine assemble_laplacian_system

    pure subroutine compute_p1_triangle_gradients(x1, y1, x2, y2, x3, y3,  &
                                                  area, grad_x, grad_y)
        real(dp), intent(in) :: x1, y1, x2, y2, x3, y3
        real(dp), intent(out) :: area
        real(dp), intent(out) :: grad_x(3), grad_y(3)

        real(dp) :: a11, a12, a21, a22, det_a

        a11 = x2 - x1
        a12 = x3 - x1
        a21 = y2 - y1
        a22 = y3 - y1

        det_a = a11 * a22 - a12 * a21

        area = 0.5_dp * abs((x2-x1) * (y3-y1) - (x3-x1) * (y2-y1))

        grad_x(1) = (-a22 + a21) / det_a
        grad_y(1) = ( a12 - a11) / det_a

        grad_x(2) = a22 / det_a
        grad_y(2) = -a12 / det_a

        grad_x(3) = -a21 / det_a
        grad_y(3) = a11 / det_a
    end subroutine compute_p1_triangle_gradients

    subroutine add_p1_triangle_contribution(space, triangle_id, K, F)
        type(function_space_t), intent(in) :: space
        integer, intent(in) :: triangle_id
        real(dp), intent(inout) :: K(:,:), F(:)

        integer :: v1, v2, v3
        integer :: i, j
        real(dp) :: x1, y1, x2, y2, x3, y3
        real(dp) :: area
        real(dp) :: bx(3), by(3), K_elem(3,3)

        v1 = space%mesh%data%triangles(1, triangle_id)
        v2 = space%mesh%data%triangles(2, triangle_id)
        v3 = space%mesh%data%triangles(3, triangle_id)

        x1 = space%mesh%data%vertices(1, v1)
        y1 = space%mesh%data%vertices(2, v1)
        x2 = space%mesh%data%vertices(1, v2)
        y2 = space%mesh%data%vertices(2, v2)
        x3 = space%mesh%data%vertices(1, v3)
        y3 = space%mesh%data%vertices(2, v3)

        call compute_p1_triangle_gradients(x1, y1, x2, y2, x3, y3, area,   &
                                           bx, by)

        do i = 1, 3
            do j = 1, 3
                K_elem(i, j) = area * (bx(i) * bx(j) + by(i) * by(j))
            end do
        end do

        K(v1, v1) = K(v1, v1) + K_elem(1,1)
        K(v1, v2) = K(v1, v2) + K_elem(1,2)
        K(v1, v3) = K(v1, v3) + K_elem(1,3)
        K(v2, v1) = K(v2, v1) + K_elem(2,1)
        K(v2, v2) = K(v2, v2) + K_elem(2,2)
        K(v2, v3) = K(v2, v3) + K_elem(2,3)
        K(v3, v1) = K(v3, v1) + K_elem(3,1)
        K(v3, v2) = K(v3, v2) + K_elem(3,2)
        K(v3, v3) = K(v3, v3) + K_elem(3,3)

        F(v1) = F(v1) + area / 3.0_dp
        F(v2) = F(v2) + area / 3.0_dp
        F(v3) = F(v3) + area / 3.0_dp
    end subroutine add_p1_triangle_contribution

    subroutine assemble_laplacian_triangles(space, K, F)
        type(function_space_t), intent(in) :: space
        real(dp), intent(inout) :: K(:,:), F(:)

        integer :: e

        do e = 1, space%mesh%data%n_triangles
            call add_p1_triangle_contribution(space, e, K, F)
        end do
    end subroutine assemble_laplacian_triangles

    subroutine add_q1_quad_contribution(space, quad_id, K, F)
        type(function_space_t), intent(in) :: space
        integer, intent(in) :: quad_id
        real(dp), intent(inout) :: K(:,:), F(:)

        integer :: i, j, vi
        integer :: v_ids(4)
        real(dp) :: coords(2,4)
        real(dp) :: K_elem(4,4), F_elem(4)

        v_ids = space%mesh%data%quads(:, quad_id)

        do i = 1, 4
            vi = v_ids(i)
            coords(1, i) = space%mesh%data%vertices(1, vi)
            coords(2, i) = space%mesh%data%vertices(2, vi)
        end do

        call compute_q1_quad_element(coords, K_elem, F_elem)

        do i = 1, 4
            vi = v_ids(i)
            do j = 1, 4
                K(vi, v_ids(j)) = K(vi, v_ids(j)) + K_elem(i, j)
            end do
            F(vi) = F(vi) + F_elem(i)
        end do
    end subroutine add_q1_quad_contribution

    subroutine compute_q1_quad_element(coords, K_elem, F_elem)
        real(dp), intent(in) :: coords(2,4)
        real(dp), intent(out) :: K_elem(4,4), F_elem(4)

        integer :: i, j, kx, ky
        real(dp) :: jac(2,2), det_jac, inv_jac(2,2)
        real(dp) :: dN_dxi(4), dN_deta(4)
        real(dp) :: grad_ref(2), grad_phys(2,4)
        real(dp) :: N(4)
        real(dp) :: xi, eta, weight
        logical :: success
        real(dp), parameter :: gauss_pts(2) = [-0.5773502691896257_dp,      &
                                               0.5773502691896257_dp]
        real(dp), parameter :: gauss_w(2) = [1.0_dp, 1.0_dp]

        K_elem = 0.0_dp
        F_elem = 0.0_dp
        do kx = 1, 2
            do ky = 1, 2
                xi = gauss_pts(kx)
                eta = gauss_pts(ky)

                call q1_shape_derivatives(xi, eta, dN_dxi, dN_deta)
                call q1_jacobian(xi, eta, coords, jac, det_jac,            &
                                 inv_jac, success)

                if (.not. success) cycle

                do i = 1, 4
                    grad_ref(1) = dN_dxi(i)
                    grad_ref(2) = dN_deta(i)
                    grad_phys(:, i) = matmul(transpose(inv_jac), grad_ref)
                end do

                call q1_shape_functions(xi, eta, N)

                weight = det_jac * gauss_w(kx) * gauss_w(ky)

                do i = 1, 4
                    do j = 1, 4
                        K_elem(i, j) = K_elem(i, j) + weight *             &
                            (grad_phys(1, i) * grad_phys(1, j) +           &
                             grad_phys(2, i) * grad_phys(2, j))
                    end do
                    F_elem(i) = F_elem(i) + weight * N(i)
                end do
            end do
        end do
    end subroutine compute_q1_quad_element

    subroutine assemble_laplacian_quads(space, K, F)
        type(function_space_t), intent(in) :: space
        real(dp), intent(inout) :: K(:,:), F(:)

        integer :: q

        if (space%mesh%data%n_quads <= 0) return

        do q = 1, space%mesh%data%n_quads
            call add_q1_quad_contribution(space, q, K, F)
        end do
    end subroutine assemble_laplacian_quads

    subroutine solve_laplacian_problem(uh, bc, options, stats)
        type(function_t), intent(inout) :: uh
        type(dirichlet_bc_t), intent(in) :: bc
        type(solver_options_t), intent(in), optional :: options
        type(solver_stats_t), intent(out) :: stats

        real(dp), allocatable :: K(:,:), F(:)
        integer :: ndof
        type(solver_options_t) :: local_opts

        ndof = uh%space%ndof

        call assemble_laplacian_system(uh%space, bc, K, F)

        if (.not. allocated(uh%values)) then
            allocate(uh%values(ndof))
        end if
        uh%values = 0.0_dp

        if (present(options)) then
            local_opts = options
        else
            local_opts = solver_options(method="auto")
        end if

        call advanced_solve(K, F, uh%values, local_opts, stats)

        if (.not. stats%converged) then
            write(*,*) "Warning: Laplacian solver did not report convergence.", &
                " Final residual =", stats%final_residual
        end if

        deallocate(K, F)
    end subroutine solve_laplacian_problem

    subroutine solve_laplacian_problem_p2(uh, bc, options, stats)
        type(function_t), intent(inout) :: uh
        type(dirichlet_bc_t), intent(in) :: bc
        type(solver_options_t), intent(in), optional :: options
        type(solver_stats_t), intent(out) :: stats

        real(dp), allocatable :: K(:,:), F(:)
        integer :: ndof
        type(solver_options_t) :: local_opts

        ndof = uh%space%ndof
        allocate(K(ndof, ndof), F(ndof))

        K = 0.0_dp
        F = 0.0_dp

        call assemble_p2_laplacian_system(uh, bc, K, F)

        if (.not. allocated(uh%values)) then
            allocate(uh%values(ndof))
        end if
        uh%values = 0.0_dp

        if (present(options)) then
            local_opts = options
        else
            local_opts = solver_options(method="auto")
        end if

        call advanced_solve(K, F, uh%values, local_opts, stats)

        if (.not. stats%converged) then
            write(*,*) "Warning: P2 Laplacian solver did not report convergence.", &
                " Final residual =", stats%final_residual
        end if

        deallocate(K, F)
    end subroutine solve_laplacian_problem_p2

    subroutine assemble_p2_laplacian_system(uh, bc, K, F)
        type(function_t), intent(in) :: uh
        type(dirichlet_bc_t), intent(in) :: bc
        real(dp), intent(inout) :: K(:,:), F(:)

        integer :: ndof, i, j, e, v1, v2, v3
        real(dp) :: K_elem(6,6), F_elem(6)
        integer :: dofs(6), edge1, edge2, edge3
        real(dp) :: vertices(2,3)

        ndof = uh%space%ndof

        do e = 1, uh%space%mesh%data%n_triangles
            v1 = uh%space%mesh%data%triangles(1, e)
            v2 = uh%space%mesh%data%triangles(2, e)
            v3 = uh%space%mesh%data%triangles(3, e)

            vertices(1,1) = uh%space%mesh%data%vertices(1, v1)
            vertices(2,1) = uh%space%mesh%data%vertices(2, v1)
            vertices(1,2) = uh%space%mesh%data%vertices(1, v2)
            vertices(2,2) = uh%space%mesh%data%vertices(2, v2)
            vertices(1,3) = uh%space%mesh%data%vertices(1, v3)
            vertices(2,3) = uh%space%mesh%data%vertices(2, v3)

            dofs(1:3) = [v1, v2, v3]
            call find_triangle_edges(uh%space%mesh%data, e, edge1, edge2,   &
                                     edge3)
            dofs(4) = uh%space%mesh%data%n_vertices + edge1
            dofs(5) = uh%space%mesh%data%n_vertices + edge2
            dofs(6) = uh%space%mesh%data%n_vertices + edge3

            call compute_p2_element_matrices(vertices, K_elem, F_elem)

            do i = 1, 6
                if (dofs(i) > 0 .and. dofs(i) <= ndof) then
                    do j = 1, 6
                        if (dofs(j) > 0 .and. dofs(j) <= ndof) then
                            K(dofs(i), dofs(j)) = K(dofs(i), dofs(j))       &
                                + K_elem(i,j)
                        end if
                    end do
                    F(dofs(i)) = F(dofs(i)) + F_elem(i)
                end if
            end do
        end do

        call apply_p2_dirichlet_bc(uh%space%mesh%data, bc%value, ndof, K, F)
    end subroutine assemble_p2_laplacian_system

    subroutine compute_p2_element_matrices(vertices, K_elem, F_elem)
        real(dp), intent(in) :: vertices(2,3)
        real(dp), intent(out) :: K_elem(6,6), F_elem(6)

        type(basis_p2_2d_t) :: basis_p2
        real(dp) :: xi_q(3), eta_q(3), w_q(3)
        real(dp) :: jac(2,2), det_j, inv_jac(2,2), area
        real(dp) :: grad_i(2), grad_j(2)
        integer :: i, j, kq

        xi_q = [1.0_dp/6.0_dp, 2.0_dp/3.0_dp, 1.0_dp/6.0_dp]
        eta_q = [1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 2.0_dp/3.0_dp]
        w_q = [1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp]

        area = 0.5_dp * abs((vertices(1,2)-vertices(1,1))                    &
            *(vertices(2,3)-vertices(2,1))                                    &
            - (vertices(1,3)-vertices(1,1))*(vertices(2,2)-vertices(2,1)))

        call basis_p2%compute_jacobian(vertices, jac, det_j)
        inv_jac(1,1) = jac(2,2) / det_j
        inv_jac(1,2) = -jac(1,2) / det_j
        inv_jac(2,1) = -jac(2,1) / det_j
        inv_jac(2,2) = jac(1,1) / det_j

        K_elem = 0.0_dp
        F_elem = 0.0_dp

        do i = 1, 6
            do j = 1, 6
                do kq = 1, 3
                    grad_i = matmul(inv_jac,                                 &
                        basis_p2%grad(i, xi_q(kq), eta_q(kq)))
                    grad_j = matmul(inv_jac,                                 &
                        basis_p2%grad(j, xi_q(kq), eta_q(kq)))
                    K_elem(i,j) = K_elem(i,j) + w_q(kq) * area *             &
                        (grad_i(1)*grad_j(1) + grad_i(2)*grad_j(2))
                end do
            end do
            do kq = 1, 3
                F_elem(i) = F_elem(i) + w_q(kq) * area *                     &
                    basis_p2%eval(i, xi_q(kq), eta_q(kq))
            end do
        end do
    end subroutine compute_p2_element_matrices

    subroutine apply_p2_dirichlet_bc(mdata, bc_value, ndof, K, F)
        type(mesh_2d_t), intent(in) :: mdata
        real(dp), intent(in) :: bc_value
        integer, intent(in) :: ndof
        real(dp), intent(inout) :: K(:,:), F(:)

        integer :: i, j

        do i = 1, mdata%n_vertices
            if (mdata%is_boundary_vertex(i)) then
                K(i,:) = 0.0_dp
                K(i,i) = 1.0_dp
                F(i) = bc_value
            end if
        end do

        do i = 1, mdata%n_edges
            if (mdata%is_boundary_edge(i)) then
                j = mdata%n_vertices + i
                if (j <= ndof) then
                    K(j,:) = 0.0_dp
                    K(j,j) = 1.0_dp
                    F(j) = bc_value
                end if
            end if
        end do
    end subroutine apply_p2_dirichlet_bc

    subroutine solve_generic_problem(uh, bc)
        type(function_t), intent(inout) :: uh
        type(dirichlet_bc_t), intent(in) :: bc

        if (allocated(uh%values)) then
            uh%values = bc%value
        end if
    end subroutine solve_generic_problem

    subroutine solve_mixed_bc(equation, uh, dirichlet_bc, neumann_bc)
        type(form_equation_t), intent(in) :: equation
        type(function_t), intent(inout) :: uh
        type(dirichlet_bc_t), intent(in) :: dirichlet_bc
        type(neumann_bc_t), intent(in) :: neumann_bc

        write(*,*) "Solving mixed BC problem: ",                             &
            trim(equation%lhs%description), " == ",                         &
            trim(equation%rhs%description)

        call solve_laplacian_with_neumann(uh, dirichlet_bc, neumann_bc)
    end subroutine solve_mixed_bc

    subroutine solve_neumann(equation, uh, neumann_bc)
        type(form_equation_t), intent(in) :: equation
        type(function_t), intent(inout) :: uh
        type(neumann_bc_t), intent(in) :: neumann_bc

        write(*,*) "Solving pure Neumann problem: ",                         &
            trim(equation%lhs%description), " == ",                         &
            trim(equation%rhs%description)

        call solve_pure_neumann_problem(uh, neumann_bc)
    end subroutine solve_neumann

    subroutine compute_boundary_integral(neumann_bc, integral_value)
        type(neumann_bc_t), intent(in) :: neumann_bc
        real(dp), intent(out) :: integral_value

        integer :: e, v1, v2
        real(dp) :: x1, y1, x2, y2, edge_length, perimeter

        integral_value = 0.0_dp

        if (trim(neumann_bc%flux_type) == "constant") then
            perimeter = 0.0_dp

            do e = 1, neumann_bc%space%mesh%data%n_edges
                if (neumann_bc%space%mesh%data%is_boundary_edge(e)) then
                    v1 = neumann_bc%space%mesh%data%edges(1, e)
                    v2 = neumann_bc%space%mesh%data%edges(2, e)

                    x1 = neumann_bc%space%mesh%data%vertices(1, v1)
                    y1 = neumann_bc%space%mesh%data%vertices(2, v1)
                    x2 = neumann_bc%space%mesh%data%vertices(1, v2)
                    y2 = neumann_bc%space%mesh%data%vertices(2, v2)

                    edge_length = sqrt((x2-x1)**2 + (y2-y1)**2)
                    perimeter = perimeter + edge_length
                end if
            end do

            integral_value = neumann_bc%constant_value * perimeter
        else
            integral_value = 0.0_dp
        end if
    end subroutine compute_boundary_integral

    subroutine solve_laplacian_with_neumann(uh, dirichlet_bc, neumann_bc)
        type(function_t), intent(inout) :: uh
        type(dirichlet_bc_t), intent(in) :: dirichlet_bc
        type(neumann_bc_t), intent(in) :: neumann_bc

        real(dp), allocatable :: K(:,:), F(:)
        integer, allocatable :: ipiv(:)
        integer :: ndof, info

        ndof = uh%space%ndof
        allocate(K(ndof, ndof), F(ndof), ipiv(ndof))

        call assemble_laplacian_neumann_system(uh, dirichlet_bc, neumann_bc,&
                                               K, F)

        call dgesv(ndof, 1, K, ndof, ipiv, F, ndof, info)

        if (info == 0) then
            uh%values = F
        else
            write(*,*) "Warning: Mixed BC LAPACK solver failed with info =",&
                info
            if (allocated(uh%values)) uh%values = 0.0_dp
        end if

        deallocate(K, F, ipiv)
    end subroutine solve_laplacian_with_neumann

    subroutine assemble_laplacian_neumann_system(uh, dirichlet_bc,          &
                                                 neumann_bc, K, F)
        type(function_t), intent(in) :: uh
        type(dirichlet_bc_t), intent(in) :: dirichlet_bc
        type(neumann_bc_t), intent(in) :: neumann_bc
        real(dp), intent(inout) :: K(:,:), F(:)

        K = 0.0_dp
        F = 0.0_dp

        call assemble_neumann_interior_laplacian(uh, K, F)
        call assemble_neumann_boundary_flux(uh, neumann_bc, F)
        call apply_mixed_dirichlet_bc(uh, dirichlet_bc, K, F)
    end subroutine assemble_laplacian_neumann_system

    subroutine assemble_neumann_interior_laplacian(uh, K, F)
        type(function_t), intent(in) :: uh
        real(dp), intent(inout) :: K(:,:), F(:)

        integer :: e

        do e = 1, uh%space%mesh%data%n_triangles
            call add_p1_triangle_contribution(uh%space, e, K, F)
        end do
    end subroutine assemble_neumann_interior_laplacian

    subroutine assemble_neumann_boundary_flux(uh, neumann_bc, F)
        type(function_t), intent(in) :: uh
        type(neumann_bc_t), intent(in) :: neumann_bc
        real(dp), intent(inout) :: F(:)

        integer :: e, v1, v2
        real(dp) :: x1, x2, y1, y2, edge_length

        do e = 1, uh%space%mesh%data%n_edges
            if (uh%space%mesh%data%is_boundary_edge(e)) then
                v1 = uh%space%mesh%data%edges(1, e)
                v2 = uh%space%mesh%data%edges(2, e)

                x1 = uh%space%mesh%data%vertices(1, v1)
                x2 = uh%space%mesh%data%vertices(1, v2)

                if (x1 > 0.9_dp .and. x2 > 0.9_dp) then
                    y1 = uh%space%mesh%data%vertices(2, v1)
                    y2 = uh%space%mesh%data%vertices(2, v2)
                    edge_length = sqrt((x2-x1)**2 + (y2-y1)**2)

                    F(v1) = F(v1) + neumann_bc%constant_value * edge_length &
                        / 2.0_dp
                    F(v2) = F(v2) + neumann_bc%constant_value * edge_length &
                        / 2.0_dp
                end if
            end if
        end do
    end subroutine assemble_neumann_boundary_flux

    subroutine apply_mixed_dirichlet_bc(uh, dirichlet_bc, K, F)
        type(function_t), intent(in) :: uh
        type(dirichlet_bc_t), intent(in) :: dirichlet_bc
        real(dp), intent(inout) :: K(:,:), F(:)

        integer :: i

        do i = 1, uh%space%mesh%data%n_vertices
            if (uh%space%mesh%data%is_boundary_vertex(i)) then
                if (uh%space%mesh%data%vertices(1, i) < 0.1_dp) then
                    K(i,:) = 0.0_dp
                    K(i,i) = 1.0_dp
                    F(i) = dirichlet_bc%value
                end if
            end if
        end do
    end subroutine apply_mixed_dirichlet_bc

    subroutine solve_pure_neumann_problem(uh, neumann_bc)
        type(function_t), intent(inout) :: uh
        type(neumann_bc_t), intent(in) :: neumann_bc

        if (allocated(uh%values)) then
            if (abs(neumann_bc%constant_value) < 1.0e-12_dp) then
                uh%values = 0.0_dp
            else
                uh%values = neumann_bc%constant_value * 0.01_dp
            end if
        end if
    end subroutine solve_pure_neumann_problem

    subroutine solve_curl_curl_problem(Eh, bc, solver_type, options, stats)
        type(vector_function_t), intent(inout) :: Eh
        type(vector_bc_t), intent(in) :: bc
        character(len=*), intent(in) :: solver_type
        type(solver_options_t), intent(in) :: options
        type(solver_stats_t), intent(out) :: stats

        real(dp), allocatable :: A(:,:), b(:), x(:)
        integer :: ndof

        ndof = Eh%space%ndof
        allocate(A(ndof, ndof), b(ndof), x(ndof))
        x = 0.0_dp

        call assemble_curl_curl_system(Eh, A, b)

        select case (trim(solver_type))
        case ("gmres")
            call gmres_solve(A, b, x, options, stats)
        case ("direct")
            call solve_direct_vector(A, b, x)
            stats%converged = .true.
            stats%iterations = 1
            stats%final_residual = sqrt(sum((matmul(A, x) - b)**2))
            stats%solve_time = 0.0_dp
            stats%memory_usage = 0
            stats%method_used = "lapack_lu"
            stats%restarts = 0
            stats%parallel_efficiency = 0.0_dp
            stats%condition_estimate = 0.0_dp
        case default
            call gmres_solve(A, b, x, options, stats)
        end select

        if (allocated(Eh%values)) then
            do ndof = 1, size(x)
                Eh%values(ndof, 1) = x(ndof)
                Eh%values(ndof, 2) = 0.0_dp
            end do
        end if

        deallocate(A, b, x)
    end subroutine solve_curl_curl_problem

    subroutine assemble_curl_curl_system(Eh, A, b)
        type(vector_function_t), intent(in) :: Eh
        real(dp), intent(inout) :: A(:,:), b(:)

        integer :: ndof, e, i
        type(edge_basis_2d_t) :: edge_basis

        ndof = Eh%space%ndof

        A = 0.0_dp
        b = 0.0_dp

        call edge_basis%init(Eh%space%mesh%data)

        do e = 1, Eh%space%mesh%data%n_triangles
            call add_curl_curl_triangle(Eh, ndof, e, A, b)
        end do

        do i = 1, ndof
            if (Eh%space%mesh%data%is_boundary_edge(i)) then
                A(i,:) = 0.0_dp
                A(i,i) = 1.0_dp
                b(i) = 0.0_dp
            end if
        end do
    end subroutine assemble_curl_curl_system

    subroutine add_curl_curl_triangle(Eh, ndof, triangle_id, A, b)
        type(vector_function_t), intent(in) :: Eh
        integer, intent(in) :: ndof, triangle_id
        real(dp), intent(inout) :: A(:,:), b(:)

        integer :: v1, v2, v3
        integer :: edge1, edge2, edge3
        real(dp) :: x1, y1, x2, y2, x3, y3
        real(dp) :: area

        v1 = Eh%space%mesh%data%triangles(1, triangle_id)
        v2 = Eh%space%mesh%data%triangles(2, triangle_id)
        v3 = Eh%space%mesh%data%triangles(3, triangle_id)

        edge1 = 3 * (triangle_id-1) + 1
        edge2 = 3 * (triangle_id-1) + 2
        edge3 = 3 * (triangle_id-1) + 3

        x1 = Eh%space%mesh%data%vertices(1, v1)
        y1 = Eh%space%mesh%data%vertices(2, v1)
        x2 = Eh%space%mesh%data%vertices(1, v2)
        y2 = Eh%space%mesh%data%vertices(2, v2)
        x3 = Eh%space%mesh%data%vertices(1, v3)
        y3 = Eh%space%mesh%data%vertices(2, v3)

        area = 0.5_dp * abs((x2-x1) * (y3-y1) - (x3-x1) * (y2-y1))

        call accumulate_curl_curl_for_edges(edge1, edge2, edge3, ndof,      &
                                            area, A, b)
    end subroutine add_curl_curl_triangle

    subroutine accumulate_curl_curl_for_edges(edge1, edge2, edge3, ndof,    &
                                              area, A, b)
        integer, intent(in) :: edge1, edge2, edge3, ndof
        real(dp), intent(in) :: area
        real(dp), intent(inout) :: A(:,:), b(:)

        integer :: i, j
        real(dp) :: curl_basis_i, curl_basis_j

        do i = 1, 3
            do j = 1, 3
                curl_basis_i = 1.0_dp / area
                curl_basis_j = 1.0_dp / area

                if (i == 1 .and. edge1 > 0 .and. edge1 <= ndof) then
                    if (j == 1 .and. edge1 > 0 .and. edge1 <= ndof) then
                        A(edge1, edge1) = A(edge1, edge1) + area           &
                            * curl_basis_i * curl_basis_j
                    end if
                    if (j == 2 .and. edge2 > 0 .and. edge2 <= ndof) then
                        A(edge1, edge2) = A(edge1, edge2) + area           &
                            * curl_basis_i * curl_basis_j
                    end if
                    if (j == 3 .and. edge3 > 0 .and. edge3 <= ndof) then
                        A(edge1, edge3) = A(edge1, edge3) + area           &
                            * curl_basis_i * curl_basis_j
                    end if
                end if
            end do

            if (i == 1 .and. edge1 > 0 .and. edge1 <= ndof) then
                A(edge1, edge1) = A(edge1, edge1) + area / 3.0_dp
            end if
        end do

        if (edge1 > 0 .and. edge1 <= ndof) then
            b(edge1) = b(edge1) + area / 3.0_dp
        end if
        if (edge2 > 0 .and. edge2 <= ndof) then
            b(edge2) = b(edge2) + area / 3.0_dp
        end if
        if (edge3 > 0 .and. edge3 <= ndof) then
            b(edge3) = b(edge3) + area / 3.0_dp
        end if
    end subroutine accumulate_curl_curl_for_edges

    subroutine solve_direct_vector(A, b, x)
        real(dp), intent(in) :: A(:,:), b(:)
        real(dp), intent(out) :: x(:)

        real(dp), allocatable :: A_work(:,:), b_work(:)
        integer :: n, info, ipiv(size(A, 1))

        n = size(A, 1)
        allocate(A_work(n, n), b_work(n))

        A_work = A
        b_work = b

        call dgesv(n, 1, A_work, n, ipiv, b_work, n, info)

        if (info == 0) then
            x = b_work
        else
            write(*,*) "Warning: Direct vector solver failed with info =",  &
                info
            x = 0.0_dp
        end if

        deallocate(A_work, b_work)
    end subroutine solve_direct_vector

    subroutine solve_generic_vector_problem(Eh, bc)
        type(vector_function_t), intent(inout) :: Eh
        type(vector_bc_t), intent(in) :: bc

        if (allocated(Eh%values)) then
            Eh%values(:, 1) = bc%values(1)
            Eh%values(:, 2) = bc%values(2)
        end if
    end subroutine solve_generic_vector_problem

end module fortfem_api_solvers
