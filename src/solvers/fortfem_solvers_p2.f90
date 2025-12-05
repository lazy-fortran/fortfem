module fortfem_solvers_p2
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_api_types, only: function_t, dirichlet_bc_t
    use fortfem_api_mesh, only: find_triangle_edges
    use fortfem_advanced_solvers, only: solver_options_t, solver_stats_t, &
        solver_options, advanced_solve => solve
    use basis_p2_2d_module, only: basis_p2_2d_t
    implicit none

    private

    public :: solve_laplacian_problem_p2
    public :: assemble_p2_laplacian_system

contains

    subroutine solve_laplacian_problem_p2(uh, bc, options, stats)
        type(function_t), intent(inout) :: uh
        type(dirichlet_bc_t), intent(in) :: bc
        type(solver_options_t), intent(in), optional :: options
        type(solver_stats_t), intent(out) :: stats

        real(dp), allocatable :: K(:, :), F(:)
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
        real(dp), intent(inout) :: K(:, :), F(:)

        integer :: ndof, i, j, e, v1, v2, v3
        real(dp) :: K_elem(6, 6), F_elem(6)
        integer :: dofs(6), edge1, edge2, edge3
        real(dp) :: vertices(2, 3)

        ndof = uh%space%ndof

        do e = 1, uh%space%mesh%data%n_triangles
            v1 = uh%space%mesh%data%triangles(1, e)
            v2 = uh%space%mesh%data%triangles(2, e)
            v3 = uh%space%mesh%data%triangles(3, e)

            vertices(1, 1) = uh%space%mesh%data%vertices(1, v1)
            vertices(2, 1) = uh%space%mesh%data%vertices(2, v1)
            vertices(1, 2) = uh%space%mesh%data%vertices(1, v2)
            vertices(2, 2) = uh%space%mesh%data%vertices(2, v2)
            vertices(1, 3) = uh%space%mesh%data%vertices(1, v3)
            vertices(2, 3) = uh%space%mesh%data%vertices(2, v3)

            dofs(1:3) = [v1, v2, v3]
            call find_triangle_edges(uh%space%mesh%data, e, edge1, edge2, &
                edge3)
            dofs(4) = uh%space%mesh%data%n_vertices + edge1
            dofs(5) = uh%space%mesh%data%n_vertices + edge2
            dofs(6) = uh%space%mesh%data%n_vertices + edge3

            call compute_p2_element_matrices(vertices, K_elem, F_elem)

            do i = 1, 6
                if (dofs(i) > 0 .and. dofs(i) <= ndof) then
                    do j = 1, 6
                        if (dofs(j) > 0 .and. dofs(j) <= ndof) then
                            K(dofs(i), dofs(j)) = K(dofs(i), dofs(j)) + &
                                K_elem(i, j)
                        end if
                    end do
                    F(dofs(i)) = F(dofs(i)) + F_elem(i)
                end if
            end do
        end do

        call apply_p2_dirichlet_bc(uh%space%mesh%data, bc%value, ndof, K, F)
    end subroutine assemble_p2_laplacian_system

    subroutine compute_p2_element_matrices(vertices, K_elem, F_elem)
        real(dp), intent(in) :: vertices(2, 3)
        real(dp), intent(out) :: K_elem(6, 6), F_elem(6)

        type(basis_p2_2d_t) :: basis_p2
        real(dp) :: xi_q(3), eta_q(3), w_q(3)
        real(dp) :: jac(2, 2), det_j, inv_jac(2, 2), area
        real(dp) :: grad_i(2), grad_j(2)
        integer :: i, j, kq

        xi_q = [1.0_dp/6.0_dp, 2.0_dp/3.0_dp, 1.0_dp/6.0_dp]
        eta_q = [1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 2.0_dp/3.0_dp]
        w_q = [1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp]

        area = 0.5_dp*abs((vertices(1, 2) - vertices(1, 1))* &
            (vertices(2, 3) - vertices(2, 1)) - &
            (vertices(1, 3) - vertices(1, 1))* &
            (vertices(2, 2) - vertices(2, 1)))

        call basis_p2%compute_jacobian(vertices, jac, det_j)
        inv_jac(1, 1) = jac(2, 2)/det_j
        inv_jac(1, 2) = -jac(1, 2)/det_j
        inv_jac(2, 1) = -jac(2, 1)/det_j
        inv_jac(2, 2) = jac(1, 1)/det_j

        K_elem = 0.0_dp
        F_elem = 0.0_dp

        do i = 1, 6
            do j = 1, 6
                do kq = 1, 3
                    grad_i = matmul(inv_jac, basis_p2%grad(i, xi_q(kq), &
                        eta_q(kq)))
                    grad_j = matmul(inv_jac, basis_p2%grad(j, xi_q(kq), &
                        eta_q(kq)))
                    K_elem(i, j) = K_elem(i, j) + w_q(kq)*area* &
                        (grad_i(1)*grad_j(1) + grad_i(2)*grad_j(2))
                end do
            end do
            do kq = 1, 3
                F_elem(i) = F_elem(i) + w_q(kq)*area* &
                    basis_p2%eval(i, xi_q(kq), eta_q(kq))
            end do
        end do
    end subroutine compute_p2_element_matrices

    subroutine apply_p2_dirichlet_bc(mdata, bc_value, ndof, K, F)
        type(mesh_2d_t), intent(in) :: mdata
        real(dp), intent(in) :: bc_value
        integer, intent(in) :: ndof
        real(dp), intent(inout) :: K(:, :), F(:)

        integer :: i, j

        do i = 1, mdata%n_vertices
            if (mdata%is_boundary_vertex(i)) then
                K(i, :) = 0.0_dp
                K(i, i) = 1.0_dp
                F(i) = bc_value
            end if
        end do

        do i = 1, mdata%n_edges
            if (mdata%is_boundary_edge(i)) then
                j = mdata%n_vertices + i
                if (j <= ndof) then
                    K(j, :) = 0.0_dp
                    K(j, j) = 1.0_dp
                    F(j) = bc_value
                end if
            end if
        end do
    end subroutine apply_p2_dirichlet_bc

end module fortfem_solvers_p2

