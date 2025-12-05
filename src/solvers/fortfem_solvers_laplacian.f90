module fortfem_solvers_laplacian
    use fortfem_kinds, only: dp
    use fortfem_api_types, only: function_space_t, function_t, dirichlet_bc_t
    use fortfem_api_forms, only: form_equation_t
    use fortfem_advanced_solvers, only: solver_options_t, solver_stats_t, &
                                        solver_options, advanced_solve => solve
    use basis_q1_quad_2d_module, only: q1_shape_functions, &
                                       q1_shape_derivatives, q1_jacobian
    use fortfem_solvers_p2, only: solve_laplacian_problem_p2
    implicit none

    private

    public :: assemble_laplacian_system
    public :: solve_scalar
    public :: solve_laplacian_problem
    public :: solve_generic_problem
    public :: add_p1_triangle_contribution

contains

    subroutine solve_scalar(equation, uh, bc, options, stats)
        type(form_equation_t), intent(in) :: equation
        type(function_t), intent(inout) :: uh
        type(dirichlet_bc_t), intent(in) :: bc
        type(solver_options_t), intent(in), optional :: options
        type(solver_stats_t), intent(out), optional :: stats

        type(solver_stats_t) :: local_stats
        logical :: have_stats

        write (*, *) "Solving: ", trim(equation%lhs%description), " == ", &
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

    subroutine assemble_laplacian_system(space, bc, K, F)
        type(function_space_t), intent(in) :: space
        type(dirichlet_bc_t), intent(in) :: bc
        real(dp), allocatable, intent(out) :: K(:, :), F(:)

        integer :: ndof, i

        ndof = space%ndof
        allocate (K(ndof, ndof), F(ndof))

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

    pure subroutine compute_p1_triangle_gradients(x1, y1, x2, y2, x3, y3, &
                                                  area, grad_x, grad_y)
        real(dp), intent(in) :: x1, y1, x2, y2, x3, y3
        real(dp), intent(out) :: area
        real(dp), intent(out) :: grad_x(3), grad_y(3)

        real(dp) :: a11, a12, a21, a22, det_a

        a11 = x2 - x1
        a12 = x3 - x1
        a21 = y2 - y1
        a22 = y3 - y1

        det_a = a11*a22 - a12*a21

        area = 0.5_dp*abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))

        grad_x(1) = (-a22 + a21)/det_a
        grad_y(1) = (a12 - a11)/det_a

        grad_x(2) = a22/det_a
        grad_y(2) = -a12/det_a

        grad_x(3) = -a21/det_a
        grad_y(3) = a11/det_a
    end subroutine compute_p1_triangle_gradients

    subroutine add_p1_triangle_contribution(space, triangle_id, K, F)
        type(function_space_t), intent(in) :: space
        integer, intent(in) :: triangle_id
        real(dp), intent(inout) :: K(:, :), F(:)

        integer :: v1, v2, v3
        integer :: i, j
        real(dp) :: x1, y1, x2, y2, x3, y3
        real(dp) :: area
        real(dp) :: bx(3), by(3), K_elem(3, 3)

        v1 = space%mesh%data%triangles(1, triangle_id)
        v2 = space%mesh%data%triangles(2, triangle_id)
        v3 = space%mesh%data%triangles(3, triangle_id)

        x1 = space%mesh%data%vertices(1, v1)
        y1 = space%mesh%data%vertices(2, v1)
        x2 = space%mesh%data%vertices(1, v2)
        y2 = space%mesh%data%vertices(2, v2)
        x3 = space%mesh%data%vertices(1, v3)
        y3 = space%mesh%data%vertices(2, v3)

        call compute_p1_triangle_gradients(x1, y1, x2, y2, x3, y3, area, bx, &
                                           by)

        do i = 1, 3
            do j = 1, 3
                K_elem(i, j) = area*(bx(i)*bx(j) + by(i)*by(j))
            end do
        end do

        K(v1, v1) = K(v1, v1) + K_elem(1, 1)
        K(v1, v2) = K(v1, v2) + K_elem(1, 2)
        K(v1, v3) = K(v1, v3) + K_elem(1, 3)
        K(v2, v1) = K(v2, v1) + K_elem(2, 1)
        K(v2, v2) = K(v2, v2) + K_elem(2, 2)
        K(v2, v3) = K(v2, v3) + K_elem(2, 3)
        K(v3, v1) = K(v3, v1) + K_elem(3, 1)
        K(v3, v2) = K(v3, v2) + K_elem(3, 2)
        K(v3, v3) = K(v3, v3) + K_elem(3, 3)

        F(v1) = F(v1) + area/3.0_dp
        F(v2) = F(v2) + area/3.0_dp
        F(v3) = F(v3) + area/3.0_dp
    end subroutine add_p1_triangle_contribution

    subroutine assemble_laplacian_triangles(space, K, F)
        type(function_space_t), intent(in) :: space
        real(dp), intent(inout) :: K(:, :), F(:)

        integer :: e

        do e = 1, space%mesh%data%n_triangles
            call add_p1_triangle_contribution(space, e, K, F)
        end do
    end subroutine assemble_laplacian_triangles

    subroutine add_q1_quad_contribution(space, quad_id, K, F)
        type(function_space_t), intent(in) :: space
        integer, intent(in) :: quad_id
        real(dp), intent(inout) :: K(:, :), F(:)

        integer :: i, j, vi
        integer :: v_ids(4)
        real(dp) :: coords(2, 4)
        real(dp) :: K_elem(4, 4), F_elem(4)

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
        real(dp), intent(in) :: coords(2, 4)
        real(dp), intent(out) :: K_elem(4, 4), F_elem(4)

        integer :: i, j, kx, ky
        real(dp) :: jac(2, 2), det_jac, inv_jac(2, 2)
        real(dp) :: dN_dxi(4), dN_deta(4)
        real(dp) :: grad_ref(2), grad_phys(2, 4)
        real(dp) :: N(4)
        real(dp) :: xi, eta, weight
        logical :: success
        real(dp), parameter :: gauss_pts(2) = [-0.5773502691896257_dp, &
                                               0.5773502691896257_dp]
        real(dp), parameter :: gauss_w(2) = [1.0_dp, 1.0_dp]

        K_elem = 0.0_dp
        F_elem = 0.0_dp
        do kx = 1, 2
            do ky = 1, 2
                xi = gauss_pts(kx)
                eta = gauss_pts(ky)

                call q1_shape_derivatives(xi, eta, dN_dxi, dN_deta)
                call q1_jacobian(xi, eta, coords, jac, det_jac, inv_jac, &
                                 success)

                if (.not. success) cycle

                do i = 1, 4
                    grad_ref(1) = dN_dxi(i)
                    grad_ref(2) = dN_deta(i)
                    grad_phys(:, i) = matmul(transpose(inv_jac), grad_ref)
                end do

                call q1_shape_functions(xi, eta, N)

                weight = det_jac*gauss_w(kx)*gauss_w(ky)

                do i = 1, 4
                    do j = 1, 4
                        K_elem(i, j) = K_elem(i, j) + weight* &
                                       (grad_phys(1, i)*grad_phys(1, j) + &
                                        grad_phys(2, i)*grad_phys(2, j))
                    end do
                    F_elem(i) = F_elem(i) + weight*N(i)
                end do
            end do
        end do
    end subroutine compute_q1_quad_element

    subroutine assemble_laplacian_quads(space, K, F)
        type(function_space_t), intent(in) :: space
        real(dp), intent(inout) :: K(:, :), F(:)

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

        real(dp), allocatable :: K(:, :), F(:)
        integer :: ndof
        type(solver_options_t) :: local_opts

        ndof = uh%space%ndof

        call assemble_laplacian_system(uh%space, bc, K, F)

        if (.not. allocated(uh%values)) then
            allocate (uh%values(ndof))
        end if
        uh%values = 0.0_dp

        if (present(options)) then
            local_opts = options
        else
            local_opts = solver_options(method="auto")
        end if

        call advanced_solve(K, F, uh%values, local_opts, stats)

        if (.not. stats%converged) then
            write (*, *) "Warning: Laplacian solver did not report convergence.", &
                " Final residual =", stats%final_residual
        end if

        deallocate (K, F)
    end subroutine solve_laplacian_problem

    subroutine solve_generic_problem(uh, bc)
        type(function_t), intent(inout) :: uh
        type(dirichlet_bc_t), intent(in) :: bc

        if (allocated(uh%values)) then
            uh%values = bc%value
        end if
    end subroutine solve_generic_problem

end module fortfem_solvers_laplacian

