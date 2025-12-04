program test_poisson_quads
    use fortfem_kinds, only: dp
    use fortfem_api, only: mesh_t, function_space_t, trial_function_t,     &
        test_function_t, function_t, dirichlet_bc_t, form_expr_t,          &
        unit_square_mesh, structured_quad_mesh, function_space,            &
        trial_function, test_function, constant, dirichlet_bc,             &
        function, inner, grad, dx, solve, plot, operator(*), operator(==)
    use check, only: check_condition, check_summary
    implicit none

    write(*,*) "Testing Poisson solver on Q1 quadrilateral meshes..."

    call test_poisson_quads_basic()
    call test_poisson_quads_vs_triangles()
    call test_poisson_quads_plot()

    call check_summary("Poisson on Quadrilateral Meshes")

contains

    subroutine test_poisson_quads_basic()
        type(mesh_t) :: quad_mesh
        type(function_space_t) :: Vh_q
        type(trial_function_t) :: u_q
        type(test_function_t) :: v_q
        type(function_t) :: f_q, uh_q
        type(dirichlet_bc_t) :: bc_q
        type(form_expr_t) :: a_q, L_q
        real(dp) :: max_val, max_interior, center_val
        real(dp) :: x, y, min_dist, dist
        integer :: i, center_node
        logical :: is_boundary

        quad_mesh = structured_quad_mesh(8, 8, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        Vh_q = function_space(quad_mesh, "Lagrange", 1)

        u_q = trial_function(Vh_q)
        v_q = test_function(Vh_q)
        f_q = constant(1.0_dp)

        a_q = inner(grad(u_q), grad(v_q)) * dx
        L_q = f_q * v_q * dx

        bc_q = dirichlet_bc(Vh_q, 0.0_dp)
        uh_q = function(Vh_q)

        call solve(a_q == L_q, uh_q, bc_q)

        ! Basic sanity checks on solution
        max_val = maxval(uh_q%values)
        max_interior = 0.0_dp
        min_dist = huge(1.0_dp)
        center_node = 1

        do i = 1, Vh_q%ndof
            is_boundary = quad_mesh%data%is_boundary_vertex(i)
            x = quad_mesh%data%vertices(1, i)
            y = quad_mesh%data%vertices(2, i)

            if (is_boundary) then
                call check_condition(abs(uh_q%values(i)) < 1.0e-12_dp, &
                    "Quad Poisson: boundary value approximately zero")
            else
                max_interior = max(max_interior, uh_q%values(i))
                dist = (x - 0.5_dp)**2 + (y - 0.5_dp)**2
                if (dist < min_dist) then
                    min_dist = dist
                    center_node = i
                end if
            end if
        end do

        center_val = uh_q%values(center_node)

        call check_condition(max_interior > 0.0_dp, &
            "Quad Poisson: interior solution positive")
        call check_condition(max_interior < 0.5_dp, &
            "Quad Poisson: interior solution bounded")
        call check_condition(abs(max_val - max_interior) < 1.0e-10_dp, &
            "Quad Poisson: maximum attained in interior")
        call check_condition(center_val > 0.8_dp * max_interior, &
            "Quad Poisson: maximum located near center")

        write(*,*) "   Quad Poisson basic: max =", max_interior
        write(*,*) "   Quad Poisson basic: center node =", center_node
        write(*,*) "   Quad Poisson basic: center value =", center_val
    end subroutine test_poisson_quads_basic

    subroutine test_poisson_quads_vs_triangles()
        type(mesh_t) :: quad_mesh, tri_mesh
        type(function_space_t) :: Vh_q, Vh_t
        type(trial_function_t) :: u_q, u_t
        type(test_function_t) :: v_q, v_t
        type(function_t) :: f_q, f_t, uh_q, uh_t
        type(dirichlet_bc_t) :: bc_q, bc_t
        type(form_expr_t) :: a_q, L_q, a_t, L_t
        real(dp) :: max_q, max_t, rel_diff

        quad_mesh = structured_quad_mesh(8, 8, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        tri_mesh = unit_square_mesh(8)

        Vh_q = function_space(quad_mesh, "Lagrange", 1)
        Vh_t = function_space(tri_mesh, "Lagrange", 1)

        u_q = trial_function(Vh_q)
        v_q = test_function(Vh_q)
        f_q = constant(1.0_dp)
        a_q = inner(grad(u_q), grad(v_q)) * dx
        L_q = f_q * v_q * dx
        bc_q = dirichlet_bc(Vh_q, 0.0_dp)
        uh_q = function(Vh_q)
        call solve(a_q == L_q, uh_q, bc_q)

        u_t = trial_function(Vh_t)
        v_t = test_function(Vh_t)
        f_t = constant(1.0_dp)
        a_t = inner(grad(u_t), grad(v_t)) * dx
        L_t = f_t * v_t * dx
        bc_t = dirichlet_bc(Vh_t, 0.0_dp)
        uh_t = function(Vh_t)
        call solve(a_t == L_t, uh_t, bc_t)

        max_q = maxval(uh_q%values)
        max_t = maxval(uh_t%values)

        if (max_t > 0.0_dp) then
            rel_diff = abs(max_q - max_t) / max_t
        else
            rel_diff = 0.0_dp
        end if

        call check_condition(max_q > 0.0_dp .and. max_t > 0.0_dp, &
            "Quad vs Tri Poisson: positive maxima")
        call check_condition(rel_diff < 0.2_dp, &
            "Quad vs Tri Poisson: maxima reasonably close")

        write(*,*) "   Quad vs Tri: max_quad =", max_q, " max_tri =", max_t
        write(*,*) "   Quad vs Tri: relative difference =", rel_diff
    end subroutine test_poisson_quads_vs_triangles

    subroutine test_poisson_quads_plot()
        type(mesh_t) :: quad_mesh
        type(function_space_t) :: Vh_q
        type(trial_function_t) :: u_q
        type(test_function_t) :: v_q
        type(function_t) :: f_q, uh_q
        type(dirichlet_bc_t) :: bc_q
        type(form_expr_t) :: a_q, L_q

        quad_mesh = structured_quad_mesh(12, 12, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        Vh_q = function_space(quad_mesh, "Lagrange", 1)

        u_q = trial_function(Vh_q)
        v_q = test_function(Vh_q)
        f_q = constant(1.0_dp)

        a_q = inner(grad(u_q), grad(v_q)) * dx
        L_q = f_q * v_q * dx

        bc_q = dirichlet_bc(Vh_q, 0.0_dp)
        uh_q = function(Vh_q)

        call solve(a_q == L_q, uh_q, bc_q)

        call plot(uh_q, filename="build/poisson_quads_solution.png", &
                  title="Poisson solution on Q1 quadrilateral mesh")

        call check_condition(any(uh_q%values /= 0.0_dp), &
            "Quad Poisson plot: non-trivial solution plotted")

        write(*,*) "   Quad Poisson plot: generated build/poisson_quads_solution.png"
    end subroutine test_poisson_quads_plot

end program test_poisson_quads
