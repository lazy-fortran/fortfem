program simple_poisson
    use fortfem_core, only: dirichlet_bc_t, function_space_t, mesh_t, &
        unit_square_mesh
    use fortfem_feec, only: constant, dirichlet_bc, dx, form_expr_t, function, &
        function_space, function_t, grad, inner, operator(*), &
        operator(==), solve, test_function, test_function_t, trial_function, &
        trial_function_t
    use fortfem_kinds, only: dp
    use fortfem_plot, only: plot
    implicit none

    type(mesh_t) :: mesh
    type(function_space_t) :: Vh
    type(trial_function_t) :: u
    type(test_function_t) :: v
    type(function_t) :: f, uh
    type(dirichlet_bc_t) :: bc
    type(form_expr_t) :: a, L

    mesh = unit_square_mesh(20)
    Vh = function_space(mesh, "Lagrange", 1)

    u = trial_function(Vh)
    v = test_function(Vh)
    f = constant(1.0_dp)

    a = inner(grad(u), grad(v))*dx
    L = f*v*dx

    bc = dirichlet_bc(Vh, 0.0_dp)
    uh = function(Vh)

    call solve(a == L, uh, bc)

    ! Put the solved field first so the gallery opens on the physical result.
    call plot(uh, filename="poisson_solution.png", &
        title="Poisson Solution: -Δu = 1", &
        colormap="viridis")
    call plot(uh, filename="primary.png", &
        title="Poisson solution: -Δu = 1", colormap="viridis")

    ! Mesh and diagnostic views follow the solution.
    call plot(mesh, filename="poisson_mesh.png", title="Poisson Mesh (20x20)")
    write(*,*) "Simple Poisson example completed!"
    write(*,*) "Generated files:"
    write(*,*) "  - Mesh: poisson_mesh.png"
    write(*,*) "  - Solution: poisson_solution.png"

end program simple_poisson
