program simple_poisson
    use fortfem_kinds
    use fortfem_api
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
    
    ! Plot mesh
    call plot(mesh, filename="poisson_mesh.png", title="Poisson Mesh (20x20)")
    
    ! Plot solution
    call plot(uh, filename="poisson_solution.png", &
              title="Poisson Solution: -Δu = 1", &
              colormap="viridis")

    write(*,*) "Simple Poisson example completed!"
    write(*,*) "Generated files:"
    write(*,*) "  - Mesh: poisson_mesh.png"
    write(*,*) "  - Solution: poisson_solution.png"

end program simple_poisson