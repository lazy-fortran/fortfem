program fenics_style_poisson
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

    write(*,*) "=== FEniCS-Style Poisson Equation ==="
    write(*,*) ""
    write(*,*) "Solving: -Δu = 1 on [0,1]²"
    write(*,*) "with u = 0 on the boundary"
    write(*,*) ""

    ! Create mesh and function space
    mesh = unit_square_mesh(10)
    Vh = function_space(mesh, "Lagrange", 1)
    
    write(*,*) "Mesh statistics:"
    write(*,*) "  Vertices:", mesh%data%n_vertices
    write(*,*) "  Elements:", mesh%data%n_triangles
    write(*,*) ""

    ! Define trial and test functions  
    u = trial_function(Vh)
    v = test_function(Vh)

    ! Define source function f = 1
    f = constant(1.0_dp)

    ! Define weak form using FEniCS syntax
    a = inner(grad(u), grad(v))*dx
    L = f*v*dx

    write(*,*) "Weak form:"
    write(*,*) "  a(u,v) = ", trim(a%description)
    write(*,*) "  L(v)   = ", trim(L%description)
    write(*,*) ""

    ! Set up boundary conditions
    bc = dirichlet_bc(Vh, 0.0_dp)

    ! Create solution function
    uh = function(Vh)

    ! Solve the system using FEniCS-style interface
    write(*,*) "Solving linear system..."
    call solve(a == L, uh, bc)
    write(*,*) ""

    ! Display results
    write(*,*) "Solution statistics:"
    write(*,*) "  Max u =", maxval(uh%values)
    write(*,*) "  Min u =", minval(uh%values)
    write(*,*) ""
    write(*,*) "Note: For -Δu = 1 on [0,1]², the exact maximum is 1/8 = 0.125"
    write(*,*) ""
    write(*,*) "FEniCS-style example completed successfully!"

end program fenics_style_poisson