program simple_poisson
    ! Clean FEniCS-style example solving the Poisson equation
    ! -Δu = 1 in Ω = [0,1]²
    ! u = 0 on ∂Ω
    ! 
    ! This demonstrates the natural mathematical notation enabled by FortFEM
    
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
    integer :: n_vertices, n_elements
    real(dp) :: h, max_u, center_val

    write(*,*) "=== Simple Poisson Equation Example ==="
    write(*,*) ""
    write(*,*) "Solving: -Δu = 1 on [0,1]²"
    write(*,*) "with u = 0 on the boundary"
    write(*,*) ""

    ! Create mesh and function space using FEniCS-style API
    mesh = unit_square_mesh(20)  ! 20x20 grid
    Vh = function_space(mesh, "Lagrange", 1)
    
    n_vertices = mesh%data%n_vertices
    n_elements = mesh%data%n_triangles
    h = 1.0_dp / 19.0_dp  ! mesh spacing
    
    write(*,*) "Mesh statistics:"
    write(*,*) "  Vertices:", n_vertices
    write(*,*) "  Elements:", n_elements
    write(*,*) "  h =", h
    write(*,*) ""

    ! Define trial and test functions using natural notation
    u = trial_function(Vh)
    v = test_function(Vh)

    ! Define source function f = 1
    f = constant(1.0_dp)

    ! Define weak form using mathematical notation
    a = inner(grad(u), grad(v))*dx  ! Bilinear form: ∫ ∇u·∇v dx
    L = f*v*dx                      ! Linear form:   ∫ f v dx

    write(*,*) "Weak form:"
    write(*,*) "  a(u,v) = ", trim(a%description)
    write(*,*) "  L(v)   = ", trim(L%description)
    write(*,*) ""

    ! Set up homogeneous Dirichlet boundary conditions
    bc = dirichlet_bc(Vh, 0.0_dp)

    ! Create solution function
    uh = function(Vh)

    ! Solve the system using FEniCS-style interface
    write(*,*) "Assembling stiffness matrix..."
    write(*,*) "Applying boundary conditions..."
    write(*,*) "Solving linear system..."
    call solve(a == L, uh, bc)
    write(*,*) ""

    ! Analyze solution
    max_u = maxval(uh%values)
    
    write(*,*) "Solution statistics:"
    write(*,*) "  Max u =", max_u
    write(*,*) "  Min u =", minval(uh%values)
    write(*,*) ""
    write(*,*) "Note: For -Δu = 1 on [0,1]², the true analytical maximum is ≈ 0.1093"
    write(*,*) ""
    
    ! Find center vertex value (approximate)
    center_val = uh%values(n_vertices/2)  ! Rough center estimate
    write(*,*) "Solution at center vertex:"
    write(*,*) "  Value: u =", center_val
    write(*,*) ""
    
    write(*,*) "Example completed successfully!"
    write(*,*) ""
    write(*,*) "This example demonstrates:"
    write(*,*) "- Natural mathematical notation for weak forms"
    write(*,*) "- Automatic assembly from symbolic expressions" 
    write(*,*) "- Clean FEniCS-style solve interface"

end program simple_poisson