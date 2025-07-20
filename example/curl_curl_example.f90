program curl_curl_example
    ! Clean FEniCS-style example solving the curl-curl equation
    ! ∇ × (∇ × E) + E = J in Ω = [0,1]²
    ! E × n = 0 on ∂Ω (tangential boundary condition)
    ! 
    ! This demonstrates FEniCS-style syntax for electromagnetic problems
    ! using edge elements (Nédélec) and iterative GMRES solver

    use fortfem_kinds
    use fortfem_api
    implicit none

    type(mesh_t) :: mesh
    type(vector_function_space_t) :: Vh
    type(vector_trial_function_t) :: E
    type(vector_test_function_t) :: F
    type(vector_function_t) :: J, Eh
    type(vector_bc_t) :: bc
    type(form_expr_t) :: a, L
    integer :: n_vertices, n_elements, n_dofs
    real(dp) :: h, max_E

    write(*,*) "=== Curl-Curl Electromagnetic Example ==="
    write(*,*) ""
    write(*,*) "Solving: ∇ × (∇ × E) + E = J on [0,1]²"
    write(*,*) "with E × n = 0 on the boundary"
    write(*,*) ""

    ! Create mesh and vector function space using FEniCS-style API
    mesh = unit_square_mesh(8)  ! 8x8 grid for efficiency
    Vh = vector_function_space(mesh, "Nedelec", 1)
    
    n_vertices = mesh%data%n_vertices
    n_elements = mesh%data%n_triangles
    n_dofs = Vh%ndof
    h = 1.0_dp / 7.0_dp  ! mesh spacing
    
    write(*,*) "Mesh statistics:"
    write(*,*) "  Vertices:", n_vertices
    write(*,*) "  Elements:", n_elements
    write(*,*) "  Vector DOFs:", n_dofs
    write(*,*) "  h =", h
    write(*,*) ""

    ! Define trial and test functions using natural notation
    E = vector_trial_function(Vh)
    F = vector_test_function(Vh)

    ! Define source current J = [1, 0] (x-directed current)
    J = vector_function(Vh)
    J%values(:, 1) = 1.0_dp  ! Jx = 1
    J%values(:, 2) = 0.0_dp  ! Jy = 0

    ! Define weak form using mathematical notation
    a = inner(curl(E), curl(F))*dx + inner(E, F)*dx  ! Bilinear form
    L = inner(J, F)*dx                                ! Linear form

    write(*,*) "Weak form:"
    write(*,*) "  a(E,F) = ", trim(a%description)
    write(*,*) "  L(F)   = ", trim(L%description)
    write(*,*) ""

    ! Set up tangential boundary conditions: E × n = 0
    bc = vector_bc(Vh, [0.0_dp, 0.0_dp], "tangential")

    ! Create solution function
    Eh = vector_function(Vh)

    ! Solve the system using GMRES iterative solver
    write(*,*) "Assembling curl-curl matrix..."
    write(*,*) "Applying tangential boundary conditions..."
    write(*,*) "Solving with GMRES iterative solver..."
    call solve(a == L, Eh, bc)
    write(*,*) ""

    ! Analyze solution
    max_E = maxval(sqrt(Eh%values(:,1)**2 + Eh%values(:,2)**2))
    
    write(*,*) "Solution statistics:"
    write(*,*) "  Max |E| =", max_E
    write(*,*) "  Max Ex =", maxval(abs(Eh%values(:,1)))
    write(*,*) "  Max Ey =", maxval(abs(Eh%values(:,2)))
    write(*,*) ""
    
    ! Plot the vector solution
    write(*,*) "Creating vector field plot..."
    call plot(Eh, filename="curl_curl_solution.png", &
              label="Curl-Curl Solution: ∇×(∇×E) + E = J", &
              plot_type="streamplot")
    write(*,*) ""
    
    write(*,*) "Example completed successfully!"
    write(*,*) ""
    write(*,*) "This example demonstrates:"
    write(*,*) "- FEniCS-style syntax for vector problems"
    write(*,*) "- Curl operator with edge elements (Nédélec)"
    write(*,*) "- Tangential boundary conditions"
    write(*,*) "- GMRES iterative solver for large systems"
    write(*,*) "- Vector field visualization with streamplots"

end program curl_curl_example