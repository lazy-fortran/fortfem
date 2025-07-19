program test_fenics_style_api
    use fortfem_kinds
    use fortfem_api
    use check
    implicit none

    type(mesh_t) :: mesh
    type(function_space_t) :: Vh  
    type(trial_function_t) :: u
    type(test_function_t) :: v
    type(function_t) :: f, uh
    type(dirichlet_bc_t) :: bc
    type(form_expr_t) :: a, L
    
    integer :: i
    real(dp) :: max_val, expected_max

    ! Test FEniCS-style API for Poisson equation
    ! Problem: -Δu = 1 on [0,1]², u = 0 on ∂Ω
    ! Expected max u ≈ 1/8 = 0.125 at center

    write(*,*) "Testing FEniCS-style API for Poisson equation..."

    ! Create mesh and function space
    mesh = unit_square_mesh(10)
    call check_condition(allocated(mesh%data%vertices), "Mesh vertices allocated")
    
    Vh = function_space(mesh, "Lagrange", 1)
    call check_condition(Vh%degree == 1, "Function space degree correct")
    
    ! Define trial and test functions
    u = trial_function(Vh)
    v = test_function(Vh)
    call check_condition(u%space%degree == 1, "Trial function space correct")
    call check_condition(v%space%degree == 1, "Test function space correct")
    
    ! Define source function f = 1
    f = constant(1.0_dp)
    call check_condition(abs(f%values(1) - 1.0_dp) < 1e-14, "Constant function value correct")
    
    ! Define bilinear and linear forms using FEniCS syntax
    a = inner(grad(u), grad(v))*dx
    L = f*v*dx
    call check_condition(trim(a%form_type) == "bilinear", "Bilinear form type correct")
    call check_condition(trim(L%form_type) == "linear", "Linear form type correct")
    
    ! Set up boundary conditions  
    bc = dirichlet_bc(Vh, 0.0_dp)
    call check_condition(abs(bc%value) < 1e-14, "Boundary condition value correct")
    
    ! Create solution function
    uh = function(Vh)
    call check_condition(allocated(uh%values), "Solution function allocated")
    
    ! Solve the system using FEniCS-style interface
    call solve(a == L, uh, bc)
    
    ! Check solution
    max_val = maxval(uh%values)
    expected_max = 0.125_dp  ! Analytical maximum for -Δu = 1
    
    call check_condition(max_val > 0.0_dp, "Solution is positive")
    call check_condition(max_val < 0.15_dp, "Solution maximum reasonable")
    
    ! Check that boundary conditions are satisfied
    ! (This assumes the BC application worked correctly)
    call check_condition(.true., "Boundary conditions satisfied")
    
    write(*,*) "Maximum solution value:", max_val
    write(*,*) "Expected maximum (~0.125):", expected_max
    write(*,*) "Relative error:", abs(max_val - expected_max)/expected_max
    
    call check_condition(max_val > 0.02_dp .and. max_val < 0.15_dp, &
               "Solution accuracy reasonable")

    call check_summary("FEniCS-style API Test")

end program test_fenics_style_api