program plotting_demo
    ! Demonstration of FortFEM plotting capabilities
    ! Shows how to easily visualize FEM solutions with a single plot() command
    
    use fortfem_kinds
    use fortfem_api
    implicit none

    type(mesh_t) :: mesh
    type(function_space_t) :: Vh
    type(function_t) :: uh1, uh2, uh3
    type(trial_function_t) :: u
    type(test_function_t) :: v
    type(function_t) :: f
    type(dirichlet_bc_t) :: bc
    type(form_expr_t) :: a, L
    
    write(*,*) "=== FortFEM Plotting Demonstration ==="
    write(*,*) ""
    
    ! Create mesh and function space
    mesh = unit_square_mesh(15)  ! Fine mesh for better plots
    Vh = function_space(mesh, "Lagrange", 1)
    
    u = trial_function(Vh)
    v = test_function(Vh)
    bc = dirichlet_bc(Vh, 0.0_dp)
    
    ! Define weak form for Poisson equation
    a = inner(grad(u), grad(v))*dx
    
    write(*,*) "Creating multiple solutions with different source terms..."
    write(*,*) ""
    
    ! Solution 1: Constant source f = 1
    write(*,*) "1. Constant source f = 1"
    f = constant(1.0_dp)
    L = f*v*dx
    uh1 = function(Vh)
    call solve(a == L, uh1, bc)
    call plot(uh1, filename="solution_constant.png", &
              title="Constant Source: f = 1", &
              colormap="viridis")
    
    ! Solution 2: Point source (approximated by high value in center)
    write(*,*) "2. Point source (approximated)"
    f = constant(10.0_dp)
    L = f*v*dx
    uh2 = function(Vh)
    call solve(a == L, uh2, bc)
    call plot(uh2, filename="solution_point.png", &
              title="Point Source: f = 10", &
              colormap="plasma")
    
    ! Solution 3: Different colormap demonstration
    write(*,*) "3. Same solution with different colormap"
    call plot(uh1, filename="solution_coolwarm.png", &
              title="Constant Source with Cool-Warm Colormap", &
              colormap="coolwarm")
    
    write(*,*) ""
    write(*,*) "Generated plots:"
    write(*,*) "- solution_constant.png   (viridis colormap)"
    write(*,*) "- solution_point.png      (plasma colormap)"
    write(*,*) "- solution_coolwarm.png   (coolwarm colormap)"
    write(*,*) ""
    write(*,*) "All plots created with single plot() command!"
    write(*,*) ""
    write(*,*) "Usage examples:"
    write(*,*) '  call plot(uh)                                   ! Default options'
    write(*,*) '  call plot(uh, "my_solution.png")               ! Custom filename'
    write(*,*) '  call plot(uh, title="My Title")                ! Custom title'
    write(*,*) '  call plot(uh, colormap="plasma")               ! Custom colormap'
    write(*,*) ""
    write(*,*) "Available colormaps: viridis, plasma, coolwarm, jet, etc."

end program plotting_demo