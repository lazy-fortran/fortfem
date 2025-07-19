program poisson_forms
    ! Example of target API with FEniCS-style forms in Fortran
    ! This shows the desired syntax - not yet implemented
    use fortfem
    implicit none
    
    ! Types
    type(mesh_t) :: mesh
    type(function_space_t) :: Vh
    type(trial_function_t) :: u
    type(test_function_t) :: v  
    type(function_t) :: uh, f, g
    type(dirichlet_bc_t), allocatable :: bcs(:)
    type(form_t) :: a, L
    
    ! Create mesh and function space
    mesh = unit_square_mesh(32, 32)
    Vh = function_space(mesh, "Lagrange", 1)
    
    ! Define functions
    u = trial_function(Vh)
    v = test_function(Vh)
    
    ! Source term: f = 2π²sin(πx)sin(πy)
    f = function(Vh)
    call f%interpolate(source_function)
    
    ! Bilinear and linear forms
    a = inner(grad(u), grad(v))*dx
    L = f*v*dx
    
    ! Boundary conditions (all edges)
    allocate(bcs(1))
    bcs(1) = dirichlet_bc(Vh, 0.0_dp, "on_boundary")
    
    ! Solve
    uh = function(Vh)
    call solve(a == L, uh, bcs)
    
    ! Error norm
    g = function(Vh)
    call g%interpolate(exact_solution)
    print *, "L2 error:", norm(uh - g, "L2")
    
    ! Save solution
    call uh%write_vtk("poisson_solution.vtk")
    
contains

    real(dp) function source_function(x, y)
        real(dp), intent(in) :: x, y
        source_function = 2.0_dp * pi**2 * sin(pi*x) * sin(pi*y)
    end function
    
    real(dp) function exact_solution(x, y)
        real(dp), intent(in) :: x, y
        exact_solution = sin(pi*x) * sin(pi*y)
    end function

end program poisson_forms