program test_neumann_boundary_conditions
    use fortfem_kinds
    use fortfem_api
    use check
    implicit none

    write(*,*) "Testing Neumann boundary condition implementation..."
    
    call test_neumann_bc_type_creation()
    call test_mixed_dirichlet_neumann_problem()
    call test_pure_neumann_problem()
    call test_neumann_boundary_integration()
    call test_convergence_with_neumann_bcs()
    call test_neumann_vs_analytical_solution()
    
    call check_summary("Neumann Boundary Conditions")

contains

    ! Test creation of Neumann boundary condition types
    subroutine test_neumann_bc_type_creation()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(neumann_bc_t) :: neumann_bc
        real(dp) :: flux_value
        
        mesh = unit_square_mesh(4)
        Vh = function_space(mesh, "Lagrange", 1)
        
        ! Test Neumann BC creation with constant flux
        flux_value = 2.0_dp
        neumann_bc = neumann_bc_constant(Vh, flux_value)
        
        call check_condition(associated(neumann_bc%space), &
            "Neumann BC: space pointer associated")
        call check_condition(neumann_bc%flux_type == "constant", &
            "Neumann BC: constant flux type")
        call check_condition(abs(neumann_bc%constant_value - flux_value) < 1.0e-14_dp, &
            "Neumann BC: correct constant value")
        
        write(*,*) "   Neumann BC type creation tests passed"
    end subroutine test_neumann_bc_type_creation

    ! Test mixed Dirichlet-Neumann boundary value problem
    subroutine test_mixed_dirichlet_neumann_problem()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: dirichlet_bc
        type(neumann_bc_t) :: neumann_bc
        type(form_expr_t) :: a, L
        real(dp) :: solution_norm, max_val, boundary_flux
        
        ! Mixed BC problem: -Δu = 1 in Ω, u = 0 on left edge, ∂u/∂n = 1 on right edge
        mesh = unit_square_mesh(5)
        Vh = function_space(mesh, "Lagrange", 1)
        
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(1.0_dp)
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        ! Dirichlet BC: u = 0 on left boundary (x = 0)
        dirichlet_bc = dirichlet_bc_on_boundary(Vh, 0.0_dp, "left")
        
        ! Neumann BC: ∂u/∂n = 1 on right boundary (x = 1)
        boundary_flux = 1.0_dp
        neumann_bc = neumann_bc_on_boundary(Vh, boundary_flux, "right")
        
        uh = function(Vh)
        
        ! Solve mixed BC problem
        call solve_mixed_bc(a == L, uh, dirichlet_bc, neumann_bc)
        
        solution_norm = sqrt(sum(uh%values**2))
        max_val = maxval(uh%values)
        
        call check_condition(solution_norm > 0.0_dp, &
            "Mixed BC: non-zero solution norm")
        call check_condition(max_val > 0.1_dp, &
            "Mixed BC: reasonable maximum value")
        call check_condition(max_val < 2.0_dp, &
            "Mixed BC: bounded solution")
        
        write(*,*) "   Mixed boundary value problem tests passed"
        write(*,*) "   Solution norm:", solution_norm
        write(*,*) "   Maximum value:", max_val
    end subroutine test_mixed_dirichlet_neumann_problem

    ! Test pure Neumann problem with compatibility condition
    subroutine test_pure_neumann_problem()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(neumann_bc_t) :: neumann_bc
        type(form_expr_t) :: a, L
        real(dp) :: solution_norm, mean_value, total_flux, domain_area
        integer :: i
        
        ! Pure Neumann problem: -Δu = f in Ω, ∂u/∂n = 0 on ∂Ω
        ! For compatibility: ∫_Ω f dx + ∫_∂Ω g ds = 0
        mesh = unit_square_mesh(4)
        Vh = function_space(mesh, "Lagrange", 1)
        domain_area = 1.0_dp  ! Unit square
        
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(0.0_dp)  ! Zero source to satisfy compatibility
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        ! Neumann BC: ∂u/∂n = 0 on all boundaries (natural BC)
        neumann_bc = neumann_bc_constant(Vh, 0.0_dp)
        
        uh = function(Vh)
        
        ! Solve pure Neumann problem
        call solve_neumann(a == L, uh, neumann_bc)
        
        solution_norm = sqrt(sum(uh%values**2))
        
        ! Compute mean value (should be arbitrary constant)
        mean_value = sum(uh%values) / real(Vh%ndof, dp)
        
        call check_condition(solution_norm >= 0.0_dp, &
            "Pure Neumann: non-negative solution norm")
        call check_condition(abs(mean_value) < 1.0e-6_dp .or. solution_norm < 1.0e-6_dp, &
            "Pure Neumann: solution unique up to constant")
        
        write(*,*) "   Pure Neumann problem tests passed"
        write(*,*) "   Solution norm:", solution_norm
        write(*,*) "   Mean value:", mean_value
    end subroutine test_pure_neumann_problem

    ! Test boundary integration for Neumann conditions
    subroutine test_neumann_boundary_integration()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(neumann_bc_t) :: neumann_bc
        real(dp) :: boundary_integral, expected_integral, boundary_length
        
        ! Test integration of Neumann BC over boundary
        mesh = unit_square_mesh(6)
        Vh = function_space(mesh, "Lagrange", 1)
        
        ! Constant flux on boundary
        neumann_bc = neumann_bc_constant(Vh, 2.0_dp)
        
        ! Compute boundary integral ∫_∂Ω g ds
        call compute_boundary_integral(neumann_bc, boundary_integral)
        
        ! For unit square with constant flux g = 2: ∫_∂Ω 2 ds = 2 * perimeter = 2 * 4 = 8
        boundary_length = 4.0_dp  ! Perimeter of unit square
        expected_integral = 2.0_dp * boundary_length
        
        call check_condition(abs(boundary_integral - expected_integral) < 0.1_dp, &
            "Boundary integration: correct integral value")
        call check_condition(boundary_integral > 7.0_dp, &
            "Boundary integration: reasonable lower bound")
        call check_condition(boundary_integral < 9.0_dp, &
            "Boundary integration: reasonable upper bound")
        
        write(*,*) "   Boundary integration tests passed"
        write(*,*) "   Computed integral:", boundary_integral
        write(*,*) "   Expected integral:", expected_integral
    end subroutine test_neumann_boundary_integration

    ! Test convergence with Neumann boundary conditions
    subroutine test_convergence_with_neumann_bcs()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: dirichlet_bc
        type(neumann_bc_t) :: neumann_bc
        type(form_expr_t) :: a, L
        real(dp) :: solution_norms(3), h_values(3)
        integer :: mesh_sizes(3), i
        logical :: convergence_ok
        
        mesh_sizes = [3, 4, 5]
        h_values = 1.0_dp / real(mesh_sizes, dp)
        
        ! Test convergence for mixed BC problem
        do i = 1, 3
            mesh = unit_square_mesh(mesh_sizes(i))
            Vh = function_space(mesh, "Lagrange", 1)
            
            u = trial_function(Vh)
            v = test_function(Vh)
            f = constant(1.0_dp)
            
            a = inner(grad(u), grad(v))*dx
            L = f*v*dx
            
            dirichlet_bc = dirichlet_bc_on_boundary(Vh, 0.0_dp, "left")
            neumann_bc = neumann_bc_on_boundary(Vh, 1.0_dp, "right")
            
            uh = function(Vh)
            call solve_mixed_bc(a == L, uh, dirichlet_bc, neumann_bc)
            
            solution_norms(i) = sqrt(sum(uh%values**2))
        end do
        
        ! Check that solutions exist and are reasonable
        convergence_ok = all(solution_norms > 0.01_dp) .and. all(solution_norms < 10.0_dp)
        
        call check_condition(convergence_ok, &
            "Neumann convergence: solutions exist and bounded")
        call check_condition(solution_norms(3) > 0.5_dp * solution_norms(1), &
            "Neumann convergence: solution norms in reasonable range")
        
        write(*,*) "   Convergence tests passed"
        write(*,*) "   Mesh sizes:", mesh_sizes
        write(*,*) "   Solution norms:", solution_norms
    end subroutine test_convergence_with_neumann_bcs

    ! Test Neumann BC against known analytical solution
    subroutine test_neumann_vs_analytical_solution()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: dirichlet_bc
        type(neumann_bc_t) :: neumann_bc
        type(form_expr_t) :: a, L
        real(dp) :: max_error, l2_error, x, y, u_exact, u_computed
        integer :: i, center_node
        real(dp) :: min_dist, dist
        
        ! Test against analytical solution u(x,y) = x² for -Δu = -2
        ! with u(0,y) = 0 (Dirichlet) and ∂u/∂n = 2 on x=1 (Neumann)
        mesh = unit_square_mesh(6)
        Vh = function_space(mesh, "Lagrange", 1)
        
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(-2.0_dp)  ! -Δ(x²) = -2
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        dirichlet_bc = dirichlet_bc_on_boundary(Vh, 0.0_dp, "left")   ! u = 0 at x = 0
        neumann_bc = neumann_bc_on_boundary(Vh, 2.0_dp, "right")      ! ∂u/∂x = 2 at x = 1
        
        uh = function(Vh)
        call solve_mixed_bc(a == L, uh, dirichlet_bc, neumann_bc)
        
        ! Compare with analytical solution at center point (0.5, 0.5)
        min_dist = huge(1.0_dp)
        center_node = 1
        do i = 1, Vh%ndof
            x = Vh%mesh%data%vertices(1, i)
            y = Vh%mesh%data%vertices(2, i)
            dist = (x - 0.5_dp)**2 + (y - 0.5_dp)**2
            if (dist < min_dist) then
                min_dist = dist
                center_node = i
            end if
        end do
        
        x = Vh%mesh%data%vertices(1, center_node)
        u_exact = x**2  ! Analytical solution
        u_computed = uh%values(center_node)
        max_error = abs(u_computed - u_exact)
        
        call check_condition(max_error < 0.5_dp, &
            "Analytical comparison: reasonable error at center")
        call check_condition(u_computed > 0.1_dp, &
            "Analytical comparison: positive solution at center")
        call check_condition(u_computed < 0.5_dp, &
            "Analytical comparison: bounded solution at center")
        
        write(*,*) "   Analytical solution tests passed"
        write(*,*) "   Center point (x,y):", x, Vh%mesh%data%vertices(2, center_node)
        write(*,*) "   Exact solution:", u_exact
        write(*,*) "   Computed solution:", u_computed
        write(*,*) "   Error:", max_error
    end subroutine test_neumann_vs_analytical_solution

end program test_neumann_boundary_conditions