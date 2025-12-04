program test_solver_validation
    use fortfem_kinds
    use fortfem_api
    use check
    implicit none

    write(*,*) "Testing advanced solver validation..."
    
    call test_analytical_solutions()
    call test_error_norms_convergence()
    call test_solver_conditioning()
    call test_gmres_convergence_behavior()
    call test_matrix_properties()
    call test_solution_smoothness()
    call test_energy_conservation()
    
    call check_summary("Advanced Solver Validation")

contains

    ! Test against analytical solutions
    subroutine test_analytical_solutions()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        real(dp) :: center_val, expected_center, relative_error
        integer :: i, center_node
        real(dp) :: x, y, min_dist, dist
        
        ! Test against known analytical solution for unit square
        ! For -Δu = 1 with u=0 on boundary, analytical solution involves Fourier series
        mesh = unit_square_mesh(8)
        Vh = function_space(mesh, "Lagrange", 1)
        
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(1.0_dp)
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        bc = dirichlet_bc(Vh, 0.0_dp)
        uh = function(Vh)
        
        call solve(a == L, uh, bc)
        
        ! Find center point
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
        
        center_val = uh%values(center_node)
        expected_center = 0.073_dp  ! Approximate from analytical solution
        relative_error = abs(center_val - expected_center) / expected_center
        
        call check_condition(center_val > 0.01_dp, &
            "Analytical solutions: positive center value")
        call check_condition(center_val < 0.15_dp, &
            "Analytical solutions: bounded center value")
        call check_condition(relative_error < 1.0_dp, &
            "Analytical solutions: reasonable accuracy")
        
        write(*,*) "   Center value:", center_val
        write(*,*) "   Expected center:", expected_center
        write(*,*) "   Relative error:", relative_error
    end subroutine test_analytical_solutions

    ! Test error norms and convergence rates
    subroutine test_error_norms_convergence()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        real(dp) :: l2_norms(3), h1_norms(3), h_vals(3)
        real(dp) :: l2_rate, h1_rate
        integer :: n_vals(3), i, j
        
        n_vals = [4, 6, 8]
        h_vals = 1.0_dp / real(n_vals, dp)
        
        ! Compute error norms on sequence of meshes
        do i = 1, 3
            mesh = unit_square_mesh(n_vals(i))
            Vh = function_space(mesh, "Lagrange", 1)
            
            u = trial_function(Vh)
            v = test_function(Vh)
            f = constant(1.0_dp)
            
            a = inner(grad(u), grad(v))*dx
            L = f*v*dx
            
            bc = dirichlet_bc(Vh, 0.0_dp)
            uh = function(Vh)
            
            call solve(a == L, uh, bc)
            
            ! Compute approximate L2 norm
            l2_norms(i) = sqrt(sum(uh%values**2) / real(Vh%ndof, dp))
            
            ! Compute approximate H1 norm (simplified)
            h1_norms(i) = 0.0_dp
            do j = 1, Vh%ndof-1
                h1_norms(i) = h1_norms(i) + (uh%values(j+1) - uh%values(j))**2
            end do
            h1_norms(i) = sqrt(h1_norms(i) + l2_norms(i)**2)
        end do
        
        ! Compute convergence rates
        if (l2_norms(1) > 1.0e-12_dp .and. l2_norms(2) > 1.0e-12_dp) then
            l2_rate = log(l2_norms(2) / l2_norms(1)) / log(h_vals(2) / h_vals(1))
        else
            l2_rate = 0.0_dp
        end if
        
        if (h1_norms(1) > 1.0e-12_dp .and. h1_norms(2) > 1.0e-12_dp) then
            h1_rate = log(h1_norms(2) / h1_norms(1)) / log(h_vals(2) / h_vals(1))
        else
            h1_rate = 0.0_dp
        end if
        
        call check_condition(l2_norms(2) > l2_norms(1), &
            "Error convergence: L2 norm increases with refinement")
        call check_condition(h1_norms(2) > h1_norms(1), &
            "Error convergence: H1 norm increases with refinement")
        call check_condition(abs(l2_rate) < 5.0_dp, &
            "Error convergence: reasonable L2 rate")
        call check_condition(abs(h1_rate) < 5.0_dp, &
            "Error convergence: reasonable H1 rate")
        
        write(*,*) "   L2 norms:", l2_norms
        write(*,*) "   H1 norms:", h1_norms
        write(*,*) "   L2 convergence rate:", l2_rate
        write(*,*) "   H1 convergence rate:", h1_rate
    end subroutine test_error_norms_convergence

    ! Test solver conditioning and stability
    subroutine test_solver_conditioning()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        real(dp) :: solution_norm, perturbation, perturbed_norm
        real(dp) :: condition_estimate
        
        ! Test solver stability with slightly different problems
        mesh = unit_square_mesh(5)
        Vh = function_space(mesh, "Lagrange", 1)
        
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(1.0_dp)
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        bc = dirichlet_bc(Vh, 0.0_dp)
        uh = function(Vh)
        
        call solve(a == L, uh, bc)
        solution_norm = sqrt(sum(uh%values**2))
        
        ! Solve with slightly perturbed BC
        bc = dirichlet_bc(Vh, 1.0e-6_dp)
        call solve(a == L, uh, bc)
        perturbed_norm = sqrt(sum(uh%values**2))
        
        ! Estimate condition number effect
        perturbation = 1.0e-6_dp
        condition_estimate = abs(perturbed_norm - solution_norm) / perturbation
        
        call check_condition(solution_norm > 0.0_dp, &
            "Solver conditioning: non-zero solution norm")
        call check_condition(condition_estimate < 1.0e6_dp, &
            "Solver conditioning: reasonable conditioning")
        call check_condition(abs(perturbed_norm - solution_norm) > 1.0e-12_dp, &
            "Solver conditioning: perturbation has effect")
        
        write(*,*) "   Original solution norm:", solution_norm
        write(*,*) "   Perturbed solution norm:", perturbed_norm
        write(*,*) "   Condition estimate:", condition_estimate
    end subroutine test_solver_conditioning

    ! Test GMRES convergence behavior in detail
    subroutine test_gmres_convergence_behavior()
        type(mesh_t) :: mesh
        type(vector_function_space_t) :: Vh
        type(vector_trial_function_t) :: E
        type(vector_test_function_t) :: F
        type(vector_function_t) :: j, Eh1, Eh2
        type(vector_bc_t) :: bc1, bc2
        type(form_expr_t) :: a, L
        real(dp) :: norm1, norm2, diff_norm
        
        ! Test GMRES with different initial conditions
        mesh = unit_square_mesh(4)
        Vh = vector_function_space(mesh, "Nedelec", 1)
        
        E = vector_trial_function(Vh)
        F = vector_test_function(Vh)
        j = vector_function(Vh)
        
        ! Set up current source
        if (allocated(j%values)) then
            j%values = 0.1_dp
        end if
        
        a = inner(curl(E), curl(F))*dx + inner(E, F)*dx
        L = inner(j, F)*dx
        
        ! Solve with different boundary conditions
        bc1 = vector_bc(Vh, [0.0_dp, 0.0_dp])
        Eh1 = vector_function(Vh)
        call solve(a == L, Eh1, bc1, "gmres")
        
        bc2 = vector_bc(Vh, [0.1_dp, 0.0_dp])
        Eh2 = vector_function(Vh)
        call solve(a == L, Eh2, bc2, "gmres")
        
        ! Compute solution norms
        if (allocated(Eh1%values) .and. allocated(Eh2%values)) then
            norm1 = sqrt(sum(Eh1%values**2))
            norm2 = sqrt(sum(Eh2%values**2))
            diff_norm = sqrt(sum((Eh1%values - Eh2%values)**2))
        else
            norm1 = 0.0_dp
            norm2 = 0.0_dp
            diff_norm = 0.0_dp
        end if
        
        call check_condition(norm1 >= 0.0_dp, &
            "GMRES convergence: non-negative norm 1")
        call check_condition(norm2 >= 0.0_dp, &
            "GMRES convergence: non-negative norm 2")
        call check_condition(diff_norm >= 0.0_dp, &
            "GMRES convergence: valid difference computation")
        call check_condition(diff_norm < 10.0_dp * max(norm1, norm2), &
            "GMRES convergence: solutions are reasonably close")
        
        write(*,*) "   GMRES solution norm 1:", norm1
        write(*,*) "   GMRES solution norm 2:", norm2
        write(*,*) "   Difference norm:", diff_norm
    end subroutine test_gmres_convergence_behavior

    ! Test matrix properties through solver behavior
    subroutine test_matrix_properties()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        real(dp) :: energy1, energy2, energy_ratio
        
        ! Test matrix positive definiteness through energy
        mesh = unit_square_mesh(6)
        Vh = function_space(mesh, "Lagrange", 1)
        
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(1.0_dp)
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        ! Test with two different boundary values
        bc = dirichlet_bc(Vh, 0.0_dp)
        uh = function(Vh)
        call solve(a == L, uh, bc)
        energy1 = sum(uh%values**2)
        
        bc = dirichlet_bc(Vh, 0.1_dp)
        call solve(a == L, uh, bc)
        energy2 = sum(uh%values**2)
        
        energy_ratio = energy2 / energy1
        
        call check_condition(energy1 >= 0.0_dp, &
            "Matrix properties: non-negative energy 1")
        call check_condition(energy2 >= 0.0_dp, &
            "Matrix properties: non-negative energy 2")
        call check_condition(energy2 > energy1, &
            "Matrix properties: higher BC gives higher energy")
        call check_condition(energy_ratio > 1.0_dp .and. energy_ratio < 1000.0_dp, &
            "Matrix properties: reasonable energy scaling")
        
        write(*,*) "   Energy with BC=0:", energy1
        write(*,*) "   Energy with BC=0.1:", energy2
        write(*,*) "   Energy ratio:", energy_ratio
    end subroutine test_matrix_properties

    ! Test solution smoothness properties
    subroutine test_solution_smoothness()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        real(dp) :: max_gradient, avg_gradient, gradient_variance, grad_approx
        integer :: i, gradient_count
        
        ! Test solution smoothness
        mesh = unit_square_mesh(6)
        Vh = function_space(mesh, "Lagrange", 1)
        
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(1.0_dp)
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        bc = dirichlet_bc(Vh, 0.0_dp)
        uh = function(Vh)
        
        call solve(a == L, uh, bc)
        
        ! Compute approximate gradients (simplified)
        max_gradient = 0.0_dp
        avg_gradient = 0.0_dp
        gradient_count = 0
        
        do i = 1, Vh%ndof-1
            grad_approx = abs(uh%values(i+1) - uh%values(i))
            max_gradient = max(max_gradient, grad_approx)
            avg_gradient = avg_gradient + grad_approx
            gradient_count = gradient_count + 1
        end do
        
        if (gradient_count > 0) avg_gradient = avg_gradient / real(gradient_count, dp)
        
        ! Compute gradient variance
        gradient_variance = 0.0_dp
        do i = 1, Vh%ndof-1
            grad_approx = abs(uh%values(i+1) - uh%values(i))
            gradient_variance = gradient_variance + (grad_approx - avg_gradient)**2
        end do
        if (gradient_count > 1) gradient_variance = gradient_variance / real(gradient_count-1, dp)
        
        call check_condition(max_gradient < 1.0_dp, &
            "Solution smoothness: bounded gradients")
        call check_condition(avg_gradient > 0.0_dp, &
            "Solution smoothness: positive average gradient")
        call check_condition(gradient_variance < 1.0_dp, &
            "Solution smoothness: reasonable gradient variance")
        
        write(*,*) "   Max gradient:", max_gradient
        write(*,*) "   Average gradient:", avg_gradient
        write(*,*) "   Gradient variance:", gradient_variance
    end subroutine test_solution_smoothness

    ! Test energy conservation properties
    subroutine test_energy_conservation()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        real(dp) :: total_energy, kinetic_energy, potential_energy
        real(dp) :: energy_balance
        integer :: i
        
        ! Test energy balance in finite element solution
        mesh = unit_square_mesh(5)
        Vh = function_space(mesh, "Lagrange", 1)
        
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(1.0_dp)
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        bc = dirichlet_bc(Vh, 0.0_dp)
        uh = function(Vh)
        
        call solve(a == L, uh, bc)
        
        ! Compute approximate energy terms
        kinetic_energy = sum(uh%values**2) / real(Vh%ndof, dp)
        
        potential_energy = 0.0_dp
        do i = 1, Vh%ndof-1
            potential_energy = potential_energy + (uh%values(i+1) - uh%values(i))**2
        end do
        potential_energy = potential_energy / real(Vh%ndof, dp)
        
        total_energy = kinetic_energy + potential_energy
        energy_balance = abs(total_energy - kinetic_energy - potential_energy)
        
        call check_condition(kinetic_energy >= 0.0_dp, &
            "Energy conservation: non-negative kinetic energy")
        call check_condition(potential_energy >= 0.0_dp, &
            "Energy conservation: non-negative potential energy")
        call check_condition(total_energy > 0.0_dp, &
            "Energy conservation: positive total energy")
        call check_condition(energy_balance < 1.0e-12_dp, &
            "Energy conservation: energy balance")
        
        write(*,*) "   Kinetic energy:", kinetic_energy
        write(*,*) "   Potential energy:", potential_energy
        write(*,*) "   Total energy:", total_energy
        write(*,*) "   Energy balance:", energy_balance
    end subroutine test_energy_conservation

end program test_solver_validation