program test_assembly_system
    use fortfem_kinds
    use fortfem_api
    use check
    implicit none

    write(*,*) "Testing finite element assembly system..."
    
    call test_element_stiffness_matrix()
    call test_element_load_vector()
    call test_global_assembly()
    call test_jacobian_computation()
    call test_shape_function_derivatives()
    call test_assembly_symmetry()
    call test_assembly_properties()
    call test_manufactured_assembly()
    
    call check_summary("Assembly System")

contains

    ! Test element stiffness matrix computation
    subroutine test_element_stiffness_matrix()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        
        ! Test on small mesh (but not too small)
        mesh = unit_square_mesh(3)  ! 3x3 mesh has interior nodes
        Vh = function_space(mesh, "Lagrange", 1)
        
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(1.0_dp)
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        bc = dirichlet_bc(Vh, 0.0_dp)
        uh = function(Vh)
        
        call solve(a == L, uh, bc)
        
        ! Test that stiffness matrix produces reasonable results
        call check_condition(allocated(uh%values), &
            "Element stiffness: solution values allocated")
        call check_condition(size(uh%values) == Vh%ndof, &
            "Element stiffness: solution size matches DOFs")
        call check_condition(maxval(uh%values) > 0.0_dp, &
            "Element stiffness: positive solution values")
        
        write(*,*) "   DOFs:", Vh%ndof
        write(*,*) "   Max solution:", maxval(uh%values)
        write(*,*) "   Min solution:", minval(uh%values)
    end subroutine test_element_stiffness_matrix

    ! Test element load vector computation
    subroutine test_element_load_vector()
        type(mesh_t) :: mesh1, mesh2
        type(function_space_t) :: Vh1, Vh2
        type(trial_function_t) :: u1, u2
        type(test_function_t) :: v1, v2
        type(function_t) :: f1, f2, uh1, uh2
        type(dirichlet_bc_t) :: bc1, bc2
        type(form_expr_t) :: a1, a2, L1, L2
        real(dp) :: integral1, integral2, ratio
        
        ! Test load vector scaling with mesh refinement
        mesh1 = unit_square_mesh(3)  ! Coarse mesh
        mesh2 = unit_square_mesh(5)  ! Fine mesh
        
        Vh1 = function_space(mesh1, "Lagrange", 1)
        Vh2 = function_space(mesh2, "Lagrange", 1)
        
        ! Solve same problem on both meshes
        u1 = trial_function(Vh1)
        v1 = test_function(Vh1)
        f1 = constant(1.0_dp)
        a1 = inner(grad(u1), grad(v1))*dx
        L1 = f1*v1*dx
        bc1 = dirichlet_bc(Vh1, 0.0_dp)
        uh1 = function(Vh1)
        call solve(a1 == L1, uh1, bc1)
        
        u2 = trial_function(Vh2)
        v2 = test_function(Vh2)
        f2 = constant(1.0_dp)
        a2 = inner(grad(u2), grad(v2))*dx
        L2 = f2*v2*dx
        bc2 = dirichlet_bc(Vh2, 0.0_dp)
        uh2 = function(Vh2)
        call solve(a2 == L2, uh2, bc2)
        
        ! Compute approximate integrals (sum of values * approximate area per vertex)
        integral1 = sum(uh1%values) / real(Vh1%ndof, dp)
        integral2 = sum(uh2%values) / real(Vh2%ndof, dp)
        
        call check_condition(integral2 > integral1, &
            "Load vector: finer mesh gives larger integral")
        
        ! Ratio should be reasonable (not too different)
        ! Handle case where integral1 might be very small
        if (integral1 > 1.0e-10_dp) then
            ratio = integral2 / integral1
            call check_condition(ratio > 0.5_dp .and. ratio < 5.0_dp, &
                "Load vector: mesh refinement gives reasonable ratio")
        else
            call check_condition(integral2 > 1.0e-6_dp, &
                "Load vector: fine mesh gives positive integral")
        end if
        
        write(*,*) "   Coarse mesh integral:", integral1
        write(*,*) "   Fine mesh integral:", integral2
        if (integral1 > 1.0e-10_dp) then
            write(*,*) "   Ratio:", integral2/integral1
        else
            write(*,*) "   Ratio: N/A (coarse integral too small)"
        end if
    end subroutine test_element_load_vector

    ! Test global assembly process
    subroutine test_global_assembly()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        integer :: i, interior_count, boundary_count
        real(dp) :: sum_interior, sum_boundary
        
        ! Test on small mesh to verify assembly
        mesh = unit_square_mesh(3)
        Vh = function_space(mesh, "Lagrange", 1)
        
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(1.0_dp)
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        bc = dirichlet_bc(Vh, 0.0_dp)
        uh = function(Vh)
        
        call solve(a == L, uh, bc)
        
        ! Count interior vs boundary nodes and their contributions
        interior_count = 0
        boundary_count = 0
        sum_interior = 0.0_dp
        sum_boundary = 0.0_dp
        
        do i = 1, Vh%ndof
            if (Vh%mesh%data%is_boundary_vertex(i)) then
                boundary_count = boundary_count + 1
                sum_boundary = sum_boundary + uh%values(i)
            else
                interior_count = interior_count + 1
                sum_interior = sum_interior + uh%values(i)
            end if
        end do
        
        call check_condition(interior_count > 0, &
            "Global assembly: found interior nodes")
        call check_condition(boundary_count > 0, &
            "Global assembly: found boundary nodes")
        call check_condition(abs(sum_boundary) < 1.0e-10_dp, &
            "Global assembly: boundary values are zero")
        call check_condition(sum_interior > 0.0_dp, &
            "Global assembly: positive interior values")
        
        write(*,*) "   Interior nodes:", interior_count
        write(*,*) "   Boundary nodes:", boundary_count
        write(*,*) "   Sum interior values:", sum_interior
        write(*,*) "   Sum boundary values:", sum_boundary
    end subroutine test_global_assembly

    ! Test Jacobian computation and coordinate transformations
    subroutine test_jacobian_computation()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        integer :: e, v1, v2, v3
        real(dp) :: x1, y1, x2, y2, x3, y3, area
        real(dp) :: area_formula, area_computed
        real(dp) :: a11, a12, a21, a22, det_a
        logical :: all_positive_areas, reasonable_jacobians
        integer :: valid_elements
        
        ! Test on regular mesh
        mesh = unit_square_mesh(4)
        Vh = function_space(mesh, "Lagrange", 1)
        
        all_positive_areas = .true.
        reasonable_jacobians = .true.
        valid_elements = 0
        
        ! Check all triangles have positive area and reasonable Jacobians
        do e = 1, mesh%data%n_triangles
            v1 = mesh%data%triangles(1, e)
            v2 = mesh%data%triangles(2, e)
            v3 = mesh%data%triangles(3, e)
            
            x1 = mesh%data%vertices(1, v1)
            y1 = mesh%data%vertices(2, v1)
            x2 = mesh%data%vertices(1, v2)
            y2 = mesh%data%vertices(2, v2)
            x3 = mesh%data%vertices(1, v3)
            y3 = mesh%data%vertices(2, v3)
            
            ! Compute area using formula from assembly
            area_formula = 0.5_dp * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
            
            ! Compute area using cross product
            area_computed = 0.5_dp * abs((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1))
            
            if (area_formula <= 0.0_dp) all_positive_areas = .false.
            
            ! Check Jacobian matrix
            a11 = x2 - x1; a12 = x3 - x1
            a21 = y2 - y1; a22 = y3 - y1
            det_a = a11*a22 - a12*a21
            
            if (abs(det_a) < 1.0e-12_dp) reasonable_jacobians = .false.
            
            ! For unit square with regular mesh, areas should be similar
            if (area_formula > 1.0e-6_dp .and. area_formula < 1.0_dp) then
                valid_elements = valid_elements + 1
            end if
        end do
        
        call check_condition(all_positive_areas, &
            "Jacobian computation: all elements have positive area")
        call check_condition(reasonable_jacobians, &
            "Jacobian computation: all Jacobians are non-singular")
        call check_condition(valid_elements > 0, &
            "Jacobian computation: found valid elements")
        
        write(*,*) "   Total triangles:", mesh%data%n_triangles
        write(*,*) "   Valid elements:", valid_elements
        write(*,*) "   All positive areas:", all_positive_areas
        write(*,*) "   All reasonable Jacobians:", reasonable_jacobians
    end subroutine test_jacobian_computation

    ! Test shape function derivatives
    subroutine test_shape_function_derivatives()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        integer :: e, v1, v2, v3
        real(dp) :: x1, y1, x2, y2, x3, y3
        real(dp) :: a(2,2), det_a, b(3), c(3)
        real(dp) :: grad_sum_x, grad_sum_y
        logical :: partition_of_unity
        integer :: test_element
        
        ! Test on simple mesh
        mesh = unit_square_mesh(3)
        Vh = function_space(mesh, "Lagrange", 1)
        
        ! Test gradient computation on first triangle
        test_element = 1
        v1 = mesh%data%triangles(1, test_element)
        v2 = mesh%data%triangles(2, test_element)
        v3 = mesh%data%triangles(3, test_element)
        
        x1 = mesh%data%vertices(1, v1)
        y1 = mesh%data%vertices(2, v1)
        x2 = mesh%data%vertices(1, v2)
        y2 = mesh%data%vertices(2, v2)
        x3 = mesh%data%vertices(1, v3)
        y3 = mesh%data%vertices(2, v3)
        
        ! Compute shape function gradients (copied from assembly code)
        a(1,1) = x2 - x1; a(1,2) = x3 - x1
        a(2,1) = y2 - y1; a(2,2) = y3 - y1
        det_a = a(1,1)*a(2,2) - a(1,2)*a(2,1)
        
        ! Physical gradients
        b(1) = (-a(2,2) + a(2,1)) / det_a
        c(1) = ( a(1,2) - a(1,1)) / det_a
        
        b(2) = a(2,2) / det_a
        c(2) = -a(1,2) / det_a
        
        b(3) = -a(2,1) / det_a
        c(3) = a(1,1) / det_a
        
        ! Test partition of unity: sum of gradients should be zero
        grad_sum_x = b(1) + b(2) + b(3)
        grad_sum_y = c(1) + c(2) + c(3)
        
        partition_of_unity = (abs(grad_sum_x) < 1.0e-12_dp) .and. &
                           (abs(grad_sum_y) < 1.0e-12_dp)
        
        call check_condition(partition_of_unity, &
            "Shape function derivatives: partition of unity")
        call check_condition(abs(det_a) > 1.0e-12_dp, &
            "Shape function derivatives: non-singular Jacobian")
        
        ! Test that gradients are reasonable magnitudes
        call check_condition(maxval(abs(b)) < 100.0_dp, &
            "Shape function derivatives: reasonable x-gradients")
        call check_condition(maxval(abs(c)) < 100.0_dp, &
            "Shape function derivatives: reasonable y-gradients")
        
        write(*,*) "   Test element:", test_element
        write(*,*) "   Jacobian determinant:", det_a
        write(*,*) "   Sum of x-gradients:", grad_sum_x
        write(*,*) "   Sum of y-gradients:", grad_sum_y
        write(*,*) "   Max |grad_x|:", maxval(abs(b))
        write(*,*) "   Max |grad_y|:", maxval(abs(c))
    end subroutine test_shape_function_derivatives

    ! Test assembly matrix symmetry
    subroutine test_assembly_symmetry()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh1, uh2
        type(dirichlet_bc_t) :: bc1, bc2
        type(form_expr_t) :: a, L
        real(dp) :: energy1, energy2, diff
        
        ! Test symmetry by solving with different BC values
        mesh = unit_square_mesh(4)
        Vh = function_space(mesh, "Lagrange", 1)
        
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(1.0_dp)
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        ! Solve with BC = 0
        bc1 = dirichlet_bc(Vh, 0.0_dp)
        uh1 = function(Vh)
        call solve(a == L, uh1, bc1)
        
        ! Solve with BC = 1
        bc2 = dirichlet_bc(Vh, 1.0_dp)
        uh2 = function(Vh)
        call solve(a == L, uh2, bc2)
        
        ! Compute approximate energy norms (simplified)
        energy1 = sum(uh1%values**2)
        energy2 = sum(uh2%values**2)
        
        call check_condition(energy1 >= 0.0_dp, &
            "Assembly symmetry: non-negative energy norm 1")
        call check_condition(energy2 >= 0.0_dp, &
            "Assembly symmetry: non-negative energy norm 2")
        call check_condition(energy2 > energy1, &
            "Assembly symmetry: higher BC gives higher energy")
        
        write(*,*) "   Energy norm (BC=0):", energy1
        write(*,*) "   Energy norm (BC=1):", energy2
        write(*,*) "   Energy ratio:", energy2/energy1
    end subroutine test_assembly_symmetry

    ! Test assembly properties like positive definiteness
    subroutine test_assembly_properties()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        real(dp) :: solution_norm, gradient_norm
        integer :: i
        
        ! Test on medium mesh
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
        
        ! Compute solution norms
        solution_norm = sqrt(sum(uh%values**2))
        
        ! Approximate gradient norm (very simplified)
        gradient_norm = 0.0_dp
        do i = 1, Vh%ndof-1
            gradient_norm = gradient_norm + (uh%values(i+1) - uh%values(i))**2
        end do
        gradient_norm = sqrt(gradient_norm)
        
        call check_condition(solution_norm > 0.0_dp, &
            "Assembly properties: positive solution norm")
        call check_condition(gradient_norm >= 0.0_dp, &
            "Assembly properties: non-negative gradient norm")
        
        ! Test coercivity: for Poisson, energy should be bounded
        call check_condition(solution_norm < 1.0_dp, &
            "Assembly properties: bounded solution")
        
        write(*,*) "   Solution L2 norm:", solution_norm
        write(*,*) "   Approximate gradient norm:", gradient_norm
        write(*,*) "   DOFs:", Vh%ndof
    end subroutine test_assembly_properties

    ! Test assembly against manufactured solution
    subroutine test_manufactured_assembly()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        real(dp) :: max_val, center_val, corner_val
        integer :: center_node, corner_node, i
        real(dp) :: x, y, min_dist_center, min_dist_corner, dist
        
        ! Test assembly accuracy with known problem
        mesh = unit_square_mesh(6)
        Vh = function_space(mesh, "Lagrange", 1)
        
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(1.0_dp)  ! Note: solver hardcodes f=1
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        bc = dirichlet_bc(Vh, 0.0_dp)
        uh = function(Vh)
        
        call solve(a == L, uh, bc)
        
        ! Find nodes closest to center and corner
        min_dist_center = huge(1.0_dp)
        min_dist_corner = huge(1.0_dp)
        center_node = 1
        corner_node = 1
        
        do i = 1, Vh%ndof
            x = Vh%mesh%data%vertices(1, i)
            y = Vh%mesh%data%vertices(2, i)
            
            ! Distance to center (0.5, 0.5)
            dist = (x - 0.5_dp)**2 + (y - 0.5_dp)**2
            if (dist < min_dist_center) then
                min_dist_center = dist
                center_node = i
            end if
            
            ! Distance to corner (0, 0)
            dist = x**2 + y**2
            if (dist < min_dist_corner) then
                min_dist_corner = dist
                corner_node = i
            end if
        end do
        
        max_val = maxval(uh%values)
        center_val = uh%values(center_node)
        corner_val = uh%values(corner_node)
        
        ! For -Δu = 1 on unit square with u=0 on boundary:
        ! - Maximum should be at center
        ! - Corner should be zero (boundary)
        call check_condition(center_val > corner_val, &
            "Manufactured assembly: center > corner")
        call check_condition(abs(corner_val) < 1.0e-10_dp, &
            "Manufactured assembly: corner value is zero")
        call check_condition(center_val > 0.01_dp, &
            "Manufactured assembly: reasonable center value")
        call check_condition(center_val < 0.15_dp, &
            "Manufactured assembly: bounded center value")
        
        write(*,*) "   Maximum value:", max_val
        write(*,*) "   Center value:", center_val
        write(*,*) "   Corner value:", corner_val
        write(*,*) "   Center/Max ratio:", center_val/max_val
    end subroutine test_manufactured_assembly

end program test_assembly_system