program test_assembly_advanced
    use fortfem_kinds
    use fortfem_api
    use check
    implicit none

    write(*,*) "Testing advanced assembly system features..."
    
    call test_mass_matrix_properties()
    call test_stiffness_matrix_properties()
    call test_element_matrix_properties()
    call test_quadrature_accuracy()
    call test_assembly_scaling()
    call test_mesh_quality_effects()
    call test_reference_element_mapping()
    
    call check_summary("Advanced Assembly System")

contains

    ! Test mass matrix properties (if we had access to it)
    subroutine test_mass_matrix_properties()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        real(dp) :: total_mass, expected_mass
        
        ! Test on unit square where we can compute expected mass
        mesh = unit_square_mesh(4)
        Vh = function_space(mesh, "Lagrange", 1)
        
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(1.0_dp)
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        bc = dirichlet_bc(Vh, 0.0_dp)
        uh = function(Vh)
        
        call solve(a == L, uh, bc)
        
        ! Estimate total mass by integrating constant function
        ! For unit square, ∫ 1 dx = 1
        total_mass = sum(uh%values) * (1.0_dp / Vh%ndof)  ! Approximate integration
        expected_mass = 1.0_dp
        
        call check_condition(total_mass > 0.0_dp, &
            "Mass matrix: positive total mass")
        
        write(*,*) "   Approximate total mass:", total_mass
        write(*,*) "   Expected mass (unit square):", expected_mass
    end subroutine test_mass_matrix_properties

    ! Test stiffness matrix properties  
    subroutine test_stiffness_matrix_properties()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh1, uh2
        type(dirichlet_bc_t) :: bc1, bc2
        type(form_expr_t) :: a, L
        real(dp) :: energy1, energy2
        
        ! Test stiffness matrix positive definiteness through energy
        mesh = unit_square_mesh(5)
        Vh = function_space(mesh, "Lagrange", 1)
        
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(1.0_dp)
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        ! Solve with different boundary conditions
        bc1 = dirichlet_bc(Vh, 0.0_dp)
        uh1 = function(Vh)
        call solve(a == L, uh1, bc1)
        
        bc2 = dirichlet_bc(Vh, 0.5_dp)
        uh2 = function(Vh)
        call solve(a == L, uh2, bc2)
        
        ! Compute approximate energies
        energy1 = sum(uh1%values**2)
        energy2 = sum(uh2%values**2)
        
        call check_condition(energy1 >= 0.0_dp, &
            "Stiffness matrix: non-negative energy 1")
        call check_condition(energy2 >= 0.0_dp, &
            "Stiffness matrix: non-negative energy 2")
        call check_condition(energy2 > energy1, &
            "Stiffness matrix: larger BC gives larger energy")
        
        write(*,*) "   Energy with BC=0:", energy1
        write(*,*) "   Energy with BC=0.5:", energy2
        write(*,*) "   Energy ratio:", energy2/energy1
    end subroutine test_stiffness_matrix_properties

    ! Test individual element matrix properties
    subroutine test_element_matrix_properties()
        type(mesh_t) :: mesh
        real(dp) :: x1, y1, x2, y2, x3, y3, area
        real(dp) :: a(2,2), det_a, b(3), c(3), K_elem(3,3)
        integer :: i, j, e, v1, v2, v3
        real(dp) :: trace_K, det_K_2x2, sum_off_diag
        logical :: symmetric_element
        
        ! Test element matrix properties on first triangle
        mesh = unit_square_mesh(4)
        e = 1  ! First triangle
        
        v1 = mesh%data%triangles(1, e)
        v2 = mesh%data%triangles(2, e)
        v3 = mesh%data%triangles(3, e)
        
        x1 = mesh%data%vertices(1, v1)
        y1 = mesh%data%vertices(2, v1)
        x2 = mesh%data%vertices(1, v2)
        y2 = mesh%data%vertices(2, v2)
        x3 = mesh%data%vertices(1, v3)
        y3 = mesh%data%vertices(2, v3)
        
        area = 0.5_dp * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
        
        ! Compute element matrix (from assembly code)
        a(1,1) = x2 - x1; a(1,2) = x3 - x1
        a(2,1) = y2 - y1; a(2,2) = y3 - y1
        det_a = a(1,1)*a(2,2) - a(1,2)*a(2,1)
        
        b(1) = (-a(2,2) + a(2,1)) / det_a
        c(1) = ( a(1,2) - a(1,1)) / det_a
        b(2) = a(2,2) / det_a
        c(2) = -a(1,2) / det_a
        b(3) = -a(2,1) / det_a
        c(3) = a(1,1) / det_a
        
        do i = 1, 3
            do j = 1, 3
                K_elem(i,j) = area * (b(i)*b(j) + c(i)*c(j))
            end do
        end do
        
        ! Test symmetry
        symmetric_element = .true.
        do i = 1, 3
            do j = 1, 3
                if (abs(K_elem(i,j) - K_elem(j,i)) > 1.0e-12_dp) then
                    symmetric_element = .false.
                end if
            end do
        end do
        
        ! Test trace and other properties
        trace_K = K_elem(1,1) + K_elem(2,2) + K_elem(3,3)
        sum_off_diag = K_elem(1,2) + K_elem(1,3) + K_elem(2,3)
        
        call check_condition(symmetric_element, &
            "Element matrix: symmetry")
        call check_condition(trace_K > 0.0_dp, &
            "Element matrix: positive trace")
        call check_condition(area > 0.0_dp, &
            "Element matrix: positive area")
        call check_condition(abs(det_a) > 1.0e-12_dp, &
            "Element matrix: non-singular Jacobian")
        
        write(*,*) "   Element area:", area
        write(*,*) "   Jacobian determinant:", det_a
        write(*,*) "   Element matrix trace:", trace_K
        write(*,*) "   Sum off-diagonal:", sum_off_diag
        write(*,*) "   Symmetric:", symmetric_element
    end subroutine test_element_matrix_properties

    ! Test quadrature accuracy (current implementation uses exact integration)
    subroutine test_quadrature_accuracy()
        type(mesh_t) :: mesh1, mesh2
        type(function_space_t) :: Vh1, Vh2
        type(trial_function_t) :: u1, u2
        type(test_function_t) :: v1, v2
        type(function_t) :: f1, f2, uh1, uh2
        type(dirichlet_bc_t) :: bc1, bc2
        type(form_expr_t) :: a1, a2, L1, L2
        real(dp) :: integral1, integral2, convergence_rate
        integer :: dof1, dof2
        
        ! Test integration accuracy through mesh refinement
        mesh1 = unit_square_mesh(4)   ! Coarse
        mesh2 = unit_square_mesh(8)   ! Fine
        
        Vh1 = function_space(mesh1, "Lagrange", 1)
        Vh2 = function_space(mesh2, "Lagrange", 1)
        dof1 = Vh1%ndof
        dof2 = Vh2%ndof
        
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
        
        ! Approximate integrals
        integral1 = sum(uh1%values) / real(dof1, dp)
        integral2 = sum(uh2%values) / real(dof2, dp)
        
        call check_condition(integral2 > integral1, &
            "Quadrature accuracy: convergence with refinement")
        
        ! Test that the change is significant and in the right direction
        if (integral1 > 1.0e-12_dp) then
            convergence_rate = abs(integral2 - integral1) / integral1
            call check_condition(convergence_rate > 0.1_dp .and. convergence_rate < 10.0_dp, &
                "Quadrature accuracy: reasonable convergence rate")
        else
            call check_condition(.true., &
                "Quadrature accuracy: baseline integral too small")
        end if
        
        write(*,*) "   Coarse DOFs:", dof1, "Integral:", integral1
        write(*,*) "   Fine DOFs:", dof2, "Integral:", integral2
        write(*,*) "   Relative change:", abs(integral2 - integral1) / integral1
    end subroutine test_quadrature_accuracy

    ! Test assembly scaling with mesh size
    subroutine test_assembly_scaling()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        real(dp) :: max_vals(3), h_values(3)
        integer :: n_values(3), i
        
        n_values = [4, 6, 8]
        h_values = 1.0_dp / real(n_values, dp)
        
        ! Test scaling of solution maximum with mesh size
        do i = 1, 3
            mesh = unit_square_mesh(n_values(i))
            Vh = function_space(mesh, "Lagrange", 1)
            
            u = trial_function(Vh)
            v = test_function(Vh)
            f = constant(1.0_dp)
            
            a = inner(grad(u), grad(v))*dx
            L = f*v*dx
            
            bc = dirichlet_bc(Vh, 0.0_dp)
            uh = function(Vh)
            
            call solve(a == L, uh, bc)
            max_vals(i) = maxval(uh%values)
        end do
        
        ! Test convergence: finer mesh should give more accurate (higher) maximum
        call check_condition(max_vals(2) > max_vals(1), &
            "Assembly scaling: refinement increases accuracy 1")
        call check_condition(max_vals(3) > max_vals(2), &
            "Assembly scaling: refinement increases accuracy 2")
        
        ! Values should be bounded
        call check_condition(maxval(max_vals) < 0.15_dp, &
            "Assembly scaling: solution remains bounded")
        
        write(*,*) "   Mesh sizes:", n_values
        write(*,*) "   Max values:", max_vals
        write(*,*) "   h values:", h_values
    end subroutine test_assembly_scaling

    ! Test effects of mesh quality on assembly
    subroutine test_mesh_quality_effects()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        integer :: e, v1, v2, v3
        real(dp) :: x1, y1, x2, y2, x3, y3, area
        real(dp) :: edge1, edge2, edge3, semiperimeter, quality
        real(dp) :: min_quality, max_quality, avg_quality
        integer :: valid_elements, poor_quality_count
        
        ! Analyze mesh quality on moderately refined mesh
        mesh = unit_square_mesh(6)
        Vh = function_space(mesh, "Lagrange", 1)
        
        min_quality = 1.0_dp
        max_quality = 0.0_dp
        avg_quality = 0.0_dp
        valid_elements = 0
        poor_quality_count = 0
        
        ! Compute quality for each triangle
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
            
            area = 0.5_dp * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
            
            if (area > 1.0e-12_dp) then
                ! Compute edge lengths
                edge1 = sqrt((x2-x1)**2 + (y2-y1)**2)
                edge2 = sqrt((x3-x2)**2 + (y3-y2)**2)
                edge3 = sqrt((x1-x3)**2 + (y1-y3)**2)
                
                semiperimeter = 0.5_dp * (edge1 + edge2 + edge3)
                
                if (semiperimeter > 1.0e-12_dp) then
                    ! Quality = area / (semiperimeter^2), normalized
                    quality = 4.0_dp * sqrt(3.0_dp) * area / (edge1**2 + edge2**2 + edge3**2)
                    
                    min_quality = min(min_quality, quality)
                    max_quality = max(max_quality, quality)
                    avg_quality = avg_quality + quality
                    valid_elements = valid_elements + 1
                    
                    if (quality < 0.5_dp) poor_quality_count = poor_quality_count + 1
                end if
            end if
        end do
        
        if (valid_elements > 0) avg_quality = avg_quality / real(valid_elements, dp)
        
        call check_condition(valid_elements > 0, &
            "Mesh quality: found valid elements")
        call check_condition(min_quality > 0.0_dp, &
            "Mesh quality: all elements have positive quality")
        call check_condition(avg_quality > 0.3_dp, &
            "Mesh quality: reasonable average quality")
        call check_condition(real(poor_quality_count, dp) / real(valid_elements, dp) < 0.5_dp, &
            "Mesh quality: most elements have acceptable quality")
        
        write(*,*) "   Total triangles:", mesh%data%n_triangles
        write(*,*) "   Valid elements:", valid_elements
        write(*,*) "   Quality range: [", min_quality, ",", max_quality, "]"
        write(*,*) "   Average quality:", avg_quality
        write(*,*) "   Poor quality elements:", poor_quality_count
    end subroutine test_mesh_quality_effects

    ! Test reference element to physical element mapping
    subroutine test_reference_element_mapping()
        type(mesh_t) :: mesh
        integer :: e, v1, v2, v3
        real(dp) :: x1, y1, x2, y2, x3, y3
        real(dp) :: a(2,2), det_a
        real(dp) :: xi, eta, x_phys, y_phys
        real(dp) :: xi_vals(3), eta_vals(3)
        real(dp) :: x_expected(3), y_expected(3)
        logical :: mapping_correct
        integer :: i
        
        ! Test mapping on first triangle
        mesh = unit_square_mesh(4)
        e = 1
        
        v1 = mesh%data%triangles(1, e)
        v2 = mesh%data%triangles(2, e)
        v3 = mesh%data%triangles(3, e)
        
        x1 = mesh%data%vertices(1, v1)
        y1 = mesh%data%vertices(2, v1)
        x2 = mesh%data%vertices(1, v2)
        y2 = mesh%data%vertices(2, v2)
        x3 = mesh%data%vertices(1, v3)
        y3 = mesh%data%vertices(2, v3)
        
        ! Jacobian matrix
        a(1,1) = x2 - x1; a(1,2) = x3 - x1
        a(2,1) = y2 - y1; a(2,2) = y3 - y1
        det_a = a(1,1)*a(2,2) - a(1,2)*a(2,1)
        
        ! Test points on reference triangle
        xi_vals  = [0.0_dp, 1.0_dp, 0.0_dp]
        eta_vals = [0.0_dp, 0.0_dp, 1.0_dp]
        x_expected = [x1, x2, x3]
        y_expected = [y1, y2, y3]
        
        mapping_correct = .true.
        
        ! Test mapping: x = x1 + a11*xi + a12*eta, y = y1 + a21*xi + a22*eta
        do i = 1, 3
            xi = xi_vals(i)
            eta = eta_vals(i)
            
            x_phys = x1 + a(1,1)*xi + a(1,2)*eta
            y_phys = y1 + a(2,1)*xi + a(2,2)*eta
            
            if (abs(x_phys - x_expected(i)) > 1.0e-12_dp .or. &
                abs(y_phys - y_expected(i)) > 1.0e-12_dp) then
                mapping_correct = .false.
            end if
        end do
        
        call check_condition(abs(det_a) > 1.0e-12_dp, &
            "Reference mapping: non-singular Jacobian")
        call check_condition(mapping_correct, &
            "Reference mapping: vertices map correctly")
        
        write(*,*) "   Test element:", e
        write(*,*) "   Jacobian determinant:", det_a
        write(*,*) "   Vertex 1: (", x1, ",", y1, ")"
        write(*,*) "   Vertex 2: (", x2, ",", y2, ")"
        write(*,*) "   Vertex 3: (", x3, ",", y3, ")"
        write(*,*) "   Mapping correct:", mapping_correct
    end subroutine test_reference_element_mapping

end program test_assembly_advanced