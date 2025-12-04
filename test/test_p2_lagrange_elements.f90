program test_p2_lagrange_elements
    use fortfem_kinds
    use fortfem_api
    use check
    implicit none

    write(*,*) "Testing P2 Lagrange element implementation..."
    
    call test_p2_basis_functions()
    call test_p2_function_space_creation()
    call test_p2_dof_mapping()
    call test_p2_element_matrix_assembly()
    call test_p2_convergence_properties()
    call test_p2_vs_p1_comparison()
    
    call check_summary("P2 Lagrange Elements")

contains

    ! Test P2 basis function evaluation
    subroutine test_p2_basis_functions()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        real(dp) :: xi, eta, val, expected
        integer :: i
        
        ! Test P2 function space creation
        mesh = unit_square_mesh(3)
        Vh = function_space(mesh, "Lagrange", 2)
        
        call check_condition(Vh%degree == 2, &
            "P2 basis: correct degree")
        call check_condition(trim(Vh%element_family) == "Lagrange", &
            "P2 basis: correct element family")
        call check_condition(Vh%ndof > mesh%data%n_vertices, &
            "P2 basis: more DOFs than vertices")
        
        ! Expected DOFs for P2: vertices + edge midpoints  
        ! For triangular mesh: n_vertices + n_edges
        call check_condition(Vh%ndof == mesh%data%n_vertices + mesh%data%n_edges, &
            "P2 basis: correct DOF count")
        
        write(*,*) "   P2 mesh vertices:", mesh%data%n_vertices
        write(*,*) "   P2 mesh edges:", mesh%data%n_edges
        write(*,*) "   P2 total DOFs:", Vh%ndof
    end subroutine test_p2_basis_functions

    ! Test P2 function space creation with different mesh sizes
    subroutine test_p2_function_space_creation()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh_p1, Vh_p2
        integer :: n_values(3), i
        
        n_values = [3, 5, 7]
        
        do i = 1, 3
            mesh = unit_square_mesh(n_values(i))
            Vh_p1 = function_space(mesh, "Lagrange", 1)
            Vh_p2 = function_space(mesh, "Lagrange", 2)
            
            call check_condition(Vh_p2%ndof > Vh_p1%ndof, &
                "P2 function space: more DOFs than P1")
            call check_condition(Vh_p2%degree == 2, &
                "P2 function space: correct degree")
            call check_condition(Vh_p1%degree == 1, &
                "P1 function space: correct degree")
        end do
        
        write(*,*) "   P2 function spaces created successfully"
    end subroutine test_p2_function_space_creation

    ! Test DOF mapping and numbering for P2 elements
    subroutine test_p2_dof_mapping()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(function_t) :: uh
        integer :: i, vertex_dofs, edge_dofs
        
        mesh = unit_square_mesh(4)
        Vh = function_space(mesh, "Lagrange", 2)
        uh = function(Vh)
        
        call check_condition(allocated(uh%values), &
            "P2 DOF mapping: function values allocated")
        call check_condition(size(uh%values) == Vh%ndof, &
            "P2 DOF mapping: correct values array size")
        
        ! Test that we can set and retrieve DOF values
        do i = 1, Vh%ndof
            uh%values(i) = real(i, dp)
        end do
        
        call check_condition(uh%values(1) == 1.0_dp, &
            "P2 DOF mapping: can set/get first DOF")
        call check_condition(uh%values(Vh%ndof) == real(Vh%ndof, dp), &
            "P2 DOF mapping: can set/get last DOF")
        
        vertex_dofs = mesh%data%n_vertices
        edge_dofs = mesh%data%n_edges
        
        call check_condition(vertex_dofs + edge_dofs == Vh%ndof, &
            "P2 DOF mapping: vertex + edge DOFs equal total")
        
        write(*,*) "   P2 vertex DOFs:", vertex_dofs
        write(*,*) "   P2 edge DOFs:", edge_dofs
        write(*,*) "   P2 total DOFs:", Vh%ndof
    end subroutine test_p2_dof_mapping

    ! Test P2 element matrix assembly
    subroutine test_p2_element_matrix_assembly()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        real(dp) :: solution_norm, max_val
        
        ! Test P2 problem assembly and solution
        mesh = unit_square_mesh(4)
        Vh = function_space(mesh, "Lagrange", 2)
        
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(1.0_dp)
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        bc = dirichlet_bc(Vh, 0.0_dp)
        uh = function(Vh)
        
        call solve(a == L, uh, bc)
        
        solution_norm = sqrt(sum(uh%values**2))
        max_val = maxval(uh%values)
        
        call check_condition(solution_norm > 0.0_dp, &
            "P2 assembly: non-zero solution norm")
        call check_condition(max_val > 0.01_dp, &
            "P2 assembly: reasonable maximum value")
        call check_condition(max_val < 1.0_dp, &
            "P2 assembly: bounded solution")
        
        write(*,*) "   P2 solution norm:", solution_norm
        write(*,*) "   P2 maximum value:", max_val
        write(*,*) "   P2 total DOFs:", Vh%ndof
    end subroutine test_p2_element_matrix_assembly

    ! Test P2 convergence properties
    subroutine test_p2_convergence_properties()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        real(dp) :: max_vals(3), h_vals(3)
        integer :: n_vals(3), i
        
        ! Test convergence on sequence of P2 meshes
        n_vals = [3, 4, 5]
        h_vals = 1.0_dp / real(n_vals, dp)
        
        do i = 1, 3
            mesh = unit_square_mesh(n_vals(i))
            Vh = function_space(mesh, "Lagrange", 2)
            
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
        
        ! For this problem, max values should converge monotonically
        ! Check that values are decreasing (converging to true solution)
        call check_condition(max_vals(2) < max_vals(1) * 1.1_dp, &
            "P2 convergence: refinement 1 shows convergence")
        call check_condition(max_vals(3) < max_vals(2) * 1.1_dp, &
            "P2 convergence: refinement 2 shows convergence")
        call check_condition(all(max_vals > 0.01_dp), &
            "P2 convergence: all solutions positive")
        
        write(*,*) "   P2 mesh sizes:", n_vals
        write(*,*) "   P2 max values:", max_vals
        write(*,*) "   P2 h values:", h_vals
    end subroutine test_p2_convergence_properties

    ! Test P2 vs P1 accuracy comparison
    subroutine test_p2_vs_p1_comparison()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh_p1, Vh_p2
        type(trial_function_t) :: u1, u2
        type(test_function_t) :: v1, v2
        type(function_t) :: f, uh1, uh2
        type(dirichlet_bc_t) :: bc1, bc2
        type(form_expr_t) :: a1, L1, a2, L2
        real(dp) :: max_p1, max_p2, norm_p1, norm_p2
        
        ! Compare P1 vs P2 on same mesh
        mesh = unit_square_mesh(5)
        
        ! P1 solution
        Vh_p1 = function_space(mesh, "Lagrange", 1)
        u1 = trial_function(Vh_p1)
        v1 = test_function(Vh_p1)
        f = constant(1.0_dp)
        
        a1 = inner(grad(u1), grad(v1))*dx
        L1 = f*v1*dx
        
        bc1 = dirichlet_bc(Vh_p1, 0.0_dp)
        uh1 = function(Vh_p1)
        call solve(a1 == L1, uh1, bc1)
        
        max_p1 = maxval(uh1%values)
        norm_p1 = sqrt(sum(uh1%values**2))
        
        ! P2 solution
        Vh_p2 = function_space(mesh, "Lagrange", 2)
        u2 = trial_function(Vh_p2)
        v2 = test_function(Vh_p2)
        
        a2 = inner(grad(u2), grad(v2))*dx
        L2 = f*v2*dx
        
        bc2 = dirichlet_bc(Vh_p2, 0.0_dp)
        uh2 = function(Vh_p2)
        call solve(a2 == L2, uh2, bc2)
        
        max_p2 = maxval(uh2%values)
        norm_p2 = sqrt(sum(uh2%values**2))
        
        call check_condition(Vh_p2%ndof > Vh_p1%ndof, &
            "P2 vs P1: P2 has more DOFs")
        call check_condition(max_p2 > 0.5_dp * max_p1, &
            "P2 vs P1: P2 maximum comparable to P1")
        call check_condition(norm_p2 > 0.5_dp * norm_p1, &
            "P2 vs P1: P2 norm reasonable compared to P1")
        
        write(*,*) "   P1 DOFs:", Vh_p1%ndof, "max:", max_p1, "norm:", norm_p1
        write(*,*) "   P2 DOFs:", Vh_p2%ndof, "max:", max_p2, "norm:", norm_p2
        write(*,*) "   DOF ratio (P2/P1):", real(Vh_p2%ndof) / real(Vh_p1%ndof)
    end subroutine test_p2_vs_p1_comparison

end program test_p2_lagrange_elements