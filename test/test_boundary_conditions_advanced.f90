program test_boundary_conditions_advanced
    use fortfem_kinds
    use fortfem_api
    use check
    implicit none

    write(*,*) "Testing advanced boundary condition scenarios..."
    
    call test_mixed_bc_different_segments()
    call test_bc_on_refined_mesh()
    call test_bc_matrix_modification()
    call test_bc_corner_constraints()
    call test_bc_convergence_study()
    
    call check_summary("Advanced Boundary Conditions")

contains

    ! Test different BC values on different boundary segments (if supported)
    subroutine test_mixed_bc_different_segments()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        integer :: i, left_count, right_count, top_count, bottom_count
        real(dp) :: x, y
        
        ! Create mesh
        mesh = unit_square_mesh(6)
        Vh = function_space(mesh, "Lagrange", 1)
        
        left_count = 0
        right_count = 0
        top_count = 0
        bottom_count = 0
        
        ! Count vertices on each boundary segment
        do i = 1, Vh%mesh%data%n_vertices
            if (Vh%mesh%data%is_boundary_vertex(i)) then
                x = Vh%mesh%data%vertices(1, i)
                y = Vh%mesh%data%vertices(2, i)
                
                if (abs(x - 0.0_dp) < 1.0e-10_dp) left_count = left_count + 1
                if (abs(x - 1.0_dp) < 1.0e-10_dp) right_count = right_count + 1
                if (abs(y - 0.0_dp) < 1.0e-10_dp) bottom_count = bottom_count + 1
                if (abs(y - 1.0_dp) < 1.0e-10_dp) top_count = top_count + 1
            end if
        end do
        
        call check_condition(left_count > 0 .and. right_count > 0 .and. &
                           top_count > 0 .and. bottom_count > 0, &
            "Mixed BC segments: all boundary segments have vertices")
        
        ! Test that each segment has reasonable number of vertices
        call check_condition(left_count >= 5 .and. left_count <= 10, &
            "Mixed BC segments: left boundary has reasonable vertex count")
        call check_condition(right_count >= 5 .and. right_count <= 10, &
            "Mixed BC segments: right boundary has reasonable vertex count")
        
        write(*,*) "   Left boundary vertices:", left_count
        write(*,*) "   Right boundary vertices:", right_count
        write(*,*) "   Top boundary vertices:", top_count
        write(*,*) "   Bottom boundary vertices:", bottom_count
    end subroutine test_mixed_bc_different_segments

    ! Test boundary conditions on refined mesh
    subroutine test_bc_on_refined_mesh()
        type(mesh_t) :: mesh_coarse, mesh_fine
        type(function_space_t) :: Vh_coarse, Vh_fine
        type(trial_function_t) :: u_coarse, u_fine
        type(test_function_t) :: v_coarse, v_fine
        type(function_t) :: f_coarse, f_fine, uh_coarse, uh_fine
        type(dirichlet_bc_t) :: bc_coarse, bc_fine
        type(form_expr_t) :: a_coarse, a_fine, L_coarse, L_fine
        real(dp) :: max_coarse, max_fine
        
        ! Test on coarse mesh
        mesh_coarse = unit_square_mesh(4)
        Vh_coarse = function_space(mesh_coarse, "Lagrange", 1)
        
        u_coarse = trial_function(Vh_coarse)
        v_coarse = test_function(Vh_coarse)
        f_coarse = constant(1.0_dp)
        
        a_coarse = inner(grad(u_coarse), grad(v_coarse))*dx
        L_coarse = f_coarse*v_coarse*dx
        
        bc_coarse = dirichlet_bc(Vh_coarse, 0.0_dp)
        uh_coarse = function(Vh_coarse)
        
        call solve(a_coarse == L_coarse, uh_coarse, bc_coarse)
        max_coarse = maxval(uh_coarse%values)
        
        ! Test on fine mesh
        mesh_fine = unit_square_mesh(8)
        Vh_fine = function_space(mesh_fine, "Lagrange", 1)
        
        u_fine = trial_function(Vh_fine)
        v_fine = test_function(Vh_fine)
        f_fine = constant(1.0_dp)
        
        a_fine = inner(grad(u_fine), grad(v_fine))*dx
        L_fine = f_fine*v_fine*dx
        
        bc_fine = dirichlet_bc(Vh_fine, 0.0_dp)
        uh_fine = function(Vh_fine)
        
        call solve(a_fine == L_fine, uh_fine, bc_fine)
        max_fine = maxval(uh_fine%values)
        
        ! Fine mesh should give more accurate solution (higher maximum)
        call check_condition(max_fine > max_coarse, &
            "Refined mesh: finer mesh gives more accurate solution")
        call check_condition(abs(max_fine - max_coarse) / max_coarse < 0.5_dp, &
            "Refined mesh: solutions are reasonably close")
        
        write(*,*) "   Coarse mesh (4x4) max:", max_coarse
        write(*,*) "   Fine mesh (8x8) max:", max_fine
        write(*,*) "   Relative difference:", abs(max_fine - max_coarse) / max_coarse
        call plot(uh_fine, filename="build/bc_refined_mesh_solution.png", &
                  title="Refined mesh BC solution (8x8)")
    end subroutine test_bc_on_refined_mesh

    ! Test that BC modification properly zeroes matrix rows and sets diagonal
    subroutine test_bc_matrix_modification()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        integer :: i, j, boundary_node, ndof
        logical :: row_zeroed, diagonal_is_one, found_boundary_node
        type(dirichlet_bc_t) :: bc
        real(dp), allocatable :: K(:,:), F(:)
        real(dp) :: max_offdiag, tol, bc_value
        
        ! Create small mesh for testing
        mesh = unit_square_mesh(3)
        Vh = function_space(mesh, "Lagrange", 1)
        ndof = Vh%ndof
        
        ! Find a boundary node
        boundary_node = -1
        do i = 1, ndof
            if (Vh%mesh%data%is_boundary_vertex(i)) then
                boundary_node = i
                exit
            end if
        end do
        
        found_boundary_node = boundary_node > 0
        call check_condition(found_boundary_node, &
            "Matrix modification: found boundary node for testing")
        
        if (found_boundary_node) then
            write(*,*) "   Testing with boundary node:", boundary_node
            write(*,*) "   Total DOFs:", ndof
            
            ! Assemble Laplacian system and inspect modified matrix row
            bc_value = 2.0_dp
            bc = dirichlet_bc(Vh, bc_value)
            call assemble_laplacian_system(Vh, bc, K, F)
            
            tol = 1.0e-10_dp
            max_offdiag = 0.0_dp
            do j = 1, ndof
                if (j == boundary_node) cycle
                max_offdiag = max(max_offdiag, abs(K(boundary_node, j)))
            end do
            row_zeroed = max_offdiag < tol
            diagonal_is_one = abs(K(boundary_node, boundary_node) - 1.0_dp) < tol
            
            call check_condition(row_zeroed, &
                "Matrix modification: boundary row off-diagonals zeroed")
            call check_condition(diagonal_is_one, &
                "Matrix modification: boundary diagonal set to one")
            call check_condition(abs(F(boundary_node) - bc_value) < tol, &
                "Matrix modification: RHS matches BC value")
            
            deallocate(K, F)
        end if
        
    end subroutine test_bc_matrix_modification

    ! Test corner node constraints (nodes that belong to multiple boundaries)
    subroutine test_bc_corner_constraints()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        integer :: i, corner_node
        real(dp) :: x, y, corner_value
        logical :: found_corner
        
        ! Create mesh
        mesh = unit_square_mesh(5)
        Vh = function_space(mesh, "Lagrange", 1)
        
        ! Solve problem
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(1.0_dp)
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        bc = dirichlet_bc(Vh, 1.5_dp)
        uh = function(Vh)
        
        call solve(a == L, uh, bc)
        
        ! Find a corner node and check its value
        found_corner = .false.
        do i = 1, Vh%ndof
            if (Vh%mesh%data%is_boundary_vertex(i)) then
                x = Vh%mesh%data%vertices(1, i)
                y = Vh%mesh%data%vertices(2, i)
                
                ! Check for bottom-left corner
                if (abs(x - 0.0_dp) < 1.0e-10_dp .and. abs(y - 0.0_dp) < 1.0e-10_dp) then
                    corner_node = i
                    corner_value = uh%values(i)
                    found_corner = .true.
                    exit
                end if
            end if
        end do
        
        call check_condition(found_corner, &
            "Corner constraints: found corner node")
        
        if (found_corner) then
            call check_condition(abs(corner_value - 1.5_dp) < 1.0e-12_dp, &
                "Corner constraints: corner value matches BC")
            write(*,*) "   Corner node:", corner_node
            write(*,*) "   Corner value:", corner_value
            write(*,*) "   Expected:", 1.5_dp
        end if
    end subroutine test_bc_corner_constraints

    ! Test convergence properties with boundary conditions
    subroutine test_bc_convergence_study()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        real(dp) :: max_val_4, max_val_8, max_val_16
        real(dp) :: conv_rate_1, conv_rate_2
        
        ! Test on 4x4 mesh
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
        max_val_4 = maxval(uh%values)
        
        ! Test on 8x8 mesh
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
        max_val_8 = maxval(uh%values)
        
        ! Test on 16x16 mesh
        mesh = unit_square_mesh(16)
        Vh = function_space(mesh, "Lagrange", 1)
        
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(1.0_dp)
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        bc = dirichlet_bc(Vh, 0.0_dp)
        uh = function(Vh)
        
        call solve(a == L, uh, bc)
        max_val_16 = maxval(uh%values)
        
        ! Check convergence: values should increase and converge
        call check_condition(max_val_8 > max_val_4, &
            "Convergence: 8x8 mesh gives higher maximum than 4x4")
        call check_condition(max_val_16 > max_val_8, &
            "Convergence: 16x16 mesh gives higher maximum than 8x8")
        
        ! Check that convergence is slowing down (approaching limit)
        conv_rate_1 = (max_val_8 - max_val_4) / max_val_4
        conv_rate_2 = (max_val_16 - max_val_8) / max_val_8
        
        call check_condition(conv_rate_2 < conv_rate_1, &
            "Convergence: convergence rate is decreasing")
        
        write(*,*) "   Max value 4x4:", max_val_4
        write(*,*) "   Max value 8x8:", max_val_8
        write(*,*) "   Max value 16x16:", max_val_16
        write(*,*) "   Convergence rate 4->8:", conv_rate_1
        write(*,*) "   Convergence rate 8->16:", conv_rate_2
    end subroutine test_bc_convergence_study

end program test_boundary_conditions_advanced
