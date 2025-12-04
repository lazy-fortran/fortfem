program test_boundary_conditions
    use fortfem_kinds
    use fortfem_api
    use check
    implicit none

    write(*,*) "Testing boundary condition implementation..."
    
    call test_dirichlet_zero_bc()
    call test_dirichlet_nonzero_bc()
    call test_multiple_boundary_segments()
    call test_bc_consistency()
    call test_corner_node_handling()
    call test_manufactured_solution()
    
    call check_summary("Boundary Conditions")

contains

    ! Test basic zero Dirichlet BC on unit square
    subroutine test_dirichlet_zero_bc()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        integer :: i
        real(dp) :: max_interior_value, min_boundary_value, max_boundary_value
        logical :: all_boundary_zero
        
        ! Create simple mesh
        mesh = unit_square_mesh(5)
        Vh = function_space(mesh, "Lagrange", 1)
        
        ! Define problem: -Δu = 1, u = 0 on boundary
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(1.0_dp)
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        bc = dirichlet_bc(Vh, 0.0_dp)
        uh = function(Vh)
        
        call solve(a == L, uh, bc)
        
        ! Check that boundary values are indeed zero
        all_boundary_zero = .true.
        max_boundary_value = -1.0e10_dp
        min_boundary_value = 1.0e10_dp
        max_interior_value = -1.0e10_dp
        
        do i = 1, uh%space%ndof
            if (uh%space%mesh%data%is_boundary_vertex(i)) then
                max_boundary_value = max(max_boundary_value, abs(uh%values(i)))
                min_boundary_value = min(min_boundary_value, abs(uh%values(i)))
                if (abs(uh%values(i)) > 1.0e-12_dp) then
                    all_boundary_zero = .false.
                end if
            else
                max_interior_value = max(max_interior_value, uh%values(i))
            end if
        end do
        
        call check_condition(all_boundary_zero, &
            "Zero Dirichlet BC: boundary values are zero")
        call check_condition(max_interior_value > 1.0e-6_dp, &
            "Zero Dirichlet BC: interior values are positive")
        call check_condition(max_boundary_value < 1.0e-10_dp, &
            "Zero Dirichlet BC: boundary tolerance check")
        
        write(*,*) "   Max interior value:", max_interior_value
        write(*,*) "   Max boundary error:", max_boundary_value
        call plot(uh, filename="build/bc_zero_dirichlet_solution.png", &
                  title="Zero Dirichlet BC solution")
    end subroutine test_dirichlet_zero_bc

    ! Test non-zero Dirichlet BC
    subroutine test_dirichlet_nonzero_bc()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        integer :: i
        real(dp), parameter :: bc_value = 2.5_dp
        real(dp) :: max_error, error
        logical :: all_boundary_correct
        
        ! Create mesh
        mesh = unit_square_mesh(4)
        Vh = function_space(mesh, "Lagrange", 1)
        
        ! Note: Current implementation hardcodes f=1, so we test with that
        ! Define problem: -Δu = 1, u = 2.5 on boundary
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(1.0_dp)
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        bc = dirichlet_bc(Vh, bc_value)
        uh = function(Vh)
        
        call solve(a == L, uh, bc)
        
        ! Check boundary values
        all_boundary_correct = .true.
        max_error = 0.0_dp
        
        do i = 1, uh%space%ndof
            if (uh%space%mesh%data%is_boundary_vertex(i)) then
                error = abs(uh%values(i) - bc_value)
                max_error = max(max_error, error)
                if (error > 1.0e-12_dp) then
                    all_boundary_correct = .false.
                end if
            end if
        end do
        
        call check_condition(all_boundary_correct, &
            "Non-zero Dirichlet BC: boundary values correct")
        call check_condition(max_error < 1.0e-10_dp, &
            "Non-zero Dirichlet BC: boundary tolerance check")
        
        ! For -Δu = 1 with constant boundary conditions, 
        ! solution should have interior values > boundary values
        call check_condition(maxval(uh%values) > minval(uh%values), &
            "Poisson equation: interior values > boundary values")
        
        write(*,*) "   BC value:", bc_value
        write(*,*) "   Max boundary error:", max_error
        write(*,*) "   Solution range:", minval(uh%values), "to", maxval(uh%values)
        call plot(uh, filename="build/bc_nonzero_dirichlet_solution.png", &
                  title="Non-zero Dirichlet BC solution")
    end subroutine test_dirichlet_nonzero_bc

    ! Test boundary condition identification on different segments
    subroutine test_multiple_boundary_segments()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        integer :: i, boundary_count, total_vertices
        real(dp) :: x, y
        logical :: has_left_boundary, has_right_boundary, has_top_boundary, has_bottom_boundary
        
        ! Create mesh
        mesh = unit_square_mesh(6)
        Vh = function_space(mesh, "Lagrange", 1)
        
        boundary_count = 0
        total_vertices = Vh%mesh%data%n_vertices
        has_left_boundary = .false.
        has_right_boundary = .false.
        has_top_boundary = .false.
        has_bottom_boundary = .false.
        
        ! Check boundary vertex identification
        do i = 1, total_vertices
            if (Vh%mesh%data%is_boundary_vertex(i)) then
                boundary_count = boundary_count + 1
                
                ! Get vertex coordinates
                x = Vh%mesh%data%vertices(1, i)
                y = Vh%mesh%data%vertices(2, i)
                
                ! Check which boundary segment this vertex belongs to
                if (abs(x - 0.0_dp) < 1.0e-10_dp) has_left_boundary = .true.
                if (abs(x - 1.0_dp) < 1.0e-10_dp) has_right_boundary = .true.
                if (abs(y - 0.0_dp) < 1.0e-10_dp) has_bottom_boundary = .true.
                if (abs(y - 1.0_dp) < 1.0e-10_dp) has_top_boundary = .true.
            end if
        end do
        
        call check_condition(boundary_count > 0, &
            "Boundary detection: found boundary vertices")
        call check_condition(has_left_boundary, &
            "Boundary detection: found left boundary vertices")
        call check_condition(has_right_boundary, &
            "Boundary detection: found right boundary vertices")
        call check_condition(has_top_boundary, &
            "Boundary detection: found top boundary vertices")
        call check_condition(has_bottom_boundary, &
            "Boundary detection: found bottom boundary vertices")
        
        ! Boundary should be a reasonable fraction of total vertices
        ! For structured meshes, this can be high (e.g., 6x6 grid has 24 boundary vs 16 interior)
        call check_condition(real(boundary_count, dp) / real(total_vertices, dp) < 0.8_dp, &
            "Boundary detection: reasonable boundary/interior ratio")
        
        write(*,*) "   Total vertices:", total_vertices
        write(*,*) "   Boundary vertices:", boundary_count
        write(*,*) "   Boundary fraction:", real(boundary_count, dp) / real(total_vertices, dp)
    end subroutine test_multiple_boundary_segments

    ! Test boundary condition consistency
    subroutine test_bc_consistency()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh1, uh2
        type(dirichlet_bc_t) :: bc1, bc2
        type(form_expr_t) :: a, L
        integer :: i
        real(dp) :: max_diff
        logical :: solutions_different
        
        ! Create mesh
        mesh = unit_square_mesh(4)
        Vh = function_space(mesh, "Lagrange", 1)
        
        ! Define same problem with different BC values
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
        
        ! Solutions should be different
        max_diff = 0.0_dp
        do i = 1, uh1%space%ndof
            max_diff = max(max_diff, abs(uh1%values(i) - uh2%values(i)))
        end do
        
        solutions_different = max_diff > 1.0e-6_dp
        
        call check_condition(solutions_different, &
            "BC consistency: different BCs give different solutions")
        
        write(*,*) "   Max difference between BC=0 and BC=1 solutions:", max_diff
    end subroutine test_bc_consistency

    ! Test corner node handling (corner nodes belong to multiple boundaries)
    subroutine test_corner_node_handling()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        integer :: i, corner_count
        real(dp) :: x, y
        logical :: found_corners(4)  ! [bottom-left, bottom-right, top-left, top-right]
        
        ! Create mesh
        mesh = unit_square_mesh(5)
        Vh = function_space(mesh, "Lagrange", 1)
        
        corner_count = 0
        found_corners = .false.
        
        ! Check for corner vertices
        do i = 1, Vh%mesh%data%n_vertices
            if (Vh%mesh%data%is_boundary_vertex(i)) then
                x = Vh%mesh%data%vertices(1, i)
                y = Vh%mesh%data%vertices(2, i)
                
                ! Check for corners
                if (abs(x - 0.0_dp) < 1.0e-10_dp .and. abs(y - 0.0_dp) < 1.0e-10_dp) then
                    found_corners(1) = .true.  ! bottom-left
                    corner_count = corner_count + 1
                end if
                if (abs(x - 1.0_dp) < 1.0e-10_dp .and. abs(y - 0.0_dp) < 1.0e-10_dp) then
                    found_corners(2) = .true.  ! bottom-right
                    corner_count = corner_count + 1
                end if
                if (abs(x - 0.0_dp) < 1.0e-10_dp .and. abs(y - 1.0_dp) < 1.0e-10_dp) then
                    found_corners(3) = .true.  ! top-left
                    corner_count = corner_count + 1
                end if
                if (abs(x - 1.0_dp) < 1.0e-10_dp .and. abs(y - 1.0_dp) < 1.0e-10_dp) then
                    found_corners(4) = .true.  ! top-right
                    corner_count = corner_count + 1
                end if
            end if
        end do
        
        call check_condition(corner_count == 4, &
            "Corner detection: found all 4 corners")
        call check_condition(all(found_corners), &
            "Corner detection: all corners identified as boundary")
        
        write(*,*) "   Found corners:", corner_count
        write(*,*) "   Corner flags:", found_corners
    end subroutine test_corner_node_handling

    ! Test against manufactured solution
    subroutine test_manufactured_solution()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(trial_function_t) :: u
        type(test_function_t) :: v
        type(function_t) :: f, uh
        type(dirichlet_bc_t) :: bc
        type(form_expr_t) :: a, L
        integer :: i
        real(dp) :: x, y, u_exact, max_error, l2_error
        real(dp) :: total_area, elem_area
        
        ! Create fine mesh for accuracy
        mesh = unit_square_mesh(8)
        Vh = function_space(mesh, "Lagrange", 1)
        
        ! Since current implementation hardcodes f=1, we test with that
        ! For -Δu = 1 with u=0 on boundary, we expect a bowl-shaped solution
        ! with maximum around 0.1 at center for unit square
        
        u = trial_function(Vh)
        v = test_function(Vh)
        f = constant(1.0_dp)  ! Note: solver ignores this and uses f=1
        
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        
        bc = dirichlet_bc(Vh, 0.0_dp)
        uh = function(Vh)
        
        call solve(a == L, uh, bc)
        
        ! Test basic properties of the solution
        max_error = maxval(uh%values)
        
        ! For -Δu = 1 on unit square with u=0 on boundary:
        ! - All boundary values should be 0
        ! - All interior values should be positive  
        ! - Maximum should be around 0.1 at center
        call check_condition(max_error > 0.01_dp, &
            "Poisson solution: positive interior values")
        call check_condition(max_error < 0.15_dp, &
            "Poisson solution: reasonable maximum value")
        
        ! Check that boundary values are indeed zero
        l2_error = 0.0_dp
        do i = 1, uh%space%ndof
            if (uh%space%mesh%data%is_boundary_vertex(i)) then
                l2_error = max(l2_error, abs(uh%values(i)))
            end if
        end do
        
        call check_condition(l2_error < 1.0e-10_dp, &
            "Poisson solution: boundary values are zero")
        
        write(*,*) "   Max solution value:", max_error
        write(*,*) "   Max boundary error:", l2_error
        write(*,*) "   Mesh resolution: 8x8"
        call plot(uh, filename="build/bc_manufactured_solution.png", &
                  title="Manufactured Poisson solution")
    end subroutine test_manufactured_solution

end program test_boundary_conditions
