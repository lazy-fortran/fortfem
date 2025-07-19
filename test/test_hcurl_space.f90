program test_hcurl_space
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    implicit none

    call test_hcurl_space_creation()
    call test_hcurl_space_evaluation()
    call test_hcurl_space_curl_evaluation()
    call test_hcurl_space_boundary_conditions()
    
    print *, "All H(curl) space tests passed!"

contains

    subroutine test_hcurl_space_creation()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        
        ! Create mesh with edge connectivity
        call create_reference_triangle(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        
        ! Create H(curl) space
        call space%init(mesh)
        
        ! Test DOF counts
        if (space%get_n_dofs() /= 3) then
            print *, "Error: expected 3 DOFs, got", space%get_n_dofs()
            stop 1
        end if
        
        if (space%get_n_interior_dofs() /= 0) then
            print *, "Error: expected 0 interior DOFs, got", space%get_n_interior_dofs()
            stop 1
        end if
        
        if (space%get_n_boundary_dofs() /= 3) then
            print *, "Error: expected 3 boundary DOFs, got", space%get_n_boundary_dofs()
            stop 1
        end if
        
        call space%destroy()
        call mesh%destroy()
        print *, "H(curl) space creation test passed"
    end subroutine
    
    subroutine test_hcurl_space_evaluation()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        real(dp) :: coeff(3), values(2)
        real(dp) :: xi, eta
        
        ! Create mesh and space
        call create_reference_triangle(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        ! Test evaluation with unit coefficients
        coeff = [1.0_dp, 0.0_dp, 0.0_dp]  ! Only first basis function
        xi = 0.33_dp
        eta = 0.33_dp
        
        call space%evaluate_at_point(1, xi, eta, coeff, values)
        
        ! Values should be non-zero for first basis function
        if (abs(values(1)) < 1e-16_dp .and. abs(values(2)) < 1e-16_dp) then
            print *, "Error: evaluation gives zero values for non-zero coefficient"
            stop 1
        end if
        
        ! Test with zero coefficients
        coeff = [0.0_dp, 0.0_dp, 0.0_dp]
        call space%evaluate_at_point(1, xi, eta, coeff, values)
        
        if (abs(values(1)) > 1e-12_dp .or. abs(values(2)) > 1e-12_dp) then
            print *, "Error: evaluation should give zero for zero coefficients"
            stop 1
        end if
        
        call space%destroy()
        call mesh%destroy()
        print *, "H(curl) space evaluation test passed"
    end subroutine
    
    subroutine test_hcurl_space_curl_evaluation()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        real(dp) :: coeff(3), curl_value
        real(dp) :: xi, eta
        
        ! Create mesh and space
        call create_reference_triangle(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        ! Test curl evaluation with unit coefficients
        coeff = [1.0_dp, 0.0_dp, 0.0_dp]  ! Only first basis function
        xi = 0.33_dp
        eta = 0.33_dp
        
        call space%evaluate_curl_at_point(1, xi, eta, coeff, curl_value)
        
        ! For RT0 on reference triangle, curl of first basis function = 0.0
        if (abs(curl_value - 0.0_dp) > 1e-12_dp) then
            print *, "Error: curl evaluation incorrect"
            print *, "Expected: 0.0, got:", curl_value
            stop 1
        end if
        
        ! Test with different coefficients
        coeff = [0.0_dp, 1.0_dp, 0.0_dp]  ! Only second basis function
        call space%evaluate_curl_at_point(1, xi, eta, coeff, curl_value)
        
        if (abs(curl_value - 2.0_dp) > 1e-12_dp) then
            print *, "Error: curl evaluation incorrect for second basis"
            print *, "Expected: 2.0, got:", curl_value
            stop 1
        end if
        
        call space%destroy()
        call mesh%destroy()
        print *, "H(curl) space curl evaluation test passed"
    end subroutine
    
    subroutine test_hcurl_space_boundary_conditions()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        integer :: i, n_essential
        
        ! Create mesh and space
        call create_reference_triangle(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        ! Apply homogeneous essential BC
        call space%apply_essential_bc(0.0_dp)
        
        ! Count essential DOFs
        n_essential = 0
        do i = 1, space%get_n_dofs()
            if (space%dof_is_essential(i) == 1) then
                n_essential = n_essential + 1
            end if
        end do
        
        ! For single triangle, all edges are boundary edges
        if (n_essential /= 3) then
            print *, "Error: expected 3 essential DOFs, got", n_essential
            stop 1
        end if
        
        ! Check that essential values are set correctly
        do i = 1, space%get_n_dofs()
            if (space%dof_is_essential(i) == 1) then
                if (abs(space%essential_values(i)) > 1e-12_dp) then
                    print *, "Error: essential value not set to zero"
                    stop 1
                end if
            end if
        end do
        
        call space%destroy()
        call mesh%destroy()
        print *, "H(curl) space boundary conditions test passed"
    end subroutine
    
    subroutine create_reference_triangle(mesh)
        type(mesh_2d_t), intent(out) :: mesh
        
        mesh%n_vertices = 3
        mesh%n_triangles = 1
        
        allocate(mesh%vertices(2, 3))
        allocate(mesh%triangles(3, 1))
        
        ! Reference triangle vertices: (0,0), (1,0), (0,1)
        mesh%vertices(1, 1) = 0.0_dp
        mesh%vertices(2, 1) = 0.0_dp
        mesh%vertices(1, 2) = 1.0_dp
        mesh%vertices(2, 2) = 0.0_dp
        mesh%vertices(1, 3) = 0.0_dp
        mesh%vertices(2, 3) = 1.0_dp
        
        ! Triangle connectivity
        mesh%triangles(1, 1) = 1
        mesh%triangles(2, 1) = 2
        mesh%triangles(3, 1) = 3
    end subroutine

end program test_hcurl_space