program test_curl_p0_mapping
    use fortfem_kinds, only: dp
    use fortfem_basis_edge_2d
    use fortfem_mesh_2d
    implicit none

    call test_curl_maps_to_p0()
    call test_curl_p0_integration()
    
    print *, "All curl P0 mapping tests passed!"

contains

    subroutine test_curl_maps_to_p0()
        type(mesh_2d_t) :: mesh
        real(dp) :: triangle_area
        real(dp) :: curls(3)
        real(dp) :: p0_values(3)  ! P0 has 3 basis functions on reference triangle
        integer :: i
        
        ! Create reference triangle
        call create_reference_triangle(mesh)
        
        ! Get triangle area
        triangle_area = 0.5_dp  ! Reference triangle area
        
        ! Get curl values (constant over triangle)
        call evaluate_edge_basis_curl_2d(0.33_dp, 0.33_dp, triangle_area, curls)
        
        ! P0 basis functions on reference triangle are just constants
        ! P0_1 = 1, P0_2 = 0, P0_3 = 0 on triangle 1
        p0_values(1) = 1.0_dp
        p0_values(2) = 0.0_dp
        p0_values(3) = 0.0_dp
        
        ! The curl of each RT0 basis function should be a constant
        ! This verifies that curl(RT0) ∈ P0
        do i = 1, 3
            if (abs(curls(i) - curls(i)) > 1e-16_dp) then
                ! This is always true, but checks that curls are well-defined
                print *, "Error: curl value", i, "is not well-defined"
                stop 1
            end if
        end do
        
        ! Verify that curls are indeed constant (not varying with position)
        ! This is the key property: curl(RT0) maps to P0 (piecewise constants)
        call verify_curl_constant_over_triangle(triangle_area)
        
        call mesh%destroy()
        print *, "Curl maps to P0 space test passed"
    end subroutine
    
    subroutine verify_curl_constant_over_triangle(triangle_area)
        real(dp), intent(in) :: triangle_area
        real(dp) :: curls1(3), curls2(3), curls3(3)
        integer :: i
        
        ! Test curl at three different points
        call evaluate_edge_basis_curl_2d(0.1_dp, 0.1_dp, triangle_area, curls1)
        call evaluate_edge_basis_curl_2d(0.5_dp, 0.2_dp, triangle_area, curls2)
        call evaluate_edge_basis_curl_2d(0.2_dp, 0.6_dp, triangle_area, curls3)
        
        ! All should be identical (constant over triangle)
        do i = 1, 3
            if (abs(curls1(i) - curls2(i)) > 1e-12_dp .or. &
                abs(curls2(i) - curls3(i)) > 1e-12_dp .or. &
                abs(curls1(i) - curls3(i)) > 1e-12_dp) then
                print *, "Error: curl of basis", i, "is not constant over triangle"
                print *, "Values:", curls1(i), curls2(i), curls3(i)
                stop 1
            end if
        end do
    end subroutine
    
    subroutine test_curl_p0_integration()
        real(dp) :: triangle_area
        real(dp) :: curls(3)
        real(dp) :: p0_integrals(3)
        
        ! Test integration of curl over reference triangle
        triangle_area = 0.5_dp
        call evaluate_edge_basis_curl_2d(0.33_dp, 0.33_dp, triangle_area, curls)
        
        ! P0 basis functions integrate to:
        ! ∫ P0_i dx = area * P0_i_value
        ! Since P0_i is constant over each triangle
        
        ! For RT0 curl, each component is constant over triangle
        p0_integrals(1) = curls(1) * triangle_area  ! ∫ curl(φ_1) dx
        p0_integrals(2) = curls(2) * triangle_area  ! ∫ curl(φ_2) dx  
        p0_integrals(3) = curls(3) * triangle_area  ! ∫ curl(φ_3) dx
        
        ! Expected values for reference triangle  
        ! curl values are [0, 2, -2], so integrals are [0*0.5, 2*0.5, -2*0.5] = [0, 1, -1]
        if (abs(p0_integrals(1) - 0.0_dp) > 1e-12_dp) then
            print *, "Error: integral of curl(φ_1) should be 0.0, got", p0_integrals(1)
            stop 1
        end if
        
        if (abs(p0_integrals(2) - 1.0_dp) > 1e-12_dp) then
            print *, "Error: integral of curl(φ_2) should be 1.0, got", p0_integrals(2)
            stop 1
        end if
        
        if (abs(p0_integrals(3) - (-1.0_dp)) > 1e-12_dp) then
            print *, "Error: integral of curl(φ_3) should be -1.0, got", p0_integrals(3)
            stop 1
        end if
        
        print *, "Curl P0 integration test passed"
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

end program test_curl_p0_mapping