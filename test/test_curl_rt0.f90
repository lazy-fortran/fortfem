program test_curl_rt0
    use fortfem_kinds, only: dp
    use fortfem_basis_edge_2d
    implicit none

    call test_curl_rt0_constants()
    call test_curl_rt0_values()
    call test_curl_integration()
    
    print *, "All curl RT0 tests passed!"

contains

    subroutine test_curl_rt0_constants()
        real(dp) :: xi, eta, triangle_area
        real(dp) :: curls(3)
        integer :: i
        
        ! Test curl values at different points in reference triangle
        triangle_area = 0.5_dp  ! Reference triangle area
        
        ! Test at multiple points - curl should be constant
        xi = 0.1_dp
        eta = 0.1_dp
        call evaluate_edge_basis_curl_2d(xi, eta, triangle_area, curls)
        
        ! Check that curls are constant (RT0 property)
        if (abs(curls(1) - 0.0_dp) > 1e-12_dp) then
            print *, "Error: curl of basis 1 should be 0.0, got", curls(1)
            stop 1
        end if
        
        if (abs(curls(2) - 2.0_dp) > 1e-12_dp) then
            print *, "Error: curl of basis 2 should be 2.0, got", curls(2)
            stop 1
        end if
        
        if (abs(curls(3) - (-2.0_dp)) > 1e-12_dp) then
            print *, "Error: curl of basis 3 should be -2.0, got", curls(3)
            stop 1
        end if
        
        ! Test at different point - should be same values
        xi = 0.8_dp
        eta = 0.1_dp
        call evaluate_edge_basis_curl_2d(xi, eta, triangle_area, curls)
        
        if (abs(curls(1) - 0.0_dp) > 1e-12_dp .or. &
            abs(curls(2) - 2.0_dp) > 1e-12_dp .or. &
            abs(curls(3) - (-2.0_dp)) > 1e-12_dp) then
            print *, "Error: curl values not constant across reference triangle"
            stop 1
        end if
        
        print *, "RT0 curl constant test passed"
    end subroutine
    
    subroutine test_curl_rt0_values()
        real(dp) :: triangle_area
        real(dp) :: curls(3)
        
        ! Test with different triangle areas
        triangle_area = 1.0_dp
        call evaluate_edge_basis_curl_2d(0.33_dp, 0.33_dp, triangle_area, curls)
        
        ! For area = 1.0, curl should be: 0, 1/area, -1/area = 0, 1.0, -1.0
        if (abs(curls(1) - 0.0_dp) > 1e-12_dp .or. &
            abs(curls(2) - 1.0_dp) > 1e-12_dp .or. &
            abs(curls(3) - (-1.0_dp)) > 1e-12_dp) then
            print *, "Error: curl values incorrect for unit area triangle"
            print *, "Expected: 0.0, 1.0, -1.0"
            print *, "Got:", curls
            stop 1
        end if
        
        ! Test with different area
        triangle_area = 0.25_dp
        call evaluate_edge_basis_curl_2d(0.33_dp, 0.33_dp, triangle_area, curls)
        
        ! For area = 0.25, curl should be: 0, 1/area, -1/area = 0, 4.0, -4.0
        if (abs(curls(1) - 0.0_dp) > 1e-12_dp .or. &
            abs(curls(2) - 4.0_dp) > 1e-12_dp .or. &
            abs(curls(3) - (-4.0_dp)) > 1e-12_dp) then
            print *, "Error: curl values incorrect for area=0.25 triangle"
            print *, "Expected: 0.0, 4.0, -4.0"
            print *, "Got:", curls
            stop 1
        end if
        
        print *, "RT0 curl values test passed"
    end subroutine
    
    subroutine test_curl_integration()
        real(dp) :: triangle_area
        real(dp) :: curls(3)
        real(dp) :: integral_sum
        
        ! Test that curl integrates properly over reference triangle
        triangle_area = 0.5_dp  ! Reference triangle area
        call evaluate_edge_basis_curl_2d(0.33_dp, 0.33_dp, triangle_area, curls)
        
        ! Integral of curl over triangle = curl * area
        integral_sum = (curls(1) + curls(2) + curls(3)) * triangle_area
        
        ! For RT0, sum of curls should be: 0 + 2 + (-2) = 0
        ! So integral should be 0 * 0.5 = 0.0
        if (abs(integral_sum - 0.0_dp) > 1e-12_dp) then
            print *, "Error: curl integral incorrect"
            print *, "Expected integral: 0.0, got:", integral_sum
            stop 1
        end if
        
        print *, "RT0 curl integration test passed"
    end subroutine

end program test_curl_rt0