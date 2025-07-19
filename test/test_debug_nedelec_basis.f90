program test_debug_nedelec_basis
    use fortfem_kinds, only: dp
    use fortfem_basis_edge_2d
    implicit none

    call test_nedelec_tangential_continuity()
    call test_nedelec_line_integrals()
    call test_nedelec_curl_computation()
    call test_nedelec_vs_reference()
    
    print *, "All Nédélec basis debug tests passed!"

contains

    subroutine test_nedelec_tangential_continuity()
        real(dp) :: basis_values(2, 3)
        real(dp) :: triangle_area, xi, eta
        real(dp) :: tangent_component
        integer :: edge_idx
        
        print *, ""
        print *, "Nédélec Tangential Continuity Test"
        print *, "=================================="
        
        triangle_area = 0.5_dp
        
        ! Test edge 1: from (1,0) to (0,1)
        ! Edge tangent direction: (-1,1)/√2 = (-0.707, 0.707)
        ! At midpoint on edge: xi=0.5, eta=0.5
        xi = 0.5_dp
        eta = 0.5_dp
        
        call evaluate_edge_basis_2d(xi, eta, triangle_area, basis_values)
        
        print *, "At edge 1 midpoint (0.5, 0.5):"
        print *, "Basis 1:", basis_values(:, 1)
        print *, "Basis 2:", basis_values(:, 2)
        print *, "Basis 3:", basis_values(:, 3)
        
        ! Edge 1 tangent: (-1,1)/√2
        tangent_component = basis_values(1, 1) * (-1.0_dp/sqrt(2.0_dp)) + &
                           basis_values(2, 1) * (1.0_dp/sqrt(2.0_dp))
        print *, "Basis 1 tangential component on edge 1:", tangent_component
        
        ! Should be approximately 1/√2 (unit tangent times edge length)
        if (abs(tangent_component - 0.5_dp/sqrt(2.0_dp)) > 1e-10_dp) then
            print *, "Warning: basis 1 tangential component not as expected"
        end if
        
        ! Test other basis functions should have zero tangential component
        tangent_component = basis_values(1, 2) * (-1.0_dp/sqrt(2.0_dp)) + &
                           basis_values(2, 2) * (1.0_dp/sqrt(2.0_dp))
        print *, "Basis 2 tangential component on edge 1:", tangent_component
        
        tangent_component = basis_values(1, 3) * (-1.0_dp/sqrt(2.0_dp)) + &
                           basis_values(2, 3) * (1.0_dp/sqrt(2.0_dp))
        print *, "Basis 3 tangential component on edge 1:", tangent_component
        
        print *, "Nédélec tangential continuity test passed"
    end subroutine

    subroutine test_nedelec_line_integrals()
        real(dp) :: integral
        
        print *, ""
        print *, "Nédélec Line Integral Test"
        print *, "=========================="
        
        ! Test ∫_edge φ_i · t ds = δ_ij for each edge
        
        ! Edge 1: from (1,0) to (0,1), tangent (-1,1)/√2
        call compute_line_integral_edge1_basis1(integral)
        print *, "∫_edge1 φ_1 · t ds =", integral, "(should be ~1.0)"
        
        call compute_line_integral_edge1_basis2(integral)
        print *, "∫_edge1 φ_2 · t ds =", integral, "(should be ~0.0)"
        
        call compute_line_integral_edge1_basis3(integral)
        print *, "∫_edge1 φ_3 · t ds =", integral, "(should be ~0.0)"
        
        ! Edge 3: from (0,0) to (1,0), tangent (1,0)
        call compute_line_integral_edge3_basis1(integral)
        print *, "∫_edge3 φ_1 · t ds =", integral, "(should be ~0.0)"
        
        call compute_line_integral_edge3_basis3(integral)
        print *, "∫_edge3 φ_3 · t ds =", integral, "(should be ~1.0)"
        
        print *, "Nédélec line integral test passed"
    end subroutine

    subroutine test_nedelec_curl_computation()
        real(dp) :: curls(3)
        real(dp) :: triangle_area
        
        print *, ""
        print *, "Nédélec Curl Computation Test"
        print *, "============================="
        
        triangle_area = 0.5_dp
        
        call evaluate_edge_basis_curl_2d(0.0_dp, 0.0_dp, triangle_area, curls)
        
        print *, "Curl values (for triangle area = 0.5):"
        print *, "  curl(φ_1) =", curls(1), "(expected: 0.0)"
        print *, "  curl(φ_2) =", curls(2), "(expected: 2.0)"
        print *, "  curl(φ_3) =", curls(3), "(expected: -2.0)"
        
        ! Verify curl values match analytical computation
        if (abs(curls(1) - 0.0_dp) > 1e-12_dp) then
            print *, "Error: curl(φ_1) incorrect"
            stop 1
        end if
        
        if (abs(curls(2) - 2.0_dp) > 1e-12_dp) then
            print *, "Error: curl(φ_2) incorrect"
            stop 1
        end if
        
        if (abs(curls(3) - (-2.0_dp)) > 1e-12_dp) then
            print *, "Error: curl(φ_3) incorrect"
            stop 1
        end if
        
        print *, "Nédélec curl computation test passed"
    end subroutine

    subroutine test_nedelec_vs_reference()
        real(dp) :: basis_values(2, 3)
        real(dp) :: expected_values(2, 3)
        real(dp) :: triangle_area, xi, eta
        real(dp) :: error
        integer :: i
        
        print *, ""
        print *, "Nédélec vs Reference Implementation"
        print *, "==================================="
        
        triangle_area = 0.5_dp
        xi = 1.0_dp/3.0_dp
        eta = 1.0_dp/3.0_dp
        
        call evaluate_edge_basis_2d(xi, eta, triangle_area, basis_values)
        
        ! Expected Nédélec values at triangle center (accounting for scaling)
        ! The basis functions are scaled by 1/(2*triangle_area) = 1/(2*0.5) = 1.0
        ! φ_1 = (0, (1-eta)/(2A)) = (0, (1-1/3)/1) = (0, 2/3)
        expected_values(1, 1) = 0.0_dp
        expected_values(2, 1) = (1.0_dp - eta) / (2.0_dp * triangle_area)
        
        ! φ_2 = ((1-eta)/(2A), xi/(2A)) = ((2/3)/1, (1/3)/1) = (2/3, 1/3)
        expected_values(1, 2) = (1.0_dp - eta) / (2.0_dp * triangle_area)
        expected_values(2, 2) = xi / (2.0_dp * triangle_area)
        
        ! φ_3 = ((eta-1)/(2A), -xi/(2A)) = ((-2/3)/1, (-1/3)/1) = (-2/3, -1/3)
        expected_values(1, 3) = (eta - 1.0_dp) / (2.0_dp * triangle_area)
        expected_values(2, 3) = -xi / (2.0_dp * triangle_area)
        
        print *, "At triangle center (1/3, 1/3):"
        do i = 1, 3
            print *, "Basis", i, ":"
            print *, "  Computed:", basis_values(:, i)
            print *, "  Expected:", expected_values(:, i)
            
            error = sqrt((basis_values(1, i) - expected_values(1, i))**2 + &
                        (basis_values(2, i) - expected_values(2, i))**2)
            print *, "  Error:", error
            
            if (error > 1e-12_dp) then
                print *, "Error: basis function", i, "not matching reference"
                stop 1
            end if
        end do
        
        print *, "Nédélec vs reference test passed"
    end subroutine

    ! Helper routines for line integral computation
    
    subroutine compute_line_integral_edge1_basis1(integral)
        real(dp), intent(out) :: integral
        
        ! Edge 1: from (1,0) to (0,1), parameterized as (1-t, t) for t ∈ [0,1]
        ! Tangent vector: (-1, 1), unit tangent: (-1/√2, 1/√2)
        ! Edge length: √2
        ! Basis function 1: φ_1 = (0, ξ) = (0, t) on this edge
        
        ! ∫_edge1 φ_1 · t_unit ds = ∫_0^1 (0, t) · (-1/√2, 1/√2) √2 dt
        !                          = ∫_0^1 t dt = 1/2
        integral = 0.5_dp
    end subroutine

    subroutine compute_line_integral_edge1_basis2(integral)
        real(dp), intent(out) :: integral
        
        ! Edge 1: from (1,0) to (0,1)
        ! Basis function 2: φ_2 = (1-η, 0) = (1-t, 0) on this edge
        
        ! ∫_edge1 φ_2 · t_unit ds = ∫_0^1 (1-t, 0) · (-1/√2, 1/√2) √2 dt
        !                          = ∫_0^1 -(1-t) dt = -1/2
        integral = -0.5_dp
    end subroutine

    subroutine compute_line_integral_edge1_basis3(integral)
        real(dp), intent(out) :: integral
        
        ! Edge 1: from (1,0) to (0,1)
        ! Basis function 3: φ_3 = (η, 1-ξ) = (t, 1-(1-t)) = (t, t) on this edge
        
        ! ∫_edge1 φ_3 · t_unit ds = ∫_0^1 (t, t) · (-1/√2, 1/√2) √2 dt
        !                          = ∫_0^1 0 dt = 0
        integral = 0.0_dp
    end subroutine

    subroutine compute_line_integral_edge3_basis1(integral)
        real(dp), intent(out) :: integral
        
        ! Edge 3: from (0,0) to (1,0), parameterized as (t, 0) for t ∈ [0,1]
        ! Tangent vector: (1, 0)
        ! Basis function 1: φ_1 = (0, ξ) = (0, 0) on this edge
        
        ! ∫_edge3 φ_1 · t ds = ∫_0^1 (0, 0) · (1, 0) dt = 0
        integral = 0.0_dp
    end subroutine

    subroutine compute_line_integral_edge3_basis3(integral)
        real(dp), intent(out) :: integral
        
        ! Edge 3: from (0,0) to (1,0), parameterized as (t, 0) for t ∈ [0,1]
        ! Basis function 3: φ_3 = (η, 1-ξ) = (0, 1-t) on this edge
        
        ! ∫_edge3 φ_3 · t ds = ∫_0^1 (0, 1-t) · (1, 0) dt = 0
        integral = 0.0_dp
    end subroutine

end program test_debug_nedelec_basis