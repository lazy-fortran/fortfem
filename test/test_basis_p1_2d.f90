program test_basis_p1_2d
    ! Unit tests for P1 2D basis functions
    use fortfem_kinds
    use basis_p1_2d_module
    implicit none
    
    logical :: all_passed = .true.
    
    print *, "=== Testing P1 2D Basis Functions ==="
    print *, ""
    
    call test_basis_values()
    call test_basis_gradients()
    call test_partition_of_unity()
    call test_linear_reproduction()
    call test_jacobian_computation()
    call test_coordinate_transformation()
    
    print *, ""
    if (all_passed) then
        print *, "✅ ALL P1 BASIS TESTS PASSED!"
    else
        print *, "❌ SOME P1 BASIS TESTS FAILED!"
        stop 1
    end if
    
contains

    subroutine test_basis_values()
        type(basis_p1_2d_t) :: basis
        real(dp) :: val
        real(dp), parameter :: tol = 1e-14_dp
        logical :: passed = .true.
        
        print *, "Testing basis function values..."
        
        ! Test at vertices
        ! At (0,0): phi_1 = 1, phi_2 = 0, phi_3 = 0
        val = basis%eval(1, 0.0_dp, 0.0_dp)
        if (abs(val - 1.0_dp) > tol) then
            print *, "  ❌ phi_1(0,0) incorrect:", val
            passed = .false.
        end if
        
        val = basis%eval(2, 0.0_dp, 0.0_dp)
        if (abs(val) > tol) then
            print *, "  ❌ phi_2(0,0) incorrect:", val
            passed = .false.
        end if
        
        val = basis%eval(3, 0.0_dp, 0.0_dp)
        if (abs(val) > tol) then
            print *, "  ❌ phi_3(0,0) incorrect:", val
            passed = .false.
        end if
        
        ! At (1,0): phi_1 = 0, phi_2 = 1, phi_3 = 0
        val = basis%eval(1, 1.0_dp, 0.0_dp)
        if (abs(val) > tol) then
            print *, "  ❌ phi_1(1,0) incorrect:", val
            passed = .false.
        end if
        
        val = basis%eval(2, 1.0_dp, 0.0_dp)
        if (abs(val - 1.0_dp) > tol) then
            print *, "  ❌ phi_2(1,0) incorrect:", val
            passed = .false.
        end if
        
        ! At (0,1): phi_1 = 0, phi_2 = 0, phi_3 = 1
        val = basis%eval(3, 0.0_dp, 1.0_dp)
        if (abs(val - 1.0_dp) > tol) then
            print *, "  ❌ phi_3(0,1) incorrect:", val
            passed = .false.
        end if
        
        ! At centroid (1/3, 1/3): all basis functions = 1/3
        val = basis%eval(1, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp)
        if (abs(val - 1.0_dp/3.0_dp) > tol) then
            print *, "  ❌ phi_1 at centroid incorrect:", val
            passed = .false.
        end if
        
        if (passed) then
            print *, "  ✅ Basis function values test passed"
        else
            all_passed = .false.
        end if
    end subroutine test_basis_values
    
    subroutine test_basis_gradients()
        type(basis_p1_2d_t) :: basis
        real(dp) :: grad(2)
        real(dp), parameter :: tol = 1e-14_dp
        logical :: passed = .true.
        
        print *, "Testing basis function gradients..."
        
        ! Gradients are constant for P1 elements
        ! grad(phi_1) = [-1, -1]
        grad = basis%grad(1, 0.5_dp, 0.3_dp)  ! Any point
        if (abs(grad(1) + 1.0_dp) > tol .or. abs(grad(2) + 1.0_dp) > tol) then
            print *, "  ❌ grad(phi_1) incorrect:", grad
            passed = .false.
        end if
        
        ! grad(phi_2) = [1, 0]
        grad = basis%grad(2, 0.2_dp, 0.7_dp)
        if (abs(grad(1) - 1.0_dp) > tol .or. abs(grad(2)) > tol) then
            print *, "  ❌ grad(phi_2) incorrect:", grad
            passed = .false.
        end if
        
        ! grad(phi_3) = [0, 1]
        grad = basis%grad(3, 0.1_dp, 0.1_dp)
        if (abs(grad(1)) > tol .or. abs(grad(2) - 1.0_dp) > tol) then
            print *, "  ❌ grad(phi_3) incorrect:", grad
            passed = .false.
        end if
        
        if (passed) then
            print *, "  ✅ Basis function gradients test passed"
        else
            all_passed = .false.
        end if
    end subroutine test_basis_gradients
    
    subroutine test_partition_of_unity()
        type(basis_p1_2d_t) :: basis
        real(dp) :: sum_phi, xi, eta
        real(dp), parameter :: tol = 1e-14_dp
        logical :: passed = .true.
        integer :: i, j
        
        print *, "Testing partition of unity..."
        
        ! Test at various points
        do i = 0, 10
            do j = 0, 10-i
                xi = real(i, dp) / 10.0_dp
                eta = real(j, dp) / 10.0_dp
                
                ! Sum of basis functions should be 1
                sum_phi = basis%eval(1, xi, eta) + &
                         basis%eval(2, xi, eta) + &
                         basis%eval(3, xi, eta)
                
                if (abs(sum_phi - 1.0_dp) > tol) then
                    print *, "  ❌ Partition of unity failed at (", xi, ",", eta, ")"
                    print *, "    Sum =", sum_phi
                    passed = .false.
                    exit
                end if
            end do
            if (.not. passed) exit
        end do
        
        if (passed) then
            print *, "  ✅ Partition of unity test passed"
        else
            all_passed = .false.
        end if
    end subroutine test_partition_of_unity
    
    subroutine test_linear_reproduction()
        type(basis_p1_2d_t) :: basis
        real(dp) :: vertices(2,3), xi, eta, x, y
        real(dp) :: f_exact, f_approx
        real(dp) :: f_nodes(3)
        real(dp), parameter :: tol = 1e-14_dp
        logical :: passed = .true.
        integer :: i
        
        print *, "Testing linear function reproduction..."
        
        ! Test triangle
        vertices(:,1) = [0.0_dp, 0.0_dp]
        vertices(:,2) = [2.0_dp, 0.0_dp]
        vertices(:,3) = [0.0_dp, 3.0_dp]
        
        ! Test function: f(x,y) = 2x + 3y + 1
        ! Nodal values
        f_nodes(1) = 1.0_dp              ! f(0,0) = 1
        f_nodes(2) = 5.0_dp              ! f(2,0) = 5
        f_nodes(3) = 10.0_dp             ! f(0,3) = 10
        
        ! Test at various points
        do i = 1, 5
            xi = 0.2_dp * real(i-1, dp)
            eta = 0.2_dp * (1.0_dp - xi)
            
            ! Get physical coordinates
            call basis%transform_to_physical(xi, eta, vertices, x, y)
            
            ! Exact value
            f_exact = 2.0_dp*x + 3.0_dp*y + 1.0_dp
            
            ! Approximated value
            f_approx = f_nodes(1) * basis%eval(1, xi, eta) + &
                      f_nodes(2) * basis%eval(2, xi, eta) + &
                      f_nodes(3) * basis%eval(3, xi, eta)
            
            if (abs(f_exact - f_approx) > tol) then
                print *, "  ❌ Linear reproduction failed at (", x, ",", y, ")"
                print *, "    Exact:", f_exact, "Approx:", f_approx
                passed = .false.
            end if
        end do
        
        if (passed) then
            print *, "  ✅ Linear function reproduction test passed"
        else
            all_passed = .false.
        end if
    end subroutine test_linear_reproduction
    
    subroutine test_jacobian_computation()
        type(basis_p1_2d_t) :: basis
        real(dp) :: vertices(2,3), jac(2,2), det_j
        real(dp) :: expected_det
        real(dp), parameter :: tol = 1e-14_dp
        logical :: passed = .true.
        
        print *, "Testing Jacobian computation..."
        
        ! Test 1: Reference triangle
        vertices(:,1) = [0.0_dp, 0.0_dp]
        vertices(:,2) = [1.0_dp, 0.0_dp]
        vertices(:,3) = [0.0_dp, 1.0_dp]
        
        call basis%compute_jacobian(vertices, jac, det_j)
        
        ! Expected Jacobian: [[1, 0], [0, 1]]
        if (abs(jac(1,1) - 1.0_dp) > tol .or. abs(jac(1,2)) > tol .or. &
            abs(jac(2,1)) > tol .or. abs(jac(2,2) - 1.0_dp) > tol) then
            print *, "  ❌ Reference triangle Jacobian incorrect!"
            print *, "    J =", jac
            passed = .false.
        end if
        
        if (abs(det_j - 1.0_dp) > tol) then
            print *, "  ❌ Reference triangle determinant incorrect:", det_j
            passed = .false.
        end if
        
        ! Test 2: Scaled triangle
        vertices(:,1) = [0.0_dp, 0.0_dp]
        vertices(:,2) = [2.0_dp, 0.0_dp]
        vertices(:,3) = [0.0_dp, 3.0_dp]
        
        call basis%compute_jacobian(vertices, jac, det_j)
        expected_det = 6.0_dp  ! 2*3
        
        if (abs(det_j - expected_det) > tol) then
            print *, "  ❌ Scaled triangle determinant incorrect!"
            print *, "    Expected:", expected_det, "Got:", det_j
            passed = .false.
        end if
        
        if (passed) then
            print *, "  ✅ Jacobian computation test passed"
        else
            all_passed = .false.
        end if
    end subroutine test_jacobian_computation
    
    subroutine test_coordinate_transformation()
        type(basis_p1_2d_t) :: basis
        real(dp) :: vertices(2,3), x, y
        real(dp), parameter :: tol = 1e-14_dp
        logical :: passed = .true.
        
        print *, "Testing coordinate transformation..."
        
        ! Test triangle
        vertices(:,1) = [1.0_dp, 2.0_dp]
        vertices(:,2) = [4.0_dp, 3.0_dp]
        vertices(:,3) = [2.0_dp, 5.0_dp]
        
        ! Test vertex mapping
        call basis%transform_to_physical(0.0_dp, 0.0_dp, vertices, x, y)
        if (abs(x - 1.0_dp) > tol .or. abs(y - 2.0_dp) > tol) then
            print *, "  ❌ Vertex 1 mapping incorrect:", x, y
            passed = .false.
        end if
        
        call basis%transform_to_physical(1.0_dp, 0.0_dp, vertices, x, y)
        if (abs(x - 4.0_dp) > tol .or. abs(y - 3.0_dp) > tol) then
            print *, "  ❌ Vertex 2 mapping incorrect:", x, y
            passed = .false.
        end if
        
        call basis%transform_to_physical(0.0_dp, 1.0_dp, vertices, x, y)
        if (abs(x - 2.0_dp) > tol .or. abs(y - 5.0_dp) > tol) then
            print *, "  ❌ Vertex 3 mapping incorrect:", x, y
            passed = .false.
        end if
        
        ! Test centroid mapping
        call basis%transform_to_physical(1.0_dp/3.0_dp, 1.0_dp/3.0_dp, vertices, x, y)
        if (abs(x - 7.0_dp/3.0_dp) > tol .or. abs(y - 10.0_dp/3.0_dp) > tol) then
            print *, "  ❌ Centroid mapping incorrect:", x, y
            passed = .false.
        end if
        
        if (passed) then
            print *, "  ✅ Coordinate transformation test passed"
        else
            all_passed = .false.
        end if
    end subroutine test_coordinate_transformation

end program test_basis_p1_2d