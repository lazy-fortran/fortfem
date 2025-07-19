program test_basis_p2_2d
    ! Unit tests for P2 2D basis functions
    use fortfem_kinds
    use basis_p2_2d_module
    implicit none
    
    logical :: all_passed = .true.
    
    print *, "=== Testing P2 2D Basis Functions ==="
    print *, ""
    
    call test_basis_values_at_nodes()
    call test_basis_gradients()
    call test_partition_of_unity()
    call test_quadratic_reproduction()
    call test_jacobian_computation()
    call test_coordinate_transformation()
    call test_hessian_symmetry()
    call test_interpolation_property()
    
    print *, ""
    if (all_passed) then
        print *, "✅ ALL P2 BASIS TESTS PASSED!"
    else
        print *, "❌ SOME P2 BASIS TESTS FAILED!"
        stop 1
    end if
    
contains

    subroutine test_basis_values_at_nodes()
        type(basis_p2_2d_t) :: basis
        real(dp) :: val
        real(dp), parameter :: tol = 1e-14_dp
        logical :: passed = .true.
        integer :: i, j
        real(dp) :: xi, eta
        
        print *, "Testing basis function values at nodes..."
        
        ! Test Lagrange property: phi_i(node_j) = delta_ij
        do i = 1, 6
            do j = 1, 6
                xi = basis%nodes(1, j)
                eta = basis%nodes(2, j)
                val = basis%eval(i, xi, eta)
                
                if (i == j) then
                    ! Should be 1
                    if (abs(val - 1.0_dp) > tol) then
                        print *, "  ❌ phi_", i, "(node_", j, ") =", val, " (expected 1)"
                        passed = .false.
                    end if
                else
                    ! Should be 0
                    if (abs(val) > tol) then
                        print *, "  ❌ phi_", i, "(node_", j, ") =", val, " (expected 0)"
                        passed = .false.
                    end if
                end if
            end do
        end do
        
        if (passed) then
            print *, "  ✅ Basis function values at nodes test passed"
        else
            all_passed = .false.
        end if
    end subroutine test_basis_values_at_nodes
    
    subroutine test_basis_gradients()
        type(basis_p2_2d_t) :: basis
        real(dp) :: grad(2), grad_fd(2), val_plus, val_minus
        real(dp), parameter :: tol = 1e-8_dp
        real(dp), parameter :: h = 1e-6_dp
        logical :: passed = .true.
        integer :: i
        real(dp) :: xi, eta
        
        print *, "Testing basis function gradients..."
        
        ! Test at an interior point
        xi = 0.3_dp
        eta = 0.2_dp
        
        do i = 1, 6
            ! Analytical gradient
            grad = basis%grad(i, xi, eta)
            
            ! Finite difference gradient
            val_plus = basis%eval(i, xi + h, eta)
            val_minus = basis%eval(i, xi - h, eta)
            grad_fd(1) = (val_plus - val_minus) / (2.0_dp * h)
            
            val_plus = basis%eval(i, xi, eta + h)
            val_minus = basis%eval(i, xi, eta - h)
            grad_fd(2) = (val_plus - val_minus) / (2.0_dp * h)
            
            ! Compare
            if (abs(grad(1) - grad_fd(1)) > tol .or. abs(grad(2) - grad_fd(2)) > tol) then
                print *, "  ❌ grad(phi_", i, ") mismatch at (", xi, ",", eta, ")"
                print *, "    Analytical:", grad
                print *, "    Finite diff:", grad_fd
                passed = .false.
            end if
        end do
        
        if (passed) then
            print *, "  ✅ Basis function gradients test passed"
        else
            all_passed = .false.
        end if
    end subroutine test_basis_gradients
    
    subroutine test_partition_of_unity()
        type(basis_p2_2d_t) :: basis
        real(dp) :: sum_phi, xi, eta
        real(dp), parameter :: tol = 1e-14_dp
        logical :: passed = .true.
        integer :: i, j, k
        
        print *, "Testing partition of unity..."
        
        ! Test at various points in the reference triangle
        do i = 0, 10
            do j = 0, 10-i
                xi = real(i, dp) / 10.0_dp
                eta = real(j, dp) / 10.0_dp
                
                ! Sum of all basis functions
                sum_phi = 0.0_dp
                do k = 1, 6
                    sum_phi = sum_phi + basis%eval(k, xi, eta)
                end do
                
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
    
    subroutine test_quadratic_reproduction()
        type(basis_p2_2d_t) :: basis
        real(dp) :: vertices(2,3), xi, eta, x, y
        real(dp) :: f_exact, f_approx
        real(dp) :: f_nodes(6)
        real(dp), parameter :: tol = 1e-14_dp
        logical :: passed = .true.
        integer :: i, k
        real(dp) :: nodes_physical(2,6)
        
        print *, "Testing quadratic function reproduction..."
        
        ! Test triangle
        vertices(:,1) = [0.0_dp, 0.0_dp]
        vertices(:,2) = [2.0_dp, 0.0_dp]
        vertices(:,3) = [0.0_dp, 3.0_dp]
        
        ! Physical nodes for P2 element
        nodes_physical(:,1:3) = vertices
        nodes_physical(:,4) = 0.5_dp * (vertices(:,1) + vertices(:,2))
        nodes_physical(:,5) = 0.5_dp * (vertices(:,2) + vertices(:,3))
        nodes_physical(:,6) = 0.5_dp * (vertices(:,3) + vertices(:,1))
        
        ! Test function: f(x,y) = x^2 + 2*x*y + y^2 - x + 3*y + 1
        ! Compute nodal values
        do i = 1, 6
            x = nodes_physical(1,i)
            y = nodes_physical(2,i)
            f_nodes(i) = x**2 + 2.0_dp*x*y + y**2 - x + 3.0_dp*y + 1.0_dp
        end do
        
        ! Test at various points
        do k = 1, 10
            xi = 0.1_dp * real(k-1, dp)
            eta = 0.1_dp * (1.0_dp - xi)
            
            ! Get physical coordinates
            call basis%transform_to_physical(xi, eta, vertices, x, y)
            
            ! Exact value
            f_exact = x**2 + 2.0_dp*x*y + y**2 - x + 3.0_dp*y + 1.0_dp
            
            ! Approximated value
            f_approx = 0.0_dp
            do i = 1, 6
                f_approx = f_approx + f_nodes(i) * basis%eval(i, xi, eta)
            end do
            
            if (abs(f_exact - f_approx) > tol) then
                print *, "  ❌ Quadratic reproduction failed at (", x, ",", y, ")"
                print *, "    Exact:", f_exact, "Approx:", f_approx
                passed = .false.
            end if
        end do
        
        if (passed) then
            print *, "  ✅ Quadratic function reproduction test passed"
        else
            all_passed = .false.
        end if
    end subroutine test_quadratic_reproduction
    
    subroutine test_jacobian_computation()
        type(basis_p2_2d_t) :: basis
        real(dp) :: vertices(2,3), jac(2,2), det_j
        real(dp) :: expected_det
        real(dp), parameter :: tol = 1e-14_dp
        logical :: passed = .true.
        
        print *, "Testing Jacobian computation..."
        
        ! Test 1: Reference triangle (should be identity for P2)
        vertices(:,1) = [0.0_dp, 0.0_dp]
        vertices(:,2) = [1.0_dp, 0.0_dp]
        vertices(:,3) = [0.0_dp, 1.0_dp]
        
        call basis%compute_jacobian(vertices, jac, det_j)
        
        ! For reference triangle, Jacobian should be [[1,0],[0,1]] at centroid
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
        type(basis_p2_2d_t) :: basis
        real(dp) :: vertices(2,3), x, y
        real(dp), parameter :: tol = 1e-14_dp
        logical :: passed = .true.
        integer :: i
        real(dp) :: xi, eta
        
        print *, "Testing coordinate transformation..."
        
        ! Test triangle
        vertices(:,1) = [1.0_dp, 2.0_dp]
        vertices(:,2) = [4.0_dp, 3.0_dp]
        vertices(:,3) = [2.0_dp, 5.0_dp]
        
        ! Test node mapping (all 6 nodes)
        do i = 1, 6
            xi = basis%nodes(1, i)
            eta = basis%nodes(2, i)
            
            call basis%transform_to_physical(xi, eta, vertices, x, y)
            
            ! Check vertices
            if (i <= 3) then
                if (abs(x - vertices(1,i)) > tol .or. abs(y - vertices(2,i)) > tol) then
                    print *, "  ❌ Vertex", i, "mapping incorrect:", x, y
                    passed = .false.
                end if
            else
                ! Check edge midpoints
                select case(i)
                case(4)  ! Edge 1-2 midpoint
                    if (abs(x - 2.5_dp) > tol .or. abs(y - 2.5_dp) > tol) then
                        print *, "  ❌ Edge 1-2 midpoint mapping incorrect:", x, y
                        passed = .false.
                    end if
                case(5)  ! Edge 2-3 midpoint
                    if (abs(x - 3.0_dp) > tol .or. abs(y - 4.0_dp) > tol) then
                        print *, "  ❌ Edge 2-3 midpoint mapping incorrect:", x, y
                        passed = .false.
                    end if
                case(6)  ! Edge 3-1 midpoint
                    if (abs(x - 1.5_dp) > tol .or. abs(y - 3.5_dp) > tol) then
                        print *, "  ❌ Edge 3-1 midpoint mapping incorrect:", x, y
                        passed = .false.
                    end if
                end select
            end if
        end do
        
        if (passed) then
            print *, "  ✅ Coordinate transformation test passed"
        else
            all_passed = .false.
        end if
    end subroutine test_coordinate_transformation
    
    subroutine test_hessian_symmetry()
        type(basis_p2_2d_t) :: basis
        real(dp) :: hess(2,2)
        real(dp), parameter :: tol = 1e-14_dp
        logical :: passed = .true.
        integer :: i
        real(dp) :: xi, eta
        
        print *, "Testing Hessian symmetry..."
        
        xi = 0.3_dp
        eta = 0.2_dp
        
        do i = 1, 6
            hess = basis%hessian(i, xi, eta)
            
            ! Check symmetry: H(1,2) = H(2,1)
            if (abs(hess(1,2) - hess(2,1)) > tol) then
                print *, "  ❌ Hessian not symmetric for phi_", i
                print *, "    H(1,2) =", hess(1,2), "H(2,1) =", hess(2,1)
                passed = .false.
            end if
        end do
        
        if (passed) then
            print *, "  ✅ Hessian symmetry test passed"
        else
            all_passed = .false.
        end if
    end subroutine test_hessian_symmetry
    
    subroutine test_interpolation_property()
        type(basis_p2_2d_t) :: basis
        real(dp) :: vertices(2,3)
        real(dp) :: f_values(6), f_interp
        real(dp) :: xi, eta, x, y
        real(dp), parameter :: tol = 1e-14_dp
        logical :: passed = .true.
        integer :: i, j
        real(dp) :: nodes_physical(2,6)
        
        print *, "Testing interpolation property..."
        
        ! Test triangle
        vertices(:,1) = [0.0_dp, 0.0_dp]
        vertices(:,2) = [1.0_dp, 0.0_dp]
        vertices(:,3) = [0.0_dp, 1.0_dp]
        
        ! Physical nodes
        nodes_physical(:,1:3) = vertices
        nodes_physical(:,4) = 0.5_dp * (vertices(:,1) + vertices(:,2))
        nodes_physical(:,5) = 0.5_dp * (vertices(:,2) + vertices(:,3))
        nodes_physical(:,6) = 0.5_dp * (vertices(:,3) + vertices(:,1))
        
        ! Function values at nodes: f(x,y) = sin(pi*x) * cos(pi*y)
        do i = 1, 6
            x = nodes_physical(1,i)
            y = nodes_physical(2,i)
            f_values(i) = sin(3.14159265359_dp * x) * cos(3.14159265359_dp * y)
        end do
        
        ! Check interpolation at nodes
        do j = 1, 6
            xi = basis%nodes(1, j)
            eta = basis%nodes(2, j)
            
            f_interp = 0.0_dp
            do i = 1, 6
                f_interp = f_interp + f_values(i) * basis%eval(i, xi, eta)
            end do
            
            if (abs(f_interp - f_values(j)) > tol) then
                print *, "  ❌ Interpolation failed at node", j
                print *, "    Expected:", f_values(j), "Got:", f_interp
                passed = .false.
            end if
        end do
        
        if (passed) then
            print *, "  ✅ Interpolation property test passed"
        else
            all_passed = .false.
        end if
    end subroutine test_interpolation_property

end program test_basis_p2_2d