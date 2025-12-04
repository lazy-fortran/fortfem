program test_p2_basis_functions
    use fortfem_kinds
    use basis_p2_2d_module
    use check
    implicit none

    write(*,*) "Testing P2 basis function properties..."
    
    call test_p2_basis_evaluation()
    call test_p2_basis_gradients()
    call test_p2_partition_of_unity()
    call test_p2_lagrange_property()
    call test_p2_hessian_computation()
    
    call check_summary("P2 Basis Functions")

contains

    ! Test P2 basis function evaluation at specific points
    subroutine test_p2_basis_evaluation()
        type(basis_p2_2d_t) :: basis
        real(dp) :: val
        
        ! Test vertex basis functions at vertices
        ! φ₁(0,0) = 1, φ₁ at other vertices = 0
        val = basis%eval(1, 0.0_dp, 0.0_dp)
        call check_condition(abs(val - 1.0_dp) < 1.0e-14_dp, &
            "P2 evaluation: φ₁(0,0) = 1")
        
        val = basis%eval(1, 1.0_dp, 0.0_dp)
        call check_condition(abs(val - 0.0_dp) < 1.0e-14_dp, &
            "P2 evaluation: φ₁(1,0) = 0")
        
        val = basis%eval(1, 0.0_dp, 1.0_dp)
        call check_condition(abs(val - 0.0_dp) < 1.0e-14_dp, &
            "P2 evaluation: φ₁(0,1) = 0")
        
        ! φ₂(1,0) = 1, φ₂ at other vertices = 0
        val = basis%eval(2, 1.0_dp, 0.0_dp)
        call check_condition(abs(val - 1.0_dp) < 1.0e-14_dp, &
            "P2 evaluation: φ₂(1,0) = 1")
        
        val = basis%eval(2, 0.0_dp, 0.0_dp)
        call check_condition(abs(val - 0.0_dp) < 1.0e-14_dp, &
            "P2 evaluation: φ₂(0,0) = 0")
        
        ! Test edge basis functions at edge midpoints
        ! φ₄(0.5,0) = 1 (edge 1-2 midpoint)
        val = basis%eval(4, 0.5_dp, 0.0_dp)
        call check_condition(abs(val - 1.0_dp) < 1.0e-14_dp, &
            "P2 evaluation: φ₄(0.5,0) = 1")
        
        val = basis%eval(4, 0.0_dp, 0.0_dp)
        call check_condition(abs(val - 0.0_dp) < 1.0e-14_dp, &
            "P2 evaluation: φ₄(0,0) = 0")
        
        write(*,*) "   P2 basis evaluation tests passed"
    end subroutine test_p2_basis_evaluation

    ! Test P2 basis function gradients
    subroutine test_p2_basis_gradients()
        type(basis_p2_2d_t) :: basis
        real(dp) :: grad(2), expected_grad(2)
        
        ! Test gradient of φ₁ at center
        grad = basis%grad(1, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp)
        ! At center, λ₁ = λ₂ = λ₃ = 1/3
        ! φ₁ = λ₁(2λ₁ - 1) = (1/3)(2/3 - 1) = (1/3)(-1/3) = -1/9
        ! ∇φ₁ = (4λ₁ - 1)∇λ₁ = (4/3 - 1)[-1, -1] = (1/3)[-1, -1] = [-1/3, -1/3]
        expected_grad = [-1.0_dp/3.0_dp, -1.0_dp/3.0_dp]
        
        call check_condition(abs(grad(1) - expected_grad(1)) < 1.0e-12_dp, &
            "P2 gradients: φ₁ x-gradient at center")
        call check_condition(abs(grad(2) - expected_grad(2)) < 1.0e-12_dp, &
            "P2 gradients: φ₁ y-gradient at center")
        
        ! Test gradient of φ₂ at center
        grad = basis%grad(2, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp)
        ! ∇φ₂ = (4λ₂ - 1)∇λ₂ = (1/3)[1, 0] = [1/3, 0]
        expected_grad = [1.0_dp/3.0_dp, 0.0_dp]
        
        call check_condition(abs(grad(1) - expected_grad(1)) < 1.0e-12_dp, &
            "P2 gradients: φ₂ x-gradient at center")
        call check_condition(abs(grad(2) - expected_grad(2)) < 1.0e-12_dp, &
            "P2 gradients: φ₂ y-gradient at center")
        
        write(*,*) "   P2 gradient computation tests passed"
    end subroutine test_p2_basis_gradients

    ! Test partition of unity property for P2 elements
    subroutine test_p2_partition_of_unity()
        type(basis_p2_2d_t) :: basis
        real(dp) :: xi, eta, sum_basis
        integer :: i, j, k
        
        ! Test at several points that sum of basis functions = 1
        do i = 0, 5
            do j = 0, 5-i
                xi = real(i, dp) / 5.0_dp
                eta = real(j, dp) / 5.0_dp
                
                if (xi + eta <= 1.0_dp) then
                    sum_basis = 0.0_dp
                    do k = 1, 6
                        sum_basis = sum_basis + basis%eval(k, xi, eta)
                    end do
                    
                    call check_condition(abs(sum_basis - 1.0_dp) < 1.0e-14_dp, &
                        "P2 partition of unity: sum = 1")
                end if
            end do
        end do
        
        write(*,*) "   P2 partition of unity tests passed"
    end subroutine test_p2_partition_of_unity

    ! Test Lagrange property: φᵢ(xⱼ) = δᵢⱼ
    subroutine test_p2_lagrange_property()
        type(basis_p2_2d_t) :: basis
        real(dp) :: nodes(2,6)
        real(dp) :: val
        integer :: i, j
        
        ! P2 nodes: vertices + edge midpoints
        nodes(:,1) = [0.0_dp, 0.0_dp]  ! Vertex 1
        nodes(:,2) = [1.0_dp, 0.0_dp]  ! Vertex 2
        nodes(:,3) = [0.0_dp, 1.0_dp]  ! Vertex 3
        nodes(:,4) = [0.5_dp, 0.0_dp]  ! Edge 1-2 midpoint
        nodes(:,5) = [0.5_dp, 0.5_dp]  ! Edge 2-3 midpoint
        nodes(:,6) = [0.0_dp, 0.5_dp]  ! Edge 3-1 midpoint
        
        do i = 1, 6
            do j = 1, 6
                val = basis%eval(i, nodes(1,j), nodes(2,j))
                if (i == j) then
                    call check_condition(abs(val - 1.0_dp) < 1.0e-14_dp, &
                        "P2 Lagrange: φᵢ(xᵢ) = 1")
                else
                    call check_condition(abs(val - 0.0_dp) < 1.0e-14_dp, &
                        "P2 Lagrange: φᵢ(xⱼ) = 0 for i≠j")
                end if
            end do
        end do
        
        write(*,*) "   P2 Lagrange property tests passed"
    end subroutine test_p2_lagrange_property

    ! Test P2 Hessian computation
    subroutine test_p2_hessian_computation()
        type(basis_p2_2d_t) :: basis
        real(dp) :: hess(2,2)
        
        ! Test Hessian of φ₁ at a point
        hess = basis%hessian(1, 0.2_dp, 0.3_dp)
        
        ! For φ₁ = λ₁(2λ₁ - 1), the Hessian should be constant
        ! Since λ₁ = 1 - ξ - η, we have ∇²φ₁ = 4∇λ₁ ⊗ ∇λ₁
        call check_condition(abs(hess(1,1) - 4.0_dp) < 1.0e-14_dp, &
            "P2 Hessian: φ₁ ∂²/∂ξ²")
        call check_condition(abs(hess(1,2) - 4.0_dp) < 1.0e-14_dp, &
            "P2 Hessian: φ₁ ∂²/∂ξ∂η")
        call check_condition(abs(hess(2,2) - 4.0_dp) < 1.0e-14_dp, &
            "P2 Hessian: φ₁ ∂²/∂η²")
        
        write(*,*) "   P2 Hessian computation tests passed"
    end subroutine test_p2_hessian_computation

end program test_p2_basis_functions