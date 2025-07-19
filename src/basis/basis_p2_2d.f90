module basis_p2_2d_module
    use fortfem_kinds
    implicit none
    private
    
    public :: basis_p2_2d_t
    
    type :: basis_p2_2d_t
        ! Reference triangle nodes for P2 elements
        ! Vertices: (0,0), (1,0), (0,1)
        ! Edge midpoints: (0.5,0), (0.5,0.5), (0,0.5)
        real(dp) :: nodes(2,6) = reshape([ &
            0.0_dp, 0.0_dp, &  ! Node 1 (vertex)
            1.0_dp, 0.0_dp, &  ! Node 2 (vertex)
            0.0_dp, 1.0_dp, &  ! Node 3 (vertex)
            0.5_dp, 0.0_dp, &  ! Node 4 (edge 1-2 midpoint)
            0.5_dp, 0.5_dp, &  ! Node 5 (edge 2-3 midpoint)
            0.0_dp, 0.5_dp  &  ! Node 6 (edge 3-1 midpoint)
        ], [2, 6])
    contains
        procedure :: eval
        procedure :: grad
        procedure :: hessian
        procedure :: transform_to_physical
        procedure :: compute_jacobian
        procedure :: get_num_dofs
    end type basis_p2_2d_t
    
contains

    pure function get_num_dofs(this) result(n)
        class(basis_p2_2d_t), intent(in) :: this
        integer :: n
        n = 6
    end function get_num_dofs

    pure function eval(this, i, xi, eta) result(val)
        class(basis_p2_2d_t), intent(in) :: this
        integer, intent(in) :: i
        real(dp), intent(in) :: xi, eta
        real(dp) :: val
        real(dp) :: lambda1, lambda2, lambda3
        
        ! Barycentric coordinates
        lambda1 = 1.0_dp - xi - eta
        lambda2 = xi
        lambda3 = eta
        
        select case (i)
        case (1)
            ! Vertex 1: phi_1 = lambda1 * (2*lambda1 - 1)
            val = lambda1 * (2.0_dp * lambda1 - 1.0_dp)
        case (2)
            ! Vertex 2: phi_2 = lambda2 * (2*lambda2 - 1)
            val = lambda2 * (2.0_dp * lambda2 - 1.0_dp)
        case (3)
            ! Vertex 3: phi_3 = lambda3 * (2*lambda3 - 1)
            val = lambda3 * (2.0_dp * lambda3 - 1.0_dp)
        case (4)
            ! Edge 1-2 midpoint: phi_4 = 4 * lambda1 * lambda2
            val = 4.0_dp * lambda1 * lambda2
        case (5)
            ! Edge 2-3 midpoint: phi_5 = 4 * lambda2 * lambda3
            val = 4.0_dp * lambda2 * lambda3
        case (6)
            ! Edge 3-1 midpoint: phi_6 = 4 * lambda3 * lambda1
            val = 4.0_dp * lambda3 * lambda1
        case default
            val = 0.0_dp
        end select
        
    end function eval
    
    pure function grad(this, i, xi, eta) result(gradient)
        class(basis_p2_2d_t), intent(in) :: this
        integer, intent(in) :: i
        real(dp), intent(in) :: xi, eta
        real(dp) :: gradient(2)
        real(dp) :: lambda1, lambda2, lambda3
        real(dp) :: d_lambda1_dxi, d_lambda1_deta
        real(dp) :: d_lambda2_dxi, d_lambda2_deta
        real(dp) :: d_lambda3_dxi, d_lambda3_deta
        
        ! Barycentric coordinates
        lambda1 = 1.0_dp - xi - eta
        lambda2 = xi
        lambda3 = eta
        
        ! Gradients of barycentric coordinates
        d_lambda1_dxi = -1.0_dp
        d_lambda1_deta = -1.0_dp
        d_lambda2_dxi = 1.0_dp
        d_lambda2_deta = 0.0_dp
        d_lambda3_dxi = 0.0_dp
        d_lambda3_deta = 1.0_dp
        
        select case (i)
        case (1)
            ! grad(phi_1) = grad(lambda1 * (2*lambda1 - 1))
            gradient(1) = d_lambda1_dxi * (4.0_dp * lambda1 - 1.0_dp)
            gradient(2) = d_lambda1_deta * (4.0_dp * lambda1 - 1.0_dp)
        case (2)
            ! grad(phi_2) = grad(lambda2 * (2*lambda2 - 1))
            gradient(1) = d_lambda2_dxi * (4.0_dp * lambda2 - 1.0_dp)
            gradient(2) = d_lambda2_deta * (4.0_dp * lambda2 - 1.0_dp)
        case (3)
            ! grad(phi_3) = grad(lambda3 * (2*lambda3 - 1))
            gradient(1) = d_lambda3_dxi * (4.0_dp * lambda3 - 1.0_dp)
            gradient(2) = d_lambda3_deta * (4.0_dp * lambda3 - 1.0_dp)
        case (4)
            ! grad(phi_4) = grad(4 * lambda1 * lambda2)
            gradient(1) = 4.0_dp * (d_lambda1_dxi * lambda2 + lambda1 * d_lambda2_dxi)
            gradient(2) = 4.0_dp * (d_lambda1_deta * lambda2 + lambda1 * d_lambda2_deta)
        case (5)
            ! grad(phi_5) = grad(4 * lambda2 * lambda3)
            gradient(1) = 4.0_dp * (d_lambda2_dxi * lambda3 + lambda2 * d_lambda3_dxi)
            gradient(2) = 4.0_dp * (d_lambda2_deta * lambda3 + lambda2 * d_lambda3_deta)
        case (6)
            ! grad(phi_6) = grad(4 * lambda3 * lambda1)
            gradient(1) = 4.0_dp * (d_lambda3_dxi * lambda1 + lambda3 * d_lambda1_dxi)
            gradient(2) = 4.0_dp * (d_lambda3_deta * lambda1 + lambda3 * d_lambda1_deta)
        case default
            gradient = [0.0_dp, 0.0_dp]
        end select
        
    end function grad
    
    pure function hessian(this, i, xi, eta) result(hess)
        class(basis_p2_2d_t), intent(in) :: this
        integer, intent(in) :: i
        real(dp), intent(in) :: xi, eta
        real(dp) :: hess(2,2)
        real(dp) :: d2_lambda1_dxi2, d2_lambda1_dxideta, d2_lambda1_deta2
        real(dp) :: d2_lambda2_dxi2, d2_lambda2_dxideta, d2_lambda2_deta2
        real(dp) :: d2_lambda3_dxi2, d2_lambda3_dxideta, d2_lambda3_deta2
        
        ! Second derivatives of barycentric coordinates (all zero for linear functions)
        d2_lambda1_dxi2 = 0.0_dp
        d2_lambda1_dxideta = 0.0_dp
        d2_lambda1_deta2 = 0.0_dp
        d2_lambda2_dxi2 = 0.0_dp
        d2_lambda2_dxideta = 0.0_dp
        d2_lambda2_deta2 = 0.0_dp
        d2_lambda3_dxi2 = 0.0_dp
        d2_lambda3_dxideta = 0.0_dp
        d2_lambda3_deta2 = 0.0_dp
        
        select case (i)
        case (1)
            ! hess(phi_1) = hess(lambda1 * (2*lambda1 - 1))
            hess(1,1) = 4.0_dp * (-1.0_dp) * (-1.0_dp)  ! d2/dxi2
            hess(1,2) = 4.0_dp * (-1.0_dp) * (-1.0_dp)  ! d2/dxideta
            hess(2,1) = hess(1,2)                        ! d2/detadxi
            hess(2,2) = 4.0_dp * (-1.0_dp) * (-1.0_dp)  ! d2/deta2
        case (2)
            ! hess(phi_2) = hess(lambda2 * (2*lambda2 - 1))
            hess(1,1) = 4.0_dp * 1.0_dp * 1.0_dp       ! d2/dxi2
            hess(1,2) = 0.0_dp                          ! d2/dxideta
            hess(2,1) = 0.0_dp                          ! d2/detadxi
            hess(2,2) = 0.0_dp                          ! d2/deta2
        case (3)
            ! hess(phi_3) = hess(lambda3 * (2*lambda3 - 1))
            hess(1,1) = 0.0_dp                          ! d2/dxi2
            hess(1,2) = 0.0_dp                          ! d2/dxideta
            hess(2,1) = 0.0_dp                          ! d2/detadxi
            hess(2,2) = 4.0_dp * 1.0_dp * 1.0_dp       ! d2/deta2
        case (4)
            ! hess(phi_4) = hess(4 * lambda1 * lambda2)
            hess(1,1) = 0.0_dp                          ! d2/dxi2
            hess(1,2) = 4.0_dp * (-1.0_dp) * 1.0_dp    ! d2/dxideta
            hess(2,1) = hess(1,2)                       ! d2/detadxi
            hess(2,2) = 0.0_dp                          ! d2/deta2
        case (5)
            ! hess(phi_5) = hess(4 * lambda2 * lambda3)
            hess(1,1) = 0.0_dp                          ! d2/dxi2
            hess(1,2) = 4.0_dp * 1.0_dp * 1.0_dp       ! d2/dxideta
            hess(2,1) = hess(1,2)                       ! d2/detadxi
            hess(2,2) = 0.0_dp                          ! d2/deta2
        case (6)
            ! hess(phi_6) = hess(4 * lambda3 * lambda1)
            hess(1,1) = 0.0_dp                          ! d2/dxi2
            hess(1,2) = 4.0_dp * 1.0_dp * (-1.0_dp)    ! d2/dxideta
            hess(2,1) = hess(1,2)                       ! d2/detadxi
            hess(2,2) = 0.0_dp                          ! d2/deta2
        case default
            hess = 0.0_dp
        end select
        
    end function hessian
    
    pure subroutine transform_to_physical(this, xi, eta, vertices, x, y)
        class(basis_p2_2d_t), intent(in) :: this
        real(dp), intent(in) :: xi, eta
        real(dp), intent(in) :: vertices(2,3)
        real(dp), intent(out) :: x, y
        integer :: i
        real(dp) :: nodes_physical(2,6)
        
        ! Set vertex nodes
        nodes_physical(:,1:3) = vertices(:,1:3)
        
        ! Compute edge midpoints
        nodes_physical(:,4) = 0.5_dp * (vertices(:,1) + vertices(:,2))  ! Edge 1-2
        nodes_physical(:,5) = 0.5_dp * (vertices(:,2) + vertices(:,3))  ! Edge 2-3
        nodes_physical(:,6) = 0.5_dp * (vertices(:,3) + vertices(:,1))  ! Edge 3-1
        
        x = 0.0_dp
        y = 0.0_dp
        
        ! Quadratic transformation using basis functions
        do i = 1, 6
            x = x + nodes_physical(1,i) * this%eval(i, xi, eta)
            y = y + nodes_physical(2,i) * this%eval(i, xi, eta)
        end do
        
    end subroutine transform_to_physical
    
    pure subroutine compute_jacobian(this, vertices, jac, det_j)
        class(basis_p2_2d_t), intent(in) :: this
        real(dp), intent(in) :: vertices(2,3)
        real(dp), intent(out) :: jac(2,2)
        real(dp), intent(out) :: det_j
        integer :: i
        real(dp) :: grad_ref(2)
        real(dp) :: nodes_physical(2,6)
        
        ! Set vertex nodes
        nodes_physical(:,1:3) = vertices(:,1:3)
        
        ! Compute edge midpoints
        nodes_physical(:,4) = 0.5_dp * (vertices(:,1) + vertices(:,2))  ! Edge 1-2
        nodes_physical(:,5) = 0.5_dp * (vertices(:,2) + vertices(:,3))  ! Edge 2-3
        nodes_physical(:,6) = 0.5_dp * (vertices(:,3) + vertices(:,1))  ! Edge 3-1
        
        ! Initialize Jacobian
        jac = 0.0_dp
        
        ! Jacobian of transformation at reference element center (1/3, 1/3)
        ! J = sum_i nodes_i * grad(phi_i)^T
        do i = 1, 6
            grad_ref = this%grad(i, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp)
            jac(1,1) = jac(1,1) + nodes_physical(1,i) * grad_ref(1)
            jac(1,2) = jac(1,2) + nodes_physical(1,i) * grad_ref(2)
            jac(2,1) = jac(2,1) + nodes_physical(2,i) * grad_ref(1)
            jac(2,2) = jac(2,2) + nodes_physical(2,i) * grad_ref(2)
        end do
        
        ! Determinant
        det_j = jac(1,1) * jac(2,2) - jac(1,2) * jac(2,1)
        
    end subroutine compute_jacobian

end module basis_p2_2d_module