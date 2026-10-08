module basis_p2_2d_module
    use fortfem_kinds
    use fortfem_generated_p2_basis, only: generated_p2_jet
    implicit none
    private

    public :: basis_p2_2d_t

    type :: basis_p2_2d_t
        ! Reference triangle nodes for P2 elements
        ! Vertices: (0,0), (1,0), (0,1)
        ! Edge midpoints: (0.5,0), (0.5,0.5), (0,0.5)
        real(dp) :: nodes(2,6) = reshape([ &
            0.0_dp, 0.0_dp, & ! Node 1 (vertex)
            1.0_dp, 0.0_dp, & ! Node 2 (vertex)
            0.0_dp, 1.0_dp, & ! Node 3 (vertex)
            0.5_dp, 0.0_dp, & ! Node 4 (edge 1-2 midpoint)
            0.5_dp, 0.5_dp, & ! Node 5 (edge 2-3 midpoint)
            0.0_dp, 0.5_dp  & ! Node 6 (edge 3-1 midpoint)
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
        real(dp) :: val, g(2), h(2, 2)
        call generated_p2_jet(i, xi, eta, val, g, h)
    end function eval

    pure function grad(this, i, xi, eta) result(gradient)
        class(basis_p2_2d_t), intent(in) :: this
        integer, intent(in) :: i
        real(dp), intent(in) :: xi, eta
        real(dp) :: gradient(2), v, h(2, 2)
        call generated_p2_jet(i, xi, eta, v, gradient, h)
    end function grad

    pure function hessian(this, i, xi, eta) result(hess)
        class(basis_p2_2d_t), intent(in) :: this
        integer, intent(in) :: i
        real(dp), intent(in) :: xi, eta
        real(dp) :: hess(2, 2), v, g(2)
        call generated_p2_jet(i, xi, eta, v, g, hess)
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
        nodes_physical(:,4) = 0.5_dp * (vertices(:,1) + vertices(:,2)) ! Edge 1-2
        nodes_physical(:,5) = 0.5_dp * (vertices(:,2) + vertices(:,3)) ! Edge 2-3
        nodes_physical(:,6) = 0.5_dp * (vertices(:,3) + vertices(:,1)) ! Edge 3-1

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
        nodes_physical(:,4) = 0.5_dp * (vertices(:,1) + vertices(:,2)) ! Edge 1-2
        nodes_physical(:,5) = 0.5_dp * (vertices(:,2) + vertices(:,3)) ! Edge 2-3
        nodes_physical(:,6) = 0.5_dp * (vertices(:,3) + vertices(:,1)) ! Edge 3-1

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
