module basis_p2_2d_module
    use fortfem_kinds
    use fortfem_generated_p2_basis, only: generated_p2_jet
    use fortfem_generated_affine_triangle_geometry, only: generated_affine_triangle_geometry
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
        real(dp) :: mapped(2), jacobian(2,2), determinant
        call generated_affine_triangle_geometry(xi, eta, vertices, mapped, jacobian, determinant)
        x=mapped(1)
        y=mapped(2)
    end subroutine transform_to_physical

    pure subroutine compute_jacobian(this, vertices, jac, det_j)
        class(basis_p2_2d_t), intent(in) :: this
        real(dp), intent(in) :: vertices(2,3)
        real(dp), intent(out) :: jac(2,2), det_j
        real(dp) :: mapped(2)
        call generated_affine_triangle_geometry(0.0_dp, 0.0_dp, vertices, mapped, jac, det_j)
    end subroutine compute_jacobian

end module basis_p2_2d_module
