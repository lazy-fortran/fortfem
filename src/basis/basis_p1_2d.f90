module basis_p1_2d_module
    use fortfem_kinds
    use fortfem_generated_p1_basis, only: generated_p1_jet
    use fortfem_generated_affine_triangle_geometry, only: generated_affine_triangle_geometry
    implicit none
    private

    public :: basis_p1_2d_t

    type :: basis_p1_2d_t
        ! Reference triangle nodes
        real(dp) :: nodes(2,3) = reshape([ &
            0.0_dp, 0.0_dp, & ! Node 1
            1.0_dp, 0.0_dp, & ! Node 2
            0.0_dp, 1.0_dp  & ! Node 3
            ], [2, 3])
    contains
        procedure :: eval
        procedure :: grad
        procedure :: transform_to_physical
        procedure :: compute_jacobian
    end type basis_p1_2d_t

contains

    pure function eval(this, i, xi, eta) result(val)
        class(basis_p1_2d_t), intent(in) :: this
        integer, intent(in) :: i
        real(dp), intent(in) :: xi, eta
        real(dp) :: val, g(2), h(2, 2)
        call generated_p1_jet(i, xi, eta, val, g, h)
    end function eval

    pure function grad(this, i, xi, eta) result(gradient)
        class(basis_p1_2d_t), intent(in) :: this
        integer, intent(in) :: i
        real(dp), intent(in) :: xi, eta
        real(dp) :: gradient(2), v, h(2, 2)
        call generated_p1_jet(i, xi, eta, v, gradient, h)
    end function grad

    pure subroutine transform_to_physical(this, xi, eta, vertices, x, y)
        class(basis_p1_2d_t), intent(in) :: this
        real(dp), intent(in) :: xi, eta
        real(dp), intent(in) :: vertices(2,3)
        real(dp), intent(out) :: x, y
        real(dp) :: mapped(2), jacobian(2,2), determinant
        call generated_affine_triangle_geometry(xi, eta, vertices, mapped, jacobian, determinant)
        x=mapped(1)
        y=mapped(2)
    end subroutine transform_to_physical

    pure subroutine compute_jacobian(this, vertices, jac, det_j)
        class(basis_p1_2d_t), intent(in) :: this
        real(dp), intent(in) :: vertices(2,3)
        real(dp), intent(out) :: jac(2,2), det_j
        real(dp) :: mapped(2)
        call generated_affine_triangle_geometry(0.0_dp, 0.0_dp, vertices, mapped, jac, det_j)
    end subroutine compute_jacobian

end module basis_p1_2d_module
