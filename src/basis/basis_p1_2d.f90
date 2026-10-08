module basis_p1_2d_module
    use fortfem_kinds
    use fortfem_generated_p1_basis, only: generated_p1_jet
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
        integer :: i

        x = 0.0_dp
        y = 0.0_dp

        ! Linear transformation using basis functions
        do i = 1, 3
            x = x + vertices(1,i) * this%eval(i, xi, eta)
            y = y + vertices(2,i) * this%eval(i, xi, eta)
        end do

    end subroutine transform_to_physical

    pure subroutine compute_jacobian(this, vertices, jac, det_j)
        class(basis_p1_2d_t), intent(in) :: this
        real(dp), intent(in) :: vertices(2,3)
        real(dp), intent(out) :: jac(2,2)
        real(dp), intent(out) :: det_j
        integer :: i
        real(dp) :: grad_ref(2)

        ! Initialize Jacobian
        jac = 0.0_dp

        ! Jacobian of transformation
        ! J = sum_i vertices_i * grad(phi_i)^T
        do i = 1, 3
            grad_ref = this%grad(i, 0.0_dp, 0.0_dp) ! Constant for P1
            jac(1,1) = jac(1,1) + vertices(1,i) * grad_ref(1)
            jac(1,2) = jac(1,2) + vertices(1,i) * grad_ref(2)
            jac(2,1) = jac(2,1) + vertices(2,i) * grad_ref(1)
            jac(2,2) = jac(2,2) + vertices(2,i) * grad_ref(2)
        end do

        ! Determinant
        det_j = jac(1,1) * jac(2,2) - jac(1,2) * jac(2,1)

    end subroutine compute_jacobian

end module basis_p1_2d_module
