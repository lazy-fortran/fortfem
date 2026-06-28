module basis_p1_2d_module
    use fortfem_kinds
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
        real(dp) :: val

        select case (i)
        case (1)
            ! phi_1 = 1 - xi - eta
            val = 1.0_dp - xi - eta
        case (2)
            ! phi_2 = xi
            val = xi
        case (3)
            ! phi_3 = eta
            val = eta
        case default
            val = 0.0_dp
        end select

    end function eval

    pure function grad(this, i, xi, eta) result(gradient)
        class(basis_p1_2d_t), intent(in) :: this
        integer, intent(in) :: i
        real(dp), intent(in) :: xi, eta
        real(dp) :: gradient(2)

        ! P1 gradients are constant
        select case (i)
        case (1)
            ! grad(phi_1) = [-1, -1]
            gradient = [-1.0_dp, -1.0_dp]
        case (2)
            ! grad(phi_2) = [1, 0]
            gradient = [1.0_dp, 0.0_dp]
        case (3)
            ! grad(phi_3) = [0, 1]
            gradient = [0.0_dp, 1.0_dp]
        case default
            gradient = [0.0_dp, 0.0_dp]
        end select

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
