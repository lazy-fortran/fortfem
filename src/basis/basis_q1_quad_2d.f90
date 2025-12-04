module basis_q1_quad_2d_module
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: q1_shape_functions
    public :: q1_shape_derivatives
    public :: q1_jacobian
    public :: q1_reference_to_physical

contains

    pure subroutine q1_shape_functions(xi, eta, N)
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: N(4)
        ! Q1 bilinear shape functions on reference square [-1,1] x [-1,1]
        N(1) = 0.25_dp * (1.0_dp - xi) * (1.0_dp - eta)
        N(2) = 0.25_dp * (1.0_dp + xi) * (1.0_dp - eta)
        N(3) = 0.25_dp * (1.0_dp + xi) * (1.0_dp + eta)
        N(4) = 0.25_dp * (1.0_dp - xi) * (1.0_dp + eta)
    end subroutine q1_shape_functions

    pure subroutine q1_shape_derivatives(xi, eta, dN_dxi, dN_deta)
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: dN_dxi(4), dN_deta(4)

        dN_dxi(1) = -0.25_dp * (1.0_dp - eta)
        dN_dxi(2) =  0.25_dp * (1.0_dp - eta)
        dN_dxi(3) =  0.25_dp * (1.0_dp + eta)
        dN_dxi(4) = -0.25_dp * (1.0_dp + eta)

        dN_deta(1) = -0.25_dp * (1.0_dp - xi)
        dN_deta(2) = -0.25_dp * (1.0_dp + xi)
        dN_deta(3) =  0.25_dp * (1.0_dp + xi)
        dN_deta(4) =  0.25_dp * (1.0_dp - xi)
    end subroutine q1_shape_derivatives

    pure subroutine q1_jacobian(xi, eta, coords, jac, det_jac, inv_jac, success)
        real(dp), intent(in) :: xi, eta
        real(dp), intent(in) :: coords(2,4)
        real(dp), intent(out) :: jac(2,2)
        real(dp), intent(out) :: det_jac
        real(dp), intent(out) :: inv_jac(2,2)
        logical, intent(out) :: success

        real(dp) :: dN_dxi(4), dN_deta(4)

        call q1_shape_derivatives(xi, eta, dN_dxi, dN_deta)

        ! Jacobian of mapping x(xi,eta) = sum_i N_i(xi,eta) * x_i
        jac(1,1) = sum(dN_dxi  * coords(1,:))   ! dx/dxi
        jac(1,2) = sum(dN_deta * coords(1,:))   ! dx/deta
        jac(2,1) = sum(dN_dxi  * coords(2,:))   ! dy/dxi
        jac(2,2) = sum(dN_deta * coords(2,:))   ! dy/deta

        det_jac = jac(1,1) * jac(2,2) - jac(1,2) * jac(2,1)

        success = abs(det_jac) > 1.0e-12_dp
        if (success) then
            inv_jac(1,1) =  jac(2,2) / det_jac
            inv_jac(1,2) = -jac(1,2) / det_jac
            inv_jac(2,1) = -jac(2,1) / det_jac
            inv_jac(2,2) =  jac(1,1) / det_jac
        else
            inv_jac = 0.0_dp
        end if
    end subroutine q1_jacobian

    pure subroutine q1_reference_to_physical(xi_ref, eta_ref, coords, x_phys, y_phys)
        real(dp), intent(in) :: xi_ref, eta_ref
        real(dp), intent(in) :: coords(2,4)
        real(dp), intent(out) :: x_phys, y_phys

        real(dp) :: N(4)

        call q1_shape_functions(xi_ref, eta_ref, N)
        x_phys = sum(N * coords(1,:))
        y_phys = sum(N * coords(2,:))
    end subroutine q1_reference_to_physical

end module basis_q1_quad_2d_module

