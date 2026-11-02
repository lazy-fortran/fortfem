module fortfem_triangle_bernstein_jets
    !! Degree dispatch for the generated Bernstein jets on the reference
    !! triangle. Index order: beta_2 (power of xi) outer, beta_3 (power of eta)
    !! inner; hessian rows are (xx, xy, yy).
    use fortfem_kinds, only: dp
    use fortfem_generated_triangle_bernstein_jet_0, only: &
        generated_triangle_bernstein_jet_0
    use fortfem_generated_triangle_bernstein_jet_1, only: &
        generated_triangle_bernstein_jet_1
    use fortfem_generated_triangle_bernstein_jet_2, only: &
        generated_triangle_bernstein_jet_2
    use fortfem_generated_triangle_bernstein_jet_3, only: &
        generated_triangle_bernstein_jet_3
    use fortfem_generated_triangle_bernstein_jet_4, only: &
        generated_triangle_bernstein_jet_4
    use fortfem_generated_triangle_bernstein_jet_5, only: &
        generated_triangle_bernstein_jet_5
    use fortfem_generated_triangle_bernstein_jet_6, only: &
        generated_triangle_bernstein_jet_6
    use fortfem_generated_triangle_bernstein_jet_7, only: &
        generated_triangle_bernstein_jet_7
    use fortfem_generated_triangle_bernstein_jet_8, only: &
        generated_triangle_bernstein_jet_8
    implicit none

    private

    public :: triangle_bernstein_jet

contains

    pure subroutine triangle_bernstein_jet(degree, xi, eta, value, gradient, &
            hessian)
        integer, intent(in) :: degree
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: value(:), gradient(:, :), hessian(:, :)

        real(dp) :: lambda_1

        lambda_1 = 1.0_dp - xi - eta
        select case (degree)
        case (0)
            call generated_triangle_bernstein_jet_0(lambda_1, xi, eta, value, &
                gradient, hessian)
        case (1)
            call generated_triangle_bernstein_jet_1(lambda_1, xi, eta, value, &
                gradient, hessian)
        case (2)
            call generated_triangle_bernstein_jet_2(lambda_1, xi, eta, value, &
                gradient, hessian)
        case (3)
            call generated_triangle_bernstein_jet_3(lambda_1, xi, eta, value, &
                gradient, hessian)
        case (4)
            call generated_triangle_bernstein_jet_4(lambda_1, xi, eta, value, &
                gradient, hessian)
        case (5)
            call generated_triangle_bernstein_jet_5(lambda_1, xi, eta, value, &
                gradient, hessian)
        case (6)
            call generated_triangle_bernstein_jet_6(lambda_1, xi, eta, value, &
                gradient, hessian)
        case (7)
            call generated_triangle_bernstein_jet_7(lambda_1, xi, eta, value, &
                gradient, hessian)
        case (8)
            call generated_triangle_bernstein_jet_8(lambda_1, xi, eta, value, &
                gradient, hessian)
        case default
            value = 0.0_dp
            gradient = 0.0_dp
            hessian = 0.0_dp
        end select
    end subroutine triangle_bernstein_jet

end module fortfem_triangle_bernstein_jets
