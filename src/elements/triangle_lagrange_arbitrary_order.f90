module fortfem_triangle_lagrange_arbitrary_order
    use fortfem_kinds, only: dp
    use fortfem_generated_p1_basis, only: generated_p1_jet
    use fortfem_generated_p2_basis, only: generated_p2_jet
    use fortfem_generated_barycentric_jet2, only: generated_barycentric_jet2
    use fortfem_generated_cardinal_product_jet1_order1, only: &
        generated_cardinal_product_jet1_order1
    use fortfem_generated_factor_product_jet2_order1, only: &
        generated_factor_product_jet2_order1
    implicit none

    private

    public :: assignment(=)
    public :: evaluate_triangle_lagrange_basis
    public :: initialize_triangle_lagrange_basis
    public :: triangle_lagrange_basis_t
    public :: triangle_lagrange_nodes

    type :: triangle_lagrange_basis_t
        private
        integer :: degree = -1
        integer :: dof_count = 0
        integer, allocatable :: barycentric_powers(:, :)
        real(dp), allocatable :: nodes(:, :)
    end type triangle_lagrange_basis_t

    interface assignment(=)
        module procedure assign_triangle_lagrange_basis
    end interface

contains

    subroutine initialize_triangle_lagrange_basis(degree, basis, status)
        integer, intent(in) :: degree
        type(triangle_lagrange_basis_t), intent(out) :: basis
        integer, intent(out) :: status

        integer :: dof, first_power, second_power, third_power

        status = 1
        if (degree < 0) return

        basis%degree = degree
        basis%dof_count = (degree + 1) * (degree + 2) / 2
        allocate(basis%barycentric_powers(3, basis%dof_count))
        allocate(basis%nodes(2, basis%dof_count))

        if (degree == 0) then
            basis%barycentric_powers(:, 1) = 0
            basis%nodes(:, 1) = [1.0_dp / 3.0_dp, 1.0_dp / 3.0_dp]
        else
            dof = 0
            do first_power = 0, degree
                do second_power = 0, degree - first_power
                    third_power = degree - first_power - second_power
                    dof = dof + 1
                    basis%barycentric_powers(1, dof) = first_power
                    basis%barycentric_powers(2, dof) = second_power
                    basis%barycentric_powers(3, dof) = third_power
                    basis%nodes(1, dof) = &
                        real(second_power, dp) / real(degree, dp)
                    basis%nodes(2, dof) = &
                        real(third_power, dp) / real(degree, dp)
                end do
            end do
        end if

        status = 0
    end subroutine initialize_triangle_lagrange_basis

    subroutine evaluate_triangle_lagrange_basis( &
            basis, xi, eta, values, gradients, status)
        type(triangle_lagrange_basis_t), intent(in) :: basis
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(:), gradients(:, :)
        integer, intent(out) :: status

        real(dp) :: barycentric(3), barycentric_gradients(2, 3)
        real(dp) :: factors(3), derivatives(3), normalization
        real(dp) :: next_value, next_derivative, hessian(2, 2)
        integer, parameter :: p1_order(3) = [3, 2, 1]
        integer, parameter :: p2_order(6) = [3, 5, 2, 6, 4, 1]
        integer :: basis_dof, component, factor

        values = 0.0_dp
        gradients = 0.0_dp
        status = 1
        if (basis%degree < 0 .or. basis%dof_count < 1) return
        if (size(values) /= basis%dof_count) return
        if (size(gradients, 1) /= 2 .or. &
            size(gradients, 2) /= basis%dof_count) return
        if (xi < -64.0_dp * epsilon(1.0_dp) .or. &
            eta < -64.0_dp * epsilon(1.0_dp)) return
        if (xi + eta > 1.0_dp + 64.0_dp * epsilon(1.0_dp)) return

        ! Dispatch the common orders to their retained generated jets. The
        ! permutations translate vertex/edge numbering to barycentric order.
        select case (basis%degree)
        case (0)
            values = 1.0_dp
            status = 0
            return
        case (1)
            do basis_dof = 1, 3
                call generated_p1_jet(p1_order(basis_dof), xi, eta, &
                    values(basis_dof), gradients(:, basis_dof), hessian)
            end do
            status = 0
            return
        case (2)
            do basis_dof = 1, 6
                call generated_p2_jet(p2_order(basis_dof), xi, eta, &
                    values(basis_dof), gradients(:, basis_dof), hessian)
            end do
            status = 0
            return
        end select

        ! Equispaced simplex nodal polynomials are products of the cardinal
        ! factors (k lambda - j)/(j + 1). Reuse the generated jet also used
        ! by tetrahedra; no modal transform or per-query allocation is needed.
        call generated_barycentric_jet2(xi, eta, &
            barycentric(1), barycentric_gradients(1, 1), barycentric_gradients(2, 1), &
            barycentric(2), barycentric_gradients(1, 2), barycentric_gradients(2, 2), &
            barycentric(3), barycentric_gradients(1, 3), barycentric_gradients(2, 3))
        do basis_dof = 1, basis%dof_count
            factors = 1.0_dp
            derivatives = 0.0_dp
            do component = 1, 3
                normalization = 1.0_dp
                do factor = 0, basis%barycentric_powers(component, basis_dof) - 1
                    normalization = normalization*real(factor + 1, dp)
                    call generated_cardinal_product_jet1_order1( &
                        factors(component), derivatives(component), &
                        real(basis%degree, dp), real(factor, dp), &
                        barycentric(component), next_value, next_derivative)
                    factors(component) = next_value
                    derivatives(component) = next_derivative
                end do
                factors(component) = factors(component)/normalization
                derivatives(component) = derivatives(component)/normalization
            end do
            call generated_factor_product_jet2_order1( &
                factors(1), derivatives(1), factors(2), derivatives(2), &
                factors(3), derivatives(3), values(basis_dof), &
                gradients(1, basis_dof), gradients(2, basis_dof))
        end do
        status = 0
    end subroutine evaluate_triangle_lagrange_basis

    subroutine triangle_lagrange_nodes(basis, nodes)
        type(triangle_lagrange_basis_t), intent(in) :: basis
        real(dp), allocatable, intent(out) :: nodes(:, :)

        allocate(nodes(2, basis%dof_count))
        nodes = basis%nodes
    end subroutine triangle_lagrange_nodes

    subroutine assign_triangle_lagrange_basis(left, right)
        type(triangle_lagrange_basis_t), intent(out) :: left
        type(triangle_lagrange_basis_t), intent(in) :: right

        left%degree = right%degree
        left%dof_count = right%dof_count
        if (allocated(right%barycentric_powers)) then
            allocate(left%barycentric_powers( &
                size(right%barycentric_powers, 1), &
                size(right%barycentric_powers, 2)))
            left%barycentric_powers = right%barycentric_powers
        end if
        if (allocated(right%nodes)) then
            allocate(left%nodes(size(right%nodes, 1), size(right%nodes, 2)))
            left%nodes = right%nodes
        end if
    end subroutine assign_triangle_lagrange_basis

end module fortfem_triangle_lagrange_arbitrary_order
