module fortfem_triangle_lagrange_arbitrary_order
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve
    use fortfem_generated_barycentric_jet2, only: generated_barycentric_jet2
    use fortfem_generated_product_jet1_order1, only: generated_product_jet1_order1
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
        real(dp), allocatable :: coefficients(:, :)
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

        real(dp), allocatable :: bernstein_values(:), gradients(:, :)
        real(dp), allocatable :: inverse(:, :), vandermonde(:, :)
        integer :: dof, first_power, info, second_power, third_power

        status = 1
        if (degree < 0) return

        basis%degree = degree
        basis%dof_count = (degree + 1) * (degree + 2) / 2
        allocate(basis%barycentric_powers(3, basis%dof_count))
        allocate(basis%coefficients(basis%dof_count, basis%dof_count))
        allocate(basis%nodes(2, basis%dof_count))
        allocate(bernstein_values(basis%dof_count))
        allocate(gradients(2, basis%dof_count))
        allocate(vandermonde(basis%dof_count, basis%dof_count))
        allocate(inverse(basis%dof_count, basis%dof_count))

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

        do dof = 1, basis%dof_count
            call evaluate_bernstein_polynomials( &
                basis, basis%nodes(1, dof), basis%nodes(2, dof), &
                bernstein_values, gradients)
            vandermonde(dof, :) = bernstein_values
        end do
        basis%coefficients = 0.0_dp
        do dof = 1, basis%dof_count
            basis%coefficients(dof, dof) = 1.0_dp
        end do
        call dense_solve(vandermonde, basis%coefficients, inverse, info)
        if (info /= 0) then
            status = 2
            return
        end if
        basis%coefficients = inverse
        status = 0
    end subroutine initialize_triangle_lagrange_basis

    subroutine evaluate_triangle_lagrange_basis( &
            basis, xi, eta, values, gradients, status)
        type(triangle_lagrange_basis_t), intent(in) :: basis
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(:), gradients(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: bernstein_gradients(:, :)
        real(dp), allocatable :: bernstein_values(:)
        integer :: basis_dof, modal_dof

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

        allocate(bernstein_values(basis%dof_count))
        allocate(bernstein_gradients(2, basis%dof_count))
        call evaluate_bernstein_polynomials( &
            basis, xi, eta, bernstein_values, bernstein_gradients)
        do basis_dof = 1, basis%dof_count
            do modal_dof = 1, basis%dof_count
                values(basis_dof) = values(basis_dof) + &
                    basis%coefficients(modal_dof, basis_dof) * &
                    bernstein_values(modal_dof)
                gradients(:, basis_dof) = gradients(:, basis_dof) + &
                    basis%coefficients(modal_dof, basis_dof) * &
                    bernstein_gradients(:, modal_dof)
            end do
        end do
        status = 0
    end subroutine evaluate_triangle_lagrange_basis

    subroutine triangle_lagrange_nodes(basis, nodes)
        type(triangle_lagrange_basis_t), intent(in) :: basis
        real(dp), allocatable, intent(out) :: nodes(:, :)

        allocate(nodes(2, basis%dof_count))
        nodes = basis%nodes
    end subroutine triangle_lagrange_nodes

    pure subroutine evaluate_bernstein_polynomials( &
            basis, xi, eta, values, gradients)
        type(triangle_lagrange_basis_t), intent(in) :: basis
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(:), gradients(:, :)

        real(dp) :: barycentric(3), barycentric_gradients(2, 3), coefficient
        real(dp) :: factors(3), derivatives(3), next_value, next_derivative
        integer :: dof, first_power, second_power, third_power, component, factor

        call generated_barycentric_jet2(xi, eta, &
            barycentric(1), barycentric_gradients(1, 1), barycentric_gradients(2, 1), &
            barycentric(2), barycentric_gradients(1, 2), barycentric_gradients(2, 2), &
            barycentric(3), barycentric_gradients(1, 3), barycentric_gradients(2, 3))
        do dof = 1, basis%dof_count
            first_power = basis%barycentric_powers(1, dof)
            second_power = basis%barycentric_powers(2, dof)
            third_power = basis%barycentric_powers(3, dof)
            coefficient = multinomial_coefficient( &
                basis%degree, first_power, second_power, third_power)
            factors = 1.0_dp
            factors(1) = coefficient
            derivatives = 0.0_dp
            do component = 1, 3
                do factor = 1, basis%barycentric_powers(component, dof)
                    call generated_product_jet1_order1(factors(component), &
                        derivatives(component), barycentric(component), 1.0_dp, &
                        next_value, next_derivative)
                    factors(component) = next_value
                    derivatives(component) = next_derivative
                end do
            end do
            call generated_factor_product_jet2_order1( &
                factors(1), derivatives(1), factors(2), derivatives(2), &
                factors(3), derivatives(3), values(dof), gradients(1, dof), &
                gradients(2, dof))
        end do
    end subroutine evaluate_bernstein_polynomials

    pure function multinomial_coefficient( &
            degree, first_power, second_power, third_power) result(value)
        integer, intent(in) :: degree, first_power, second_power, third_power
        real(dp) :: value

        value = factorial(degree) / ( &
            factorial(first_power) * factorial(second_power) * &
            factorial(third_power))
    end function multinomial_coefficient

    pure function factorial(n) result(value)
        integer, intent(in) :: n
        real(dp) :: value

        integer :: factor

        value = 1.0_dp
        do factor = 2, n
            value = value * real(factor, dp)
        end do
    end function factorial

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
        if (allocated(right%coefficients)) then
            allocate(left%coefficients( &
                size(right%coefficients, 1), size(right%coefficients, 2)))
            left%coefficients = right%coefficients
        end if
        if (allocated(right%nodes)) then
            allocate(left%nodes(size(right%nodes, 1), size(right%nodes, 2)))
            left%nodes = right%nodes
        end if
    end subroutine assign_triangle_lagrange_basis

end module fortfem_triangle_lagrange_arbitrary_order
