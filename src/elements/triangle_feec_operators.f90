module fortfem_triangle_feec_operators
    use fortfem_kinds, only: dp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        evaluate_triangle_lagrange_basis, initialize_triangle_lagrange_basis, &
        triangle_lagrange_basis_t, triangle_lagrange_nodes
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    private

    public :: build_triangle_discrete_gradient

contains

    subroutine build_triangle_discrete_gradient(order, matrix, status)
        integer, intent(in) :: order
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        type(triangle_lagrange_basis_t) :: lagrange_basis
        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: eta(:), gradients(:, :), nodes(:, :)
        real(dp), allocatable :: triangle_weights(:), values(:), xi(:)
        real(dp) :: edge_point(2), polynomial, tangent(2)
        integer :: component, edge, exponent, node, row
        integer :: scalar_dof, scalar_dof_count, total_degree
        integer :: x_degree, y_degree

        status = 1
        if (order < 1) return
        call initialize_triangle_lagrange_basis( &
            order, lagrange_basis, status)
        if (status /= 0) return
        call triangle_lagrange_nodes(lagrange_basis, nodes)
        scalar_dof_count = size(nodes, 2)
        allocate(matrix(order * (order + 2), scalar_dof_count))
        allocate(values(scalar_dof_count), gradients(2, scalar_dof_count))
        allocate(edge_nodes(order + 2), edge_weights(order + 2))
        call gauss_legendre_ab( &
            order + 2, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        matrix = 0.0_dp

        row = 0
        do edge = 1, 3
            do exponent = 0, order - 1
                row = row + 1
                do node = 1, size(edge_nodes)
                    call reference_edge( &
                        edge, edge_nodes(node), edge_point, tangent)
                    call evaluate_triangle_lagrange_basis( &
                        lagrange_basis, edge_point(1), edge_point(2), &
                        values, gradients, status)
                    if (status /= 0) return
                    polynomial = shifted_legendre(exponent, edge_nodes(node))
                    do scalar_dof = 1, scalar_dof_count
                        matrix(row, scalar_dof) = matrix(row, scalar_dof) + &
                            edge_weights(node) * polynomial * &
                            dot_product(gradients(:, scalar_dof), tangent)
                    end do
                end do
            end do
        end do

        call triangle_duffy_quadrature( &
            2 * order, xi, eta, triangle_weights, status)
        if (status /= 0) return
        do component = 1, 2
            do total_degree = 0, order - 2
                do x_degree = 0, total_degree
                    y_degree = total_degree - x_degree
                    row = row + 1
                    do node = 1, size(xi)
                        call evaluate_triangle_lagrange_basis( &
                            lagrange_basis, xi(node), eta(node), &
                            values, gradients, status)
                        if (status /= 0) return
                        polynomial = xi(node)**x_degree * eta(node)**y_degree
                        do scalar_dof = 1, scalar_dof_count
                            matrix(row, scalar_dof) = &
                                matrix(row, scalar_dof) + &
                                triangle_weights(node) * polynomial * &
                                gradients(component, scalar_dof)
                        end do
                    end do
                end do
            end do
        end do
        if (row /= size(matrix, 1)) then
            status = 2
            return
        end if
        status = 0
    end subroutine build_triangle_discrete_gradient

    pure subroutine reference_edge(edge, parameter, point, tangent)
        integer, intent(in) :: edge
        real(dp), intent(in) :: parameter
        real(dp), intent(out) :: point(2), tangent(2)

        select case (edge)
        case (1)
            point(1) = parameter
            point(2) = 0.0_dp
            tangent(1) = 1.0_dp
            tangent(2) = 0.0_dp
        case (2)
            point(1) = 1.0_dp - parameter
            point(2) = parameter
            tangent(1) = -1.0_dp
            tangent(2) = 1.0_dp
        case (3)
            point(1) = 0.0_dp
            point(2) = 1.0_dp - parameter
            tangent(1) = 0.0_dp
            tangent(2) = -1.0_dp
        end select
    end subroutine reference_edge

    pure function shifted_legendre(degree, parameter) result(value)
        integer, intent(in) :: degree
        real(dp), intent(in) :: parameter
        real(dp) :: value

        real(dp) :: current, previous, coordinate
        integer :: polynomial_degree

        coordinate = 2.0_dp * parameter - 1.0_dp
        if (degree == 0) then
            value = 1.0_dp
            return
        end if
        previous = 1.0_dp
        current = coordinate
        do polynomial_degree = 1, degree - 1
            value = ( &
                real(2 * polynomial_degree + 1, dp) * coordinate * current - &
                real(polynomial_degree, dp) * previous) / &
                real(polynomial_degree + 1, dp)
            previous = current
            current = value
        end do
        value = current
    end function shifted_legendre

end module fortfem_triangle_feec_operators
