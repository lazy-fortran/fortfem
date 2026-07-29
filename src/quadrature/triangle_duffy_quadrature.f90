module fortfem_triangle_duffy_quadrature
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    private

    public :: triangle_duffy_quadrature

contains

    subroutine triangle_duffy_quadrature( &
            polynomial_degree, xi, eta, weights, status)
        integer, intent(in) :: polynomial_degree
        real(dp), allocatable, intent(out) :: xi(:), eta(:), weights(:)
        integer, intent(out) :: status

        real(dp), allocatable :: line_nodes(:), line_weights(:)
        integer :: first_node, node_count, point, second_node

        status = 1
        if (polynomial_degree < 0) return

        node_count = (polynomial_degree + 3) / 2
        allocate(line_nodes(node_count), line_weights(node_count))
        call gauss_legendre_ab( &
            node_count, 0.0_dp, 1.0_dp, line_nodes, line_weights)
        allocate(xi(node_count**2), eta(node_count**2), weights(node_count**2))

        point = 0
        do first_node = 1, node_count
            do second_node = 1, node_count
                point = point + 1
                xi(point) = line_nodes(first_node)
                eta(point) = (1.0_dp - line_nodes(first_node)) * &
                    line_nodes(second_node)
                weights(point) = line_weights(first_node) * &
                    line_weights(second_node) * &
                    (1.0_dp - line_nodes(first_node))
            end do
        end do
        status = 0
    end subroutine triangle_duffy_quadrature

end module fortfem_triangle_duffy_quadrature
