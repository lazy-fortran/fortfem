module fortfem_tetra_duffy_quadrature
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    private

    public :: tetra_duffy_quadrature

contains

    subroutine tetra_duffy_quadrature( &
            polynomial_degree, x, y, z, weights, status)
        integer, intent(in) :: polynomial_degree
        real(dp), allocatable, intent(out) :: x(:), y(:), z(:), weights(:)
        integer, intent(out) :: status

        real(dp), allocatable :: line_nodes(:), line_weights(:)
        real(dp) :: first_complement, second_complement
        integer :: first_node, node_count, point, second_node, third_node

        status = 1
        if (polynomial_degree < 0) return

        node_count = (polynomial_degree + 4) / 2
        allocate(line_nodes(node_count), line_weights(node_count))
        call gauss_legendre_ab( &
            node_count, 0.0_dp, 1.0_dp, line_nodes, line_weights)
        allocate( &
            x(node_count**3), y(node_count**3), z(node_count**3), &
            weights(node_count**3))

        point = 0
        do first_node = 1, node_count
            first_complement = 1.0_dp - line_nodes(first_node)
            do second_node = 1, node_count
                second_complement = 1.0_dp - line_nodes(second_node)
                do third_node = 1, node_count
                    point = point + 1
                    x(point) = line_nodes(first_node)
                    y(point) = first_complement * line_nodes(second_node)
                    z(point) = first_complement * second_complement * &
                        line_nodes(third_node)
                    weights(point) = line_weights(first_node) * &
                        line_weights(second_node) * line_weights(third_node) * &
                        first_complement**2 * second_complement
                end do
            end do
        end do
        status = 0
    end subroutine tetra_duffy_quadrature

end module fortfem_tetra_duffy_quadrature
