module fortfem_tetra_nedelec_first_order
    use fortfem_kinds, only: dp
    implicit none

    private

    public :: evaluate_tetra_nedelec_first_order

contains

    pure subroutine evaluate_tetra_nedelec_first_order( &
            point, values, curls, status)
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: values(3, 6), curls(3, 6)
        integer, intent(out) :: status

        real(dp) :: gradients(3, 4), lambda(4), tolerance
        integer :: edge, edge_vertices(2, 6), first, second

        values = 0.0_dp
        curls = 0.0_dp
        status = 1
        tolerance = 64.0_dp * epsilon(1.0_dp)
        lambda = [ &
            1.0_dp - point(1) - point(2) - point(3), &
            point(1), point(2), point(3)]
        if (any(lambda < -tolerance)) return
        if (any(lambda > 1.0_dp + tolerance)) return

        gradients(:, 1) = [-1.0_dp, -1.0_dp, -1.0_dp]
        gradients(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
        gradients(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
        gradients(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
        edge_vertices(:, 1) = [1, 2]
        edge_vertices(:, 2) = [1, 3]
        edge_vertices(:, 3) = [1, 4]
        edge_vertices(:, 4) = [2, 3]
        edge_vertices(:, 5) = [2, 4]
        edge_vertices(:, 6) = [3, 4]

        do edge = 1, 6
            first = edge_vertices(1, edge)
            second = edge_vertices(2, edge)
            values(:, edge) = &
                lambda(first) * gradients(:, second) - &
                lambda(second) * gradients(:, first)
            curls(:, edge) = 2.0_dp * cross_product( &
                gradients(:, first), gradients(:, second))
        end do
        status = 0
    end subroutine evaluate_tetra_nedelec_first_order

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product(1) = first(2) * second(3) - first(3) * second(2)
        product(2) = first(3) * second(1) - first(1) * second(3)
        product(3) = first(1) * second(2) - first(2) * second(1)
    end function cross_product

end module fortfem_tetra_nedelec_first_order
