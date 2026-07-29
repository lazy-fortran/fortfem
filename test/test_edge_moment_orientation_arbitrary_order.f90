program test_edge_moment_orientation_arbitrary_order
    use check, only: check_condition, check_summary
    use fortfem_api, only: apply_edge_moment_orientation
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    integer, parameter :: order = 5
    real(dp) :: forward(order), recovered(order), reversed(order)
    real(dp) :: direct_reversed(order)
    real(dp) :: nodes(10), weights(10), trace
    integer :: moment, node, status
    logical :: all_passed

    all_passed = .true.
    call gauss_legendre_ab(10, 0.0_dp, 1.0_dp, nodes, weights)
    forward = 0.0_dp
    direct_reversed = 0.0_dp
    do moment = 0, order - 1
        do node = 1, size(nodes)
            trace = 1.0_dp + 2.0_dp * nodes(node) - &
                3.0_dp * nodes(node)**2 + 0.5_dp * nodes(node)**4
            forward(moment + 1) = forward(moment + 1) + weights(node) * &
                trace * shifted_legendre(moment, nodes(node))
            trace = 1.0_dp + 2.0_dp * (1.0_dp - nodes(node)) - &
                3.0_dp * (1.0_dp - nodes(node))**2 + &
                0.5_dp * (1.0_dp - nodes(node))**4
            direct_reversed(moment + 1) = &
                direct_reversed(moment + 1) - weights(node) * trace * &
                shifted_legendre(moment, nodes(node))
        end do
    end do

    call apply_edge_moment_orientation( &
        order, -1, forward, reversed, status)
    call record_condition(status == 0 .and. &
        maxval(abs(reversed - direct_reversed)) < 2.0e-14_dp, &
        "Alternating edge transform matches direct reversed moments")
    call apply_edge_moment_orientation( &
        order, -1, reversed, recovered, status)
    call record_condition(maxval(abs(recovered - forward)) < 2.0e-14_dp, &
        "Edge reversal transform is an involution")
    call apply_edge_moment_orientation( &
        order, 1, forward, recovered, status)
    call record_condition(maxval(abs(recovered - forward)) < 1.0e-15_dp, &
        "Aligned edge orientation leaves all moments unchanged")
    call apply_edge_moment_orientation( &
        order, 0, forward, recovered, status)
    call record_condition(status /= 0, &
        "Edge moment transform rejects an invalid orientation")

    call check_summary("Arbitrary-order edge moment orientation")
    if (.not. all_passed) error stop 1

contains

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

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_edge_moment_orientation_arbitrary_order
