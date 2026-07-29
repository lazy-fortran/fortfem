program test_triangle_nedelec_commuting_projection
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_triangle_discrete_gradient, &
        initialize_triangle_lagrange_basis, interpolate_triangle_nedelec, &
        triangle_lagrange_basis_t, triangle_lagrange_nodes
    use fortfem_kinds, only: dp
    implicit none

    type(triangle_lagrange_basis_t) :: lagrange_basis
    real(dp), allocatable :: direct_dofs(:), gradient_matrix(:, :)
    real(dp), allocatable :: nodal_values(:), nodes(:, :), projected_dofs(:)
    real(dp) :: physical_point(2), vertices(2, 3)
    integer :: active_order, node, order, status
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.2_dp, -0.1_dp]
    vertices(:, 2) = [1.7_dp, 0.3_dp]
    vertices(:, 3) = [-0.2_dp, 1.4_dp]

    do order = 1, 4
        active_order = order
        call initialize_triangle_lagrange_basis( &
            order, lagrange_basis, status)
        call triangle_lagrange_nodes(lagrange_basis, nodes)
        allocate(nodal_values(size(nodes, 2)))
        do node = 1, size(nodes, 2)
            physical_point = vertices(:, 1) + &
                (vertices(:, 2) - vertices(:, 1)) * nodes(1, node) + &
                (vertices(:, 3) - vertices(:, 1)) * nodes(2, node)
            nodal_values(node) = manufactured_potential( &
                physical_point(1), physical_point(2), order)
        end do
        call build_triangle_discrete_gradient( &
            order, gradient_matrix, status)
        allocate(projected_dofs(size(gradient_matrix, 1)))
        projected_dofs = matmul(gradient_matrix, nodal_values)
        call interpolate_triangle_nedelec( &
            vertices, order, 2 * order, manufactured_gradient, &
            direct_dofs, status)
        call record_condition(status == 0 .and. &
            size(direct_dofs) == order * (order + 2), &
            "Physical Nedelec interpolation has the exact trimmed dimension")
        call record_condition(maxval(abs( &
            direct_dofs - projected_dofs)) < 2.0e-9_dp, &
            "Nedelec projection commutes with the affine discrete gradient")
        deallocate( &
            direct_dofs, gradient_matrix, nodal_values, nodes, projected_dofs)
    end do

    call interpolate_triangle_nedelec( &
        vertices, 0, 2, manufactured_gradient, direct_dofs, status)
    call record_condition(status /= 0, &
        "Physical Nedelec interpolation rejects order zero")

    call check_summary("Triangle Nedelec commuting projection")
    if (.not. all_passed) error stop 1

contains

    pure function manufactured_potential(x, y, order) result(value)
        real(dp), intent(in) :: x, y
        integer, intent(in) :: order
        real(dp) :: value

        value = x**order + x * y**(order - 1)
    end function manufactured_potential

    subroutine manufactured_gradient(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        value(1) = real(active_order, dp) * x**(active_order - 1) + &
            y**(active_order - 1)
        value(2) = 0.0_dp
        if (active_order > 1) then
            value(2) = real(active_order - 1, dp) * &
                x * y**(active_order - 2)
        end if
    end subroutine manufactured_gradient

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_triangle_nedelec_commuting_projection
