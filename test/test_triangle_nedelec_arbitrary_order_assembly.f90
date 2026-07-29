program test_triangle_nedelec_arbitrary_order_assembly
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_triangle_nedelec_curl_mass_element, &
        build_triangle_discrete_gradient, initialize_triangle_lagrange_basis, &
        triangle_lagrange_basis_t, triangle_lagrange_nodes
    use fortfem_kinds, only: dp
    implicit none

    type(triangle_lagrange_basis_t) :: lagrange_basis
    real(dp), allocatable :: dofs(:), gradient_matrix(:, :), matrix(:, :)
    real(dp), allocatable :: matrix_times_dofs(:), nodal_values(:), nodes(:, :)
    real(dp) :: energy, exact_energy
    real(dp) :: reference_vertices(2, 3), stretched_vertices(2, 3)
    integer :: order, status
    logical :: all_passed

    all_passed = .true.
    reference_vertices(:, 1) = [0.0_dp, 0.0_dp]
    reference_vertices(:, 2) = [1.0_dp, 0.0_dp]
    reference_vertices(:, 3) = [0.0_dp, 1.0_dp]
    stretched_vertices(:, 1) = [0.0_dp, 0.0_dp]
    stretched_vertices(:, 2) = [2.0_dp, 0.0_dp]
    stretched_vertices(:, 3) = [0.0_dp, 3.0_dp]

    call assemble_triangle_nedelec_curl_mass_element( &
        reference_vertices, 1, 2, matrix, status)
    call record_condition(status == 0, &
        "Order-one curl-plus-mass element assembly succeeds")
    call record_condition(abs(matrix(2, 2) - 13.0_dp / 6.0_dp) < 1.0e-13_dp, &
        "Element energy matches the exact field v=(-y,x)")
    deallocate(matrix)

    do order = 1, 4
        call initialize_triangle_lagrange_basis(order, lagrange_basis, status)
        call triangle_lagrange_nodes(lagrange_basis, nodes)
        allocate(nodal_values(size(nodes, 2)))
        nodal_values = nodes(1, :)**order
        call build_triangle_discrete_gradient(order, gradient_matrix, status)
        allocate(dofs(size(gradient_matrix, 1)))
        dofs = matmul(gradient_matrix, nodal_values)

        call assemble_triangle_nedelec_curl_mass_element( &
            reference_vertices, order, 2 * order, matrix, status)
        allocate(matrix_times_dofs(size(dofs)))
        matrix_times_dofs = matmul(matrix, dofs)
        energy = dot_product(dofs, matrix_times_dofs)
        exact_energy = real(order, dp) / real(2 * (2 * order - 1), dp)
        call record_condition(status == 0 .and. &
            size(matrix, 1) == order * (order + 2), &
            "Arbitrary-order element matrix has the exact trimmed dimension")
        call record_condition(maxval(abs(matrix - transpose(matrix))) < &
            2.0e-12_dp, "Arbitrary-order element matrix is symmetric")
        call record_condition(abs(energy - exact_energy) < 2.0e-9_dp, &
            "Element matrix integrates the exact gradient of x^order")

        deallocate( &
            dofs, gradient_matrix, matrix, matrix_times_dofs, nodal_values, &
            nodes)
    end do

    order = 4
    call initialize_triangle_lagrange_basis(order, lagrange_basis, status)
    call triangle_lagrange_nodes(lagrange_basis, nodes)
    allocate(nodal_values(size(nodes, 2)))
    nodal_values = (2.0_dp * nodes(1, :))**order
    call build_triangle_discrete_gradient(order, gradient_matrix, status)
    allocate(dofs(size(gradient_matrix, 1)))
    dofs = matmul(gradient_matrix, nodal_values)
    call assemble_triangle_nedelec_curl_mass_element( &
        stretched_vertices, order, 2 * order, matrix, status)
    allocate(matrix_times_dofs(size(dofs)))
    matrix_times_dofs = matmul(matrix, dofs)
    energy = dot_product(dofs, matrix_times_dofs)
    exact_energy = 6.0_dp * real(order * order, dp) * &
        2.0_dp**(2 * order - 2) / real((2 * order - 1) * (2 * order), dp)
    call record_condition(abs(energy - exact_energy) < 2.0e-7_dp, &
        "Covariant Piola assembly gives exact energy on a stretched triangle")
    deallocate( &
        dofs, gradient_matrix, matrix, matrix_times_dofs, nodal_values, nodes)

    call assemble_triangle_nedelec_curl_mass_element( &
        reference_vertices, 0, 2, matrix, status)
    call record_condition(status /= 0, &
        "Arbitrary-order element assembly rejects order zero")

    call check_summary("Arbitrary-order triangle Nedelec element assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_triangle_nedelec_arbitrary_order_assembly
