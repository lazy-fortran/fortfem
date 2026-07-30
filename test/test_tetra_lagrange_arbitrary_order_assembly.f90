program test_tetra_lagrange_arbitrary_order_assembly
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_tetra_lagrange_stiffness_csc, &
        assemble_tetra_lagrange_stiffness_element, &
        build_tetra_lagrange_dof_map, initialize_tetra_lagrange, &
        tetra_lagrange_nodes, tetra_lagrange_t
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    type(tetra_lagrange_t) :: basis
    type(csc_t) :: matrix
    type(fortsparse_status_t) :: sparse_status
    integer, allocatable :: global_dofs(:, :)
    real(dp), allocatable :: coefficients(:), element_matrix(:, :)
    real(dp), allocatable :: nodes(:, :), product(:)
    real(dp) :: mesh_vertices(3, 5), point(3), vertices(3, 4)
    integer :: cell, degree, global_count, local_dof, status
    integer :: tetrahedra(4, 2)
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.2_dp, -0.3_dp, 0.1_dp]
    vertices(:, 2) = [2.1_dp, -0.5_dp, 0.4_dp]
    vertices(:, 3) = [0.6_dp, 1.4_dp, -0.2_dp]
    vertices(:, 4) = [-0.1_dp, 0.2_dp, 1.6_dp]
    do degree = 1, 4
        call assemble_tetra_lagrange_stiffness_element( &
            vertices, degree, 2*degree + 2, element_matrix, status)
        call record_condition(status == 0 .and. &
            maxval(abs(element_matrix - transpose(element_matrix))) < &
            2.0e-11_dp, &
            "Tetrahedral H1 element assembly is symmetric at every order")
    end do

    mesh_vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    mesh_vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    mesh_vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
    mesh_vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
    mesh_vertices(:, 5) = [0.0_dp, 0.0_dp, -1.0_dp]
    tetrahedra(:, 1) = [1, 2, 3, 4]
    tetrahedra(:, 2) = [1, 3, 2, 5]
    call build_tetra_lagrange_dof_map( &
        1, tetrahedra, global_dofs, global_count, status)
    call initialize_tetra_lagrange(1, basis, status)
    call tetra_lagrange_nodes(basis, nodes)
    allocate(coefficients(global_count))
    coefficients = 0.0_dp
    do cell = 1, 2
        vertices = mesh_vertices(:, tetrahedra(:, cell))
        do local_dof = 1, size(nodes, 2)
            point = vertices(:, 1) + &
                nodes(1, local_dof)*(vertices(:, 2) - vertices(:, 1)) + &
                nodes(2, local_dof)*(vertices(:, 3) - vertices(:, 1)) + &
                nodes(3, local_dof)*(vertices(:, 4) - vertices(:, 1))
            coefficients(global_dofs(local_dof, cell)) = point(1)
        end do
    end do
    call assemble_tetra_lagrange_stiffness_csc( &
        mesh_vertices, tetrahedra, 1, 4, matrix, sparse_status)
    product = csc_matvec(matrix, coefficients)
    call record_condition(sparse_status%code == 0 .and. &
        abs(dot_product(coefficients, product) - 1.0_dp/3.0_dp) < &
        2.0e-12_dp, &
        "Global tetrahedral H1 assembly gives exact linear-field energy")

    call check_summary("Tetrahedral H1 arbitrary-order assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_lagrange_arbitrary_order_assembly
