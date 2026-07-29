program test_tetra_nedelec_sparse_assembly
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_tetra_nedelec_curl_mass_csc, &
        build_tetra_edge_dof_map
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    type(csc_t) :: matrix
    type(fortsparse_status_t) :: sparse_status
    integer, allocatable :: edges(:, :), global_dofs(:, :), orientations(:, :)
    integer :: edge, status, tetrahedra(4, 2)
    real(dp), allocatable :: dofs(:), matrix_times_dofs(:)
    real(dp) :: constant_field(3), edge_midpoint(3), edge_vector(3), energy
    real(dp) :: curl_field(3), vertices(3, 5)
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
    vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
    vertices(:, 5) = [0.0_dp, 0.0_dp, -1.0_dp]
    tetrahedra(:, 1) = [1, 2, 3, 4]
    tetrahedra(:, 2) = [1, 3, 2, 5]
    call build_tetra_edge_dof_map( &
        tetrahedra, edges, global_dofs, orientations, status)
    allocate(dofs(size(edges, 2)), matrix_times_dofs(size(edges, 2)))

    constant_field = [1.0_dp, 2.0_dp, -1.0_dp]
    do edge = 1, size(edges, 2)
        edge_vector = vertices(:, edges(2, edge)) - &
            vertices(:, edges(1, edge))
        dofs(edge) = dot_product(constant_field, edge_vector)
    end do
    call assemble_tetra_nedelec_curl_mass_csc( &
        vertices, tetrahedra, matrix, sparse_status, 4.0_dp, 1.5_dp)
    matrix_times_dofs = csc_matvec(matrix, dofs)
    energy = dot_product(dofs, matrix_times_dofs)
    call record_condition(sparse_status%code == 0 .and. &
        abs(energy - 3.0_dp) < 2.0e-12_dp, &
        "Tetrahedral mass assembly has exact constant-field energy")

    curl_field = [1.0_dp, -2.0_dp, 3.0_dp]
    do edge = 1, size(edges, 2)
        edge_vector = vertices(:, edges(2, edge)) - &
            vertices(:, edges(1, edge))
        edge_midpoint = 0.5_dp * ( &
            vertices(:, edges(1, edge)) + vertices(:, edges(2, edge)))
        dofs(edge) = dot_product( &
            0.5_dp * cross_product(curl_field, edge_midpoint), edge_vector)
    end do
    call assemble_tetra_nedelec_curl_mass_csc( &
        vertices, tetrahedra, matrix, sparse_status, 2.0_dp, 0.0_dp)
    matrix_times_dofs = csc_matvec(matrix, dofs)
    energy = dot_product(dofs, matrix_times_dofs)
    call record_condition(abs(energy - 28.0_dp / 3.0_dp) < 3.0e-12_dp, &
        "Tetrahedral curl assembly has exact rotational-field energy")

    call check_summary("Sparse tetrahedral Nedelec assembly")
    if (.not. all_passed) error stop 1

contains

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product(1) = first(2) * second(3) - first(3) * second(2)
        product(2) = first(3) * second(1) - first(1) * second(3)
        product(3) = first(1) * second(2) - first(2) * second(1)
    end function cross_product

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_nedelec_sparse_assembly
