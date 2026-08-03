program test_tetra_nedelec_pml_sparse
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_tetra_nedelec_pml_csc, &
        build_tetra_edge_dof_map
    use fortfem_boundary, only: cartesian_curl_curl_pml_coefficients
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_matvec, csc_z_t, fortsparse_status_t
    implicit none

    type(csc_z_t) :: matrix
    type(fortsparse_status_t) :: sparse_status
    integer, allocatable :: edges(:, :), global_dofs(:, :), orientations(:, :)
    integer :: edge, order, status
    integer :: tetrahedra(4, 2)
    real(dp) :: vertices(3, 5), edge_vector(3), field(3)
    complex(dp), allocatable :: dofs(:), matrix_times_dofs(:)
    complex(dp) :: curl_coefficient(3), energy, expected
    complex(dp) :: mass_coefficient(3), stretch(3, 2)
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
    vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
    vertices(:, 5) = [0.0_dp, 0.0_dp, -1.0_dp]
    tetrahedra(:, 1) = [1, 2, 3, 4]
    tetrahedra(:, 2) = [1, 3, 2, 5]
    stretch(:, 1) = [ &
        cmplx(1.0_dp, 0.6_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
        cmplx(1.0_dp, 0.0_dp, dp)]
    stretch(:, 2) = stretch(:, 1)
    call build_tetra_edge_dof_map( &
        tetrahedra, edges, global_dofs, orientations, status)
    allocate(dofs(size(edges, 2)), matrix_times_dofs(size(edges, 2)))
    field = [1.0_dp, 2.0_dp, -1.0_dp]
    do edge = 1, size(edges, 2)
        edge_vector = vertices(:, edges(2, edge)) - &
            vertices(:, edges(1, edge))
        dofs(edge) = cmplx(dot_product(field, edge_vector), 0.0_dp, dp)
    end do
    call assemble_tetra_nedelec_pml_csc( &
        vertices, tetrahedra, 1, stretch, 1.3_dp, matrix, sparse_status)
    matrix_times_dofs = csc_matvec(matrix, dofs)
    energy = sum(dofs*matrix_times_dofs)
    call cartesian_curl_curl_pml_coefficients( &
        stretch(:, 1), curl_coefficient, mass_coefficient, status)
    expected = -1.3_dp**2/3.0_dp*sum(mass_coefficient*field**2)
    call record_condition(sparse_status%code == 0 .and. &
        abs(energy - expected) < 3.0e-12_dp, &
        "global PML mass has exact constant-field energy")

    do order = 2, 4
        call assemble_tetra_nedelec_pml_csc( &
            vertices, tetrahedra, order, stretch, 1.3_dp, matrix, sparse_status)
        call record_condition(sparse_status%code == 0 .and. matrix%nrow > 0, &
            "global complex PML assembly supports higher Nedelec order")
    end do
    call check_summary("Sparse arbitrary-order tetrahedral Nedelec PML")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_nedelec_pml_sparse
