program test_nedelec_sparse_assembly
    use check, only: check_condition, check_summary
    use fortsparse, only: csc_t, csc_matvec, fortsparse_status_t
    use fortfem_assembly_nedelec_2d, only: &
        assemble_nedelec_axisymmetric_fourier_csc
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    implicit none

    type(mesh_2d_t) :: mesh
    type(csc_t) :: matrix
    type(fortsparse_status_t) :: status
    real(dp), allocatable :: dofs(:), matrix_times_dofs(:)
    real(dp) :: energy
    integer :: dof_id, edge_id, vertex_a, vertex_b
    logical :: all_passed

    all_passed = .true.
    mesh%n_vertices = 4
    mesh%n_triangles = 2
    mesh%has_triangles = .true.
    allocate(mesh%vertices(2, 4), mesh%triangles(3, 2))
    mesh%vertices(:, 1) = [1.0_dp, 0.0_dp]
    mesh%vertices(:, 2) = [2.0_dp, 0.0_dp]
    mesh%vertices(:, 3) = [2.0_dp, 1.0_dp]
    mesh%vertices(:, 4) = [1.0_dp, 1.0_dp]
    mesh%triangles(:, 1) = [1, 2, 3]
    mesh%triangles(:, 2) = [1, 3, 4]
    call mesh%build_edge_connectivity()
    call mesh%build_edge_dof_numbering()

    call assemble_nedelec_axisymmetric_fourier_csc( &
        mesh, 0, 2, matrix, status)
    call record_condition(status%code == 0, &
        "Sparse weighted Nedelec assembly succeeds")
    call record_condition(matrix%nrow == 5 .and. matrix%ncol == 5, &
        "Sparse weighted Nedelec matrix uses global edge dimensions")
    call record_condition(matrix%nnz == 17, &
        "Sparse weighted Nedelec assembly compresses shared-edge entries")

    allocate(dofs(mesh%n_edges), matrix_times_dofs(mesh%n_edges))
    do edge_id = 1, mesh%n_edges
        vertex_a = mesh%edges(1, edge_id)
        vertex_b = mesh%edges(2, edge_id)
        dof_id = mesh%edge_to_dof(edge_id) + 1
        dofs(dof_id) = 0.25_dp * ( &
            -(mesh%vertices(2, vertex_a) + mesh%vertices(2, vertex_b)) * &
            (mesh%vertices(1, vertex_b) - mesh%vertices(1, vertex_a)) + &
            (mesh%vertices(1, vertex_a) + mesh%vertices(1, vertex_b)) * &
            (mesh%vertices(2, vertex_b) - mesh%vertices(2, vertex_a)))
    end do
    matrix_times_dofs = csc_matvec(matrix, dofs)
    energy = dot_product(dofs, matrix_times_dofs)
    call record_condition(abs(energy - 1.5_dp) < 1.0e-13_dp, &
        "Sparse operator reproduces exact integral of R curl(A)^2")

    call check_summary("Sparse weighted Nedelec assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_nedelec_sparse_assembly
