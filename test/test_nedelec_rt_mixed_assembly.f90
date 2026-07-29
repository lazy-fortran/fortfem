program test_nedelec_rt_mixed_assembly
    use check, only: check_condition, check_summary
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    use fortfem_assembly_mixed_2d, only: assemble_nedelec_rt_mass_csc
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    implicit none

    type(mesh_2d_t) :: mesh
    type(csc_t) :: matrix
    type(fortsparse_status_t) :: status
    real(dp), allocatable :: nedelec_x(:), nedelec_y(:), rt_field(:)
    real(dp), allocatable :: load(:)
    real(dp) :: edge_vector(2), pairing_x, pairing_y
    integer :: dof_id, edge_id, vertex_a, vertex_b
    logical :: all_passed

    all_passed = .true.
    mesh%n_vertices = 4
    mesh%n_triangles = 2
    mesh%has_triangles = .true.
    allocate(mesh%vertices(2, 4), mesh%triangles(3, 2))
    mesh%vertices(:, 1) = [0.0_dp, 0.0_dp]
    mesh%vertices(:, 2) = [1.0_dp, 0.0_dp]
    mesh%vertices(:, 3) = [1.0_dp, 1.0_dp]
    mesh%vertices(:, 4) = [0.0_dp, 1.0_dp]
    mesh%triangles(:, 1) = [1, 2, 3]
    mesh%triangles(:, 2) = [1, 3, 4]
    call mesh%build_edge_connectivity()
    call mesh%build_edge_dof_numbering()

    call assemble_nedelec_rt_mass_csc(mesh, 2, matrix, status)
    call record_condition(status%code == 0 .and. matrix%nnz == 17, &
        "Mixed Nedelec-RT0 assembly compresses shared-edge entries")

    allocate(nedelec_x(mesh%n_edges), nedelec_y(mesh%n_edges))
    allocate(rt_field(mesh%n_edges), load(mesh%n_edges))
    do edge_id = 1, mesh%n_edges
        vertex_a = mesh%edges(1, edge_id)
        vertex_b = mesh%edges(2, edge_id)
        edge_vector = mesh%vertices(:, vertex_b) - &
            mesh%vertices(:, vertex_a)
        dof_id = mesh%edge_to_dof(edge_id) + 1
        nedelec_x(dof_id) = edge_vector(1)
        nedelec_y(dof_id) = edge_vector(2)
        rt_field(dof_id) = &
            2.0_dp * edge_vector(2) - 3.0_dp * edge_vector(1)
    end do

    load = csc_matvec(matrix, rt_field)
    pairing_x = dot_product(nedelec_x, load)
    pairing_y = dot_product(nedelec_y, load)
    call record_condition(abs(pairing_x - 2.0_dp) < 1.0e-13_dp, &
        "Mixed operator integrates the exact x component of a constant field")
    call record_condition(abs(pairing_y - 3.0_dp) < 1.0e-13_dp, &
        "Mixed operator integrates the exact y component of a constant field")

    call check_summary("Mixed Nedelec-RT0 assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_nedelec_rt_mixed_assembly
