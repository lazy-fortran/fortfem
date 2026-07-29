program test_capi_mixed_assembly
    use check, only: check_condition, check_summary
    use iso_c_binding, only: c_double, c_int
    implicit none

    real(c_double), parameter :: vertices(8) = [ &
        0.0_c_double, 0.0_c_double, 1.0_c_double, 0.0_c_double, &
        1.0_c_double, 1.0_c_double, 0.0_c_double, 1.0_c_double]
    integer(c_int), parameter :: triangles(6) = [ &
        0_c_int, 1_c_int, 2_c_int, 0_c_int, 2_c_int, 3_c_int]
    integer(c_int) :: col_ptr(6), edges(10), row_ind(17)
    integer(c_int) :: triangle_edge_dofs(6), triangle_edge_signs(6)
    integer(c_int) :: handle, n_dofs, n_edges, nnz, status
    real(c_double) :: load(5), nedelec_x(5), nedelec_y(5)
    real(c_double) :: rt_field(5), values(17)
    real(c_double) :: edge_R, edge_Z
    integer :: column, edge, entry, row
    logical :: all_passed

    interface
        subroutine fortfem_triangle_mesh_create( &
                n_vertices, vertices, n_triangles, triangles, &
                handle, n_edges, status) &
                bind(c, name="fortfem_triangle_mesh_create")
            import c_double, c_int
            integer(c_int), value :: n_vertices, n_triangles
            real(c_double), intent(in) :: vertices(*)
            integer(c_int), intent(in) :: triangles(*)
            integer(c_int), intent(out) :: handle, n_edges, status
        end subroutine fortfem_triangle_mesh_create

        subroutine fortfem_triangle_mesh_edges( &
                handle, edge_capacity, n_edges, edges, &
                triangle_edge_dofs, triangle_edge_signs, status) &
                bind(c, name="fortfem_triangle_mesh_edges")
            import c_int
            integer(c_int), value :: handle, edge_capacity
            integer(c_int), intent(out) :: n_edges, edges(*)
            integer(c_int), intent(out) :: triangle_edge_dofs(*)
            integer(c_int), intent(out) :: triangle_edge_signs(*)
            integer(c_int), intent(out) :: status
        end subroutine fortfem_triangle_mesh_edges

        subroutine fortfem_nedelec_rt0_mass_csc( &
                handle, quadrature_degree, nnz_capacity, &
                n_dofs, nnz, col_ptr, row_ind, values, status) &
                bind(c, name="fortfem_nedelec_rt0_mass_csc")
            import c_double, c_int
            integer(c_int), value :: handle, quadrature_degree
            integer(c_int), value :: nnz_capacity
            integer(c_int), intent(out) :: n_dofs, nnz
            integer(c_int), intent(out) :: col_ptr(*), row_ind(*)
            real(c_double), intent(out) :: values(*)
            integer(c_int), intent(out) :: status
        end subroutine fortfem_nedelec_rt0_mass_csc

        subroutine fortfem_triangle_mesh_free(handle, status) &
                bind(c, name="fortfem_triangle_mesh_free")
            import c_int
            integer(c_int), value :: handle
            integer(c_int), intent(out) :: status
        end subroutine fortfem_triangle_mesh_free
    end interface

    all_passed = .true.
    call fortfem_triangle_mesh_create( &
        4_c_int, vertices, 2_c_int, triangles, handle, n_edges, status)
    call fortfem_triangle_mesh_edges( &
        handle, 5_c_int, n_edges, edges, triangle_edge_dofs, &
        triangle_edge_signs, status)
    call fortfem_nedelec_rt0_mass_csc( &
        handle, 2_c_int, 17_c_int, n_dofs, nnz, &
        col_ptr, row_ind, values, status)
    call record_condition(status == 0_c_int .and. nnz == 17_c_int, &
        "C mixed assembly exports compressed CSC storage")

    do edge = 1, 5
        edge_R = vertices(2 * edges(2 * edge) + 1) - &
            vertices(2 * edges(2 * edge - 1) + 1)
        edge_Z = vertices(2 * edges(2 * edge) + 2) - &
            vertices(2 * edges(2 * edge - 1) + 2)
        nedelec_x(edge) = edge_R
        nedelec_y(edge) = edge_Z
        rt_field(edge) = 2.0_c_double * edge_Z - 3.0_c_double * edge_R
    end do
    load = 0.0_c_double
    do column = 1, n_dofs
        do entry = col_ptr(column) + 1, col_ptr(column + 1)
            row = row_ind(entry) + 1
            load(row) = load(row) + values(entry) * rt_field(column)
        end do
    end do
    call record_condition( &
        abs(dot_product(nedelec_x, load) - 2.0_c_double) < 1.0e-13_c_double, &
        "C mixed operator integrates an exact constant x pairing")
    call record_condition( &
        abs(dot_product(nedelec_y, load) - 3.0_c_double) < 1.0e-13_c_double, &
        "C mixed operator integrates an exact constant y pairing")

    call fortfem_triangle_mesh_free(handle, status)
    call check_summary("C mixed Nedelec-RT0 assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_capi_mixed_assembly
