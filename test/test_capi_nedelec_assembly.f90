program test_capi_nedelec_assembly
    use check, only: check_condition, check_summary
    use iso_c_binding, only: c_double, c_int
    implicit none

    real(c_double), parameter :: vertices(8) = [ &
        1.0_c_double, 0.0_c_double, 2.0_c_double, 0.0_c_double, &
        2.0_c_double, 1.0_c_double, 1.0_c_double, 1.0_c_double]
    integer(c_int), parameter :: triangles(6) = [ &
        0_c_int, 1_c_int, 2_c_int, 0_c_int, 2_c_int, 3_c_int]
    integer(c_int) :: col_ptr(6), edges(10), row_ind(17)
    integer(c_int) :: triangle_edge_dofs(6), triangle_edge_signs(6)
    integer(c_int) :: handle, n_dofs, n_edges, nnz, status
    real(c_double) :: dofs(5), matrix_times_dofs(5), values(17)
    real(c_double) :: energy
    integer :: column, edge, entry, row, vertex_a, vertex_b
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

        subroutine fortfem_nedelec_axisymmetric_fourier_csc( &
                handle, mode, quadrature_degree, nnz_capacity, &
                n_dofs, nnz, col_ptr, row_ind, values, status) &
                bind(c, name="fortfem_nedelec_axisymmetric_fourier_csc")
            import c_double, c_int
            integer(c_int), value :: handle, mode, quadrature_degree
            integer(c_int), value :: nnz_capacity
            integer(c_int), intent(out) :: n_dofs, nnz
            integer(c_int), intent(out) :: col_ptr(*), row_ind(*)
            real(c_double), intent(out) :: values(*)
            integer(c_int), intent(out) :: status
        end subroutine fortfem_nedelec_axisymmetric_fourier_csc

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

    call fortfem_nedelec_axisymmetric_fourier_csc( &
        handle, 0_c_int, 2_c_int, 16_c_int, n_dofs, nnz, &
        col_ptr, row_ind, values, status)
    call record_condition(status == -2_c_int .and. &
        n_dofs == 5_c_int .and. nnz == 17_c_int, &
        "C sparse assembly reports the required capacity")

    call fortfem_nedelec_axisymmetric_fourier_csc( &
        handle, 0_c_int, 2_c_int, 17_c_int, n_dofs, nnz, &
        col_ptr, row_ind, values, status)
    call record_condition(status == 0_c_int .and. &
        col_ptr(1) == 0_c_int .and. col_ptr(6) == 17_c_int, &
        "C sparse assembly exports zero-based CSC storage")

    do edge = 1, 5
        vertex_a = edges(2 * edge - 1)
        vertex_b = edges(2 * edge)
        dofs(edge) = 0.25_c_double * ( &
            -(vertices(2 * vertex_a + 2) + vertices(2 * vertex_b + 2)) * &
            (vertices(2 * vertex_b + 1) - vertices(2 * vertex_a + 1)) + &
            (vertices(2 * vertex_a + 1) + vertices(2 * vertex_b + 1)) * &
            (vertices(2 * vertex_b + 2) - vertices(2 * vertex_a + 2)))
    end do
    matrix_times_dofs = 0.0_c_double
    do column = 1, n_dofs
        do entry = col_ptr(column) + 1, col_ptr(column + 1)
            row = row_ind(entry) + 1
            matrix_times_dofs(row) = matrix_times_dofs(row) + &
                values(entry) * dofs(column)
        end do
    end do
    energy = dot_product(dofs, matrix_times_dofs)
    call record_condition(abs(energy - 1.5_c_double) < 1.0e-13_c_double, &
        "C sparse operator reproduces exact integral of R curl(A)^2")

    call fortfem_triangle_mesh_free(handle, status)
    call check_summary("C weighted Nedelec assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_capi_nedelec_assembly
