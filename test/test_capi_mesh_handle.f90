program test_capi_mesh_handle
    use iso_c_binding, only: c_double, c_int
    use check, only: check_condition, check_summary
    implicit none

    real(c_double), parameter :: vertices(8) = [ &
        0.0_c_double, 0.0_c_double, 1.0_c_double, 0.0_c_double, &
        1.0_c_double, 1.0_c_double, 0.0_c_double, 1.0_c_double]
    integer(c_int), parameter :: triangles(6) = [ &
        0_c_int, 1_c_int, 2_c_int, 0_c_int, 2_c_int, 3_c_int]
    integer(c_int) :: edges(10), triangle_edge_dofs(6)
    integer(c_int) :: triangle_edge_signs(6)
    integer(c_int) :: handle, n_edges, status
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
    call record_condition(status == 0_c_int .and. handle > 0_c_int, &
        "C triangle mesh creation returns a live handle")
    call record_condition(n_edges == 5_c_int, &
        "C triangle mesh handle retains the Euler edge count")

    call fortfem_triangle_mesh_edges( &
        handle, 5_c_int, n_edges, edges, triangle_edge_dofs, &
        triangle_edge_signs, status)
    call record_condition(status == 0_c_int, &
        "C triangle mesh handle exports its oriented edge map")
    call record_condition(count(triangle_edge_signs == -1_c_int) > 0 .and. &
        count(triangle_edge_signs == 1_c_int) > 0, &
        "C triangle mesh handle preserves both edge orientations")

    call fortfem_triangle_mesh_free(handle, status)
    call record_condition(status == 0_c_int, &
        "C triangle mesh handle releases its resources")
    call fortfem_triangle_mesh_edges( &
        handle, 5_c_int, n_edges, edges, triangle_edge_dofs, &
        triangle_edge_signs, status)
    call record_condition(status /= 0_c_int, &
        "C mesh API rejects a released handle")

    call check_summary("C triangle mesh handle")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_capi_mesh_handle
