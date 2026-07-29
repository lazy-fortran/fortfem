program test_capi_mesh_extension
    use check, only: check_condition, check_summary
    use iso_c_binding, only: c_double, c_int
    implicit none

    real(c_double), parameter :: core_vertices(8) = [ &
        2.0_c_double, 0.0_c_double, 3.0_c_double, 0.0_c_double, &
        3.0_c_double, 1.0_c_double, 2.0_c_double, 1.0_c_double]
    integer(c_int), parameter :: core_triangles(6) = [ &
        0_c_int, 1_c_int, 2_c_int, 0_c_int, 2_c_int, 3_c_int]
    real(c_double), parameter :: outer_vertices(8) = [ &
        1.0_c_double, -1.0_c_double, 4.0_c_double, -1.0_c_double, &
        4.0_c_double, 2.0_c_double, 1.0_c_double, 2.0_c_double]
    integer(c_int) :: edges(34), triangle_edge_dofs(30)
    integer(c_int) :: triangle_edge_signs(30)
    integer(c_int) :: core_handle, extended_handle, first_boundary_dof
    integer(c_int) :: n_edges, n_triangles, n_vertices, status
    integer :: dof
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

        subroutine fortfem_triangle_mesh_extend( &
                core_handle, n_outer_vertices, outer_vertices, &
                handle, n_vertices, n_triangles, n_edges, &
                first_boundary_dof, status) &
                bind(c, name="fortfem_triangle_mesh_extend")
            import c_double, c_int
            integer(c_int), value :: core_handle, n_outer_vertices
            real(c_double), intent(in) :: outer_vertices(*)
            integer(c_int), intent(out) :: handle, n_vertices, n_triangles
            integer(c_int), intent(out) :: n_edges, first_boundary_dof, status
        end subroutine fortfem_triangle_mesh_extend

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
        4_c_int, core_vertices, 2_c_int, core_triangles, &
        core_handle, n_edges, status)
    call fortfem_triangle_mesh_extend( &
        core_handle, 4_c_int, outer_vertices, extended_handle, &
        n_vertices, n_triangles, n_edges, first_boundary_dof, status)
    call record_condition(status == 0_c_int .and. extended_handle > 0_c_int, &
        "C mesh extension returns a retained mesh handle")
    call record_condition(n_vertices == 8_c_int .and. &
        n_triangles == 10_c_int .and. n_edges == 17_c_int, &
        "C mesh extension reports exact nested-square topology")
    call record_condition(first_boundary_dof == 13_c_int, &
        "C mesh extension identifies the trailing outer-boundary DOFs")

    call fortfem_triangle_mesh_edges( &
        extended_handle, 17_c_int, n_edges, edges, &
        triangle_edge_dofs, triangle_edge_signs, status)
    call record_condition(status == 0_c_int, &
        "Extended C mesh exports its oriented edge map")
    do dof = 14, 17
        call record_condition( &
            edges(2 * dof - 1) >= 4_c_int .and. &
            edges(2 * dof) >= 4_c_int, &
            "Trailing boundary DOF joins supplied outer vertices")
    end do

    call fortfem_triangle_mesh_free(extended_handle, status)
    call fortfem_triangle_mesh_free(core_handle, status)
    call check_summary("C triangle mesh extension")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_capi_mesh_extension
