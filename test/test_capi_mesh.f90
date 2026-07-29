program test_capi_mesh
    use iso_c_binding, only: c_double, c_int
    use check, only: check_condition, check_summary
    implicit none

    integer(c_int), parameter :: n_vertices = 4_c_int
    integer(c_int), parameter :: n_triangles = 2_c_int
    integer(c_int), parameter :: edge_capacity = 5_c_int
    real(c_double), parameter :: vertices(8) = [ &
        0.0_c_double, 0.0_c_double, 1.0_c_double, 0.0_c_double, &
        1.0_c_double, 1.0_c_double, 0.0_c_double, 1.0_c_double]
    integer(c_int), parameter :: triangles(6) = [ &
        0_c_int, 1_c_int, 2_c_int, 0_c_int, 2_c_int, 3_c_int]
    integer(c_int) :: edges(2 * edge_capacity)
    integer(c_int) :: triangle_edge_dofs(3 * n_triangles)
    integer(c_int) :: triangle_edge_signs(3 * n_triangles)
    integer(c_int) :: n_edges, status
    integer :: triangle_id, local_edge, offset, dof
    integer(c_int) :: local_start, local_end, edge_start, edge_end
    logical :: all_passed

    interface
        subroutine fortfem_triangle_edge_map( &
                n_vertices, vertices, n_triangles, triangles, &
                edge_capacity, n_edges, edges, triangle_edge_dofs, &
                triangle_edge_signs, status) &
                bind(c, name="fortfem_triangle_edge_map")
            import c_double, c_int
            integer(c_int), value :: n_vertices, n_triangles, edge_capacity
            real(c_double), intent(in) :: vertices(*)
            integer(c_int), intent(in) :: triangles(*)
            integer(c_int), intent(out) :: n_edges, edges(*)
            integer(c_int), intent(out) :: triangle_edge_dofs(*)
            integer(c_int), intent(out) :: triangle_edge_signs(*)
            integer(c_int), intent(out) :: status
        end subroutine fortfem_triangle_edge_map
    end interface

    all_passed = .true.
    call fortfem_triangle_edge_map( &
        n_vertices, vertices, n_triangles, triangles, edge_capacity, &
        n_edges, edges, triangle_edge_dofs, triangle_edge_signs, status)

    call record_condition(status == 0_c_int, &
        "C mesh topology call reports success")
    call record_condition(n_edges == 5_c_int, &
        "C mesh topology returns the Euler edge count")

    do triangle_id = 0, int(n_triangles) - 1
        do local_edge = 0, 2
            offset = 3 * triangle_id + local_edge
            local_start = triangles(offset + 1)
            local_end = triangles(3 * triangle_id + mod(local_edge + 1, 3) + 1)
            dof = int(triangle_edge_dofs(offset + 1))
            edge_start = edges(2 * dof + 1)
            edge_end = edges(2 * dof + 2)

            call record_condition(dof >= 0 .and. dof < int(n_edges), &
                "C triangle edge DOF is in range")
            call record_condition( &
                min(local_start, local_end) == edge_start .and. &
                max(local_start, local_end) == edge_end, &
                "C triangle edge DOF references the matching global edge")
            if (local_start == edge_start) then
                call record_condition(triangle_edge_signs(offset + 1) == 1, &
                    "C edge sign preserves the local direction")
            else
                call record_condition(triangle_edge_signs(offset + 1) == -1, &
                    "C edge sign reverses the local direction")
            end if
        end do
    end do

    call check_summary("C mesh topology API")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_capi_mesh
