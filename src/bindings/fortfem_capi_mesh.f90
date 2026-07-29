module fortfem_capi_mesh
    use iso_c_binding, only: c_double, c_int
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    implicit none
    private

    public :: fortfem_triangle_edge_map

contains

    subroutine fortfem_triangle_edge_map( &
            n_vertices, vertices, n_triangles, triangles, &
            edge_capacity, n_edges, edges, triangle_edge_dofs, &
            triangle_edge_signs, status) &
            bind(c, name="fortfem_triangle_edge_map")
        integer(c_int), value :: n_vertices, n_triangles, edge_capacity
        real(c_double), intent(in) :: vertices(*)
        integer(c_int), intent(in) :: triangles(*)
        integer(c_int), intent(out) :: n_edges, edges(*)
        integer(c_int), intent(out) :: triangle_edge_dofs(*)
        integer(c_int), intent(out) :: triangle_edge_signs(*)
        integer(c_int), intent(out) :: status

        type(mesh_2d_t) :: mesh
        integer :: triangle_id, vertex_id, component, local_edge
        integer :: edge_id, dof_id, edge_dofs(3), edge_signs(3)

        n_edges = 0_c_int
        status = -1_c_int
        if (n_vertices < 3_c_int .or. n_triangles < 1_c_int .or. &
            edge_capacity < 1_c_int) return

        mesh%n_vertices = int(n_vertices)
        mesh%n_triangles = int(n_triangles)
        mesh%has_triangles = .true.
        allocate(mesh%vertices(2, mesh%n_vertices))
        allocate(mesh%triangles(3, mesh%n_triangles))

        do vertex_id = 1, mesh%n_vertices
            do component = 1, 2
                mesh%vertices(component, vertex_id) = real( &
                    vertices(2 * (vertex_id - 1) + component), dp)
            end do
        end do

        do triangle_id = 1, mesh%n_triangles
            do local_edge = 1, 3
                vertex_id = int(triangles(3 * (triangle_id - 1) + local_edge))
                if (vertex_id < 0 .or. vertex_id >= mesh%n_vertices) then
                    call mesh%destroy()
                    return
                end if
                mesh%triangles(local_edge, triangle_id) = vertex_id + 1
            end do
        end do

        call mesh%build_edge_connectivity()
        if (mesh%n_edges > int(edge_capacity)) then
            n_edges = int(mesh%n_edges, c_int)
            status = -2_c_int
            call mesh%destroy()
            return
        end if
        call mesh%build_edge_dof_numbering()
        n_edges = int(mesh%n_edges, c_int)

        do dof_id = 0, mesh%n_edges - 1
            edge_id = mesh%dof_to_edge(dof_id + 1)
            edges(2 * dof_id + 1) = int( &
                mesh%edges(1, edge_id) - 1, c_int)
            edges(2 * dof_id + 2) = int( &
                mesh%edges(2, edge_id) - 1, c_int)
        end do

        do triangle_id = 1, mesh%n_triangles
            call mesh%get_triangle_edge_dofs( &
                triangle_id, edge_dofs, edge_signs)
            do local_edge = 1, 3
                triangle_edge_dofs(3 * (triangle_id - 1) + local_edge) = &
                    int(edge_dofs(local_edge), c_int)
                triangle_edge_signs(3 * (triangle_id - 1) + local_edge) = &
                    int(edge_signs(local_edge), c_int)
            end do
        end do

        status = 0_c_int
        call mesh%destroy()
    end subroutine fortfem_triangle_edge_map

end module fortfem_capi_mesh
