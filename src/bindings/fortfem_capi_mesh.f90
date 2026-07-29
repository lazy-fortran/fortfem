module fortfem_capi_mesh
    use iso_c_binding, only: c_double, c_int
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    implicit none
    private

    integer, parameter :: max_meshes = 32
    type(mesh_2d_t), save :: meshes(max_meshes)
    logical, save :: mesh_occupied(max_meshes) = .false.

    public :: fortfem_triangle_edge_map
    public :: fortfem_triangle_mesh_create
    public :: fortfem_triangle_mesh_edges
    public :: fortfem_triangle_mesh_free

contains

    subroutine fortfem_triangle_mesh_create( &
            n_vertices, vertices, n_triangles, triangles, &
            handle, n_edges, status) &
            bind(c, name="fortfem_triangle_mesh_create")
        integer(c_int), value :: n_vertices, n_triangles
        real(c_double), intent(in) :: vertices(*)
        integer(c_int), intent(in) :: triangles(*)
        integer(c_int), intent(out) :: handle, n_edges, status

        integer :: slot

        handle = 0_c_int
        n_edges = 0_c_int
        status = -1_c_int
        if (n_vertices < 3_c_int .or. n_triangles < 1_c_int) return

        slot = find_free_mesh_slot()
        if (slot == 0) then
            status = -3_c_int
            return
        end if
        call initialize_mesh( &
            meshes(slot), n_vertices, vertices, n_triangles, triangles, status)
        if (status /= 0_c_int) return

        mesh_occupied(slot) = .true.
        handle = int(slot, c_int)
        n_edges = int(meshes(slot)%n_edges, c_int)
    end subroutine fortfem_triangle_mesh_create

    subroutine fortfem_triangle_mesh_edges( &
            handle, edge_capacity, n_edges, edges, &
            triangle_edge_dofs, triangle_edge_signs, status) &
            bind(c, name="fortfem_triangle_mesh_edges")
        integer(c_int), value :: handle, edge_capacity
        integer(c_int), intent(out) :: n_edges, edges(*)
        integer(c_int), intent(out) :: triangle_edge_dofs(*)
        integer(c_int), intent(out) :: triangle_edge_signs(*)
        integer(c_int), intent(out) :: status

        integer :: slot

        n_edges = 0_c_int
        status = -3_c_int
        slot = int(handle)
        if (.not. valid_mesh_handle(slot)) return
        call export_edge_map(meshes(slot), edge_capacity, n_edges, edges, &
            triangle_edge_dofs, triangle_edge_signs, status)
    end subroutine fortfem_triangle_mesh_edges

    subroutine fortfem_triangle_mesh_free(handle, status) &
            bind(c, name="fortfem_triangle_mesh_free")
        integer(c_int), value :: handle
        integer(c_int), intent(out) :: status

        integer :: slot

        status = -3_c_int
        slot = int(handle)
        if (.not. valid_mesh_handle(slot)) return

        call meshes(slot)%destroy()
        mesh_occupied(slot) = .false.
        status = 0_c_int
    end subroutine fortfem_triangle_mesh_free

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

        integer(c_int) :: handle, free_status

        if (edge_capacity < 1_c_int) then
            n_edges = 0_c_int
            status = -1_c_int
            return
        end if
        call fortfem_triangle_mesh_create( &
            n_vertices, vertices, n_triangles, triangles, &
            handle, n_edges, status)
        if (status /= 0_c_int) return
        call fortfem_triangle_mesh_edges(handle, edge_capacity, n_edges, &
            edges, triangle_edge_dofs, triangle_edge_signs, status)
        call fortfem_triangle_mesh_free(handle, free_status)
    end subroutine fortfem_triangle_edge_map

    subroutine initialize_mesh( &
            mesh, n_vertices, vertices, n_triangles, triangles, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer(c_int), intent(in) :: n_vertices, n_triangles
        real(c_double), intent(in) :: vertices(*)
        integer(c_int), intent(in) :: triangles(*)
        integer(c_int), intent(out) :: status

        integer :: triangle_id, vertex_id, component, local_vertex

        status = -1_c_int
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
            do local_vertex = 1, 3
                vertex_id = int( &
                    triangles(3 * (triangle_id - 1) + local_vertex))
                if (vertex_id < 0 .or. vertex_id >= mesh%n_vertices) then
                    call mesh%destroy()
                    return
                end if
                mesh%triangles(local_vertex, triangle_id) = vertex_id + 1
            end do
        end do

        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        status = 0_c_int
    end subroutine initialize_mesh

    subroutine export_edge_map( &
            mesh, edge_capacity, n_edges, edges, &
            triangle_edge_dofs, triangle_edge_signs, status)
        type(mesh_2d_t), intent(in) :: mesh
        integer(c_int), intent(in) :: edge_capacity
        integer(c_int), intent(out) :: n_edges, edges(*)
        integer(c_int), intent(out) :: triangle_edge_dofs(*)
        integer(c_int), intent(out) :: triangle_edge_signs(*)
        integer(c_int), intent(out) :: status

        integer :: triangle_id, local_edge, edge_id, dof_id
        integer :: edge_dofs(3), edge_signs(3)

        n_edges = int(mesh%n_edges, c_int)
        if (mesh%n_edges > int(edge_capacity)) then
            status = -2_c_int
            return
        end if

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
    end subroutine export_edge_map

    integer function find_free_mesh_slot() result(slot)
        integer :: candidate

        slot = 0
        do candidate = 1, max_meshes
            if (.not. mesh_occupied(candidate)) then
                slot = candidate
                return
            end if
        end do
    end function find_free_mesh_slot

    logical function valid_mesh_handle(slot) result(valid)
        integer, intent(in) :: slot

        valid = slot >= 1 .and. slot <= max_meshes
        if (valid) valid = mesh_occupied(slot)
    end function valid_mesh_handle

end module fortfem_capi_mesh
