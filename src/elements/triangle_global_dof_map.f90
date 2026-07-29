module fortfem_triangle_global_dof_map
    use fortfem_mesh_2d, only: mesh_2d_t
    implicit none

    private

    public :: build_triangle_trimmed_dof_map
    public :: build_triangle_full_vector_dof_map

contains

    subroutine build_triangle_full_vector_dof_map( &
            mesh, degree, global_dofs, transforms, global_dof_count, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree
        integer, allocatable, intent(out) :: global_dofs(:, :), transforms(:, :)
        integer, intent(out) :: global_dof_count, status

        integer :: cell_dof, cell_dof_count, edge, edge_dof_count
        integer :: edge_dofs(3), edge_moment_count, edge_orientations(3)
        integer :: local_dof, local_dof_count, moment, triangle

        global_dof_count = 0
        status = 1
        if (degree < 1 .or. mesh%n_triangles < 1) return

        if (.not. allocated(mesh%edges)) then
            call mesh%build_edge_connectivity()
        end if
        if (.not. allocated(mesh%edges)) return
        if (.not. allocated(mesh%boundary_edges)) then
            call mesh%find_boundary()
        end if
        if (.not. allocated(mesh%edge_to_dof)) then
            call mesh%build_edge_dof_numbering()
        end if
        if (.not. allocated(mesh%edge_to_dof)) return

        edge_moment_count = degree + 1
        edge_dof_count = mesh%n_edges * edge_moment_count
        cell_dof_count = degree * degree - 1
        local_dof_count = (degree + 1) * (degree + 2)
        global_dof_count = edge_dof_count + &
            mesh%n_triangles * cell_dof_count
        allocate(global_dofs(local_dof_count, mesh%n_triangles))
        allocate(transforms(local_dof_count, mesh%n_triangles))
        global_dofs = 0
        transforms = 1

        do triangle = 1, mesh%n_triangles
            call mesh%get_triangle_edge_dofs( &
                triangle, edge_dofs, edge_orientations)
            do edge = 1, 3
                do moment = 1, edge_moment_count
                    local_dof = (edge - 1) * edge_moment_count + moment
                    global_dofs(local_dof, triangle) = &
                        edge_dofs(edge) * edge_moment_count + moment
                    if (edge_orientations(edge) == -1) then
                        if (mod(moment, 2) == 1) then
                            transforms(local_dof, triangle) = -1
                        end if
                    end if
                end do
            end do
            do cell_dof = 1, cell_dof_count
                local_dof = 3 * edge_moment_count + cell_dof
                global_dofs(local_dof, triangle) = edge_dof_count + &
                    (triangle - 1) * cell_dof_count + cell_dof
            end do
        end do
        status = 0
    end subroutine build_triangle_full_vector_dof_map

    subroutine build_triangle_trimmed_dof_map( &
            mesh, order, global_dofs, transforms, global_dof_count, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: order
        integer, allocatable, intent(out) :: global_dofs(:, :), transforms(:, :)
        integer, intent(out) :: global_dof_count, status

        integer :: cell_dof_count, cell_dof, edge, edge_dof_count
        integer :: edge_dofs(3), edge_orientations(3)
        integer :: local_dof, local_dof_count, moment, triangle

        global_dof_count = 0
        status = 1
        if (order < 1 .or. mesh%n_triangles < 1) return

        if (.not. allocated(mesh%edges)) then
            call mesh%build_edge_connectivity()
        end if
        if (.not. allocated(mesh%edges)) return
        if (.not. allocated(mesh%boundary_edges)) then
            call mesh%find_boundary()
        end if
        if (.not. allocated(mesh%edge_to_dof)) then
            call mesh%build_edge_dof_numbering()
        end if
        if (.not. allocated(mesh%edge_to_dof)) return

        edge_dof_count = mesh%n_edges * order
        cell_dof_count = order * (order - 1)
        local_dof_count = order * (order + 2)
        global_dof_count = edge_dof_count + &
            mesh%n_triangles * cell_dof_count
        allocate(global_dofs(local_dof_count, mesh%n_triangles))
        allocate(transforms(local_dof_count, mesh%n_triangles))
        global_dofs = 0
        transforms = 1

        do triangle = 1, mesh%n_triangles
            call mesh%get_triangle_edge_dofs( &
                triangle, edge_dofs, edge_orientations)
            do edge = 1, 3
                do moment = 1, order
                    local_dof = (edge - 1) * order + moment
                    global_dofs(local_dof, triangle) = &
                        edge_dofs(edge) * order + moment
                    if (edge_orientations(edge) == -1) then
                        if (mod(moment, 2) == 1) then
                            transforms(local_dof, triangle) = -1
                        end if
                    end if
                end do
            end do
            do cell_dof = 1, cell_dof_count
                local_dof = 3 * order + cell_dof
                global_dofs(local_dof, triangle) = edge_dof_count + &
                    (triangle - 1) * cell_dof_count + cell_dof
            end do
        end do
        status = 0
    end subroutine build_triangle_trimmed_dof_map

end module fortfem_triangle_global_dof_map
