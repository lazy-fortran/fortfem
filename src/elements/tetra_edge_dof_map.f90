module fortfem_tetra_edge_dof_map
    implicit none

    private

    public :: build_tetra_edge_dof_map

contains

    subroutine build_tetra_edge_dof_map( &
            tetrahedra, edges, global_dofs, orientations, status)
        integer, intent(in) :: tetrahedra(:, :)
        integer, allocatable, intent(out) :: edges(:, :)
        integer, allocatable, intent(out) :: global_dofs(:, :)
        integer, allocatable, intent(out) :: orientations(:, :)
        integer, intent(out) :: status

        integer, allocatable :: temporary_edges(:, :)
        integer :: edge, edge_count, first, first_vertex, vertex
        integer :: local_edge, local_vertices(2, 6), second, second_vertex
        integer :: tetrahedron, other_vertex
        logical :: found

        status = 1
        if (size(tetrahedra, 1) /= 4) return
        if (size(tetrahedra, 2) < 1) return
        if (any(tetrahedra < 1)) return
        local_vertices(:, 1) = [1, 2]
        local_vertices(:, 2) = [1, 3]
        local_vertices(:, 3) = [1, 4]
        local_vertices(:, 4) = [2, 3]
        local_vertices(:, 5) = [2, 4]
        local_vertices(:, 6) = [3, 4]
        allocate(global_dofs(6, size(tetrahedra, 2)))
        allocate(orientations(6, size(tetrahedra, 2)))
        allocate(temporary_edges(2, 6 * size(tetrahedra, 2)))

        edge_count = 0
        do tetrahedron = 1, size(tetrahedra, 2)
            do vertex = 1, 3
                do other_vertex = vertex + 1, 4
                    if (tetrahedra(vertex, tetrahedron) == &
                        tetrahedra(other_vertex, tetrahedron)) return
                end do
            end do
            do local_edge = 1, 6
                first_vertex = tetrahedra( &
                    local_vertices(1, local_edge), tetrahedron)
                second_vertex = tetrahedra( &
                    local_vertices(2, local_edge), tetrahedron)
                first = min(first_vertex, second_vertex)
                second = max(first_vertex, second_vertex)
                found = .false.
                do edge = 1, edge_count
                    if (temporary_edges(1, edge) == first .and. &
                        temporary_edges(2, edge) == second) then
                        found = .true.
                        exit
                    end if
                end do
                if (.not. found) then
                    edge_count = edge_count + 1
                    edge = edge_count
                    temporary_edges(:, edge) = [first, second]
                end if
                global_dofs(local_edge, tetrahedron) = edge
                if (first_vertex < second_vertex) then
                    orientations(local_edge, tetrahedron) = 1
                else
                    orientations(local_edge, tetrahedron) = -1
                end if
            end do
        end do
        allocate(edges(2, edge_count))
        edges = temporary_edges(:, :edge_count)
        status = 0
    end subroutine build_tetra_edge_dof_map

end module fortfem_tetra_edge_dof_map
