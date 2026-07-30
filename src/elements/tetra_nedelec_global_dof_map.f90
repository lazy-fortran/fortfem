module fortfem_tetra_nedelec_global_dof_map
    use fortfem_tetra_face_moment_transforms, only: &
        build_tetra_face_basis_to_local_matrix
    use fortfem_kinds, only: dp
    use fortfem_tetra_edge_dof_map, only: build_tetra_edge_dof_map
    implicit none

    private

    public :: build_tetra_nedelec_basis_transform
    public :: build_tetra_nedelec_dof_map

contains

    subroutine build_tetra_nedelec_dof_map( &
            order, tetrahedra, edges, faces, global_dofs, edge_orientations, &
            face_permutations, status)
        integer, intent(in) :: order, tetrahedra(:, :)
        integer, allocatable, intent(out) :: edges(:, :), faces(:, :)
        integer, allocatable, intent(out) :: global_dofs(:, :)
        integer, allocatable, intent(out) :: edge_orientations(:, :)
        integer, allocatable, intent(out) :: face_permutations(:, :, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_indices(:, :), temporary_faces(:, :)
        integer :: canonical_vertices(3), cell_dof_count, dof
        integer :: dof_count, edge, edge_dof_count, face, face_count
        integer :: face_dof_count, face_start, local_face
        integer :: local_face_vertices(3, 4), local_vertices(3)
        integer :: tetrahedron, vertex
        logical :: found

        status = 1
        if (order < 1) return
        call build_tetra_edge_dof_map( &
            tetrahedra, edges, edge_indices, edge_orientations, status)
        if (status /= 0) return

        dof_count = order * (order + 2) * (order + 3) / 2
        face_dof_count = order * (order - 1)
        cell_dof_count = order * (order - 1) * (order - 2) / 2
        edge_dof_count = order * size(edges, 2)
        allocate(global_dofs(dof_count, size(tetrahedra, 2)))
        allocate(face_permutations(3, 4, size(tetrahedra, 2)))
        allocate(temporary_faces(3, 4 * size(tetrahedra, 2)))
        local_face_vertices(:, 1) = [1, 2, 3]
        local_face_vertices(:, 2) = [1, 2, 4]
        local_face_vertices(:, 3) = [1, 3, 4]
        local_face_vertices(:, 4) = [2, 3, 4]

        do tetrahedron = 1, size(tetrahedra, 2)
            do edge = 1, 6
                do dof = 1, order
                    global_dofs((edge - 1) * order + dof, tetrahedron) = &
                        (edge_indices(edge, tetrahedron) - 1) * order + dof
                end do
            end do
        end do

        face_count = 0
        do tetrahedron = 1, size(tetrahedra, 2)
            do local_face = 1, 4
                local_vertices = tetrahedra( &
                    local_face_vertices(:, local_face), tetrahedron)
                canonical_vertices = local_vertices
                call sort_three(canonical_vertices)
                found = .false.
                do face = 1, face_count
                    if (all(temporary_faces(:, face) == &
                        canonical_vertices)) then
                        found = .true.
                        exit
                    end if
                end do
                if (.not. found) then
                    face_count = face_count + 1
                    face = face_count
                    temporary_faces(:, face) = canonical_vertices
                end if
                do vertex = 1, 3
                    face_permutations(vertex, local_face, tetrahedron) = &
                        canonical_position( &
                        local_vertices(vertex), canonical_vertices)
                end do
                face_start = 6 * order + &
                    (local_face - 1) * face_dof_count
                do dof = 1, face_dof_count
                    global_dofs(face_start + dof, tetrahedron) = &
                        edge_dof_count + (face - 1) * face_dof_count + dof
                end do
            end do
        end do

        face_start = 6 * order + 4 * face_dof_count
        do tetrahedron = 1, size(tetrahedra, 2)
            do dof = 1, cell_dof_count
                global_dofs(face_start + dof, tetrahedron) = &
                    edge_dof_count + face_count * face_dof_count + &
                    (tetrahedron - 1) * cell_dof_count + dof
            end do
        end do
        allocate(faces(3, face_count))
        faces = temporary_faces(:, :face_count)
        status = 0
    end subroutine build_tetra_nedelec_dof_map

    subroutine build_tetra_nedelec_basis_transform( &
            order, edge_orientations, face_permutations, transform, status)
        integer, intent(in) :: order, edge_orientations(6)
        integer, intent(in) :: face_permutations(3, 4)
        real(dp), intent(out) :: transform(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: face_transform(:, :)
        integer :: dof, dof_count, edge, face, face_dof_count, start

        transform = 0.0_dp
        status = 1
        dof_count = order * (order + 2) * (order + 3) / 2
        if (order < 1) return
        if (size(transform, 1) /= dof_count .or. &
            size(transform, 2) /= dof_count) return
        if (any(abs(edge_orientations) /= 1)) return

        do dof = 1, dof_count
            transform(dof, dof) = 1.0_dp
        end do
        do edge = 1, 6
            if (edge_orientations(edge) == 1) cycle
            do dof = 1, order, 2
                start = (edge - 1) * order + dof
                transform(start, start) = -1.0_dp
            end do
        end do

        face_dof_count = order * (order - 1)
        if (face_dof_count == 0) then
            status = 0
            return
        end if
        allocate(face_transform(face_dof_count, face_dof_count))
        do face = 1, 4
            start = 6 * order + (face - 1) * face_dof_count
            transform( &
                start + 1:start + face_dof_count, &
                start + 1:start + face_dof_count) = 0.0_dp
            call build_tetra_face_basis_to_local_matrix( &
                order, face_permutations(:, face), face_transform, status)
            if (status /= 0) return
            transform( &
                start + 1:start + face_dof_count, &
                start + 1:start + face_dof_count) = face_transform
        end do
        status = 0
    end subroutine build_tetra_nedelec_basis_transform

    pure function canonical_position(vertex, canonical_vertices) result(position)
        integer, intent(in) :: vertex, canonical_vertices(3)
        integer :: position

        do position = 1, 3
            if (canonical_vertices(position) == vertex) return
        end do
        position = 0
    end function canonical_position

    pure subroutine sort_three(values)
        integer, intent(inout) :: values(3)

        integer :: first, second, temporary

        do first = 1, 2
            do second = first + 1, 3
                if (values(first) <= values(second)) cycle
                temporary = values(first)
                values(first) = values(second)
                values(second) = temporary
            end do
        end do
    end subroutine sort_three

end module fortfem_tetra_nedelec_global_dof_map
