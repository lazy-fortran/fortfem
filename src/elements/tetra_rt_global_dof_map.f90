module fortfem_tetra_rt_global_dof_map
    use fortfem_tetra_face_moment_transforms, only: &
        build_tetra_rt_face_basis_to_local_matrix
    use fortfem_kinds, only: dp
    implicit none

    private

    public :: build_tetra_rt_basis_transform
    public :: build_tetra_rt_dof_map

contains

    subroutine build_tetra_rt_dof_map( &
            degree, tetrahedra, faces, global_dofs, face_orientations, &
            face_permutations, status)
        integer, intent(in) :: degree, tetrahedra(:, :)
        integer, allocatable, intent(out) :: faces(:, :)
        integer, allocatable, intent(out) :: global_dofs(:, :)
        integer, allocatable, intent(out) :: face_orientations(:, :)
        integer, allocatable, intent(out) :: face_permutations(:, :, :)
        integer, intent(out) :: status

        integer, allocatable :: temporary_faces(:, :)
        integer :: canonical_vertices(3), cell_dof_count, dof
        integer :: dof_count, face, face_count, face_dof_count
        integer :: local_face, local_face_vertices(3, 4), local_vertices(3)
        integer :: start, tetrahedron, vertex
        logical :: found

        status = 1
        if (degree < 0) return
        if (size(tetrahedra, 1) /= 4) return
        do tetrahedron = 1, size(tetrahedra, 2)
            if (any(tetrahedra(:, tetrahedron) < 1)) return
            do vertex = 1, 4
                if (count( &
                    tetrahedra(:, tetrahedron) == &
                    tetrahedra(vertex, tetrahedron)) /= 1) return
            end do
        end do

        dof_count = (degree + 1)*(degree + 2)*(degree + 4)/2
        face_dof_count = (degree + 1)*(degree + 2)/2
        cell_dof_count = dof_count - 4*face_dof_count
        allocate(global_dofs(dof_count, size(tetrahedra, 2)))
        allocate(face_orientations(4, size(tetrahedra, 2)))
        allocate(face_permutations(3, 4, size(tetrahedra, 2)))
        allocate(temporary_faces(3, 4*size(tetrahedra, 2)))
        local_face_vertices(:, 1) = [1, 4, 3]
        local_face_vertices(:, 2) = [1, 2, 4]
        local_face_vertices(:, 3) = [1, 3, 2]
        local_face_vertices(:, 4) = [2, 3, 4]

        face_count = 0
        do tetrahedron = 1, size(tetrahedra, 2)
            do local_face = 1, 4
                local_vertices = tetrahedra( &
                    local_face_vertices(:, local_face), tetrahedron)
                canonical_vertices = local_vertices
                call sort_three(canonical_vertices)
                found = .false.
                do face = 1, face_count
                    if (all( &
                        temporary_faces(:, face) == &
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
                face_orientations(local_face, tetrahedron) = &
                    permutation_sign( &
                    face_permutations(:, local_face, tetrahedron))
                start = (local_face - 1)*face_dof_count
                do dof = 1, face_dof_count
                    global_dofs(start + dof, tetrahedron) = &
                        (face - 1)*face_dof_count + dof
                end do
            end do
        end do

        start = 4*face_dof_count
        do tetrahedron = 1, size(tetrahedra, 2)
            do dof = 1, cell_dof_count
                global_dofs(start + dof, tetrahedron) = &
                    face_count*face_dof_count + &
                    (tetrahedron - 1)*cell_dof_count + dof
            end do
        end do
        allocate(faces(3, face_count))
        faces = temporary_faces(:, :face_count)
        status = 0
    end subroutine build_tetra_rt_dof_map

    subroutine build_tetra_rt_basis_transform( &
            degree, face_orientations, face_permutations, transform, status)
        integer, intent(in) :: degree, face_orientations(4)
        integer, intent(in) :: face_permutations(3, 4)
        real(dp), intent(out) :: transform(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: face_transform(:, :)
        integer :: dof, dof_count, face, face_dof_count, start

        transform = 0.0_dp
        status = 1
        if (degree < 0) return
        dof_count = (degree + 1)*(degree + 2)*(degree + 4)/2
        face_dof_count = (degree + 1)*(degree + 2)/2
        if (size(transform, 1) /= dof_count) return
        if (size(transform, 2) /= dof_count) return
        if (any(abs(face_orientations) /= 1)) return

        do dof = 1, dof_count
            transform(dof, dof) = 1.0_dp
        end do
        allocate(face_transform(face_dof_count, face_dof_count))
        do face = 1, 4
            start = (face - 1)*face_dof_count
            call build_tetra_rt_face_basis_to_local_matrix( &
                degree, face_permutations(:, face), face_transform, status)
            if (status /= 0) return
            transform( &
                start + 1:start + face_dof_count, &
                start + 1:start + face_dof_count) = &
                face_orientations(face)*face_transform
        end do
        status = 0
    end subroutine build_tetra_rt_basis_transform

    pure function canonical_position(vertex, canonical_vertices) &
            result(position)
        integer, intent(in) :: vertex, canonical_vertices(3)
        integer :: position

        do position = 1, 3
            if (canonical_vertices(position) == vertex) return
        end do
        position = 0
    end function canonical_position

    pure function permutation_sign(permutation) result(sign_)
        integer, intent(in) :: permutation(3)
        integer :: inversions, first, second, sign_

        inversions = 0
        do first = 1, 2
            do second = first + 1, 3
                if (permutation(first) > permutation(second)) then
                    inversions = inversions + 1
                end if
            end do
        end do
        sign_ = merge(1, -1, modulo(inversions, 2) == 0)
    end function permutation_sign

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

end module fortfem_tetra_rt_global_dof_map
