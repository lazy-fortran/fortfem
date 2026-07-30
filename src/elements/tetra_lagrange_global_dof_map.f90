module fortfem_tetra_lagrange_global_dof_map
    use fortfem_tetra_lagrange_arbitrary_order, only: &
        initialize_tetra_lagrange, tetra_lagrange_barycentric_indices, &
        tetra_lagrange_dof_count, tetra_lagrange_t
    implicit none

    private

    public :: build_tetra_lagrange_dof_map

contains

    subroutine build_tetra_lagrange_dof_map( &
            degree, tetrahedra, global_dofs, global_dof_count, status)
        integer, intent(in) :: degree, tetrahedra(:, :)
        integer, allocatable, intent(out) :: global_dofs(:, :)
        integer, intent(out) :: global_dof_count, status

        type(tetra_lagrange_t) :: basis
        integer, allocatable :: indices(:, :)
        integer, allocatable :: key_vertices(:, :), key_weights(:, :)
        integer :: cell, existing, local_dof, local_dof_count
        integer :: local_vertices(4), local_weights(4), vertex
        logical :: found

        global_dof_count = 0
        status = 1
        if (degree < 1 .or. degree > 4) return
        if (size(tetrahedra, 1) /= 4) return
        do cell = 1, size(tetrahedra, 2)
            if (any(tetrahedra(:, cell) < 1)) return
            do vertex = 1, 4
                if (count( &
                    tetrahedra(:, cell) == tetrahedra(vertex, cell)) /= &
                    1) return
            end do
        end do
        call initialize_tetra_lagrange(degree, basis, status)
        if (status /= 0) return
        call tetra_lagrange_barycentric_indices(basis, indices)
        local_dof_count = tetra_lagrange_dof_count(basis)
        allocate(global_dofs(local_dof_count, size(tetrahedra, 2)))
        allocate( &
            key_vertices(4, local_dof_count*size(tetrahedra, 2)), &
            key_weights(4, local_dof_count*size(tetrahedra, 2)))
        key_vertices = 0
        key_weights = 0

        do cell = 1, size(tetrahedra, 2)
            do local_dof = 1, local_dof_count
                local_vertices = 0
                local_weights = 0
                existing = 0
                do vertex = 1, 4
                    if (indices(vertex, local_dof) == 0) cycle
                    existing = existing + 1
                    local_vertices(existing) = tetrahedra(vertex, cell)
                    local_weights(existing) = indices(vertex, local_dof)
                end do
                call sort_key(local_vertices, local_weights)
                found = .false.
                do existing = 1, global_dof_count
                    if (all(key_vertices(:, existing) == local_vertices) .and. &
                        all(key_weights(:, existing) == local_weights)) then
                        found = .true.
                        exit
                    end if
                end do
                if (.not. found) then
                    global_dof_count = global_dof_count + 1
                    existing = global_dof_count
                    key_vertices(:, existing) = local_vertices
                    key_weights(:, existing) = local_weights
                end if
                global_dofs(local_dof, cell) = existing
            end do
        end do
        status = 0
    end subroutine build_tetra_lagrange_dof_map

    pure subroutine sort_key(vertices, weights)
        integer, intent(inout) :: vertices(4), weights(4)

        integer :: first, second, temporary

        do first = 1, 3
            do second = first + 1, 4
                if (vertices(second) == 0) cycle
                if (vertices(first) /= 0 .and. &
                    vertices(first) <= vertices(second)) cycle
                temporary = vertices(first)
                vertices(first) = vertices(second)
                vertices(second) = temporary
                temporary = weights(first)
                weights(first) = weights(second)
                weights(second) = temporary
            end do
        end do
    end subroutine sort_key

end module fortfem_tetra_lagrange_global_dof_map
