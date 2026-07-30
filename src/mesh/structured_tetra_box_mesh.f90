module fortfem_structured_tetra_box_mesh
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det3
    implicit none

    private

    public :: generate_structured_tetra_box_mesh

contains

    subroutine generate_structured_tetra_box_mesh( &
            bounds, counts, vertices, tetrahedra, status)
        real(dp), intent(in) :: bounds(3, 2)
        integer, intent(in) :: counts(3)
        real(dp), allocatable, intent(out) :: vertices(:, :)
        integer, allocatable, intent(out) :: tetrahedra(:, :)
        integer, intent(out) :: status

        integer :: brick(8), cell, i, j, k, temporary, vertex
        real(dp) :: coordinate(3), jacobian(3, 3)

        status = 1
        if (any(counts < 1)) return
        if (any(bounds(:, 2) <= bounds(:, 1))) return
        allocate(vertices(3, product(counts + 1)))
        do k = 0, counts(3)
            do j = 0, counts(2)
                do i = 0, counts(1)
                    vertex = vertex_index(i, j, k, counts)
                    coordinate = real([i, j, k], dp)/real(counts, dp)
                    vertices(:, vertex) = bounds(:, 1) + &
                        coordinate*(bounds(:, 2) - bounds(:, 1))
                end do
            end do
        end do

        allocate(tetrahedra(4, 6*product(counts)))
        cell = 0
        do k = 0, counts(3) - 1
            do j = 0, counts(2) - 1
                do i = 0, counts(1) - 1
                    brick = [ &
                        vertex_index(i, j, k, counts), &
                        vertex_index(i + 1, j, k, counts), &
                        vertex_index(i, j + 1, k, counts), &
                        vertex_index(i + 1, j + 1, k, counts), &
                        vertex_index(i, j, k + 1, counts), &
                        vertex_index(i + 1, j, k + 1, counts), &
                        vertex_index(i, j + 1, k + 1, counts), &
                        vertex_index(i + 1, j + 1, k + 1, counts)]
                    call add_brick_tetrahedra(brick, tetrahedra, cell)
                end do
            end do
        end do
        do cell = 1, size(tetrahedra, 2)
            jacobian(:, 1) = vertices(:, tetrahedra(2, cell)) - &
                vertices(:, tetrahedra(1, cell))
            jacobian(:, 2) = vertices(:, tetrahedra(3, cell)) - &
                vertices(:, tetrahedra(1, cell))
            jacobian(:, 3) = vertices(:, tetrahedra(4, cell)) - &
                vertices(:, tetrahedra(1, cell))
            if (det3(jacobian) < 0.0_dp) then
                temporary = tetrahedra(3, cell)
                tetrahedra(3, cell) = tetrahedra(4, cell)
                tetrahedra(4, cell) = temporary
            end if
        end do
        status = 0
    end subroutine generate_structured_tetra_box_mesh

    pure integer function vertex_index(i, j, k, counts) result(index)
        integer, intent(in) :: i, j, k, counts(3)

        index = 1 + i + (counts(1) + 1)*(j + (counts(2) + 1)*k)
    end function vertex_index

    pure subroutine add_brick_tetrahedra(brick, tetrahedra, cell_count)
        integer, intent(in) :: brick(8)
        integer, intent(inout) :: tetrahedra(:, :), cell_count

        integer :: local_cells(4, 6), local_cell

        local_cells(:, 1) = brick([1, 2, 4, 8])
        local_cells(:, 2) = brick([1, 2, 6, 8])
        local_cells(:, 3) = brick([1, 3, 4, 8])
        local_cells(:, 4) = brick([1, 3, 7, 8])
        local_cells(:, 5) = brick([1, 5, 6, 8])
        local_cells(:, 6) = brick([1, 5, 7, 8])
        do local_cell = 1, 6
            cell_count = cell_count + 1
            tetrahedra(:, cell_count) = local_cells(:, local_cell)
        end do
    end subroutine add_brick_tetrahedra

end module fortfem_structured_tetra_box_mesh
