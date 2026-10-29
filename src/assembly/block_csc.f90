module fortfem_block_csc
    !! Dense local blocks of several groups to one square CSC matrix.
    !!
    !! A group is a connectivity conn(nloc, nblocks) of global dofs; block e couples
    !! conn(:, e) with itself, and an entry with dof <= 0 is dropped. Groups may
    !! differ in nloc (vector P1 triangles, boundary edges, a dense interface block).
    !! The pattern holds, per group, the CSC value position of every block entry and
    !! the inverse map from a value position to its block entries in ascending
    !! order. Adding is then a gather-sum per position, parallel over positions and
    !! bitwise independent of the thread count; gather is parallel over blocks.
    !!
    !! Assembly is linear in the blocks: the JVP of block_csc_add is block_csc_add of
    !! the block tangents and its VJP is block_csc_gather (real 0/1 coefficients, so
    !! the same holds for complex values with conjugate-linear cotangents).
    use, intrinsic :: iso_fortran_env, only: int64
    use fortfem_element_to_csr, only: csr_sort_row
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_t, csc_z_t
    implicit none
    private

    public :: block_csc_pattern_t
    public :: block_csc_init
    public :: block_csc_add_group
    public :: block_csc_finalize
    public :: block_csc_new_matrix
    public :: block_csc_add
    public :: block_csc_gather

    type :: block_csc_group_t
        integer :: nloc = 0
        integer :: nblocks = 0
        integer, allocatable :: conn(:, :) ! kept until finalize
        integer, allocatable :: position(:, :, :) ! (nloc, nloc, nblocks); 0 dropped
        integer, allocatable :: source_ptr(:) ! (nnz + 1)
        integer, allocatable :: source(:) ! a + nloc*(b - 1 + nloc*(e - 1))
    end type block_csc_group_t

    type :: block_csc_pattern_t
        integer :: n = 0
        integer :: nnz = 0
        logical :: finalized = .false.
        integer, allocatable :: col_ptr(:) ! (n + 1)
        integer, allocatable :: row_idx(:) ! (nnz), sorted within each column
        type(block_csc_group_t), allocatable :: groups(:)
    end type block_csc_pattern_t

    interface block_csc_new_matrix
        module procedure block_csc_new_matrix_real, block_csc_new_matrix_complex
    end interface block_csc_new_matrix

    interface block_csc_add
        module procedure block_csc_add_real, block_csc_add_complex
    end interface block_csc_add

    interface block_csc_gather
        module procedure block_csc_gather_real, block_csc_gather_complex
    end interface block_csc_gather

contains

    subroutine block_csc_init(pattern, n)
        !! Empty pattern of an n x n matrix; add groups, then finalize.
        type(block_csc_pattern_t), intent(out) :: pattern
        integer, intent(in) :: n

        if (n < 1) error stop "block_csc_init: n must be positive"
        pattern%n = n
        allocate(pattern%groups(0))
    end subroutine block_csc_init

    subroutine block_csc_add_group(pattern, conn, group)
        !! Register conn(nloc, nblocks); group is its index for add and gather.
        type(block_csc_pattern_t), intent(inout) :: pattern
        integer, intent(in) :: conn(:, :)
        integer, intent(out) :: group

        type(block_csc_group_t), allocatable :: groups(:)

        if (pattern%finalized) error stop "block_csc_add_group: pattern is final"
        if (.not. allocated(pattern%groups)) error stop "block_csc_add_group: no init"
        if (any(conn > pattern%n)) error stop "block_csc_add_group: dof exceeds n"
        if (int(size(conn, 1), int64)**2*size(conn, 2) > int(huge(1), int64)) &
            error stop "block_csc_add_group: group exceeds default integer range"
        group = size(pattern%groups) + 1
        allocate(groups(group))
        groups(1:group - 1) = pattern%groups
        groups(group)%nloc = size(conn, 1)
        groups(group)%nblocks = size(conn, 2)
        groups(group)%conn = conn
        call move_alloc(groups, pattern%groups)
    end subroutine block_csc_add_group

    subroutine block_csc_finalize(pattern)
        !! Column structure, per-group positions and inverse maps.
        type(block_csc_pattern_t), intent(inout) :: pattern

        integer, allocatable :: dof_ptr(:), dof_group(:), dof_block(:), marker(:)
        integer, allocatable :: slot(:)
        integer :: column, count_column, e, g, k, n, ngroups, nloc, p, row, start
        integer(int64) :: total

        if (pattern%finalized) return
        if (.not. allocated(pattern%groups)) error stop "block_csc_finalize: no init"
        n = pattern%n
        ngroups = size(pattern%groups)
        call dof_to_block(pattern, n, dof_ptr, dof_group, dof_block)

        ! Column lengths: the union of the dofs of every block holding the column.
        allocate(pattern%col_ptr(n + 1))
        !$omp parallel default(none) &
        !$omp shared(pattern, n, dof_ptr, dof_group, dof_block) &
        !$omp private(marker, column, count_column, k, g, e, nloc, p, row)
        allocate(marker(n))
        marker = 0
        !$omp do schedule(static)
        do column = 1, n
            count_column = 0
            do k = dof_ptr(column), dof_ptr(column + 1) - 1
                g = dof_group(k)
                e = dof_block(k)
                nloc = pattern%groups(g)%nloc
                do p = 1, nloc
                    row = pattern%groups(g)%conn(p, e)
                    if (row <= 0) cycle
                    if (marker(row) == column) cycle
                    marker(row) = column
                    count_column = count_column + 1
                end do
            end do
            pattern%col_ptr(column + 1) = count_column
        end do
        !$omp end do
        deallocate(marker)
        !$omp end parallel
        total = 1
        pattern%col_ptr(1) = 1
        do column = 1, n
            total = total + pattern%col_ptr(column + 1)
            if (total > int(huge(1), int64)) &
                error stop "block_csc_finalize: nnz exceeds default integer range"
            pattern%col_ptr(column + 1) = int(total)
        end do
        pattern%nnz = pattern%col_ptr(n + 1) - 1
        allocate(pattern%row_idx(pattern%nnz))
        do g = 1, ngroups
            associate (group => pattern%groups(g))
                allocate(group%position(group%nloc, group%nloc, group%nblocks))
                allocate(group%source_ptr(pattern%nnz + 1))
                group%position = 0
            end associate
        end do

        ! Rows sorted per column, then positions and per-position entry counts.
        ! A position (row, column) only receives entries whose column dof is
        ! `column`, so each column owner writes disjoint data.
        !$omp parallel default(none) &
        !$omp shared(pattern, n, dof_ptr, dof_group, dof_block) &
        !$omp private(marker, slot, column, start, p, k, g, e)
        allocate(marker(n), slot(n))
        marker = 0
        !$omp do schedule(static)
        do column = 1, n
            start = pattern%col_ptr(column)
            p = start - 1
            do k = dof_ptr(column), dof_ptr(column + 1) - 1
                g = dof_group(k)
                e = dof_block(k)
                call collect_rows(pattern%groups(g)%conn(:, e), column, marker, &
                    pattern%row_idx, p)
            end do
            call csr_sort_row(pattern%row_idx, start, p)
            do p = start, pattern%col_ptr(column + 1) - 1
                slot(pattern%row_idx(p)) = p
                do g = 1, size(pattern%groups)
                    pattern%groups(g)%source_ptr(p + 1) = 0
                end do
            end do
            do k = dof_ptr(column), dof_ptr(column + 1) - 1
                call set_positions(pattern%groups(dof_group(k)), dof_block(k), &
                    column, slot)
            end do
        end do
        !$omp end do
        deallocate(marker, slot)
        !$omp end parallel

        do g = 1, ngroups
            call build_sources(pattern%groups(g), pattern%nnz, &
                pattern%col_ptr, dof_ptr, dof_group, dof_block, g)
            deallocate(pattern%groups(g)%conn)
        end do
        pattern%finalized = .true.
    end subroutine block_csc_finalize

    subroutine dof_to_block(pattern, n, dof_ptr, dof_group, dof_block)
        !! Counting sort of (group, block) incidences per dof, ascending in group
        !! and block. A dof repeated within a block is listed once for it.
        type(block_csc_pattern_t), intent(in) :: pattern
        integer, intent(in) :: n
        integer, allocatable, intent(out) :: dof_ptr(:), dof_group(:), dof_block(:)

        integer, allocatable :: fill(:), stamp(:)
        integer :: a, dof, e, g, label

        allocate(dof_ptr(n + 1), stamp(n), fill(n))
        dof_ptr = 0
        stamp = 0
        label = 0
        do g = 1, size(pattern%groups)
            do e = 1, pattern%groups(g)%nblocks
                label = label + 1
                do a = 1, pattern%groups(g)%nloc
                    dof = pattern%groups(g)%conn(a, e)
                    if (dof <= 0) cycle
                    if (stamp(dof) == label) cycle
                    stamp(dof) = label
                    dof_ptr(dof + 1) = dof_ptr(dof + 1) + 1
                end do
            end do
        end do
        dof_ptr(1) = 1
        do dof = 1, n
            dof_ptr(dof + 1) = dof_ptr(dof + 1) + dof_ptr(dof)
        end do
        allocate(dof_group(dof_ptr(n + 1) - 1), dof_block(dof_ptr(n + 1) - 1))
        fill = dof_ptr(1:n)
        stamp = 0
        label = 0
        do g = 1, size(pattern%groups)
            do e = 1, pattern%groups(g)%nblocks
                label = label + 1
                do a = 1, pattern%groups(g)%nloc
                    dof = pattern%groups(g)%conn(a, e)
                    if (dof <= 0) cycle
                    if (stamp(dof) == label) cycle
                    stamp(dof) = label
                    dof_group(fill(dof)) = g
                    dof_block(fill(dof)) = e
                    fill(dof) = fill(dof) + 1
                end do
            end do
        end do
    end subroutine dof_to_block

    pure subroutine collect_rows(dofs, column, marker, row_idx, last)
        integer, intent(in) :: dofs(:), column
        integer, intent(inout) :: marker(:), row_idx(:), last
        integer :: a

        do a = 1, size(dofs)
            if (dofs(a) <= 0) cycle
            if (marker(dofs(a)) == column) cycle
            marker(dofs(a)) = column
            last = last + 1
            row_idx(last) = dofs(a)
        end do
    end subroutine collect_rows

    subroutine set_positions(group, e, column, slot)
        !! Positions of column `column` of block e (slot maps row -> position)
        !! and the entry count of each position.
        type(block_csc_group_t), intent(inout) :: group
        integer, intent(in) :: e, column, slot(:)
        integer :: a, b, p

        do b = 1, group%nloc
            if (group%conn(b, e) /= column) cycle
            do a = 1, group%nloc
                if (group%conn(a, e) <= 0) cycle
                p = slot(group%conn(a, e))
                group%position(a, b, e) = p
                group%source_ptr(p + 1) = group%source_ptr(p + 1) + 1
            end do
        end do
    end subroutine set_positions

    subroutine build_sources(group, nnz, col_ptr, dof_ptr, dof_group, dof_block, g)
        !! Inverse map: entries of each position in ascending (block, column, row)
        !! order, which is the fixed summation order of block_csc_add.
        type(block_csc_group_t), intent(inout) :: group
        integer, intent(in) :: nnz, col_ptr(:), dof_ptr(:), dof_group(:)
        integer, intent(in) :: dof_block(:), g

        integer, allocatable :: fill(:)
        integer :: a, b, column, e, k, nloc, p

        group%source_ptr(1) = 1
        do p = 1, nnz
            group%source_ptr(p + 1) = group%source_ptr(p + 1) + group%source_ptr(p)
        end do
        allocate(group%source(group%source_ptr(nnz + 1) - 1), fill(nnz))
        nloc = group%nloc
        !$omp parallel do default(none) schedule(static) &
        !$omp shared(group, col_ptr, dof_ptr, dof_group, dof_block, g, fill, nloc) &
        !$omp private(column, k, e, b, a, p)
        do column = 1, size(col_ptr) - 1
            do p = col_ptr(column), col_ptr(column + 1) - 1
                fill(p) = group%source_ptr(p)
            end do
            do k = dof_ptr(column), dof_ptr(column + 1) - 1
                if (dof_group(k) /= g) cycle
                e = dof_block(k)
                do b = 1, nloc
                    if (group%conn(b, e) /= column) cycle
                    do a = 1, nloc
                        p = group%position(a, b, e)
                        if (p == 0) cycle
                        group%source(fill(p)) = a + nloc*(b - 1 + nloc*(e - 1))
                        fill(p) = fill(p) + 1
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine build_sources

    subroutine block_csc_new_matrix_real(pattern, matrix)
        !! Real CSC matrix with the pattern and zero values.
        type(block_csc_pattern_t), intent(in) :: pattern
        type(csc_t), intent(out) :: matrix

        if (.not. pattern%finalized) error stop "block_csc_new_matrix: not final"
        matrix%nrow = pattern%n
        matrix%ncol = pattern%n
        matrix%nnz = pattern%nnz
        matrix%col_ptr = pattern%col_ptr
        matrix%row_idx = pattern%row_idx
        allocate(matrix%val(pattern%nnz))
        call zero_real(pattern, matrix%val)
    end subroutine block_csc_new_matrix_real

    subroutine block_csc_new_matrix_complex(pattern, matrix)
        !! Complex CSC matrix with the pattern and zero values.
        type(block_csc_pattern_t), intent(in) :: pattern
        type(csc_z_t), intent(out) :: matrix

        if (.not. pattern%finalized) error stop "block_csc_new_matrix: not final"
        matrix%nrow = pattern%n
        matrix%ncol = pattern%n
        matrix%nnz = pattern%nnz
        matrix%col_ptr = pattern%col_ptr
        matrix%row_idx = pattern%row_idx
        allocate(matrix%val(pattern%nnz))
        call zero_complex(pattern, matrix%val)
    end subroutine block_csc_new_matrix_complex

    subroutine zero_real(pattern, values)
        type(block_csc_pattern_t), intent(in) :: pattern
        real(dp), intent(out) :: values(:)
        integer :: column

        ! First touch by the column owner of block_csc_add.
        !$omp parallel do default(none) shared(pattern, values) private(column) &
        !$omp schedule(static)
        do column = 1, pattern%n
            values(pattern%col_ptr(column):pattern%col_ptr(column + 1) - 1) = 0.0_dp
        end do
        !$omp end parallel do
    end subroutine zero_real

    subroutine zero_complex(pattern, values)
        type(block_csc_pattern_t), intent(in) :: pattern
        complex(dp), intent(out) :: values(:)
        integer :: column

        !$omp parallel do default(none) shared(pattern, values) private(column) &
        !$omp schedule(static)
        do column = 1, pattern%n
            values(pattern%col_ptr(column):pattern%col_ptr(column + 1) - 1) = &
                (0.0_dp, 0.0_dp)
        end do
        !$omp end parallel do
    end subroutine zero_complex

    subroutine block_csc_add_real(pattern, group, blocks, values)
        !! values += sum of blocks(a, b, e) at (conn(a, e), conn(b, e)).
        type(block_csc_pattern_t), intent(in) :: pattern
        integer, intent(in) :: group
        real(dp), intent(in), contiguous :: blocks(:, :, :)
        real(dp), intent(inout) :: values(:)

        call check_shapes(pattern, group, shape(blocks), size(values), "add")
        call add_real(pattern%groups(group), pattern%nnz, blocks, values)
    end subroutine block_csc_add_real

    subroutine add_real(group, nnz, blocks, values)
        type(block_csc_group_t), intent(in) :: group
        integer, intent(in) :: nnz
        real(dp), intent(in) :: blocks(*)
        real(dp), intent(inout) :: values(:)
        real(dp) :: total
        integer :: k, p

        !$omp parallel do default(none) shared(group, nnz, blocks, values) &
        !$omp private(p, k, total) schedule(static)
        do p = 1, nnz
            if (group%source_ptr(p) == group%source_ptr(p + 1)) cycle
            total = 0.0_dp
            do k = group%source_ptr(p), group%source_ptr(p + 1) - 1
                total = total + blocks(group%source(k))
            end do
            values(p) = values(p) + total
        end do
        !$omp end parallel do
    end subroutine add_real

    subroutine block_csc_add_complex(pattern, group, blocks, values)
        !! values += sum of blocks(a, b, e) at (conn(a, e), conn(b, e)).
        type(block_csc_pattern_t), intent(in) :: pattern
        integer, intent(in) :: group
        complex(dp), intent(in), contiguous :: blocks(:, :, :)
        complex(dp), intent(inout) :: values(:)

        call check_shapes(pattern, group, shape(blocks), size(values), "add")
        call add_complex(pattern%groups(group), pattern%nnz, blocks, values)
    end subroutine block_csc_add_complex

    subroutine add_complex(group, nnz, blocks, values)
        type(block_csc_group_t), intent(in) :: group
        integer, intent(in) :: nnz
        complex(dp), intent(in) :: blocks(*)
        complex(dp), intent(inout) :: values(:)
        complex(dp) :: total
        integer :: k, p

        !$omp parallel do default(none) shared(group, nnz, blocks, values) &
        !$omp private(p, k, total) schedule(static)
        do p = 1, nnz
            if (group%source_ptr(p) == group%source_ptr(p + 1)) cycle
            total = (0.0_dp, 0.0_dp)
            do k = group%source_ptr(p), group%source_ptr(p + 1) - 1
                total = total + blocks(group%source(k))
            end do
            values(p) = values(p) + total
        end do
        !$omp end parallel do
    end subroutine add_complex

    subroutine block_csc_gather_real(pattern, group, values_bar, blocks_bar)
        !! Adjoint of block_csc_add: blocks_bar(a, b, e) = values_bar at the
        !! position of entry (a, b, e), zero for dropped entries.
        type(block_csc_pattern_t), intent(in) :: pattern
        integer, intent(in) :: group
        real(dp), intent(in) :: values_bar(:)
        real(dp), intent(out) :: blocks_bar(:, :, :)

        call check_shapes(pattern, group, shape(blocks_bar), size(values_bar), &
            "gather")
        call gather_real(pattern%groups(group)%position, values_bar, blocks_bar)
    end subroutine block_csc_gather_real

    subroutine gather_real(position, values_bar, blocks_bar)
        integer, intent(in) :: position(:, :, :)
        real(dp), intent(in) :: values_bar(:)
        real(dp), intent(out) :: blocks_bar(:, :, :)
        integer :: a, b, e, p

        !$omp parallel do default(none) shared(position, values_bar, blocks_bar) &
        !$omp private(e, b, a, p) schedule(static)
        do e = 1, size(position, 3)
            do b = 1, size(position, 2)
                do a = 1, size(position, 1)
                    p = position(a, b, e)
                    if (p == 0) then
                        blocks_bar(a, b, e) = 0.0_dp
                    else
                        blocks_bar(a, b, e) = values_bar(p)
                    end if
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine gather_real

    subroutine block_csc_gather_complex(pattern, group, values_bar, blocks_bar)
        !! Adjoint of block_csc_add: blocks_bar(a, b, e) = values_bar at the
        !! position of entry (a, b, e), zero for dropped entries.
        type(block_csc_pattern_t), intent(in) :: pattern
        integer, intent(in) :: group
        complex(dp), intent(in) :: values_bar(:)
        complex(dp), intent(out) :: blocks_bar(:, :, :)

        call check_shapes(pattern, group, shape(blocks_bar), size(values_bar), &
            "gather")
        call gather_complex(pattern%groups(group)%position, values_bar, blocks_bar)
    end subroutine block_csc_gather_complex

    subroutine gather_complex(position, values_bar, blocks_bar)
        integer, intent(in) :: position(:, :, :)
        complex(dp), intent(in) :: values_bar(:)
        complex(dp), intent(out) :: blocks_bar(:, :, :)
        integer :: a, b, e, p

        !$omp parallel do default(none) shared(position, values_bar, blocks_bar) &
        !$omp private(e, b, a, p) schedule(static)
        do e = 1, size(position, 3)
            do b = 1, size(position, 2)
                do a = 1, size(position, 1)
                    p = position(a, b, e)
                    if (p == 0) then
                        blocks_bar(a, b, e) = (0.0_dp, 0.0_dp)
                    else
                        blocks_bar(a, b, e) = values_bar(p)
                    end if
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine gather_complex

    subroutine check_shapes(pattern, group, block_shape, value_count, operation)
        type(block_csc_pattern_t), intent(in) :: pattern
        integer, intent(in) :: group, block_shape(3), value_count
        character(*), intent(in) :: operation

        if (.not. pattern%finalized) &
            error stop "block_csc: pattern used before finalize"
        if (group < 1 .or. group > size(pattern%groups)) &
            error stop "block_csc: group index out of range"
        if (value_count /= pattern%nnz) then
            write (*, "(a)") "block_csc_"//operation//": values size /= nnz"
            error stop "block_csc: value array size mismatch"
        end if
        associate (g => pattern%groups(group))
            if (any(block_shape /= [g%nloc, g%nloc, g%nblocks])) then
                write (*, "(a)") "block_csc_"//operation//": block shape mismatch"
                error stop "block_csc: block shape mismatch"
            end if
        end associate
    end subroutine check_shapes

end module fortfem_block_csc
