module fortfem_element_to_csr
    !! Element blocks to a row-compressed (CSR) matrix.
    !!
    !! The pattern is symmetric: row i holds column j whenever dofs i and j share an
    !! element. Entries of a row are sorted by column. Two storage forms exist.
    !!
    !! full form      : the diagonal is stored and (A x)_i = sum_j a_ij x_j.
    !! exchange form  : only off-diagonals are stored and the diagonal is implied by
    !!                  a zero row sum, (A x)_i = sum_{j/=i} a_ij (x_j - x_i).
    !!
    !! The node-to-element map is kept in the matrix. Every pass over rows is
    !! row-parallel (OpenMP) and sums each row in a fixed order (ascending element,
    !! then local indices), so the values are bitwise independent of the thread count.
    use, intrinsic :: iso_fortran_env, only: real32
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: csr_matrix_t
    public :: csr_build_pattern
    public :: csr_assemble_blocks
    public :: csr_assemble_batched
    public :: csr_find_entry
    public :: csr_element_colouring
    public :: build_node_to_element

    type :: csr_matrix_t
        integer :: n = 0
        logical :: exchange = .false.
        integer :: nloc = 0
        integer, allocatable :: rowptr(:)   ! (n+1)
        integer, allocatable :: colidx(:)   ! (nnz), sorted within each row
        integer, allocatable :: tpos(:)     ! (nnz), position of the entry (j,i)
        real(dp), allocatable :: val(:)     ! (nnz)
        integer, allocatable :: elem_ptr(:) ! (n+1), node-to-element map
        integer, allocatable :: elem_idx(:) ! elements of dof i, ascending
        integer, allocatable :: slot(:, :, :) ! (nloc, nloc, nel): slot(b, a, e) is the
                                              ! colidx position of (conn(a, e),
                                              ! conn(b, e)), contiguous along the
                                              ! row a; 0 if dropped
    end type csr_matrix_t

    interface csr_assemble_batched
        module procedure csr_assemble_batched_fp64, csr_assemble_batched_fp32
    end interface

contains

    subroutine build_node_to_element(n_dofs, conn, elem_ptr, elem_idx)
        !! Counting sort of the element connectivity conn(nloc, nel). Elements of a
        !! dof come out in ascending order. A dof repeated within one element is
        !! listed once for that element.
        integer, intent(in) :: n_dofs, conn(:, :)
        integer, allocatable, intent(out) :: elem_ptr(:), elem_idx(:)

        integer, allocatable :: fill(:)
        logical(1), allocatable :: repeated(:, :)
        integer :: a, b, dof, e, nel, nloc

        nloc = size(conn, 1)
        nel = size(conn, 2)
        if (any(conn < 1) .or. any(conn > n_dofs)) then
            error stop "build_node_to_element: connectivity outside 1:n_dofs"
        end if
        allocate(elem_ptr(n_dofs + 1), fill(n_dofs), repeated(nloc, nel))
        ! The repeated-dof test is quadratic in nloc: evaluate it once, in
        ! parallel, instead of in both serial counting passes.
        !$omp parallel do default(shared) private(e, a) schedule(static)
        do e = 1, nel
            do a = 1, nloc
                repeated(a, e) = is_repeated(conn(:, e), a)
            end do
        end do
        !$omp end parallel do
        elem_ptr = 0
        do e = 1, nel
            do a = 1, nloc
                if (repeated(a, e)) cycle
                elem_ptr(conn(a, e) + 1) = elem_ptr(conn(a, e) + 1) + 1
            end do
        end do
        elem_ptr(1) = 1
        do dof = 1, n_dofs
            elem_ptr(dof + 1) = elem_ptr(dof + 1) + elem_ptr(dof)
        end do
        allocate(elem_idx(elem_ptr(n_dofs + 1) - 1))
        fill = elem_ptr(1:n_dofs)
        do e = 1, nel
            do a = 1, nloc
                if (repeated(a, e)) cycle
                b = conn(a, e)
                elem_idx(fill(b)) = e
                fill(b) = fill(b) + 1
            end do
        end do
    contains
        pure logical function is_repeated(column, index)
            integer, intent(in) :: column(:), index
            integer :: previous

            is_repeated = .false.
            do previous = 1, index - 1
                if (column(previous) == column(index)) then
                    is_repeated = .true.
                    return
                end if
            end do
        end function is_repeated
    end subroutine build_node_to_element

    subroutine csr_build_pattern(n_dofs, conn, exchange, matrix, with_slots)
        !! Symmetric sparsity pattern from conn(nloc, nel). The values are zeroed.
        !! with_slots=.true. also caches the element-to-CSR position map matrix%slot,
        !! which csr_assemble_blocks then uses instead of a search per entry.
        integer, intent(in) :: n_dofs, conn(:, :)
        logical, intent(in) :: exchange
        type(csr_matrix_t), intent(out) :: matrix
        logical, intent(in), optional :: with_slots

        integer, allocatable :: marker(:)
        integer :: a, b, count_row, e, i, j, k, p, position, start
        integer(8) :: total

        matrix%n = n_dofs
        matrix%exchange = exchange
        matrix%nloc = size(conn, 1)
        call build_node_to_element(n_dofs, conn, matrix%elem_ptr, matrix%elem_idx)
        allocate(matrix%rowptr(n_dofs + 1))

        ! Pass 1: row lengths.
        !$omp parallel default(none) shared(n_dofs, conn, exchange, matrix) &
        !$omp private(marker, i, k, e, b, j, count_row)
        allocate(marker(n_dofs))
        marker = 0
        !$omp do schedule(static)
        do i = 1, n_dofs
            count_row = 0
            if (.not. exchange) then
                marker(i) = i
                count_row = 1
            end if
            do k = matrix%elem_ptr(i), matrix%elem_ptr(i + 1) - 1
                e = matrix%elem_idx(k)
                do b = 1, size(conn, 1)
                    j = conn(b, e)
                    if (marker(j) == i) cycle
                    if (exchange .and. j == i) cycle
                    marker(j) = i
                    count_row = count_row + 1
                end do
            end do
            matrix%rowptr(i + 1) = count_row
        end do
        !$omp end do
        deallocate(marker)
        !$omp end parallel

        total = 1
        matrix%rowptr(1) = 1
        do i = 1, n_dofs
            total = total + matrix%rowptr(i + 1)
            if (total > int(huge(1), 8)) then
                error stop "csr_build_pattern: nnz exceeds default integer range"
            end if
            matrix%rowptr(i + 1) = int(total)
        end do
        allocate(matrix%colidx(matrix%rowptr(n_dofs + 1) - 1))
        allocate(matrix%tpos(matrix%rowptr(n_dofs + 1) - 1))
        allocate(matrix%val(matrix%rowptr(n_dofs + 1) - 1))

        ! Pass 2: columns, sorted per row. First touch of val by the row owner.
        !$omp parallel default(none) shared(n_dofs, conn, exchange, matrix) &
        !$omp private(marker, i, k, e, b, j, position, start, p)
        allocate(marker(n_dofs))
        marker = 0
        !$omp do schedule(static)
        do i = 1, n_dofs
            start = matrix%rowptr(i)
            position = start - 1
            if (.not. exchange) then
                marker(i) = i
                position = position + 1
                matrix%colidx(position) = i
            end if
            do k = matrix%elem_ptr(i), matrix%elem_ptr(i + 1) - 1
                e = matrix%elem_idx(k)
                do b = 1, size(conn, 1)
                    j = conn(b, e)
                    if (marker(j) == i) cycle
                    if (exchange .and. j == i) cycle
                    marker(j) = i
                    position = position + 1
                    matrix%colidx(position) = j
                end do
            end do
            call sort_row(matrix%colidx, start, position)
            do p = start, position
                matrix%val(p) = 0.0_dp
            end do
        end do
        !$omp end do
        deallocate(marker)
        !$omp end parallel

        ! Transposed entry positions (the pattern is symmetric).
        !$omp parallel do default(none) shared(n_dofs, matrix) private(i, k, p) &
        !$omp schedule(static)
        do i = 1, n_dofs
            do k = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                p = csr_find_entry(matrix, matrix%colidx(k), i)
                matrix%tpos(k) = p
            end do
        end do
        !$omp end parallel do
        if (any(matrix%tpos < 1)) error stop "csr_build_pattern: pattern not symmetric"
        if (present(with_slots)) then
            if (with_slots) then
                ! Row-parallel: marker holds the position of each column of row
                ! i, so a slot is a lookup instead of a search. Every element
                ! entry of row i is a column of row i.
                allocate(matrix%slot(matrix%nloc, matrix%nloc, size(conn, 2)))
                !$omp parallel default(none) shared(n_dofs, conn, matrix, exchange) &
                !$omp private(marker, i, k, e, a, b, j)
                allocate(marker(n_dofs))
                !$omp do schedule(static)
                do i = 1, n_dofs
                    do k = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                        marker(matrix%colidx(k)) = k
                    end do
                    do k = matrix%elem_ptr(i), matrix%elem_ptr(i + 1) - 1
                        e = matrix%elem_idx(k)
                        do a = 1, matrix%nloc
                            if (conn(a, e) /= i) cycle
                            do b = 1, matrix%nloc
                                j = conn(b, e)
                                if (exchange .and. j == i) then
                                    matrix%slot(b, a, e) = 0
                                else
                                    matrix%slot(b, a, e) = marker(j)
                                end if
                            end do
                        end do
                    end do
                end do
                !$omp end do
                deallocate(marker)
                !$omp end parallel
            end if
        end if
    end subroutine csr_build_pattern

    pure subroutine sort_row(column, first, last)
        !! Ascending sort of column(first:last): insertion sort for short rows,
        !! heapsort (O(n log n)) for the long rows of high-order or Fourier-coupled
        !! elements, where the quadratic insertion sort dominated the pattern build.
        integer, intent(inout) :: column(:)
        integer, intent(in) :: first, last

        integer, parameter :: short_row = 32
        integer :: i, j, key

        if (last - first + 1 > short_row) then
            call heap_sort(column(first:last))
            return
        end if
        do i = first + 1, last
            key = column(i)
            j = i - 1
            do while (j >= first)
                if (column(j) <= key) exit
                column(j + 1) = column(j)
                j = j - 1
            end do
            column(j + 1) = key
        end do
    end subroutine sort_row

    pure subroutine heap_sort(v)
        integer, intent(inout) :: v(:)

        integer :: n, last, t

        n = size(v)
        do last = n/2, 1, -1
            call sift_down(v, last, n)
        end do
        do last = n, 2, -1
            t = v(1); v(1) = v(last); v(last) = t
            call sift_down(v, 1, last - 1)
        end do
    contains
        pure subroutine sift_down(v, root, bound)
            integer, intent(inout) :: v(:)
            integer, intent(in) :: root, bound
            integer :: child, parent, key

            parent = root
            key = v(parent)
            do
                child = 2*parent
                if (child > bound) exit
                if (child < bound) then
                    if (v(child + 1) > v(child)) child = child + 1
                end if
                if (v(child) <= key) exit
                v(parent) = v(child)
                parent = child
            end do
            v(parent) = key
        end subroutine sift_down
    end subroutine heap_sort

    pure function csr_find_entry(matrix, row, column) result(position)
        !! Position of entry (row, column) in colidx, or 0 when it is not stored.
        type(csr_matrix_t), intent(in) :: matrix
        integer, intent(in) :: row, column
        integer :: position

        integer :: high, low, middle

        position = 0
        low = matrix%rowptr(row)
        high = matrix%rowptr(row + 1) - 1
        do while (low <= high)
            middle = low + (high - low)/2
            if (matrix%colidx(middle) == column) then
                position = middle
                return
            else if (matrix%colidx(middle) < column) then
                low = middle + 1
            else
                high = middle - 1
            end if
        end do
    end function csr_find_entry

    subroutine csr_assemble_blocks(matrix, conn, blocks, row_sum_defect, symmetric)
        !! Overwrite the values with the sum of the element blocks
        !! blocks(a, b, e) = K_ab of element e (row a, column b). The pattern must
        !! come from csr_build_pattern with the same conn. In exchange form the
        !! diagonal is dropped; row_sum_defect (optional) returns the largest
        !! |sum_b K_ab| summed over a row relative to the largest absolute entry,
        !! which is the information lost by dropping it.
        !! symmetric=.true. (with the slots of csr_build_pattern) declares every
        !! block exactly symmetric, blocks(a, b, e) == blocks(b, a, e), and reads
        !! row a of a block as its contiguous column a. The sums run in the same
        !! order, so the values are bitwise those of the default path.
        type(csr_matrix_t), intent(inout) :: matrix
        integer, intent(in) :: conn(:, :)
        real(dp), intent(in) :: blocks(:, :, :)
        real(dp), intent(out), optional :: row_sum_defect
        logical, intent(in), optional :: symmetric

        real(dp) :: defect, row_sum, scale
        integer :: a, b, e, i, j, k, p
        logical :: columns

        if (size(blocks, 1) /= matrix%nloc .or. size(blocks, 2) /= matrix%nloc .or. &
            size(blocks, 3) /= size(conn, 2) .or. size(conn, 1) /= matrix%nloc) then
            error stop "csr_assemble_blocks: block or connectivity shape mismatch"
        end if
        columns = .false.
        if (present(symmetric)) columns = symmetric .and. allocated(matrix%slot) &
            .and. .not. present(row_sum_defect)
        if (columns) then
            !$omp parallel do default(none) shared(matrix, conn, blocks) &
            !$omp private(i, k, e, a, b, p) schedule(static)
            do i = 1, matrix%n
                matrix%val(matrix%rowptr(i):matrix%rowptr(i + 1) - 1) = 0.0_dp
                do k = matrix%elem_ptr(i), matrix%elem_ptr(i + 1) - 1
                    e = matrix%elem_idx(k)
                    do a = 1, matrix%nloc
                        if (conn(a, e) /= i) cycle
                        do b = 1, matrix%nloc
                            p = matrix%slot(b, a, e)
                            if (p == 0) cycle
                            matrix%val(p) = matrix%val(p) + blocks(b, a, e)
                        end do
                    end do
                end do
            end do
            !$omp end parallel do
            return
        end if
        defect = 0.0_dp
        scale = 0.0_dp
        !$omp parallel do default(none) shared(matrix, conn, blocks) &
        !$omp private(i, k, e, a, b, j, p, row_sum) schedule(static) &
        !$omp reduction(max: defect, scale)
        do i = 1, matrix%n
            matrix%val(matrix%rowptr(i):matrix%rowptr(i + 1) - 1) = 0.0_dp
            row_sum = 0.0_dp
            do k = matrix%elem_ptr(i), matrix%elem_ptr(i + 1) - 1
                e = matrix%elem_idx(k)
                do a = 1, matrix%nloc
                    if (conn(a, e) /= i) cycle
                    do b = 1, matrix%nloc
                        j = conn(b, e)
                        row_sum = row_sum + blocks(a, b, e)
                        scale = max(scale, abs(blocks(a, b, e)))
                        if (matrix%exchange .and. j == i) cycle
                        if (allocated(matrix%slot)) then
                            p = matrix%slot(b, a, e)
                        else
                            p = csr_find_entry(matrix, i, j)
                        end if
                        matrix%val(p) = matrix%val(p) + blocks(a, b, e)
                    end do
                end do
            end do
            defect = max(defect, abs(row_sum))
        end do
        !$omp end parallel do
        if (present(row_sum_defect)) then
            row_sum_defect = 0.0_dp
            if (scale > 0.0_dp) row_sum_defect = defect/scale
        end if
    end subroutine csr_assemble_blocks

    subroutine csr_assemble_batched_fp64(matrix, conn, blocks)
        type(csr_matrix_t), intent(inout) :: matrix
        integer, intent(in) :: conn(:, :)
        real(dp), intent(in) :: blocks(:, :, :)
        ! Same connectivity contract as csr_assemble_blocks. SoA cell storage
        ! stays narrow; conversion occurs in registers before fp64 accumulation.
        integer :: a, b, e, i, j, k, p

        if (size(blocks,1) /= size(conn,2) .or. &
            size(blocks,2) /= matrix%nloc .or. &
            size(blocks,3) /= matrix%nloc .or. &
            size(conn,1) /= matrix%nloc) &
            error stop 'csr_assemble_batched: shape mismatch'
        !$omp parallel do default(none) shared(matrix, conn, blocks) &
        !$omp private(i, k, e, a, b, j, p) schedule(static)
        do i = 1, matrix%n
            matrix%val(matrix%rowptr(i):matrix%rowptr(i+1)-1) = 0.0_dp
            do k = matrix%elem_ptr(i), matrix%elem_ptr(i+1)-1
                e = matrix%elem_idx(k)
                do a = 1, matrix%nloc
                    if (conn(a,e) /= i) cycle
                    do b = 1, matrix%nloc
                        j = conn(b,e)
                        if (matrix%exchange .and. j == i) cycle
                        p = csr_find_entry(matrix,i,j)
                        matrix%val(p) = matrix%val(p) + real(blocks(e,a,b),dp)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine csr_assemble_batched_fp64

    subroutine csr_assemble_batched_fp32(matrix, conn, blocks)
        type(csr_matrix_t), intent(inout) :: matrix
        integer, intent(in) :: conn(:, :)
        real(real32), intent(in) :: blocks(:, :, :)
        ! Same connectivity contract as csr_assemble_blocks. SoA cell storage
        ! stays narrow; conversion occurs in registers before fp64 accumulation.
        integer :: a, b, e, i, j, k, p

        if (size(blocks,1) /= size(conn,2) .or. &
            size(blocks,2) /= matrix%nloc .or. &
            size(blocks,3) /= matrix%nloc .or. &
            size(conn,1) /= matrix%nloc) &
            error stop 'csr_assemble_batched: shape mismatch'
        !$omp parallel do default(none) shared(matrix, conn, blocks) &
        !$omp private(i, k, e, a, b, j, p) schedule(static)
        do i = 1, matrix%n
            matrix%val(matrix%rowptr(i):matrix%rowptr(i+1)-1) = 0.0_dp
            do k = matrix%elem_ptr(i), matrix%elem_ptr(i+1)-1
                e = matrix%elem_idx(k)
                do a = 1, matrix%nloc
                    if (conn(a,e) /= i) cycle
                    do b = 1, matrix%nloc
                        j = conn(b,e)
                        if (matrix%exchange .and. j == i) cycle
                        p = csr_find_entry(matrix,i,j)
                        matrix%val(p) = matrix%val(p) + real(blocks(e,a,b),dp)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine csr_assemble_batched_fp32

    subroutine csr_element_colouring(matrix, conn, colour, n_colours)
        !! Greedy element colouring: elements sharing a dof get different colours.
        !! Serial and deterministic (ascending element order).
        type(csr_matrix_t), intent(in) :: matrix
        integer, intent(in) :: conn(:, :)
        integer, allocatable, intent(out) :: colour(:)
        integer, intent(out) :: n_colours

        integer, allocatable :: stamp(:)
        integer :: a, c, e, f, k, nel

        nel = size(conn, 2)
        allocate(colour(nel), stamp(0:nel + 1))
        colour = 0
        stamp = 0
        n_colours = 0
        do e = 1, nel
            do a = 1, size(conn, 1)
                do k = matrix%elem_ptr(conn(a, e)), matrix%elem_ptr(conn(a, e) + 1) - 1
                    f = matrix%elem_idx(k)
                    if (f /= e .and. colour(f) > 0) stamp(colour(f)) = e
                end do
            end do
            c = 1
            do while (c <= n_colours)
                if (stamp(c) /= e) exit
                c = c + 1
            end do
            colour(e) = c
            n_colours = max(n_colours, c)
        end do
    end subroutine csr_element_colouring

end module fortfem_element_to_csr
