program test_block_csc
    !! Block groups to one CSC matrix: dense reference for mixed groups (vector P1
    !! triangles, boundary edges, one dense interface block, dropped dofs), real and
    !! complex; gather as the exact adjoint of add; linearity of add (its own JVP);
    !! values independent of the OpenMP thread count.
    use check, only: check_condition, check_summary
    use fortfem_block_csc, only: block_csc_add, block_csc_add_group, &
        block_csc_finalize, block_csc_gather, block_csc_init, &
        block_csc_new_matrix, block_csc_pattern_t
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_is_valid, csc_t, csc_z_t
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    implicit none

    call check_mixed_groups()
    call check_thread_invariance()
    call check_summary("block_csc")

contains

    subroutine check_mixed_groups()
        type(block_csc_pattern_t) :: pattern
        type(csc_t) :: real_matrix
        type(csc_z_t) :: complex_matrix
        integer, allocatable :: cells(:, :), edges(:, :), interface(:, :)
        real(dp), allocatable :: cell_blocks(:, :, :), edge_blocks(:, :, :)
        real(dp), allocatable :: interface_blocks(:, :, :), dense(:, :)
        real(dp), allocatable :: values_bar(:), cell_bar(:, :, :), tangent(:, :, :)
        complex(dp), allocatable :: zcell(:, :, :), zedge(:, :, :)
        complex(dp), allocatable :: zinterface(:, :, :), zdense(:, :), zbar(:)
        complex(dp), allocatable :: zcell_bar(:, :, :), zedge_bar(:, :, :)
        complex(dp), allocatable :: zinterface_bar(:, :, :)
        real(dp), allocatable :: base(:), shifted(:), linear(:)
        complex(dp) :: left, right
        real(dp) :: real_left, real_right
        integer :: cell_group, edge_group, interface_group, n, nx

        nx = 5
        call vector_p1_grid(nx, cells, edges, interface)
        n = 2*nx*nx
        ! Drop the y component of vertex 1 and both components of vertex 2 (as
        ! a caller removing constrained dofs would).
        where (cells == 2 .or. cells == 3 .or. cells == 4) cells = 0
        where (edges == 2 .or. edges == 3 .or. edges == 4) edges = 0

        call block_csc_init(pattern, n)
        call block_csc_add_group(pattern, cells, cell_group)
        call block_csc_add_group(pattern, edges, edge_group)
        call block_csc_add_group(pattern, interface, interface_group)
        call block_csc_finalize(pattern)
        call check_condition(cell_group == 1 .and. edge_group == 2 .and. &
            interface_group == 3, "groups are numbered in the order added")

        call random_blocks(cells, 11, cell_blocks)
        call random_blocks(edges, 12, edge_blocks)
        call random_blocks(interface, 13, interface_blocks)
        call block_csc_new_matrix(pattern, real_matrix)
        call check_condition(csc_is_valid(real_matrix) .and. &
            real_matrix%nrow == n .and. real_matrix%ncol == n, &
            "real matrix is a valid n x n CSC")
        call check_condition(rows_sorted(real_matrix%col_ptr, real_matrix%row_idx), &
            "row indices sorted within each column")
        call block_csc_add(pattern, cell_group, cell_blocks, real_matrix%val)
        call block_csc_add(pattern, edge_group, edge_blocks, real_matrix%val)
        call block_csc_add(pattern, interface_group, interface_blocks, &
            real_matrix%val)
        allocate(dense(n, n))
        dense = 0.0_dp
        call dense_add(cells, cell_blocks, dense)
        call dense_add(edges, edge_blocks, dense)
        call dense_add(interface, interface_blocks, dense)
        call check_condition(maxval(abs(expand(real_matrix%col_ptr, &
            real_matrix%row_idx, real_matrix%val, n) - dense)) <= &
            1.0e-14_dp*maxval(abs(dense)), "real: sparse equals dense reference")
        call check_condition(pattern_is_coupling(real_matrix%col_ptr, &
            real_matrix%row_idx, n, [cells_flat(cells), cells_flat(edges), &
            cells_flat(interface)]), "pattern is exactly the coupled dof pairs")

        ! Assembly is linear, so its JVP is add of the tangents.
        allocate(tangent, mold=cell_blocks)
        call random_blocks(cells, 14, tangent)
        base = real_matrix%val
        real_matrix%val = 0.0_dp
        call block_csc_add(pattern, cell_group, tangent, real_matrix%val)
        linear = real_matrix%val
        real_matrix%val = base
        call block_csc_add(pattern, cell_group, tangent, real_matrix%val)
        shifted = real_matrix%val
        call check_condition(maxval(abs(shifted - base - linear)) <= &
            1.0e-14_dp*maxval(abs(shifted)), "add accumulates and is linear")

        allocate(values_bar(real_matrix%nnz))
        allocate(cell_bar, mold=cell_blocks)
        call random_number(values_bar)
        call block_csc_gather(pattern, cell_group, values_bar, cell_bar)
        real_left = dot_product(values_bar, linear)
        real_right = sum(cell_bar*tangent)
        call check_condition(abs(real_left - real_right) <= &
            1.0e-13_dp*abs(real_left), "real: gather is the adjoint of add")
        call check_condition(all(dropped_zero(cells, cell_bar)), &
            "gather gives zero for dropped entries")

        call complex_blocks(cells, 21, zcell)
        call complex_blocks(edges, 22, zedge)
        call complex_blocks(interface, 23, zinterface)
        call block_csc_new_matrix(pattern, complex_matrix)
        call check_condition(csc_is_valid(complex_matrix) .and. &
            all(complex_matrix%row_idx == real_matrix%row_idx), &
            "complex matrix shares the pattern")
        call block_csc_add(pattern, cell_group, zcell, complex_matrix%val)
        call block_csc_add(pattern, edge_group, zedge, complex_matrix%val)
        call block_csc_add(pattern, interface_group, zinterface, &
            complex_matrix%val)
        allocate(zdense(n, n))
        zdense = (0.0_dp, 0.0_dp)
        call dense_add_z(cells, zcell, zdense)
        call dense_add_z(edges, zedge, zdense)
        call dense_add_z(interface, zinterface, zdense)
        call check_condition(maxval(abs(expand_z(complex_matrix, n) - zdense)) <= &
            1.0e-14_dp*maxval(abs(zdense)), "complex: sparse equals dense reference")

        allocate(zbar(complex_matrix%nnz))
        call random_number(values_bar)
        zbar = cmplx(values_bar, 1.0_dp - 2.0_dp*values_bar(size(values_bar):1:-1), &
            dp)
        allocate(zcell_bar, mold=zcell)
        allocate(zedge_bar, mold=zedge)
        allocate(zinterface_bar, mold=zinterface)
        call block_csc_gather(pattern, cell_group, zbar, zcell_bar)
        call block_csc_gather(pattern, edge_group, zbar, zedge_bar)
        call block_csc_gather(pattern, interface_group, zbar, zinterface_bar)
        left = sum(conjg(zbar)*complex_matrix%val)
        right = sum(conjg(zcell_bar)*zcell) + sum(conjg(zedge_bar)*zedge) + &
            sum(conjg(zinterface_bar)*zinterface)
        call check_condition(abs(left - right) <= 1.0e-13_dp*abs(left), &
            "complex: <V, add(B)> = <gather(V), B> over all groups")
    end subroutine check_mixed_groups

    subroutine check_thread_invariance()
        type(csc_z_t) :: one, many
        integer, allocatable :: cells(:, :), edges(:, :), interface(:, :)
        complex(dp), allocatable :: zcell(:, :, :), zedge(:, :, :)
        complex(dp), allocatable :: zinterface(:, :, :)
        integer :: nx, threads

        threads = 1
        !$ threads = omp_get_max_threads()
        nx = 60
        call vector_p1_grid(nx, cells, edges, interface)
        call complex_blocks(cells, 31, zcell)
        call complex_blocks(edges, 32, zedge)
        call complex_blocks(interface, 33, zinterface)
        !$ call omp_set_num_threads(1)
        call assemble_groups(2*nx*nx, cells, edges, interface, zcell, zedge, &
            zinterface, one)
        !$ call omp_set_num_threads(max(4, threads))
        call assemble_groups(2*nx*nx, cells, edges, interface, zcell, zedge, &
            zinterface, many)
        !$ call omp_set_num_threads(threads)
        call check_condition(all(one%row_idx == many%row_idx) .and. &
            all(one%col_ptr == many%col_ptr), "pattern independent of threads")
        call check_condition(all(one%val == many%val), &
            "values bitwise independent of threads")
    end subroutine check_thread_invariance

    subroutine assemble_groups(n, cells, edges, interface, zcell, zedge, &
            zinterface, matrix)
        integer, intent(in) :: n, cells(:, :), edges(:, :), interface(:, :)
        complex(dp), intent(in) :: zcell(:, :, :), zedge(:, :, :)
        complex(dp), intent(in) :: zinterface(:, :, :)
        type(csc_z_t), intent(out) :: matrix
        type(block_csc_pattern_t) :: pattern
        integer :: cell_group, edge_group, interface_group

        call block_csc_init(pattern, n)
        call block_csc_add_group(pattern, cells, cell_group)
        call block_csc_add_group(pattern, edges, edge_group)
        call block_csc_add_group(pattern, interface, interface_group)
        call block_csc_finalize(pattern)
        call block_csc_new_matrix(pattern, matrix)
        call block_csc_add(pattern, cell_group, zcell, matrix%val)
        call block_csc_add(pattern, edge_group, zedge, matrix%val)
        call block_csc_add(pattern, interface_group, zinterface, matrix%val)
    end subroutine assemble_groups

    subroutine vector_p1_grid(nx, cells, edges, interface)
        !! Vector dofs (2*node - 1, 2*node) of a structured nx x nx triangle grid:
        !! triangles (6 dofs), bottom boundary edges (4 dofs) and one dense block
        !! coupling every dof of the top row (2*nx dofs).
        integer, intent(in) :: nx
        integer, allocatable, intent(out) :: cells(:, :), edges(:, :)
        integer, allocatable, intent(out) :: interface(:, :)
        integer :: i, j, k, node(4), triangle(3, 2), t

        allocate(cells(6, 2*(nx - 1)**2), edges(4, nx - 1), interface(2*nx, 1))
        k = 0
        do j = 1, nx - 1
            do i = 1, nx - 1
                node = [i + (j - 1)*nx, i + 1 + (j - 1)*nx, i + 1 + j*nx, i + j*nx]
                triangle(:, 1) = [node(1), node(2), node(3)]
                triangle(:, 2) = [node(1), node(3), node(4)]
                do t = 1, 2
                    k = k + 1
                    cells(:, k) = vector_dofs(triangle(:, t))
                end do
            end do
        end do
        do i = 1, nx - 1
            edges(:, i) = vector_dofs([i, i + 1])
        end do
        interface(:, 1) = vector_dofs([(i + (nx - 1)*nx, i=1, nx)])
    end subroutine vector_p1_grid

    pure function vector_dofs(nodes) result(dofs)
        integer, intent(in) :: nodes(:)
        integer :: dofs(2*size(nodes))

        dofs(1::2) = 2*nodes - 1
        dofs(2::2) = 2*nodes
    end function vector_dofs

    subroutine random_blocks(conn, seed, blocks)
        integer, intent(in) :: conn(:, :), seed
        real(dp), allocatable, intent(inout) :: blocks(:, :, :)

        call set_seed(seed)
        if (.not. allocated(blocks)) &
            allocate(blocks(size(conn, 1), size(conn, 1), size(conn, 2)))
        call random_number(blocks)
        blocks = 2.0_dp*blocks - 1.0_dp
    end subroutine random_blocks

    subroutine complex_blocks(conn, seed, blocks)
        integer, intent(in) :: conn(:, :), seed
        complex(dp), allocatable, intent(out) :: blocks(:, :, :)
        real(dp), allocatable :: re(:, :, :), im(:, :, :)

        call random_blocks(conn, seed, re)
        call random_blocks(conn, seed + 100, im)
        blocks = cmplx(re, im, dp)
    end subroutine complex_blocks

    subroutine set_seed(seed)
        integer, intent(in) :: seed
        integer, allocatable :: state(:)
        integer :: k, size_state

        call random_seed(size=size_state)
        allocate(state(size_state))
        state = seed + 37*[(k, k=1, size_state)]
        call random_seed(put=state)
    end subroutine set_seed

    pure subroutine dense_add(conn, blocks, dense)
        integer, intent(in) :: conn(:, :)
        real(dp), intent(in) :: blocks(:, :, :)
        real(dp), intent(inout) :: dense(:, :)
        integer :: a, b, e

        do e = 1, size(conn, 2)
            do b = 1, size(conn, 1)
                if (conn(b, e) <= 0) cycle
                do a = 1, size(conn, 1)
                    if (conn(a, e) <= 0) cycle
                    dense(conn(a, e), conn(b, e)) = dense(conn(a, e), conn(b, e)) &
                        + blocks(a, b, e)
                end do
            end do
        end do
    end subroutine dense_add

    pure subroutine dense_add_z(conn, blocks, dense)
        integer, intent(in) :: conn(:, :)
        complex(dp), intent(in) :: blocks(:, :, :)
        complex(dp), intent(inout) :: dense(:, :)
        integer :: a, b, e

        do e = 1, size(conn, 2)
            do b = 1, size(conn, 1)
                if (conn(b, e) <= 0) cycle
                do a = 1, size(conn, 1)
                    if (conn(a, e) <= 0) cycle
                    dense(conn(a, e), conn(b, e)) = dense(conn(a, e), conn(b, e)) &
                        + blocks(a, b, e)
                end do
            end do
        end do
    end subroutine dense_add_z

    pure function expand(col_ptr, row_idx, values, n) result(dense)
        integer, intent(in) :: col_ptr(:), row_idx(:), n
        real(dp), intent(in) :: values(:)
        real(dp) :: dense(n, n)
        integer :: column, entry

        dense = 0.0_dp
        do column = 1, n
            do entry = col_ptr(column), col_ptr(column + 1) - 1
                dense(row_idx(entry), column) = values(entry)
            end do
        end do
    end function expand

    pure function expand_z(matrix, n) result(dense)
        type(csc_z_t), intent(in) :: matrix
        integer, intent(in) :: n
        complex(dp) :: dense(n, n)
        integer :: column, entry

        dense = (0.0_dp, 0.0_dp)
        do column = 1, n
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                dense(matrix%row_idx(entry), column) = matrix%val(entry)
            end do
        end do
    end function expand_z

    pure logical function rows_sorted(col_ptr, row_idx)
        integer, intent(in) :: col_ptr(:), row_idx(:)
        integer :: column, entry

        rows_sorted = .true.
        do column = 1, size(col_ptr) - 1
            do entry = col_ptr(column) + 1, col_ptr(column + 1) - 1
                if (row_idx(entry) <= row_idx(entry - 1)) rows_sorted = .false.
            end do
        end do
    end function rows_sorted

    function cells_flat(conn) result(pairs)
        !! Coupled (row, column) pairs of one group, encoded row + n*column.
        integer, intent(in) :: conn(:, :)
        integer, allocatable :: pairs(:)
        integer :: a, b, e, k

        allocate(pairs(size(conn, 1)**2*size(conn, 2)))
        k = 0
        do e = 1, size(conn, 2)
            do b = 1, size(conn, 1)
                do a = 1, size(conn, 1)
                    k = k + 1
                    pairs(k) = -1
                    if (conn(a, e) > 0 .and. conn(b, e) > 0) &
                        pairs(k) = conn(a, e) + 100000*conn(b, e)
                end do
            end do
        end do
    end function cells_flat

    logical function pattern_is_coupling(col_ptr, row_idx, n, pairs)
        integer, intent(in) :: col_ptr(:), row_idx(:), n, pairs(:)
        logical, allocatable :: coupled(:, :), stored(:, :)
        integer :: column, entry, k

        allocate(coupled(n, n), stored(n, n))
        coupled = .false.
        stored = .false.
        do k = 1, size(pairs)
            if (pairs(k) < 0) cycle
            coupled(mod(pairs(k), 100000), pairs(k)/100000) = .true.
        end do
        do column = 1, n
            do entry = col_ptr(column), col_ptr(column + 1) - 1
                stored(row_idx(entry), column) = .true.
            end do
        end do
        pattern_is_coupling = all(coupled .eqv. stored)
    end function pattern_is_coupling

    pure function dropped_zero(conn, blocks) result(ok)
        integer, intent(in) :: conn(:, :)
        real(dp), intent(in) :: blocks(:, :, :)
        logical :: ok(size(conn, 2))
        integer :: a, b, e

        ok = .true.
        do e = 1, size(conn, 2)
            do b = 1, size(conn, 1)
                do a = 1, size(conn, 1)
                    if (conn(a, e) > 0 .and. conn(b, e) > 0) cycle
                    if (blocks(a, b, e) /= 0.0_dp) ok(e) = .false.
                end do
            end do
        end do
    end function dropped_zero

end program test_block_csc
