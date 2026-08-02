module fortfem_tree_cotree_iga_parity
    !! Orientation parity for a fixed tree-cotree gauge on IGA DOFs.
    !!
    !! Two signed local-to-global maps that differ only by a global sign for
    !! each equivalence class describe the same physical H(curl) basis.  This
    !! contract assembles both maps from identical local data, applies the
    !! high-order/IGA tree-cotree selector, and checks the direct system,
    !! solution, and oriented period under that signed change of basis.
    use fortfem_kinds, only: dp
    use fortfem_tree_cotree_gauge, only: build_tree_cotree_dof_map, &
        build_tree_cotree_gauge, tree_cotree_gauge_t
    implicit none
    private

    type, public :: tree_cotree_iga_parity_t
        integer :: global_dof_count = 0
        integer :: free_dof_count = 0
        real(dp) :: matrix_error = huge(1.0_dp)
        real(dp) :: rhs_error = huge(1.0_dp)
        real(dp) :: solution_error = huge(1.0_dp)
        real(dp) :: period_error = huge(1.0_dp)
        real(dp) :: period_value = 0.0_dp
    end type tree_cotree_iga_parity_t

    public :: evaluate_tree_cotree_iga_parity

contains

    subroutine evaluate_tree_cotree_iga_parity( &
            incidence, control_edge_local, signed_map_a, signed_map_b, &
            local_matrix, local_rhs, period_weights, diagnostics, status)
        integer, intent(in) :: incidence(:, :), control_edge_local(:)
        integer, intent(in) :: signed_map_a(:), signed_map_b(:)
        real(dp), intent(in) :: local_matrix(:, :), local_rhs(:)
        real(dp), intent(in) :: period_weights(:)
        type(tree_cotree_iga_parity_t), intent(out) :: diagnostics
        integer, intent(out) :: status

        type(tree_cotree_gauge_t) :: gauge
        logical, allocatable :: constrained(:)
        integer, allocatable :: free_dofs(:), global_sign(:)
        real(dp), allocatable :: matrix_a(:, :), matrix_b(:, :)
        real(dp), allocatable :: rhs_a(:), rhs_b(:)
        real(dp), allocatable :: reduced_a(:, :), reduced_b(:, :)
        real(dp), allocatable :: reduced_rhs_a(:), reduced_rhs_b(:)
        real(dp), allocatable :: solution_a(:), solution_b(:)
        real(dp), allocatable :: full_solution_a(:), full_solution_b(:)
        integer :: local_count, global_count, edge_count, edge, local
        integer :: i, j, row, column, id_a, id_b, sign_a, sign_b
        integer :: gauge_status, solve_status
        real(dp) :: expected, period_a, period_b

        diagnostics = tree_cotree_iga_parity_t()
        status = 1
        local_count = size(signed_map_a)
        if (local_count < 1 .or. size(signed_map_b) /= local_count) return
        if (size(local_matrix, 1) /= local_count .or. &
            size(local_matrix, 2) /= local_count .or. &
            size(local_rhs) /= local_count) return
        if (size(incidence, 2) /= size(control_edge_local)) return
        if (size(period_weights) < 1) return
        if (any(signed_map_a == 0) .or. any(signed_map_b == 0)) return
        if (any(abs(signed_map_a) /= abs(signed_map_b))) return
        global_count = maxval(abs(signed_map_a))
        if (global_count < 1 .or. size(period_weights) /= global_count) return
        diagnostics%global_dof_count = global_count

        allocate(global_sign(global_count))
        global_sign = 0
        do local = 1, local_count
            id_a = abs(signed_map_a(local))
            id_b = abs(signed_map_b(local))
            if (id_a /= id_b) return
            sign_a = sign_of(signed_map_a(local))
            sign_b = sign_of(signed_map_b(local))
            if (global_sign(id_a) == 0) then
                global_sign(id_a) = sign_a*sign_b
            else if (global_sign(id_a) /= sign_a*sign_b) then
                status = 2
                return
            end if
        end do
        if (any(global_sign == 0)) return

        edge_count = size(incidence, 2)
        do edge = 1, edge_count
            if (control_edge_local(edge) < 1 .or. &
                control_edge_local(edge) > local_count) return
        end do
        call build_tree_cotree_gauge(incidence, gauge, gauge_status)
        if (gauge_status /= 0) return
        allocate(matrix_a(global_count, global_count), &
            matrix_b(global_count, global_count), rhs_a(global_count), &
            rhs_b(global_count))
        matrix_a = 0.0_dp
        matrix_b = 0.0_dp
        rhs_a = 0.0_dp
        rhs_b = 0.0_dp
        do local = 1, local_count
            id_a = abs(signed_map_a(local))
            id_b = abs(signed_map_b(local))
            sign_a = sign_of(signed_map_a(local))
            sign_b = sign_of(signed_map_b(local))
            rhs_a(id_a) = rhs_a(id_a) + real(sign_a, dp)*local_rhs(local)
            rhs_b(id_b) = rhs_b(id_b) + real(sign_b, dp)*local_rhs(local)
            do j = 1, local_count
                id_b = abs(signed_map_b(j))
                sign_b = sign_of(signed_map_b(j))
                matrix_a(id_a, abs(signed_map_a(j))) = &
                    matrix_a(id_a, abs(signed_map_a(j))) + &
                    real(sign_a*sign_of(signed_map_a(j)), dp)* &
                    local_matrix(local, j)
                matrix_b(id_b, abs(signed_map_b(j))) = &
                    matrix_b(id_b, abs(signed_map_b(j))) + &
                    real(sign_b*sign_of(signed_map_b(j)), dp)* &
                    local_matrix(local, j)
            end do
        end do

        call build_tree_cotree_dof_map( &
            gauge, abs(signed_map_a(control_edge_local)), global_count, &
            constrained, free_dofs, gauge_status)
        if (gauge_status /= 0 .or. size(free_dofs) < 1) return
        diagnostics%free_dof_count = size(free_dofs)
        allocate(reduced_a(size(free_dofs), size(free_dofs)), &
            reduced_b(size(free_dofs), size(free_dofs)), &
            reduced_rhs_a(size(free_dofs)), reduced_rhs_b(size(free_dofs)))
        do row = 1, size(free_dofs)
            reduced_rhs_a(row) = rhs_a(free_dofs(row))
            reduced_rhs_b(row) = rhs_b(free_dofs(row))
            do column = 1, size(free_dofs)
                reduced_a(row, column) = &
                    matrix_a(free_dofs(row), free_dofs(column))
                reduced_b(row, column) = &
                    matrix_b(free_dofs(row), free_dofs(column))
            end do
        end do
        diagnostics%matrix_error = 0.0_dp
        diagnostics%rhs_error = 0.0_dp
        do row = 1, size(free_dofs)
            do column = 1, size(free_dofs)
                expected = real(global_sign(free_dofs(row))* &
                    global_sign(free_dofs(column)), dp)*reduced_a(row, column)
                diagnostics%matrix_error = max(diagnostics%matrix_error, &
                    abs(reduced_b(row, column) - expected))
            end do
            expected = real(global_sign(free_dofs(row)), dp)*reduced_rhs_a(row)
            diagnostics%rhs_error = max(diagnostics%rhs_error, &
                abs(reduced_rhs_b(row) - expected))
        end do

        allocate(solution_a(size(free_dofs)), solution_b(size(free_dofs)))
        call solve_dense_system(reduced_a, reduced_rhs_a, solution_a, solve_status)
        if (solve_status /= 0) then
            status = 3
            return
        end if
        call solve_dense_system(reduced_b, reduced_rhs_b, solution_b, solve_status)
        if (solve_status /= 0) then
            status = 3
            return
        end if
        diagnostics%solution_error = 0.0_dp
        do row = 1, size(free_dofs)
            expected = real(global_sign(free_dofs(row)), dp)*solution_a(row)
            diagnostics%solution_error = max(diagnostics%solution_error, &
                abs(solution_b(row) - expected))
        end do
        allocate(full_solution_a(global_count), full_solution_b(global_count))
        full_solution_a = 0.0_dp
        full_solution_b = 0.0_dp
        do row = 1, size(free_dofs)
            full_solution_a(free_dofs(row)) = solution_a(row)
            full_solution_b(free_dofs(row)) = solution_b(row)
        end do
        period_a = dot_product(period_weights, full_solution_a)
        period_b = dot_product(real(global_sign, dp)*period_weights, &
            full_solution_b)
        diagnostics%period_value = period_a
        diagnostics%period_error = abs(period_b - period_a)
        status = 0
    end subroutine evaluate_tree_cotree_iga_parity

    pure integer function sign_of(value) result(sign)
        integer, intent(in) :: value

        if (value < 0) then
            sign = -1
        else
            sign = 1
        end if
    end function sign_of

    subroutine solve_dense_system(matrix, rhs, solution, status)
        real(dp), intent(in) :: matrix(:, :), rhs(:)
        real(dp), intent(out) :: solution(:)
        integer, intent(out) :: status
        real(dp), allocatable :: work(:, :), work_rhs(:), row_buffer(:)
        real(dp) :: scale, pivot_value, factor
        integer :: n, row, column, pivot, pivot_row

        n = size(rhs)
        status = 1
        solution = 0.0_dp
        if (size(matrix, 1) /= n .or. size(matrix, 2) /= n) return
        if (n < 1) return
        allocate(work(n, n), work_rhs(n), row_buffer(n))
        work = matrix
        work_rhs = rhs
        scale = max(1.0_dp, maxval(abs(work)))
        do column = 1, n
            pivot_row = column
            do row = column + 1, n
                if (abs(work(row, column)) > abs(work(pivot_row, column))) &
                    pivot_row = row
            end do
            pivot_value = work(pivot_row, column)
            if (abs(pivot_value) <= 100.0_dp*epsilon(1.0_dp)*scale) return
            if (pivot_row /= column) then
                row_buffer = work(column, :)
                work(column, :) = work(pivot_row, :)
                work(pivot_row, :) = row_buffer
                factor = work_rhs(column)
                work_rhs(column) = work_rhs(pivot_row)
                work_rhs(pivot_row) = factor
            end if
            do row = column + 1, n
                factor = work(row, column)/work(column, column)
                work(row, column:n) = work(row, column:n) - &
                    factor*work(column, column:n)
                work_rhs(row) = work_rhs(row) - factor*work_rhs(column)
            end do
        end do
        do row = n, 1, -1
            solution(row) = work_rhs(row)
            do pivot = row + 1, n
                solution(row) = solution(row) - work(row, pivot)*solution(pivot)
            end do
            solution(row) = solution(row)/work(row, row)
        end do
        status = 0
    end subroutine solve_dense_system

end module fortfem_tree_cotree_iga_parity
