module fortfem_scalar_helmholtz_planar_dtn_2d
    use fortfem_kinds, only: dp
    use fortfem_planar_helmholtz_dtn, only: &
        assemble_planar_helmholtz_dtn_form, &
        assemble_planar_helmholtz_dtn_form_jvp, &
        assemble_planar_helmholtz_dtn_form_vjp
    use fortfem_sparse_direct, only: sparse_direct_factor_adjoint_csc, &
        sparse_direct_factor_csc, sparse_direct_factor_t, sparse_direct_free, &
        sparse_direct_solve_factored, sparse_direct_solve_factored_jvp, &
        sparse_direct_solve_factored_vjp
    use fortsparse, only: csc_from_triplet, csc_z_t, fortsparse_status_t
    implicit none
    private

    public :: solve_scalar_helmholtz_planar_dtn_p1
    public :: solve_scalar_helmholtz_planar_dtn_p1_jvp
    public :: solve_scalar_helmholtz_planar_dtn_p1_vjp

contains

    subroutine solve_scalar_helmholtz_planar_dtn_p1( &
            vertices, triangles, dtn_nodes, wavenumber, period, volume_load, &
            dirichlet_nodes, dirichlet_values, solution, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), dtn_nodes(:)
        real(dp), intent(in) :: wavenumber, period
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        complex(dp), intent(out) :: solution(:)
        integer, intent(out) :: status

        type(csc_z_t) :: sparse_matrix
        type(fortsparse_status_t) :: sparse_status
        complex(dp), allocatable :: dtn_form(:, :), matrix(:, :), rhs(:)
        integer :: vertex_count

        solution = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        vertex_count = size(vertices, 2)
        if (.not. valid_solver_inputs( &
            vertices, triangles, dtn_nodes, wavenumber, period, volume_load, &
            dirichlet_nodes, dirichlet_values, size(solution))) return

        allocate(matrix(vertex_count, vertex_count), rhs(vertex_count))
        call assemble_volume_matrix( &
            vertices, triangles, wavenumber, matrix, status)
        if (status /= 0) return
        allocate(dtn_form(size(dtn_nodes), size(dtn_nodes)))
        call assemble_planar_helmholtz_dtn_form( &
            size(dtn_nodes), wavenumber, period, dtn_form, status)
        if (status /= 0) return
        matrix(dtn_nodes, dtn_nodes) = &
            matrix(dtn_nodes, dtn_nodes) - dtn_form
        rhs = volume_load
        call apply_dirichlet_elimination( &
            matrix, rhs, dirichlet_nodes, dirichlet_values)
        call dense_to_csc(matrix, sparse_matrix, sparse_status)
        if (sparse_status%code /= 0) return
        call solve_csc_once(sparse_matrix, rhs, solution, status)
        if (status /= 0) return
        status = 0
    end subroutine solve_scalar_helmholtz_planar_dtn_p1

    subroutine solve_scalar_helmholtz_planar_dtn_p1_jvp( &
            vertices, triangles, dtn_nodes, wavenumber, period, volume_load, &
            dirichlet_nodes, dirichlet_values, wavenumber_dot, period_dot, &
            volume_load_dot, dirichlet_values_dot, solution_dot, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), dtn_nodes(:)
        real(dp), intent(in) :: wavenumber, period
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        real(dp), intent(in) :: wavenumber_dot, period_dot
        complex(dp), intent(in) :: volume_load_dot(:)
        complex(dp), intent(in) :: dirichlet_values_dot(:)
        complex(dp), intent(out) :: solution_dot(:)
        integer, intent(out) :: status

        type(csc_z_t) :: matrix_csc, matrix_dot_csc
        type(fortsparse_status_t) :: sparse_status
        type(sparse_direct_factor_t) :: factor
        complex(dp), allocatable :: dtn_form(:, :), dtn_form_dot(:, :)
        complex(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
        complex(dp), allocatable :: rhs(:), rhs_dot(:), solution(:)
        integer :: vertex_count

        solution_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        vertex_count = size(vertices, 2)
        if (.not. valid_solver_inputs( &
            vertices, triangles, dtn_nodes, wavenumber, period, volume_load, &
            dirichlet_nodes, dirichlet_values, size(solution_dot))) return
        if (size(volume_load_dot) /= vertex_count) return
        if (size(dirichlet_values_dot) /= size(dirichlet_values)) return

        allocate ( &
            matrix(vertex_count, vertex_count), &
            matrix_dot(vertex_count, vertex_count), rhs(vertex_count), &
            rhs_dot(vertex_count), solution(vertex_count))
        call assemble_volume_matrix( &
            vertices, triangles, wavenumber, matrix, status)
        if (status /= 0) return
        call assemble_volume_matrix_jvp( &
            vertices, triangles, wavenumber, wavenumber_dot, matrix_dot, status)
        if (status /= 0) return
        allocate ( &
            dtn_form(size(dtn_nodes), size(dtn_nodes)), &
            dtn_form_dot(size(dtn_nodes), size(dtn_nodes)))
        call assemble_planar_helmholtz_dtn_form( &
            size(dtn_nodes), wavenumber, period, dtn_form, status)
        if (status /= 0) return
        call assemble_planar_helmholtz_dtn_form_jvp( &
            size(dtn_nodes), wavenumber, period, wavenumber_dot, period_dot, &
            dtn_form_dot, status)
        if (status /= 0) return
        matrix(dtn_nodes, dtn_nodes) = &
            matrix(dtn_nodes, dtn_nodes) - dtn_form
        matrix_dot(dtn_nodes, dtn_nodes) = &
            matrix_dot(dtn_nodes, dtn_nodes) - dtn_form_dot
        rhs = volume_load
        rhs_dot = volume_load_dot
        call apply_dirichlet_elimination_jvp( &
            matrix, matrix_dot, rhs, rhs_dot, dirichlet_nodes, &
            dirichlet_values, dirichlet_values_dot)
        call dense_to_csc(matrix, matrix_csc, sparse_status)
        if (sparse_status%code /= 0) return
        call dense_to_csc(matrix_dot, matrix_dot_csc, sparse_status)
        if (sparse_status%code /= 0) return
        call sparse_direct_factor_csc( &
            factor, vertex_count, matrix_csc%col_ptr, matrix_csc%row_idx, &
            matrix_csc%val, status)
        if (status /= 0) return
        call sparse_direct_solve_factored(factor, rhs, solution, status)
        if (status /= 0) then
            call sparse_direct_free(factor)
            return
        end if
        call sparse_direct_solve_factored_jvp( &
            factor, vertex_count, matrix_dot_csc%col_ptr, &
            matrix_dot_csc%row_idx, matrix_dot_csc%val, solution, rhs_dot, &
            solution_dot, status)
        call sparse_direct_free(factor)
    end subroutine solve_scalar_helmholtz_planar_dtn_p1_jvp

    subroutine solve_scalar_helmholtz_planar_dtn_p1_vjp( &
            vertices, triangles, dtn_nodes, wavenumber, period, volume_load, &
            dirichlet_nodes, dirichlet_values, solution, solution_bar, &
            volume_load_bar, dirichlet_values_bar, wavenumber_bar, period_bar, &
            status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), dtn_nodes(:)
        real(dp), intent(in) :: wavenumber, period
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        complex(dp), intent(in) :: solution(:), solution_bar(:)
        complex(dp), intent(out) :: volume_load_bar(:)
        complex(dp), intent(out) :: dirichlet_values_bar(:)
        real(dp), intent(out) :: wavenumber_bar, period_bar
        integer, intent(out) :: status

        type(csc_z_t) :: matrix_csc
        type(fortsparse_status_t) :: sparse_status
        type(sparse_direct_factor_t) :: adjoint_factor
        complex(dp), allocatable :: dtn_form(:, :), dtn_form_bar(:, :)
        complex(dp), allocatable :: matrix(:, :), matrix_bar(:, :)
        complex(dp), allocatable :: matrix_pre(:, :), matrix_pre_bar(:, :)
        complex(dp), allocatable :: matrix_values_bar(:), rhs(:), rhs_bar(:)
        real(dp) :: volume_wavenumber_bar
        integer :: vertex_count

        volume_load_bar = cmplx(0.0_dp, 0.0_dp, dp)
        dirichlet_values_bar = cmplx(0.0_dp, 0.0_dp, dp)
        wavenumber_bar = 0.0_dp
        period_bar = 0.0_dp
        status = 1
        vertex_count = size(vertices, 2)
        if (.not. valid_solver_inputs( &
            vertices, triangles, dtn_nodes, wavenumber, period, volume_load, &
            dirichlet_nodes, dirichlet_values, size(solution))) return
        if (size(solution_bar) /= vertex_count) return
        if (size(volume_load_bar) /= vertex_count) return
        if (size(dirichlet_values_bar) /= size(dirichlet_values)) return

        allocate ( &
            matrix_pre(vertex_count, vertex_count), &
            matrix(vertex_count, vertex_count), rhs(vertex_count), &
            rhs_bar(vertex_count))
        call assemble_volume_matrix( &
            vertices, triangles, wavenumber, matrix_pre, status)
        if (status /= 0) return
        allocate (dtn_form(size(dtn_nodes), size(dtn_nodes)))
        call assemble_planar_helmholtz_dtn_form( &
            size(dtn_nodes), wavenumber, period, dtn_form, status)
        if (status /= 0) return
        matrix_pre(dtn_nodes, dtn_nodes) = &
            matrix_pre(dtn_nodes, dtn_nodes) - dtn_form
        matrix = matrix_pre
        rhs = volume_load
        call apply_dirichlet_elimination( &
            matrix, rhs, dirichlet_nodes, dirichlet_values)
        call dense_to_csc( &
            matrix, matrix_csc, sparse_status, retain_zeros=.true.)
        if (sparse_status%code /= 0) return
        allocate (matrix_values_bar(matrix_csc%nnz))
        call sparse_direct_factor_adjoint_csc( &
            adjoint_factor, vertex_count, matrix_csc%col_ptr, &
            matrix_csc%row_idx, matrix_csc%val, status)
        if (status /= 0) return
        call sparse_direct_solve_factored_vjp( &
            adjoint_factor, vertex_count, matrix_csc%col_ptr, &
            matrix_csc%row_idx, solution, solution_bar, rhs_bar, &
            matrix_values_bar, status)
        call sparse_direct_free(adjoint_factor)
        if (status /= 0) return
        call csc_values_to_dense( &
            matrix_csc, matrix_values_bar, matrix_bar)
        call reverse_dirichlet_elimination( &
            matrix_pre, matrix_bar, rhs_bar, dirichlet_nodes, &
            dirichlet_values, matrix_pre_bar, volume_load_bar, &
            dirichlet_values_bar)

        allocate (dtn_form_bar(size(dtn_nodes), size(dtn_nodes)))
        dtn_form_bar = -matrix_pre_bar(dtn_nodes, dtn_nodes)
        call assemble_planar_helmholtz_dtn_form_vjp( &
            size(dtn_nodes), wavenumber, period, dtn_form_bar, &
            wavenumber_bar, period_bar, status)
        if (status /= 0) return
        call volume_wavenumber_vjp( &
            vertices, triangles, wavenumber, matrix_pre_bar, &
            volume_wavenumber_bar, status)
        if (status /= 0) return
        wavenumber_bar = wavenumber_bar + volume_wavenumber_bar
        status = 0
    end subroutine solve_scalar_helmholtz_planar_dtn_p1_vjp

    pure logical function valid_solver_inputs( &
            vertices, triangles, dtn_nodes, wavenumber, period, volume_load, &
            dirichlet_nodes, dirichlet_values, solution_size) result(valid)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), dtn_nodes(:)
        real(dp), intent(in) :: wavenumber, period
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        integer, intent(in) :: solution_size

        integer :: vertex_count

        valid = .false.
        vertex_count = size(vertices, 2)
        if (size(vertices, 1) /= 2 .or. size(triangles, 1) /= 3) return
        if (vertex_count < 1 .or. size(triangles, 2) < 1) return
        if (wavenumber <= 0.0_dp .or. period <= 0.0_dp) return
        if (size(dtn_nodes) < 3) return
        if (size(volume_load) /= vertex_count) return
        if (solution_size /= vertex_count) return
        if (size(dirichlet_nodes) /= size(dirichlet_values)) return
        if (any(triangles < 1) .or. any(triangles > vertex_count)) return
        if (any(dtn_nodes < 1) .or. any(dtn_nodes > vertex_count)) return
        if (any(dirichlet_nodes < 1) .or. &
            any(dirichlet_nodes > vertex_count)) return
        valid = .true.
    end function valid_solver_inputs

    subroutine assemble_volume_matrix( &
            vertices, triangles, wavenumber, volume_matrix, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: wavenumber
        complex(dp), intent(out) :: volume_matrix(:, :)
        integer, intent(out) :: status

        real(dp) :: area, determinant
        real(dp) :: gradients(2, 3), local_mass(3, 3)
        real(dp) :: x1, x2, x3, y1, y2, y3
        integer :: element, first, local_nodes(3), second

        volume_matrix = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        do element = 1, size(triangles, 2)
            local_nodes = triangles(:, element)
            x1 = vertices(1, local_nodes(1))
            y1 = vertices(2, local_nodes(1))
            x2 = vertices(1, local_nodes(2))
            y2 = vertices(2, local_nodes(2))
            x3 = vertices(1, local_nodes(3))
            y3 = vertices(2, local_nodes(3))
            determinant = (x2 - x1)*(y3 - y1) - &
                (x3 - x1)*(y2 - y1)
            area = 0.5_dp*abs(determinant)
            if (area <= tiny(1.0_dp)) return
            gradients(:, 1) = [y2 - y3, x3 - x2]/determinant
            gradients(:, 2) = [y3 - y1, x1 - x3]/determinant
            gradients(:, 3) = [y1 - y2, x2 - x1]/determinant
            call local_triangle_mass(area, local_mass)
            do first = 1, 3
                do second = 1, 3
                    volume_matrix(local_nodes(first), local_nodes(second)) = &
                        volume_matrix(local_nodes(first), local_nodes(second)) + &
                        area*dot_product( &
                        gradients(:, first), gradients(:, second)) - &
                        wavenumber**2*local_mass(first, second)
                end do
            end do
        end do
        status = 0
    end subroutine assemble_volume_matrix

    subroutine assemble_volume_matrix_jvp( &
            vertices, triangles, wavenumber, wavenumber_dot, &
            volume_matrix_dot, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: wavenumber, wavenumber_dot
        complex(dp), intent(out) :: volume_matrix_dot(:, :)
        integer, intent(out) :: status

        real(dp) :: area, determinant, local_mass(3, 3)
        integer :: element, first, local_nodes(3), second

        volume_matrix_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        do element = 1, size(triangles, 2)
            local_nodes = triangles(:, element)
            determinant = &
                (vertices(1, local_nodes(2)) - vertices(1, local_nodes(1)))* &
                (vertices(2, local_nodes(3)) - vertices(2, local_nodes(1))) - &
                (vertices(1, local_nodes(3)) - vertices(1, local_nodes(1)))* &
                (vertices(2, local_nodes(2)) - vertices(2, local_nodes(1)))
            area = 0.5_dp*abs(determinant)
            if (area <= tiny(1.0_dp)) return
            call local_triangle_mass(area, local_mass)
            do first = 1, 3
                do second = 1, 3
                    volume_matrix_dot( &
                        local_nodes(first), local_nodes(second)) = &
                        volume_matrix_dot( &
                        local_nodes(first), local_nodes(second)) - &
                        2.0_dp*wavenumber*wavenumber_dot* &
                        local_mass(first, second)
                end do
            end do
        end do
        status = 0
    end subroutine assemble_volume_matrix_jvp

    pure subroutine local_triangle_mass(area, mass)
        real(dp), intent(in) :: area
        real(dp), intent(out) :: mass(3, 3)

        mass = area/12.0_dp
        mass(1, 1) = area/6.0_dp
        mass(2, 2) = area/6.0_dp
        mass(3, 3) = area/6.0_dp
    end subroutine local_triangle_mass

    subroutine apply_dirichlet_elimination( &
            matrix, rhs, nodes, values)
        complex(dp), intent(inout) :: matrix(:, :), rhs(:)
        integer, intent(in) :: nodes(:)
        complex(dp), intent(in) :: values(:)

        integer :: entry, node

        do entry = 1, size(nodes)
            node = nodes(entry)
            rhs = rhs - matrix(:, node)*values(entry)
            matrix(:, node) = cmplx(0.0_dp, 0.0_dp, dp)
            matrix(node, :) = cmplx(0.0_dp, 0.0_dp, dp)
            matrix(node, node) = cmplx(1.0_dp, 0.0_dp, dp)
            rhs(node) = values(entry)
        end do
    end subroutine apply_dirichlet_elimination

    subroutine apply_dirichlet_elimination_jvp( &
            matrix, matrix_dot, rhs, rhs_dot, nodes, values, values_dot)
        complex(dp), intent(inout) :: matrix(:, :), matrix_dot(:, :)
        complex(dp), intent(inout) :: rhs(:), rhs_dot(:)
        integer, intent(in) :: nodes(:)
        complex(dp), intent(in) :: values(:), values_dot(:)

        integer :: entry, node

        do entry = 1, size(nodes)
            node = nodes(entry)
            rhs_dot = rhs_dot - matrix_dot(:, node)*values(entry) - &
                matrix(:, node)*values_dot(entry)
            rhs = rhs - matrix(:, node)*values(entry)
            matrix(:, node) = cmplx(0.0_dp, 0.0_dp, dp)
            matrix(node, :) = cmplx(0.0_dp, 0.0_dp, dp)
            matrix(node, node) = cmplx(1.0_dp, 0.0_dp, dp)
            matrix_dot(:, node) = cmplx(0.0_dp, 0.0_dp, dp)
            matrix_dot(node, :) = cmplx(0.0_dp, 0.0_dp, dp)
            rhs(node) = values(entry)
            rhs_dot(node) = values_dot(entry)
        end do
    end subroutine apply_dirichlet_elimination_jvp

    subroutine reverse_dirichlet_elimination( &
            matrix_pre, matrix_bar, rhs_bar, nodes, values, matrix_pre_bar, &
            load_bar, values_bar)
        complex(dp), intent(in) :: matrix_pre(:, :), matrix_bar(:, :)
        complex(dp), intent(in) :: rhs_bar(:)
        integer, intent(in) :: nodes(:)
        complex(dp), intent(in) :: values(:)
        complex(dp), allocatable, intent(out) :: matrix_pre_bar(:, :)
        complex(dp), intent(out) :: load_bar(:), values_bar(:)

        logical, allocatable :: free(:)
        integer :: column, entry, node, row

        allocate (free(size(rhs_bar)), source=.true.)
        free(nodes) = .false.
        allocate (matrix_pre_bar(size(matrix_pre, 1), size(matrix_pre, 2)), &
            source=cmplx(0.0_dp, 0.0_dp, dp))
        load_bar = cmplx(0.0_dp, 0.0_dp, dp)
        where (free) load_bar = rhs_bar
        values_bar = cmplx(0.0_dp, 0.0_dp, dp)
        do column = 1, size(matrix_pre, 2)
            if (.not. free(column)) cycle
            do row = 1, size(matrix_pre, 1)
                if (free(row)) matrix_pre_bar(row, column) = &
                    matrix_bar(row, column)
            end do
        end do
        do entry = 1, size(nodes)
            node = nodes(entry)
            values_bar(entry) = rhs_bar(node)
            do row = 1, size(rhs_bar)
                if (.not. free(row)) cycle
                values_bar(entry) = values_bar(entry) - &
                    conjg(matrix_pre(row, node))*rhs_bar(row)
                matrix_pre_bar(row, node) = matrix_pre_bar(row, node) - &
                    rhs_bar(row)*conjg(values(entry))
            end do
        end do
    end subroutine reverse_dirichlet_elimination

    subroutine volume_wavenumber_vjp( &
            vertices, triangles, wavenumber, matrix_bar, wavenumber_bar, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: wavenumber
        complex(dp), intent(in) :: matrix_bar(:, :)
        real(dp), intent(out) :: wavenumber_bar
        integer, intent(out) :: status

        complex(dp), allocatable :: matrix_dot(:, :)

        allocate (matrix_dot(size(matrix_bar, 1), size(matrix_bar, 2)))
        call assemble_volume_matrix_jvp( &
            vertices, triangles, wavenumber, 1.0_dp, matrix_dot, status)
        if (status /= 0) return
        wavenumber_bar = real(sum(conjg(matrix_bar)*matrix_dot), dp)
    end subroutine volume_wavenumber_vjp

    subroutine dense_to_csc(matrix, sparse_matrix, sparse_status, retain_zeros)
        complex(dp), intent(in) :: matrix(:, :)
        type(csc_z_t), intent(out) :: sparse_matrix
        type(fortsparse_status_t), intent(out) :: sparse_status
        logical, intent(in), optional :: retain_zeros

        complex(dp), allocatable :: values(:)
        integer, allocatable :: columns(:), rows(:)
        logical :: keep_zeros
        integer :: column, entry, row

        keep_zeros = .false.
        if (present(retain_zeros)) keep_zeros = retain_zeros
        allocate ( &
            rows(size(matrix)), columns(size(matrix)), values(size(matrix)))
        entry = 0
        do column = 1, size(matrix, 2)
            do row = 1, size(matrix, 1)
                if (.not. keep_zeros) then
                    if (abs(matrix(row, column)) <= tiny(1.0_dp)) cycle
                end if
                entry = entry + 1
                rows(entry) = row
                columns(entry) = column
                values(entry) = matrix(row, column)
            end do
        end do
        call csc_from_triplet( &
            size(matrix, 1), size(matrix, 2), rows(:entry), columns(:entry), &
            values(:entry), sparse_matrix, sparse_status)
    end subroutine dense_to_csc

    subroutine csc_values_to_dense(sparse_matrix, values, matrix)
        type(csc_z_t), intent(in) :: sparse_matrix
        complex(dp), intent(in) :: values(:)
        complex(dp), allocatable, intent(out) :: matrix(:, :)

        integer :: column, entry

        allocate ( &
            matrix(sparse_matrix%nrow, sparse_matrix%ncol), &
            source=cmplx(0.0_dp, 0.0_dp, dp))
        do column = 1, sparse_matrix%ncol
            do entry = sparse_matrix%col_ptr(column), &
                    sparse_matrix%col_ptr(column + 1) - 1
                matrix(sparse_matrix%row_idx(entry), column) = values(entry)
            end do
        end do
    end subroutine csc_values_to_dense

    subroutine solve_csc_once(matrix, rhs, solution, status)
        type(csc_z_t), intent(in) :: matrix
        complex(dp), intent(in) :: rhs(:)
        complex(dp), intent(out) :: solution(:)
        integer, intent(out) :: status

        type(sparse_direct_factor_t) :: factor

        call sparse_direct_factor_csc( &
            factor, matrix%nrow, matrix%col_ptr, matrix%row_idx, matrix%val, &
            status)
        if (status /= 0) return
        call sparse_direct_solve_factored(factor, rhs, solution, status)
        call sparse_direct_free(factor)
    end subroutine solve_csc_once

end module fortfem_scalar_helmholtz_planar_dtn_2d
