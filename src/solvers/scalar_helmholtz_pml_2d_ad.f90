module fortfem_scalar_helmholtz_pml_2d_ad
    !! Analytical products for the complex P1 scalar Helmholtz PML system.
    use fortfem_cartesian_helmholtz_pml, only: &
        cartesian_scalar_helmholtz_pml_coefficients, &
        cartesian_scalar_helmholtz_pml_coefficients_jvp, &
        cartesian_scalar_helmholtz_pml_coefficients_vjp
    use fortfem_kinds, only: dp
    use fortfem_scalar_helmholtz_pml_2d, only: &
        solve_scalar_helmholtz_pml_p1_2d
    use fortfem_sparse_direct, only: &
        sparse_direct_factor_adjoint_csc, sparse_direct_factor_csc, &
        sparse_direct_factor_t, sparse_direct_free, &
        sparse_direct_solve_factored_jvp, sparse_direct_solve_factored_vjp, &
        sparse_direct_solve_factored
    use fortsparse, only: csc_from_triplet, csc_z_t, fortsparse_status_t
    implicit none
    private

    public :: solve_scalar_helmholtz_pml_p1_2d_jvp
    public :: solve_scalar_helmholtz_pml_p1_2d_vjp

contains

    subroutine solve_scalar_helmholtz_pml_p1_2d_jvp( &
            vertices, triangles, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, vertices_dot, stretch_dot, &
            wave_number_dot, volume_load_dot, dirichlet_values_dot, &
            solution_dot, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        complex(dp), intent(in) :: stretch(:, :)
        real(dp), intent(in) :: wave_number
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        real(dp), intent(in) :: vertices_dot(:, :)
        complex(dp), intent(in) :: stretch_dot(:, :)
        real(dp), intent(in) :: wave_number_dot
        complex(dp), intent(in) :: volume_load_dot(:)
        complex(dp), intent(in) :: dirichlet_values_dot(:)
        complex(dp), allocatable, intent(out) :: solution_dot(:)
        integer, intent(out) :: status

        type(csc_z_t) :: matrix, matrix_dot
        type(fortsparse_status_t) :: sparse_status
        type(sparse_direct_factor_t) :: factor
        complex(dp), allocatable :: solution(:), rhs(:), rhs_dot(:)
        complex(dp), allocatable :: values(:), values_dot(:), values_pre(:)
        integer, allocatable :: rows(:), columns(:)
        integer :: vertex_count

        status = 1
        if (allocated(solution_dot)) deallocate(solution_dot)
        vertex_count = size(vertices, 2)
        if (.not. valid_inputs_2d( &
            vertices, triangles, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, vertex_count)) return
        if (any(shape(vertices_dot) /= shape(vertices))) return
        if (any(shape(stretch_dot) /= shape(stretch))) return
        if (size(volume_load_dot) /= vertex_count) return
        if (size(dirichlet_values_dot) /= size(dirichlet_values)) return

        call solve_scalar_helmholtz_pml_p1_2d( &
            vertices, triangles, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, solution, status)
        if (status /= 0) return
        call assemble_2d_raw( &
            vertices, triangles, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, vertices_dot, stretch_dot, &
            wave_number_dot, volume_load_dot, dirichlet_values_dot, rows, &
            columns, values, values_dot, values_pre, rhs, rhs_dot, status)
        if (status /= 0) return
        call csc_from_triplet( &
            vertex_count, vertex_count, rows, columns, values, matrix, &
            sparse_status)
        if (sparse_status%code /= 0) then
            status = sparse_status%code
            return
        end if
        call csc_from_triplet( &
            vertex_count, vertex_count, rows, columns, values_dot, matrix_dot, &
            sparse_status)
        if (sparse_status%code /= 0) then
            status = sparse_status%code
            return
        end if
        allocate(solution_dot(vertex_count))
        call sparse_direct_factor_csc( &
            factor, vertex_count, matrix%col_ptr, matrix%row_idx, matrix%val, &
            status)
        if (status /= 0) return
        call sparse_direct_solve_factored_jvp( &
            factor, vertex_count, matrix_dot%col_ptr, matrix_dot%row_idx, &
            matrix_dot%val, solution, rhs_dot, solution_dot, status)
        call sparse_direct_free(factor)
    end subroutine solve_scalar_helmholtz_pml_p1_2d_jvp

    subroutine solve_scalar_helmholtz_pml_p1_2d_vjp( &
            vertices, triangles, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, solution, solution_bar, &
            vertices_bar, stretch_bar, wave_number_bar, volume_load_bar, &
            dirichlet_values_bar, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        complex(dp), intent(in) :: stretch(:, :)
        real(dp), intent(in) :: wave_number
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        complex(dp), intent(in) :: solution(:), solution_bar(:)
        real(dp), intent(out) :: vertices_bar(:, :)
        complex(dp), intent(out) :: stretch_bar(:, :)
        real(dp), intent(out) :: wave_number_bar
        complex(dp), intent(out) :: volume_load_bar(:)
        complex(dp), intent(out) :: dirichlet_values_bar(:)
        integer, intent(out) :: status

        type(csc_z_t) :: matrix
        type(fortsparse_status_t) :: sparse_status
        type(sparse_direct_factor_t) :: adjoint_factor
        complex(dp), allocatable :: rhs(:), rhs_dot(:), rhs_bar(:)
        complex(dp), allocatable :: values(:), values_dot(:), values_pre(:)
        complex(dp), allocatable :: matrix_values_bar(:)
        integer, allocatable :: rows(:), columns(:)
        real(dp), allocatable :: vertices_dot(:,:)
        complex(dp), allocatable :: stretch_dot(:,:), volume_load_dot(:)
        complex(dp), allocatable :: dirichlet_values_dot(:)
        logical, allocatable :: constrained(:)
        integer, allocatable :: constraint_index(:)
        integer :: vertex_count, element, local_column, local_row, entry
        integer :: row, column, local_nodes(3), constraint
        complex(dp) :: local_matrix_bar(3, 3), local_value_bar
        real(dp) :: local_vertices_bar(2, 3)
        complex(dp) :: local_stretch_bar(2)

        vertices_bar = 0.0_dp
        stretch_bar = cmplx(0.0_dp, 0.0_dp, dp)
        wave_number_bar = 0.0_dp
        volume_load_bar = cmplx(0.0_dp, 0.0_dp, dp)
        dirichlet_values_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        vertex_count = size(vertices, 2)
        if (.not. valid_inputs_2d( &
            vertices, triangles, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, vertex_count)) return
        if (any(shape(vertices_bar) /= shape(vertices))) return
        if (any(shape(stretch_bar) /= shape(stretch))) return
        if (size(solution) /= vertex_count .or. &
            size(solution_bar) /= vertex_count) return
        if (size(volume_load_bar) /= vertex_count .or. &
            size(dirichlet_values_bar) /= size(dirichlet_values)) return

        allocate(vertices_dot(2, vertex_count), stretch_dot(2, size(triangles, 2)))
        allocate(volume_load_dot(vertex_count), &
            dirichlet_values_dot(size(dirichlet_values)))
        vertices_dot = 0.0_dp
        stretch_dot = cmplx(0.0_dp, 0.0_dp, dp)
        volume_load_dot = cmplx(0.0_dp, 0.0_dp, dp)
        dirichlet_values_dot = cmplx(0.0_dp, 0.0_dp, dp)
        call assemble_2d_raw( &
            vertices, triangles, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, vertices_dot, stretch_dot, &
            0.0_dp, volume_load_dot, dirichlet_values_dot, rows, columns, &
            values, values_dot, values_pre, rhs, rhs_dot, status)
        if (status /= 0) return
        call csc_from_triplet( &
            vertex_count, vertex_count, rows, columns, values, matrix, &
            sparse_status)
        if (sparse_status%code /= 0) then
            status = sparse_status%code
            return
        end if
        allocate(rhs_bar(vertex_count), matrix_values_bar(matrix%nnz))
        call sparse_direct_factor_adjoint_csc( &
            adjoint_factor, vertex_count, matrix%col_ptr, matrix%row_idx, &
            matrix%val, status)
        if (status /= 0) return
        call sparse_direct_solve_factored_vjp( &
            adjoint_factor, vertex_count, matrix%col_ptr, matrix%row_idx, &
            solution, solution_bar, rhs_bar, matrix_values_bar, status)
        call sparse_direct_free(adjoint_factor)
        if (status /= 0) return

        allocate(constrained(vertex_count), constraint_index(vertex_count))
        call build_constraint_map( &
            vertex_count, dirichlet_nodes, constrained, constraint_index, status)
        if (status /= 0) return
        do row = 1, vertex_count
            if (constrained(row)) cycle
            volume_load_bar(row) = rhs_bar(row)
        end do
        do constraint = 1, size(dirichlet_nodes)
            dirichlet_values_bar(constraint) = &
                rhs_bar(dirichlet_nodes(constraint))
        end do

        do element = 1, size(triangles, 2)
            local_nodes = triangles(:, element)
            local_matrix_bar = cmplx(0.0_dp, 0.0_dp, dp)
            do local_column = 1, 3
                column = local_nodes(local_column)
                do local_row = 1, 3
                    row = local_nodes(local_row)
                    entry = (element - 1)*9 + (local_column - 1)*3 + local_row
                    local_value_bar = cmplx(0.0_dp, 0.0_dp, dp)
                    if (.not. constrained(row)) then
                        if (constrained(column)) then
                            constraint = constraint_index(column)
                            local_value_bar = -rhs_bar(row)*conjg( &
                                dirichlet_values(constraint))
                            dirichlet_values_bar(constraint) = &
                                dirichlet_values_bar(constraint) - &
                                rhs_bar(row)*conjg(values_pre(entry))
                        else
                            local_value_bar = csc_value_bar( &
                                matrix, matrix_values_bar, row, column)
                        end if
                    end if
                    local_matrix_bar(local_row, local_column) = &
                        local_value_bar
                end do
            end do
            local_vertices_bar = 0.0_dp
            local_stretch_bar = cmplx(0.0_dp, 0.0_dp, dp)
            call accumulate_triangle_vjp( &
                vertices, local_nodes, stretch(:, element), wave_number, &
                local_matrix_bar, local_vertices_bar, local_stretch_bar, &
                wave_number_bar, status)
            if (status /= 0) return
            do local_row = 1, 3
                vertices_bar(:, local_nodes(local_row)) = &
                    vertices_bar(:, local_nodes(local_row)) + &
                    local_vertices_bar(:, local_row)
            end do
            stretch_bar(:, element) = stretch_bar(:, element) + &
                local_stretch_bar
        end do
        status = 0
    end subroutine solve_scalar_helmholtz_pml_p1_2d_vjp

    subroutine assemble_2d_raw( &
            vertices, triangles, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, vertices_dot, stretch_dot, &
            wave_number_dot, volume_load_dot, dirichlet_values_dot, rows, &
            columns, values, values_dot, values_pre, rhs, rhs_dot, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        complex(dp), intent(in) :: stretch(:, :)
        real(dp), intent(in) :: wave_number
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        real(dp), intent(in) :: vertices_dot(:, :)
        complex(dp), intent(in) :: stretch_dot(:, :)
        real(dp), intent(in) :: wave_number_dot
        complex(dp), intent(in) :: volume_load_dot(:)
        complex(dp), intent(in) :: dirichlet_values_dot(:)
        integer, allocatable, intent(out) :: rows(:), columns(:)
        complex(dp), allocatable, intent(out) :: values(:), values_dot(:)
        complex(dp), allocatable, intent(out) :: values_pre(:)
        complex(dp), allocatable, intent(out) :: rhs(:), rhs_dot(:)
        integer, intent(out) :: status

        logical, allocatable :: constrained(:)
        integer, allocatable :: constraint_index(:)
        complex(dp) :: gradient_coefficient(3), gradient_dot(3)
        complex(dp) :: mass_coefficient, mass_dot, stretch_3d(3)
        complex(dp) :: stretch_3d_dot(3), element_matrix(3, 3)
        complex(dp) :: element_matrix_dot(3, 3)
        real(dp) :: area, area_dot, determinant, determinant_dot
        real(dp) :: gradients(2, 3), gradients_dot(2, 3)
        real(dp) :: local_mass, local_mass_dot
        integer :: element, entry, local_column, local_row, node, constraint
        integer :: vertex_count, local_nodes(3), local_status

        status = 1
        if (.not. valid_inputs_2d( &
            vertices, triangles, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, size(vertices, 2))) return
        if (any(shape(vertices_dot) /= shape(vertices))) return
        if (any(shape(stretch_dot) /= shape(stretch))) return
        if (size(volume_load_dot) /= size(volume_load)) return
        if (size(dirichlet_values_dot) /= size(dirichlet_values)) return
        vertex_count = size(vertices, 2)
        call build_constraint_map( &
            vertex_count, dirichlet_nodes, constrained, constraint_index, status)
        if (status /= 0) return
        entry = 9*size(triangles, 2) + size(dirichlet_nodes)
        allocate(rows(entry), columns(entry), values(entry), values_dot(entry), &
            values_pre(entry), rhs(vertex_count), rhs_dot(vertex_count))
        values = cmplx(0.0_dp, 0.0_dp, dp)
        values_dot = cmplx(0.0_dp, 0.0_dp, dp)
        values_pre = cmplx(0.0_dp, 0.0_dp, dp)
        rhs = volume_load
        rhs_dot = volume_load_dot
        entry = 0
        do element = 1, size(triangles, 2)
            local_nodes = triangles(:, element)
            call triangle_geometry_jvp( &
                vertices(:, local_nodes), vertices_dot(:, local_nodes), area, &
                area_dot, determinant, determinant_dot, gradients, &
                gradients_dot, local_status)
            if (local_status /= 0) return
            stretch_3d = [stretch(:, element), cmplx(1.0_dp, 0.0_dp, dp)]
            stretch_3d_dot = [ &
                stretch_dot(:, element), cmplx(0.0_dp, 0.0_dp, dp)]
            call cartesian_scalar_helmholtz_pml_coefficients( &
                stretch_3d, gradient_coefficient, mass_coefficient, &
                local_status)
            if (local_status /= 0) return
            call cartesian_scalar_helmholtz_pml_coefficients_jvp( &
                stretch_3d, stretch_3d_dot, gradient_dot, mass_dot, &
                local_status)
            if (local_status /= 0) return
            do local_column = 1, 3
                do local_row = 1, 3
                    local_mass = area/12.0_dp
                    local_mass_dot = area_dot/12.0_dp
                    if (local_row == local_column) then
                        local_mass = area/6.0_dp
                        local_mass_dot = area_dot/6.0_dp
                    end if
                    element_matrix(local_row, local_column) = area*( &
                        gradient_coefficient(1)*gradients(1, local_row)* &
                        gradients(1, local_column) + &
                        gradient_coefficient(2)*gradients(2, local_row)* &
                        gradients(2, local_column)) - &
                        wave_number**2*mass_coefficient*local_mass
                    element_matrix_dot(local_row, local_column) = &
                        area_dot*(gradient_coefficient(1)* &
                        gradients(1, local_row)*gradients(1, local_column) + &
                        gradient_coefficient(2)*gradients(2, local_row)* &
                        gradients(2, local_column)) + area*( &
                        gradient_dot(1)*gradients(1, local_row)* &
                        gradients(1, local_column) + &
                        gradient_coefficient(1)*(gradients_dot(1, local_row)* &
                        gradients(1, local_column) + gradients(1, local_row)* &
                        gradients_dot(1, local_column)) + &
                        gradient_dot(2)*gradients(2, local_row)* &
                        gradients(2, local_column) + &
                        gradient_coefficient(2)*(gradients_dot(2, local_row)* &
                        gradients(2, local_column) + gradients(2, local_row)* &
                        gradients_dot(2, local_column))) - &
                        (2.0_dp*wave_number*wave_number_dot*mass_coefficient* &
                        local_mass + wave_number**2*mass_dot*local_mass + &
                        wave_number**2*mass_coefficient*local_mass_dot)
                    entry = entry + 1
                    node = local_nodes(local_row)
                    rows(entry) = node
                    columns(entry) = local_nodes(local_column)
                    values_pre(entry) = element_matrix(local_row, local_column)
                    values(entry) = values_pre(entry)
                    values_dot(entry) = element_matrix_dot(local_row, local_column)
                    if (constrained(local_nodes(local_column))) then
                        constraint = constraint_index(local_nodes(local_column))
                        rhs(node) = rhs(node) - values_pre(entry)* &
                            dirichlet_values(constraint)
                        rhs_dot(node) = rhs_dot(node) - &
                            values_dot(entry)*dirichlet_values(constraint) - &
                            values_pre(entry)*dirichlet_values_dot(constraint)
                    end if
                    if (constrained(node) .or. &
                        constrained(local_nodes(local_column))) then
                        values(entry) = cmplx(0.0_dp, 0.0_dp, dp)
                        values_dot(entry) = cmplx(0.0_dp, 0.0_dp, dp)
                    end if
                end do
            end do
        end do
        do constraint = 1, size(dirichlet_nodes)
            entry = entry + 1
            rows(entry) = dirichlet_nodes(constraint)
            columns(entry) = dirichlet_nodes(constraint)
            values(entry) = cmplx(1.0_dp, 0.0_dp, dp)
            values_dot(entry) = cmplx(0.0_dp, 0.0_dp, dp)
            rhs(dirichlet_nodes(constraint)) = dirichlet_values(constraint)
            rhs_dot(dirichlet_nodes(constraint)) = &
                dirichlet_values_dot(constraint)
        end do
        status = 0
    end subroutine assemble_2d_raw

    subroutine accumulate_triangle_vjp( &
            vertices, local_nodes, stretch, wave_number, matrix_bar, &
            vertices_bar, stretch_bar, wave_number_bar, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: local_nodes(3)
        complex(dp), intent(in) :: stretch(2)
        real(dp), intent(in) :: wave_number
        complex(dp), intent(in) :: matrix_bar(3, 3)
        real(dp), intent(inout) :: vertices_bar(:, :)
        complex(dp), intent(inout) :: stretch_bar(2)
        real(dp), intent(inout) :: wave_number_bar
        integer, intent(out) :: status

        complex(dp) :: gradient_coefficient(3), gradient_bar(3)
        complex(dp) :: mass_coefficient, mass_bar, stretch_3d(3)
        complex(dp) :: stretch_3d_bar(3)
        real(dp) :: area, determinant, gradients(2, 3), area_bar
        real(dp) :: gradient_bar_real(2, 3), darea, ddet, dgrad(2, 3)
        complex(dp) :: dmatrix_area
        real(dp) :: direction(2, 3), direction_area, direction_det
        real(dp) :: direction_grad(2, 3)
        integer :: local_column, local_row, coordinate, node, local_status

        direction = 0.0_dp
        call triangle_geometry_jvp( &
            vertices(:, local_nodes), direction, area, darea, determinant, &
            ddet, gradients, dgrad, local_status)
        if (local_status /= 0) then
            status = local_status
            return
        end if
        stretch_3d = [stretch, cmplx(1.0_dp, 0.0_dp, dp)]
        call cartesian_scalar_helmholtz_pml_coefficients( &
            stretch_3d, gradient_coefficient, mass_coefficient, local_status)
        if (local_status /= 0) then
            status = local_status
            return
        end if
        gradient_bar = cmplx(0.0_dp, 0.0_dp, dp)
        mass_bar = cmplx(0.0_dp, 0.0_dp, dp)
        gradient_bar_real = 0.0_dp
        area_bar = 0.0_dp
        do local_column = 1, 3
            do local_row = 1, 3
                dmatrix_area = gradient_coefficient(1)* &
                    gradients(1, local_row)*gradients(1, local_column) + &
                    gradient_coefficient(2)*gradients(2, local_row)* &
                    gradients(2, local_column) - wave_number**2* &
                    mass_coefficient*merge(1.0_dp/6.0_dp, 1.0_dp/12.0_dp, &
                    local_row == local_column)
                area_bar = area_bar + real(conjg(matrix_bar(local_row, &
                    local_column))*dmatrix_area, dp)
                mass_bar = mass_bar + matrix_bar(local_row, local_column)* &
                    (-wave_number**2*area*merge(1.0_dp/6.0_dp, &
                    1.0_dp/12.0_dp, local_row == local_column))
                gradient_bar(1) = gradient_bar(1) + &
                    matrix_bar(local_row, local_column)*area* &
                    gradients(1, local_row)*gradients(1, local_column)
                gradient_bar(2) = gradient_bar(2) + &
                    matrix_bar(local_row, local_column)*area* &
                    gradients(2, local_row)*gradients(2, local_column)
                gradient_bar_real(1, local_row) = &
                    gradient_bar_real(1, local_row) + real(conjg( &
                    matrix_bar(local_row, local_column))*area* &
                    gradient_coefficient(1)*gradients(1, local_column), dp)
                gradient_bar_real(1, local_column) = &
                    gradient_bar_real(1, local_column) + real(conjg( &
                    matrix_bar(local_row, local_column))*area* &
                    gradient_coefficient(1)*gradients(1, local_row), dp)
                gradient_bar_real(2, local_row) = &
                    gradient_bar_real(2, local_row) + real(conjg( &
                    matrix_bar(local_row, local_column))*area* &
                    gradient_coefficient(2)*gradients(2, local_column), dp)
                gradient_bar_real(2, local_column) = &
                    gradient_bar_real(2, local_column) + real(conjg( &
                    matrix_bar(local_row, local_column))*area* &
                    gradient_coefficient(2)*gradients(2, local_row), dp)
                wave_number_bar = wave_number_bar + real(conjg( &
                    matrix_bar(local_row, local_column))* &
                    (-2.0_dp*wave_number*mass_coefficient*area* &
                    merge(1.0_dp/6.0_dp, 1.0_dp/12.0_dp, &
                    local_row == local_column)), dp)
            end do
        end do
        call cartesian_scalar_helmholtz_pml_coefficients_vjp( &
            stretch_3d, gradient_bar, mass_bar, stretch_3d_bar, local_status)
        if (local_status /= 0) then
            status = local_status
            return
        end if
        stretch_bar = stretch_bar + stretch_3d_bar(1:2)

        do local_row = 1, 3
            do coordinate = 1, 2
                direction = 0.0_dp
                direction(coordinate, local_row) = 1.0_dp
                call triangle_geometry_jvp( &
                    vertices(:, local_nodes), direction, area, direction_area, &
                    determinant, direction_det, gradients, direction_grad, &
                    local_status)
                if (local_status /= 0) then
                    status = local_status
                    return
                end if
                darea = area_bar*direction_area + sum( &
                    gradient_bar_real*direction_grad)
                node = local_row
                vertices_bar(coordinate, node) = &
                    vertices_bar(coordinate, node) + darea
            end do
        end do
        status = 0
    end subroutine accumulate_triangle_vjp

    subroutine triangle_geometry_jvp( &
            points, points_dot, area, area_dot, determinant, determinant_dot, &
            gradients, gradients_dot, status)
        real(dp), intent(in) :: points(2, 3), points_dot(2, 3)
        real(dp), intent(out) :: area, area_dot, determinant, determinant_dot
        real(dp), intent(out) :: gradients(2, 3), gradients_dot(2, 3)
        integer, intent(out) :: status

        real(dp) :: numerator(2, 3), numerator_dot(2, 3)
        real(dp) :: a(2), b(2), a_dot(2), b_dot(2)
        integer :: node

        a = points(:, 2) - points(:, 1)
        b = points(:, 3) - points(:, 1)
        a_dot = points_dot(:, 2) - points_dot(:, 1)
        b_dot = points_dot(:, 3) - points_dot(:, 1)
        determinant = a(1)*b(2) - b(1)*a(2)
        determinant_dot = a_dot(1)*b(2) + a(1)*b_dot(2) - &
            b_dot(1)*a(2) - b(1)*a_dot(2)
        status = 1
        if (abs(determinant) <= tiny(1.0_dp)) return
        area = 0.5_dp*abs(determinant)
        area_dot = 0.5_dp*sign(1.0_dp, determinant)*determinant_dot
        numerator(:, 1) = [points(2, 2) - points(2, 3), &
            points(1, 3) - points(1, 2)]
        numerator(:, 2) = [points(2, 3) - points(2, 1), &
            points(1, 1) - points(1, 3)]
        numerator(:, 3) = [points(2, 1) - points(2, 2), &
            points(1, 2) - points(1, 1)]
        numerator_dot(:, 1) = [points_dot(2, 2) - points_dot(2, 3), &
            points_dot(1, 3) - points_dot(1, 2)]
        numerator_dot(:, 2) = [points_dot(2, 3) - points_dot(2, 1), &
            points_dot(1, 1) - points_dot(1, 3)]
        numerator_dot(:, 3) = [points_dot(2, 1) - points_dot(2, 2), &
            points_dot(1, 2) - points_dot(1, 1)]
        do node = 1, 3
            gradients(:, node) = numerator(:, node)/determinant
            gradients_dot(:, node) = numerator_dot(:, node)/determinant - &
                numerator(:, node)*determinant_dot/determinant**2
        end do
        status = 0
    end subroutine triangle_geometry_jvp

    subroutine build_constraint_map( &
            vertex_count, dirichlet_nodes, constrained, constraint_index, status)
        integer, intent(in) :: vertex_count, dirichlet_nodes(:)
        logical, allocatable, intent(out) :: constrained(:)
        integer, allocatable, intent(out) :: constraint_index(:)
        integer, intent(out) :: status

        integer :: constraint, node

        allocate(constrained(vertex_count), constraint_index(vertex_count))
        constrained = .false.
        constraint_index = 0
        status = 1
        do constraint = 1, size(dirichlet_nodes)
            node = dirichlet_nodes(constraint)
            if (node < 1 .or. node > vertex_count .or. constrained(node)) return
            constrained(node) = .true.
            constraint_index(node) = constraint
        end do
        status = 0
    end subroutine build_constraint_map

    function csc_value_bar(matrix, values_bar, row, column) result(value)
        type(csc_z_t), intent(in) :: matrix
        complex(dp), intent(in) :: values_bar(:)
        integer, intent(in) :: row, column
        complex(dp) :: value
        integer :: entry

        value = cmplx(0.0_dp, 0.0_dp, dp)
        do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
            if (matrix%row_idx(entry) == row) then
                value = values_bar(entry)
                return
            end if
        end do
    end function csc_value_bar

    pure logical function valid_inputs_2d( &
            vertices, triangles, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, solution_size) result(valid)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), dirichlet_nodes(:)
        complex(dp), intent(in) :: stretch(:, :), volume_load(:)
        real(dp), intent(in) :: wave_number
        complex(dp), intent(in) :: dirichlet_values(:)
        integer, intent(in) :: solution_size

        valid = size(vertices, 1) == 2 .and. size(vertices, 2) >= 3
        valid = valid .and. size(triangles, 1) == 3
        valid = valid .and. size(triangles, 2) >= 1
        valid = valid .and. size(stretch, 1) == 2
        valid = valid .and. size(stretch, 2) == size(triangles, 2)
        valid = valid .and. size(volume_load) == solution_size
        valid = valid .and. size(dirichlet_nodes) == size(dirichlet_values)
        valid = valid .and. solution_size == size(vertices, 2)
        valid = valid .and. wave_number > 0.0_dp
        valid = valid .and. all(triangles >= 1) .and. &
            all(triangles <= solution_size)
        valid = valid .and. all(dirichlet_nodes >= 1) .and. &
            all(dirichlet_nodes <= solution_size)
    end function valid_inputs_2d

end module fortfem_scalar_helmholtz_pml_2d_ad
