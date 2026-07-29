module fortfem_laplace_symmetric_coupling_2d
    use fortfem_advanced_solvers, only: solve_dense => solve, &
        solver_options, solver_options_t, solver_stats_t
    use fortfem_kinds, only: dp
    use fortfem_laplace_boundary_operators_2d, only: &
        assemble_laplace_double_layer_mixed_linear, &
        assemble_laplace_hypersingular_linear, &
        assemble_laplace_single_layer_constant
    use fortfem_helmholtz_boundary_operators_2d, only: &
        assemble_helmholtz_double_layer_mixed_linear, &
        assemble_helmholtz_hypersingular_linear, &
        assemble_helmholtz_single_layer_constant
    use fortsparse, only: csc_from_triplet, csc_z_t, fortsparse_status_t, &
        sparse_solve_once
    implicit none

    private

    public :: assemble_laplace_symmetric_coupling_p1_p0
    public :: assemble_helmholtz_symmetric_coupling_p1_p0
    public :: solve_laplace_symmetric_coupling_p1_p0
    public :: solve_helmholtz_symmetric_coupling_p1_p0

contains

    subroutine assemble_helmholtz_symmetric_coupling_p1_p0( &
            vertices, triangles, panel_start, panel_end, panel_nodes, &
            wavenumber, quadrature_order, matrix, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: panel_nodes(:, :)
        real(dp), intent(in) :: wavenumber
        integer, intent(in) :: quadrature_order
        complex(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: double_layer(:, :)
        complex(dp), allocatable :: hypersingular(:, :)
        complex(dp), allocatable :: single_layer(:, :)
        real(dp), allocatable :: mass(:, :), volume_mass(:, :)
        real(dp), allocatable :: stiffness(:, :)
        integer :: operator_status, panel_count, total_dofs, vertex_count

        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        vertex_count = size(vertices, 2)
        panel_count = size(panel_start, 2)
        total_dofs = vertex_count + panel_count
        if (.not. valid_discretization_shapes( &
            vertices, triangles, panel_start, panel_end, panel_nodes, &
            quadrature_order)) return
        if (size(matrix, 1) /= total_dofs .or. &
            size(matrix, 2) /= total_dofs) return
        if (wavenumber <= 0.0_dp) return
        if (any(triangles < 1) .or. any(triangles > vertex_count)) return
        if (any(panel_nodes < 1) .or. any(panel_nodes > vertex_count)) return
        if (.not. panels_match_trace_nodes( &
            vertices, panel_start, panel_end, panel_nodes)) return

        allocate(stiffness(vertex_count, vertex_count))
        allocate(volume_mass(vertex_count, vertex_count))
        call assemble_p1_stiffness( &
            vertices, triangles, stiffness, operator_status)
        if (operator_status /= 0) return
        call assemble_p1_mass(vertices, triangles, volume_mass, operator_status)
        if (operator_status /= 0) return

        allocate(mass(panel_count, vertex_count))
        call assemble_boundary_mass( &
            panel_start, panel_end, panel_nodes, mass, operator_status)
        if (operator_status /= 0) return
        allocate(double_layer(panel_count, vertex_count))
        call assemble_helmholtz_double_layer_mixed_linear( &
            panel_start, panel_end, panel_nodes, vertex_count, wavenumber, &
            quadrature_order, double_layer, operator_status)
        if (operator_status /= 0) return
        allocate(hypersingular(vertex_count, vertex_count))
        call assemble_helmholtz_hypersingular_linear( &
            panel_start, panel_end, panel_nodes, vertex_count, wavenumber, &
            quadrature_order, hypersingular, operator_status)
        if (operator_status /= 0) return
        allocate(single_layer(panel_count, panel_count))
        call assemble_helmholtz_single_layer_constant( &
            panel_start, panel_end, wavenumber, quadrature_order, &
            single_layer, operator_status)
        if (operator_status /= 0) return

        matrix(1:vertex_count, 1:vertex_count) = &
            stiffness - wavenumber**2*volume_mass + hypersingular
        matrix(vertex_count + 1:, vertex_count + 1:) = single_layer
        matrix(vertex_count + 1:, 1:vertex_count) = &
            0.5_dp*mass - double_layer
        matrix(1:vertex_count, vertex_count + 1:) = &
            -transpose(matrix(vertex_count + 1:, 1:vertex_count))
        status = 0
    end subroutine assemble_helmholtz_symmetric_coupling_p1_p0

    subroutine assemble_laplace_symmetric_coupling_p1_p0( &
            vertices, triangles, panel_start, panel_end, panel_nodes, &
            quadrature_order, matrix, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: panel_nodes(:, :)
        integer, intent(in) :: quadrature_order
        real(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: double_layer(:, :), hypersingular(:, :)
        real(dp), allocatable :: mass(:, :), single_layer(:, :)
        real(dp), allocatable :: stiffness(:, :)
        integer :: operator_status, panel_count, vertex_count

        matrix = 0.0_dp
        status = 1
        vertex_count = size(vertices, 2)
        panel_count = size(panel_start, 2)
        if (.not. valid_discretization_shapes( &
            vertices, triangles, panel_start, panel_end, panel_nodes, &
            quadrature_order)) return
        if (size(matrix, 1) /= vertex_count + panel_count .or. &
            size(matrix, 2) /= vertex_count + panel_count) return
        if (any(triangles < 1) .or. any(triangles > vertex_count)) return
        if (any(panel_nodes < 1) .or. any(panel_nodes > vertex_count)) return
        if (.not. panels_match_trace_nodes( &
            vertices, panel_start, panel_end, panel_nodes)) return

        allocate(stiffness(vertex_count, vertex_count))
        call assemble_p1_stiffness( &
            vertices, triangles, stiffness, operator_status)
        if (operator_status /= 0) return

        allocate(mass(panel_count, vertex_count))
        call assemble_boundary_mass( &
            panel_start, panel_end, panel_nodes, mass, operator_status)
        if (operator_status /= 0) return

        allocate(double_layer(panel_count, vertex_count))
        call assemble_laplace_double_layer_mixed_linear( &
            panel_start, panel_end, panel_nodes, vertex_count, &
            quadrature_order, double_layer, operator_status)
        if (operator_status /= 0) return

        allocate(hypersingular(vertex_count, vertex_count))
        call assemble_laplace_hypersingular_linear( &
            panel_start, panel_end, panel_nodes, vertex_count, &
            quadrature_order, hypersingular, operator_status)
        if (operator_status /= 0) return

        allocate(single_layer(panel_count, panel_count))
        call assemble_laplace_single_layer_constant( &
            panel_start, panel_end, quadrature_order, single_layer, &
            operator_status)
        if (operator_status /= 0) return

        matrix(1:vertex_count, 1:vertex_count) = stiffness + hypersingular
        matrix(vertex_count + 1:, vertex_count + 1:) = single_layer
        matrix(vertex_count + 1:, 1:vertex_count) = &
            0.5_dp * mass - double_layer
        matrix(1:vertex_count, vertex_count + 1:) = &
            -transpose(matrix(vertex_count + 1:, 1:vertex_count))
        status = 0
    end subroutine assemble_laplace_symmetric_coupling_p1_p0

    subroutine solve_laplace_symmetric_coupling_p1_p0( &
            vertices, triangles, panel_start, panel_end, panel_nodes, &
            quadrature_order, volume_load, dirichlet_jump, neumann_jump, &
            interior_solution, exterior_flux, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: panel_nodes(:, :)
        integer, intent(in) :: quadrature_order
        real(dp), intent(in) :: volume_load(:), dirichlet_jump(:)
        real(dp), intent(in) :: neumann_jump(:)
        real(dp), intent(out) :: interior_solution(:), exterior_flux(:)
        integer, intent(out) :: status

        real(dp), allocatable :: hypersingular(:, :), mass(:, :)
        real(dp), allocatable :: matrix(:, :), right_hand_side(:), solution(:)
        type(solver_options_t) :: options
        type(solver_stats_t) :: statistics
        integer :: operator_status, panel_count, total_dofs, vertex_count

        interior_solution = 0.0_dp
        exterior_flux = 0.0_dp
        status = 1
        vertex_count = size(vertices, 2)
        panel_count = size(panel_start, 2)
        total_dofs = vertex_count + panel_count
        if (size(volume_load) /= vertex_count .or. &
            size(dirichlet_jump) /= vertex_count) return
        if (size(neumann_jump) /= panel_count) return
        if (size(interior_solution) /= vertex_count .or. &
            size(exterior_flux) /= panel_count) return

        allocate(matrix(total_dofs, total_dofs))
        call assemble_laplace_symmetric_coupling_p1_p0( &
            vertices, triangles, panel_start, panel_end, panel_nodes, &
            quadrature_order, matrix, operator_status)
        if (operator_status /= 0) return

        allocate(mass(panel_count, vertex_count))
        call assemble_boundary_mass( &
            panel_start, panel_end, panel_nodes, mass, operator_status)
        if (operator_status /= 0) return
        allocate(hypersingular(vertex_count, vertex_count))
        call assemble_laplace_hypersingular_linear( &
            panel_start, panel_end, panel_nodes, vertex_count, &
            quadrature_order, hypersingular, operator_status)
        if (operator_status /= 0) return

        allocate(right_hand_side(total_dofs), solution(total_dofs))
        right_hand_side(1:vertex_count) = volume_load + &
            matmul(transpose(mass), neumann_jump) + &
            matmul(hypersingular, dirichlet_jump)
        right_hand_side(vertex_count + 1:) = matmul( &
            matrix(vertex_count + 1:, 1:vertex_count), dirichlet_jump)
        solution = 0.0_dp
        options = solver_options(method="lapack_lu")
        call solve_dense(matrix, right_hand_side, solution, options, statistics)
        if (.not. statistics%converged) return

        interior_solution = solution(1:vertex_count)
        exterior_flux = solution(vertex_count + 1:)
        status = 0
    end subroutine solve_laplace_symmetric_coupling_p1_p0

    subroutine solve_helmholtz_symmetric_coupling_p1_p0( &
            vertices, triangles, panel_start, panel_end, panel_nodes, &
            wavenumber, quadrature_order, volume_load, dirichlet_jump, &
            neumann_jump, interior_solution, exterior_flux, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: panel_nodes(:, :)
        real(dp), intent(in) :: wavenumber
        integer, intent(in) :: quadrature_order
        complex(dp), intent(in) :: volume_load(:), dirichlet_jump(:)
        complex(dp), intent(in) :: neumann_jump(:)
        complex(dp), intent(out) :: interior_solution(:), exterior_flux(:)
        integer, intent(out) :: status

        type(csc_z_t) :: sparse_matrix
        type(fortsparse_status_t) :: sparse_status
        complex(dp), allocatable :: hypersingular(:, :), matrix(:, :)
        complex(dp), allocatable :: right_hand_side(:), solution(:)
        complex(dp), allocatable :: triplet_values(:)
        real(dp), allocatable :: mass(:, :)
        integer, allocatable :: columns(:), rows(:)
        integer :: column, entry, operator_status, panel_count, row
        integer :: total_dofs, vertex_count

        interior_solution = cmplx(0.0_dp, 0.0_dp, dp)
        exterior_flux = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        vertex_count = size(vertices, 2)
        panel_count = size(panel_start, 2)
        total_dofs = vertex_count + panel_count
        if (size(volume_load) /= vertex_count .or. &
            size(dirichlet_jump) /= vertex_count) return
        if (size(neumann_jump) /= panel_count) return
        if (size(interior_solution) /= vertex_count .or. &
            size(exterior_flux) /= panel_count) return

        allocate(matrix(total_dofs, total_dofs))
        call assemble_helmholtz_symmetric_coupling_p1_p0( &
            vertices, triangles, panel_start, panel_end, panel_nodes, &
            wavenumber, quadrature_order, matrix, operator_status)
        if (operator_status /= 0) return

        allocate(mass(panel_count, vertex_count))
        call assemble_boundary_mass( &
            panel_start, panel_end, panel_nodes, mass, operator_status)
        if (operator_status /= 0) return
        allocate(hypersingular(vertex_count, vertex_count))
        call assemble_helmholtz_hypersingular_linear( &
            panel_start, panel_end, panel_nodes, vertex_count, wavenumber, &
            quadrature_order, hypersingular, operator_status)
        if (operator_status /= 0) return

        allocate(right_hand_side(total_dofs), solution(total_dofs))
        right_hand_side(1:vertex_count) = volume_load + &
            matmul(transpose(mass), neumann_jump) + &
            matmul(hypersingular, dirichlet_jump)
        right_hand_side(vertex_count + 1:) = matmul( &
            matrix(vertex_count + 1:, 1:vertex_count), dirichlet_jump)

        allocate(rows(total_dofs**2), columns(total_dofs**2))
        allocate(triplet_values(total_dofs**2))
        entry = 0
        do column = 1, total_dofs
            do row = 1, total_dofs
                entry = entry + 1
                rows(entry) = row
                columns(entry) = column
                triplet_values(entry) = matrix(row, column)
            end do
        end do
        call csc_from_triplet( &
            total_dofs, total_dofs, rows, columns, triplet_values, &
            sparse_matrix, sparse_status)
        if (sparse_status%code /= 0) return
        call sparse_solve_once( &
            sparse_matrix, right_hand_side, solution, sparse_status)
        if (sparse_status%code /= 0) return

        interior_solution = solution(1:vertex_count)
        exterior_flux = solution(vertex_count + 1:)
        status = 0
    end subroutine solve_helmholtz_symmetric_coupling_p1_p0

    pure logical function valid_discretization_shapes( &
            vertices, triangles, panel_start, panel_end, panel_nodes, &
            quadrature_order) result(valid)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: panel_nodes(:, :)
        integer, intent(in) :: quadrature_order

        integer :: panel_count, vertex_count

        valid = .false.
        vertex_count = size(vertices, 2)
        panel_count = size(panel_start, 2)
        if (size(vertices, 1) /= 2 .or. vertex_count < 3) return
        if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
        if (size(panel_start, 1) /= 2 .or. panel_count < 1) return
        if (size(panel_end, 1) /= 2 .or. &
            size(panel_end, 2) /= panel_count) return
        if (size(panel_nodes, 1) /= 2 .or. &
            size(panel_nodes, 2) /= panel_count) return
        if (quadrature_order < 1) return
        valid = .true.
    end function valid_discretization_shapes

    pure logical function panels_match_trace_nodes( &
            vertices, panel_start, panel_end, panel_nodes) result(matches)
        real(dp), intent(in) :: vertices(:, :)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: panel_nodes(:, :)

        real(dp) :: scale, tolerance
        integer :: panel

        matches = .false.
        do panel = 1, size(panel_start, 2)
            scale = max(1.0_dp, norm2(panel_start(:, panel)), &
                norm2(panel_end(:, panel)))
            tolerance = 64.0_dp * epsilon(1.0_dp) * scale
            if (norm2(panel_start(:, panel) - &
                vertices(:, panel_nodes(1, panel))) > tolerance) return
            if (norm2(panel_end(:, panel) - &
                vertices(:, panel_nodes(2, panel))) > tolerance) return
        end do
        matches = .true.
    end function panels_match_trace_nodes

    pure subroutine assemble_p1_stiffness( &
            vertices, triangles, stiffness, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(out) :: stiffness(:, :)
        integer, intent(out) :: status

        real(dp) :: area, determinant, gradients(2, 3), local_matrix(3, 3)
        real(dp) :: x(3), y(3)
        integer :: first_basis, first_node, second_basis, second_node
        integer :: triangle

        stiffness = 0.0_dp
        status = 2
        do triangle = 1, size(triangles, 2)
            x = vertices(1, triangles(:, triangle))
            y = vertices(2, triangles(:, triangle))
            determinant = (x(2) - x(1)) * (y(3) - y(1)) - &
                (x(3) - x(1)) * (y(2) - y(1))
            if (abs(determinant) <= 64.0_dp * epsilon(1.0_dp)) return
            area = 0.5_dp * abs(determinant)
            gradients(1, :) = [y(2) - y(3), y(3) - y(1), y(1) - y(2)] / &
                determinant
            gradients(2, :) = [x(3) - x(2), x(1) - x(3), x(2) - x(1)] / &
                determinant
            local_matrix = area * matmul(transpose(gradients), gradients)
            do second_basis = 1, 3
                second_node = triangles(second_basis, triangle)
                do first_basis = 1, 3
                    first_node = triangles(first_basis, triangle)
                    stiffness(first_node, second_node) = &
                        stiffness(first_node, second_node) + &
                        local_matrix(first_basis, second_basis)
                end do
            end do
        end do
        status = 0
    end subroutine assemble_p1_stiffness

    pure subroutine assemble_p1_mass(vertices, triangles, mass, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(out) :: mass(:, :)
        integer, intent(out) :: status

        real(dp) :: area, determinant, local_matrix(3, 3), x(3), y(3)
        integer :: first_basis, first_node, second_basis, second_node, triangle

        mass = 0.0_dp
        status = 2
        do triangle = 1, size(triangles, 2)
            x = vertices(1, triangles(:, triangle))
            y = vertices(2, triangles(:, triangle))
            determinant = (x(2) - x(1))*(y(3) - y(1)) - &
                (x(3) - x(1))*(y(2) - y(1))
            if (abs(determinant) <= 64.0_dp*epsilon(1.0_dp)) return
            area = 0.5_dp*abs(determinant)
            local_matrix = area/12.0_dp
            do first_basis = 1, 3
                local_matrix(first_basis, first_basis) = area/6.0_dp
            end do
            do second_basis = 1, 3
                second_node = triangles(second_basis, triangle)
                do first_basis = 1, 3
                    first_node = triangles(first_basis, triangle)
                    mass(first_node, second_node) = &
                        mass(first_node, second_node) + &
                        local_matrix(first_basis, second_basis)
                end do
            end do
        end do
        status = 0
    end subroutine assemble_p1_mass

    pure subroutine assemble_boundary_mass( &
            panel_start, panel_end, panel_nodes, mass, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: panel_nodes(:, :)
        real(dp), intent(out) :: mass(:, :)
        integer, intent(out) :: status

        real(dp) :: length
        integer :: endpoint, panel

        mass = 0.0_dp
        status = 2
        do panel = 1, size(panel_start, 2)
            length = norm2(panel_end(:, panel) - panel_start(:, panel))
            if (length <= 0.0_dp) return
            do endpoint = 1, 2
                mass(panel, panel_nodes(endpoint, panel)) = &
                    mass(panel, panel_nodes(endpoint, panel)) + 0.5_dp * length
            end do
        end do
        status = 0
    end subroutine assemble_boundary_mass

end module fortfem_laplace_symmetric_coupling_2d
