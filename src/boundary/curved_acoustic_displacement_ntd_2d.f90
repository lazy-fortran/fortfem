module fortfem_curved_acoustic_displacement_ntd_2d
    !! Acoustic displacement-to-pressure map on a closed polygonal interface.
    !!
    !! For the time convention exp(-i omega t), the fluid momentum equation
    !! gives q = d_n p = rho omega^2 u_n.  The outgoing exterior Calderon
    !! equation
    !!
    !!     (1/2 M - K) p + V q = 0
    !!
    !! is discretized with continuous P1 pressure and panelwise P0 normal
    !! displacement.  Panel orientation must be counter-clockwise, so the
    !! right normal points from the bounded FEM domain into the exterior.
    use fortfem_helmholtz_boundary_operators_2d, only: &
        assemble_helmholtz_double_layer_mixed_linear, &
        assemble_helmholtz_hypersingular_linear, &
        assemble_helmholtz_single_layer_constant
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve
    use fortsparse, only: csc_from_triplet, csc_z_t, fortsparse_status_t, &
        sparse_destroy, sparse_factor, sparse_solve, sparse_solver_t
    implicit none

    private

    public :: assemble_curved_acoustic_displacement_ntd_form_2d
    public :: apply_curved_acoustic_displacement_ntd_2d

contains

    subroutine apply_curved_acoustic_displacement_ntd_2d( &
            panel_start, panel_end, panel_nodes, wavenumber, density, &
            angular_frequency, quadrature_order, normal_displacement, &
            pressure, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: panel_nodes(:, :)
        real(dp), intent(in) :: wavenumber, density, angular_frequency
        integer, intent(in) :: quadrature_order
        complex(dp), intent(in) :: normal_displacement(:)
        complex(dp), intent(out) :: pressure(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: matrix(:, :), neumann_operator(:, :)
        complex(dp), allocatable :: neumann_trace(:), right_hand_side(:)
        integer :: node_count, operator_status, panel_count

        pressure = (0.0_dp, 0.0_dp)
        status = 1
        panel_count = size(panel_start, 2)
        node_count = size(pressure)
        if (panel_count < 1 .or. node_count /= panel_count) return
        if (size(panel_start, 1) /= 2 .or. size(panel_end, 1) /= 2) return
        if (size(panel_end, 2) /= panel_count) return
        if (size(panel_nodes, 1) /= 2 .or. &
            size(panel_nodes, 2) /= panel_count) return
        if (any(panel_nodes < 1) .or. any(panel_nodes > node_count)) return
        if (size(normal_displacement) /= panel_count) return
        if (wavenumber <= 0.0_dp .or. density <= 0.0_dp .or. &
            angular_frequency <= 0.0_dp .or. quadrature_order < 1) return

        allocate(matrix(node_count, node_count))
        allocate(neumann_operator(node_count, panel_count))
        call assemble_combined_calderon_system( &
            panel_start, panel_end, panel_nodes, wavenumber, &
            quadrature_order, matrix, neumann_operator, operator_status)
        if (operator_status /= 0) return

        allocate(neumann_trace(panel_count), right_hand_side(panel_count))
        neumann_trace = density*angular_frequency**2*normal_displacement
        right_hand_side = matmul(neumann_operator, neumann_trace)
        call dense_solve(matrix, right_hand_side, pressure, operator_status)
        if (operator_status /= 0) then
            status = 2
            pressure = (0.0_dp, 0.0_dp)
            return
        end if
        status = 0
    end subroutine apply_curved_acoustic_displacement_ntd_2d

    subroutine assemble_curved_acoustic_displacement_ntd_form_2d( &
            panel_start, panel_end, panel_nodes, wavenumber, density, &
            angular_frequency, quadrature_order, form, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: panel_nodes(:, :)
        real(dp), intent(in) :: wavenumber, density, angular_frequency
        integer, intent(in) :: quadrature_order
        complex(dp), intent(out) :: form(:, :)
        integer, intent(out) :: status

        type(csc_z_t) :: sparse_matrix
        type(fortsparse_status_t) :: sparse_status
        type(sparse_solver_t) :: factor
        complex(dp), allocatable :: matrix(:, :), neumann_operator(:, :)
        complex(dp), allocatable :: pressure_operator(:, :), values(:)
        complex(dp), allocatable :: right_hand_side(:), solution(:)
        real(dp), allocatable :: displacement_trace(:, :), traction_mass(:, :)
        integer, allocatable :: columns(:), rows(:)
        real(dp) :: length, normal(2)
        integer :: column, component, endpoint, entry, node_count
        integer :: operator_status, panel, panel_count

        form = (0.0_dp, 0.0_dp)
        status = 1
        panel_count = size(panel_start, 2)
        node_count = panel_count
        if (size(form, 1) /= 2*node_count .or. &
            size(form, 2) /= 2*node_count) return
        if (.not. valid_boundary_inputs( &
            panel_start, panel_end, panel_nodes, node_count, wavenumber, &
            density, angular_frequency, quadrature_order)) return

        allocate(matrix(node_count, node_count))
        allocate(neumann_operator(node_count, panel_count))
        call assemble_combined_calderon_system( &
            panel_start, panel_end, panel_nodes, wavenumber, &
            quadrature_order, matrix, neumann_operator, operator_status)
        if (operator_status /= 0) return

        allocate(rows(node_count**2), columns(node_count**2))
        allocate(values(node_count**2))
        entry = 0
        do column = 1, node_count
            do panel = 1, node_count
                entry = entry + 1
                rows(entry) = panel
                columns(entry) = column
                values(entry) = matrix(panel, column)
            end do
        end do
        call csc_from_triplet( &
            node_count, node_count, rows, columns, values, sparse_matrix, &
            sparse_status)
        if (sparse_status%code /= 0) return
        call sparse_factor(factor, sparse_matrix, sparse_status)
        if (sparse_status%code /= 0) return

        allocate(pressure_operator(node_count, panel_count))
        allocate(right_hand_side(node_count), solution(node_count))
        do column = 1, panel_count
            right_hand_side = density*angular_frequency**2* &
                neumann_operator(:, column)
            call sparse_solve( &
                factor, right_hand_side, solution, sparse_status)
            if (sparse_status%code /= 0) then
                call sparse_destroy(factor)
                status = 2
                return
            end if
            pressure_operator(:, column) = solution
        end do
        call sparse_destroy(factor)

        allocate(displacement_trace(panel_count, 2*node_count))
        allocate(traction_mass(2*node_count, node_count))
        displacement_trace = 0.0_dp
        traction_mass = 0.0_dp
        do panel = 1, panel_count
            length = norm2(panel_end(:, panel) - panel_start(:, panel))
            normal = [ &
                panel_end(2, panel) - panel_start(2, panel), &
                panel_start(1, panel) - panel_end(1, panel)]/length
            do endpoint = 1, 2
                do component = 1, 2
                    displacement_trace(panel, &
                        2*panel_nodes(endpoint, panel) - 2 + component) = &
                        0.5_dp*normal(component)
                end do
            end do
            call add_panel_traction_mass( &
                panel, length, normal, panel_nodes, traction_mass)
        end do
        form = matmul(traction_mass, &
            matmul(pressure_operator, displacement_trace))
        status = 0
    end subroutine assemble_curved_acoustic_displacement_ntd_form_2d

    subroutine assemble_combined_calderon_system( &
            panel_start, panel_end, panel_nodes, wavenumber, quadrature_order, &
            matrix, neumann_operator, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: panel_nodes(:, :)
        real(dp), intent(in) :: wavenumber
        integer, intent(in) :: quadrature_order
        complex(dp), intent(out) :: matrix(:, :), neumann_operator(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: double_layer(:, :), hypersingular(:, :)
        complex(dp), allocatable :: single_layer(:, :)
        real(dp), allocatable :: mass(:, :)
        integer :: node_count, operator_status, panel_count

        matrix = (0.0_dp, 0.0_dp)
        neumann_operator = (0.0_dp, 0.0_dp)
        status = 1
        panel_count = size(panel_start, 2)
        node_count = size(matrix, 1)
        allocate(mass(panel_count, node_count))
        call assemble_boundary_mass( &
            panel_start, panel_end, panel_nodes, mass, operator_status)
        if (operator_status /= 0) return
        allocate(double_layer(panel_count, node_count))
        call assemble_helmholtz_double_layer_mixed_linear( &
            panel_start, panel_end, panel_nodes, node_count, wavenumber, &
            quadrature_order, double_layer, operator_status)
        if (operator_status /= 0) return
        allocate(single_layer(panel_count, panel_count))
        call assemble_helmholtz_single_layer_constant( &
            panel_start, panel_end, wavenumber, quadrature_order, &
            single_layer, operator_status)
        if (operator_status /= 0) return
        allocate(hypersingular(node_count, node_count))
        call assemble_helmholtz_hypersingular_linear( &
            panel_start, panel_end, panel_nodes, node_count, wavenumber, &
            quadrature_order, hypersingular, operator_status)
        if (operator_status /= 0) return

        ! Burton--Miller combination of both exterior Calderon equations.
        matrix = hypersingular + cmplx(0.0_dp, wavenumber, dp)*( &
            0.5_dp*mass - double_layer)
        neumann_operator = -( &
            0.5_dp*transpose(mass) + transpose(double_layer)) - &
            cmplx(0.0_dp, wavenumber, dp)*single_layer
        status = 0
    end subroutine assemble_combined_calderon_system

    pure logical function valid_boundary_inputs( &
            panel_start, panel_end, panel_nodes, node_count, wavenumber, &
            density, angular_frequency, quadrature_order)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: panel_nodes(:, :), node_count
        real(dp), intent(in) :: wavenumber, density, angular_frequency
        integer, intent(in) :: quadrature_order

        integer :: panel_count

        valid_boundary_inputs = .false.
        panel_count = size(panel_start, 2)
        if (panel_count < 1 .or. node_count /= panel_count) return
        if (size(panel_start, 1) /= 2 .or. size(panel_end, 1) /= 2) return
        if (size(panel_end, 2) /= panel_count) return
        if (size(panel_nodes, 1) /= 2 .or. &
            size(panel_nodes, 2) /= panel_count) return
        if (any(panel_nodes < 1) .or. any(panel_nodes > node_count)) return
        if (wavenumber <= 0.0_dp .or. density <= 0.0_dp .or. &
            angular_frequency <= 0.0_dp .or. quadrature_order < 1) return
        valid_boundary_inputs = .true.
    end function valid_boundary_inputs

    pure subroutine add_panel_traction_mass( &
            panel, length, normal, panel_nodes, traction_mass)
        integer, intent(in) :: panel
        real(dp), intent(in) :: length, normal(2)
        integer, intent(in) :: panel_nodes(:, :)
        real(dp), intent(inout) :: traction_mass(:, :)

        integer :: component, first, second

        do first = 1, 2
            do second = 1, 2
                do component = 1, 2
                    traction_mass( &
                        2*panel_nodes(first, panel) - 2 + component, &
                        panel_nodes(second, panel)) = traction_mass( &
                        2*panel_nodes(first, panel) - 2 + component, &
                        panel_nodes(second, panel)) + normal(component)* &
                        length*merge(2.0_dp, 1.0_dp, first == second)/6.0_dp
                end do
            end do
        end do
    end subroutine add_panel_traction_mass

    pure subroutine assemble_boundary_mass( &
            panel_start, panel_end, panel_nodes, mass, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: panel_nodes(:, :)
        real(dp), intent(out) :: mass(:, :)
        integer, intent(out) :: status

        real(dp) :: length
        integer :: endpoint, panel

        mass = 0.0_dp
        status = 1
        do panel = 1, size(panel_start, 2)
            length = norm2(panel_end(:, panel) - panel_start(:, panel))
            if (length <= 0.0_dp) return
            do endpoint = 1, 2
                mass(panel, panel_nodes(endpoint, panel)) = &
                    mass(panel, panel_nodes(endpoint, panel)) + 0.5_dp*length
            end do
        end do
        status = 0
    end subroutine assemble_boundary_mass

end module fortfem_curved_acoustic_displacement_ntd_2d
