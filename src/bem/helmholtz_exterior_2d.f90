module fortfem_helmholtz_exterior_2d
    use fortfem_kinds, only: dp
    use fortfem_helmholtz_boundary_operators_2d, only: &
        assemble_helmholtz_double_layer_constant, &
        assemble_helmholtz_single_layer_constant
    use fortnum_linalg, only: dense_solve
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_special_complex_bessel, only: hankel_h1_real
    use fortnum_status, only: fortnum_status_t
    implicit none

    private

    public :: evaluate_helmholtz_combined_potential_adaptive_constant
    public :: evaluate_helmholtz_combined_potential_constant
    public :: solve_helmholtz_cfie_constant

contains

    subroutine evaluate_helmholtz_combined_potential_adaptive_constant( &
            panel_start, panel_end, wavenumber, coupling, density, points, &
            quadrature_order, relative_tolerance, max_depth, field, &
            error_estimate, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        real(dp), intent(in) :: wavenumber, coupling
        complex(dp), intent(in) :: density(:)
        real(dp), intent(in) :: points(:, :)
        integer, intent(in) :: quadrature_order
        real(dp), intent(in) :: relative_tolerance
        integer, intent(in) :: max_depth
        complex(dp), intent(out) :: field(:)
        real(dp), intent(out) :: error_estimate(:)
        integer, intent(out) :: status

        real(dp), allocatable :: nodes(:), weights(:)
        complex(dp) :: contribution
        real(dp) :: contribution_error, length, normal(2), tolerance
        real(dp) :: panel_end_local(2), panel_start_local(2), target(2)
        integer :: local_status, panel, panel_count, point

        field = (0.0_dp, 0.0_dp)
        error_estimate = 0.0_dp
        status = 1
        panel_count = size(panel_start, 2)
        if (size(panel_start, 1) /= 2 .or. size(panel_end, 1) /= 2) return
        if (size(panel_end, 2) /= panel_count .or. panel_count < 1) return
        if (size(density) /= panel_count .or. size(points, 1) /= 2) return
        if (size(field) /= size(points, 2) .or. &
            size(error_estimate) /= size(points, 2)) return
        if (wavenumber <= 0.0_dp .or. coupling <= 0.0_dp) return
        if (quadrature_order < 1 .or. relative_tolerance <= 0.0_dp) return
        if (max_depth < 1) return

        allocate(nodes(quadrature_order), weights(quadrature_order))
        call gauss_legendre_ab( &
            quadrature_order, 0.0_dp, 1.0_dp, nodes, weights)
        tolerance = relative_tolerance / real(panel_count, dp)
        do panel = 1, panel_count
            panel_start_local = panel_start(:, panel)
            panel_end_local = panel_end(:, panel)
            length = norm2(panel_end_local - panel_start_local)
            if (length <= 0.0_dp) then
                status = 2
                return
            end if
            normal(1) = (panel_end_local(2) - panel_start_local(2)) / length
            normal(2) = (panel_start_local(1) - panel_end_local(1)) / length
            do point = 1, size(points, 2)
                target = points(:, point)
                if (point_segment_distance( &
                    target, panel_start_local, panel_end_local) <= &
                    64.0_dp * epsilon(1.0_dp) * max(1.0_dp, length)) then
                    status = 2
                    return
                end if
                call adaptive_panel_potential( &
                    panel_start_local, panel_end_local, normal, length, &
                    target, wavenumber, coupling, nodes, weights, &
                    0.0_dp, 1.0_dp, tolerance, 0, max_depth, contribution, &
                    contribution_error, local_status)
                if (local_status /= 0) then
                    status = local_status
                    return
                end if
                field(point) = field(point) + density(panel) * contribution
                error_estimate(point) = error_estimate(point) + &
                    abs(density(panel)) * contribution_error
            end do
        end do
        status = 0
    end subroutine evaluate_helmholtz_combined_potential_adaptive_constant

    subroutine solve_helmholtz_cfie_constant( &
            panel_start, panel_end, wavenumber, coupling, dirichlet_trace, &
            quadrature_order, density, relative_residual, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        real(dp), intent(in) :: wavenumber, coupling
        complex(dp), intent(in) :: dirichlet_trace(:)
        integer, intent(in) :: quadrature_order
        complex(dp), intent(out) :: density(:)
        real(dp), intent(out) :: relative_residual
        integer, intent(out) :: status

        complex(dp), allocatable :: matrix(:, :), original_matrix(:, :)
        complex(dp), allocatable :: right_hand_side(:)
        real(dp), allocatable :: lengths(:)
        integer :: info, panel, panel_count

        density = (0.0_dp, 0.0_dp)
        relative_residual = huge(1.0_dp)
        status = 1
        panel_count = size(panel_start, 2)
        if (size(panel_start, 1) /= 2 .or. size(panel_end, 1) /= 2) return
        if (size(panel_end, 2) /= panel_count .or. panel_count < 1) return
        if (size(dirichlet_trace) /= panel_count .or. &
            size(density) /= panel_count) return
        if (wavenumber <= 0.0_dp .or. coupling <= 0.0_dp) return
        if (quadrature_order < 1) return

        allocate(matrix(panel_count, panel_count))
        allocate(original_matrix(panel_count, panel_count))
        allocate(right_hand_side(panel_count), lengths(panel_count))
        do panel = 1, panel_count
            lengths(panel) = norm2( &
                panel_end(:, panel) - panel_start(:, panel))
        end do
        if (any(lengths <= 0.0_dp)) return

        call assemble_helmholtz_double_layer_constant( &
            panel_start, panel_end, wavenumber, quadrature_order, &
            matrix, status)
        if (status /= 0) return
        call add_single_layer_cfie_term( &
            panel_start, panel_end, wavenumber, coupling, quadrature_order, &
            matrix, status)
        if (status /= 0) return
        do panel = 1, panel_count
            matrix(panel, panel) = matrix(panel, panel) + 0.5_dp * lengths(panel)
        end do

        right_hand_side = -lengths * dirichlet_trace
        original_matrix = matrix
        call dense_solve(matrix, right_hand_side, density, info)
        if (info /= 0) then
            status = 2
            return
        end if
        relative_residual = dense_relative_residual( &
            original_matrix, density, right_hand_side)
        status = 0
    end subroutine solve_helmholtz_cfie_constant

    subroutine add_single_layer_cfie_term( &
            panel_start, panel_end, wavenumber, coupling, quadrature_order, &
            matrix, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        real(dp), intent(in) :: wavenumber, coupling
        integer, intent(in) :: quadrature_order
        complex(dp), intent(inout) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: single_layer(:, :)
        integer :: panel_count

        panel_count = size(panel_start, 2)
        allocate(single_layer(panel_count, panel_count))
        call assemble_helmholtz_single_layer_constant( &
            panel_start, panel_end, wavenumber, quadrature_order, &
            single_layer, status)
        if (status /= 0) return
        matrix = matrix - cmplx(0.0_dp, coupling, dp) * single_layer
    end subroutine add_single_layer_cfie_term

    subroutine evaluate_helmholtz_combined_potential_constant( &
            panel_start, panel_end, wavenumber, coupling, density, points, &
            quadrature_order, field, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        real(dp), intent(in) :: wavenumber, coupling
        complex(dp), intent(in) :: density(:)
        real(dp), intent(in) :: points(:, :)
        integer, intent(in) :: quadrature_order
        complex(dp), intent(out) :: field(:)
        integer, intent(out) :: status

        real(dp), allocatable :: nodes(:), weights(:)
        real(dp) :: displacement(2), distance, length, normal(2)
        real(dp) :: source(2)
        complex(dp) :: double_kernel, hankel0, hankel1, single_kernel
        type(fortnum_status_t) :: special_status
        integer :: node, panel, panel_count, point

        field = (0.0_dp, 0.0_dp)
        status = 1
        panel_count = size(panel_start, 2)
        if (size(panel_start, 1) /= 2 .or. size(panel_end, 1) /= 2) return
        if (size(panel_end, 2) /= panel_count .or. panel_count < 1) return
        if (size(density) /= panel_count .or. size(points, 1) /= 2) return
        if (size(field) /= size(points, 2)) return
        if (wavenumber <= 0.0_dp .or. coupling <= 0.0_dp) return
        if (quadrature_order < 1) return

        allocate(nodes(quadrature_order), weights(quadrature_order))
        call gauss_legendre_ab( &
            quadrature_order, 0.0_dp, 1.0_dp, nodes, weights)
        status = 2
        do panel = 1, panel_count
            length = norm2(panel_end(:, panel) - panel_start(:, panel))
            if (length <= 0.0_dp) return
            normal(1) = (panel_end(2, panel) - panel_start(2, panel)) / length
            normal(2) = (panel_start(1, panel) - panel_end(1, panel)) / length
            do node = 1, quadrature_order
                source = panel_start(:, panel) + nodes(node) * &
                    (panel_end(:, panel) - panel_start(:, panel))
                do point = 1, size(points, 2)
                    displacement = points(:, point) - source
                    distance = norm2(displacement)
                    if (distance <= 0.0_dp) return
                    call hankel_h1_real( &
                        0, wavenumber * distance, hankel0, special_status)
                    if (special_status%code /= 0) return
                    call hankel_h1_real( &
                        1, wavenumber * distance, hankel1, special_status)
                    if (special_status%code /= 0) return
                    single_kernel = cmplx(0.0_dp, 0.25_dp, dp) * hankel0
                    double_kernel = cmplx(0.0_dp, 0.25_dp, dp) * &
                        wavenumber * hankel1 * &
                        dot_product(displacement, normal) / distance
                    field(point) = field(point) + density(panel) * length * &
                        weights(node) * (double_kernel - &
                        cmplx(0.0_dp, coupling, dp) * single_kernel)
                end do
            end do
        end do
        status = 0
    end subroutine evaluate_helmholtz_combined_potential_constant

    recursive subroutine adaptive_panel_potential( &
            panel_start, panel_end, normal, panel_length, target, wavenumber, &
            coupling, nodes, weights, interval_start, interval_end, tolerance, &
            depth, max_depth, value, error_estimate, status)
        real(dp), intent(in) :: panel_start(2), panel_end(2), normal(2)
        real(dp), intent(in) :: panel_length, target(2), wavenumber, coupling
        real(dp), intent(in) :: nodes(:), weights(:)
        real(dp), intent(in) :: interval_start, interval_end, tolerance
        integer, intent(in) :: depth, max_depth
        complex(dp), intent(out) :: value
        real(dp), intent(out) :: error_estimate
        integer, intent(out) :: status

        complex(dp) :: left_value, parent_value, right_value
        real(dp) :: left_error, midpoint, right_error
        integer :: local_status

        midpoint = 0.5_dp * (interval_start + interval_end)
        call panel_interval_potential( &
            panel_start, panel_end, normal, panel_length, target, wavenumber, &
            coupling, nodes, weights, interval_start, interval_end, &
            parent_value, status)
        if (status /= 0) return
        call panel_interval_potential( &
            panel_start, panel_end, normal, panel_length, target, wavenumber, &
            coupling, nodes, weights, interval_start, midpoint, &
            left_value, status)
        if (status /= 0) return
        call panel_interval_potential( &
            panel_start, panel_end, normal, panel_length, target, wavenumber, &
            coupling, nodes, weights, midpoint, interval_end, &
            right_value, status)
        if (status /= 0) return

        value = left_value + right_value
        error_estimate = abs(value - parent_value)
        if (error_estimate <= tolerance * max(1.0_dp, abs(value))) then
            status = 0
            return
        end if
        if (depth >= max_depth) then
            status = 3
            return
        end if

        call adaptive_panel_potential( &
            panel_start, panel_end, normal, panel_length, target, wavenumber, &
            coupling, nodes, weights, interval_start, midpoint, &
            0.5_dp * tolerance, depth + 1, max_depth, left_value, left_error, &
            local_status)
        if (local_status /= 0) then
            status = local_status
            return
        end if
        call adaptive_panel_potential( &
            panel_start, panel_end, normal, panel_length, target, wavenumber, &
            coupling, nodes, weights, midpoint, interval_end, &
            0.5_dp * tolerance, depth + 1, max_depth, right_value, right_error, &
            local_status)
        if (local_status /= 0) then
            status = local_status
            return
        end if
        value = left_value + right_value
        error_estimate = left_error + right_error
        status = 0
    end subroutine adaptive_panel_potential

    subroutine panel_interval_potential( &
            panel_start, panel_end, normal, panel_length, target, wavenumber, &
            coupling, nodes, weights, interval_start, interval_end, &
            value, status)
        real(dp), intent(in) :: panel_start(2), panel_end(2), normal(2)
        real(dp), intent(in) :: panel_length, target(2), wavenumber, coupling
        real(dp), intent(in) :: nodes(:), weights(:)
        real(dp), intent(in) :: interval_start, interval_end
        complex(dp), intent(out) :: value
        integer, intent(out) :: status

        complex(dp) :: double_kernel, hankel0, hankel1, single_kernel
        real(dp) :: displacement(2), distance, parameter, source(2)
        type(fortnum_status_t) :: special_status
        integer :: node

        value = (0.0_dp, 0.0_dp)
        status = 2
        do node = 1, size(nodes)
            parameter = interval_start + &
                (interval_end - interval_start) * nodes(node)
            source = panel_start + parameter * (panel_end - panel_start)
            displacement = target - source
            distance = norm2(displacement)
            if (distance <= 0.0_dp) return
            call hankel_h1_real( &
                0, wavenumber * distance, hankel0, special_status)
            if (special_status%code /= 0) return
            call hankel_h1_real( &
                1, wavenumber * distance, hankel1, special_status)
            if (special_status%code /= 0) return
            single_kernel = cmplx(0.0_dp, 0.25_dp, dp) * hankel0
            double_kernel = cmplx(0.0_dp, 0.25_dp, dp) * &
                wavenumber * hankel1 * dot_product(displacement, normal) / &
                distance
            value = value + weights(node) * (double_kernel - &
                cmplx(0.0_dp, coupling, dp) * single_kernel)
        end do
        value = panel_length * (interval_end - interval_start) * value
        status = 0
    end subroutine panel_interval_potential

    pure function point_segment_distance(point, start, finish) result(distance)
        real(dp), intent(in) :: point(2), start(2), finish(2)
        real(dp) :: distance

        real(dp) :: direction(2), parameter

        direction = finish - start
        parameter = dot_product(point - start, direction) / sum(direction**2)
        parameter = max(0.0_dp, min(1.0_dp, parameter))
        distance = norm2(point - (start + parameter * direction))
    end function point_segment_distance

    pure function dense_relative_residual(matrix, solution, rhs) result(residual)
        complex(dp), intent(in) :: matrix(:, :), solution(:), rhs(:)
        real(dp) :: residual

        complex(dp) :: row_value
        real(dp) :: residual_squared, rhs_squared
        integer :: column, row

        residual_squared = 0.0_dp
        rhs_squared = sum(abs(rhs)**2)
        do row = 1, size(solution)
            row_value = (0.0_dp, 0.0_dp)
            do column = 1, size(solution)
                row_value = row_value + matrix(row, column) * solution(column)
            end do
            residual_squared = residual_squared + abs(row_value - rhs(row))**2
        end do
        residual = sqrt(residual_squared / max(rhs_squared, tiny(1.0_dp)))
    end function dense_relative_residual

end module fortfem_helmholtz_exterior_2d
