module fortfem_laplace_boundary_operators_2d
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    private

    public :: assemble_laplace_double_layer_constant
    public :: assemble_laplace_single_layer_constant

contains

    subroutine assemble_laplace_double_layer_constant( &
            panel_start, panel_end, quadrature_order, matrix, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: quadrature_order
        real(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: lengths(:), nodes(:), weights(:)
        real(dp) :: first_end(2), first_start(2), integral
        real(dp) :: second_end(2), second_start(2), source_normal(2)
        integer :: first_panel, panel_count, second_panel

        matrix = 0.0_dp
        status = 1
        panel_count = size(panel_start, 2)
        if (size(panel_start, 1) /= 2 .or. size(panel_end, 1) /= 2) return
        if (size(panel_end, 2) /= panel_count .or. panel_count < 1) return
        if (size(matrix, 1) /= panel_count .or. &
            size(matrix, 2) /= panel_count) return
        if (quadrature_order < 1) return

        allocate(lengths(panel_count))
        do first_panel = 1, panel_count
            lengths(first_panel) = norm2( &
                panel_end(:, first_panel) - panel_start(:, first_panel))
        end do
        if (any(lengths <= 0.0_dp)) return

        allocate(nodes(quadrature_order), weights(quadrature_order))
        call gauss_legendre_ab( &
            quadrature_order, 0.0_dp, 1.0_dp, nodes, weights)

        do first_panel = 1, panel_count
            first_start = panel_start(:, first_panel)
            first_end = panel_end(:, first_panel)
            do second_panel = 1, panel_count
                if (first_panel == second_panel) cycle
                second_start = panel_start(:, second_panel)
                second_end = panel_end(:, second_panel)
                source_normal(1) = (second_end(2) - second_start(2)) / &
                    lengths(second_panel)
                source_normal(2) = (second_start(1) - second_end(1)) / &
                    lengths(second_panel)
                call double_layer_pair_integral( &
                    first_start, first_end, second_start, second_end, &
                    source_normal, lengths(first_panel), lengths(second_panel), &
                    nodes, weights, integral, status)
                if (status /= 0) return
                matrix(first_panel, second_panel) = integral
            end do
        end do
        status = 0
    end subroutine assemble_laplace_double_layer_constant

    subroutine assemble_laplace_single_layer_constant( &
            panel_start, panel_end, quadrature_order, matrix, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: quadrature_order
        real(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: nodes(:), weights(:)
        real(dp), allocatable :: lengths(:)
        real(dp) :: first_end(2), first_start(2)
        real(dp) :: logarithmic_integral
        real(dp) :: second_end(2), second_start(2)
        integer :: first_panel, panel_count, second_panel

        matrix = 0.0_dp
        status = 1
        panel_count = size(panel_start, 2)
        if (size(panel_start, 1) /= 2 .or. size(panel_end, 1) /= 2) return
        if (size(panel_end, 2) /= panel_count .or. panel_count < 1) return
        if (size(matrix, 1) /= panel_count .or. &
            size(matrix, 2) /= panel_count) return
        if (quadrature_order < 1) return

        allocate(lengths(panel_count))
        do first_panel = 1, panel_count
            lengths(first_panel) = norm2( &
                panel_end(:, first_panel) - panel_start(:, first_panel))
        end do
        if (any(lengths <= 0.0_dp)) return

        allocate(nodes(quadrature_order), weights(quadrature_order))
        call gauss_legendre_ab( &
            quadrature_order, 0.0_dp, 1.0_dp, nodes, weights)

        do first_panel = 1, panel_count
            logarithmic_integral = lengths(first_panel)**2 * &
                (log(lengths(first_panel)) - 1.5_dp)
            matrix(first_panel, first_panel) = &
                -logarithmic_integral / (2.0_dp * acos(-1.0_dp))

            do second_panel = 1, first_panel - 1
                first_start = panel_start(:, first_panel)
                first_end = panel_end(:, first_panel)
                second_start = panel_start(:, second_panel)
                second_end = panel_end(:, second_panel)
                call panel_pair_logarithmic_integral( &
                    first_start, first_end, second_start, second_end, &
                    lengths(first_panel), lengths(second_panel), nodes, weights, &
                    logarithmic_integral, status)
                if (status /= 0) return
                matrix(first_panel, second_panel) = &
                    -logarithmic_integral / (2.0_dp * acos(-1.0_dp))
                matrix(second_panel, first_panel) = matrix(first_panel, second_panel)
            end do
        end do
        status = 0
    end subroutine assemble_laplace_single_layer_constant

    subroutine double_layer_pair_integral( &
            first_start, first_end, second_start, second_end, source_normal, &
            first_length, second_length, nodes, weights, integral, status)
        real(dp), intent(in) :: first_start(2), first_end(2)
        real(dp), intent(in) :: second_start(2), second_end(2)
        real(dp), intent(in) :: source_normal(2)
        real(dp), intent(in) :: first_length, second_length
        real(dp), intent(in) :: nodes(:), weights(:)
        real(dp), intent(out) :: integral
        integer, intent(out) :: status

        real(dp) :: first_away(2), second_away(2)
        integer :: shared_endpoint_count

        call shared_endpoint_vectors( &
            first_start, first_end, second_start, second_end, &
            first_length, second_length, first_away, second_away, &
            shared_endpoint_count)
        select case (shared_endpoint_count)
        case (0)
            call regular_double_layer_pair_integral( &
                first_start, first_end, second_start, second_end, &
                source_normal, first_length, second_length, nodes, weights, &
                integral, status)
        case (1)
            call endpoint_singular_double_layer_pair_integral( &
                first_away, second_away, source_normal, first_length, &
                second_length, nodes, weights, integral, status)
        case default
            integral = 0.0_dp
            status = 2
        end select
    end subroutine double_layer_pair_integral

    subroutine panel_pair_logarithmic_integral( &
            first_start, first_end, second_start, second_end, &
            first_length, second_length, nodes, weights, integral, status)
        real(dp), intent(in) :: first_start(2), first_end(2)
        real(dp), intent(in) :: second_start(2), second_end(2)
        real(dp), intent(in) :: first_length, second_length
        real(dp), intent(in) :: nodes(:), weights(:)
        real(dp), intent(out) :: integral
        integer, intent(out) :: status

        real(dp) :: first_away(2), second_away(2)
        integer :: shared_endpoint_count

        call shared_endpoint_vectors( &
            first_start, first_end, second_start, second_end, &
            first_length, second_length, first_away, second_away, &
            shared_endpoint_count)
        select case (shared_endpoint_count)
        case (0)
            call regular_pair_integral( &
                first_start, first_end, second_start, second_end, &
                first_length, second_length, nodes, weights, integral)
            status = 0
        case (1)
            call endpoint_singular_pair_integral( &
                first_away, second_away, first_length, second_length, &
                nodes, weights, integral, status)
        case default
            integral = 0.0_dp
            status = 2
        end select
    end subroutine panel_pair_logarithmic_integral

    subroutine shared_endpoint_vectors( &
            first_start, first_end, second_start, second_end, &
            first_length, second_length, first_away, second_away, count)
        real(dp), intent(in) :: first_start(2), first_end(2)
        real(dp), intent(in) :: second_start(2), second_end(2)
        real(dp), intent(in) :: first_length, second_length
        real(dp), intent(out) :: first_away(2), second_away(2)
        integer, intent(out) :: count

        real(dp) :: tolerance

        count = 0
        first_away = 0.0_dp
        second_away = 0.0_dp
        tolerance = 64.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp, first_length, second_length)

        if (norm2(first_start - second_start) <= tolerance) then
            count = count + 1
            first_away = first_end - first_start
            second_away = second_end - second_start
        end if
        if (norm2(first_start - second_end) <= tolerance) then
            count = count + 1
            first_away = first_end - first_start
            second_away = second_start - second_end
        end if
        if (norm2(first_end - second_start) <= tolerance) then
            count = count + 1
            first_away = first_start - first_end
            second_away = second_end - second_start
        end if
        if (norm2(first_end - second_end) <= tolerance) then
            count = count + 1
            first_away = first_start - first_end
            second_away = second_start - second_end
        end if
    end subroutine shared_endpoint_vectors

    subroutine endpoint_singular_pair_integral( &
            first_away, second_away, first_length, second_length, &
            nodes, weights, integral, status)
        real(dp), intent(in) :: first_away(2), second_away(2)
        real(dp), intent(in) :: first_length, second_length
        real(dp), intent(in) :: nodes(:), weights(:)
        real(dp), intent(out) :: integral
        integer, intent(out) :: status

        real(dp) :: first_distance, second_distance
        integer :: node

        integral = 0.0_dp
        status = 2
        do node = 1, size(nodes)
            first_distance = norm2(first_away - nodes(node) * second_away)
            second_distance = norm2(nodes(node) * first_away - second_away)
            if (first_distance <= 0.0_dp .or. second_distance <= 0.0_dp) return
            integral = integral + weights(node) * ( &
                -0.5_dp + 0.5_dp * &
                (log(first_distance) + log(second_distance)))
        end do
        integral = first_length * second_length * integral
        status = 0
    end subroutine endpoint_singular_pair_integral

    subroutine endpoint_singular_double_layer_pair_integral( &
            first_away, second_away, source_normal, first_length, &
            second_length, nodes, weights, integral, status)
        real(dp), intent(in) :: first_away(2), second_away(2)
        real(dp), intent(in) :: source_normal(2)
        real(dp), intent(in) :: first_length, second_length
        real(dp), intent(in) :: nodes(:), weights(:)
        real(dp), intent(out) :: integral
        integer, intent(out) :: status

        real(dp) :: first_distance_squared, normal_projection
        real(dp) :: second_distance_squared
        integer :: node

        integral = 0.0_dp
        status = 2
        normal_projection = dot_product(first_away, source_normal)
        do node = 1, size(nodes)
            first_distance_squared = sum( &
                (first_away - nodes(node) * second_away)**2)
            second_distance_squared = sum( &
                (nodes(node) * first_away - second_away)**2)
            if (first_distance_squared <= 0.0_dp .or. &
                second_distance_squared <= 0.0_dp) return
            integral = integral + weights(node) * normal_projection * ( &
                1.0_dp / first_distance_squared + &
                nodes(node) / second_distance_squared)
        end do
        integral = first_length * second_length * integral / &
            (2.0_dp * acos(-1.0_dp))
        status = 0
    end subroutine endpoint_singular_double_layer_pair_integral

    subroutine regular_pair_integral( &
            first_start, first_end, second_start, second_end, &
            first_length, second_length, nodes, weights, integral)
        real(dp), intent(in) :: first_start(2), first_end(2)
        real(dp), intent(in) :: second_start(2), second_end(2)
        real(dp), intent(in) :: first_length, second_length
        real(dp), intent(in) :: nodes(:), weights(:)
        real(dp), intent(out) :: integral

        real(dp) :: first_point(2), second_point(2)
        integer :: first_node, second_node

        integral = 0.0_dp
        do first_node = 1, size(nodes)
            first_point = first_start + &
                nodes(first_node) * (first_end - first_start)
            do second_node = 1, size(nodes)
                second_point = second_start + &
                    nodes(second_node) * (second_end - second_start)
                integral = integral + weights(first_node) * &
                    weights(second_node) * log(norm2(first_point - second_point))
            end do
        end do
        integral = first_length * second_length * integral
    end subroutine regular_pair_integral

    subroutine regular_double_layer_pair_integral( &
            first_start, first_end, second_start, second_end, source_normal, &
            first_length, second_length, nodes, weights, integral, status)
        real(dp), intent(in) :: first_start(2), first_end(2)
        real(dp), intent(in) :: second_start(2), second_end(2)
        real(dp), intent(in) :: source_normal(2)
        real(dp), intent(in) :: first_length, second_length
        real(dp), intent(in) :: nodes(:), weights(:)
        real(dp), intent(out) :: integral
        integer, intent(out) :: status

        real(dp) :: displacement(2), distance_squared
        real(dp) :: first_point(2), second_point(2)
        integer :: first_node, second_node

        integral = 0.0_dp
        status = 2
        do first_node = 1, size(nodes)
            first_point = first_start + &
                nodes(first_node) * (first_end - first_start)
            do second_node = 1, size(nodes)
                second_point = second_start + &
                    nodes(second_node) * (second_end - second_start)
                displacement = first_point - second_point
                distance_squared = sum(displacement**2)
                if (distance_squared <= 0.0_dp) return
                integral = integral + weights(first_node) * &
                    weights(second_node) * &
                    dot_product(displacement, source_normal) / distance_squared
            end do
        end do
        integral = first_length * second_length * integral / &
            (2.0_dp * acos(-1.0_dp))
        status = 0
    end subroutine regular_double_layer_pair_integral

end module fortfem_laplace_boundary_operators_2d
