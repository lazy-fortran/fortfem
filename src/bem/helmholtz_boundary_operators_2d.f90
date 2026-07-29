module fortfem_helmholtz_boundary_operators_2d
    use fortfem_kinds, only: dp
    use fortfem_laplace_boundary_operators_2d, only: &
        assemble_laplace_double_layer_constant, &
        assemble_laplace_double_layer_mixed_linear, &
        assemble_laplace_single_layer_constant
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_special_complex_bessel, only: hankel_h1_real
    use fortnum_status, only: fortnum_status_t
    implicit none

    private

    public :: assemble_helmholtz_adjoint_double_layer_constant
    public :: assemble_helmholtz_double_layer_constant
    public :: assemble_helmholtz_double_layer_mixed_linear
    public :: assemble_helmholtz_hypersingular_linear
    public :: assemble_helmholtz_single_layer_constant
    public :: assemble_helmholtz_single_layer_linear

contains

    subroutine assemble_helmholtz_double_layer_mixed_linear( &
            panel_start, panel_end, panel_nodes, node_count, wavenumber, &
            quadrature_order, matrix, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: panel_nodes(:, :)
        integer, intent(in) :: node_count
        real(dp), intent(in) :: wavenumber
        integer, intent(in) :: quadrature_order
        complex(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp) :: correction(2)
        real(dp), allocatable :: laplace(:, :), lengths(:), nodes(:), weights(:)
        real(dp) :: source_normal(2)
        integer :: endpoint, source_panel, target_panel

        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (node_count < 1 .or. size(panel_start, 1) /= 2 .or. &
            size(panel_end, 1) /= 2) return
        if (size(panel_end, 2) /= size(panel_start, 2) .or. &
            size(panel_start, 2) < 1) return
        if (size(panel_nodes, 1) /= 2 .or. &
            size(panel_nodes, 2) /= size(panel_start, 2)) return
        if (any(panel_nodes < 1) .or. any(panel_nodes > node_count)) return
        if (size(matrix, 1) /= size(panel_start, 2) .or. &
            size(matrix, 2) /= node_count) return
        if (wavenumber <= 0.0_dp .or. quadrature_order < 1) return

        allocate(laplace(size(matrix, 1), size(matrix, 2)))
        call assemble_laplace_double_layer_mixed_linear( &
            panel_start, panel_end, panel_nodes, node_count, quadrature_order, &
            laplace, status)
        if (status /= 0) return
        matrix = cmplx(laplace, 0.0_dp, dp)

        allocate(lengths(size(panel_start, 2)))
        allocate(nodes(quadrature_order), weights(quadrature_order))
        do source_panel = 1, size(panel_start, 2)
            lengths(source_panel) = norm2( &
                panel_end(:, source_panel) - panel_start(:, source_panel))
            if (lengths(source_panel) <= 0.0_dp) return
        end do
        call gauss_legendre_ab( &
            quadrature_order, 0.0_dp, 1.0_dp, nodes, weights)

        do target_panel = 1, size(panel_start, 2)
            do source_panel = 1, size(panel_start, 2)
                if (target_panel == source_panel) cycle
                source_normal = [ &
                    panel_end(2, source_panel) - panel_start(2, source_panel), &
                    panel_start(1, source_panel) - panel_end(1, source_panel)] / &
                    lengths(source_panel)
                call regular_double_layer_remainder_mixed_moment( &
                    panel_start(:, target_panel), panel_end(:, target_panel), &
                    panel_start(:, source_panel), panel_end(:, source_panel), &
                    source_normal, lengths(target_panel), lengths(source_panel), &
                    wavenumber, nodes, weights, correction, status)
                if (status /= 0) return
                do endpoint = 1, 2
                    matrix(target_panel, panel_nodes(endpoint, source_panel)) = &
                        matrix(target_panel, &
                        panel_nodes(endpoint, source_panel)) + &
                        correction(endpoint)
                end do
            end do
        end do
        status = 0
    end subroutine assemble_helmholtz_double_layer_mixed_linear

    subroutine assemble_helmholtz_hypersingular_linear( &
            panel_start, panel_end, panel_nodes, node_count, wavenumber, &
            quadrature_order, matrix, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: panel_nodes(:, :)
        integer, intent(in) :: node_count
        real(dp), intent(in) :: wavenumber
        integer, intent(in) :: quadrature_order
        complex(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp), parameter :: endpoint_sign(2) = [-1.0_dp, 1.0_dp]
        complex(dp), allocatable :: constant_kernel(:, :), linear_kernel(:, :)
        integer, allocatable :: discontinuous_nodes(:, :)
        real(dp), allocatable :: lengths(:), normals(:, :)
        real(dp) :: first_derivative, normal_product, second_derivative
        integer :: first_endpoint, first_local_node, first_node, first_panel
        integer :: panel_count, second_endpoint, second_local_node
        integer :: second_node, second_panel

        matrix = (0.0_dp, 0.0_dp)
        status = 1
        panel_count = size(panel_start, 2)
        if (node_count < 1 .or. size(panel_nodes, 1) /= 2) return
        if (size(panel_start, 1) /= 2 .or. size(panel_end, 1) /= 2) return
        if (size(panel_end, 2) /= panel_count .or. panel_count < 1) return
        if (size(panel_nodes, 2) /= panel_count) return
        if (size(matrix, 1) /= node_count .or. &
            size(matrix, 2) /= node_count) return
        if (any(panel_nodes < 1) .or. any(panel_nodes > node_count)) return
        if (wavenumber <= 0.0_dp .or. quadrature_order < 1) return

        allocate(constant_kernel(panel_count, panel_count))
        allocate(discontinuous_nodes(2, panel_count))
        allocate(lengths(panel_count), normals(2, panel_count))
        allocate(linear_kernel(2 * panel_count, 2 * panel_count))
        do first_panel = 1, panel_count
            lengths(first_panel) = norm2( &
                panel_end(:, first_panel) - panel_start(:, first_panel))
            if (lengths(first_panel) <= 0.0_dp) return
            normals(1, first_panel) = ( &
                panel_end(2, first_panel) - panel_start(2, first_panel)) / &
                lengths(first_panel)
            normals(2, first_panel) = ( &
                panel_start(1, first_panel) - panel_end(1, first_panel)) / &
                lengths(first_panel)
            discontinuous_nodes(1, first_panel) = 2 * first_panel - 1
            discontinuous_nodes(2, first_panel) = 2 * first_panel
        end do

        call assemble_helmholtz_single_layer_constant( &
            panel_start, panel_end, wavenumber, quadrature_order, &
            constant_kernel, status)
        if (status /= 0) return
        call assemble_helmholtz_single_layer_linear( &
            panel_start, panel_end, discontinuous_nodes, 2 * panel_count, &
            wavenumber, quadrature_order, linear_kernel, status)
        if (status /= 0) return

        do first_panel = 1, panel_count
            do second_panel = 1, panel_count
                normal_product = dot_product( &
                    normals(:, first_panel), normals(:, second_panel))
                do first_endpoint = 1, 2
                    first_node = panel_nodes(first_endpoint, first_panel)
                    first_local_node = 2 * (first_panel - 1) + first_endpoint
                    first_derivative = endpoint_sign(first_endpoint) / &
                        lengths(first_panel)
                    do second_endpoint = 1, 2
                        second_node = panel_nodes( &
                            second_endpoint, second_panel)
                        second_local_node = &
                            2 * (second_panel - 1) + second_endpoint
                        second_derivative = endpoint_sign(second_endpoint) / &
                            lengths(second_panel)
                        matrix(first_node, second_node) = &
                            matrix(first_node, second_node) + &
                            first_derivative * &
                            constant_kernel(first_panel, second_panel) * &
                            second_derivative - wavenumber**2 * &
                            normal_product * &
                            linear_kernel(first_local_node, second_local_node)
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_helmholtz_hypersingular_linear

    subroutine assemble_helmholtz_single_layer_linear( &
            panel_start, panel_end, panel_nodes, node_count, wavenumber, &
            quadrature_order, matrix, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: panel_nodes(:, :)
        integer, intent(in) :: node_count
        real(dp), intent(in) :: wavenumber
        integer, intent(in) :: quadrature_order
        complex(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp) :: local_matrix(2, 2)
        real(dp), allocatable :: lengths(:), nodes(:), weights(:)
        real(dp) :: first_away(2), second_away(2)
        real(dp) :: first_end_local(2), first_start_local(2)
        real(dp) :: second_end_local(2), second_start_local(2)
        integer :: first_endpoint, first_node, first_panel, first_shared
        integer :: panel_count, second_endpoint, second_node, second_panel
        integer :: second_shared, shared_count

        matrix = (0.0_dp, 0.0_dp)
        status = 1
        panel_count = size(panel_start, 2)
        if (node_count < 1 .or. size(panel_nodes, 1) /= 2) return
        if (size(panel_start, 1) /= 2 .or. size(panel_end, 1) /= 2) return
        if (size(panel_end, 2) /= panel_count .or. panel_count < 1) return
        if (size(panel_nodes, 2) /= panel_count) return
        if (size(matrix, 1) /= node_count .or. &
            size(matrix, 2) /= node_count) return
        if (any(panel_nodes < 1) .or. any(panel_nodes > node_count)) return
        if (wavenumber <= 0.0_dp .or. quadrature_order < 1) return

        allocate(lengths(panel_count))
        allocate(nodes(quadrature_order), weights(quadrature_order))
        do first_panel = 1, panel_count
            lengths(first_panel) = norm2( &
                panel_end(:, first_panel) - panel_start(:, first_panel))
        end do
        if (any(lengths <= 0.0_dp)) return
        call gauss_legendre_ab( &
            quadrature_order, 0.0_dp, 1.0_dp, nodes, weights)

        do first_panel = 1, panel_count
            first_start_local = panel_start(:, first_panel)
            first_end_local = panel_end(:, first_panel)
            do second_panel = 1, panel_count
                second_start_local = panel_start(:, second_panel)
                second_end_local = panel_end(:, second_panel)
                if (first_panel == second_panel) then
                    call self_linear_single_layer_moment( &
                        lengths(first_panel), wavenumber, nodes, weights, &
                        local_matrix, status)
                else
                    call shared_linear_panel_endpoint( &
                        first_start_local, first_end_local, &
                        second_start_local, second_end_local, &
                        lengths(first_panel), lengths(second_panel), &
                        first_away, second_away, first_shared, second_shared, &
                        shared_count)
                    select case (shared_count)
                    case (0)
                        call regular_linear_single_layer_moment( &
                            first_start_local, first_end_local, &
                            second_start_local, second_end_local, &
                            lengths(first_panel), lengths(second_panel), &
                            wavenumber, nodes, weights, local_matrix, status)
                    case (1)
                        call adjacent_linear_single_layer_moment( &
                            first_start_local, first_end_local, &
                            second_start_local, second_end_local, &
                            first_away, second_away, first_shared, second_shared, &
                            lengths(first_panel), lengths(second_panel), &
                            wavenumber, nodes, weights, local_matrix, status)
                    case default
                        status = 2
                    end select
                end if
                if (status /= 0) return

                do first_endpoint = 1, 2
                    first_node = panel_nodes(first_endpoint, first_panel)
                    do second_endpoint = 1, 2
                        second_node = panel_nodes(second_endpoint, second_panel)
                        matrix(first_node, second_node) = &
                            matrix(first_node, second_node) + &
                            local_matrix(first_endpoint, second_endpoint)
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_helmholtz_single_layer_linear

    subroutine assemble_helmholtz_adjoint_double_layer_constant( &
            panel_start, panel_end, wavenumber, quadrature_order, matrix, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        real(dp), intent(in) :: wavenumber
        integer, intent(in) :: quadrature_order
        complex(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: double_layer(:, :)
        integer :: panel_count

        matrix = (0.0_dp, 0.0_dp)
        status = 1
        panel_count = size(panel_start, 2)
        if (size(matrix, 1) /= panel_count .or. &
            size(matrix, 2) /= panel_count) return

        allocate(double_layer(panel_count, panel_count))
        call assemble_helmholtz_double_layer_constant( &
            panel_start, panel_end, wavenumber, quadrature_order, &
            double_layer, status)
        if (status /= 0) return
        matrix = transpose(double_layer)
    end subroutine assemble_helmholtz_adjoint_double_layer_constant

    subroutine assemble_helmholtz_double_layer_constant( &
            panel_start, panel_end, wavenumber, quadrature_order, matrix, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        real(dp), intent(in) :: wavenumber
        integer, intent(in) :: quadrature_order
        complex(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp) :: correction
        real(dp), allocatable :: laplace(:, :), lengths(:)
        real(dp), allocatable :: nodes(:), weights(:)
        real(dp) :: first_end(2), first_start(2)
        real(dp) :: second_end(2), second_start(2), source_normal(2)
        integer :: first_panel, panel_count, second_panel

        matrix = (0.0_dp, 0.0_dp)
        status = 1
        panel_count = size(panel_start, 2)
        if (size(panel_start, 1) /= 2 .or. size(panel_end, 1) /= 2) return
        if (size(panel_end, 2) /= panel_count .or. panel_count < 1) return
        if (size(matrix, 1) /= panel_count .or. &
            size(matrix, 2) /= panel_count) return
        if (wavenumber <= 0.0_dp .or. quadrature_order < 1) return

        allocate(laplace(panel_count, panel_count), lengths(panel_count))
        allocate(nodes(quadrature_order), weights(quadrature_order))
        call assemble_laplace_double_layer_constant( &
            panel_start, panel_end, quadrature_order, laplace, status)
        if (status /= 0) return
        matrix = cmplx(laplace, 0.0_dp, dp)

        do first_panel = 1, panel_count
            lengths(first_panel) = norm2( &
                panel_end(:, first_panel) - panel_start(:, first_panel))
        end do
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
                call regular_double_layer_remainder_integral( &
                    first_start, first_end, second_start, second_end, &
                    source_normal, lengths(first_panel), lengths(second_panel), &
                    wavenumber, nodes, weights, correction, status)
                if (status /= 0) return
                matrix(first_panel, second_panel) = &
                    matrix(first_panel, second_panel) + correction
            end do
        end do
        status = 0
    end subroutine assemble_helmholtz_double_layer_constant

    subroutine assemble_helmholtz_single_layer_constant( &
            panel_start, panel_end, wavenumber, quadrature_order, matrix, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        real(dp), intent(in) :: wavenumber
        integer, intent(in) :: quadrature_order
        complex(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp) :: correction
        real(dp), allocatable :: laplace(:, :), lengths(:)
        real(dp), allocatable :: nodes(:), weights(:)
        real(dp) :: first_end(2), first_start(2)
        real(dp) :: second_end(2), second_start(2)
        integer :: first_panel, panel_count, second_panel

        matrix = (0.0_dp, 0.0_dp)
        status = 1
        panel_count = size(panel_start, 2)
        if (size(panel_start, 1) /= 2 .or. size(panel_end, 1) /= 2) return
        if (size(panel_end, 2) /= panel_count .or. panel_count < 1) return
        if (size(matrix, 1) /= panel_count .or. &
            size(matrix, 2) /= panel_count) return
        if (wavenumber <= 0.0_dp .or. quadrature_order < 1) return

        allocate(laplace(panel_count, panel_count), lengths(panel_count))
        allocate(nodes(quadrature_order), weights(quadrature_order))
        call assemble_laplace_single_layer_constant( &
            panel_start, panel_end, quadrature_order, laplace, status)
        if (status /= 0) return
        matrix = cmplx(laplace, 0.0_dp, dp)

        do first_panel = 1, panel_count
            lengths(first_panel) = norm2( &
                panel_end(:, first_panel) - panel_start(:, first_panel))
        end do
        call gauss_legendre_ab( &
            quadrature_order, 0.0_dp, 1.0_dp, nodes, weights)

        do first_panel = 1, panel_count
            call self_remainder_integral( &
                lengths(first_panel), wavenumber, nodes, weights, &
                correction, status)
            if (status /= 0) return
            matrix(first_panel, first_panel) = &
                matrix(first_panel, first_panel) + correction

            first_start = panel_start(:, first_panel)
            first_end = panel_end(:, first_panel)
            do second_panel = 1, first_panel - 1
                second_start = panel_start(:, second_panel)
                second_end = panel_end(:, second_panel)
                call regular_remainder_integral( &
                    first_start, first_end, second_start, second_end, &
                    lengths(first_panel), lengths(second_panel), wavenumber, &
                    nodes, weights, correction, status)
                if (status /= 0) return
                matrix(first_panel, second_panel) = &
                    matrix(first_panel, second_panel) + correction
                matrix(second_panel, first_panel) = matrix(first_panel, second_panel)
            end do
        end do
        status = 0
    end subroutine assemble_helmholtz_single_layer_constant

    subroutine self_remainder_integral( &
            panel_length, wavenumber, nodes, weights, integral, status)
        real(dp), intent(in) :: panel_length, wavenumber
        real(dp), intent(in) :: nodes(:), weights(:)
        complex(dp), intent(out) :: integral
        integer, intent(out) :: status

        complex(dp) :: remainder
        integer :: node

        integral = (0.0_dp, 0.0_dp)
        do node = 1, size(nodes)
            call helmholtz_laplace_remainder( &
                panel_length * nodes(node), wavenumber, remainder, status)
            if (status /= 0) return
            integral = integral + &
                2.0_dp * weights(node) * (1.0_dp - nodes(node)) * remainder
        end do
        integral = panel_length**2 * integral
        status = 0
    end subroutine self_remainder_integral

    subroutine self_linear_single_layer_moment( &
            panel_length, wavenumber, nodes, weights, matrix, status)
        real(dp), intent(in) :: panel_length, wavenumber
        real(dp), intent(in) :: nodes(:), weights(:)
        complex(dp), intent(out) :: matrix(2, 2)
        integer, intent(out) :: status

        complex(dp) :: diagonal_correction, off_diagonal_correction, remainder
        real(dp) :: diagonal_laplace, off_diagonal_laplace
        real(dp) :: diagonal_weight, off_diagonal_weight, parameter
        integer :: node

        diagonal_laplace = panel_length**2 / (2.0_dp * acos(-1.0_dp)) * &
            (7.0_dp / 16.0_dp - 0.25_dp * log(panel_length))
        off_diagonal_laplace = panel_length**2 / &
            (2.0_dp * acos(-1.0_dp)) * &
            (5.0_dp / 16.0_dp - 0.25_dp * log(panel_length))
        diagonal_correction = (0.0_dp, 0.0_dp)
        off_diagonal_correction = (0.0_dp, 0.0_dp)

        do node = 1, size(nodes)
            parameter = nodes(node)
            diagonal_weight = parameter**3 / 3.0_dp - parameter + &
                2.0_dp / 3.0_dp
            off_diagonal_weight = 1.0_dp / 3.0_dp - parameter**3 / 3.0_dp
            call helmholtz_laplace_remainder( &
                panel_length * parameter, wavenumber, remainder, status)
            if (status /= 0) return
            diagonal_correction = diagonal_correction + &
                weights(node) * diagonal_weight * remainder
            off_diagonal_correction = off_diagonal_correction + &
                weights(node) * off_diagonal_weight * remainder
        end do
        diagonal_correction = panel_length**2 * diagonal_correction
        off_diagonal_correction = panel_length**2 * off_diagonal_correction

        matrix(1, 1) = diagonal_laplace + diagonal_correction
        matrix(2, 2) = matrix(1, 1)
        matrix(1, 2) = off_diagonal_laplace + off_diagonal_correction
        matrix(2, 1) = matrix(1, 2)
        status = 0
    end subroutine self_linear_single_layer_moment

    subroutine shared_linear_panel_endpoint( &
            first_start, first_end, second_start, second_end, &
            first_length, second_length, first_away, second_away, &
            first_shared, second_shared, count)
        real(dp), intent(in) :: first_start(2), first_end(2)
        real(dp), intent(in) :: second_start(2), second_end(2)
        real(dp), intent(in) :: first_length, second_length
        real(dp), intent(out) :: first_away(2), second_away(2)
        integer, intent(out) :: first_shared, second_shared, count

        real(dp) :: tolerance

        count = 0
        first_away = 0.0_dp
        second_away = 0.0_dp
        first_shared = 0
        second_shared = 0
        tolerance = 64.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp, first_length, second_length)

        if (norm2(first_start - second_start) <= tolerance) then
            count = count + 1
            first_away = first_end - first_start
            second_away = second_end - second_start
            first_shared = 1
            second_shared = 1
        end if
        if (norm2(first_start - second_end) <= tolerance) then
            count = count + 1
            first_away = first_end - first_start
            second_away = second_start - second_end
            first_shared = 1
            second_shared = 2
        end if
        if (norm2(first_end - second_start) <= tolerance) then
            count = count + 1
            first_away = first_start - first_end
            second_away = second_end - second_start
            first_shared = 2
            second_shared = 1
        end if
        if (norm2(first_end - second_end) <= tolerance) then
            count = count + 1
            first_away = first_start - first_end
            second_away = second_start - second_end
            first_shared = 2
            second_shared = 2
        end if
    end subroutine shared_linear_panel_endpoint

    subroutine adjacent_linear_single_layer_moment( &
            first_start, first_end, second_start, second_end, &
            first_away, second_away, first_shared, second_shared, &
            first_length, second_length, wavenumber, nodes, weights, &
            matrix, status)
        real(dp), intent(in) :: first_start(2), first_end(2)
        real(dp), intent(in) :: second_start(2), second_end(2)
        real(dp), intent(in) :: first_away(2), second_away(2)
        integer, intent(in) :: first_shared, second_shared
        real(dp), intent(in) :: first_length, second_length, wavenumber
        real(dp), intent(in) :: nodes(:), weights(:)
        complex(dp), intent(out) :: matrix(2, 2)
        integer, intent(out) :: status

        complex(dp) :: remainder_matrix(2, 2)
        real(dp) :: laplace_matrix(2, 2)

        call adjacent_linear_laplace_moment( &
            first_away, second_away, first_shared, second_shared, &
            first_length, second_length, nodes, weights, laplace_matrix, status)
        if (status /= 0) return
        call regular_linear_remainder_moment( &
            first_start, first_end, second_start, second_end, &
            first_length, second_length, wavenumber, nodes, weights, &
            remainder_matrix, status)
        if (status /= 0) return
        matrix = cmplx(laplace_matrix, 0.0_dp, dp) + remainder_matrix
    end subroutine adjacent_linear_single_layer_moment

    subroutine adjacent_linear_laplace_moment( &
            first_away, second_away, first_shared, second_shared, &
            first_length, second_length, nodes, weights, matrix, status)
        real(dp), intent(in) :: first_away(2), second_away(2)
        integer, intent(in) :: first_shared, second_shared
        real(dp), intent(in) :: first_length, second_length
        real(dp), intent(in) :: nodes(:), weights(:)
        real(dp), intent(out) :: matrix(2, 2)
        integer, intent(out) :: status

        real(dp) :: first_a, first_b, second_a, second_b
        real(dp) :: distance_first, distance_second, integral
        real(dp) :: parameter
        integer :: first_endpoint, node, second_endpoint

        matrix = 0.0_dp
        status = 2
        do node = 1, size(nodes)
            parameter = nodes(node)
            distance_first = norm2(first_away - parameter * second_away)
            distance_second = norm2(parameter * first_away - second_away)
            if (distance_first <= 0.0_dp .or. distance_second <= 0.0_dp) return

            do first_endpoint = 1, 2
                call away_basis_coefficients( &
                    first_endpoint, first_shared, first_a, first_b)
                do second_endpoint = 1, 2
                    call away_basis_coefficients( &
                        second_endpoint, second_shared, second_a, second_b)
                    integral = integrated_duffy_log_product( &
                        first_a, first_b, second_a, &
                        parameter * second_b, log(distance_first))
                    integral = integral + integrated_duffy_log_product( &
                        first_a, parameter * first_b, second_a, &
                        second_b, log(distance_second))
                    matrix(first_endpoint, second_endpoint) = &
                        matrix(first_endpoint, second_endpoint) - &
                        weights(node) * integral * first_length * &
                        second_length / (2.0_dp * acos(-1.0_dp))
                end do
            end do
        end do
        status = 0
    end subroutine adjacent_linear_laplace_moment

    pure subroutine away_basis_coefficients( &
            endpoint, shared_endpoint, constant, slope)
        integer, intent(in) :: endpoint, shared_endpoint
        real(dp), intent(out) :: constant, slope

        if (endpoint == shared_endpoint) then
            constant = 1.0_dp
            slope = -1.0_dp
        else
            constant = 0.0_dp
            slope = 1.0_dp
        end if
    end subroutine away_basis_coefficients

    pure function integrated_duffy_log_product( &
            first_a, first_b, second_a, second_b, log_distance) result(value)
        real(dp), intent(in) :: first_a, first_b, second_a, second_b
        real(dp), intent(in) :: log_distance
        real(dp) :: value

        real(dp) :: linear, quadratic, cubic

        linear = first_a * second_a
        quadratic = first_a * second_b + first_b * second_a
        cubic = first_b * second_b
        value = -linear / 4.0_dp - quadratic / 9.0_dp - cubic / 16.0_dp + &
            log_distance * ( &
            linear / 2.0_dp + quadratic / 3.0_dp + cubic / 4.0_dp)
    end function integrated_duffy_log_product

    subroutine regular_linear_remainder_moment( &
            first_start, first_end, second_start, second_end, &
            first_length, second_length, wavenumber, nodes, weights, &
            matrix, status)
        real(dp), intent(in) :: first_start(2), first_end(2)
        real(dp), intent(in) :: second_start(2), second_end(2)
        real(dp), intent(in) :: first_length, second_length, wavenumber
        real(dp), intent(in) :: nodes(:), weights(:)
        complex(dp), intent(out) :: matrix(2, 2)
        integer, intent(out) :: status

        complex(dp) :: remainder
        real(dp) :: distance, first_basis(2), first_point(2), second_basis(2)
        real(dp) :: second_point(2)
        integer :: first_endpoint, first_node, second_endpoint, second_node

        matrix = (0.0_dp, 0.0_dp)
        status = 2
        do first_node = 1, size(nodes)
            first_point = first_start + &
                nodes(first_node) * (first_end - first_start)
            first_basis(1) = 1.0_dp - nodes(first_node)
            first_basis(2) = nodes(first_node)
            do second_node = 1, size(nodes)
                second_point = second_start + &
                    nodes(second_node) * (second_end - second_start)
                second_basis(1) = 1.0_dp - nodes(second_node)
                second_basis(2) = nodes(second_node)
                distance = norm2(first_point - second_point)
                if (distance <= 0.0_dp) return
                call helmholtz_laplace_remainder( &
                    distance, wavenumber, remainder, status)
                if (status /= 0) return
                do first_endpoint = 1, 2
                    do second_endpoint = 1, 2
                        matrix(first_endpoint, second_endpoint) = &
                            matrix(first_endpoint, second_endpoint) + &
                            weights(first_node) * weights(second_node) * &
                            first_basis(first_endpoint) * &
                            second_basis(second_endpoint) * remainder
                    end do
                end do
            end do
        end do
        matrix = first_length * second_length * matrix
        status = 0
    end subroutine regular_linear_remainder_moment

    subroutine regular_linear_single_layer_moment( &
            first_start, first_end, second_start, second_end, &
            first_length, second_length, wavenumber, nodes, weights, &
            matrix, status)
        real(dp), intent(in) :: first_start(2), first_end(2)
        real(dp), intent(in) :: second_start(2), second_end(2)
        real(dp), intent(in) :: first_length, second_length, wavenumber
        real(dp), intent(in) :: nodes(:), weights(:)
        complex(dp), intent(out) :: matrix(2, 2)
        integer, intent(out) :: status

        complex(dp) :: hankel, kernel
        real(dp) :: distance, first_basis(2), first_point(2), second_basis(2)
        real(dp) :: second_point(2)
        type(fortnum_status_t) :: special_status
        integer :: first_endpoint, first_node, second_endpoint, second_node

        matrix = (0.0_dp, 0.0_dp)
        status = 2
        do first_node = 1, size(nodes)
            first_point = first_start + &
                nodes(first_node) * (first_end - first_start)
            first_basis(1) = 1.0_dp - nodes(first_node)
            first_basis(2) = nodes(first_node)
            do second_node = 1, size(nodes)
                second_point = second_start + &
                    nodes(second_node) * (second_end - second_start)
                second_basis(1) = 1.0_dp - nodes(second_node)
                second_basis(2) = nodes(second_node)
                distance = norm2(first_point - second_point)
                if (distance <= 0.0_dp) return
                call hankel_h1_real( &
                    0, wavenumber * distance, hankel, special_status)
                if (special_status%code /= 0) return
                kernel = cmplx(0.0_dp, 0.25_dp, dp) * hankel
                do first_endpoint = 1, 2
                    do second_endpoint = 1, 2
                        matrix(first_endpoint, second_endpoint) = &
                            matrix(first_endpoint, second_endpoint) + &
                            weights(first_node) * weights(second_node) * &
                            first_basis(first_endpoint) * &
                            second_basis(second_endpoint) * kernel
                    end do
                end do
            end do
        end do
        matrix = first_length * second_length * matrix
        status = 0
    end subroutine regular_linear_single_layer_moment

    subroutine regular_remainder_integral( &
            first_start, first_end, second_start, second_end, &
            first_length, second_length, wavenumber, nodes, weights, &
            integral, status)
        real(dp), intent(in) :: first_start(2), first_end(2)
        real(dp), intent(in) :: second_start(2), second_end(2)
        real(dp), intent(in) :: first_length, second_length, wavenumber
        real(dp), intent(in) :: nodes(:), weights(:)
        complex(dp), intent(out) :: integral
        integer, intent(out) :: status

        complex(dp) :: remainder
        real(dp) :: distance, first_point(2), second_point(2)
        integer :: first_node, second_node

        integral = (0.0_dp, 0.0_dp)
        status = 2
        do first_node = 1, size(nodes)
            first_point = first_start + &
                nodes(first_node) * (first_end - first_start)
            do second_node = 1, size(nodes)
                second_point = second_start + &
                    nodes(second_node) * (second_end - second_start)
                distance = norm2(first_point - second_point)
                if (distance <= 0.0_dp) return
                call helmholtz_laplace_remainder( &
                    distance, wavenumber, remainder, status)
                if (status /= 0) return
                integral = integral + weights(first_node) * &
                    weights(second_node) * remainder
            end do
        end do
        integral = first_length * second_length * integral
        status = 0
    end subroutine regular_remainder_integral

    subroutine regular_double_layer_remainder_integral( &
            first_start, first_end, second_start, second_end, source_normal, &
            first_length, second_length, wavenumber, nodes, weights, &
            integral, status)
        real(dp), intent(in) :: first_start(2), first_end(2)
        real(dp), intent(in) :: second_start(2), second_end(2)
        real(dp), intent(in) :: source_normal(2)
        real(dp), intent(in) :: first_length, second_length, wavenumber
        real(dp), intent(in) :: nodes(:), weights(:)
        complex(dp), intent(out) :: integral
        integer, intent(out) :: status

        complex(dp) :: hankel, helmholtz_kernel
        real(dp) :: displacement(2), distance, laplace_kernel, projection
        real(dp) :: first_point(2), second_point(2)
        type(fortnum_status_t) :: special_status
        integer :: first_node, second_node

        integral = (0.0_dp, 0.0_dp)
        status = 2
        do first_node = 1, size(nodes)
            first_point = first_start + &
                nodes(first_node) * (first_end - first_start)
            do second_node = 1, size(nodes)
                second_point = second_start + &
                    nodes(second_node) * (second_end - second_start)
                displacement = first_point - second_point
                distance = norm2(displacement)
                if (distance <= 0.0_dp) return
                projection = dot_product(displacement, source_normal)
                call hankel_h1_real( &
                    1, wavenumber * distance, hankel, special_status)
                if (special_status%code /= 0) return
                helmholtz_kernel = cmplx(0.0_dp, 0.25_dp, dp) * &
                    wavenumber * hankel * projection / distance
                laplace_kernel = projection / &
                    (2.0_dp * acos(-1.0_dp) * distance**2)
                integral = integral + weights(first_node) * &
                    weights(second_node) * &
                    (helmholtz_kernel - laplace_kernel)
            end do
        end do
        integral = first_length * second_length * integral
        status = 0
    end subroutine regular_double_layer_remainder_integral

    subroutine regular_double_layer_remainder_mixed_moment( &
            first_start, first_end, second_start, second_end, source_normal, &
            first_length, second_length, wavenumber, nodes, weights, &
            integral, status)
        real(dp), intent(in) :: first_start(2), first_end(2)
        real(dp), intent(in) :: second_start(2), second_end(2)
        real(dp), intent(in) :: source_normal(2)
        real(dp), intent(in) :: first_length, second_length, wavenumber
        real(dp), intent(in) :: nodes(:), weights(:)
        complex(dp), intent(out) :: integral(2)
        integer, intent(out) :: status

        complex(dp) :: hankel, kernel
        real(dp) :: basis(2), displacement(2), distance, projection
        real(dp) :: first_point(2), laplace_kernel, second_point(2)
        type(fortnum_status_t) :: special_status
        integer :: endpoint, first_node, second_node

        integral = cmplx(0.0_dp, 0.0_dp, dp)
        status = 2
        do first_node = 1, size(nodes)
            first_point = first_start + &
                nodes(first_node)*(first_end - first_start)
            do second_node = 1, size(nodes)
                second_point = second_start + &
                    nodes(second_node)*(second_end - second_start)
                basis = [1.0_dp - nodes(second_node), nodes(second_node)]
                displacement = first_point - second_point
                distance = norm2(displacement)
                if (distance <= 0.0_dp) return
                projection = dot_product(displacement, source_normal)
                call hankel_h1_real( &
                    1, wavenumber*distance, hankel, special_status)
                if (special_status%code /= 0) return
                kernel = cmplx(0.0_dp, 0.25_dp, dp)*wavenumber*hankel* &
                    projection/distance
                laplace_kernel = projection/( &
                    2.0_dp*acos(-1.0_dp)*distance**2)
                do endpoint = 1, 2
                    integral(endpoint) = integral(endpoint) + &
                        weights(first_node)*weights(second_node)* &
                        basis(endpoint)*(kernel - laplace_kernel)
                end do
            end do
        end do
        integral = first_length*second_length*integral
        status = 0
    end subroutine regular_double_layer_remainder_mixed_moment

    subroutine helmholtz_laplace_remainder( &
            distance, wavenumber, remainder, status)
        real(dp), intent(in) :: distance, wavenumber
        complex(dp), intent(out) :: remainder
        integer, intent(out) :: status

        complex(dp) :: hankel
        type(fortnum_status_t) :: special_status

        remainder = (0.0_dp, 0.0_dp)
        status = 2
        if (distance <= 0.0_dp .or. wavenumber <= 0.0_dp) return
        call hankel_h1_real(0, wavenumber * distance, hankel, special_status)
        if (special_status%code /= 0) return
        remainder = cmplx(0.0_dp, 0.25_dp, dp) * hankel + &
            cmplx(log(distance) / (2.0_dp * acos(-1.0_dp)), 0.0_dp, dp)
        status = 0
    end subroutine helmholtz_laplace_remainder

end module fortfem_helmholtz_boundary_operators_2d
