module fortfem_helmholtz_boundary_operators_2d
    use fortfem_kinds, only: dp
    use fortfem_laplace_boundary_operators_2d, only: &
        assemble_laplace_double_layer_constant, &
        assemble_laplace_single_layer_constant
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_special_complex_bessel, only: hankel_h1_real
    use fortnum_status, only: fortnum_status_t
    implicit none

    private

    public :: assemble_helmholtz_adjoint_double_layer_constant
    public :: assemble_helmholtz_double_layer_constant
    public :: assemble_helmholtz_single_layer_constant

contains

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
