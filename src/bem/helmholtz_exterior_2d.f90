module fortfem_helmholtz_exterior_2d
    use fortfem_kinds, only: dp
    use fortfem_helmholtz_boundary_operators_2d, only: &
        assemble_helmholtz_double_layer_constant, &
        assemble_helmholtz_single_layer_constant
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_special_complex_bessel, only: hankel_h1_real
    use fortnum_status, only: fortnum_status_t
    implicit none

    private

    public :: evaluate_helmholtz_combined_potential_constant
    public :: solve_helmholtz_cfie_constant

    interface
        subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            complex(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*)
            complex(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine zgesv
    end interface

contains

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
        complex(dp), allocatable :: solution_matrix(:, :)
        real(dp), allocatable :: lengths(:)
        integer, allocatable :: pivots(:)
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
        allocate(solution_matrix(panel_count, 1))
        allocate(pivots(panel_count))
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
        solution_matrix(:, 1) = right_hand_side
        call zgesv( &
            panel_count, 1, matrix, panel_count, pivots, solution_matrix, &
            panel_count, info)
        if (info /= 0) then
            status = 2
            return
        end if
        density = solution_matrix(:, 1)
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
