module fortfem_maxwell_efie_rwg_3d
    !! Electric-field integral equation in a lowest-order RWG space.
    !!
    !! The scalar-potential block is a Galerkin product D^T V D because RWG
    !! surface divergences are panelwise constant. The vector-potential block
    !! uses centroid moments of the affine RWG functions with the adaptively
    !! integrated Helmholtz panel operator. Keeping the two blocks separate
    !! exposes the low-frequency scaling and permits compatible preconditioning.
    use fortfem_helmholtz_galerkin_3d, only: &
        assemble_helmholtz_single_layer_p0_adaptive_3d
    use fortfem_kinds, only: dp
    use fortfem_maxwell_rwg_surface, only: &
        build_maxwell_rwg_surface_space, evaluate_maxwell_rwg_basis
    implicit none
    private

    public :: assemble_maxwell_efie_rwg_3d
    public :: assemble_maxwell_rwg_potential_operators_3d

contains

    subroutine assemble_maxwell_efie_rwg_3d( &
            vertices, triangles, wave_number, impedance, quadrature_degree, &
            tolerance, max_depth, matrix, status)
        real(dp), intent(in) :: vertices(:, :), wave_number, impedance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        real(dp), intent(in) :: tolerance
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: scalar_potential(:, :)
        complex(dp), allocatable :: vector_potential(:, :)

        status = 1
        if (wave_number <= 0.0_dp .or. impedance <= 0.0_dp) return
        call assemble_maxwell_rwg_potential_operators_3d( &
            vertices, triangles, wave_number, quadrature_degree, tolerance, &
            max_depth, vector_potential, scalar_potential, status)
        if (status /= 0) return
        matrix = cmplx(0.0_dp, wave_number*impedance, dp)*vector_potential - &
            cmplx(0.0_dp, impedance/wave_number, dp)*scalar_potential
        status = 0
    end subroutine assemble_maxwell_efie_rwg_3d

    subroutine assemble_maxwell_rwg_potential_operators_3d( &
            vertices, triangles, wave_number, quadrature_degree, tolerance, &
            max_depth, vector_potential, scalar_potential, status)
        real(dp), intent(in) :: vertices(:, :), wave_number, tolerance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: vector_potential(:, :)
        complex(dp), allocatable, intent(out) :: scalar_potential(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: panel_operator(:, :)
        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp) :: first_divergence, first_value(3)
        real(dp) :: second_divergence, second_value(3)
        real(dp) :: first_point(3), second_point(3)
        integer :: first, first_panel, first_slot, second, second_panel
        integer :: second_slot

        status = 1
        if (wave_number < 0.0_dp) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0 .or. size(edge_vertices, 2) == 0) return
        call assemble_helmholtz_single_layer_p0_adaptive_3d( &
            vertices, triangles, wave_number, quadrature_degree, tolerance, &
            max_depth, panel_operator, status)
        if (status /= 0) return
        allocate( &
            vector_potential(size(edge_vertices, 2), size(edge_vertices, 2)), &
            scalar_potential(size(edge_vertices, 2), size(edge_vertices, 2)))
        vector_potential = cmplx(0.0_dp, 0.0_dp, dp)
        scalar_potential = cmplx(0.0_dp, 0.0_dp, dp)
        do first = 1, size(edge_vertices, 2)
            do second = 1, first
                do first_slot = 1, 2
                    first_panel = edge_triangles(first_slot, first)
                    first_point = panel_centroid( &
                        vertices(:, triangles(:, first_panel)))
                    call evaluate_maxwell_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        first, first_panel, first_point, first_value, &
                        first_divergence, status)
                    if (status /= 0) return
                    do second_slot = 1, 2
                        second_panel = edge_triangles(second_slot, second)
                        second_point = panel_centroid( &
                            vertices(:, triangles(:, second_panel)))
                        call evaluate_maxwell_rwg_basis( &
                            vertices, triangles, edge_vertices, edge_triangles, &
                            second, second_panel, second_point, second_value, &
                            second_divergence, status)
                        if (status /= 0) return
                        vector_potential(first, second) = &
                            vector_potential(first, second) + &
                            dot_product(first_value, second_value)* &
                            panel_operator(first_panel, second_panel)
                        scalar_potential(first, second) = &
                            scalar_potential(first, second) + &
                            first_divergence*second_divergence* &
                            panel_operator(first_panel, second_panel)
                    end do
                end do
                vector_potential(second, first) = vector_potential(first, second)
                scalar_potential(second, first) = scalar_potential(first, second)
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_rwg_potential_operators_3d

    pure function panel_centroid(vertices) result(point)
        real(dp), intent(in) :: vertices(3, 3)
        real(dp) :: point(3)

        point = sum(vertices, dim=2)/3.0_dp
    end function panel_centroid

end module fortfem_maxwell_efie_rwg_3d
