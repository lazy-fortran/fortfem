module fortfem_maxwell_efie_bc_3d
    !! Electric-field operators with Buffa--Christiansen trial and test traces.
    use fortfem_helmholtz_galerkin_3d, only: &
        assemble_helmholtz_single_layer_p0_adaptive_3d
    use fortfem_kinds, only: dp
    use fortfem_maxwell_bc_surface, only: build_maxwell_bc_transformation
    use fortfem_maxwell_localized_rwg_surface, only: &
        evaluate_maxwell_localized_rwg_basis
    implicit none
    private

    public :: assemble_maxwell_bc_scalar_potential_3d
    public :: build_maxwell_bc_panel_divergence

contains

    subroutine assemble_maxwell_bc_scalar_potential_3d( &
            vertices, triangles, wave_number, quadrature_degree, tolerance, &
            max_depth, matrix, status)
        real(dp), intent(in) :: vertices(:, :), wave_number, tolerance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: panel_operator(:, :)
        real(dp), allocatable :: divergence(:, :)
        real(dp), allocatable :: refined_vertices(:, :), transformation(:, :)
        integer, allocatable :: refined_triangles(:, :)

        status = 1
        if (wave_number < 0.0_dp) return
        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status)
        if (status /= 0) return
        call build_maxwell_bc_panel_divergence( &
            vertices, triangles, divergence, status)
        if (status /= 0) return
        if (size(divergence, 2) /= size(transformation, 2)) return
        call assemble_helmholtz_single_layer_p0_adaptive_3d( &
            refined_vertices, refined_triangles, wave_number, &
            quadrature_degree, tolerance, max_depth, panel_operator, status)
        if (status /= 0) return
        matrix = matmul( &
            transpose(cmplx(divergence, 0.0_dp, dp)), &
            matmul(panel_operator, cmplx(divergence, 0.0_dp, dp)))
        status = 0
    end subroutine assemble_maxwell_bc_scalar_potential_3d

    subroutine build_maxwell_bc_panel_divergence( &
            vertices, triangles, divergence, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), allocatable, intent(out) :: divergence(:, :)
        integer, intent(out) :: status

        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: refined_vertices(:, :)
        real(dp), allocatable :: transformation(:, :)
        real(dp) :: local_divergence, local_value(3), point(3)
        integer :: basis, local_edge, refined_panel, row

        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status)
        if (status /= 0) return
        allocate(divergence( &
            size(refined_triangles, 2), size(transformation, 2)))
        divergence = 0.0_dp
        do refined_panel = 1, size(refined_triangles, 2)
            point = sum( &
                refined_vertices(:, refined_triangles(:, refined_panel)), &
                dim=2)/3.0_dp
            do local_edge = 1, 3
                call evaluate_maxwell_localized_rwg_basis( &
                    refined_vertices, refined_triangles, refined_panel, &
                    local_edge, point, local_value, local_divergence, status)
                if (status /= 0) return
                row = 3*(refined_panel - 1) + local_edge
                do basis = 1, size(transformation, 2)
                    divergence(refined_panel, basis) = &
                        divergence(refined_panel, basis) + &
                        transformation(row, basis)*local_divergence
                end do
            end do
        end do
        status = 0
    end subroutine build_maxwell_bc_panel_divergence

end module fortfem_maxwell_efie_bc_3d
