module fortfem_helmholtz_bem_state_ad_3d
    use fortfem_helmholtz_galerkin_3d, only: &
        assemble_helmholtz_single_layer_p0_3d
    use fortfem_helmholtz_panel_pair_3d, only: &
        assemble_helmholtz_single_layer_p0_3d_jvp, &
        assemble_helmholtz_single_layer_p0_3d_vjp
    use fortfem_kinds, only: dp
    use fortfem_surface_triangle_areas_3d, only: &
        assemble_surface_triangle_areas_3d, &
        assemble_surface_triangle_areas_3d_jvp, &
        assemble_surface_triangle_areas_3d_vjp
    use fortnum_linalg, only: dense_solve, linear_solve_complex_jvp, &
        linear_solve_complex_vjp
    implicit none

    private

    public :: solve_helmholtz_dirichlet_p0_3d_jvp
    public :: solve_helmholtz_dirichlet_p0_3d_vjp

contains

    subroutine solve_helmholtz_dirichlet_p0_3d_jvp( &
            vertices, triangles, boundary_value, wave_number, &
            quadrature_degree, vertices_dot, boundary_value_dot, &
            wave_number_dot, density_dot, status)
        real(dp), intent(in) :: vertices(:, :), wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(in) :: boundary_value, boundary_value_dot
        real(dp), intent(in) :: vertices_dot(:, :), wave_number_dot
        complex(dp), allocatable, intent(out) :: density_dot(:)
        integer, intent(out) :: status

        real(dp), allocatable :: areas(:), areas_dot(:)
        complex(dp), allocatable :: density(:), matrix(:, :), matrix_dot(:, :)
        complex(dp), allocatable :: rhs(:), rhs_dot(:)
        integer :: info

        status = 1
        call assemble_helmholtz_single_layer_p0_3d( &
            vertices, triangles, wave_number, quadrature_degree, matrix, status)
        if (status /= 0) return
        call assemble_helmholtz_single_layer_p0_3d_jvp( &
            vertices, triangles, wave_number, quadrature_degree, vertices_dot, &
            wave_number_dot, matrix_dot, status)
        if (status /= 0) return
        call assemble_surface_triangle_areas_3d_jvp( &
            vertices, triangles, vertices_dot, areas, areas_dot, status)
        if (status /= 0) return
        rhs = boundary_value*areas
        rhs_dot = boundary_value_dot*areas + boundary_value*areas_dot
        allocate(density(size(rhs)), density_dot(size(rhs)))
        call dense_solve(matrix, rhs, density, info)
        if (info /= 0) then
            status = 2
            return
        end if
        call linear_solve_complex_jvp( &
            size(density), matrix, density, matrix_dot, rhs_dot, density_dot, &
            info)
        if (info /= 0) then
            status = 2
            return
        end if
        status = 0
    end subroutine solve_helmholtz_dirichlet_p0_3d_jvp

    subroutine solve_helmholtz_dirichlet_p0_3d_vjp( &
            vertices, triangles, boundary_value, wave_number, &
            quadrature_degree, density_bar, vertices_bar, boundary_value_bar, &
            wave_number_bar, status)
        real(dp), intent(in) :: vertices(:, :), wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(in) :: boundary_value, density_bar(:)
        real(dp), intent(out) :: vertices_bar(:, :), wave_number_bar
        complex(dp), intent(out) :: boundary_value_bar
        integer, intent(out) :: status

        real(dp), allocatable :: areas(:), areas_bar(:)
        complex(dp), allocatable :: b_bar(:), density(:), matrix(:, :)
        complex(dp), allocatable :: matrix_bar(:, :), rhs(:)
        real(dp) :: assembly_vertices_bar(size(vertices, 1), size(vertices, 2))
        real(dp) :: area_vertices_bar(size(vertices, 1), size(vertices, 2))
        integer :: info

        vertices_bar = 0.0_dp
        boundary_value_bar = cmplx(0.0_dp, 0.0_dp, dp)
        wave_number_bar = 0.0_dp
        status = 1
        if (any(shape(vertices_bar) /= shape(vertices))) return
        if (size(density_bar) /= size(triangles, 2)) return
        call assemble_helmholtz_single_layer_p0_3d( &
            vertices, triangles, wave_number, quadrature_degree, matrix, status)
        if (status /= 0) return
        call assemble_surface_triangle_areas_3d( &
            vertices, triangles, areas, status)
        if (status /= 0) return
        rhs = boundary_value*areas
        allocate(density(size(rhs)), b_bar(size(rhs)))
        allocate(matrix_bar(size(matrix, 1), size(matrix, 2)))
        call dense_solve(matrix, rhs, density, info)
        if (info /= 0) then
            status = 2
            return
        end if
        call linear_solve_complex_vjp( &
            size(density), matrix, density, density_bar, matrix_bar, b_bar, &
            info)
        if (info /= 0) then
            status = 2
            return
        end if
        call assemble_helmholtz_single_layer_p0_3d_vjp( &
            vertices, triangles, wave_number, quadrature_degree, matrix_bar, &
            assembly_vertices_bar, wave_number_bar, status)
        if (status /= 0) return
        areas_bar = real(conjg(b_bar)*boundary_value, dp)
        call assemble_surface_triangle_areas_3d_vjp( &
            vertices, triangles, areas_bar, area_vertices_bar, status)
        if (status /= 0) return
        vertices_bar = assembly_vertices_bar + area_vertices_bar
        boundary_value_bar = sum(b_bar*areas)
        status = 0
    end subroutine solve_helmholtz_dirichlet_p0_3d_vjp

end module fortfem_helmholtz_bem_state_ad_3d
