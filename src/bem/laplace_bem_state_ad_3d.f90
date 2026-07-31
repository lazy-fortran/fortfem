module fortfem_laplace_bem_state_ad_3d
    use fortfem_kinds, only: dp
    use fortfem_laplace_galerkin_3d, only: &
        assemble_laplace_single_layer_p0_3d
    use fortfem_laplace_panel_pair_3d, only: &
        assemble_laplace_single_layer_p0_3d_jvp, &
        assemble_laplace_single_layer_p0_3d_vjp
    use fortfem_surface_triangle_areas_3d, only: &
        assemble_surface_triangle_areas_3d, &
        assemble_surface_triangle_areas_3d_jvp, &
        assemble_surface_triangle_areas_3d_vjp
    use fortnum_linalg, only: dense_solve, linear_solve_jvp, linear_solve_vjp
    implicit none

    private

    public :: solve_laplace_dirichlet_p0_3d_jvp
    public :: solve_laplace_dirichlet_p0_3d_vjp

contains

    subroutine solve_laplace_dirichlet_p0_3d_jvp( &
            vertices, triangles, boundary_value, quadrature_degree, &
            vertices_dot, boundary_value_dot, density_dot, capacity_dot, &
            status)
        real(dp), intent(in) :: vertices(:, :), boundary_value
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), intent(in) :: vertices_dot(:, :), boundary_value_dot
        real(dp), allocatable, intent(out) :: density_dot(:)
        real(dp), intent(out) :: capacity_dot
        integer, intent(out) :: status

        real(dp), allocatable :: areas(:), areas_dot(:), density(:)
        real(dp), allocatable :: matrix(:, :), matrix_dot(:, :), rhs(:)
        real(dp), allocatable :: rhs_dot(:)
        integer :: info

        capacity_dot = 0.0_dp
        status = 1
        call assemble_laplace_single_layer_p0_3d( &
            vertices, triangles, quadrature_degree, matrix, status)
        if (status /= 0) return
        call assemble_laplace_single_layer_p0_3d_jvp( &
            vertices, triangles, quadrature_degree, vertices_dot, matrix_dot, &
            status)
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
        call linear_solve_jvp( &
            size(density), matrix, density, matrix_dot, rhs_dot, density_dot, &
            info)
        if (info /= 0) then
            status = 2
            return
        end if
        capacity_dot = dot_product(density_dot, areas) + &
            dot_product(density, areas_dot)
        status = 0
    end subroutine solve_laplace_dirichlet_p0_3d_jvp

    subroutine solve_laplace_dirichlet_p0_3d_vjp( &
            vertices, triangles, boundary_value, quadrature_degree, &
            density_bar, capacity_bar, vertices_bar, boundary_value_bar, &
            status)
        real(dp), intent(in) :: vertices(:, :), boundary_value
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), intent(in) :: density_bar(:), capacity_bar
        real(dp), intent(out) :: vertices_bar(:, :), boundary_value_bar
        integer, intent(out) :: status

        real(dp), allocatable :: adjoint(:), areas(:), areas_bar(:)
        real(dp), allocatable :: density(:), matrix(:, :), matrix_bar(:, :)
        real(dp), allocatable :: rhs(:), state_bar(:)
        real(dp) :: assembly_vertices_bar(size(vertices, 1), size(vertices, 2))
        real(dp) :: area_vertices_bar(size(vertices, 1), size(vertices, 2))
        integer :: info

        vertices_bar = 0.0_dp
        boundary_value_bar = 0.0_dp
        status = 1
        if (any(shape(vertices_bar) /= shape(vertices))) return
        if (size(density_bar) /= size(triangles, 2)) return
        call assemble_laplace_single_layer_p0_3d( &
            vertices, triangles, quadrature_degree, matrix, status)
        if (status /= 0) return
        call assemble_surface_triangle_areas_3d( &
            vertices, triangles, areas, status)
        if (status /= 0) return
        rhs = boundary_value*areas
        allocate(density(size(rhs)), adjoint(size(rhs)))
        call dense_solve(matrix, rhs, density, info)
        if (info /= 0) then
            status = 2
            return
        end if
        state_bar = density_bar + capacity_bar*areas
        allocate(matrix_bar(size(matrix, 1), size(matrix, 2)))
        call linear_solve_vjp( &
            size(density), matrix, density, state_bar, matrix_bar, adjoint, &
            info)
        if (info /= 0) then
            status = 2
            return
        end if
        call assemble_laplace_single_layer_p0_3d_vjp( &
            vertices, triangles, quadrature_degree, matrix_bar, &
            assembly_vertices_bar, status)
        if (status /= 0) return
        areas_bar = capacity_bar*density + boundary_value*adjoint
        call assemble_surface_triangle_areas_3d_vjp( &
            vertices, triangles, areas_bar, area_vertices_bar, status)
        if (status /= 0) return
        vertices_bar = assembly_vertices_bar + area_vertices_bar
        boundary_value_bar = dot_product(adjoint, areas)
        status = 0
    end subroutine solve_laplace_dirichlet_p0_3d_vjp

end module fortfem_laplace_bem_state_ad_3d
