module fortfem_surface_triangle_areas_3d
    use fortfem_kinds, only: dp
    use fortfem_surface_triangle_geometry_3d, only: &
        evaluate_surface_triangle_geometry_3d, &
        evaluate_surface_triangle_geometry_3d_jvp, &
        evaluate_surface_triangle_geometry_3d_vjp
    implicit none

    private

    public :: assemble_surface_triangle_areas_3d
    public :: assemble_surface_triangle_areas_3d_jvp
    public :: assemble_surface_triangle_areas_3d_vjp

contains

    subroutine assemble_surface_triangle_areas_3d( &
            vertices, triangles, areas, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), allocatable, intent(out) :: areas(:)
        integer, intent(out) :: status

        real(dp) :: jacobian, normal(3), point(3)
        integer :: triangle

        status = 1
        if (.not. valid_surface(vertices, triangles)) return
        allocate(areas(size(triangles, 2)))
        do triangle = 1, size(triangles, 2)
            call evaluate_surface_triangle_geometry_3d( &
                vertices(:, triangles(:, triangle)), 0.0_dp, 0.0_dp, point, &
                jacobian, normal, status)
            if (status /= 0) return
            areas(triangle) = 0.5_dp*jacobian
        end do
        status = 0
    end subroutine assemble_surface_triangle_areas_3d

    subroutine assemble_surface_triangle_areas_3d_jvp( &
            vertices, triangles, vertices_dot, areas, areas_dot, status)
        real(dp), intent(in) :: vertices(:, :), vertices_dot(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), allocatable, intent(out) :: areas(:), areas_dot(:)
        integer, intent(out) :: status

        real(dp) :: jacobian, jacobian_dot, normal(3), normal_dot(3)
        real(dp) :: point(3), point_dot(3)
        integer :: triangle

        status = 1
        if (.not. valid_surface(vertices, triangles)) return
        if (any(shape(vertices_dot) /= shape(vertices))) return
        allocate(areas(size(triangles, 2)), areas_dot(size(triangles, 2)))
        do triangle = 1, size(triangles, 2)
            call evaluate_surface_triangle_geometry_3d( &
                vertices(:, triangles(:, triangle)), 0.0_dp, 0.0_dp, point, &
                jacobian, normal, status)
            if (status /= 0) return
            call evaluate_surface_triangle_geometry_3d_jvp( &
                vertices(:, triangles(:, triangle)), 0.0_dp, 0.0_dp, &
                vertices_dot(:, triangles(:, triangle)), 0.0_dp, 0.0_dp, &
                point_dot, jacobian_dot, normal_dot, status)
            if (status /= 0) return
            areas(triangle) = 0.5_dp*jacobian
            areas_dot(triangle) = 0.5_dp*jacobian_dot
        end do
        status = 0
    end subroutine assemble_surface_triangle_areas_3d_jvp

    subroutine assemble_surface_triangle_areas_3d_vjp( &
            vertices, triangles, areas_bar, vertices_bar, status)
        real(dp), intent(in) :: vertices(:, :), areas_bar(:)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(out) :: vertices_bar(:, :)
        integer, intent(out) :: status

        real(dp) :: local_vertices_bar(3, 3), xi_bar, eta_bar
        integer :: local, triangle

        vertices_bar = 0.0_dp
        status = 1
        if (.not. valid_surface(vertices, triangles)) return
        if (size(areas_bar) /= size(triangles, 2)) return
        if (any(shape(vertices_bar) /= shape(vertices))) return
        do triangle = 1, size(triangles, 2)
            call evaluate_surface_triangle_geometry_3d_vjp( &
                vertices(:, triangles(:, triangle)), 0.0_dp, 0.0_dp, &
                [0.0_dp, 0.0_dp, 0.0_dp], 0.5_dp*areas_bar(triangle), &
                [0.0_dp, 0.0_dp, 0.0_dp], local_vertices_bar, xi_bar, &
                eta_bar, status)
            if (status /= 0) return
            do local = 1, 3
                vertices_bar(:, triangles(local, triangle)) = &
                    vertices_bar(:, triangles(local, triangle)) + &
                    local_vertices_bar(:, local)
            end do
        end do
        status = 0
    end subroutine assemble_surface_triangle_areas_3d_vjp

    pure function valid_surface(vertices, triangles) result(valid)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        logical :: valid

        valid = size(vertices, 1) == 3 .and. size(vertices, 2) >= 3 .and. &
            size(triangles, 1) == 3 .and. size(triangles, 2) >= 1
        if (.not. valid) return
        valid = all(triangles >= 1) .and. &
            all(triangles <= size(vertices, 2))
    end function valid_surface

end module fortfem_surface_triangle_areas_3d
