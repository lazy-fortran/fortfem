module fortfem_surface_triangle_geometry_3d
    use fortfem_generated_surface_triangle_geometry_3d, only: &
        generated_surface_triangle_geometry_3d
    use fortfem_generated_surface_triangle_geometry_3d_jvp, only: &
        generated_surface_triangle_geometry_3d_jvp
    use fortfem_generated_surface_triangle_geometry_3d_vjp, only: &
        generated_surface_triangle_geometry_3d_vjp
    use fortfem_kinds, only: dp
    implicit none

    private

    public :: evaluate_surface_triangle_geometry_3d
    public :: evaluate_surface_triangle_geometry_3d_jvp
    public :: evaluate_surface_triangle_geometry_3d_vjp

contains

    pure subroutine evaluate_surface_triangle_geometry_3d( &
            vertices, xi, eta, point, jacobian, normal, status)
        real(dp), intent(in) :: vertices(3, 3), xi, eta
        real(dp), intent(out) :: point(3), jacobian, normal(3)
        integer, intent(out) :: status

        point = 0.0_dp
        jacobian = surface_jacobian(vertices)
        normal = 0.0_dp
        status = 1
        if (jacobian <= 64.0_dp*epsilon(1.0_dp)) return
        call generated_surface_triangle_geometry_3d( &
            vertices(1, 1), vertices(2, 1), vertices(3, 1), &
            vertices(1, 2), vertices(2, 2), vertices(3, 2), &
            vertices(1, 3), vertices(2, 3), vertices(3, 3), xi, eta, &
            point(1), point(2), point(3), jacobian, normal(1), normal(2), &
            normal(3))
        status = 0
    end subroutine evaluate_surface_triangle_geometry_3d

    pure subroutine evaluate_surface_triangle_geometry_3d_jvp( &
            vertices, xi, eta, vertices_dot, xi_dot, eta_dot, point_dot, &
            jacobian_dot, normal_dot, status)
        real(dp), intent(in) :: vertices(3, 3), xi, eta
        real(dp), intent(in) :: vertices_dot(3, 3), xi_dot, eta_dot
        real(dp), intent(out) :: point_dot(3), jacobian_dot, normal_dot(3)
        integer, intent(out) :: status

        point_dot = 0.0_dp
        jacobian_dot = 0.0_dp
        normal_dot = 0.0_dp
        status = 1
        if (surface_jacobian(vertices) <= 64.0_dp*epsilon(1.0_dp)) return
        call generated_surface_triangle_geometry_3d_jvp( &
            vertices(1, 1), vertices(2, 1), vertices(3, 1), &
            vertices(1, 2), vertices(2, 2), vertices(3, 2), &
            vertices(1, 3), vertices(2, 3), vertices(3, 3), xi, eta, &
            vertices_dot(1, 1), vertices_dot(2, 1), vertices_dot(3, 1), &
            vertices_dot(1, 2), vertices_dot(2, 2), vertices_dot(3, 2), &
            vertices_dot(1, 3), vertices_dot(2, 3), vertices_dot(3, 3), &
            xi_dot, eta_dot, point_dot(1), point_dot(2), point_dot(3), &
            jacobian_dot, normal_dot(1), normal_dot(2), normal_dot(3))
        status = 0
    end subroutine evaluate_surface_triangle_geometry_3d_jvp

    pure subroutine evaluate_surface_triangle_geometry_3d_vjp( &
            vertices, xi, eta, point_bar, jacobian_bar, normal_bar, &
            vertices_bar, xi_bar, eta_bar, status)
        real(dp), intent(in) :: vertices(3, 3), xi, eta
        real(dp), intent(in) :: point_bar(3), jacobian_bar, normal_bar(3)
        real(dp), intent(out) :: vertices_bar(3, 3), xi_bar, eta_bar
        integer, intent(out) :: status

        vertices_bar = 0.0_dp
        xi_bar = 0.0_dp
        eta_bar = 0.0_dp
        status = 1
        if (surface_jacobian(vertices) <= 64.0_dp*epsilon(1.0_dp)) return
        call generated_surface_triangle_geometry_3d_vjp( &
            vertices(1, 1), vertices(2, 1), vertices(3, 1), &
            vertices(1, 2), vertices(2, 2), vertices(3, 2), &
            vertices(1, 3), vertices(2, 3), vertices(3, 3), xi, eta, &
            point_bar(1), point_bar(2), point_bar(3), jacobian_bar, &
            normal_bar(1), normal_bar(2), normal_bar(3), &
            vertices_bar(1, 1), vertices_bar(2, 1), vertices_bar(3, 1), &
            vertices_bar(1, 2), vertices_bar(2, 2), vertices_bar(3, 2), &
            vertices_bar(1, 3), vertices_bar(2, 3), vertices_bar(3, 3), &
            xi_bar, eta_bar)
        status = 0
    end subroutine evaluate_surface_triangle_geometry_3d_vjp

    pure function surface_jacobian(vertices) result(jacobian)
        real(dp), intent(in) :: vertices(3, 3)
        real(dp) :: jacobian

        real(dp) :: first(3), second(3), area_vector(3)

        first = vertices(:, 2) - vertices(:, 1)
        second = vertices(:, 3) - vertices(:, 1)
        area_vector = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
        jacobian = norm2(area_vector)
    end function surface_jacobian

end module fortfem_surface_triangle_geometry_3d
