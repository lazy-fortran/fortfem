module fortfem_laplace_sphere_panel_pair_ad_3d
    use fortfem_generated_laplace_single_layer_integrand, only: &
        generated_laplace_single_layer_integrand
    use fortfem_generated_laplace_single_layer_integrand_jvp, only: &
        generated_laplace_single_layer_integrand_jvp
    use fortfem_generated_laplace_single_layer_integrand_vjp, only: &
        generated_laplace_single_layer_integrand_vjp
    use fortfem_kinds, only: dp
    use fortfem_sphere_curved_panel, only: &
        evaluate_sphere_curved_panel, &
        evaluate_sphere_curved_panel_jvp, &
        evaluate_sphere_curved_panel_vjp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    private

    public :: integrate_laplace_sphere_panel_p0_3d
    public :: integrate_laplace_sphere_panel_p0_3d_jvp
    public :: integrate_laplace_sphere_panel_p0_3d_vjp

contains

    subroutine integrate_laplace_sphere_panel_p0_3d( &
            first_vertices, second_vertices, radius, quadrature_degree, &
            value, status)
        real(dp), intent(in) :: first_vertices(3, 3), second_vertices(3, 3)
        real(dp), intent(in) :: radius
        integer, intent(in) :: quadrature_degree
        real(dp), intent(out) :: value
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: first_jacobian, first_point(3), kernel_scale
        real(dp) :: second_jacobian, second_point(3), term
        real(dp) :: tangent_eta(3), tangent_xi(3)
        integer :: first_node, quadrature_status, second_node

        value = 0.0_dp
        status = 1
        if (radius <= 0.0_dp .or. quadrature_degree < 0) return
        if (all(first_vertices == second_vertices)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        do first_node = 1, size(weights)
            call evaluate_sphere_curved_panel( &
                first_vertices, radius, xi(first_node), eta(first_node), &
                first_point, tangent_xi, tangent_eta, first_jacobian, status)
            if (status /= 0) return
            do second_node = 1, size(weights)
                call evaluate_sphere_curved_panel( &
                    second_vertices, radius, xi(second_node), eta(second_node), &
                    second_point, tangent_xi, tangent_eta, second_jacobian, &
                    status)
                if (status /= 0) return
                if (norm2(first_point - second_point) <= &
                    64.0_dp*epsilon(1.0_dp)) return
                kernel_scale = weights(first_node)*weights(second_node)/ &
                    (4.0_dp*acos(-1.0_dp))
                call generated_laplace_single_layer_integrand( &
                    first_point(1), first_point(2), first_point(3), &
                    second_point(1), second_point(2), second_point(3), &
                    first_jacobian, second_jacobian, kernel_scale, term)
                value = value + term
            end do
        end do
        status = 0
    end subroutine integrate_laplace_sphere_panel_p0_3d

    subroutine integrate_laplace_sphere_panel_p0_3d_jvp( &
            first_vertices, second_vertices, radius, quadrature_degree, &
            first_vertices_dot, second_vertices_dot, radius_dot, value_dot, &
            status)
        real(dp), intent(in) :: first_vertices(3, 3), second_vertices(3, 3)
        real(dp), intent(in) :: radius
        integer, intent(in) :: quadrature_degree
        real(dp), intent(in) :: first_vertices_dot(3, 3)
        real(dp), intent(in) :: second_vertices_dot(3, 3), radius_dot
        real(dp), intent(out) :: value_dot
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: first_jacobian, first_jacobian_dot, first_point(3)
        real(dp) :: first_point_dot(3), first_tangent_eta(3)
        real(dp) :: first_tangent_eta_dot(3), first_tangent_xi(3)
        real(dp) :: first_tangent_xi_dot(3), kernel_scale
        real(dp) :: second_jacobian, second_jacobian_dot, second_point(3)
        real(dp) :: second_point_dot(3), second_tangent_eta(3)
        real(dp) :: second_tangent_eta_dot(3), second_tangent_xi(3)
        real(dp) :: second_tangent_xi_dot(3), term_dot
        integer :: first_node, quadrature_status, second_node

        value_dot = 0.0_dp
        status = 1
        if (radius <= 0.0_dp .or. quadrature_degree < 0) return
        if (all(first_vertices == second_vertices)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        do first_node = 1, size(weights)
            call evaluate_sphere_curved_panel( &
                first_vertices, radius, xi(first_node), eta(first_node), &
                first_point, first_tangent_xi, first_tangent_eta, &
                first_jacobian, status)
            if (status /= 0) return
            call evaluate_sphere_curved_panel_jvp( &
                first_vertices, radius, xi(first_node), eta(first_node), &
                first_vertices_dot, radius_dot, 0.0_dp, 0.0_dp, &
                first_point_dot, first_tangent_xi_dot, first_tangent_eta_dot, &
                first_jacobian_dot, status)
            if (status /= 0) return
            do second_node = 1, size(weights)
                call evaluate_sphere_curved_panel( &
                    second_vertices, radius, xi(second_node), eta(second_node), &
                    second_point, second_tangent_xi, second_tangent_eta, &
                    second_jacobian, status)
                if (status /= 0) return
                call evaluate_sphere_curved_panel_jvp( &
                    second_vertices, radius, xi(second_node), eta(second_node), &
                    second_vertices_dot, radius_dot, 0.0_dp, 0.0_dp, &
                    second_point_dot, second_tangent_xi_dot, &
                    second_tangent_eta_dot, second_jacobian_dot, status)
                if (status /= 0) return
                if (norm2(first_point - second_point) <= &
                    64.0_dp*epsilon(1.0_dp)) return
                kernel_scale = weights(first_node)*weights(second_node)/ &
                    (4.0_dp*acos(-1.0_dp))
                call generated_laplace_single_layer_integrand_jvp( &
                    first_point(1), first_point(2), first_point(3), &
                    second_point(1), second_point(2), second_point(3), &
                    first_jacobian, second_jacobian, kernel_scale, &
                    first_point_dot(1), first_point_dot(2), first_point_dot(3), &
                    second_point_dot(1), second_point_dot(2), &
                    second_point_dot(3), first_jacobian_dot, &
                    second_jacobian_dot, term_dot)
                value_dot = value_dot + term_dot
            end do
        end do
        status = 0
    end subroutine integrate_laplace_sphere_panel_p0_3d_jvp

    subroutine integrate_laplace_sphere_panel_p0_3d_vjp( &
            first_vertices, second_vertices, radius, quadrature_degree, &
            value_bar, first_vertices_bar, second_vertices_bar, radius_bar, &
            status)
        real(dp), intent(in) :: first_vertices(3, 3), second_vertices(3, 3)
        real(dp), intent(in) :: radius
        integer, intent(in) :: quadrature_degree
        real(dp), intent(in) :: value_bar
        real(dp), intent(out) :: first_vertices_bar(3, 3)
        real(dp), intent(out) :: second_vertices_bar(3, 3), radius_bar
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: first_jacobian, first_jacobian_bar, first_point(3)
        real(dp) :: first_point_bar(3), first_tangent_eta(3)
        real(dp) :: first_tangent_xi(3), second_jacobian
        real(dp) :: second_jacobian_bar, second_point(3)
        real(dp) :: second_point_bar(3), second_tangent_eta(3)
        real(dp) :: second_tangent_xi(3), kernel_scale
        real(dp) :: first_local_bar(3, 3), second_local_bar(3, 3)
        real(dp) :: first_radius_bar, second_radius_bar
        real(dp) :: xi_bar, eta_bar, term
        integer :: first_node, quadrature_status, second_node

        first_vertices_bar = 0.0_dp
        second_vertices_bar = 0.0_dp
        radius_bar = 0.0_dp
        status = 1
        if (radius <= 0.0_dp .or. quadrature_degree < 0) return
        if (all(first_vertices == second_vertices)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        do first_node = 1, size(weights)
            call evaluate_sphere_curved_panel( &
                first_vertices, radius, xi(first_node), eta(first_node), &
                first_point, first_tangent_xi, first_tangent_eta, &
                first_jacobian, status)
            if (status /= 0) return
            do second_node = 1, size(weights)
                call evaluate_sphere_curved_panel( &
                    second_vertices, radius, xi(second_node), eta(second_node), &
                    second_point, second_tangent_xi, second_tangent_eta, &
                    second_jacobian, status)
                if (status /= 0) return
                if (norm2(first_point - second_point) <= &
                    64.0_dp*epsilon(1.0_dp)) return
                kernel_scale = weights(first_node)*weights(second_node)/ &
                    (4.0_dp*acos(-1.0_dp))
                call generated_laplace_single_layer_integrand( &
                    first_point(1), first_point(2), first_point(3), &
                    second_point(1), second_point(2), second_point(3), &
                    first_jacobian, second_jacobian, kernel_scale, term)
                call generated_laplace_single_layer_integrand_vjp( &
                    first_point(1), first_point(2), first_point(3), &
                    second_point(1), second_point(2), second_point(3), &
                    first_jacobian, second_jacobian, kernel_scale, value_bar, &
                    first_point_bar(1), first_point_bar(2), first_point_bar(3), &
                    second_point_bar(1), second_point_bar(2), &
                    second_point_bar(3), first_jacobian_bar, &
                    second_jacobian_bar)
                call evaluate_sphere_curved_panel_vjp( &
                    first_vertices, radius, xi(first_node), eta(first_node), &
                    first_point_bar, [0.0_dp, 0.0_dp, 0.0_dp], &
                    [0.0_dp, 0.0_dp, 0.0_dp], first_jacobian_bar, &
                    first_local_bar, first_radius_bar, xi_bar, eta_bar, status)
                if (status /= 0) return
                call evaluate_sphere_curved_panel_vjp( &
                    second_vertices, radius, xi(second_node), eta(second_node), &
                    second_point_bar, [0.0_dp, 0.0_dp, 0.0_dp], &
                    [0.0_dp, 0.0_dp, 0.0_dp], second_jacobian_bar, &
                    second_local_bar, second_radius_bar, xi_bar, eta_bar, status)
                if (status /= 0) return
                first_vertices_bar = first_vertices_bar + first_local_bar
                second_vertices_bar = second_vertices_bar + second_local_bar
                radius_bar = radius_bar + first_radius_bar + second_radius_bar
            end do
        end do
        status = 0
    end subroutine integrate_laplace_sphere_panel_p0_3d_vjp

end module fortfem_laplace_sphere_panel_pair_ad_3d
