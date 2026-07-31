module fortfem_laplace_torus_panel_pair_ad_3d
    use fortfem_generated_laplace_single_layer_integrand, only: &
        generated_laplace_single_layer_integrand
    use fortfem_generated_laplace_single_layer_integrand_jvp, only: &
        generated_laplace_single_layer_integrand_jvp
    use fortfem_generated_laplace_single_layer_integrand_vjp, only: &
        generated_laplace_single_layer_integrand_vjp
    use fortfem_kinds, only: dp
    use fortfem_torus_curved_panel, only: &
        evaluate_torus_curved_panel, evaluate_torus_curved_panel_jvp, &
        evaluate_torus_curved_panel_vjp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none

    private

    public :: integrate_laplace_torus_panel_p0_3d
    public :: integrate_laplace_torus_panel_p0_3d_jvp
    public :: integrate_laplace_torus_panel_p0_3d_vjp

contains

    subroutine integrate_laplace_torus_panel_p0_3d( &
            first_parameters, second_parameters, major_radius, minor_radius, &
            quadrature_degree, value, status)
        real(dp), intent(in) :: first_parameters(2, 3), second_parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: quadrature_degree
        real(dp), intent(out) :: value
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: first_jacobian, first_point(3), kernel_scale
        real(dp) :: second_jacobian, second_point(3), term
        real(dp) :: first_tangent_xi(3), first_tangent_eta(3)
        real(dp) :: second_tangent_xi(3), second_tangent_eta(3)
        integer :: first_node, quadrature_status, second_node

        value = 0.0_dp
        status = 1
        if (major_radius <= minor_radius) return
        if (minor_radius <= 0.0_dp) return
        if (all(first_parameters == second_parameters)) return
        if (quadrature_degree < 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        do first_node = 1, size(weights)
            call evaluate_torus_curved_panel( &
                first_parameters, major_radius, minor_radius, xi(first_node), &
                eta(first_node), first_point, first_tangent_xi, &
                first_tangent_eta, first_jacobian, status)
            if (status /= 0) return
            do second_node = 1, size(weights)
                call evaluate_torus_curved_panel( &
                    second_parameters, major_radius, minor_radius, &
                    xi(second_node), eta(second_node), second_point, &
                    second_tangent_xi, second_tangent_eta, second_jacobian, &
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
    end subroutine integrate_laplace_torus_panel_p0_3d

    subroutine integrate_laplace_torus_panel_p0_3d_jvp( &
            first_parameters, second_parameters, major_radius, minor_radius, &
            quadrature_degree, first_parameters_dot, second_parameters_dot, &
            major_radius_dot, minor_radius_dot, value_dot, status)
        real(dp), intent(in) :: first_parameters(2, 3), second_parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: quadrature_degree
        real(dp), intent(in) :: first_parameters_dot(2, 3)
        real(dp), intent(in) :: second_parameters_dot(2, 3)
        real(dp), intent(in) :: major_radius_dot, minor_radius_dot
        real(dp), intent(out) :: value_dot
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: first_jacobian, first_jacobian_dot, first_point(3)
        real(dp) :: first_point_dot(3), kernel_scale, second_jacobian
        real(dp) :: second_jacobian_dot, second_point(3), second_point_dot(3)
        real(dp) :: first_tangent_xi(3), first_tangent_eta(3)
        real(dp) :: first_tangent_xi_dot(3), first_tangent_eta_dot(3)
        real(dp) :: second_tangent_xi(3), second_tangent_eta(3)
        real(dp) :: second_tangent_xi_dot(3), second_tangent_eta_dot(3)
        real(dp) :: term_dot
        integer :: first_node, quadrature_status, second_node

        value_dot = 0.0_dp
        status = 1
        if (major_radius <= minor_radius) return
        if (minor_radius <= 0.0_dp) return
        if (all(first_parameters == second_parameters)) return
        if (quadrature_degree < 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        do first_node = 1, size(weights)
            call evaluate_torus_curved_panel( &
                first_parameters, major_radius, minor_radius, xi(first_node), &
                eta(first_node), first_point, first_tangent_xi, &
                first_tangent_eta, first_jacobian, status)
            if (status /= 0) return
            call evaluate_torus_curved_panel_jvp( &
                first_parameters, major_radius, minor_radius, xi(first_node), &
                eta(first_node), first_parameters_dot, major_radius_dot, &
                minor_radius_dot, 0.0_dp, 0.0_dp, first_point_dot, &
                first_tangent_xi_dot, first_tangent_eta_dot, &
                first_jacobian_dot, status)
            if (status /= 0) return
            do second_node = 1, size(weights)
                call evaluate_torus_curved_panel( &
                    second_parameters, major_radius, minor_radius, &
                    xi(second_node), eta(second_node), second_point, &
                    second_tangent_xi, second_tangent_eta, second_jacobian, &
                    status)
                if (status /= 0) return
                call evaluate_torus_curved_panel_jvp( &
                    second_parameters, major_radius, minor_radius, &
                    xi(second_node), eta(second_node), second_parameters_dot, &
                    major_radius_dot, minor_radius_dot, 0.0_dp, 0.0_dp, &
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
    end subroutine integrate_laplace_torus_panel_p0_3d_jvp

    subroutine integrate_laplace_torus_panel_p0_3d_vjp( &
            first_parameters, second_parameters, major_radius, minor_radius, &
            quadrature_degree, value_bar, first_parameters_bar, &
            second_parameters_bar, major_radius_bar, minor_radius_bar, status)
        real(dp), intent(in) :: first_parameters(2, 3), second_parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius, value_bar
        integer, intent(in) :: quadrature_degree
        real(dp), intent(out) :: first_parameters_bar(2, 3)
        real(dp), intent(out) :: second_parameters_bar(2, 3)
        real(dp), intent(out) :: major_radius_bar, minor_radius_bar
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), first_jacobian_bar(:)
        real(dp), allocatable :: first_point_bar(:, :), second_jacobian_bar(:)
        real(dp), allocatable :: second_point_bar(:, :), weights(:), xi(:)
        real(dp) :: first_jacobian, first_point(3), kernel_scale
        real(dp) :: second_jacobian, second_point(3)
        real(dp) :: local_first_jacobian_bar, local_major_radius_bar
        real(dp) :: local_first_point_bar(3), local_first_parameters_bar(2, 3)
        real(dp) :: local_minor_radius_bar, local_second_jacobian_bar
        real(dp) :: local_second_point_bar(3), local_second_parameters_bar(2, 3)
        real(dp) :: first_tangent_xi(3), first_tangent_eta(3)
        real(dp) :: second_tangent_xi(3), second_tangent_eta(3)
        real(dp) :: xi_bar, eta_bar
        integer :: first_node, quadrature_status, second_node

        first_parameters_bar = 0.0_dp
        second_parameters_bar = 0.0_dp
        major_radius_bar = 0.0_dp
        minor_radius_bar = 0.0_dp
        status = 1
        if (major_radius <= minor_radius) return
        if (minor_radius <= 0.0_dp) return
        if (all(first_parameters == second_parameters)) return
        if (quadrature_degree < 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        allocate(first_point_bar(3, size(weights)), source=0.0_dp)
        allocate(second_point_bar(3, size(weights)), source=0.0_dp)
        allocate(first_jacobian_bar(size(weights)), source=0.0_dp)
        allocate(second_jacobian_bar(size(weights)), source=0.0_dp)
        do first_node = 1, size(weights)
            call evaluate_torus_curved_panel( &
                first_parameters, major_radius, minor_radius, xi(first_node), &
                eta(first_node), first_point, first_tangent_xi, &
                first_tangent_eta, first_jacobian, status)
            if (status /= 0) return
            do second_node = 1, size(weights)
                call evaluate_torus_curved_panel( &
                    second_parameters, major_radius, minor_radius, &
                    xi(second_node), eta(second_node), second_point, &
                    second_tangent_xi, second_tangent_eta, second_jacobian, &
                    status)
                if (status /= 0) return
                if (norm2(first_point - second_point) <= &
                    64.0_dp*epsilon(1.0_dp)) return
                kernel_scale = weights(first_node)*weights(second_node)/ &
                    (4.0_dp*acos(-1.0_dp))
                call generated_laplace_single_layer_integrand_vjp( &
                    first_point(1), first_point(2), first_point(3), &
                    second_point(1), second_point(2), second_point(3), &
                    first_jacobian, second_jacobian, kernel_scale, value_bar, &
                    local_first_point_bar(1), local_first_point_bar(2), &
                    local_first_point_bar(3), local_second_point_bar(1), &
                    local_second_point_bar(2), local_second_point_bar(3), &
                    local_first_jacobian_bar, local_second_jacobian_bar)
                first_point_bar(:, first_node) = &
                    first_point_bar(:, first_node) + local_first_point_bar
                second_point_bar(:, second_node) = &
                    second_point_bar(:, second_node) + local_second_point_bar
                first_jacobian_bar(first_node) = &
                    first_jacobian_bar(first_node) + local_first_jacobian_bar
                second_jacobian_bar(second_node) = &
                    second_jacobian_bar(second_node) + local_second_jacobian_bar
            end do
        end do
        do first_node = 1, size(weights)
            call evaluate_torus_curved_panel_vjp( &
                first_parameters, major_radius, minor_radius, xi(first_node), &
                eta(first_node), first_point_bar(:, first_node), &
                [0.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 0.0_dp, 0.0_dp], &
                first_jacobian_bar(first_node), local_first_parameters_bar, &
                local_major_radius_bar, local_minor_radius_bar, xi_bar, &
                eta_bar, status)
            if (status /= 0) return
            first_parameters_bar = first_parameters_bar + &
                local_first_parameters_bar
            major_radius_bar = major_radius_bar + local_major_radius_bar
            minor_radius_bar = minor_radius_bar + local_minor_radius_bar
        end do
        do second_node = 1, size(weights)
            call evaluate_torus_curved_panel_vjp( &
                second_parameters, major_radius, minor_radius, xi(second_node), &
                eta(second_node), second_point_bar(:, second_node), &
                [0.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 0.0_dp, 0.0_dp], &
                second_jacobian_bar(second_node), local_second_parameters_bar, &
                local_major_radius_bar, local_minor_radius_bar, xi_bar, &
                eta_bar, status)
            if (status /= 0) return
            second_parameters_bar = second_parameters_bar + &
                local_second_parameters_bar
            major_radius_bar = major_radius_bar + local_major_radius_bar
            minor_radius_bar = minor_radius_bar + local_minor_radius_bar
        end do
        status = 0
    end subroutine integrate_laplace_torus_panel_p0_3d_vjp

end module fortfem_laplace_torus_panel_pair_ad_3d
