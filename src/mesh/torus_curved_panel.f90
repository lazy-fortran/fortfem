module fortfem_torus_curved_panel
    !! Exact parametric map from a reference triangle in angle space to a torus.
    use fortfem_generated_torus_curved_panel, only: &
        generated_torus_curved_panel
    use fortfem_generated_torus_curved_panel_jvp, only: &
        generated_torus_curved_panel_jvp
    use fortfem_generated_torus_curved_panel_vjp, only: &
        generated_torus_curved_panel_vjp
    use fortfem_kinds, only: dp
    implicit none

    private

    public :: evaluate_torus_curved_panel
    public :: evaluate_torus_curved_panel_jvp
    public :: evaluate_torus_curved_panel_vjp

contains

    pure subroutine evaluate_torus_curved_panel( &
            parameters, major_radius, minor_radius, xi, eta, point, &
            tangent_xi, tangent_eta, surface_jacobian, status)
        real(dp), intent(in) :: parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius, xi, eta
        real(dp), intent(out) :: point(3), tangent_xi(3), tangent_eta(3)
        real(dp), intent(out) :: surface_jacobian
        integer, intent(out) :: status

        real(dp) :: local_parameters(2, 3)

        point = 0.0_dp
        tangent_xi = 0.0_dp
        tangent_eta = 0.0_dp
        surface_jacobian = 0.0_dp
        status = 1
        if (major_radius <= minor_radius) return
        if (minor_radius <= 0.0_dp) return
        if (xi < 0.0_dp) return
        if (eta < 0.0_dp) return
        if (xi + eta > 1.0_dp) return
        local_parameters = parameters
        call unwrap_parameters(local_parameters)
        call generated_torus_curved_panel( &
            local_parameters(1, 1), local_parameters(2, 1), &
            local_parameters(1, 2), local_parameters(2, 2), &
            local_parameters(1, 3), local_parameters(2, 3), &
            major_radius, minor_radius, xi, eta, point(1), point(2), &
            point(3), tangent_xi(1), tangent_xi(2), tangent_xi(3), &
            tangent_eta(1), tangent_eta(2), tangent_eta(3), surface_jacobian)
        if (surface_jacobian <= tiny(1.0_dp)) return
        status = 0
    end subroutine evaluate_torus_curved_panel

    pure subroutine evaluate_torus_curved_panel_jvp( &
            parameters, major_radius, minor_radius, xi, eta, parameters_dot, &
            major_radius_dot, minor_radius_dot, xi_dot, eta_dot, point_dot, &
            tangent_xi_dot, tangent_eta_dot, surface_jacobian_dot, status)
        real(dp), intent(in) :: parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius, xi, eta
        real(dp), intent(in) :: parameters_dot(2, 3)
        real(dp), intent(in) :: major_radius_dot, minor_radius_dot
        real(dp), intent(in) :: xi_dot, eta_dot
        real(dp), intent(out) :: point_dot(3), tangent_xi_dot(3)
        real(dp), intent(out) :: tangent_eta_dot(3), surface_jacobian_dot
        integer, intent(out) :: status

        real(dp) :: local_parameters(2, 3)
        real(dp) :: point(3), tangent_xi(3), tangent_eta(3), jacobian

        point_dot = 0.0_dp
        tangent_xi_dot = 0.0_dp
        tangent_eta_dot = 0.0_dp
        surface_jacobian_dot = 0.0_dp
        call evaluate_torus_curved_panel( &
            parameters, major_radius, minor_radius, xi, eta, point, &
            tangent_xi, tangent_eta, jacobian, status)
        if (status /= 0) return
        local_parameters = parameters
        call unwrap_parameters(local_parameters)
        call generated_torus_curved_panel_jvp( &
            local_parameters(1, 1), local_parameters(2, 1), &
            local_parameters(1, 2), local_parameters(2, 2), &
            local_parameters(1, 3), local_parameters(2, 3), &
            major_radius, minor_radius, xi, eta, &
            parameters_dot(1, 1), parameters_dot(2, 1), &
            parameters_dot(1, 2), parameters_dot(2, 2), &
            parameters_dot(1, 3), parameters_dot(2, 3), &
            major_radius_dot, minor_radius_dot, xi_dot, eta_dot, &
            point_dot(1), point_dot(2), point_dot(3), &
            tangent_xi_dot(1), tangent_xi_dot(2), tangent_xi_dot(3), &
            tangent_eta_dot(1), tangent_eta_dot(2), tangent_eta_dot(3), &
            surface_jacobian_dot)
    end subroutine evaluate_torus_curved_panel_jvp

    pure subroutine evaluate_torus_curved_panel_vjp( &
            parameters, major_radius, minor_radius, xi, eta, point_bar, &
            tangent_xi_bar, tangent_eta_bar, surface_jacobian_bar, &
            parameters_bar, major_radius_bar, minor_radius_bar, xi_bar, &
            eta_bar, status)
        real(dp), intent(in) :: parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius, xi, eta
        real(dp), intent(in) :: point_bar(3), tangent_xi_bar(3)
        real(dp), intent(in) :: tangent_eta_bar(3), surface_jacobian_bar
        real(dp), intent(out) :: parameters_bar(2, 3)
        real(dp), intent(out) :: major_radius_bar, minor_radius_bar
        real(dp), intent(out) :: xi_bar, eta_bar
        integer, intent(out) :: status

        real(dp) :: local_parameters(2, 3)
        real(dp) :: point(3), tangent_xi(3), tangent_eta(3), jacobian

        parameters_bar = 0.0_dp
        major_radius_bar = 0.0_dp
        minor_radius_bar = 0.0_dp
        xi_bar = 0.0_dp
        eta_bar = 0.0_dp
        call evaluate_torus_curved_panel( &
            parameters, major_radius, minor_radius, xi, eta, point, &
            tangent_xi, tangent_eta, jacobian, status)
        if (status /= 0) return
        local_parameters = parameters
        call unwrap_parameters(local_parameters)
        call generated_torus_curved_panel_vjp( &
            local_parameters(1, 1), local_parameters(2, 1), &
            local_parameters(1, 2), local_parameters(2, 2), &
            local_parameters(1, 3), local_parameters(2, 3), &
            major_radius, minor_radius, xi, eta, point_bar(1), point_bar(2), &
            point_bar(3), tangent_xi_bar(1), tangent_xi_bar(2), &
            tangent_xi_bar(3), tangent_eta_bar(1), tangent_eta_bar(2), &
            tangent_eta_bar(3), surface_jacobian_bar, &
            parameters_bar(1, 1), parameters_bar(2, 1), &
            parameters_bar(1, 2), parameters_bar(2, 2), &
            parameters_bar(1, 3), parameters_bar(2, 3), major_radius_bar, &
            minor_radius_bar, xi_bar, eta_bar)
    end subroutine evaluate_torus_curved_panel_vjp

    pure subroutine unwrap_parameters(parameters)
        real(dp), intent(inout) :: parameters(2, 3)

        real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp
        integer :: component, vertex

        do vertex = 2, 3
            do component = 1, 2
                do while (parameters(component, vertex) - &
                        parameters(component, 1) > pi)
                    parameters(component, vertex) = &
                        parameters(component, vertex) - 2.0_dp*pi
                end do
                do while (parameters(component, vertex) - &
                        parameters(component, 1) < -pi)
                    parameters(component, vertex) = &
                        parameters(component, vertex) + 2.0_dp*pi
                end do
            end do
        end do
    end subroutine unwrap_parameters

end module fortfem_torus_curved_panel
