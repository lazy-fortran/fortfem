module fortfem_torus_curved_panel
    !! Exact parametric map from a reference triangle in angle space to a torus.
    use fortfem_kinds, only: dp
    implicit none

    private

    public :: evaluate_torus_curved_panel

contains

    pure subroutine evaluate_torus_curved_panel( &
            parameters, major_radius, minor_radius, xi, eta, point, &
            tangent_xi, tangent_eta, surface_jacobian, status)
        real(dp), intent(in) :: parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius, xi, eta
        real(dp), intent(out) :: point(3), tangent_xi(3), tangent_eta(3)
        real(dp), intent(out) :: surface_jacobian
        integer, intent(out) :: status

        real(dp) :: local_parameters(2, 3), parameter_eta(2)
        real(dp) :: parameter_xi(2), phi, theta
        real(dp) :: derivative_phi(3), derivative_theta(3)
        integer :: vertex

        point = 0.0_dp
        tangent_xi = 0.0_dp
        tangent_eta = 0.0_dp
        surface_jacobian = 0.0_dp
        status = 1
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        if (xi < 0.0_dp .or. eta < 0.0_dp .or. xi + eta > 1.0_dp) return
        local_parameters = parameters
        do vertex = 2, 3
            do while ( &
                    local_parameters(1, vertex) - local_parameters(1, 1) > &
                    acos(-1.0_dp))
                local_parameters(1, vertex) = &
                    local_parameters(1, vertex) - 2.0_dp*acos(-1.0_dp)
            end do
            do while ( &
                    local_parameters(1, vertex) - local_parameters(1, 1) < &
                    -acos(-1.0_dp))
                local_parameters(1, vertex) = &
                    local_parameters(1, vertex) + 2.0_dp*acos(-1.0_dp)
            end do
            do while ( &
                    local_parameters(2, vertex) - local_parameters(2, 1) > &
                    acos(-1.0_dp))
                local_parameters(2, vertex) = &
                    local_parameters(2, vertex) - 2.0_dp*acos(-1.0_dp)
            end do
            do while ( &
                    local_parameters(2, vertex) - local_parameters(2, 1) < &
                    -acos(-1.0_dp))
                local_parameters(2, vertex) = &
                    local_parameters(2, vertex) + 2.0_dp*acos(-1.0_dp)
            end do
        end do
        parameter_xi = local_parameters(:, 2) - local_parameters(:, 1)
        parameter_eta = local_parameters(:, 3) - local_parameters(:, 1)
        theta = local_parameters(1, 1) + &
            xi*parameter_xi(1) + eta*parameter_eta(1)
        phi = local_parameters(2, 1) + &
            xi*parameter_xi(2) + eta*parameter_eta(2)
        point = [ &
            (major_radius + minor_radius*cos(theta))*cos(phi), &
            (major_radius + minor_radius*cos(theta))*sin(phi), &
            minor_radius*sin(theta)]
        derivative_theta = [ &
            -minor_radius*sin(theta)*cos(phi), &
            -minor_radius*sin(theta)*sin(phi), minor_radius*cos(theta)]
        derivative_phi = [ &
            -(major_radius + minor_radius*cos(theta))*sin(phi), &
            (major_radius + minor_radius*cos(theta))*cos(phi), 0.0_dp]
        tangent_xi = &
            parameter_xi(1)*derivative_theta + &
            parameter_xi(2)*derivative_phi
        tangent_eta = &
            parameter_eta(1)*derivative_theta + &
            parameter_eta(2)*derivative_phi
        surface_jacobian = norm2(cross_product(tangent_xi, tangent_eta))
        if (surface_jacobian <= tiny(1.0_dp)) return
        status = 0
    end subroutine evaluate_torus_curved_panel

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_torus_curved_panel
