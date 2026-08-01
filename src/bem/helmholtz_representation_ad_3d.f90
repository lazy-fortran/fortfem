module fortfem_helmholtz_representation_ad_3d
    use fortfem_kinds, only: dp
    use fortfem_helmholtz_representation_3d, only: &
        evaluate_helmholtz_representation_torus_curved_3d
    use fortfem_torus_curved_panel, only: &
        evaluate_torus_curved_panel, evaluate_torus_curved_panel_jvp, &
        evaluate_torus_curved_panel_vjp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    private

    public :: evaluate_helmholtz_representation_torus_curved_3d_geometry_jvp
    public :: evaluate_helmholtz_representation_torus_curved_3d_geometry_vjp

contains

    subroutine evaluate_helmholtz_representation_torus_curved_3d_geometry_jvp( &
            parameters, triangles, major_radius, minor_radius, &
            dirichlet_values, neumann_values, target, wave_number, &
            quadrature_degree, parameters_dot, major_radius_dot, &
            minor_radius_dot, target_dot, value_dot, status)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, target(3)
        real(dp), intent(in) :: wave_number
        complex(dp), intent(in) :: dirichlet_values(:), neumann_values(:)
        integer, intent(in) :: quadrature_degree
        real(dp), intent(in) :: parameters_dot(:, :)
        real(dp), intent(in) :: major_radius_dot, minor_radius_dot
        real(dp), intent(in) :: target_dot(3)
        complex(dp), intent(out) :: value_dot
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: barycentric(3), displacement(3), displacement_dot(3)
        real(dp) :: jacobian, jacobian_dot, normal(3), normal_dot(3)
        real(dp) :: point(3), point_dot(3), radius, radius_dot
        real(dp) :: tangent_eta(3), tangent_eta_dot(3)
        real(dp) :: tangent_xi(3), tangent_xi_dot(3)
        real(dp) :: area_cross(3), area_cross_dot(3), dot_normal
        complex(dp) :: green, green_dot, green_normal_derivative
        complex(dp) :: green_normal_derivative_dot, factor, factor_dot
        complex(dp) :: trace_value, neumann_value
        integer :: element, node, quadrature_status

        value_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_representation_inputs( &
            parameters, triangles, dirichlet_values, neumann_values, &
            wave_number, quadrature_degree)) return
        if (any(shape(parameters_dot) /= shape(parameters))) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        do element = 1, size(triangles, 2)
            do node = 1, size(weights)
                barycentric = [ &
                    1.0_dp - xi(node) - eta(node), xi(node), eta(node)]
                call evaluate_torus_curved_panel( &
                    parameters(:, triangles(:, element)), major_radius, &
                    minor_radius, xi(node), eta(node), point, tangent_xi, &
                    tangent_eta, jacobian, status)
                if (status /= 0) return
                call evaluate_torus_curved_panel_jvp( &
                    parameters(:, triangles(:, element)), major_radius, &
                    minor_radius, xi(node), eta(node), &
                    parameters_dot(:, triangles(:, element)), &
                    major_radius_dot, minor_radius_dot, 0.0_dp, 0.0_dp, &
                    point_dot, tangent_xi_dot, tangent_eta_dot, jacobian_dot, &
                    status)
                if (status /= 0) return
                area_cross = cross_product(tangent_xi, tangent_eta)
                area_cross_dot = cross_product(tangent_xi_dot, tangent_eta) + &
                    cross_product(tangent_xi, tangent_eta_dot)
                normal = area_cross/jacobian
                normal_dot = (area_cross_dot - normal*jacobian_dot)/jacobian
                displacement = target - point
                displacement_dot = target_dot - point_dot
                radius = norm2(displacement)
                if (radius <= tiny(1.0_dp)) return
                radius_dot = dot_product(displacement, displacement_dot)/radius
                green = exp(cmplx(0.0_dp, wave_number*radius, dp))/ &
                    (4.0_dp*acos(-1.0_dp)*radius)
                green_dot = green*cmplx( &
                    -1.0_dp/radius, wave_number, dp)*radius_dot
                dot_normal = dot_product(displacement, normal)
                factor = cmplx(1.0_dp, -wave_number*radius, dp)/radius**2
                factor_dot = (cmplx(0.0_dp, -wave_number, dp)/radius**2 - &
                    2.0_dp*cmplx(1.0_dp, -wave_number*radius, dp)/ &
                    radius**3)*radius_dot
                green_normal_derivative = green*factor*dot_normal
                green_normal_derivative_dot = &
                    (green_dot*factor + green*factor_dot)*dot_normal + &
                    green*factor*(dot_product(displacement_dot, normal) + &
                    dot_product(displacement, normal_dot))
                trace_value = dot_product( &
                    dirichlet_values(triangles(:, element)), barycentric)
                if (size(neumann_values) == size(parameters, 2)) then
                    neumann_value = dot_product( &
                        neumann_values(triangles(:, element)), barycentric)
                else
                    neumann_value = neumann_values(element)
                end if
                value_dot = value_dot + weights(node)*( &
                    jacobian_dot*(trace_value*green_normal_derivative - &
                    green*neumann_value) + jacobian*( &
                    trace_value*green_normal_derivative_dot - &
                    green_dot*neumann_value))
            end do
        end do
        status = 0
    end subroutine evaluate_helmholtz_representation_torus_curved_3d_geometry_jvp

    subroutine evaluate_helmholtz_representation_torus_curved_3d_geometry_vjp( &
            parameters, triangles, major_radius, minor_radius, &
            dirichlet_values, neumann_values, target, wave_number, &
            quadrature_degree, value_bar, parameters_bar, major_radius_bar, &
            minor_radius_bar, target_bar, status)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, target(3)
        real(dp), intent(in) :: wave_number
        complex(dp), intent(in) :: dirichlet_values(:), neumann_values(:)
        integer, intent(in) :: quadrature_degree
        complex(dp), intent(in) :: value_bar
        real(dp), intent(out) :: parameters_bar(:, :)
        real(dp), intent(out) :: major_radius_bar, minor_radius_bar
        real(dp), intent(out) :: target_bar(3)
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: barycentric(3), displacement(3)
        real(dp) :: jacobian, normal(3), point(3), radius
        real(dp) :: tangent_eta(3), tangent_xi(3), area_cross(3)
        real(dp) :: point_bar(3), tangent_eta_bar(3), tangent_xi_bar(3)
        real(dp) :: jacobian_bar, normal_bar(3), radius_bar, a_bar
        real(dp) :: displacement_bar(3), scale, dot_normal
        real(dp) :: local_major_radius_bar, local_minor_radius_bar
        real(dp) :: local_parameters_bar(2, 3), xi_bar, eta_bar
        complex(dp) :: green, factor, factor_derivative
        complex(dp) :: green_bar, factor_bar, trace_value, neumann_value
        complex(dp) :: green_normal_derivative, sum_bar
        integer :: element, local_node, node, quadrature_status

        parameters_bar = 0.0_dp
        major_radius_bar = 0.0_dp
        minor_radius_bar = 0.0_dp
        target_bar = 0.0_dp
        status = 1
        if (.not. valid_representation_inputs( &
            parameters, triangles, dirichlet_values, neumann_values, &
            wave_number, quadrature_degree)) return
        if (any(shape(parameters_bar) /= shape(parameters))) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        do element = 1, size(triangles, 2)
            do node = 1, size(weights)
                barycentric = [ &
                    1.0_dp - xi(node) - eta(node), xi(node), eta(node)]
                call evaluate_torus_curved_panel( &
                    parameters(:, triangles(:, element)), major_radius, &
                    minor_radius, xi(node), eta(node), point, tangent_xi, &
                    tangent_eta, jacobian, status)
                if (status /= 0) return
                normal = cross_product(tangent_xi, tangent_eta)/jacobian
                area_cross = cross_product(tangent_xi, tangent_eta)
                displacement = target - point
                radius = norm2(displacement)
                if (radius <= tiny(1.0_dp)) return
                green = exp(cmplx(0.0_dp, wave_number*radius, dp))/ &
                    (4.0_dp*acos(-1.0_dp)*radius)
                factor = cmplx(1.0_dp, -wave_number*radius, dp)/radius**2
                dot_normal = dot_product(displacement, normal)
                green_normal_derivative = green*factor*dot_normal
                trace_value = dot_product( &
                    dirichlet_values(triangles(:, element)), barycentric)
                if (size(neumann_values) == size(parameters, 2)) then
                    neumann_value = dot_product( &
                        neumann_values(triangles(:, element)), barycentric)
                else
                    neumann_value = neumann_values(element)
                end if
                scale = weights(node)
                sum_bar = scale*jacobian*value_bar
                jacobian_bar = scale*real( &
                    conjg(value_bar)*(trace_value*green_normal_derivative - &
                    green*neumann_value), dp)
                factor_bar = conjg(green)*dot_normal* &
                    (conjg(trace_value)*sum_bar)
                green_bar = -conjg(neumann_value)*sum_bar + &
                    conjg(factor)*dot_normal*(conjg(trace_value)*sum_bar)
                factor_derivative = cmplx(0.0_dp, -wave_number, dp)/ &
                    radius**2 - 2.0_dp*cmplx(1.0_dp, -wave_number*radius, dp) &
                    /radius**3
                radius_bar = real(conjg(green_bar)*green*cmplx( &
                    -1.0_dp/radius, wave_number, dp), dp) + &
                    real(conjg(factor_bar)*factor_derivative, dp)
                a_bar = real(conjg(conjg(trace_value)*sum_bar)* &
                    green*factor, dp)
                displacement_bar = radius_bar*displacement/radius + &
                    a_bar*normal
                normal_bar = a_bar*displacement
                point_bar = -displacement_bar
                target_bar = target_bar + displacement_bar
                tangent_xi_bar = cross_product(tangent_eta, normal_bar)/jacobian
                tangent_eta_bar = cross_product(normal_bar, tangent_xi)/jacobian
                jacobian_bar = jacobian_bar - dot_product(normal_bar, &
                    area_cross)/(jacobian**2)
                call evaluate_torus_curved_panel_vjp( &
                    parameters(:, triangles(:, element)), major_radius, &
                    minor_radius, xi(node), eta(node), point_bar, &
                    tangent_xi_bar, tangent_eta_bar, jacobian_bar, &
                    local_parameters_bar, local_major_radius_bar, &
                    local_minor_radius_bar, xi_bar, eta_bar, status)
                if (status /= 0) return
                do local_node = 1, 3
                    parameters_bar(:, triangles(local_node, element)) = &
                        parameters_bar(:, triangles(local_node, element)) + &
                        local_parameters_bar(:, local_node)
                end do
                major_radius_bar = major_radius_bar + local_major_radius_bar
                minor_radius_bar = minor_radius_bar + local_minor_radius_bar
            end do
        end do
        status = 0
    end subroutine evaluate_helmholtz_representation_torus_curved_3d_geometry_vjp

    logical function valid_representation_inputs( &
            parameters, triangles, dirichlet_values, neumann_values, &
            wave_number, quadrature_degree)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        complex(dp), intent(in) :: dirichlet_values(:), neumann_values(:)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: quadrature_degree

        valid_representation_inputs = .false.
        if (size(parameters, 1) /= 2 .or. size(parameters, 2) < 3) return
        if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
        if (any(triangles < 1) .or. &
            any(triangles > size(parameters, 2))) return
        if (size(dirichlet_values) /= size(parameters, 2)) return
        if (size(neumann_values) /= size(triangles, 2) .and. &
            size(neumann_values) /= size(parameters, 2)) return
        if (wave_number < 0.0_dp .or. quadrature_degree < 0) return
        valid_representation_inputs = .true.
    end function valid_representation_inputs

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_helmholtz_representation_ad_3d
