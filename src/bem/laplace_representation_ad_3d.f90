module fortfem_laplace_representation_ad_3d
    use fortfem_kinds, only: dp
    use fortfem_laplace_representation_3d, only: &
        evaluate_laplace_representation_torus_curved_3d
    use fortfem_torus_curved_panel, only: &
        evaluate_torus_curved_panel, evaluate_torus_curved_panel_jvp, &
        evaluate_torus_curved_panel_vjp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    private

    public :: evaluate_laplace_representation_torus_curved_3d_jvp
    public :: evaluate_laplace_representation_torus_curved_3d_vjp
    public :: evaluate_laplace_representation_torus_curved_3d_geometry_jvp
    public :: evaluate_laplace_representation_torus_curved_3d_geometry_vjp

contains

    subroutine evaluate_laplace_representation_torus_curved_3d_jvp( &
            parameters, triangles, major_radius, minor_radius, &
            dirichlet_values, neumann_values, target, quadrature_degree, &
            dirichlet_values_dot, neumann_values_dot, target_dot, value_dot, &
            status)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: dirichlet_values(:), neumann_values(:)
        real(dp), intent(in) :: target(3)
        integer, intent(in) :: quadrature_degree
        real(dp), intent(in) :: dirichlet_values_dot(:)
        real(dp), intent(in) :: neumann_values_dot(:), target_dot(3)
        real(dp), intent(out) :: value_dot
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: barycentric(3), displacement(3), green
        real(dp) :: green_normal_derivative, jacobian, normal(3)
        real(dp) :: point(3), tangent_eta(3), tangent_xi(3)
        real(dp) :: trace_dot, neumann_dot, gradient(3), base_value
        integer :: element, node, quadrature_status

        value_dot = 0.0_dp
        status = 1
        if (.not. valid_representation_inputs( &
            parameters, triangles, dirichlet_values, neumann_values, &
            quadrature_degree)) return
        if (size(dirichlet_values_dot) /= size(dirichlet_values)) return
        if (size(neumann_values_dot) /= size(neumann_values)) return
        call evaluate_laplace_representation_torus_curved_3d( &
            parameters, triangles, major_radius, minor_radius, &
            dirichlet_values, neumann_values, target, quadrature_degree, &
            base_value, status, gradient)
        if (status /= 0) return
        value_dot = dot_product(gradient, target_dot)
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) then
            status = 1
            return
        end if
        do element = 1, size(triangles, 2)
            do node = 1, size(weights)
                call evaluate_torus_curved_panel( &
                    parameters(:, triangles(:, element)), major_radius, &
                    minor_radius, xi(node), eta(node), point, tangent_xi, &
                    tangent_eta, jacobian, status)
                if (status /= 0) return
                normal = cross_product(tangent_xi, tangent_eta)/jacobian
                barycentric = [ &
                    1.0_dp - xi(node) - eta(node), xi(node), eta(node)]
                displacement = target - point
                if (norm2(displacement) <= 0.0_dp) return
                green = 1.0_dp/(4.0_dp*acos(-1.0_dp)*norm2(displacement))
                green_normal_derivative = dot_product(displacement, normal)/ &
                    (4.0_dp*acos(-1.0_dp)*norm2(displacement)**3)
                trace_dot = dot_product( &
                    dirichlet_values_dot(triangles(:, element)), barycentric)
                if (size(neumann_values_dot) == size(parameters, 2)) then
                    neumann_dot = dot_product( &
                        neumann_values_dot(triangles(:, element)), barycentric)
                else
                    neumann_dot = neumann_values_dot(element)
                end if
                value_dot = value_dot + jacobian*weights(node)*( &
                    trace_dot*green_normal_derivative - green*neumann_dot)
            end do
        end do
        status = 0
    end subroutine evaluate_laplace_representation_torus_curved_3d_jvp

    subroutine evaluate_laplace_representation_torus_curved_3d_vjp( &
            parameters, triangles, major_radius, minor_radius, &
            dirichlet_values, neumann_values, target, quadrature_degree, &
            value_bar, dirichlet_values_bar, neumann_values_bar, target_bar, &
            status)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: dirichlet_values(:), neumann_values(:)
        real(dp), intent(in) :: target(3)
        integer, intent(in) :: quadrature_degree
        real(dp), intent(in) :: value_bar
        real(dp), intent(out) :: dirichlet_values_bar(:)
        real(dp), intent(out) :: neumann_values_bar(:), target_bar(3)
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: barycentric(3), displacement(3), green
        real(dp) :: green_normal_derivative, jacobian, normal(3)
        real(dp) :: normal_weight, point(3), tangent_eta(3), tangent_xi(3)
        real(dp) :: gradient(3), base_value
        integer :: element, local_node, node, quadrature_status

        dirichlet_values_bar = 0.0_dp
        neumann_values_bar = 0.0_dp
        target_bar = 0.0_dp
        status = 1
        if (.not. valid_representation_inputs( &
            parameters, triangles, dirichlet_values, neumann_values, &
            quadrature_degree)) return
        if (size(dirichlet_values_bar) /= size(dirichlet_values)) return
        if (size(neumann_values_bar) /= size(neumann_values)) return
        call evaluate_laplace_representation_torus_curved_3d( &
            parameters, triangles, major_radius, minor_radius, &
            dirichlet_values, neumann_values, target, quadrature_degree, &
            base_value, status, gradient)
        if (status /= 0) return
        target_bar = value_bar*gradient
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) then
            status = 1
            return
        end if
        do element = 1, size(triangles, 2)
            do node = 1, size(weights)
                call evaluate_torus_curved_panel( &
                    parameters(:, triangles(:, element)), major_radius, &
                    minor_radius, xi(node), eta(node), point, tangent_xi, &
                    tangent_eta, jacobian, status)
                if (status /= 0) return
                normal = cross_product(tangent_xi, tangent_eta)/jacobian
                barycentric = [ &
                    1.0_dp - xi(node) - eta(node), xi(node), eta(node)]
                displacement = target - point
                if (norm2(displacement) <= 0.0_dp) return
                green = 1.0_dp/(4.0_dp*acos(-1.0_dp)*norm2(displacement))
                green_normal_derivative = dot_product(displacement, normal)/ &
                    (4.0_dp*acos(-1.0_dp)*norm2(displacement)**3)
                normal_weight = value_bar*jacobian*weights(node)
                do local_node = 1, 3
                    dirichlet_values_bar(triangles(local_node, element)) = &
                        dirichlet_values_bar(triangles(local_node, element)) + &
                        normal_weight*barycentric(local_node)* &
                        green_normal_derivative
                    if (size(neumann_values) == size(parameters, 2)) then
                        neumann_values_bar(triangles(local_node, element)) = &
                            neumann_values_bar(triangles(local_node, element)) - &
                            normal_weight*barycentric(local_node)*green
                    end if
                end do
                if (size(neumann_values) == size(triangles, 2)) then
                    neumann_values_bar(element) = neumann_values_bar(element) - &
                        normal_weight*green
                end if
            end do
        end do
        status = 0
    end subroutine evaluate_laplace_representation_torus_curved_3d_vjp

    subroutine evaluate_laplace_representation_torus_curved_3d_geometry_jvp( &
            parameters, triangles, major_radius, minor_radius, &
            dirichlet_values, neumann_values, target, quadrature_degree, &
            parameters_dot, major_radius_dot, minor_radius_dot, target_dot, &
            value_dot, status)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: dirichlet_values(:), neumann_values(:)
        real(dp), intent(in) :: target(3)
        integer, intent(in) :: quadrature_degree
        real(dp), intent(in) :: parameters_dot(:, :)
        real(dp), intent(in) :: major_radius_dot, minor_radius_dot
        real(dp), intent(in) :: target_dot(3)
        real(dp), intent(out) :: value_dot
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: barycentric(3), displacement(3), displacement_dot(3)
        real(dp) :: green, green_dot, green_normal_derivative
        real(dp) :: green_normal_derivative_dot, jacobian, jacobian_dot
        real(dp) :: normal(3), normal_dot(3), point(3), point_dot(3)
        real(dp) :: tangent_eta(3), tangent_eta_dot(3)
        real(dp) :: tangent_xi(3), tangent_xi_dot(3)
        real(dp) :: trace_value, neumann_value, radius, radius_dot
        real(dp) :: area_cross(3), area_cross_dot(3)
        integer :: element, node, quadrature_status

        value_dot = 0.0_dp
        status = 1
        if (.not. valid_representation_inputs( &
            parameters, triangles, dirichlet_values, neumann_values, &
            quadrature_degree)) return
        if (any(shape(parameters_dot) /= shape(parameters))) return
        if (size(target_dot) /= 3) return
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
                if (radius <= 0.0_dp) return
                radius_dot = dot_product(displacement, displacement_dot)/radius
                green = 1.0_dp/(4.0_dp*acos(-1.0_dp)*radius)
                green_dot = -green*radius_dot/radius
                green_normal_derivative = &
                    dot_product(displacement, normal)/ &
                    (4.0_dp*acos(-1.0_dp)*radius**3)
                green_normal_derivative_dot = &
                    (dot_product(displacement_dot, normal) + &
                    dot_product(displacement, normal_dot) - &
                    3.0_dp*dot_product(displacement, normal)*radius_dot/radius) &
                    /(4.0_dp*acos(-1.0_dp)*radius**3)
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
    end subroutine evaluate_laplace_representation_torus_curved_3d_geometry_jvp

    subroutine evaluate_laplace_representation_torus_curved_3d_geometry_vjp( &
            parameters, triangles, major_radius, minor_radius, &
            dirichlet_values, neumann_values, target, quadrature_degree, &
            value_bar, parameters_bar, major_radius_bar, minor_radius_bar, &
            target_bar, status)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: dirichlet_values(:), neumann_values(:)
        real(dp), intent(in) :: target(3)
        integer, intent(in) :: quadrature_degree
        real(dp), intent(in) :: value_bar
        real(dp), intent(out) :: parameters_bar(:, :)
        real(dp), intent(out) :: major_radius_bar, minor_radius_bar
        real(dp), intent(out) :: target_bar(3)
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: barycentric(3), displacement(3), green
        real(dp) :: green_normal_derivative, jacobian, normal(3)
        real(dp) :: point(3), tangent_eta(3), tangent_xi(3)
        real(dp) :: point_bar(3), tangent_eta_bar(3), tangent_xi_bar(3)
        real(dp) :: jacobian_bar, normal_bar(3), area_cross(3)
        real(dp) :: trace_value, neumann_value, radius, radius_bar
        real(dp) :: displacement_bar(3), green_bar, normal_derivative_bar
        real(dp) :: scale, local_major_radius_bar, local_minor_radius_bar
        real(dp) :: local_parameters_bar(2, 3), xi_bar, eta_bar
        integer :: element, local_node, node, quadrature_status

        parameters_bar = 0.0_dp
        major_radius_bar = 0.0_dp
        minor_radius_bar = 0.0_dp
        target_bar = 0.0_dp
        status = 1
        if (.not. valid_representation_inputs( &
            parameters, triangles, dirichlet_values, neumann_values, &
            quadrature_degree)) return
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
                if (radius <= 0.0_dp) return
                green = 1.0_dp/(4.0_dp*acos(-1.0_dp)*radius)
                green_normal_derivative = &
                    dot_product(displacement, normal)/ &
                    (4.0_dp*acos(-1.0_dp)*radius**3)
                trace_value = dot_product( &
                    dirichlet_values(triangles(:, element)), barycentric)
                if (size(neumann_values) == size(parameters, 2)) then
                    neumann_value = dot_product( &
                        neumann_values(triangles(:, element)), barycentric)
                else
                    neumann_value = neumann_values(element)
                end if
                scale = value_bar*weights(node)
                jacobian_bar = scale*(trace_value*green_normal_derivative - &
                    green*neumann_value)
                normal_derivative_bar = scale*jacobian*trace_value
                green_bar = -scale*jacobian*neumann_value
                normal_bar = normal_derivative_bar*displacement/( &
                    4.0_dp*acos(-1.0_dp)*radius**3)
                radius_bar = normal_derivative_bar*dot_product( &
                    displacement, normal)*(-3.0_dp)/( &
                    4.0_dp*acos(-1.0_dp)*radius**4)
                radius_bar = radius_bar - green_bar*green/radius
                displacement_bar = radius_bar*displacement/radius + &
                    normal_derivative_bar*normal/( &
                    4.0_dp*acos(-1.0_dp)*radius**3)
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
    end subroutine evaluate_laplace_representation_torus_curved_3d_geometry_vjp

    logical function valid_representation_inputs( &
            parameters, triangles, dirichlet_values, neumann_values, &
            quadrature_degree)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: dirichlet_values(:), neumann_values(:)
        integer, intent(in) :: quadrature_degree

        valid_representation_inputs = .false.
        if (size(parameters, 1) /= 2 .or. size(parameters, 2) < 3) return
        if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
        if (any(triangles < 1) .or. &
            any(triangles > size(parameters, 2))) return
        if (size(dirichlet_values) /= size(parameters, 2)) return
        if (size(neumann_values) /= size(triangles, 2) .and. &
            size(neumann_values) /= size(parameters, 2)) return
        if (quadrature_degree < 0) return
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

end module fortfem_laplace_representation_ad_3d
