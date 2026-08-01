module fortfem_laplace_torus_curved_bem_ad_3d
    use fortfem_kinds, only: dp
    use fortfem_torus_curved_panel, only: &
        evaluate_torus_curved_panel, evaluate_torus_curved_panel_jvp, &
        evaluate_torus_curved_panel_vjp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: reference_vertices(2, 3) = reshape( &
        [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])

    public :: assemble_laplace_torus_curved_dtn_3d_geometry_jvp
    public :: assemble_laplace_torus_curved_dtn_3d_geometry_vjp

contains

    subroutine assemble_laplace_torus_curved_dtn_3d_geometry_jvp( &
            parameters, triangles, major_radius, minor_radius, &
            quadrature_degree, parameters_dot, major_radius_dot, &
            minor_radius_dot, single_layer_dot, double_layer_dot, mass_dot, &
            status)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: quadrature_degree
        real(dp), intent(in) :: parameters_dot(:, :)
        real(dp), intent(in) :: major_radius_dot, minor_radius_dot
        real(dp), allocatable, intent(out) :: single_layer_dot(:, :)
        real(dp), allocatable, intent(out) :: double_layer_dot(:, :)
        real(dp), allocatable, intent(out) :: mass_dot(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), line_nodes(:), line_weights(:)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp) :: local_double_dot(3), local_mass_dot(3)
        real(dp) :: local_single_dot
        integer :: first_panel, line_count, quadrature_status, second_panel

        status = 1
        if (.not. valid_dtn_inputs( &
            parameters, triangles, major_radius, minor_radius, &
            quadrature_degree)) return
        if (any(shape(parameters_dot) /= shape(parameters))) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        line_count = max(2, (quadrature_degree + 3)/2)
        allocate(line_nodes(line_count), line_weights(line_count))
        call gauss_legendre_ab( &
            line_count, 0.0_dp, 1.0_dp, line_nodes, line_weights)
        allocate( &
            single_layer_dot(size(triangles, 2), size(triangles, 2)), &
            double_layer_dot(size(triangles, 2), size(parameters, 2)), &
            mass_dot(size(triangles, 2), size(parameters, 2)))
        single_layer_dot = 0.0_dp
        double_layer_dot = 0.0_dp
        mass_dot = 0.0_dp
        do first_panel = 1, size(triangles, 2)
            call integrate_mass_jvp( &
                parameters(:, triangles(:, first_panel)), &
                parameters_dot(:, triangles(:, first_panel)), major_radius, &
                minor_radius, xi, eta, weights, major_radius_dot, &
                minor_radius_dot, local_mass_dot, status)
            if (status /= 0) return
            mass_dot(first_panel, triangles(:, first_panel)) = local_mass_dot
            do second_panel = 1, size(triangles, 2)
                if (first_panel == second_panel) then
                    call integrate_coincident_pair_jvp( &
                        parameters(:, triangles(:, first_panel)), &
                        parameters_dot(:, triangles(:, first_panel)), &
                        major_radius, minor_radius, xi, eta, weights, &
                        line_nodes, line_weights, major_radius_dot, &
                        minor_radius_dot, local_single_dot, local_double_dot, &
                        status)
                else
                    call integrate_regular_pair_jvp( &
                        parameters(:, triangles(:, first_panel)), &
                        parameters(:, triangles(:, second_panel)), &
                        parameters_dot(:, triangles(:, first_panel)), &
                        parameters_dot(:, triangles(:, second_panel)), &
                        major_radius, minor_radius, xi, eta, weights, &
                        major_radius_dot, minor_radius_dot, local_single_dot, &
                        local_double_dot, status)
                end if
                if (status /= 0) return
                single_layer_dot(first_panel, second_panel) = local_single_dot
                double_layer_dot(first_panel, triangles(:, second_panel)) = &
                    double_layer_dot(first_panel, triangles(:, second_panel)) + &
                    local_double_dot
            end do
        end do
        status = 0
    end subroutine assemble_laplace_torus_curved_dtn_3d_geometry_jvp

    subroutine assemble_laplace_torus_curved_dtn_3d_geometry_vjp( &
            parameters, triangles, major_radius, minor_radius, &
            quadrature_degree, single_layer_bar, double_layer_bar, mass_bar, &
            parameters_bar, major_radius_bar, minor_radius_bar, status)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: quadrature_degree
        real(dp), intent(in) :: single_layer_bar(:, :)
        real(dp), intent(in) :: double_layer_bar(:, :), mass_bar(:, :)
        real(dp), intent(out) :: parameters_bar(:, :)
        real(dp), intent(out) :: major_radius_bar, minor_radius_bar
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), line_nodes(:), line_weights(:)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp) :: local_double_bar(3), local_mass_bar(3)
        real(dp) :: local_parameters_bar(2, 3), local_major_radius_bar
        real(dp) :: local_minor_radius_bar, local_single_bar
        real(dp) :: local_parameters_bar_second(2, 3)
        real(dp) :: local_major_radius_bar_second, local_minor_radius_bar_second
        integer :: first_panel, line_count, local_node
        integer :: quadrature_status, second_panel

        parameters_bar = 0.0_dp
        major_radius_bar = 0.0_dp
        minor_radius_bar = 0.0_dp
        status = 1
        if (.not. valid_dtn_inputs( &
            parameters, triangles, major_radius, minor_radius, &
            quadrature_degree)) return
        if (any(shape(parameters_bar) /= shape(parameters))) return
        if (size(single_layer_bar, 1) /= size(triangles, 2) .or. &
            size(single_layer_bar, 2) /= size(triangles, 2)) return
        if (size(double_layer_bar, 1) /= size(triangles, 2) .or. &
            size(double_layer_bar, 2) /= size(parameters, 2)) return
        if (size(mass_bar, 1) /= size(triangles, 2) .or. &
            size(mass_bar, 2) /= size(parameters, 2)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        line_count = max(2, (quadrature_degree + 3)/2)
        allocate(line_nodes(line_count), line_weights(line_count))
        call gauss_legendre_ab( &
            line_count, 0.0_dp, 1.0_dp, line_nodes, line_weights)
        do first_panel = 1, size(triangles, 2)
            local_mass_bar = mass_bar(first_panel, triangles(:, first_panel))
            call integrate_mass_vjp( &
                parameters(:, triangles(:, first_panel)), major_radius, &
                minor_radius, xi, eta, weights, local_mass_bar, &
                local_parameters_bar, local_major_radius_bar, &
                local_minor_radius_bar, status)
            if (status /= 0) return
            do local_node = 1, 3
                parameters_bar(:, triangles(local_node, first_panel)) = &
                    parameters_bar(:, triangles(local_node, first_panel)) + &
                    local_parameters_bar(:, local_node)
            end do
            major_radius_bar = major_radius_bar + local_major_radius_bar
            minor_radius_bar = minor_radius_bar + local_minor_radius_bar
            do second_panel = 1, size(triangles, 2)
                local_single_bar = single_layer_bar(first_panel, second_panel)
                local_double_bar = double_layer_bar( &
                    first_panel, triangles(:, second_panel))
                if (first_panel == second_panel) then
                    call integrate_coincident_pair_vjp( &
                        parameters(:, triangles(:, first_panel)), major_radius, &
                        minor_radius, xi, eta, weights, line_nodes, &
                        line_weights, local_single_bar, local_double_bar, &
                        local_parameters_bar, local_major_radius_bar, &
                        local_minor_radius_bar, status)
                else
                    call integrate_regular_pair_vjp( &
                        parameters(:, triangles(:, first_panel)), &
                        parameters(:, triangles(:, second_panel)), major_radius, &
                        minor_radius, xi, eta, weights, local_single_bar, &
                        local_double_bar, local_parameters_bar, &
                        local_major_radius_bar, local_minor_radius_bar, &
                        local_parameters_bar_second, local_major_radius_bar_second, &
                        local_minor_radius_bar_second, status)
                end if
                if (status /= 0) return
                do local_node = 1, 3
                    parameters_bar(:, triangles(local_node, first_panel)) = &
                        parameters_bar(:, triangles(local_node, first_panel)) + &
                        local_parameters_bar(:, local_node)
                    if (first_panel /= second_panel) then
                        parameters_bar(:, triangles(local_node, second_panel)) = &
                            parameters_bar(:, triangles(local_node, second_panel)) + &
                            local_parameters_bar_second(:, local_node)
                    end if
                end do
                major_radius_bar = major_radius_bar + local_major_radius_bar
                minor_radius_bar = minor_radius_bar + local_minor_radius_bar
                if (first_panel /= second_panel) then
                    major_radius_bar = major_radius_bar + &
                        local_major_radius_bar_second
                    minor_radius_bar = minor_radius_bar + &
                        local_minor_radius_bar_second
                end if
            end do
        end do
        status = 0
    end subroutine assemble_laplace_torus_curved_dtn_3d_geometry_vjp

    subroutine integrate_mass_jvp( &
            panel_parameters, panel_parameters_dot, major_radius, minor_radius, &
            xi, eta, weights, major_radius_dot, minor_radius_dot, mass_dot, &
            status)
        real(dp), intent(in) :: panel_parameters(2, 3), panel_parameters_dot(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        real(dp), intent(in) :: major_radius_dot, minor_radius_dot
        real(dp), intent(out) :: mass_dot(3)
        integer, intent(out) :: status

        real(dp) :: point(3), point_dot(3), tangent_eta(3), tangent_eta_dot(3)
        real(dp) :: tangent_xi(3), tangent_xi_dot(3), jacobian
        real(dp) :: jacobian_dot, barycentric(3)
        integer :: node

        mass_dot = 0.0_dp
        status = 1
        do node = 1, size(weights)
            call evaluate_torus_curved_panel( &
                panel_parameters, major_radius, minor_radius, xi(node), &
                eta(node), point, tangent_xi, tangent_eta, jacobian, status)
            if (status /= 0) return
            call evaluate_torus_curved_panel_jvp( &
                panel_parameters, major_radius, minor_radius, xi(node), &
                eta(node), panel_parameters_dot, major_radius_dot, &
                minor_radius_dot, 0.0_dp, 0.0_dp, point_dot, tangent_xi_dot, &
                tangent_eta_dot, jacobian_dot, status)
            if (status /= 0) return
            barycentric = [1.0_dp - xi(node) - eta(node), xi(node), eta(node)]
            mass_dot = mass_dot + weights(node)*jacobian_dot*barycentric
        end do
        status = 0
    end subroutine integrate_mass_jvp

    subroutine integrate_mass_vjp( &
            panel_parameters, major_radius, minor_radius, xi, eta, weights, &
            mass_bar, parameters_bar, major_radius_bar, minor_radius_bar, &
            status)
        real(dp), intent(in) :: panel_parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: xi(:), eta(:), weights(:), mass_bar(3)
        real(dp), intent(out) :: parameters_bar(2, 3)
        real(dp), intent(out) :: major_radius_bar, minor_radius_bar
        integer, intent(out) :: status

        real(dp) :: point(3), tangent_eta(3), tangent_xi(3), jacobian
        real(dp) :: jacobian_bar, barycentric(3), xi_bar, eta_bar
        real(dp) :: local_parameters_bar(2, 3), local_major_radius_bar
        real(dp) :: local_minor_radius_bar
        integer :: node

        parameters_bar = 0.0_dp
        major_radius_bar = 0.0_dp
        minor_radius_bar = 0.0_dp
        status = 1
        do node = 1, size(weights)
            call evaluate_torus_curved_panel( &
                panel_parameters, major_radius, minor_radius, xi(node), &
                eta(node), point, tangent_xi, tangent_eta, jacobian, status)
            if (status /= 0) return
            barycentric = [1.0_dp - xi(node) - eta(node), xi(node), eta(node)]
            jacobian_bar = weights(node)*dot_product(mass_bar, barycentric)
            call evaluate_torus_curved_panel_vjp( &
                panel_parameters, major_radius, minor_radius, xi(node), &
                eta(node), [0.0_dp, 0.0_dp, 0.0_dp], &
                [0.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 0.0_dp, 0.0_dp], &
                jacobian_bar, local_parameters_bar, local_major_radius_bar, &
                local_minor_radius_bar, xi_bar, eta_bar, status)
            if (status /= 0) return
            parameters_bar = parameters_bar + local_parameters_bar
            major_radius_bar = major_radius_bar + local_major_radius_bar
            minor_radius_bar = minor_radius_bar + local_minor_radius_bar
        end do
        status = 0
    end subroutine integrate_mass_vjp

    subroutine integrate_regular_pair_jvp( &
            target_parameters, source_parameters, target_parameters_dot, &
            source_parameters_dot, major_radius, minor_radius, xi, eta, &
            weights, major_radius_dot, minor_radius_dot, single_dot, &
            double_dot, status)
        real(dp), intent(in) :: target_parameters(2, 3), source_parameters(2, 3)
        real(dp), intent(in) :: target_parameters_dot(2, 3)
        real(dp), intent(in) :: source_parameters_dot(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        real(dp), intent(in) :: major_radius_dot, minor_radius_dot
        real(dp), intent(out) :: single_dot, double_dot(3)
        integer, intent(out) :: status

        real(dp) :: target_point(3), target_point_dot(3)
        real(dp) :: source_point(3), source_point_dot(3)
        real(dp) :: target_tangent_xi(3), target_tangent_eta(3)
        real(dp) :: source_tangent_xi(3), source_tangent_eta(3)
        real(dp) :: target_tangent_xi_dot(3), target_tangent_eta_dot(3)
        real(dp) :: source_tangent_xi_dot(3), source_tangent_eta_dot(3)
        real(dp) :: target_jacobian, source_jacobian
        real(dp) :: target_jacobian_dot, source_jacobian_dot
        real(dp) :: source_barycentric(3)
        integer :: source_node, target_node

        single_dot = 0.0_dp
        double_dot = 0.0_dp
        status = 1
        do target_node = 1, size(weights)
            call evaluate_torus_curved_panel( &
                target_parameters, major_radius, minor_radius, &
                xi(target_node), eta(target_node), target_point, &
                target_tangent_xi, target_tangent_eta, target_jacobian, status)
            if (status /= 0) return
            call evaluate_torus_curved_panel_jvp( &
                target_parameters, major_radius, minor_radius, &
                xi(target_node), eta(target_node), target_parameters_dot, &
                major_radius_dot, minor_radius_dot, 0.0_dp, 0.0_dp, &
                target_point_dot, target_tangent_xi_dot, target_tangent_eta_dot, &
                target_jacobian_dot, status)
            if (status /= 0) return
            do source_node = 1, size(weights)
                call evaluate_torus_curved_panel( &
                    source_parameters, major_radius, minor_radius, &
                    xi(source_node), eta(source_node), source_point, &
                    source_tangent_xi, source_tangent_eta, source_jacobian, &
                    status)
                if (status /= 0) return
                call evaluate_torus_curved_panel_jvp( &
                    source_parameters, major_radius, minor_radius, &
                    xi(source_node), eta(source_node), source_parameters_dot, &
                    major_radius_dot, minor_radius_dot, 0.0_dp, 0.0_dp, &
                    source_point_dot, source_tangent_xi_dot, &
                    source_tangent_eta_dot, source_jacobian_dot, status)
                if (status /= 0) return
                source_barycentric = [ &
                    1.0_dp - xi(source_node) - eta(source_node), &
                    xi(source_node), eta(source_node)]
                call pair_jvp_term( &
                    target_point, target_point_dot, target_tangent_xi, &
                    target_tangent_eta, target_jacobian, target_tangent_xi_dot, &
                    target_tangent_eta_dot, target_jacobian_dot, source_point, &
                    source_point_dot, source_tangent_xi, source_tangent_eta, &
                    source_jacobian, source_tangent_xi_dot, &
                    source_tangent_eta_dot, source_jacobian_dot, &
                    weights(target_node)*weights(source_node)/(4.0_dp*pi), &
                    source_barycentric, single_dot, double_dot, status)
                if (status /= 0) return
            end do
        end do
        status = 0
    end subroutine integrate_regular_pair_jvp

    subroutine integrate_coincident_pair_jvp( &
            panel_parameters, panel_parameters_dot, major_radius, minor_radius, &
            xi, eta, weights, line_nodes, line_weights, major_radius_dot, &
            minor_radius_dot, single_dot, double_dot, status)
        real(dp), intent(in) :: panel_parameters(2, 3), panel_parameters_dot(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: xi(:), eta(:), weights(:), line_nodes(:)
        real(dp), intent(in) :: line_weights(:), major_radius_dot, minor_radius_dot
        real(dp), intent(out) :: single_dot, double_dot(3)
        integer, intent(out) :: status

        real(dp) :: direction(2), first_reference(2), second_reference(2)
        real(dp) :: wedge_first(2), wedge_second(2), wedge_jacobian, rho, t
        real(dp) :: target_point(3), target_point_dot(3), source_point(3)
        real(dp) :: source_point_dot(3), target_tangent_xi(3)
        real(dp) :: target_tangent_eta(3), source_tangent_xi(3)
        real(dp) :: source_tangent_eta(3), target_tangent_xi_dot(3)
        real(dp) :: target_tangent_eta_dot(3), source_tangent_xi_dot(3)
        real(dp) :: source_tangent_eta_dot(3), target_jacobian
        real(dp) :: source_jacobian, target_jacobian_dot, source_jacobian_dot
        real(dp) :: source_barycentric(3), reference_weight
        integer :: radial_node, tangential_node, target_node, wedge

        single_dot = 0.0_dp
        double_dot = 0.0_dp
        status = 1
        do target_node = 1, size(weights)
            first_reference = [xi(target_node), eta(target_node)]
            call evaluate_torus_curved_panel( &
                panel_parameters, major_radius, minor_radius, &
                first_reference(1), first_reference(2), target_point, &
                target_tangent_xi, target_tangent_eta, target_jacobian, status)
            if (status /= 0) return
            call evaluate_torus_curved_panel_jvp( &
                panel_parameters, major_radius, minor_radius, &
                first_reference(1), first_reference(2), panel_parameters_dot, &
                major_radius_dot, minor_radius_dot, 0.0_dp, 0.0_dp, &
                target_point_dot, target_tangent_xi_dot, &
                target_tangent_eta_dot, target_jacobian_dot, status)
            if (status /= 0) return
            do wedge = 1, 3
                wedge_first = reference_vertices(:, wedge) - first_reference
                wedge_second = reference_vertices(:, modulo(wedge, 3) + 1) - &
                    first_reference
                wedge_jacobian = abs( &
                    wedge_first(1)*wedge_second(2) - &
                    wedge_first(2)*wedge_second(1))
                do radial_node = 1, size(line_nodes)
                    rho = line_nodes(radial_node)
                    do tangential_node = 1, size(line_nodes)
                        t = line_nodes(tangential_node)
                        direction = (1.0_dp - t)*wedge_first + t*wedge_second
                        second_reference = first_reference + rho*direction
                        call evaluate_torus_curved_panel( &
                            panel_parameters, major_radius, minor_radius, &
                            second_reference(1), second_reference(2), &
                            source_point, source_tangent_xi, &
                            source_tangent_eta, source_jacobian, status)
                        if (status /= 0) return
                        call evaluate_torus_curved_panel_jvp( &
                            panel_parameters, major_radius, minor_radius, &
                            second_reference(1), second_reference(2), &
                            panel_parameters_dot, major_radius_dot, &
                            minor_radius_dot, 0.0_dp, 0.0_dp, source_point_dot, &
                            source_tangent_xi_dot, source_tangent_eta_dot, &
                            source_jacobian_dot, status)
                        if (status /= 0) return
                        source_barycentric = [ &
                            1.0_dp - sum(second_reference), second_reference(1), &
                            second_reference(2)]
                        reference_weight = weights(target_node)* &
                            line_weights(radial_node)*line_weights(tangential_node)* &
                            rho*wedge_jacobian/(4.0_dp*pi)
                        call pair_jvp_term( &
                            target_point, target_point_dot, target_tangent_xi, &
                            target_tangent_eta, target_jacobian, &
                            target_tangent_xi_dot, target_tangent_eta_dot, &
                            target_jacobian_dot, source_point, source_point_dot, &
                            source_tangent_xi, source_tangent_eta, source_jacobian, &
                            source_tangent_xi_dot, source_tangent_eta_dot, &
                            source_jacobian_dot, reference_weight, &
                            source_barycentric, single_dot, double_dot, status)
                        if (status /= 0) return
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine integrate_coincident_pair_jvp

    subroutine pair_jvp_term( &
            target_point, target_point_dot, target_tangent_xi, &
            target_tangent_eta, target_jacobian, target_tangent_xi_dot, &
            target_tangent_eta_dot, target_jacobian_dot, source_point, &
            source_point_dot, source_tangent_xi, source_tangent_eta, &
            source_jacobian, source_tangent_xi_dot, source_tangent_eta_dot, &
            source_jacobian_dot, reference_weight, source_barycentric, &
            single_dot, double_dot, status)
        real(dp), intent(in) :: target_point(3), target_point_dot(3)
        real(dp), intent(in) :: target_tangent_xi(3), target_tangent_eta(3)
        real(dp), intent(in) :: target_jacobian, target_tangent_xi_dot(3)
        real(dp), intent(in) :: target_tangent_eta_dot(3), target_jacobian_dot
        real(dp), intent(in) :: source_point(3), source_point_dot(3)
        real(dp), intent(in) :: source_tangent_xi(3), source_tangent_eta(3)
        real(dp), intent(in) :: source_jacobian, source_tangent_xi_dot(3)
        real(dp), intent(in) :: source_tangent_eta_dot(3), source_jacobian_dot
        real(dp), intent(in) :: reference_weight, source_barycentric(3)
        real(dp), intent(inout) :: single_dot, double_dot(3)
        integer, intent(out) :: status

        real(dp) :: displacement(3), displacement_dot(3), radius, radius_dot
        real(dp) :: source_normal(3), source_normal_dot(3), area_cross(3)
        real(dp) :: area_cross_dot(3), dot_normal, dot_normal_dot
        real(dp) :: weight, weight_dot, kernel, kernel_dot

        status = 1
        displacement = target_point - source_point
        displacement_dot = target_point_dot - source_point_dot
        radius = norm2(displacement)
        if (radius <= tiny(1.0_dp)) return
        radius_dot = dot_product(displacement, displacement_dot)/radius
        area_cross = cross_product(source_tangent_xi, source_tangent_eta)
        area_cross_dot = cross_product( &
            source_tangent_xi_dot, source_tangent_eta) + &
            cross_product(source_tangent_xi, source_tangent_eta_dot)
        source_normal = area_cross/source_jacobian
        source_normal_dot = (area_cross_dot - &
            source_normal*source_jacobian_dot)/source_jacobian
        dot_normal = dot_product(displacement, source_normal)
        dot_normal_dot = dot_product(displacement_dot, source_normal) + &
            dot_product(displacement, source_normal_dot)
        weight = reference_weight*target_jacobian*source_jacobian
        weight_dot = reference_weight*(target_jacobian_dot*source_jacobian + &
            target_jacobian*source_jacobian_dot)
        kernel = dot_normal/radius**3
        kernel_dot = dot_normal_dot/radius**3 - &
            3.0_dp*dot_normal*radius_dot/radius**4
        single_dot = single_dot + weight_dot/radius - &
            weight*radius_dot/radius**2
        double_dot = double_dot + (weight_dot*kernel + weight*kernel_dot)* &
            source_barycentric
        status = 0
    end subroutine pair_jvp_term

    subroutine integrate_regular_pair_vjp( &
            target_parameters, source_parameters, major_radius, minor_radius, &
            xi, eta, weights, single_bar, double_bar, target_parameters_bar, &
            target_major_radius_bar, target_minor_radius_bar, &
            source_parameters_bar, source_major_radius_bar, &
            source_minor_radius_bar, status)
        real(dp), intent(in) :: target_parameters(2, 3), source_parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        real(dp), intent(in) :: single_bar, double_bar(3)
        real(dp), intent(out) :: target_parameters_bar(2, 3)
        real(dp), intent(out) :: target_major_radius_bar, target_minor_radius_bar
        real(dp), intent(out) :: source_parameters_bar(2, 3)
        real(dp), intent(out) :: source_major_radius_bar, source_minor_radius_bar
        integer, intent(out) :: status

        real(dp) :: target_point(3), source_point(3), target_tangent_xi(3)
        real(dp) :: target_tangent_eta(3), source_tangent_xi(3)
        real(dp) :: source_tangent_eta(3), target_jacobian, source_jacobian
        real(dp) :: source_barycentric(3), reference_weight
        real(dp) :: target_point_bar(3), source_point_bar(3)
        real(dp) :: target_tangent_xi_bar(3), target_tangent_eta_bar(3)
        real(dp) :: source_tangent_xi_bar(3), source_tangent_eta_bar(3)
        real(dp) :: target_jacobian_bar, source_jacobian_bar
        real(dp) :: local_parameters_bar(2, 3), local_major_radius_bar
        real(dp) :: local_minor_radius_bar, xi_bar, eta_bar
        integer :: source_node, target_node

        target_parameters_bar = 0.0_dp
        source_parameters_bar = 0.0_dp
        target_major_radius_bar = 0.0_dp
        target_minor_radius_bar = 0.0_dp
        source_major_radius_bar = 0.0_dp
        source_minor_radius_bar = 0.0_dp
        status = 1
        do target_node = 1, size(weights)
            call evaluate_torus_curved_panel( &
                target_parameters, major_radius, minor_radius, &
                xi(target_node), eta(target_node), target_point, &
                target_tangent_xi, target_tangent_eta, target_jacobian, status)
            if (status /= 0) return
            do source_node = 1, size(weights)
                call evaluate_torus_curved_panel( &
                    source_parameters, major_radius, minor_radius, &
                    xi(source_node), eta(source_node), source_point, &
                    source_tangent_xi, source_tangent_eta, source_jacobian, &
                    status)
                if (status /= 0) return
                source_barycentric = [ &
                    1.0_dp - xi(source_node) - eta(source_node), &
                    xi(source_node), eta(source_node)]
                reference_weight = weights(target_node)*weights(source_node)/ &
                    (4.0_dp*pi)
                call pair_vjp_term( &
                    target_point, target_tangent_xi, target_tangent_eta, &
                    target_jacobian, source_point, source_tangent_xi, &
                    source_tangent_eta, source_jacobian, reference_weight, &
                    source_barycentric, single_bar, double_bar, &
                    target_point_bar, target_tangent_xi_bar, &
                    target_tangent_eta_bar, target_jacobian_bar, &
                    source_point_bar, source_tangent_xi_bar, &
                    source_tangent_eta_bar, source_jacobian_bar, status)
                if (status /= 0) return
                call evaluate_torus_curved_panel_vjp( &
                    target_parameters, major_radius, minor_radius, &
                    xi(target_node), eta(target_node), target_point_bar, &
                    target_tangent_xi_bar, target_tangent_eta_bar, &
                    target_jacobian_bar, local_parameters_bar, &
                    local_major_radius_bar, local_minor_radius_bar, xi_bar, &
                    eta_bar, status)
                if (status /= 0) return
                target_parameters_bar = target_parameters_bar + local_parameters_bar
                target_major_radius_bar = target_major_radius_bar + &
                    local_major_radius_bar
                target_minor_radius_bar = target_minor_radius_bar + &
                    local_minor_radius_bar
                call evaluate_torus_curved_panel_vjp( &
                    source_parameters, major_radius, minor_radius, &
                    xi(source_node), eta(source_node), source_point_bar, &
                    source_tangent_xi_bar, source_tangent_eta_bar, &
                    source_jacobian_bar, local_parameters_bar, &
                    local_major_radius_bar, local_minor_radius_bar, xi_bar, &
                    eta_bar, status)
                if (status /= 0) return
                source_parameters_bar = source_parameters_bar + local_parameters_bar
                source_major_radius_bar = source_major_radius_bar + &
                    local_major_radius_bar
                source_minor_radius_bar = source_minor_radius_bar + &
                    local_minor_radius_bar
            end do
        end do
        status = 0
    end subroutine integrate_regular_pair_vjp

    subroutine integrate_coincident_pair_vjp( &
            panel_parameters, major_radius, minor_radius, xi, eta, weights, &
            line_nodes, line_weights, single_bar, double_bar, parameters_bar, &
            major_radius_bar, minor_radius_bar, status)
        real(dp), intent(in) :: panel_parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: xi(:), eta(:), weights(:), line_nodes(:)
        real(dp), intent(in) :: line_weights(:), single_bar, double_bar(3)
        real(dp), intent(out) :: parameters_bar(2, 3)
        real(dp), intent(out) :: major_radius_bar, minor_radius_bar
        integer, intent(out) :: status

        real(dp) :: direction(2), first_reference(2), second_reference(2)
        real(dp) :: wedge_first(2), wedge_second(2), wedge_jacobian, rho, t
        real(dp) :: target_point(3), source_point(3), target_tangent_xi(3)
        real(dp) :: target_tangent_eta(3), source_tangent_xi(3)
        real(dp) :: source_tangent_eta(3), target_jacobian, source_jacobian
        real(dp) :: source_barycentric(3), reference_weight
        real(dp) :: target_point_bar(3), source_point_bar(3)
        real(dp) :: target_tangent_xi_bar(3), target_tangent_eta_bar(3)
        real(dp) :: source_tangent_xi_bar(3), source_tangent_eta_bar(3)
        real(dp) :: target_jacobian_bar, source_jacobian_bar
        real(dp) :: local_parameters_bar(2, 3), local_major_radius_bar
        real(dp) :: local_minor_radius_bar, xi_bar, eta_bar
        integer :: radial_node, tangential_node, target_node, wedge

        parameters_bar = 0.0_dp
        major_radius_bar = 0.0_dp
        minor_radius_bar = 0.0_dp
        status = 1
        do target_node = 1, size(weights)
            first_reference = [xi(target_node), eta(target_node)]
            call evaluate_torus_curved_panel( &
                panel_parameters, major_radius, minor_radius, &
                first_reference(1), first_reference(2), target_point, &
                target_tangent_xi, target_tangent_eta, target_jacobian, &
                status)
            if (status /= 0) return
            do wedge = 1, 3
                wedge_first = reference_vertices(:, wedge) - first_reference
                wedge_second = reference_vertices(:, modulo(wedge, 3) + 1) - &
                    first_reference
                wedge_jacobian = abs( &
                    wedge_first(1)*wedge_second(2) - &
                    wedge_first(2)*wedge_second(1))
                do radial_node = 1, size(line_nodes)
                    rho = line_nodes(radial_node)
                    do tangential_node = 1, size(line_nodes)
                        t = line_nodes(tangential_node)
                        direction = (1.0_dp - t)*wedge_first + t*wedge_second
                        second_reference = first_reference + rho*direction
                        call evaluate_torus_curved_panel( &
                            panel_parameters, major_radius, minor_radius, &
                            second_reference(1), second_reference(2), &
                            source_point, source_tangent_xi, &
                            source_tangent_eta, source_jacobian, status)
                        if (status /= 0) return
                        source_barycentric = [ &
                            1.0_dp - sum(second_reference), second_reference(1), &
                            second_reference(2)]
                        reference_weight = weights(target_node)* &
                            line_weights(radial_node)*line_weights(tangential_node)* &
                            rho*wedge_jacobian/(4.0_dp*pi)
                        call pair_vjp_term( &
                            target_point, target_tangent_xi, target_tangent_eta, &
                            target_jacobian, source_point, source_tangent_xi, &
                            source_tangent_eta, source_jacobian, reference_weight, &
                            source_barycentric, single_bar, double_bar, &
                            target_point_bar, target_tangent_xi_bar, &
                            target_tangent_eta_bar, target_jacobian_bar, &
                            source_point_bar, source_tangent_xi_bar, &
                            source_tangent_eta_bar, source_jacobian_bar, status)
                        if (status /= 0) return
                        call evaluate_torus_curved_panel_vjp( &
                            panel_parameters, major_radius, minor_radius, &
                            first_reference(1), first_reference(2), &
                            target_point_bar, target_tangent_xi_bar, &
                            target_tangent_eta_bar, target_jacobian_bar, &
                            local_parameters_bar, local_major_radius_bar, &
                            local_minor_radius_bar, xi_bar, eta_bar, status)
                        if (status /= 0) return
                        parameters_bar = parameters_bar + local_parameters_bar
                        major_radius_bar = major_radius_bar + local_major_radius_bar
                        minor_radius_bar = minor_radius_bar + local_minor_radius_bar
                        call evaluate_torus_curved_panel_vjp( &
                            panel_parameters, major_radius, minor_radius, &
                            second_reference(1), second_reference(2), &
                            source_point_bar, source_tangent_xi_bar, &
                            source_tangent_eta_bar, source_jacobian_bar, &
                            local_parameters_bar, local_major_radius_bar, &
                            local_minor_radius_bar, xi_bar, eta_bar, status)
                        if (status /= 0) return
                        parameters_bar = parameters_bar + local_parameters_bar
                        major_radius_bar = major_radius_bar + local_major_radius_bar
                        minor_radius_bar = minor_radius_bar + local_minor_radius_bar
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine integrate_coincident_pair_vjp

    subroutine pair_vjp_term( &
            target_point, target_tangent_xi, target_tangent_eta, &
            target_jacobian, source_point, source_tangent_xi, &
            source_tangent_eta, source_jacobian, reference_weight, &
            source_barycentric, single_bar, double_bar, target_point_bar, &
            target_tangent_xi_bar, target_tangent_eta_bar, &
            target_jacobian_bar, source_point_bar, source_tangent_xi_bar, &
            source_tangent_eta_bar, source_jacobian_bar, status)
        real(dp), intent(in) :: target_point(3), target_tangent_xi(3)
        real(dp), intent(in) :: target_tangent_eta(3), target_jacobian
        real(dp), intent(in) :: source_point(3), source_tangent_xi(3)
        real(dp), intent(in) :: source_tangent_eta(3), source_jacobian
        real(dp), intent(in) :: reference_weight, source_barycentric(3)
        real(dp), intent(in) :: single_bar, double_bar(3)
        real(dp), intent(out) :: target_point_bar(3)
        real(dp), intent(out) :: target_tangent_xi_bar(3)
        real(dp), intent(out) :: target_tangent_eta_bar(3)
        real(dp), intent(out) :: target_jacobian_bar, source_point_bar(3)
        real(dp), intent(out) :: source_tangent_xi_bar(3)
        real(dp), intent(out) :: source_tangent_eta_bar(3)
        real(dp), intent(out) :: source_jacobian_bar
        integer, intent(out) :: status

        real(dp) :: displacement(3), radius, source_normal(3), area_cross(3)
        real(dp) :: dot_normal, double_weight, weight, weight_bar, radius_bar
        real(dp) :: a_bar, displacement_bar(3), normal_bar(3)
        real(dp) :: cross_bar(3)

        target_point_bar = 0.0_dp
        target_tangent_xi_bar = 0.0_dp
        target_tangent_eta_bar = 0.0_dp
        target_jacobian_bar = 0.0_dp
        source_point_bar = 0.0_dp
        source_tangent_xi_bar = 0.0_dp
        source_tangent_eta_bar = 0.0_dp
        source_jacobian_bar = 0.0_dp
        displacement = target_point - source_point
        radius = norm2(displacement)
        status = 1
        if (radius <= tiny(1.0_dp)) return
        area_cross = cross_product(source_tangent_xi, source_tangent_eta)
        source_normal = area_cross/source_jacobian
        dot_normal = dot_product(displacement, source_normal)
        weight = reference_weight*target_jacobian*source_jacobian
        double_weight = dot_product(double_bar, source_barycentric)
        weight_bar = single_bar/radius + double_weight*dot_normal/radius**3
        radius_bar = -single_bar*weight/radius**2 - &
            3.0_dp*double_weight*weight*dot_normal/radius**4
        a_bar = double_weight*weight/radius**3
        displacement_bar = radius_bar*displacement/radius + a_bar*source_normal
        normal_bar = a_bar*displacement
        target_point_bar = displacement_bar
        source_point_bar = -displacement_bar
        target_jacobian_bar = weight_bar*reference_weight*source_jacobian
        source_jacobian_bar = weight_bar*reference_weight*target_jacobian
        cross_bar = normal_bar/source_jacobian
        source_jacobian_bar = source_jacobian_bar - dot_product(normal_bar, &
            area_cross)/(source_jacobian**2)
        source_tangent_xi_bar = cross_product(source_tangent_eta, cross_bar)
        source_tangent_eta_bar = cross_product(cross_bar, source_tangent_xi)
        status = 0
    end subroutine pair_vjp_term

    logical function valid_dtn_inputs( &
            parameters, triangles, major_radius, minor_radius, quadrature_degree)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: quadrature_degree

        valid_dtn_inputs = .false.
        if (size(parameters, 1) /= 2 .or. size(parameters, 2) < 3) return
        if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
        if (any(triangles < 1) .or. &
            any(triangles > size(parameters, 2))) return
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        if (quadrature_degree < 0) return
        valid_dtn_inputs = .true.
    end function valid_dtn_inputs

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_laplace_torus_curved_bem_ad_3d
