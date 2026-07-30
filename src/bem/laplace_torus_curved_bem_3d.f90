module fortfem_laplace_torus_curved_bem_3d
    use fortfem_kinds, only: dp
    use fortfem_torus_curved_panel, only: evaluate_torus_curved_panel
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_linalg, only: dense_solve
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private

    public :: assemble_laplace_torus_curved_dtn_3d
    public :: assemble_laplace_torus_curved_calderon_3d
    public :: solve_laplace_bem_dtn_torus_curved_3d

contains

    subroutine assemble_laplace_torus_curved_calderon_3d( &
            parameters, triangles, major_radius, minor_radius, &
            quadrature_degree, single_layer, double_layer, &
            adjoint_double_layer, hypersingular, mass, status)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: quadrature_degree
        real(dp), allocatable, intent(out) :: single_layer(:, :)
        real(dp), allocatable, intent(out) :: double_layer(:, :)
        real(dp), allocatable, intent(out) :: adjoint_double_layer(:, :)
        real(dp), allocatable, intent(out) :: hypersingular(:, :)
        real(dp), allocatable, intent(out) :: mass(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), line_nodes(:), line_weights(:)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp) :: local_hypersingular(3, 3)
        integer :: first_panel, line_count, quadrature_status, second_panel

        status = 1
        if (allocated(adjoint_double_layer)) deallocate(adjoint_double_layer)
        if (allocated(hypersingular)) deallocate(hypersingular)
        call assemble_laplace_torus_curved_dtn_3d( &
            parameters, triangles, major_radius, minor_radius, &
            quadrature_degree, single_layer, double_layer, mass, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) then
            status = 1
            return
        end if
        line_count = max(2, (quadrature_degree + 3)/2)
        allocate(line_nodes(line_count), line_weights(line_count))
        call gauss_legendre_ab( &
            line_count, 0.0_dp, 1.0_dp, line_nodes, line_weights)
        allocate(adjoint_double_layer(size(parameters, 2), size(triangles, 2)))
        allocate(hypersingular(size(parameters, 2), size(parameters, 2)))
        adjoint_double_layer = transpose(double_layer)
        hypersingular = 0.0_dp
        do first_panel = 1, size(triangles, 2)
            do second_panel = 1, size(triangles, 2)
                if (first_panel == second_panel) then
                    call integrate_coincident_hypersingular_pair( &
                        parameters(:, triangles(:, first_panel)), &
                        major_radius, minor_radius, xi, eta, weights, &
                        line_nodes, line_weights, local_hypersingular, status)
                else
                    call integrate_regular_hypersingular_pair( &
                        parameters(:, triangles(:, first_panel)), &
                        parameters(:, triangles(:, second_panel)), &
                        major_radius, minor_radius, xi, eta, weights, &
                        local_hypersingular, status)
                end if
                if (status /= 0) return
                hypersingular( &
                    triangles(:, first_panel), triangles(:, second_panel)) = &
                    hypersingular( &
                    triangles(:, first_panel), triangles(:, second_panel)) + &
                    local_hypersingular
            end do
        end do
        hypersingular = 0.5_dp*(hypersingular + transpose(hypersingular))
        status = 0
    end subroutine assemble_laplace_torus_curved_calderon_3d

    subroutine solve_laplace_bem_dtn_torus_curved_3d( &
            parameters, triangles, major_radius, minor_radius, &
            dirichlet_trace, quadrature_degree, neumann_trace, status)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: dirichlet_trace(:)
        integer, intent(in) :: quadrature_degree
        real(dp), allocatable, intent(out) :: neumann_trace(:)
        integer, intent(out) :: status

        real(dp), allocatable :: double_layer(:, :), mass(:, :)
        real(dp), allocatable :: right_hand_side(:), single_layer(:, :)
        integer :: info, panel_count

        status = 1
        if (allocated(neumann_trace)) deallocate(neumann_trace)
        call assemble_laplace_torus_curved_dtn_3d( &
            parameters, triangles, major_radius, minor_radius, &
            quadrature_degree, single_layer, double_layer, mass, status)
        if (status /= 0) return
        if (size(dirichlet_trace) /= size(parameters, 2)) then
            status = 1
            return
        end if
        panel_count = size(triangles, 2)
        allocate(right_hand_side(panel_count), neumann_trace(panel_count))
        right_hand_side = matmul( &
            double_layer - 0.5_dp*mass, dirichlet_trace)
        call dense_solve( &
            single_layer, right_hand_side, neumann_trace, info)
        if (info /= 0) then
            deallocate(neumann_trace)
            status = 2
            return
        end if
        status = 0
    end subroutine solve_laplace_bem_dtn_torus_curved_3d

    subroutine assemble_laplace_torus_curved_dtn_3d( &
            parameters, triangles, major_radius, minor_radius, &
            quadrature_degree, single_layer, double_layer, mass, status)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: quadrature_degree
        real(dp), allocatable, intent(out) :: single_layer(:, :)
        real(dp), allocatable, intent(out) :: double_layer(:, :)
        real(dp), allocatable, intent(out) :: mass(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), line_nodes(:), line_weights(:)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp) :: local_double(3), local_mass(3), local_single
        integer :: first_panel, line_count, quadrature_status, second_panel

        status = 1
        if (allocated(single_layer)) deallocate(single_layer)
        if (allocated(double_layer)) deallocate(double_layer)
        if (allocated(mass)) deallocate(mass)
        if (size(parameters, 1) /= 2 .or. size(parameters, 2) < 3) return
        if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
        if (any(triangles < 1)) return
        if (any(triangles > size(parameters, 2))) return
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        if (quadrature_degree < 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        line_count = max(2, (quadrature_degree + 3)/2)
        allocate(line_nodes(line_count), line_weights(line_count))
        call gauss_legendre_ab( &
            line_count, 0.0_dp, 1.0_dp, line_nodes, line_weights)
        allocate( &
            single_layer(size(triangles, 2), size(triangles, 2)), &
            double_layer(size(triangles, 2), size(parameters, 2)), &
            mass(size(triangles, 2), size(parameters, 2)))
        single_layer = 0.0_dp
        double_layer = 0.0_dp
        mass = 0.0_dp
        do first_panel = 1, size(triangles, 2)
            call integrate_trace_mass( &
                parameters(:, triangles(:, first_panel)), major_radius, &
                minor_radius, xi, eta, weights, local_mass, status)
            if (status /= 0) return
            mass(first_panel, triangles(:, first_panel)) = local_mass
            do second_panel = 1, size(triangles, 2)
                if (first_panel == second_panel) then
                    call integrate_coincident_pair( &
                        parameters(:, triangles(:, first_panel)), &
                        major_radius, minor_radius, xi, eta, weights, &
                        line_nodes, line_weights, local_single, local_double, &
                        status)
                else
                    call integrate_regular_pair( &
                        parameters(:, triangles(:, first_panel)), &
                        parameters(:, triangles(:, second_panel)), &
                        major_radius, minor_radius, xi, eta, weights, &
                        local_single, local_double, status)
                end if
                if (status /= 0) return
                single_layer(first_panel, second_panel) = local_single
                double_layer(first_panel, triangles(:, second_panel)) = &
                    double_layer(first_panel, triangles(:, second_panel)) + &
                    local_double
            end do
        end do
        status = 0
    end subroutine assemble_laplace_torus_curved_dtn_3d

    subroutine integrate_trace_mass( &
            panel_parameters, major_radius, minor_radius, xi, eta, weights, &
            mass, status)
        real(dp), intent(in) :: panel_parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        real(dp), intent(out) :: mass(3)
        integer, intent(out) :: status

        real(dp) :: jacobian, point(3), tangent_eta(3), tangent_xi(3)
        integer :: node

        mass = 0.0_dp
        do node = 1, size(weights)
            call evaluate_torus_curved_panel( &
                panel_parameters, major_radius, minor_radius, xi(node), &
                eta(node), point, tangent_xi, tangent_eta, jacobian, status)
            if (status /= 0) return
            mass = mass + weights(node)*jacobian* &
                [1.0_dp - xi(node) - eta(node), xi(node), eta(node)]
        end do
        status = 0
    end subroutine integrate_trace_mass

    subroutine integrate_regular_pair( &
            target_parameters, source_parameters, major_radius, minor_radius, &
            xi, eta, weights, single_layer, double_layer, status)
        real(dp), intent(in) :: target_parameters(2, 3)
        real(dp), intent(in) :: source_parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        real(dp), intent(out) :: single_layer, double_layer(3)
        integer, intent(out) :: status

        real(dp) :: displacement(3), jacobian_source, jacobian_target
        real(dp) :: normal_source(3), point_source(3), point_target(3)
        real(dp) :: radius, source_barycentric(3)
        real(dp) :: source_tangent_eta(3), source_tangent_xi(3)
        real(dp) :: target_tangent_eta(3), target_tangent_xi(3)
        real(dp) :: weighted_kernel
        integer :: source_node, target_node

        single_layer = 0.0_dp
        double_layer = 0.0_dp
        do target_node = 1, size(weights)
            call evaluate_torus_curved_panel( &
                target_parameters, major_radius, minor_radius, &
                xi(target_node), eta(target_node), point_target, &
                target_tangent_xi, target_tangent_eta, jacobian_target, status)
            if (status /= 0) return
            do source_node = 1, size(weights)
                call evaluate_torus_curved_panel( &
                    source_parameters, major_radius, minor_radius, &
                    xi(source_node), eta(source_node), point_source, &
                    source_tangent_xi, source_tangent_eta, jacobian_source, &
                    status)
                if (status /= 0) return
                normal_source = cross_product( &
                    source_tangent_xi, source_tangent_eta)/jacobian_source
                displacement = point_target - point_source
                radius = norm2(displacement)
                if (radius <= tiny(1.0_dp)) then
                    status = 1
                    return
                end if
                weighted_kernel = weights(target_node)*weights(source_node)* &
                    jacobian_target*jacobian_source/(4.0_dp*acos(-1.0_dp))
                single_layer = single_layer + weighted_kernel/radius
                source_barycentric = [ &
                    1.0_dp - xi(source_node) - eta(source_node), &
                    xi(source_node), eta(source_node)]
                double_layer = double_layer + weighted_kernel* &
                    dot_product(displacement, normal_source)/radius**3* &
                    source_barycentric
            end do
        end do
        status = 0
    end subroutine integrate_regular_pair

    subroutine integrate_coincident_pair( &
            panel_parameters, major_radius, minor_radius, xi, eta, weights, &
            line_nodes, line_weights, single_layer, double_layer, status)
        real(dp), intent(in) :: panel_parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        real(dp), intent(in) :: line_nodes(:), line_weights(:)
        real(dp), intent(out) :: single_layer, double_layer(3)
        integer, intent(out) :: status

        real(dp), parameter :: reference_vertices(2, 3) = reshape( &
            [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        real(dp) :: direction(2), displacement(3), first_reference(2)
        real(dp) :: jacobian_source, jacobian_target, normal_source(3)
        real(dp) :: point_source(3), point_target(3), radius, rho
        real(dp) :: second_reference(2), source_barycentric(3), t
        real(dp) :: source_tangent_eta(3), source_tangent_xi(3)
        real(dp) :: target_tangent_eta(3), target_tangent_xi(3)
        real(dp) :: wedge_first(2), wedge_jacobian, wedge_second(2)
        real(dp) :: weighted_kernel
        integer :: radial_node, tangential_node, target_node, wedge

        single_layer = 0.0_dp
        double_layer = 0.0_dp
        do target_node = 1, size(weights)
            first_reference = [xi(target_node), eta(target_node)]
            call evaluate_torus_curved_panel( &
                panel_parameters, major_radius, minor_radius, &
                first_reference(1), first_reference(2), point_target, &
                target_tangent_xi, target_tangent_eta, jacobian_target, status)
            if (status /= 0) return
            do wedge = 1, 3
                wedge_first = reference_vertices(:, wedge) - first_reference
                wedge_second = &
                    reference_vertices(:, modulo(wedge, 3) + 1) - &
                    first_reference
                wedge_jacobian = abs( &
                    wedge_first(1)*wedge_second(2) - &
                    wedge_first(2)*wedge_second(1))
                do radial_node = 1, size(line_nodes)
                    rho = line_nodes(radial_node)
                    do tangential_node = 1, size(line_nodes)
                        t = line_nodes(tangential_node)
                        direction = (1.0_dp - t)*wedge_first + &
                            t*wedge_second
                        second_reference = first_reference + rho*direction
                        call evaluate_torus_curved_panel( &
                            panel_parameters, major_radius, minor_radius, &
                            second_reference(1), second_reference(2), &
                            point_source, source_tangent_xi, &
                            source_tangent_eta, jacobian_source, status)
                        if (status /= 0) return
                        normal_source = cross_product( &
                            source_tangent_xi, source_tangent_eta)/ &
                            jacobian_source
                        displacement = point_target - point_source
                        radius = norm2(displacement)
                        if (radius <= tiny(1.0_dp)) then
                            status = 1
                            return
                        end if
                        weighted_kernel = weights(target_node)* &
                            line_weights(radial_node)* &
                            line_weights(tangential_node)*rho*wedge_jacobian* &
                            jacobian_target*jacobian_source/ &
                            (4.0_dp*acos(-1.0_dp))
                        single_layer = single_layer + weighted_kernel/radius
                        source_barycentric = [ &
                            1.0_dp - sum(second_reference), &
                            second_reference(1), second_reference(2)]
                        double_layer = double_layer + weighted_kernel* &
                            dot_product(displacement, normal_source)/ &
                            radius**3*source_barycentric
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine integrate_coincident_pair

    subroutine integrate_regular_hypersingular_pair( &
            target_parameters, source_parameters, major_radius, minor_radius, &
            xi, eta, weights, hypersingular, status)
        real(dp), intent(in) :: target_parameters(2, 3)
        real(dp), intent(in) :: source_parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        real(dp), intent(out) :: hypersingular(3, 3)
        integer, intent(out) :: status

        real(dp) :: curls_source(3, 3), curls_target(3, 3)
        real(dp) :: jacobian_source, jacobian_target, point_source(3)
        real(dp) :: point_target(3), radius, source_tangent_eta(3)
        real(dp) :: source_tangent_xi(3), target_tangent_eta(3)
        real(dp) :: target_tangent_xi(3), weighted_kernel
        integer :: source_node, target_node

        hypersingular = 0.0_dp
        do target_node = 1, size(weights)
            call evaluate_torus_curved_panel( &
                target_parameters, major_radius, minor_radius, &
                xi(target_node), eta(target_node), point_target, &
                target_tangent_xi, target_tangent_eta, jacobian_target, status)
            if (status /= 0) return
            call surface_curls( &
                target_tangent_xi, target_tangent_eta, jacobian_target, &
                curls_target, status)
            if (status /= 0) return
            do source_node = 1, size(weights)
                call evaluate_torus_curved_panel( &
                    source_parameters, major_radius, minor_radius, &
                    xi(source_node), eta(source_node), point_source, &
                    source_tangent_xi, source_tangent_eta, jacobian_source, &
                    status)
                if (status /= 0) return
                call surface_curls( &
                    source_tangent_xi, source_tangent_eta, jacobian_source, &
                    curls_source, status)
                if (status /= 0) return
                radius = norm2(point_target - point_source)
                if (radius <= tiny(1.0_dp)) then
                    status = 1
                    return
                end if
                weighted_kernel = weights(target_node)*weights(source_node)* &
                    jacobian_target*jacobian_source/ &
                    (4.0_dp*acos(-1.0_dp)*radius)
                hypersingular = hypersingular + weighted_kernel* &
                    matmul(transpose(curls_target), curls_source)
            end do
        end do
        status = 0
    end subroutine integrate_regular_hypersingular_pair

    subroutine integrate_coincident_hypersingular_pair( &
            panel_parameters, major_radius, minor_radius, xi, eta, weights, &
            line_nodes, line_weights, hypersingular, status)
        real(dp), intent(in) :: panel_parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        real(dp), intent(in) :: line_nodes(:), line_weights(:)
        real(dp), intent(out) :: hypersingular(3, 3)
        integer, intent(out) :: status

        real(dp), parameter :: reference_vertices(2, 3) = reshape( &
            [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        real(dp) :: curls_source(3, 3), curls_target(3, 3), direction(2)
        real(dp) :: first_reference(2), jacobian_source, jacobian_target
        real(dp) :: point_source(3), point_target(3), radius, rho
        real(dp) :: second_reference(2), source_tangent_eta(3)
        real(dp) :: source_tangent_xi(3), t, target_tangent_eta(3)
        real(dp) :: target_tangent_xi(3), wedge_first(2), wedge_jacobian
        real(dp) :: wedge_second(2), weighted_kernel
        integer :: radial_node, tangential_node, target_node, wedge

        hypersingular = 0.0_dp
        do target_node = 1, size(weights)
            first_reference = [xi(target_node), eta(target_node)]
            call evaluate_torus_curved_panel( &
                panel_parameters, major_radius, minor_radius, &
                first_reference(1), first_reference(2), point_target, &
                target_tangent_xi, target_tangent_eta, jacobian_target, status)
            if (status /= 0) return
            call surface_curls( &
                target_tangent_xi, target_tangent_eta, jacobian_target, &
                curls_target, status)
            if (status /= 0) return
            do wedge = 1, 3
                wedge_first = reference_vertices(:, wedge) - first_reference
                wedge_second = &
                    reference_vertices(:, modulo(wedge, 3) + 1) - &
                    first_reference
                wedge_jacobian = abs( &
                    wedge_first(1)*wedge_second(2) - &
                    wedge_first(2)*wedge_second(1))
                do radial_node = 1, size(line_nodes)
                    rho = line_nodes(radial_node)
                    do tangential_node = 1, size(line_nodes)
                        t = line_nodes(tangential_node)
                        direction = (1.0_dp - t)*wedge_first + &
                            t*wedge_second
                        second_reference = first_reference + rho*direction
                        call evaluate_torus_curved_panel( &
                            panel_parameters, major_radius, minor_radius, &
                            second_reference(1), second_reference(2), &
                            point_source, source_tangent_xi, &
                            source_tangent_eta, jacobian_source, status)
                        if (status /= 0) return
                        call surface_curls( &
                            source_tangent_xi, source_tangent_eta, &
                            jacobian_source, curls_source, status)
                        if (status /= 0) return
                        radius = norm2(point_target - point_source)
                        if (radius <= tiny(1.0_dp)) then
                            status = 1
                            return
                        end if
                        weighted_kernel = weights(target_node)* &
                            line_weights(radial_node)* &
                            line_weights(tangential_node)*rho*wedge_jacobian* &
                            jacobian_target*jacobian_source/ &
                            (4.0_dp*acos(-1.0_dp)*radius)
                        hypersingular = hypersingular + weighted_kernel* &
                            matmul(transpose(curls_target), curls_source)
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine integrate_coincident_hypersingular_pair

    pure subroutine surface_curls( &
            tangent_xi, tangent_eta, jacobian, curls, status)
        real(dp), intent(in) :: tangent_xi(3), tangent_eta(3), jacobian
        real(dp), intent(out) :: curls(3, 3)
        integer, intent(out) :: status

        real(dp), parameter :: reference_gradients(2, 3) = reshape( &
            [-1.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        real(dp) :: coefficient_eta, coefficient_xi, determinant
        real(dp) :: gradient(3), metric_eta_eta, metric_xi_eta, metric_xi_xi
        real(dp) :: normal(3)
        integer :: basis

        curls = 0.0_dp
        status = 1
        metric_xi_xi = dot_product(tangent_xi, tangent_xi)
        metric_xi_eta = dot_product(tangent_xi, tangent_eta)
        metric_eta_eta = dot_product(tangent_eta, tangent_eta)
        determinant = &
            metric_xi_xi*metric_eta_eta - metric_xi_eta**2
        if (determinant <= tiny(1.0_dp)) return
        if (jacobian <= tiny(1.0_dp)) return
        normal = cross_product(tangent_xi, tangent_eta)/jacobian
        do basis = 1, 3
            coefficient_xi = ( &
                metric_eta_eta*reference_gradients(1, basis) - &
                metric_xi_eta*reference_gradients(2, basis))/determinant
            coefficient_eta = ( &
                -metric_xi_eta*reference_gradients(1, basis) + &
                metric_xi_xi*reference_gradients(2, basis))/determinant
            gradient = &
                coefficient_xi*tangent_xi + coefficient_eta*tangent_eta
            curls(:, basis) = cross_product(normal, gradient)
        end do
        status = 0
    end subroutine surface_curls

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_laplace_torus_curved_bem_3d
