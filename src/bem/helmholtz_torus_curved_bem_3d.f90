module fortfem_helmholtz_torus_curved_bem_3d
    !! Curved-torus Helmholtz operators by Laplace singularity subtraction.
    !!
    !! The convention G=exp(i*k*r)/(4*pi*r) is outgoing for exp(-i*omega*t).
    !! The singular Laplace part uses radial Duffy quadrature; only the smooth
    !! Helmholtz-minus-Laplace correction is integrated by product quadrature.
    !! The hypersingular block uses the surface-curl regularization of Dolz,
    !! Harbrecht, and Multerer, Engineering with Computers 40 (2024),
    !! https://doi.org/10.1007/s00366-024-02013-y.
    use fortfem_kinds, only: dp
    use fortfem_laplace_torus_curved_bem_3d, only: &
        assemble_laplace_torus_curved_calderon_3d, &
        assemble_laplace_torus_curved_dtn_3d, &
        evaluate_torus_curved_p1_surface_curls
    use fortfem_torus_curved_panel, only: evaluate_torus_curved_panel
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_linalg, only: dense_solve
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private

    public :: assemble_helmholtz_torus_curved_dtn_3d
    public :: assemble_helmholtz_torus_curved_calderon_3d
    public :: solve_helmholtz_bem_dtn_torus_curved_3d

contains

    subroutine assemble_helmholtz_torus_curved_calderon_3d( &
            parameters, triangles, major_radius, minor_radius, wave_number, &
            quadrature_degree, single_layer, double_layer, &
            adjoint_double_layer, hypersingular, mass, status)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        integer, intent(in) :: quadrature_degree
        complex(dp), allocatable, intent(out) :: single_layer(:, :)
        complex(dp), allocatable, intent(out) :: double_layer(:, :)
        complex(dp), allocatable, intent(out) :: adjoint_double_layer(:, :)
        complex(dp), allocatable, intent(out) :: hypersingular(:, :)
        real(dp), allocatable, intent(out) :: mass(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), laplace_adjoint(:, :)
        real(dp), allocatable :: laplace_double(:, :), laplace_hyper(:, :)
        real(dp), allocatable :: laplace_single(:, :), line_nodes(:)
        real(dp), allocatable :: line_weights(:), weights(:), xi(:)
        complex(dp) :: local_double(3), local_hyper(3, 3), local_single
        integer :: first_panel, line_count, quadrature_status, second_panel

        status = 1
        if (allocated(single_layer)) deallocate(single_layer)
        if (allocated(double_layer)) deallocate(double_layer)
        if (allocated(adjoint_double_layer)) deallocate(adjoint_double_layer)
        if (allocated(hypersingular)) deallocate(hypersingular)
        if (allocated(mass)) deallocate(mass)
        if (wave_number <= 0.0_dp) return
        call assemble_laplace_torus_curved_calderon_3d( &
            parameters, triangles, major_radius, minor_radius, &
            quadrature_degree, laplace_single, laplace_double, &
            laplace_adjoint, laplace_hyper, mass, status)
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
        allocate( &
            single_layer(size(laplace_single, 1), size(laplace_single, 2)), &
            double_layer(size(laplace_double, 1), size(laplace_double, 2)), &
            adjoint_double_layer( &
            size(laplace_adjoint, 1), size(laplace_adjoint, 2)), &
            hypersingular(size(laplace_hyper, 1), size(laplace_hyper, 2)))
        single_layer = cmplx(laplace_single, 0.0_dp, dp)
        double_layer = cmplx(laplace_double, 0.0_dp, dp)
        hypersingular = cmplx(laplace_hyper, 0.0_dp, dp)
        do first_panel = 1, size(triangles, 2)
            do second_panel = 1, size(triangles, 2)
                call integrate_helmholtz_correction( &
                    parameters(:, triangles(:, first_panel)), &
                    parameters(:, triangles(:, second_panel)), &
                    major_radius, minor_radius, wave_number, xi, eta, &
                    weights, local_single, local_double, status)
                if (status /= 0) return
                if (first_panel == second_panel) then
                    call integrate_coincident_hypersingular_correction( &
                        parameters(:, triangles(:, first_panel)), &
                        major_radius, minor_radius, wave_number, xi, eta, &
                        weights, line_nodes, line_weights, local_hyper, status)
                else
                    call integrate_regular_hypersingular_correction( &
                        parameters(:, triangles(:, first_panel)), &
                        parameters(:, triangles(:, second_panel)), &
                        major_radius, minor_radius, wave_number, xi, eta, &
                        weights, local_hyper, status)
                end if
                if (status /= 0) return
                single_layer(first_panel, second_panel) = &
                    single_layer(first_panel, second_panel) + local_single
                double_layer(first_panel, triangles(:, second_panel)) = &
                    double_layer(first_panel, triangles(:, second_panel)) + &
                    local_double
                hypersingular( &
                    triangles(:, first_panel), triangles(:, second_panel)) = &
                    hypersingular( &
                    triangles(:, first_panel), triangles(:, second_panel)) + &
                    local_hyper
            end do
        end do
        adjoint_double_layer = transpose(double_layer)
        hypersingular = 0.5_dp*(hypersingular + transpose(hypersingular))
        status = 0
    end subroutine assemble_helmholtz_torus_curved_calderon_3d

    subroutine solve_helmholtz_bem_dtn_torus_curved_3d( &
            parameters, triangles, major_radius, minor_radius, wave_number, &
            dirichlet_trace, quadrature_degree, neumann_trace, status)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        complex(dp), intent(in) :: dirichlet_trace(:)
        integer, intent(in) :: quadrature_degree
        complex(dp), allocatable, intent(out) :: neumann_trace(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: double_layer(:, :)
        complex(dp), allocatable :: right_hand_side(:), single_layer(:, :)
        real(dp), allocatable :: mass(:, :)
        integer :: info, panel_count

        status = 1
        if (allocated(neumann_trace)) deallocate(neumann_trace)
        call assemble_helmholtz_torus_curved_dtn_3d( &
            parameters, triangles, major_radius, minor_radius, wave_number, &
            quadrature_degree, single_layer, double_layer, mass, status)
        if (status /= 0) return
        if (size(dirichlet_trace) /= size(parameters, 2)) return
        panel_count = size(triangles, 2)
        allocate(right_hand_side(panel_count), neumann_trace(panel_count))
        right_hand_side = matmul( &
            double_layer - 0.5_dp*cmplx(mass, 0.0_dp, dp), dirichlet_trace)
        call dense_solve( &
            single_layer, right_hand_side, neumann_trace, info)
        if (info /= 0) then
            deallocate(neumann_trace)
            status = 2
            return
        end if
        status = 0
    end subroutine solve_helmholtz_bem_dtn_torus_curved_3d

    subroutine assemble_helmholtz_torus_curved_dtn_3d( &
            parameters, triangles, major_radius, minor_radius, wave_number, &
            quadrature_degree, single_layer, double_layer, mass, status)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        integer, intent(in) :: quadrature_degree
        complex(dp), allocatable, intent(out) :: single_layer(:, :)
        complex(dp), allocatable, intent(out) :: double_layer(:, :)
        real(dp), allocatable, intent(out) :: mass(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), laplace_double(:, :)
        real(dp), allocatable :: laplace_single(:, :), weights(:), xi(:)
        complex(dp) :: local_double(3), local_single
        integer :: first_panel, quadrature_status, second_panel

        status = 1
        if (allocated(single_layer)) deallocate(single_layer)
        if (allocated(double_layer)) deallocate(double_layer)
        if (allocated(mass)) deallocate(mass)
        if (wave_number <= 0.0_dp) return
        call assemble_laplace_torus_curved_dtn_3d( &
            parameters, triangles, major_radius, minor_radius, &
            quadrature_degree, laplace_single, laplace_double, mass, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) then
            status = 1
            return
        end if
        allocate( &
            single_layer(size(laplace_single, 1), size(laplace_single, 2)), &
            double_layer(size(laplace_double, 1), size(laplace_double, 2)))
        single_layer = cmplx(laplace_single, 0.0_dp, dp)
        double_layer = cmplx(laplace_double, 0.0_dp, dp)
        do first_panel = 1, size(triangles, 2)
            do second_panel = 1, size(triangles, 2)
                call integrate_helmholtz_correction( &
                    parameters(:, triangles(:, first_panel)), &
                    parameters(:, triangles(:, second_panel)), &
                    major_radius, minor_radius, wave_number, xi, eta, &
                    weights, local_single, local_double, status)
                if (status /= 0) return
                single_layer(first_panel, second_panel) = &
                    single_layer(first_panel, second_panel) + local_single
                double_layer(first_panel, triangles(:, second_panel)) = &
                    double_layer(first_panel, triangles(:, second_panel)) + &
                    local_double
            end do
        end do
        status = 0
    end subroutine assemble_helmholtz_torus_curved_dtn_3d

    subroutine integrate_helmholtz_correction( &
            target_parameters, source_parameters, major_radius, minor_radius, &
            wave_number, xi, eta, weights, single_layer, double_layer, status)
        real(dp), intent(in) :: target_parameters(2, 3)
        real(dp), intent(in) :: source_parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        complex(dp), intent(out) :: single_layer, double_layer(3)
        integer, intent(out) :: status

        real(dp) :: displacement(3), jacobian_source, jacobian_target
        real(dp) :: normal_source(3), point_source(3), point_target(3), radius
        real(dp) :: source_barycentric(3), source_tangent_eta(3)
        real(dp) :: source_tangent_xi(3), target_tangent_eta(3)
        real(dp) :: target_tangent_xi(3), weight
        complex(dp) :: exponential, green_correction, normal_correction
        integer :: source_node, target_node

        single_layer = cmplx(0.0_dp, 0.0_dp, dp)
        double_layer = cmplx(0.0_dp, 0.0_dp, dp)
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
                weight = weights(target_node)*weights(source_node)* &
                    jacobian_target*jacobian_source
                if (radius <= sqrt(tiny(1.0_dp))) then
                    green_correction = &
                        cmplx(0.0_dp, wave_number/(4.0_dp*acos(-1.0_dp)), dp)
                    normal_correction = cmplx(0.0_dp, 0.0_dp, dp)
                else
                    exponential = &
                        exp(cmplx(0.0_dp, wave_number*radius, dp))
                    green_correction = (exponential - 1.0_dp)/ &
                        (4.0_dp*acos(-1.0_dp)*radius)
                    normal_correction = &
                        (exponential* &
                        cmplx(1.0_dp, -wave_number*radius, dp) - 1.0_dp)* &
                        dot_product(displacement, normal_source)/ &
                        (4.0_dp*acos(-1.0_dp)*radius**3)
                end if
                source_barycentric = [ &
                    1.0_dp - xi(source_node) - eta(source_node), &
                    xi(source_node), eta(source_node)]
                single_layer = single_layer + weight*green_correction
                double_layer = double_layer + &
                    weight*normal_correction*source_barycentric
            end do
        end do
        status = 0
    end subroutine integrate_helmholtz_correction

    subroutine integrate_regular_hypersingular_correction( &
            target_parameters, source_parameters, major_radius, minor_radius, &
            wave_number, xi, eta, weights, hypersingular, status)
        real(dp), intent(in) :: target_parameters(2, 3)
        real(dp), intent(in) :: source_parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        complex(dp), intent(out) :: hypersingular(3, 3)
        integer, intent(out) :: status

        real(dp) :: barycentric_source(3), barycentric_target(3)
        real(dp) :: curls_source(3, 3), curls_target(3, 3)
        real(dp) :: jacobian_source, jacobian_target, normal_source(3)
        real(dp) :: normal_target(3), point_source(3), point_target(3), radius
        real(dp) :: source_tangent_eta(3), source_tangent_xi(3)
        real(dp) :: target_tangent_eta(3), target_tangent_xi(3), weight
        complex(dp) :: exponential, green, green_correction
        integer :: source_basis, source_node, target_basis, target_node

        hypersingular = cmplx(0.0_dp, 0.0_dp, dp)
        do target_node = 1, size(weights)
            call evaluate_torus_curved_panel( &
                target_parameters, major_radius, minor_radius, &
                xi(target_node), eta(target_node), point_target, &
                target_tangent_xi, target_tangent_eta, jacobian_target, status)
            if (status /= 0) return
            normal_target = cross_product( &
                target_tangent_xi, target_tangent_eta)/jacobian_target
            call evaluate_torus_curved_p1_surface_curls( &
                target_tangent_xi, target_tangent_eta, jacobian_target, &
                curls_target, status)
            if (status /= 0) return
            barycentric_target = [ &
                1.0_dp - xi(target_node) - eta(target_node), &
                xi(target_node), eta(target_node)]
            do source_node = 1, size(weights)
                call evaluate_torus_curved_panel( &
                    source_parameters, major_radius, minor_radius, &
                    xi(source_node), eta(source_node), point_source, &
                    source_tangent_xi, source_tangent_eta, jacobian_source, &
                    status)
                if (status /= 0) return
                normal_source = cross_product( &
                    source_tangent_xi, source_tangent_eta)/jacobian_source
                call evaluate_torus_curved_p1_surface_curls( &
                    source_tangent_xi, source_tangent_eta, jacobian_source, &
                    curls_source, status)
                if (status /= 0) return
                radius = norm2(point_target - point_source)
                if (radius <= tiny(1.0_dp)) then
                    status = 1
                    return
                end if
                exponential = exp(cmplx(0.0_dp, wave_number*radius, dp))
                green = exponential/(4.0_dp*acos(-1.0_dp)*radius)
                green_correction = &
                    (exponential - 1.0_dp)/(4.0_dp*acos(-1.0_dp)*radius)
                weight = weights(target_node)*weights(source_node)* &
                    jacobian_target*jacobian_source
                barycentric_source = [ &
                    1.0_dp - xi(source_node) - eta(source_node), &
                    xi(source_node), eta(source_node)]
                hypersingular = hypersingular + weight*green_correction* &
                    matmul(transpose(curls_target), curls_source)
                do target_basis = 1, 3
                    do source_basis = 1, 3
                        hypersingular(target_basis, source_basis) = &
                            hypersingular(target_basis, source_basis) - &
                            weight*wave_number**2*green* &
                            dot_product(normal_target, normal_source)* &
                            barycentric_target(target_basis)* &
                            barycentric_source(source_basis)
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine integrate_regular_hypersingular_correction

    subroutine integrate_coincident_hypersingular_correction( &
            panel_parameters, major_radius, minor_radius, wave_number, &
            xi, eta, weights, line_nodes, line_weights, hypersingular, status)
        real(dp), intent(in) :: panel_parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        real(dp), intent(in) :: line_nodes(:), line_weights(:)
        complex(dp), intent(out) :: hypersingular(3, 3)
        integer, intent(out) :: status

        real(dp), parameter :: reference_vertices(2, 3) = reshape( &
            [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        real(dp) :: barycentric_source(3), barycentric_target(3)
        real(dp) :: curls_source(3, 3), curls_target(3, 3), direction(2)
        real(dp) :: first_reference(2), jacobian_source, jacobian_target
        real(dp) :: normal_source(3), normal_target(3), point_source(3)
        real(dp) :: point_target(3), radius, rho, second_reference(2)
        real(dp) :: source_tangent_eta(3), source_tangent_xi(3), t
        real(dp) :: target_tangent_eta(3), target_tangent_xi(3)
        real(dp) :: wedge_first(2), wedge_jacobian, wedge_second(2), weight
        complex(dp) :: exponential, green, green_correction
        integer :: radial_node, source_basis, tangential_node, target_basis
        integer :: target_node, wedge

        hypersingular = cmplx(0.0_dp, 0.0_dp, dp)
        do target_node = 1, size(weights)
            first_reference = [xi(target_node), eta(target_node)]
            barycentric_target = [ &
                1.0_dp - sum(first_reference), first_reference]
            call evaluate_torus_curved_panel( &
                panel_parameters, major_radius, minor_radius, &
                first_reference(1), first_reference(2), point_target, &
                target_tangent_xi, target_tangent_eta, jacobian_target, status)
            if (status /= 0) return
            normal_target = cross_product( &
                target_tangent_xi, target_tangent_eta)/jacobian_target
            call evaluate_torus_curved_p1_surface_curls( &
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
                        barycentric_source = [ &
                            1.0_dp - sum(second_reference), second_reference]
                        call evaluate_torus_curved_panel( &
                            panel_parameters, major_radius, minor_radius, &
                            second_reference(1), second_reference(2), &
                            point_source, source_tangent_xi, &
                            source_tangent_eta, jacobian_source, status)
                        if (status /= 0) return
                        normal_source = cross_product( &
                            source_tangent_xi, source_tangent_eta)/ &
                            jacobian_source
                        call evaluate_torus_curved_p1_surface_curls( &
                            source_tangent_xi, source_tangent_eta, &
                            jacobian_source, curls_source, status)
                        if (status /= 0) return
                        radius = norm2(point_target - point_source)
                        if (radius <= tiny(1.0_dp)) then
                            status = 1
                            return
                        end if
                        exponential = &
                            exp(cmplx(0.0_dp, wave_number*radius, dp))
                        green = exponential/ &
                            (4.0_dp*acos(-1.0_dp)*radius)
                        green_correction = (exponential - 1.0_dp)/ &
                            (4.0_dp*acos(-1.0_dp)*radius)
                        weight = weights(target_node)* &
                            line_weights(radial_node)* &
                            line_weights(tangential_node)*rho*wedge_jacobian* &
                            jacobian_target*jacobian_source
                        hypersingular = hypersingular + &
                            weight*green_correction* &
                            matmul(transpose(curls_target), curls_source)
                        do target_basis = 1, 3
                            do source_basis = 1, 3
                                hypersingular(target_basis, source_basis) = &
                                    hypersingular(target_basis, source_basis) - &
                                    weight*wave_number**2*green* &
                                    dot_product(normal_target, normal_source)* &
                                    barycentric_target(target_basis)* &
                                    barycentric_source(source_basis)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine integrate_coincident_hypersingular_correction

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_helmholtz_torus_curved_bem_3d
