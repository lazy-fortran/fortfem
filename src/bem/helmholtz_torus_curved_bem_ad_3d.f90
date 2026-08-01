module fortfem_helmholtz_torus_curved_bem_ad_3d
    !! Geometry products for the curved-torus Helmholtz DtN blocks.
    !!
    !! The Laplace singular blocks are differentiated by the existing
    !! FortSym-backed products.  This module differentiates the smooth
    !! Helmholtz-minus-Laplace correction at the same product quadrature
    !! points as the primal assembly, including the coincident-point limit.
    use fortfem_kinds, only: dp
    use fortfem_laplace_torus_curved_bem_ad_3d, only: &
        assemble_laplace_torus_curved_dtn_3d_geometry_jvp, &
        assemble_laplace_torus_curved_dtn_3d_geometry_vjp
    use fortfem_torus_curved_panel, only: &
        evaluate_torus_curved_panel, evaluate_torus_curved_panel_jvp, &
        evaluate_torus_curved_panel_vjp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: four_pi = 4.0_dp*pi

    public :: assemble_helmholtz_torus_curved_dtn_3d_geometry_jvp
    public :: assemble_helmholtz_torus_curved_dtn_3d_geometry_vjp

contains

    subroutine assemble_helmholtz_torus_curved_dtn_3d_geometry_jvp( &
            parameters, triangles, major_radius, minor_radius, wave_number, &
            quadrature_degree, parameters_dot, major_radius_dot, &
            minor_radius_dot, wave_number_dot, single_layer_dot, &
            double_layer_dot, mass_dot, status)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        integer, intent(in) :: quadrature_degree
        real(dp), intent(in) :: parameters_dot(:, :)
        real(dp), intent(in) :: major_radius_dot, minor_radius_dot
        real(dp), intent(in) :: wave_number_dot
        complex(dp), allocatable, intent(out) :: single_layer_dot(:, :)
        complex(dp), allocatable, intent(out) :: double_layer_dot(:, :)
        real(dp), allocatable, intent(out) :: mass_dot(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp), allocatable :: laplace_single_dot(:, :)
        real(dp), allocatable :: laplace_double_dot(:, :)
        real(dp), allocatable :: laplace_mass_dot(:, :)
        complex(dp) :: local_double_dot(3), local_single_dot
        integer :: first_panel, quadrature_status, second_panel

        status = 1
        if (allocated(single_layer_dot)) deallocate(single_layer_dot)
        if (allocated(double_layer_dot)) deallocate(double_layer_dot)
        if (allocated(mass_dot)) deallocate(mass_dot)
        if (.not. valid_inputs( &
            parameters, triangles, major_radius, minor_radius, wave_number, &
            quadrature_degree)) return
        if (any(shape(parameters_dot) /= shape(parameters))) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        call assemble_laplace_torus_curved_dtn_3d_geometry_jvp( &
            parameters, triangles, major_radius, minor_radius, &
            quadrature_degree, parameters_dot, major_radius_dot, &
            minor_radius_dot, laplace_single_dot, laplace_double_dot, &
            laplace_mass_dot, status)
        if (status /= 0) return
        allocate( &
            single_layer_dot(size(triangles, 2), size(triangles, 2)), &
            double_layer_dot(size(triangles, 2), size(parameters, 2)), &
            mass_dot(size(triangles, 2), size(parameters, 2)))
        single_layer_dot = cmplx(laplace_single_dot, 0.0_dp, dp)
        double_layer_dot = cmplx(laplace_double_dot, 0.0_dp, dp)
        mass_dot = laplace_mass_dot
        do first_panel = 1, size(triangles, 2)
            do second_panel = 1, size(triangles, 2)
                call integrate_correction_jvp( &
                    parameters(:, triangles(:, first_panel)), &
                    parameters(:, triangles(:, second_panel)), &
                    parameters_dot(:, triangles(:, first_panel)), &
                    parameters_dot(:, triangles(:, second_panel)), &
                    major_radius, minor_radius, wave_number, &
                    major_radius_dot, minor_radius_dot, wave_number_dot, &
                    xi, eta, weights, local_single_dot, local_double_dot, &
                    status)
                if (status /= 0) return
                single_layer_dot(first_panel, second_panel) = &
                    single_layer_dot(first_panel, second_panel) + &
                    local_single_dot
                double_layer_dot(first_panel, triangles(:, second_panel)) = &
                    double_layer_dot(first_panel, triangles(:, second_panel)) + &
                    local_double_dot
            end do
        end do
        status = 0
    end subroutine assemble_helmholtz_torus_curved_dtn_3d_geometry_jvp

    subroutine assemble_helmholtz_torus_curved_dtn_3d_geometry_vjp( &
            parameters, triangles, major_radius, minor_radius, wave_number, &
            quadrature_degree, single_layer_bar, double_layer_bar, mass_bar, &
            parameters_bar, major_radius_bar, minor_radius_bar, &
            wave_number_bar, status)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        integer, intent(in) :: quadrature_degree
        complex(dp), intent(in) :: single_layer_bar(:, :)
        complex(dp), intent(in) :: double_layer_bar(:, :)
        real(dp), intent(in) :: mass_bar(:, :)
        real(dp), intent(out) :: parameters_bar(:, :)
        real(dp), intent(out) :: major_radius_bar, minor_radius_bar
        real(dp), intent(out) :: wave_number_bar
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp), allocatable :: laplace_single_bar(:, :)
        real(dp), allocatable :: laplace_double_bar(:, :), laplace_mass_bar(:, :)
        real(dp) :: local_major_bar, local_minor_bar, local_wave_bar
        real(dp) :: local_parameters_bar(2, 3), local_parameters_bar_second(2, 3)
        complex(dp) :: local_single_bar, local_double_bar(3)
        integer :: first_panel, quadrature_status, second_panel, local_node

        parameters_bar = 0.0_dp
        major_radius_bar = 0.0_dp
        minor_radius_bar = 0.0_dp
        wave_number_bar = 0.0_dp
        status = 1
        if (.not. valid_inputs( &
            parameters, triangles, major_radius, minor_radius, wave_number, &
            quadrature_degree)) return
        if (size(single_layer_bar, 1) /= size(triangles, 2) .or. &
            size(single_layer_bar, 2) /= size(triangles, 2)) return
        if (size(double_layer_bar, 1) /= size(triangles, 2) .or. &
            size(double_layer_bar, 2) /= size(parameters, 2)) return
        if (size(mass_bar, 1) /= size(triangles, 2) .or. &
            size(mass_bar, 2) /= size(parameters, 2)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        allocate( &
            laplace_single_bar(size(triangles, 2), size(triangles, 2)), &
            laplace_double_bar(size(triangles, 2), size(parameters, 2)), &
            laplace_mass_bar(size(triangles, 2), size(parameters, 2)))
        laplace_single_bar = real(single_layer_bar, dp)
        laplace_double_bar = real(double_layer_bar, dp)
        laplace_mass_bar = mass_bar
        call assemble_laplace_torus_curved_dtn_3d_geometry_vjp( &
            parameters, triangles, major_radius, minor_radius, quadrature_degree, &
            laplace_single_bar, laplace_double_bar, laplace_mass_bar, &
            parameters_bar, major_radius_bar, minor_radius_bar, status)
        if (status /= 0) return
        do first_panel = 1, size(triangles, 2)
            do second_panel = 1, size(triangles, 2)
                local_single_bar = single_layer_bar(first_panel, second_panel)
                local_double_bar = double_layer_bar( &
                    first_panel, triangles(:, second_panel))
                call integrate_correction_vjp( &
                    parameters(:, triangles(:, first_panel)), &
                    parameters(:, triangles(:, second_panel)), major_radius, &
                    minor_radius, wave_number, xi, eta, weights, &
                    local_single_bar, local_double_bar, local_parameters_bar, &
                    local_parameters_bar_second, local_major_bar, &
                    local_minor_bar, local_wave_bar, status)
                if (status /= 0) return
                do local_node = 1, 3
                    parameters_bar(:, triangles(local_node, first_panel)) = &
                        parameters_bar(:, triangles(local_node, first_panel)) + &
                        local_parameters_bar(:, local_node)
                    parameters_bar(:, triangles(local_node, second_panel)) = &
                        parameters_bar(:, triangles(local_node, second_panel)) + &
                        local_parameters_bar_second(:, local_node)
                end do
                major_radius_bar = major_radius_bar + local_major_bar
                minor_radius_bar = minor_radius_bar + local_minor_bar
                wave_number_bar = wave_number_bar + local_wave_bar
            end do
        end do
        status = 0
    end subroutine assemble_helmholtz_torus_curved_dtn_3d_geometry_vjp

    subroutine integrate_correction_jvp( &
            first_parameters, second_parameters, first_parameters_dot, &
            second_parameters_dot, major_radius, minor_radius, wave_number, &
            major_radius_dot, minor_radius_dot, wave_number_dot, xi, eta, &
            weights, single_layer_dot, double_layer_dot, status)
        real(dp), intent(in) :: first_parameters(2, 3), second_parameters(2, 3)
        real(dp), intent(in) :: first_parameters_dot(2, 3)
        real(dp), intent(in) :: second_parameters_dot(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        real(dp), intent(in) :: major_radius_dot, minor_radius_dot, wave_number_dot
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        complex(dp), intent(out) :: single_layer_dot, double_layer_dot(3)
        integer, intent(out) :: status

        real(dp) :: first_jacobian, first_jacobian_dot, first_point(3)
        real(dp) :: first_point_dot(3), second_jacobian, second_jacobian_dot
        real(dp) :: second_point(3), second_point_dot(3), displacement(3)
        real(dp) :: displacement_dot(3), first_tangent_xi(3)
        real(dp) :: first_tangent_eta(3), first_tangent_xi_dot(3)
        real(dp) :: first_tangent_eta_dot(3), second_tangent_xi(3)
        real(dp) :: second_tangent_eta(3), second_tangent_xi_dot(3)
        real(dp) :: second_tangent_eta_dot(3), normal_source(3)
        real(dp) :: normal_source_dot(3), source_cross(3), source_cross_dot(3)
        real(dp) :: radius, radius_dot, weight, weight_dot, dot_normal
        real(dp) :: dot_normal_dot, phase_dot, source_barycentric(3)
        complex(dp) :: exponential, numerator, numerator_dot
        complex(dp) :: green_correction, green_correction_dot
        complex(dp) :: normal_correction, normal_correction_dot
        integer :: first_node, second_node

        single_layer_dot = cmplx(0.0_dp, 0.0_dp, dp)
        double_layer_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (size(xi) == 0 .or. size(eta) /= size(xi) .or. &
            size(weights) /= size(xi)) return
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
                source_cross = cross_product( &
                    second_tangent_xi, second_tangent_eta)
                source_cross_dot = cross_product( &
                    second_tangent_xi_dot, second_tangent_eta) + &
                    cross_product(second_tangent_xi, second_tangent_eta_dot)
                normal_source = source_cross/second_jacobian
                normal_source_dot = (source_cross_dot - &
                    normal_source*second_jacobian_dot)/second_jacobian
                displacement = first_point - second_point
                displacement_dot = first_point_dot - second_point_dot
                radius = norm2(displacement)
                weight = weights(first_node)*weights(second_node)* &
                    first_jacobian*second_jacobian
                weight_dot = weights(first_node)*weights(second_node)* ( &
                    first_jacobian_dot*second_jacobian + &
                    first_jacobian*second_jacobian_dot)
                if (radius <= sqrt(tiny(1.0_dp))) then
                    green_correction = cmplx(0.0_dp, &
                        wave_number/four_pi, dp)
                    green_correction_dot = cmplx(0.0_dp, &
                        wave_number_dot/four_pi, dp)
                    normal_correction = cmplx(0.0_dp, 0.0_dp, dp)
                    normal_correction_dot = cmplx(0.0_dp, 0.0_dp, dp)
                else
                    radius_dot = dot_product(displacement, displacement_dot)/ &
                        radius
                    phase_dot = wave_number_dot*radius + &
                        wave_number*radius_dot
                    exponential = exp(cmplx(0.0_dp, wave_number*radius, dp))
                    numerator = exponential - 1.0_dp
                    green_correction = numerator/(four_pi*radius)
                    green_correction_dot = ( &
                        exponential*cmplx(0.0_dp, phase_dot, dp)*radius - &
                        numerator*radius_dot)/(four_pi*radius**2)
                    dot_normal = dot_product(displacement, normal_source)
                    dot_normal_dot = dot_product(displacement_dot, normal_source) + &
                        dot_product(displacement, normal_source_dot)
                    numerator = exponential*cmplx(1.0_dp, &
                        -wave_number*radius, dp) - 1.0_dp
                    numerator_dot = exponential*cmplx(0.0_dp, phase_dot, dp)* &
                        cmplx(1.0_dp, -wave_number*radius, dp) + &
                        exponential*cmplx(0.0_dp, -phase_dot, dp)
                    normal_correction = numerator*dot_normal/ &
                        (four_pi*radius**3)
                    normal_correction_dot = (numerator_dot*dot_normal + &
                        numerator*dot_normal_dot)/(four_pi*radius**3) - &
                        3.0_dp*numerator*dot_normal*radius_dot/ &
                        (four_pi*radius**4)
                end if
                source_barycentric = [ &
                    1.0_dp - xi(second_node) - eta(second_node), &
                    xi(second_node), eta(second_node)]
                single_layer_dot = single_layer_dot + weight_dot* &
                    green_correction + weight*green_correction_dot
                double_layer_dot = double_layer_dot + (weight_dot* &
                    normal_correction + weight*normal_correction_dot)* &
                    source_barycentric
            end do
        end do
        status = 0
    end subroutine integrate_correction_jvp

    subroutine integrate_correction_vjp( &
            first_parameters, second_parameters, major_radius, minor_radius, &
            wave_number, xi, eta, weights, single_layer_bar, double_layer_bar, &
            first_parameters_bar, second_parameters_bar, major_radius_bar, &
            minor_radius_bar, wave_number_bar, status)
        real(dp), intent(in) :: first_parameters(2, 3), second_parameters(2, 3)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        complex(dp), intent(in) :: single_layer_bar, double_layer_bar(3)
        real(dp), intent(out) :: first_parameters_bar(2, 3)
        real(dp), intent(out) :: second_parameters_bar(2, 3)
        real(dp), intent(out) :: major_radius_bar, minor_radius_bar
        real(dp), intent(out) :: wave_number_bar
        integer, intent(out) :: status

        real(dp) :: first_jacobian, first_point(3), second_jacobian
        real(dp) :: second_point(3), displacement(3), normal_source(3)
        real(dp) :: source_cross(3), source_barycentric(3)
        real(dp) :: first_point_bar(3), second_point_bar(3)
        real(dp) :: first_tangent_xi(3), first_tangent_eta(3)
        real(dp) :: second_tangent_xi(3), second_tangent_eta(3)
        real(dp) :: first_tangent_xi_bar(3), first_tangent_eta_bar(3)
        real(dp) :: second_tangent_xi_bar(3), second_tangent_eta_bar(3)
        real(dp) :: first_jacobian_bar, second_jacobian_bar
        real(dp) :: radius, radius_bar, weight, weight_bar, dot_normal
        real(dp) :: dot_normal_bar, normal_source_bar(3), displacement_bar(3)
        real(dp) :: cross_bar(3)
        real(dp) :: local_major_first, local_minor_first
        real(dp) :: local_major_second, local_minor_second
        real(dp) :: local_first_parameters_bar(2, 3)
        real(dp) :: local_second_parameters_bar(2, 3)
        real(dp) :: xi_bar, eta_bar
        complex(dp) :: green_correction, normal_correction, green_bar
        complex(dp) :: normal_correction_bar, green_radius, green_wave
        complex(dp) :: normal_radius, normal_wave, normal_factor
        complex(dp) :: double_projection_bar
        integer :: first_node, second_node, local_node

        first_parameters_bar = 0.0_dp
        second_parameters_bar = 0.0_dp
        major_radius_bar = 0.0_dp
        minor_radius_bar = 0.0_dp
        wave_number_bar = 0.0_dp
        status = 1
        if (size(xi) == 0 .or. size(eta) /= size(xi) .or. &
            size(weights) /= size(xi)) return
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
                source_cross = cross_product( &
                    second_tangent_xi, second_tangent_eta)
                normal_source = source_cross/second_jacobian
                displacement = first_point - second_point
                radius = norm2(displacement)
                weight = weights(first_node)*weights(second_node)* &
                    first_jacobian*second_jacobian
                source_barycentric = [ &
                    1.0_dp - xi(second_node) - eta(second_node), &
                    xi(second_node), eta(second_node)]
                call correction_kernel_values( &
                    displacement, normal_source, radius, wave_number, &
                    green_correction, normal_correction, green_radius, &
                    green_wave, normal_radius, normal_wave, normal_factor)
                green_bar = weight*single_layer_bar
                double_projection_bar = cmplx(0.0_dp, 0.0_dp, dp)
                do local_node = 1, 3
                    double_projection_bar = double_projection_bar + &
                        double_layer_bar(local_node)* &
                        source_barycentric(local_node)
                end do
                normal_correction_bar = weight*double_projection_bar
                weight_bar = real(conjg(single_layer_bar)*green_correction + &
                    conjg(double_projection_bar)*normal_correction, dp)
                if (radius <= sqrt(tiny(1.0_dp))) then
                    wave_number_bar = wave_number_bar + real(conjg(green_bar)* &
                        cmplx(0.0_dp, 1.0_dp/four_pi, dp), dp)
                    displacement_bar = 0.0_dp
                    normal_source_bar = 0.0_dp
                else
                    radius_bar = real(conjg(green_bar)*green_radius, dp) + &
                        real(conjg(normal_correction_bar*dot_product( &
                        displacement, normal_source))*normal_radius, dp)
                    wave_number_bar = wave_number_bar + &
                        real(conjg(green_bar)*green_wave, dp) + &
                        real(conjg(normal_correction_bar*dot_product( &
                        displacement, normal_source))*normal_wave, dp)
                    dot_normal = dot_product(displacement, normal_source)
                    dot_normal_bar = real(conjg(normal_correction_bar)* &
                        normal_factor, dp)
                    displacement_bar = radius_bar*displacement/radius + &
                        dot_normal_bar*normal_source
                    normal_source_bar = dot_normal_bar*displacement
                end if
                first_point_bar = displacement_bar
                second_point_bar = -displacement_bar
                first_tangent_xi_bar = 0.0_dp
                first_tangent_eta_bar = 0.0_dp
                second_tangent_xi_bar = 0.0_dp
                second_tangent_eta_bar = 0.0_dp
                first_jacobian_bar = weight_bar*weights(first_node)* &
                    weights(second_node)*second_jacobian
                second_jacobian_bar = weight_bar*weights(first_node)* &
                    weights(second_node)*first_jacobian
                cross_bar = normal_source_bar/second_jacobian
                second_jacobian_bar = second_jacobian_bar - &
                    dot_product(normal_source_bar, normal_source)/second_jacobian
                second_tangent_xi_bar = cross_product( &
                    second_tangent_eta, cross_bar)
                second_tangent_eta_bar = cross_product( &
                    cross_bar, second_tangent_xi)
                call evaluate_torus_curved_panel_vjp( &
                    first_parameters, major_radius, minor_radius, xi(first_node), &
                    eta(first_node), first_point_bar, first_tangent_xi_bar, &
                    first_tangent_eta_bar, first_jacobian_bar, &
                    local_first_parameters_bar, local_major_first, local_minor_first, &
                    xi_bar, eta_bar, status)
                if (status /= 0) return
                call evaluate_torus_curved_panel_vjp( &
                    second_parameters, major_radius, minor_radius, &
                    xi(second_node), eta(second_node), second_point_bar, &
                    second_tangent_xi_bar, second_tangent_eta_bar, &
                    second_jacobian_bar, local_second_parameters_bar, &
                    local_major_second, local_minor_second, xi_bar, eta_bar, &
                    status)
                if (status /= 0) return
                first_parameters_bar = first_parameters_bar + &
                    local_first_parameters_bar
                second_parameters_bar = second_parameters_bar + &
                    local_second_parameters_bar
                major_radius_bar = major_radius_bar + local_major_first + &
                    local_major_second
                minor_radius_bar = minor_radius_bar + local_minor_first + &
                    local_minor_second
            end do
        end do
        status = 0
    end subroutine integrate_correction_vjp

    subroutine correction_kernel_values( &
            displacement, normal_source, radius, wave_number, green, normal, &
            green_radius, green_wave, normal_radius, normal_wave, normal_factor)
        real(dp), intent(in) :: displacement(3), normal_source(3), radius
        real(dp), intent(in) :: wave_number
        complex(dp), intent(out) :: green, normal, green_radius, green_wave
        complex(dp), intent(out) :: normal_radius, normal_wave, normal_factor

        complex(dp) :: exponential, numerator, numerator_radius, numerator_wave
        real(dp) :: dot_normal

        if (radius <= sqrt(tiny(1.0_dp))) then
            green = cmplx(0.0_dp, wave_number/four_pi, dp)
            normal = cmplx(0.0_dp, 0.0_dp, dp)
            green_radius = cmplx(0.0_dp, 0.0_dp, dp)
            green_wave = cmplx(0.0_dp, 1.0_dp/four_pi, dp)
            normal_radius = cmplx(0.0_dp, 0.0_dp, dp)
            normal_wave = cmplx(0.0_dp, 0.0_dp, dp)
            normal_factor = cmplx(0.0_dp, 0.0_dp, dp)
            return
        end if
        exponential = exp(cmplx(0.0_dp, wave_number*radius, dp))
        numerator = exponential - 1.0_dp
        green = numerator/(four_pi*radius)
        green_radius = (exponential*cmplx(0.0_dp, wave_number, dp)*radius - &
            numerator)/(four_pi*radius**2)
        green_wave = exponential*cmplx(0.0_dp, radius, dp)/ &
            (four_pi*radius)
        dot_normal = dot_product(displacement, normal_source)
        numerator = exponential*cmplx(1.0_dp, -wave_number*radius, dp) - 1.0_dp
        numerator_radius = exponential*wave_number**2*radius
        numerator_wave = exponential*wave_number*radius**2
        normal_factor = numerator/(four_pi*radius**3)
        normal = normal_factor*dot_normal
        normal_radius = numerator_radius/(four_pi*radius**3) - &
            3.0_dp*numerator/(four_pi*radius**4)
        normal_wave = numerator_wave/(four_pi*radius**3)
    end subroutine correction_kernel_values

    pure logical function valid_inputs( &
            parameters, triangles, major_radius, minor_radius, wave_number, &
            quadrature_degree)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        integer, intent(in) :: quadrature_degree

        valid_inputs = size(parameters, 1) == 2 .and. &
            size(parameters, 2) >= 3 .and. size(triangles, 1) == 3 .and. &
            size(triangles, 2) >= 1 .and. all(triangles >= 1) .and. &
            all(triangles <= size(parameters, 2)) .and. &
            major_radius > minor_radius .and. minor_radius > 0.0_dp .and. &
            wave_number > 0.0_dp .and. quadrature_degree >= 0
    end function valid_inputs

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_helmholtz_torus_curved_bem_ad_3d
