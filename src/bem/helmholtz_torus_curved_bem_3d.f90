module fortfem_helmholtz_torus_curved_bem_3d
    !! Curved-torus Helmholtz operators by Laplace singularity subtraction.
    !!
    !! The convention G=exp(i*k*r)/(4*pi*r) is outgoing for exp(-i*omega*t).
    !! The singular Laplace part uses radial Duffy quadrature; only the smooth
    !! Helmholtz-minus-Laplace correction is integrated by product quadrature.
    use fortfem_kinds, only: dp
    use fortfem_laplace_torus_curved_bem_3d, only: &
        assemble_laplace_torus_curved_dtn_3d
    use fortfem_torus_curved_panel, only: evaluate_torus_curved_panel
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_linalg, only: dense_solve
    implicit none
    private

    public :: assemble_helmholtz_torus_curved_dtn_3d
    public :: solve_helmholtz_bem_dtn_torus_curved_3d

contains

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

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_helmholtz_torus_curved_bem_3d
