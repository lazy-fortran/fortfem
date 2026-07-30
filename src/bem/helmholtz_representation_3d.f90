module fortfem_helmholtz_representation_3d
    !! Outgoing three-dimensional Helmholtz Green representation on triangles.
    !!
    !! The convention G=exp(i*k*r)/(4*pi*r) corresponds to outgoing waves for
    !! exp(-i*omega*t) time dependence. See Colton and Kress, Integral Equation
    !! Methods in Scattering Theory, DOI:10.1137/1.9781611973167.
    use fortfem_kinds, only: dp
    use fortfem_torus_curved_panel, only: evaluate_torus_curved_panel
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    private

    public :: evaluate_helmholtz_representation_triangles_3d
    public :: evaluate_helmholtz_representation_torus_curved_3d

contains

    subroutine evaluate_helmholtz_representation_torus_curved_3d( &
            parameters, triangles, major_radius, minor_radius, &
            dirichlet_values, neumann_values, target, wave_number, &
            quadrature_degree, value, status)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, target(3)
        real(dp), intent(in) :: wave_number
        complex(dp), intent(in) :: dirichlet_values(:), neumann_values(:)
        integer, intent(in) :: quadrature_degree
        complex(dp), intent(out) :: value
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: barycentric(3), displacement(3), jacobian, normal(3)
        real(dp) :: point(3), radius, tangent_eta(3), tangent_xi(3)
        complex(dp) :: green, green_normal_derivative, neumann_value
        complex(dp) :: trace_value
        integer :: element, node, quadrature_status

        value = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (size(parameters, 1) /= 2 .or. size(parameters, 2) < 3) return
        if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
        if (any(triangles < 1)) return
        if (any(triangles > size(parameters, 2))) return
        if (size(dirichlet_values) /= size(parameters, 2)) return
        if (size(neumann_values) /= size(triangles, 2) .and. &
            size(neumann_values) /= size(parameters, 2)) return
        if (wave_number < 0.0_dp .or. quadrature_degree < 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
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
                radius = norm2(displacement)
                if (radius <= tiny(1.0_dp)) return
                green = exp(cmplx(0.0_dp, wave_number*radius, dp))/ &
                    (4.0_dp*acos(-1.0_dp)*radius)
                green_normal_derivative = green* &
                    cmplx(1.0_dp, -wave_number*radius, dp)* &
                    dot_product(displacement, normal)/radius**2
                trace_value = dot_product( &
                    dirichlet_values(triangles(:, element)), barycentric)
                if (size(neumann_values) == size(parameters, 2)) then
                    neumann_value = dot_product( &
                        neumann_values(triangles(:, element)), barycentric)
                else
                    neumann_value = neumann_values(element)
                end if
                value = value + jacobian*weights(node)*( &
                    trace_value*green_normal_derivative - &
                    green*neumann_value)
            end do
        end do
        status = 0
    end subroutine evaluate_helmholtz_representation_torus_curved_3d

    subroutine evaluate_helmholtz_representation_triangles_3d( &
            vertices, triangles, dirichlet_values, neumann_values, target, &
            wave_number, quadrature_degree, value, status)
        real(dp), intent(in) :: vertices(:, :), target(3), wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(in) :: dirichlet_values(:), neumann_values(:)
        complex(dp), intent(out) :: value
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: barycentric(3), displacement(3), first_edge(3)
        real(dp) :: jacobian, normal(3), point(3), radius, second_edge(3)
        complex(dp) :: green, green_normal_derivative, trace_value
        integer :: element, node, quadrature_status

        value = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (size(vertices, 1) /= 3 .or. size(vertices, 2) < 3) return
        if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
        if (any(triangles < 1) .or. any(triangles > size(vertices, 2))) return
        if (size(dirichlet_values) /= size(vertices, 2)) return
        if (size(neumann_values) /= size(triangles, 2)) return
        if (wave_number < 0.0_dp .or. quadrature_degree < 0) return

        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        do element = 1, size(triangles, 2)
            first_edge = vertices(:, triangles(2, element)) - &
                vertices(:, triangles(1, element))
            second_edge = vertices(:, triangles(3, element)) - &
                vertices(:, triangles(1, element))
            normal = cross_product(first_edge, second_edge)
            jacobian = norm2(normal)
            if (jacobian <= 0.0_dp) return
            normal = normal/jacobian
            do node = 1, size(weights)
                barycentric = [1.0_dp - xi(node) - eta(node), &
                    xi(node), eta(node)]
                point = matmul( &
                    vertices(:, triangles(:, element)), barycentric)
                displacement = target - point
                radius = norm2(displacement)
                if (radius <= 0.0_dp) return
                green = exp(cmplx(0.0_dp, wave_number*radius, dp))/ &
                    (4.0_dp*acos(-1.0_dp)*radius)
                green_normal_derivative = green* &
                    cmplx(1.0_dp, -wave_number*radius, dp)* &
                    dot_product(displacement, normal)/radius**2
                trace_value = dot_product( &
                    dirichlet_values(triangles(:, element)), barycentric)
                value = value + jacobian*weights(node)*( &
                    trace_value*green_normal_derivative - &
                    green*neumann_values(element))
            end do
        end do
        status = 0
    end subroutine evaluate_helmholtz_representation_triangles_3d

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_helmholtz_representation_3d
