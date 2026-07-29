module fortfem_laplace_representation_3d
    use fortfem_kinds, only: dp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    private

    public :: evaluate_laplace_representation_triangles_3d

contains

    subroutine evaluate_laplace_representation_triangles_3d( &
            vertices, triangles, dirichlet_values, neumann_values, target, &
            quadrature_degree, value, status, gradient)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: dirichlet_values(:), neumann_values(:)
        real(dp), intent(in) :: target(3)
        integer, intent(in) :: quadrature_degree
        real(dp), intent(out) :: value
        integer, intent(out) :: status
        real(dp), intent(out), optional :: gradient(3)

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: barycentric(3), displacement(3), first_edge(3)
        real(dp) :: green, green_normal_derivative, jacobian, normal(3)
        real(dp) :: green_gradient(3), normal_gradient(3)
        real(dp) :: point(3), second_edge(3), trace_value
        integer :: element, node, quadrature_status

        value = 0.0_dp
        if (present(gradient)) gradient = 0.0_dp
        status = 1
        if (size(vertices, 1) /= 3 .or. size(vertices, 2) < 3) return
        if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
        if (any(triangles < 1) .or. any(triangles > size(vertices, 2))) return
        if (size(dirichlet_values) /= size(vertices, 2)) return
        if (size(neumann_values) /= size(triangles, 2)) return
        if (quadrature_degree < 0) return

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
                if (norm2(displacement) <= 0.0_dp) return
                green = 1.0_dp/(4.0_dp*acos(-1.0_dp)*norm2(displacement))
                green_normal_derivative = dot_product(displacement, normal)/ &
                    (4.0_dp*acos(-1.0_dp)*norm2(displacement)**3)
                trace_value = dot_product( &
                    dirichlet_values(triangles(:, element)), barycentric)
                value = value + jacobian*weights(node)*( &
                    trace_value*green_normal_derivative - &
                    green*neumann_values(element))
                if (present(gradient)) then
                    green_gradient = -displacement/( &
                        4.0_dp*acos(-1.0_dp)*norm2(displacement)**3)
                    normal_gradient = normal/( &
                        4.0_dp*acos(-1.0_dp)*norm2(displacement)**3) - &
                        3.0_dp*dot_product(displacement, normal)* &
                        displacement/( &
                        4.0_dp*acos(-1.0_dp)*norm2(displacement)**5)
                    gradient = gradient + jacobian*weights(node)*( &
                        trace_value*normal_gradient - &
                        neumann_values(element)*green_gradient)
                end if
            end do
        end do
        status = 0
    end subroutine evaluate_laplace_representation_triangles_3d

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_laplace_representation_3d
