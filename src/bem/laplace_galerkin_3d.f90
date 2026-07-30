module fortfem_laplace_galerkin_3d
    use fortfem_kinds, only: dp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none

    private

    public :: assemble_laplace_single_layer_p0_3d
    public :: assemble_laplace_calderon_p1_p0_3d
    public :: solve_laplace_dirichlet_p0_3d

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dgesv
    end interface

contains

    subroutine assemble_laplace_calderon_p1_p0_3d( &
            vertices, triangles, quadrature_degree, single_layer, &
            double_layer, adjoint_double_layer, hypersingular, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), allocatable, intent(out) :: single_layer(:, :)
        real(dp), allocatable, intent(out) :: double_layer(:, :)
        real(dp), allocatable, intent(out) :: adjoint_double_layer(:, :)
        real(dp), allocatable, intent(out) :: hypersingular(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: barycentric(3), displacement(3)
        real(dp) :: first_gradients(3, 3), first_jacobian
        real(dp) :: first_point_value(3), first_vertices(3, 3)
        real(dp) :: normal(3), radius, second_gradients(3, 3)
        real(dp) :: second_jacobian, second_vertices(3, 3)
        integer :: first, first_local, first_point, quadrature_status
        integer :: second, second_local, second_point

        status = 1
        call assemble_laplace_single_layer_p0_3d( &
            vertices, triangles, quadrature_degree, single_layer, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        allocate(double_layer(size(triangles, 2), size(vertices, 2)))
        allocate(adjoint_double_layer(size(vertices, 2), size(triangles, 2)))
        allocate(hypersingular(size(vertices, 2), size(vertices, 2)))
        double_layer = 0.0_dp
        hypersingular = 0.0_dp

        do second = 1, size(triangles, 2)
            second_vertices = vertices(:, triangles(:, second))
            call triangle_geometry( &
                second_vertices, second_jacobian, normal, second_gradients)
            if (second_jacobian <= 0.0_dp) return
            do first = 1, size(triangles, 2)
                first_vertices = vertices(:, triangles(:, first))
                call triangle_geometry( &
                    first_vertices, first_jacobian, normal, first_gradients)
                if (first_jacobian <= 0.0_dp) return
                do first_local = 1, 3
                    do second_local = 1, 3
                        hypersingular( &
                            triangles(first_local, first), &
                            triangles(second_local, second)) = &
                            hypersingular( &
                            triangles(first_local, first), &
                            triangles(second_local, second)) + &
                            dot_product( &
                            first_gradients(:, first_local), &
                            second_gradients(:, second_local))* &
                            single_layer(first, second)
                    end do
                end do
                if (first == second) cycle
                call triangle_geometry( &
                    second_vertices, second_jacobian, normal, &
                    second_gradients)
                do first_point = 1, size(weights)
                    first_point_value = triangle_point( &
                        first_vertices, xi(first_point), eta(first_point))
                    do second_point = 1, size(weights)
                        displacement = first_point_value - triangle_point( &
                            second_vertices, xi(second_point), &
                            eta(second_point))
                        radius = norm2(displacement)
                        barycentric = [ &
                            1.0_dp - xi(second_point) - eta(second_point), &
                            xi(second_point), eta(second_point)]
                        do second_local = 1, 3
                            double_layer( &
                                first, triangles(second_local, second)) = &
                                double_layer( &
                                first, triangles(second_local, second)) + &
                                first_jacobian*second_jacobian* &
                                weights(first_point)*weights(second_point)* &
                                barycentric(second_local)* &
                                dot_product(displacement, normal)/ &
                                (4.0_dp*acos(-1.0_dp)*radius**3)
                        end do
                    end do
                end do
            end do
        end do
        adjoint_double_layer = transpose(double_layer)
        status = 0
    end subroutine assemble_laplace_calderon_p1_p0_3d

    subroutine solve_laplace_dirichlet_p0_3d( &
            vertices, triangles, boundary_value, quadrature_degree, density, &
            capacity, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: boundary_value
        integer, intent(in) :: quadrature_degree
        real(dp), allocatable, intent(out) :: density(:)
        real(dp), intent(out) :: capacity
        integer, intent(out) :: status

        real(dp), allocatable :: matrix(:, :), right_hand_side(:, :)
        integer, allocatable :: pivots(:)
        real(dp) :: area
        integer :: element, info

        capacity = 0.0_dp
        status = 1
        call assemble_laplace_single_layer_p0_3d( &
            vertices, triangles, quadrature_degree, matrix, status)
        if (status /= 0) return
        allocate( &
            density(size(triangles, 2)), &
            right_hand_side(size(triangles, 2), 1), &
            pivots(size(triangles, 2)))
        do element = 1, size(triangles, 2)
            area = triangle_area( &
                vertices(:, triangles(:, element)))
            right_hand_side(element, 1) = boundary_value*area
        end do
        call dgesv( &
            size(triangles, 2), 1, matrix, size(triangles, 2), pivots, &
            right_hand_side, size(triangles, 2), info)
        if (info /= 0) then
            status = 2
            return
        end if
        density = right_hand_side(:, 1)
        do element = 1, size(triangles, 2)
            capacity = capacity + density(element)*triangle_area( &
                vertices(:, triangles(:, element)))
        end do
        status = 0
    end subroutine solve_laplace_dirichlet_p0_3d

    subroutine assemble_laplace_single_layer_p0_3d( &
            vertices, triangles, quadrature_degree, matrix, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        integer, intent(in) :: quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: first_jacobian, point(3), second_jacobian
        real(dp) :: first_vertices(3, 3), second_vertices(3, 3)
        integer :: first, first_point, quadrature_status, second
        integer :: second_point

        status = 1
        if (.not. valid_surface(vertices, triangles)) return
        if (quadrature_degree < 1) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        allocate(matrix(size(triangles, 2), size(triangles, 2)))
        matrix = 0.0_dp
        do first = 1, size(triangles, 2)
            first_vertices = vertices(:, triangles(:, first))
            first_jacobian = 2.0_dp*triangle_area(first_vertices)
            if (first_jacobian <= 0.0_dp) return
            do second = 1, first
                second_vertices = vertices(:, triangles(:, second))
                second_jacobian = 2.0_dp*triangle_area(second_vertices)
                if (second_jacobian <= 0.0_dp) return
                if (first == second) then
                    do first_point = 1, size(weights)
                        point = triangle_point( &
                            first_vertices, xi(first_point), eta(first_point))
                        matrix(first, second) = matrix(first, second) + &
                            first_jacobian*weights(first_point)* &
                            singular_triangle_potential( &
                            point, second_vertices)/(4.0_dp*acos(-1.0_dp))
                    end do
                else
                    do first_point = 1, size(weights)
                        point = triangle_point( &
                            first_vertices, xi(first_point), eta(first_point))
                        do second_point = 1, size(weights)
                            matrix(first, second) = matrix(first, second) + &
                                first_jacobian*second_jacobian* &
                                weights(first_point)*weights(second_point)/ &
                                (4.0_dp*acos(-1.0_dp)*norm2( &
                                point - triangle_point( &
                                second_vertices, xi(second_point), &
                                eta(second_point))))
                        end do
                    end do
                end if
                matrix(second, first) = matrix(first, second)
            end do
        end do
        status = 0
    end subroutine assemble_laplace_single_layer_p0_3d

    pure function singular_triangle_potential(point, vertices) result(value)
        real(dp), intent(in) :: point(3), vertices(3, 3)
        real(dp) :: value

        real(dp) :: first(3), height, length, numerator, second(3)
        real(dp) :: denominator
        integer :: edge, next

        value = 0.0_dp
        do edge = 1, 3
            next = modulo(edge, 3) + 1
            first = vertices(:, edge) - point
            second = vertices(:, next) - point
            length = norm2(second - first)
            height = norm2(cross_product(first, second))/length
            if (height <= tiny(1.0_dp)) cycle
            numerator = norm2(first) + norm2(second) + length
            denominator = norm2(first) + norm2(second) - length
            value = value + height*log(numerator/denominator)
        end do
    end function singular_triangle_potential

    pure function triangle_point(vertices, xi, eta) result(point)
        real(dp), intent(in) :: vertices(3, 3), xi, eta
        real(dp) :: point(3)

        point = vertices(:, 1) + xi*(vertices(:, 2) - vertices(:, 1)) + &
            eta*(vertices(:, 3) - vertices(:, 1))
    end function triangle_point

    pure function triangle_area(vertices) result(area)
        real(dp), intent(in) :: vertices(3, 3)
        real(dp) :: area

        area = 0.5_dp*norm2(cross_product( &
            vertices(:, 2) - vertices(:, 1), &
            vertices(:, 3) - vertices(:, 1)))
    end function triangle_area

    pure subroutine triangle_geometry( &
            vertices, jacobian, normal, gradients)
        real(dp), intent(in) :: vertices(3, 3)
        real(dp), intent(out) :: jacobian, normal(3), gradients(3, 3)

        normal = cross_product( &
            vertices(:, 2) - vertices(:, 1), &
            vertices(:, 3) - vertices(:, 1))
        jacobian = norm2(normal)
        if (jacobian <= 0.0_dp) then
            gradients = 0.0_dp
            return
        end if
        normal = normal/jacobian
        gradients(:, 1) = &
            cross_product(normal, vertices(:, 3) - vertices(:, 2))/jacobian
        gradients(:, 2) = &
            cross_product(normal, vertices(:, 1) - vertices(:, 3))/jacobian
        gradients(:, 3) = &
            cross_product(normal, vertices(:, 2) - vertices(:, 1))/jacobian
    end subroutine triangle_geometry

    pure function valid_surface(vertices, triangles) result(valid)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        logical :: valid

        valid = size(vertices, 1) == 3 .and. size(vertices, 2) >= 3 .and. &
            size(triangles, 1) == 3 .and. size(triangles, 2) >= 1
        if (.not. valid) return
        valid = all(triangles >= 1) .and. &
            all(triangles <= size(vertices, 2))
    end function valid_surface

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_laplace_galerkin_3d
