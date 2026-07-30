module fortfem_helmholtz_galerkin_3d
    !! P0 Galerkin single-layer BEM for the outgoing three-dimensional
    !! Helmholtz kernel. The singularity subtraction
    !! exp(i*k*r)/r = 1/r + (exp(i*k*r)-1)/r leaves the singular Laplace
    !! contribution to the analytical diagonal treatment.
    !!
    !! The P1/P0 Calderon assembly uses the surface-curl weak regularization
    !! of the hypersingular operator described by Dolz, Harbrecht, and
    !! Multerer, Engineering with Computers 40 (2024),
    !! https://doi.org/10.1007/s00366-024-02013-y. Coincident P1 panel pairs
    !! use a radial Duffy map, so only weakly singular kernels are integrated.
    use fortfem_kinds, only: dp
    use fortfem_laplace_galerkin_3d, only: &
        assemble_laplace_calderon_p1_p0_3d, &
        assemble_laplace_single_layer_p0_3d, &
        assemble_laplace_single_layer_p0_adaptive_3d
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    private

    public :: assemble_helmholtz_single_layer_p0_3d
    public :: assemble_helmholtz_single_layer_p0_adaptive_3d
    public :: assemble_helmholtz_calderon_p1_p0_3d
    public :: assemble_helmholtz_double_layer_p0_3d
    public :: evaluate_helmholtz_cfie_p0_3d
    public :: solve_helmholtz_cfie_p0_3d
    public :: solve_helmholtz_dirichlet_p0_3d

    interface
        subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            complex(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*)
            complex(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine zgesv
    end interface

contains

    subroutine assemble_helmholtz_calderon_p1_p0_3d( &
            vertices, triangles, wave_number, quadrature_degree, single_layer, &
            double_layer, adjoint_double_layer, hypersingular, status)
        real(dp), intent(in) :: vertices(:, :), wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), allocatable, intent(out) :: single_layer(:, :)
        complex(dp), allocatable, intent(out) :: double_layer(:, :)
        complex(dp), allocatable, intent(out) :: adjoint_double_layer(:, :)
        complex(dp), allocatable, intent(out) :: hypersingular(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: laplace_adjoint(:, :)
        real(dp), allocatable :: laplace_double(:, :), laplace_hyper(:, :)
        real(dp), allocatable :: laplace_single(:, :)
        complex(dp), allocatable :: linear_single(:, :)
        real(dp), allocatable :: eta(:), line_nodes(:), line_weights(:)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp) :: barycentric(3), correction, displacement(3)
        real(dp) :: first_curls(3, 3), first_gradients(3, 3)
        real(dp) :: first_jacobian, first_normal(3), first_point_value(3)
        real(dp) :: first_vertices(3, 3), radius
        real(dp) :: second_curls(3, 3), second_gradients(3, 3)
        real(dp) :: second_jacobian, second_normal(3)
        real(dp) :: second_vertices(3, 3)
        complex(dp) :: smooth_derivative
        integer :: first, first_local, first_point, quadrature_status
        integer :: second, second_local, second_point

        status = 1
        if (wave_number < 0.0_dp .or. quadrature_degree < 1) return
        call assemble_laplace_calderon_p1_p0_3d( &
            vertices, triangles, quadrature_degree, laplace_single, &
            laplace_double, laplace_adjoint, laplace_hyper, status)
        if (status /= 0) return
        call assemble_helmholtz_single_layer_p0_3d( &
            vertices, triangles, wave_number, quadrature_degree, &
            single_layer, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        allocate(line_nodes(quadrature_degree), line_weights(quadrature_degree))
        call gauss_legendre_ab( &
            quadrature_degree, 0.0_dp, 1.0_dp, line_nodes, line_weights)
        allocate(double_layer(size(laplace_double, 1), &
            size(laplace_double, 2)))
        allocate(adjoint_double_layer(size(laplace_adjoint, 1), &
            size(laplace_adjoint, 2)))
        allocate(hypersingular(size(laplace_hyper, 1), &
            size(laplace_hyper, 2)))
        double_layer = cmplx(laplace_double, 0.0_dp, dp)
        hypersingular = cmplx(0.0_dp, 0.0_dp, dp)

        do second = 1, size(triangles, 2)
            second_vertices = vertices(:, triangles(:, second))
            call triangle_geometry( &
                second_vertices, second_jacobian, second_normal, &
                second_gradients)
            if (second_jacobian <= 0.0_dp) return
            do first = 1, size(triangles, 2)
                first_vertices = vertices(:, triangles(:, first))
                call triangle_geometry( &
                    first_vertices, first_jacobian, first_normal, &
                    first_gradients)
                if (first_jacobian <= 0.0_dp) return
                do first_local = 1, 3
                    first_curls(:, first_local) = cross_product( &
                        first_normal, first_gradients(:, first_local))
                    do second_local = 1, 3
                        second_curls(:, second_local) = cross_product( &
                            second_normal, second_gradients(:, second_local))
                        hypersingular( &
                            triangles(first_local, first), &
                            triangles(second_local, second)) = &
                            hypersingular( &
                            triangles(first_local, first), &
                            triangles(second_local, second)) + &
                            dot_product( &
                            first_curls(:, first_local), &
                            second_curls(:, second_local))* &
                            single_layer(first, second)
                    end do
                end do
                do first_point = 1, size(weights)
                    first_point_value = triangle_point( &
                        first_vertices, xi(first_point), eta(first_point))
                    do second_point = 1, size(weights)
                        displacement = first_point_value - triangle_point( &
                            second_vertices, xi(second_point), &
                            eta(second_point))
                        radius = norm2(displacement)
                        if (radius <= 64.0_dp*epsilon(1.0_dp)) cycle
                        barycentric = [ &
                            1.0_dp - xi(second_point) - eta(second_point), &
                            xi(second_point), eta(second_point)]
                        smooth_derivative = &
                            helmholtz_double_layer_remainder( &
                            wave_number, radius, &
                            dot_product(displacement, second_normal))
                        do second_local = 1, 3
                            double_layer( &
                                first, triangles(second_local, second)) = &
                                double_layer( &
                                first, triangles(second_local, second)) + &
                                first_jacobian*second_jacobian* &
                                weights(first_point)*weights(second_point)* &
                                barycentric(second_local)*smooth_derivative
                        end do
                    end do
                end do
            end do
        end do
        call assemble_helmholtz_single_layer_p1_3d( &
            vertices, triangles, wave_number, xi, eta, weights, line_nodes, &
            line_weights, linear_single, status)
        if (status /= 0) return
        do second = 1, size(triangles, 2)
            second_vertices = vertices(:, triangles(:, second))
            call triangle_geometry( &
                second_vertices, second_jacobian, second_normal, &
                second_gradients)
            do first = 1, size(triangles, 2)
                first_vertices = vertices(:, triangles(:, first))
                call triangle_geometry( &
                    first_vertices, first_jacobian, first_normal, &
                    first_gradients)
                correction = -wave_number**2* &
                    dot_product(first_normal, second_normal)
                do first_local = 1, 3
                    do second_local = 1, 3
                        hypersingular( &
                            triangles(first_local, first), &
                            triangles(second_local, second)) = &
                            hypersingular( &
                            triangles(first_local, first), &
                            triangles(second_local, second)) + &
                            correction*linear_single( &
                            3*(first - 1) + first_local, &
                            3*(second - 1) + second_local)
                    end do
                end do
            end do
        end do
        adjoint_double_layer = transpose(double_layer)
        status = 0
    end subroutine assemble_helmholtz_calderon_p1_p0_3d

    subroutine assemble_helmholtz_single_layer_p1_3d( &
            vertices, triangles, wave_number, xi, eta, weights, line_nodes, &
            line_weights, matrix, status)
        real(dp), intent(in) :: vertices(:, :), wave_number
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        real(dp), intent(in) :: line_nodes(:), line_weights(:)
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp) :: local_matrix(3, 3)
        real(dp) :: first_vertices(3, 3), second_vertices(3, 3)
        integer :: first, first_local, second, second_local

        allocate(matrix(3*size(triangles, 2), 3*size(triangles, 2)))
        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        do first = 1, size(triangles, 2)
            first_vertices = vertices(:, triangles(:, first))
            do second = 1, first
                second_vertices = vertices(:, triangles(:, second))
                if (first == second) then
                    call integrate_coincident_p1_pair( &
                        first_vertices, wave_number, xi, eta, weights, &
                        line_nodes, line_weights, local_matrix)
                    local_matrix = 0.5_dp* &
                        (local_matrix + transpose(local_matrix))
                else
                    call integrate_regular_p1_pair( &
                        first_vertices, second_vertices, wave_number, xi, eta, &
                        weights, local_matrix)
                end if
                do first_local = 1, 3
                    do second_local = 1, 3
                        matrix(3*(first - 1) + first_local, &
                            3*(second - 1) + second_local) = &
                            local_matrix(first_local, second_local)
                        matrix(3*(second - 1) + second_local, &
                            3*(first - 1) + first_local) = &
                            local_matrix(first_local, second_local)
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_helmholtz_single_layer_p1_3d

    subroutine integrate_regular_p1_pair( &
            first_vertices, second_vertices, wave_number, xi, eta, weights, &
            matrix)
        real(dp), intent(in) :: first_vertices(3, 3), second_vertices(3, 3)
        real(dp), intent(in) :: wave_number, xi(:), eta(:), weights(:)
        complex(dp), intent(out) :: matrix(3, 3)

        real(dp) :: first_barycentric(3), first_jacobian, first_point(3)
        real(dp) :: radius, second_barycentric(3), second_jacobian
        real(dp) :: second_point(3)
        integer :: first_local, first_node, second_local, second_node

        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        first_jacobian = 2.0_dp*triangle_area(first_vertices)
        second_jacobian = 2.0_dp*triangle_area(second_vertices)
        do first_node = 1, size(weights)
            first_barycentric = [ &
                1.0_dp - xi(first_node) - eta(first_node), &
                xi(first_node), eta(first_node)]
            first_point = triangle_point( &
                first_vertices, xi(first_node), eta(first_node))
            do second_node = 1, size(weights)
                second_barycentric = [ &
                    1.0_dp - xi(second_node) - eta(second_node), &
                    xi(second_node), eta(second_node)]
                second_point = triangle_point( &
                    second_vertices, xi(second_node), eta(second_node))
                radius = norm2(first_point - second_point)
                do first_local = 1, 3
                    do second_local = 1, 3
                        matrix(first_local, second_local) = &
                            matrix(first_local, second_local) + &
                            first_jacobian*second_jacobian* &
                            weights(first_node)*weights(second_node)* &
                            first_barycentric(first_local)* &
                            second_barycentric(second_local)* &
                            helmholtz_green(wave_number, radius)
                    end do
                end do
            end do
        end do
    end subroutine integrate_regular_p1_pair

    subroutine integrate_coincident_p1_pair( &
            panel_vertices, wave_number, xi, eta, weights, line_nodes, &
            line_weights, matrix)
        real(dp), intent(in) :: panel_vertices(3, 3), wave_number
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        real(dp), intent(in) :: line_nodes(:), line_weights(:)
        complex(dp), intent(out) :: matrix(3, 3)

        real(dp) :: direction(3), first_barycentric(3), first_point(3)
        real(dp) :: jacobian, radius, rho, second_barycentric(3)
        real(dp) :: second_point(3), t, wedge_first(3), wedge_jacobian
        real(dp) :: wedge_second(3)
        integer :: first_local, first_node, radial_node, second_local
        integer :: tangential_node, wedge

        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        jacobian = 2.0_dp*triangle_area(panel_vertices)
        do first_node = 1, size(weights)
            first_barycentric = [ &
                1.0_dp - xi(first_node) - eta(first_node), &
                xi(first_node), eta(first_node)]
            first_point = triangle_point( &
                panel_vertices, xi(first_node), eta(first_node))
            do wedge = 1, 3
                wedge_first = panel_vertices(:, wedge) - first_point
                wedge_second = panel_vertices(:, modulo(wedge, 3) + 1) - &
                    first_point
                wedge_jacobian = norm2( &
                    cross_product(wedge_first, wedge_second))
                do radial_node = 1, size(line_nodes)
                    rho = line_nodes(radial_node)
                    do tangential_node = 1, size(line_nodes)
                        t = line_nodes(tangential_node)
                        direction = (1.0_dp - t)*wedge_first + t*wedge_second
                        second_point = first_point + rho*direction
                        radius = rho*norm2(direction)
                        call triangle_barycentric( &
                            panel_vertices, second_point, second_barycentric)
                        do first_local = 1, 3
                            do second_local = 1, 3
                                matrix(first_local, second_local) = &
                                    matrix(first_local, second_local) + &
                                    jacobian*weights(first_node)* &
                                    line_weights(radial_node)* &
                                    line_weights(tangential_node)*rho* &
                                    wedge_jacobian* &
                                    first_barycentric(first_local)* &
                                    second_barycentric(second_local)* &
                                    helmholtz_green(wave_number, radius)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine integrate_coincident_p1_pair

    subroutine assemble_helmholtz_single_layer_p0_adaptive_3d( &
            vertices, triangles, wave_number, quadrature_degree, tolerance, &
            max_depth, matrix, status)
        real(dp), intent(in) :: vertices(:, :), wave_number, tolerance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: laplace_matrix(:, :)

        status = 1
        if (wave_number < 0.0_dp) return
        call assemble_laplace_single_layer_p0_adaptive_3d( &
            vertices, triangles, tolerance, max_depth, laplace_matrix, status)
        if (status /= 0) return
        allocate(matrix(size(laplace_matrix, 1), size(laplace_matrix, 2)))
        matrix = cmplx(laplace_matrix, 0.0_dp, dp)
        call add_single_layer_smooth_remainder( &
            vertices, triangles, wave_number, quadrature_degree, matrix, status)
    end subroutine assemble_helmholtz_single_layer_p0_adaptive_3d

    subroutine solve_helmholtz_cfie_p0_3d( &
            vertices, triangles, boundary_value, wave_number, coupling, &
            quadrature_degree, density, status)
        real(dp), intent(in) :: vertices(:, :), wave_number, coupling
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(in) :: boundary_value
        complex(dp), allocatable, intent(out) :: density(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: double_layer(:, :), matrix(:, :)
        complex(dp), allocatable :: right_hand_side(:, :)
        complex(dp), allocatable :: single_layer(:, :)
        integer, allocatable :: pivots(:)
        real(dp) :: area
        integer :: element, info

        status = 1
        if (coupling <= 0.0_dp) return
        call assemble_helmholtz_single_layer_p0_3d( &
            vertices, triangles, wave_number, quadrature_degree, &
            single_layer, status)
        if (status /= 0) return
        call assemble_helmholtz_double_layer_p0_3d( &
            vertices, triangles, wave_number, quadrature_degree, &
            double_layer, status)
        if (status /= 0) return
        matrix = double_layer - cmplx(0.0_dp, coupling, dp)*single_layer
        allocate( &
            density(size(triangles, 2)), &
            right_hand_side(size(triangles, 2), 1), &
            pivots(size(triangles, 2)))
        do element = 1, size(triangles, 2)
            area = triangle_area(vertices(:, triangles(:, element)))
            matrix(element, element) = matrix(element, element) + 0.5_dp*area
            right_hand_side(element, 1) = boundary_value*area
        end do
        call zgesv( &
            size(triangles, 2), 1, matrix, size(triangles, 2), pivots, &
            right_hand_side, size(triangles, 2), info)
        if (info /= 0) then
            status = 2
            return
        end if
        density = right_hand_side(:, 1)
        status = 0
    end subroutine solve_helmholtz_cfie_p0_3d

    subroutine evaluate_helmholtz_cfie_p0_3d( &
            vertices, triangles, density, target, wave_number, coupling, &
            quadrature_degree, value, status)
        real(dp), intent(in) :: vertices(:, :), target(3)
        real(dp), intent(in) :: wave_number, coupling
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(in) :: density(:)
        complex(dp), intent(out) :: value
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: displacement(3), jacobian, normal(3), point(3), radius
        real(dp) :: triangle_vertices(3, 3)
        complex(dp) :: green, normal_derivative
        integer :: element, node, quadrature_status

        value = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (size(density) /= size(triangles, 2)) return
        if (wave_number < 0.0_dp .or. coupling <= 0.0_dp) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        do element = 1, size(triangles, 2)
            triangle_vertices = vertices(:, triangles(:, element))
            normal = cross_product( &
                triangle_vertices(:, 2) - triangle_vertices(:, 1), &
                triangle_vertices(:, 3) - triangle_vertices(:, 1))
            jacobian = norm2(normal)
            if (jacobian <= 0.0_dp) return
            normal = normal/jacobian
            do node = 1, size(weights)
                point = triangle_point( &
                    triangle_vertices, xi(node), eta(node))
                displacement = target - point
                radius = norm2(displacement)
                if (radius <= 0.0_dp) return
                green = exp(cmplx(0.0_dp, wave_number*radius, dp))/ &
                    (4.0_dp*acos(-1.0_dp)*radius)
                normal_derivative = green* &
                    cmplx(1.0_dp, -wave_number*radius, dp)* &
                    dot_product(displacement, normal)/radius**2
                value = value + jacobian*weights(node)*density(element)* &
                    (normal_derivative - cmplx(0.0_dp, coupling, dp)*green)
            end do
        end do
        status = 0
    end subroutine evaluate_helmholtz_cfie_p0_3d

    subroutine assemble_helmholtz_double_layer_p0_3d( &
            vertices, triangles, wave_number, quadrature_degree, matrix, &
            status)
        real(dp), intent(in) :: vertices(:, :), wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: displacement(3), first_jacobian, first_point_value(3)
        real(dp) :: normal(3), radius, second_jacobian
        real(dp) :: first_vertices(3, 3), second_vertices(3, 3)
        complex(dp) :: green
        integer :: first, first_point, quadrature_status, second, second_point

        status = 1
        if (wave_number < 0.0_dp) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        allocate(matrix(size(triangles, 2), size(triangles, 2)))
        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        do second = 1, size(triangles, 2)
            second_vertices = vertices(:, triangles(:, second))
            normal = cross_product( &
                second_vertices(:, 2) - second_vertices(:, 1), &
                second_vertices(:, 3) - second_vertices(:, 1))
            second_jacobian = norm2(normal)
            if (second_jacobian <= 0.0_dp) return
            normal = normal/second_jacobian
            do first = 1, size(triangles, 2)
                if (first == second) cycle
                first_vertices = vertices(:, triangles(:, first))
                first_jacobian = 2.0_dp*triangle_area(first_vertices)
                do first_point = 1, size(weights)
                    first_point_value = triangle_point( &
                        first_vertices, xi(first_point), eta(first_point))
                    do second_point = 1, size(weights)
                        displacement = first_point_value - triangle_point( &
                            second_vertices, xi(second_point), &
                            eta(second_point))
                        radius = norm2(displacement)
                        green = exp(cmplx( &
                            0.0_dp, wave_number*radius, dp))/ &
                            (4.0_dp*acos(-1.0_dp)*radius)
                        matrix(first, second) = matrix(first, second) + &
                            first_jacobian*second_jacobian* &
                            weights(first_point)*weights(second_point)*green* &
                            cmplx(1.0_dp, -wave_number*radius, dp)* &
                            dot_product(displacement, normal)/radius**2
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_helmholtz_double_layer_p0_3d

    subroutine solve_helmholtz_dirichlet_p0_3d( &
            vertices, triangles, boundary_value, wave_number, &
            quadrature_degree, density, status)
        real(dp), intent(in) :: vertices(:, :), wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(in) :: boundary_value
        complex(dp), allocatable, intent(out) :: density(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: matrix(:, :), right_hand_side(:, :)
        integer, allocatable :: pivots(:)
        integer :: element, info

        status = 1
        call assemble_helmholtz_single_layer_p0_3d( &
            vertices, triangles, wave_number, quadrature_degree, matrix, &
            status)
        if (status /= 0) return
        allocate( &
            density(size(triangles, 2)), &
            right_hand_side(size(triangles, 2), 1), &
            pivots(size(triangles, 2)))
        do element = 1, size(triangles, 2)
            right_hand_side(element, 1) = boundary_value*triangle_area( &
                vertices(:, triangles(:, element)))
        end do
        call zgesv( &
            size(triangles, 2), 1, matrix, size(triangles, 2), pivots, &
            right_hand_side, size(triangles, 2), info)
        if (info /= 0) then
            status = 2
            return
        end if
        density = right_hand_side(:, 1)
        status = 0
    end subroutine solve_helmholtz_dirichlet_p0_3d

    subroutine assemble_helmholtz_single_layer_p0_3d( &
            vertices, triangles, wave_number, quadrature_degree, matrix, &
            status)
        real(dp), intent(in) :: vertices(:, :), wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: laplace_matrix(:, :)
        status = 1
        if (wave_number < 0.0_dp) return
        call assemble_laplace_single_layer_p0_3d( &
            vertices, triangles, quadrature_degree, laplace_matrix, status)
        if (status /= 0) return
        allocate(matrix(size(laplace_matrix, 1), size(laplace_matrix, 2)))
        matrix = cmplx(laplace_matrix, 0.0_dp, dp)
        call add_single_layer_smooth_remainder( &
            vertices, triangles, wave_number, quadrature_degree, matrix, status)
    end subroutine assemble_helmholtz_single_layer_p0_3d

    subroutine add_single_layer_smooth_remainder( &
            vertices, triangles, wave_number, quadrature_degree, matrix, status)
        real(dp), intent(in) :: vertices(:, :), wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(inout) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: first_jacobian, first_point_value(3), radius
        real(dp) :: second_jacobian
        real(dp) :: first_vertices(3, 3), second_vertices(3, 3)
        integer :: first, first_point, quadrature_status, second
        integer :: second_point

        status = 1
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        do first = 1, size(triangles, 2)
            first_vertices = vertices(:, triangles(:, first))
            first_jacobian = 2.0_dp*triangle_area(first_vertices)
            do second = 1, first
                second_vertices = vertices(:, triangles(:, second))
                second_jacobian = 2.0_dp*triangle_area(second_vertices)
                do first_point = 1, size(weights)
                    first_point_value = triangle_point( &
                        first_vertices, xi(first_point), eta(first_point))
                    do second_point = 1, size(weights)
                        radius = norm2(first_point_value - triangle_point( &
                            second_vertices, xi(second_point), &
                            eta(second_point)))
                        matrix(first, second) = matrix(first, second) + &
                            first_jacobian*second_jacobian* &
                            weights(first_point)*weights(second_point)* &
                            smooth_remainder(wave_number, radius)
                    end do
                end do
                matrix(second, first) = matrix(first, second)
            end do
        end do
        status = 0
    end subroutine add_single_layer_smooth_remainder

    pure function smooth_remainder(wave_number, radius) result(value)
        real(dp), intent(in) :: wave_number, radius
        complex(dp) :: value

        if (radius <= 64.0_dp*epsilon(1.0_dp)) then
            value = cmplx(0.0_dp, wave_number, dp)/(4.0_dp*acos(-1.0_dp))
        else
            value = (exp(cmplx(0.0_dp, wave_number*radius, dp)) - 1.0_dp)/ &
                (4.0_dp*acos(-1.0_dp)*radius)
        end if
    end function smooth_remainder

    pure function helmholtz_double_layer_remainder( &
            wave_number, radius, normal_displacement) result(value)
        real(dp), intent(in) :: wave_number, radius, normal_displacement
        complex(dp) :: value

        complex(dp) :: phase

        if (radius <= 64.0_dp*epsilon(1.0_dp)) then
            value = cmplx(0.0_dp, 0.0_dp, dp)
            return
        end if
        phase = exp(cmplx(0.0_dp, wave_number*radius, dp))
        value = ( &
            phase*cmplx(1.0_dp, -wave_number*radius, dp) - 1.0_dp)* &
            normal_displacement/(4.0_dp*acos(-1.0_dp)*radius**3)
    end function helmholtz_double_layer_remainder

    pure function helmholtz_green(wave_number, radius) result(value)
        real(dp), intent(in) :: wave_number, radius
        complex(dp) :: value

        value = exp(cmplx(0.0_dp, wave_number*radius, dp))/ &
            (4.0_dp*acos(-1.0_dp)*radius)
    end function helmholtz_green

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

    pure subroutine triangle_barycentric(vertices, point, barycentric)
        real(dp), intent(in) :: vertices(3, 3), point(3)
        real(dp), intent(out) :: barycentric(3)

        real(dp) :: d00, d01, d11, d20, d21, denominator
        real(dp) :: first(3), offset(3), second(3)

        first = vertices(:, 2) - vertices(:, 1)
        second = vertices(:, 3) - vertices(:, 1)
        offset = point - vertices(:, 1)
        d00 = dot_product(first, first)
        d01 = dot_product(first, second)
        d11 = dot_product(second, second)
        d20 = dot_product(offset, first)
        d21 = dot_product(offset, second)
        denominator = d00*d11 - d01*d01
        barycentric(2) = (d11*d20 - d01*d21)/denominator
        barycentric(3) = (d00*d21 - d01*d20)/denominator
        barycentric(1) = 1.0_dp - barycentric(2) - barycentric(3)
    end subroutine triangle_barycentric

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

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_helmholtz_galerkin_3d
