module fortfem_laplace_panel_pair_3d
    use fortfem_generated_laplace_singular_edge_potential, only: &
        generated_laplace_singular_edge_potential
    use fortfem_generated_laplace_singular_edge_potential_jvp, only: &
        generated_laplace_singular_edge_potential_jvp
    use fortfem_generated_laplace_singular_edge_potential_vjp, only: &
        generated_laplace_singular_edge_potential_vjp
    use fortfem_generated_laplace_single_layer_integrand, only: &
        generated_laplace_single_layer_integrand
    use fortfem_generated_laplace_single_layer_integrand_jvp, only: &
        generated_laplace_single_layer_integrand_jvp
    use fortfem_generated_laplace_single_layer_integrand_vjp, only: &
        generated_laplace_single_layer_integrand_vjp
    use fortfem_kinds, only: dp
    use fortfem_surface_triangle_geometry_3d, only: &
        evaluate_surface_triangle_geometry_3d, &
        evaluate_surface_triangle_geometry_3d_jvp, &
        evaluate_surface_triangle_geometry_3d_vjp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none

    private

    public :: integrate_laplace_single_layer_regular_panel_pair_p0_3d
    public :: integrate_laplace_single_layer_regular_panel_pair_p0_3d_jvp
    public :: integrate_laplace_single_layer_regular_panel_pair_p0_3d_vjp
    public :: assemble_laplace_single_layer_p0_3d_jvp
    public :: assemble_laplace_single_layer_p0_3d_vjp

contains

    subroutine assemble_laplace_single_layer_p0_3d_jvp( &
            vertices, triangles, quadrature_degree, vertices_dot, matrix_dot, &
            status)
        real(dp), intent(in) :: vertices(:, :), vertices_dot(:, :)
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), allocatable, intent(out) :: matrix_dot(:, :)
        integer, intent(out) :: status

        real(dp) :: first_vertices(3, 3), first_vertices_dot(3, 3)
        real(dp) :: second_vertices(3, 3), second_vertices_dot(3, 3)
        real(dp) :: value_dot
        integer :: first, second

        status = 1
        if (.not. valid_surface_direction(vertices, triangles, vertices_dot)) &
            return
        allocate(matrix_dot(size(triangles, 2), size(triangles, 2)))
        matrix_dot = 0.0_dp
        do first = 1, size(triangles, 2)
            first_vertices = vertices(:, triangles(:, first))
            first_vertices_dot = vertices_dot(:, triangles(:, first))
            call integrate_laplace_single_layer_self_panel_p0_3d_jvp( &
                first_vertices, quadrature_degree, first_vertices_dot, &
                value_dot, status)
            if (status /= 0) return
            matrix_dot(first, first) = value_dot
            do second = 1, first - 1
                second_vertices = vertices(:, triangles(:, second))
                second_vertices_dot = vertices_dot(:, triangles(:, second))
                call &
                    integrate_laplace_single_layer_regular_panel_pair_p0_3d_jvp( &
                    first_vertices, second_vertices, quadrature_degree, &
                    first_vertices_dot, second_vertices_dot, value_dot, status)
                if (status /= 0) return
                matrix_dot(first, second) = value_dot
                matrix_dot(second, first) = value_dot
            end do
        end do
        status = 0
    end subroutine assemble_laplace_single_layer_p0_3d_jvp

    subroutine assemble_laplace_single_layer_p0_3d_vjp( &
            vertices, triangles, quadrature_degree, matrix_bar, vertices_bar, &
            status)
        real(dp), intent(in) :: vertices(:, :), matrix_bar(:, :)
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), intent(out) :: vertices_bar(:, :)
        integer, intent(out) :: status

        real(dp) :: first_vertices(3, 3), first_vertices_bar(3, 3)
        real(dp) :: second_vertices(3, 3), second_vertices_bar(3, 3)
        real(dp) :: value_bar
        integer :: first, local, second

        status = 1
        if (.not. valid_surface_direction(vertices, triangles, vertices)) &
            return
        if (any(shape(vertices_bar) /= shape(vertices))) return
        if (size(matrix_bar, 1) /= size(triangles, 2) .or. &
            size(matrix_bar, 2) /= size(triangles, 2)) return
        vertices_bar = 0.0_dp
        do first = 1, size(triangles, 2)
            first_vertices = vertices(:, triangles(:, first))
            call integrate_laplace_single_layer_self_panel_p0_3d_vjp( &
                first_vertices, quadrature_degree, matrix_bar(first, first), &
                first_vertices_bar, status)
            if (status /= 0) return
            do local = 1, 3
                vertices_bar(:, triangles(local, first)) = &
                    vertices_bar(:, triangles(local, first)) + &
                    first_vertices_bar(:, local)
            end do
            do second = 1, first - 1
                second_vertices = vertices(:, triangles(:, second))
                value_bar = matrix_bar(first, second) + &
                    matrix_bar(second, first)
                call &
                    integrate_laplace_single_layer_regular_panel_pair_p0_3d_vjp( &
                    first_vertices, second_vertices, quadrature_degree, &
                    value_bar, first_vertices_bar, second_vertices_bar, status)
                if (status /= 0) return
                do local = 1, 3
                    vertices_bar(:, triangles(local, first)) = &
                        vertices_bar(:, triangles(local, first)) + &
                        first_vertices_bar(:, local)
                    vertices_bar(:, triangles(local, second)) = &
                        vertices_bar(:, triangles(local, second)) + &
                        second_vertices_bar(:, local)
                end do
            end do
        end do
        status = 0
    end subroutine assemble_laplace_single_layer_p0_3d_vjp

    subroutine integrate_laplace_single_layer_self_panel_p0_3d_jvp( &
            vertices, quadrature_degree, vertices_dot, value_dot, status)
        real(dp), intent(in) :: vertices(3, 3), vertices_dot(3, 3)
        integer, intent(in) :: quadrature_degree
        real(dp), intent(out) :: value_dot
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: jacobian, jacobian_dot, normal(3), normal_dot(3)
        real(dp) :: point(3), point_dot(3), potential, potential_dot
        real(dp) :: edge_value, edge_value_dot
        integer :: edge, next, node, quadrature_status

        value_dot = 0.0_dp
        status = 1
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        do node = 1, size(weights)
            call evaluate_surface_triangle_geometry_3d( &
                vertices, xi(node), eta(node), point, jacobian, normal, status)
            if (status /= 0) return
            call evaluate_surface_triangle_geometry_3d_jvp( &
                vertices, xi(node), eta(node), vertices_dot, 0.0_dp, 0.0_dp, &
                point_dot, jacobian_dot, normal_dot, status)
            if (status /= 0) return
            potential = 0.0_dp
            potential_dot = 0.0_dp
            do edge = 1, 3
                next = modulo(edge, 3) + 1
                call generated_laplace_singular_edge_potential( &
                    point(1), point(2), point(3), vertices(1, edge), &
                    vertices(2, edge), vertices(3, edge), vertices(1, next), &
                    vertices(2, next), vertices(3, next), edge_value)
                call generated_laplace_singular_edge_potential_jvp( &
                    point(1), point(2), point(3), vertices(1, edge), &
                    vertices(2, edge), vertices(3, edge), vertices(1, next), &
                    vertices(2, next), vertices(3, next), point_dot(1), &
                    point_dot(2), point_dot(3), vertices_dot(1, edge), &
                    vertices_dot(2, edge), vertices_dot(3, edge), &
                    vertices_dot(1, next), vertices_dot(2, next), &
                    vertices_dot(3, next), edge_value_dot)
                potential = potential + edge_value
                potential_dot = potential_dot + edge_value_dot
            end do
            value_dot = value_dot + weights(node)*( &
                jacobian_dot*potential + jacobian*potential_dot)/ &
                (4.0_dp*acos(-1.0_dp))
        end do
        status = 0
    end subroutine integrate_laplace_single_layer_self_panel_p0_3d_jvp

    subroutine integrate_laplace_single_layer_self_panel_p0_3d_vjp( &
            vertices, quadrature_degree, value_bar, vertices_bar, status)
        real(dp), intent(in) :: vertices(3, 3), value_bar
        integer, intent(in) :: quadrature_degree
        real(dp), intent(out) :: vertices_bar(3, 3)
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: jacobian, jacobian_bar, normal(3), point(3)
        real(dp) :: point_bar(3), edge_point_bar(3)
        real(dp) :: first_bar(3), second_bar(3), edge_value, edge_value_bar
        real(dp) :: geometry_vertices_bar(3, 3)
        real(dp) :: local_vertices_bar(3, 3), potential, scale
        real(dp) :: xi_bar, eta_bar
        integer :: edge, next, node, quadrature_status

        vertices_bar = 0.0_dp
        status = 1
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        do node = 1, size(weights)
            call evaluate_surface_triangle_geometry_3d( &
                vertices, xi(node), eta(node), point, jacobian, normal, status)
            if (status /= 0) return
            scale = value_bar*weights(node)/(4.0_dp*acos(-1.0_dp))
            potential = 0.0_dp
            point_bar = 0.0_dp
            local_vertices_bar = 0.0_dp
            do edge = 1, 3
                next = modulo(edge, 3) + 1
                call generated_laplace_singular_edge_potential( &
                    point(1), point(2), point(3), vertices(1, edge), &
                    vertices(2, edge), vertices(3, edge), vertices(1, next), &
                    vertices(2, next), vertices(3, next), edge_value)
                edge_value_bar = scale*jacobian
                call generated_laplace_singular_edge_potential_vjp( &
                    point(1), point(2), point(3), vertices(1, edge), &
                    vertices(2, edge), vertices(3, edge), vertices(1, next), &
                    vertices(2, next), vertices(3, next), edge_value_bar, &
                    edge_point_bar(1), edge_point_bar(2), edge_point_bar(3), &
                    first_bar(1), first_bar(2), first_bar(3), second_bar(1), &
                    second_bar(2), second_bar(3))
                potential = potential + edge_value
                point_bar = point_bar + edge_point_bar
                local_vertices_bar(:, edge) = &
                    local_vertices_bar(:, edge) + first_bar
                local_vertices_bar(:, next) = &
                    local_vertices_bar(:, next) + second_bar
            end do
            jacobian_bar = scale*potential
            call evaluate_surface_triangle_geometry_3d_vjp( &
                vertices, xi(node), eta(node), point_bar, jacobian_bar, &
                [0.0_dp, 0.0_dp, 0.0_dp], geometry_vertices_bar, xi_bar, &
                eta_bar, status)
            if (status /= 0) return
            vertices_bar = vertices_bar + local_vertices_bar + &
                geometry_vertices_bar
        end do
        status = 0
    end subroutine integrate_laplace_single_layer_self_panel_p0_3d_vjp

    subroutine integrate_laplace_single_layer_regular_panel_pair_p0_3d( &
            first_vertices, second_vertices, quadrature_degree, value, status)
        real(dp), intent(in) :: first_vertices(3, 3), second_vertices(3, 3)
        integer, intent(in) :: quadrature_degree
        real(dp), intent(out) :: value
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: first_jacobian, first_point(3), kernel_scale
        real(dp) :: normal(3), second_jacobian, second_point(3), term
        integer :: first_node, quadrature_status, second_node

        value = 0.0_dp
        status = 1
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        do first_node = 1, size(weights)
            call evaluate_surface_triangle_geometry_3d( &
                first_vertices, xi(first_node), eta(first_node), first_point, &
                first_jacobian, normal, status)
            if (status /= 0) return
            do second_node = 1, size(weights)
                call evaluate_surface_triangle_geometry_3d( &
                    second_vertices, xi(second_node), eta(second_node), &
                    second_point, second_jacobian, normal, status)
                if (status /= 0) return
                if (norm2(first_point - second_point) <= &
                    64.0_dp*epsilon(1.0_dp)) return
                kernel_scale = weights(first_node)*weights(second_node)/ &
                    (4.0_dp*acos(-1.0_dp))
                call generated_laplace_single_layer_integrand( &
                    first_point(1), first_point(2), first_point(3), &
                    second_point(1), second_point(2), second_point(3), &
                    first_jacobian, second_jacobian, kernel_scale, term)
                value = value + term
            end do
        end do
        status = 0
    end subroutine integrate_laplace_single_layer_regular_panel_pair_p0_3d

    subroutine integrate_laplace_single_layer_regular_panel_pair_p0_3d_jvp( &
            first_vertices, second_vertices, quadrature_degree, &
            first_vertices_dot, second_vertices_dot, value_dot, status)
        real(dp), intent(in) :: first_vertices(3, 3), second_vertices(3, 3)
        integer, intent(in) :: quadrature_degree
        real(dp), intent(in) :: first_vertices_dot(3, 3)
        real(dp), intent(in) :: second_vertices_dot(3, 3)
        real(dp), intent(out) :: value_dot
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: first_jacobian, first_jacobian_dot, first_point(3)
        real(dp) :: first_point_dot(3), kernel_scale, normal(3), normal_dot(3)
        real(dp) :: second_jacobian, second_jacobian_dot, second_point(3)
        real(dp) :: second_point_dot(3), term_dot
        integer :: first_node, quadrature_status, second_node

        value_dot = 0.0_dp
        status = 1
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        do first_node = 1, size(weights)
            call evaluate_surface_triangle_geometry_3d( &
                first_vertices, xi(first_node), eta(first_node), first_point, &
                first_jacobian, normal, status)
            if (status /= 0) return
            call evaluate_surface_triangle_geometry_3d_jvp( &
                first_vertices, xi(first_node), eta(first_node), &
                first_vertices_dot, 0.0_dp, 0.0_dp, first_point_dot, &
                first_jacobian_dot, normal_dot, status)
            if (status /= 0) return
            do second_node = 1, size(weights)
                call evaluate_surface_triangle_geometry_3d( &
                    second_vertices, xi(second_node), eta(second_node), &
                    second_point, second_jacobian, normal, status)
                if (status /= 0) return
                call evaluate_surface_triangle_geometry_3d_jvp( &
                    second_vertices, xi(second_node), eta(second_node), &
                    second_vertices_dot, 0.0_dp, 0.0_dp, second_point_dot, &
                    second_jacobian_dot, normal_dot, status)
                if (status /= 0) return
                if (norm2(first_point - second_point) <= &
                    64.0_dp*epsilon(1.0_dp)) return
                kernel_scale = weights(first_node)*weights(second_node)/ &
                    (4.0_dp*acos(-1.0_dp))
                call generated_laplace_single_layer_integrand_jvp( &
                    first_point(1), first_point(2), first_point(3), &
                    second_point(1), second_point(2), second_point(3), &
                    first_jacobian, second_jacobian, kernel_scale, &
                    first_point_dot(1), first_point_dot(2), &
                    first_point_dot(3), second_point_dot(1), &
                    second_point_dot(2), second_point_dot(3), &
                    first_jacobian_dot, second_jacobian_dot, term_dot)
                value_dot = value_dot + term_dot
            end do
        end do
        status = 0
    end subroutine integrate_laplace_single_layer_regular_panel_pair_p0_3d_jvp

    subroutine integrate_laplace_single_layer_regular_panel_pair_p0_3d_vjp( &
            first_vertices, second_vertices, quadrature_degree, value_bar, &
            first_vertices_bar, second_vertices_bar, status)
        real(dp), intent(in) :: first_vertices(3, 3), second_vertices(3, 3)
        integer, intent(in) :: quadrature_degree
        real(dp), intent(in) :: value_bar
        real(dp), intent(out) :: first_vertices_bar(3, 3)
        real(dp), intent(out) :: second_vertices_bar(3, 3)
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), first_jacobian_bar(:)
        real(dp), allocatable :: first_point_bar(:, :)
        real(dp), allocatable :: second_jacobian_bar(:)
        real(dp), allocatable :: second_point_bar(:, :), weights(:), xi(:)
        real(dp) :: first_jacobian, first_point(3), kernel_scale
        real(dp) :: normal(3), second_jacobian, second_point(3)
        real(dp) :: local_first_jacobian_bar, local_first_point_bar(3)
        real(dp) :: local_second_jacobian_bar, local_second_point_bar(3)
        real(dp) :: local_vertices_bar(3, 3), xi_bar, eta_bar
        integer :: first_node, quadrature_status, second_node

        first_vertices_bar = 0.0_dp
        second_vertices_bar = 0.0_dp
        status = 1
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        allocate(first_point_bar(3, size(weights)), source=0.0_dp)
        allocate(second_point_bar(3, size(weights)), source=0.0_dp)
        allocate(first_jacobian_bar(size(weights)), source=0.0_dp)
        allocate(second_jacobian_bar(size(weights)), source=0.0_dp)
        do first_node = 1, size(weights)
            call evaluate_surface_triangle_geometry_3d( &
                first_vertices, xi(first_node), eta(first_node), first_point, &
                first_jacobian, normal, status)
            if (status /= 0) return
            do second_node = 1, size(weights)
                call evaluate_surface_triangle_geometry_3d( &
                    second_vertices, xi(second_node), eta(second_node), &
                    second_point, second_jacobian, normal, status)
                if (status /= 0) return
                if (norm2(first_point - second_point) <= &
                    64.0_dp*epsilon(1.0_dp)) return
                kernel_scale = weights(first_node)*weights(second_node)/ &
                    (4.0_dp*acos(-1.0_dp))
                call generated_laplace_single_layer_integrand_vjp( &
                    first_point(1), first_point(2), first_point(3), &
                    second_point(1), second_point(2), second_point(3), &
                    first_jacobian, second_jacobian, kernel_scale, value_bar, &
                    local_first_point_bar(1), local_first_point_bar(2), &
                    local_first_point_bar(3), local_second_point_bar(1), &
                    local_second_point_bar(2), local_second_point_bar(3), &
                    local_first_jacobian_bar, local_second_jacobian_bar)
                first_point_bar(:, first_node) = &
                    first_point_bar(:, first_node) + local_first_point_bar
                second_point_bar(:, second_node) = &
                    second_point_bar(:, second_node) + local_second_point_bar
                first_jacobian_bar(first_node) = &
                    first_jacobian_bar(first_node) + local_first_jacobian_bar
                second_jacobian_bar(second_node) = &
                    second_jacobian_bar(second_node) + local_second_jacobian_bar
            end do
        end do
        do first_node = 1, size(weights)
            call evaluate_surface_triangle_geometry_3d_vjp( &
                first_vertices, xi(first_node), eta(first_node), &
                first_point_bar(:, first_node), &
                first_jacobian_bar(first_node), [0.0_dp, 0.0_dp, 0.0_dp], &
                local_vertices_bar, xi_bar, eta_bar, status)
            if (status /= 0) return
            first_vertices_bar = first_vertices_bar + local_vertices_bar
            call evaluate_surface_triangle_geometry_3d_vjp( &
                second_vertices, xi(first_node), eta(first_node), &
                second_point_bar(:, first_node), &
                second_jacobian_bar(first_node), [0.0_dp, 0.0_dp, 0.0_dp], &
                local_vertices_bar, xi_bar, eta_bar, status)
            if (status /= 0) return
            second_vertices_bar = second_vertices_bar + local_vertices_bar
        end do
        status = 0
    end subroutine integrate_laplace_single_layer_regular_panel_pair_p0_3d_vjp

    pure function valid_surface_direction( &
            vertices, triangles, vertices_dot) result(valid)
        real(dp), intent(in) :: vertices(:, :), vertices_dot(:, :)
        integer, intent(in) :: triangles(:, :)
        logical :: valid

        valid = size(vertices, 1) == 3 .and. size(vertices, 2) >= 3 .and. &
            all(shape(vertices_dot) == shape(vertices)) .and. &
            size(triangles, 1) == 3 .and. size(triangles, 2) >= 1
        if (.not. valid) return
        valid = all(triangles >= 1) .and. &
            all(triangles <= size(vertices, 2))
    end function valid_surface_direction

end module fortfem_laplace_panel_pair_3d
