program test_laplace_representation_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_laplace_representation_triangles_3d, &
        evaluate_toroidal_harmonic_p, generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: cells = 30
    integer, parameter :: triangle_count = 12*cells*cells
    integer, parameter :: vertex_count = 3*triangle_count
    real(dp) :: dirichlet(vertex_count), neumann(triangle_count)
    real(dp) :: target(3), value, vertices(3, vertex_count)
    integer :: triangles(3, triangle_count), status
    logical :: all_passed

    all_passed = .true.
    call build_cube_surface(vertices, triangles, dirichlet, neumann)
    target = [3.0_dp, 0.2_dp, -0.1_dp]
    call evaluate_laplace_representation_triangles_3d( &
        vertices, triangles, dirichlet, neumann, target, 12, value, status)
    call record_condition(status == 0, &
        "Three-dimensional Laplace representation evaluates")
    call record_condition(abs(value - 1.0_dp/norm2(target)) < 7.0e-5_dp, &
        "Cube Cauchy data reproduce the exact exterior monopole")
    call test_toroidal_harmonic_representation()

    call check_summary("Three-dimensional Laplace representation")
    if (.not. all_passed) error stop 1

contains

    subroutine test_toroidal_harmonic_representation()
        real(dp), parameter :: major_radius = 2.0_dp
        real(dp), parameter :: minor_radius = 0.6_dp
        real(dp), parameter :: scale = &
            sqrt(major_radius**2 - minor_radius**2)
        real(dp), parameter :: boundary_eta = &
            acosh(major_radius/minor_radius)
        real(dp), allocatable :: torus_vertices(:, :)
        integer, allocatable :: torus_triangles(:, :)
        real(dp), allocatable :: torus_trace(:), torus_flux(:)
        real(dp) :: centroid(3), exact, numerical, phi, theta
        real(dp) :: target_eta, target_phi, target_theta
        real(dp) :: denominator, eta_at_point, first_edge(3), normal(3)
        real(dp) :: second_edge(3)
        integer :: element, point, torus_status

        call generate_torus_surface_mesh( &
            major_radius, minor_radius, 144, 112, &
            torus_vertices, torus_triangles)
        allocate(torus_trace(size(torus_vertices, 2)))
        allocate(torus_flux(size(torus_triangles, 2)))
        do point = 1, size(torus_vertices, 2)
            call inverse_toroidal_angles( &
                torus_vertices(:, point), scale, eta_at_point, theta, phi)
            call evaluate_toroidal_harmonic_p( &
                2, 1, boundary_eta, theta, phi, &
                torus_trace(point), torus_status)
            if (torus_status /= 0) error stop "toroidal trace failed"
        end do
        do element = 1, size(torus_triangles, 2)
            centroid = sum( &
                torus_vertices(:, torus_triangles(:, element)), dim=2)/3.0_dp
            first_edge = torus_vertices(:, torus_triangles(2, element)) - &
                torus_vertices(:, torus_triangles(1, element))
            second_edge = torus_vertices(:, torus_triangles(3, element)) - &
                torus_vertices(:, torus_triangles(1, element))
            normal = cross_product(first_edge, second_edge)
            normal = normal/norm2(normal)
            torus_flux(element) = &
                toroidal_normal_difference(centroid, normal, scale)
        end do

        target_eta = 0.35_dp
        target_theta = 0.8_dp
        target_phi = 0.3_dp
        denominator = cosh(target_eta) - cos(target_theta)
        target = [ &
            scale*sinh(target_eta)*cos(target_phi)/denominator, &
            scale*sinh(target_eta)*sin(target_phi)/denominator, &
            scale*sin(target_theta)/denominator]
        call evaluate_toroidal_harmonic_p( &
            2, 1, target_eta, target_theta, target_phi, exact, torus_status)
        call evaluate_laplace_representation_triangles_3d( &
            torus_vertices, torus_triangles, torus_trace, torus_flux, &
            target, 10, numerical, torus_status)
        call record_condition(torus_status == 0 .and. &
            abs(numerical - exact) < 2.0e-2_dp*abs(exact), &
            "Toroidal BEM representation matches the half-integer harmonic")
    end subroutine test_toroidal_harmonic_representation

    function toroidal_normal_difference( &
            point, normal, toroidal_scale) result(derivative)
        real(dp), intent(in) :: point(3), normal(3), toroidal_scale
        real(dp) :: derivative

        real(dp), parameter :: step = 1.0e-5_dp
        real(dp) :: eta_local, minus_value, phi_local, plus_value, theta_local
        integer :: local_status

        call inverse_toroidal_angles( &
            point + step*normal, toroidal_scale, &
            eta_local, theta_local, phi_local)
        call evaluate_toroidal_harmonic_p( &
            2, 1, eta_local, theta_local, phi_local, plus_value, local_status)
        if (local_status /= 0) error stop "positive toroidal difference failed"
        call inverse_toroidal_angles( &
            point - step*normal, toroidal_scale, &
            eta_local, theta_local, phi_local)
        call evaluate_toroidal_harmonic_p( &
            2, 1, eta_local, theta_local, phi_local, minus_value, local_status)
        if (local_status /= 0) error stop "negative toroidal difference failed"
        derivative = (plus_value - minus_value)/(2.0_dp*step)
    end function toroidal_normal_difference

    pure subroutine inverse_toroidal_angles(point, scale, eta, theta, phi)
        real(dp), intent(in) :: point(3), scale
        real(dp), intent(out) :: eta, theta, phi

        real(dp) :: cylindrical_radius

        cylindrical_radius = norm2(point(1:2))
        eta = 0.5_dp*log( &
            ((cylindrical_radius + scale)**2 + point(3)**2)/ &
            ((cylindrical_radius - scale)**2 + point(3)**2))
        theta = atan2( &
            2.0_dp*scale*point(3), &
            cylindrical_radius**2 + point(3)**2 - scale**2)
        phi = atan2(point(2), point(1))
    end subroutine inverse_toroidal_angles

    subroutine build_cube_surface( &
            points, elements, boundary_value, boundary_flux)
        real(dp), intent(out) :: points(:, :)
        integer, intent(out) :: elements(:, :)
        real(dp), intent(out) :: boundary_value(:), boundary_flux(:)

        real(dp), parameter :: normals(3, 6) = reshape([ &
            1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp], [3, 6])
        real(dp) :: corners(3, 4), centroid(3), first(3), second(3)
        real(dp) :: lower, upper
        integer :: cell_x, cell_y, face, local_triangle, triangle, vertex

        triangle = 0
        vertex = 0
        do face = 1, 6
            do cell_y = 0, cells - 1
                do cell_x = 0, cells - 1
                    lower = -1.0_dp + 2.0_dp*real(cell_x, dp)/real(cells, dp)
                    upper = -1.0_dp + 2.0_dp*real(cell_y, dp)/real(cells, dp)
                    call face_cell(face, lower, upper, corners)
                    do local_triangle = 1, 2
                        triangle = triangle + 1
                        if (local_triangle == 1) then
                            first = corners(:, 1)
                            second = corners(:, 2)
                            call add_triangle( &
                                first, second, corners(:, 3), &
                                normals(:, face), triangle, vertex, &
                                points, elements, boundary_value)
                        else
                            first = corners(:, 3)
                            second = corners(:, 2)
                            call add_triangle( &
                                first, second, corners(:, 4), &
                                normals(:, face), triangle, vertex, &
                                points, elements, boundary_value)
                        end if
                        centroid = sum(points(:, elements(:, triangle)), dim=2)/ &
                            3.0_dp
                        boundary_flux(triangle) = &
                            -dot_product(normals(:, face), centroid)/norm2(centroid)**3
                    end do
                end do
            end do
        end do

    end subroutine build_cube_surface

    subroutine add_triangle( &
            a, b, c, outward_normal, triangle, vertex, &
            points, elements, boundary_value)
        real(dp), intent(in) :: a(3), b(3), c(3), outward_normal(3)
        integer, intent(in) :: triangle
        integer, intent(inout) :: vertex
        real(dp), intent(inout) :: points(:, :), boundary_value(:)
        integer, intent(inout) :: elements(:, :)

        real(dp) :: oriented(3, 3), cross(3)

        oriented = reshape([a, b, c], [3, 3])
        cross = cross_product(b - a, c - a)
        if (dot_product(cross, outward_normal) < 0.0_dp) then
            oriented(:, 2) = c
            oriented(:, 3) = b
        end if
        elements(:, triangle) = [vertex + 1, vertex + 2, vertex + 3]
        points(:, vertex + 1:vertex + 3) = oriented
        boundary_value(vertex + 1:vertex + 3) = [ &
            1.0_dp/norm2(oriented(:, 1)), &
            1.0_dp/norm2(oriented(:, 2)), &
            1.0_dp/norm2(oriented(:, 3))]
        vertex = vertex + 3
    end subroutine add_triangle

    subroutine face_cell(face, first_coordinate, second_coordinate, corners)
        integer, intent(in) :: face
        real(dp), intent(in) :: first_coordinate, second_coordinate
        real(dp), intent(out) :: corners(3, 4)

        real(dp) :: step

        step = 2.0_dp/real(cells, dp)
        select case (face)
        case (1, 2)
            corners = reshape([ &
                merge(1.0_dp, -1.0_dp, face == 1), first_coordinate, second_coordinate, &
                merge(1.0_dp, -1.0_dp, face == 1), first_coordinate + step, second_coordinate, &
                merge(1.0_dp, -1.0_dp, face == 1), first_coordinate, second_coordinate + step, &
                merge(1.0_dp, -1.0_dp, face == 1), first_coordinate + step, second_coordinate + step], [3, 4])
        case (3, 4)
            corners = reshape([ &
                first_coordinate, merge(1.0_dp, -1.0_dp, face == 3), second_coordinate, &
                first_coordinate + step, merge(1.0_dp, -1.0_dp, face == 3), second_coordinate, &
                first_coordinate, merge(1.0_dp, -1.0_dp, face == 3), second_coordinate + step, &
                first_coordinate + step, merge(1.0_dp, -1.0_dp, face == 3), second_coordinate + step], [3, 4])
        case default
            corners = reshape([ &
                first_coordinate, second_coordinate, merge(1.0_dp, -1.0_dp, face == 5), &
                first_coordinate + step, second_coordinate, merge(1.0_dp, -1.0_dp, face == 5), &
                first_coordinate, second_coordinate + step, merge(1.0_dp, -1.0_dp, face == 5), &
                first_coordinate + step, second_coordinate + step, merge(1.0_dp, -1.0_dp, face == 5)], [3, 4])
        end select
    end subroutine face_cell

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        call check_condition(condition, description)
        all_passed = all_passed .and. condition
    end subroutine record_condition
end program test_laplace_representation_3d
