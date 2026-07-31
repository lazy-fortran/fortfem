module fortfem_fci_terminal_triangle_3d
    !! Exact first-hit search for a straight FCI trace in a 3D triangle set.
    !!
    !! The trace is parameterized as `x(t) = start + t*(finish-start)` with
    !! `0 < t <= 1`.  Triangles are oriented by their vertex ordering and the
    !! returned normal is the corresponding right-hand unit normal.  A valid
    !! search with no intersection returns status zero and hit_triangle zero.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none

    private

    real(dp), parameter :: topology_tolerance = 64.0_dp*epsilon(1.0_dp)

    public :: find_fci_first_hit_triangle_3d
    public :: find_fci_first_hit_triangle_3d_jvp

contains

    subroutine find_fci_first_hit_triangle_3d( &
            start_point, finish_point, surface_vertices, surface_triangles, &
            hit_point, hit_parameter, hit_triangle, normal, status)
        !! Find the nearest transverse triangle intersection on a trace.
        real(dp), intent(in) :: start_point(3), finish_point(3)
        real(dp), intent(in) :: surface_vertices(:, :)
        integer, intent(in) :: surface_triangles(:, :)
        real(dp), intent(out) :: hit_point(3), hit_parameter, normal(3)
        integer, intent(out) :: hit_triangle, status

        integer :: triangle, first, second, third, vertex_count, triangle_count
        real(dp) :: direction(3), edge_first(3), edge_second(3), area_vector(3)
        real(dp) :: direction_norm, area_norm, scale, denominator, numerator
        real(dp) :: parameter, u, v, inverse_denominator
        real(dp) :: offset(3), p_vector(3), q_vector(3), best_parameter

        hit_point = 0.0_dp
        hit_parameter = 0.0_dp
        hit_triangle = 0
        normal = 0.0_dp
        status = 1
        vertex_count = size(surface_vertices, 2)
        triangle_count = size(surface_triangles, 2)
        if (vertex_count < 3 .or. triangle_count < 1) return
        if (size(surface_vertices, 1) /= 3 .or. &
            size(surface_triangles, 1) /= 3) return
        if (any(.not. ieee_is_finite(start_point)) .or. &
            any(.not. ieee_is_finite(finish_point)) .or. &
            any(.not. ieee_is_finite(surface_vertices))) return
        if (any(surface_triangles < 1) .or. &
            any(surface_triangles > vertex_count)) return

        direction = finish_point - start_point
        direction_norm = sqrt(dot_product(direction, direction))
        scale = max(1.0_dp, maxval(abs(surface_vertices)), maxval(abs(start_point)), &
            maxval(abs(finish_point)))
        if (direction_norm <= topology_tolerance*scale) return

        do triangle = 1, triangle_count
            first = surface_triangles(1, triangle)
            second = surface_triangles(2, triangle)
            third = surface_triangles(3, triangle)
            edge_first = surface_vertices(:, second) - surface_vertices(:, first)
            edge_second = surface_vertices(:, third) - surface_vertices(:, first)
            area_vector = cross_3d(edge_first, edge_second)
            area_norm = sqrt(dot_product(area_vector, area_vector))
            if (area_norm <= topology_tolerance*scale*scale) return
        end do

        best_parameter = huge(1.0_dp)
        do triangle = 1, triangle_count
            first = surface_triangles(1, triangle)
            second = surface_triangles(2, triangle)
            third = surface_triangles(3, triangle)
            edge_first = surface_vertices(:, second) - surface_vertices(:, first)
            edge_second = surface_vertices(:, third) - surface_vertices(:, first)
            area_vector = cross_3d(edge_first, edge_second)
            area_norm = sqrt(dot_product(area_vector, area_vector))
            denominator = dot_product(direction, area_vector)
            if (abs(denominator) <= topology_tolerance*direction_norm*area_norm) cycle

            offset = start_point - surface_vertices(:, first)
            numerator = dot_product(-offset, area_vector)
            parameter = numerator/denominator
            if (parameter <= topology_tolerance .or. &
                parameter > 1.0_dp + topology_tolerance) cycle

            p_vector = cross_3d(direction, edge_second)
            ! Moller--Trumbore's determinant is -d.n for this orientation.
            inverse_denominator = -1.0_dp/denominator
            u = dot_product(offset, p_vector)*inverse_denominator
            q_vector = cross_3d(offset, edge_first)
            v = dot_product(direction, q_vector)*inverse_denominator
            if (u < -topology_tolerance .or. v < -topology_tolerance .or. &
                u + v > 1.0_dp + topology_tolerance) cycle
            if (parameter < best_parameter - topology_tolerance) then
                best_parameter = parameter
                hit_triangle = triangle
                hit_parameter = min(1.0_dp, max(0.0_dp, parameter))
                hit_point = start_point + hit_parameter*direction
                normal = area_vector/area_norm
            end if
        end do
        status = 0
    end subroutine find_fci_first_hit_triangle_3d

    subroutine find_fci_first_hit_triangle_3d_jvp( &
            start_point, finish_point, surface_vertices, surface_triangles, &
            start_point_dot, finish_point_dot, surface_vertices_dot, hit_point_dot, &
            hit_parameter_dot, normal_dot, status)
        !! Apply the fixed-topology JVP of the first-hit geometry.
        real(dp), intent(in) :: start_point(3), finish_point(3)
        real(dp), intent(in) :: surface_vertices(:, :)
        integer, intent(in) :: surface_triangles(:, :)
        real(dp), intent(in) :: start_point_dot(3), finish_point_dot(3)
        real(dp), intent(in) :: surface_vertices_dot(:, :)
        real(dp), intent(out) :: hit_point_dot(3), hit_parameter_dot, normal_dot(3)
        integer, intent(out) :: status

        real(dp) :: hit_point(3), hit_parameter, normal(3)
        real(dp) :: direction(3), direction_dot(3), edge_first(3), edge_second(3)
        real(dp) :: edge_first_dot(3), edge_second_dot(3), area_vector(3)
        real(dp) :: area_vector_dot(3), area_norm, area_norm_dot
        real(dp) :: offset(3), offset_dot(3), denominator, denominator_dot
        real(dp) :: numerator, numerator_dot
        integer :: hit_triangle, first, second, third, primal_status

        hit_point_dot = 0.0_dp
        hit_parameter_dot = 0.0_dp
        normal_dot = 0.0_dp
        status = 1
        if (size(surface_vertices, 1) /= 3 .or. &
            size(surface_vertices_dot, 1) /= 3 .or. &
            any(shape(surface_vertices_dot) /= shape(surface_vertices)) .or. &
            size(surface_triangles, 1) /= 3) return
        if (any(.not. ieee_is_finite(start_point_dot)) .or. &
            any(.not. ieee_is_finite(finish_point_dot)) .or. &
            any(.not. ieee_is_finite(surface_vertices_dot))) return

        call find_fci_first_hit_triangle_3d( &
            start_point, finish_point, surface_vertices, surface_triangles, &
            hit_point, hit_parameter, hit_triangle, normal, primal_status)
        if (primal_status /= 0) return
        if (hit_triangle == 0) then
            status = 0
            return
        end if

        first = surface_triangles(1, hit_triangle)
        second = surface_triangles(2, hit_triangle)
        third = surface_triangles(3, hit_triangle)
        direction = finish_point - start_point
        direction_dot = finish_point_dot - start_point_dot
        edge_first = surface_vertices(:, second) - surface_vertices(:, first)
        edge_second = surface_vertices(:, third) - surface_vertices(:, first)
        edge_first_dot = surface_vertices_dot(:, second) - &
            surface_vertices_dot(:, first)
        edge_second_dot = surface_vertices_dot(:, third) - &
            surface_vertices_dot(:, first)
        area_vector = cross_3d(edge_first, edge_second)
        area_vector_dot = cross_3d(edge_first_dot, edge_second) + &
            cross_3d(edge_first, edge_second_dot)
        area_norm = sqrt(dot_product(area_vector, area_vector))
        area_norm_dot = dot_product(area_vector, area_vector_dot)/area_norm
        offset = surface_vertices(:, first) - start_point
        offset_dot = surface_vertices_dot(:, first) - start_point_dot
        denominator = dot_product(direction, area_vector)
        denominator_dot = dot_product(direction_dot, area_vector) + &
            dot_product(direction, area_vector_dot)
        numerator = dot_product(offset, area_vector)
        numerator_dot = dot_product(offset_dot, area_vector) + &
            dot_product(offset, area_vector_dot)
        if (abs(denominator) <= topology_tolerance* &
            sqrt(dot_product(direction, direction))*area_norm) return
        hit_parameter_dot = (numerator_dot*denominator - numerator*denominator_dot)/ &
            denominator**2
        hit_point_dot = start_point_dot + hit_parameter_dot*direction + &
            hit_parameter*direction_dot
        normal_dot = area_vector_dot/area_norm - &
            normal*area_norm_dot/area_norm
        status = 0
    end subroutine find_fci_first_hit_triangle_3d_jvp

    pure function cross_3d(first, second) result(cross_product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: cross_product(3)

        cross_product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_3d

end module fortfem_fci_terminal_triangle_3d
