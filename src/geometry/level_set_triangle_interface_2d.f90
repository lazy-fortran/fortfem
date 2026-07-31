module fortfem_level_set_triangle_interface_2d
    !! Linear level-set intersection of a physical triangle.
    !!
    !! The level set is interpolated linearly from its three nodal values.  A
    !! proper cut returns the two edge intersections, their physical length,
    !! and the unit normal in the direction of the physical level-set gradient.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none

    private

    real(dp), parameter :: topology_tolerance = 64.0_dp*epsilon(1.0_dp)

    public :: evaluate_level_set_triangle_interface_2d
    public :: evaluate_level_set_triangle_cut_areas_2d
    public :: evaluate_level_set_triangle_cut_quadrature_2d

contains

    subroutine evaluate_level_set_triangle_interface_2d( &
            vertices, level_values, points, length, normal, status)
        real(dp), intent(in) :: vertices(2, 3), level_values(3)
        real(dp), intent(out) :: points(2, 2), length, normal(2)
        integer, intent(out) :: status

        integer, parameter :: edge_start(3) = [1, 1, 2]
        integer, parameter :: edge_end(3) = [2, 3, 3]
        integer :: edge, first, second, point_count
        real(dp) :: first_value, second_value, fraction
        real(dp) :: edge_point(2), edge_vector(2), first_edge(2), second_edge(2)
        real(dp) :: gradient(2), determinant, gradient_norm

        points = 0.0_dp
        length = 0.0_dp
        normal = 0.0_dp
        status = 1
        if (any(.not. ieee_is_finite(vertices)) .or. &
            any(.not. ieee_is_finite(level_values))) return
        first_edge = vertices(:, 2) - vertices(:, 1)
        second_edge = vertices(:, 3) - vertices(:, 1)
        determinant = first_edge(1)*second_edge(2) - &
            first_edge(2)*second_edge(1)
        if (abs(determinant) <= topology_tolerance) return

        point_count = 0
        do edge = 1, 3
            first = edge_start(edge)
            second = edge_end(edge)
            first_value = level_values(first)
            second_value = level_values(second)
            if (abs(first_value) <= topology_tolerance) then
                call add_unique_point( &
                    vertices(:, first), points, point_count)
            end if
            if (first_value*second_value < 0.0_dp) then
                fraction = first_value/(first_value - second_value)
                edge_vector = vertices(:, second) - vertices(:, first)
                edge_point = vertices(:, first) + fraction*edge_vector
                call add_unique_point(edge_point, points, point_count)
            end if
            if (abs(second_value) <= topology_tolerance) then
                call add_unique_point( &
                    vertices(:, second), points, point_count)
            end if
        end do
        if (point_count /= 2) return

        gradient(1) = ((level_values(2) - level_values(1))*second_edge(2) - &
            (level_values(3) - level_values(1))*first_edge(2))/determinant
        gradient(2) = (first_edge(1)*(level_values(3) - level_values(1)) - &
            second_edge(1)*(level_values(2) - level_values(1)))/determinant
        gradient_norm = sqrt(dot_product(gradient, gradient))
        if (gradient_norm <= topology_tolerance) return
        edge_vector = points(:, 2) - points(:, 1)
        length = sqrt(dot_product(edge_vector, edge_vector))
        if (length <= topology_tolerance) then
            length = 0.0_dp
            return
        end if
        normal = gradient/gradient_norm
        status = 0
    end subroutine evaluate_level_set_triangle_interface_2d

    subroutine evaluate_level_set_triangle_cut_areas_2d( &
            vertices, level_values, positive_area, negative_area, &
            interface_length, status)
        !! Return exact subcell areas for a linearly interpolated level set.
        real(dp), intent(in) :: vertices(2, 3), level_values(3)
        real(dp), intent(out) :: positive_area, negative_area
        real(dp), intent(out) :: interface_length
        integer, intent(out) :: status

        real(dp) :: points(2, 2), normal(2), determinant
        integer :: interface_status

        positive_area = 0.0_dp
        negative_area = 0.0_dp
        interface_length = 0.0_dp
        status = 1
        if (any(.not. ieee_is_finite(vertices)) .or. &
            any(.not. ieee_is_finite(level_values))) return
        determinant = (vertices(2, 2) - vertices(2, 1))* &
            (vertices(1, 3) - vertices(1, 1)) - &
            (vertices(1, 2) - vertices(1, 1))* &
            (vertices(2, 3) - vertices(2, 1))
        if (abs(determinant) <= topology_tolerance) return

        positive_area = level_set_side_area(vertices, level_values, .true.)
        negative_area = level_set_side_area(vertices, level_values, .false.)
        call evaluate_level_set_triangle_interface_2d( &
            vertices, level_values, points, interface_length, normal, &
            interface_status)
        if (interface_status /= 0) then
            if (.not. all(level_values >= -topology_tolerance) .and. &
                .not. all(level_values <= topology_tolerance)) then
                positive_area = 0.0_dp
                negative_area = 0.0_dp
                return
            end if
            interface_length = 0.0_dp
        end if
        status = 0
    end subroutine evaluate_level_set_triangle_cut_areas_2d

    subroutine evaluate_level_set_triangle_cut_quadrature_2d( &
            vertices, level_values, positive_area, positive_centroid, &
            negative_area, negative_centroid, interface_length, normal, status)
        !! Return exact degree-one quadrature data for a linear level-set cut.
        !!
        !! The positive and negative regions are clipped polygons.  Their
        !! areas and centroids therefore integrate constants and affine fields
        !! exactly.  A zero-area side has a zero centroid.  The interface
        !! normal follows the same physical-gradient orientation as the
        !! interface primitive; no normal is returned for an uncut triangle.
        real(dp), intent(in) :: vertices(2, 3), level_values(3)
        real(dp), intent(out) :: positive_area, positive_centroid(2)
        real(dp), intent(out) :: negative_area, negative_centroid(2)
        real(dp), intent(out) :: interface_length, normal(2)
        integer, intent(out) :: status

        real(dp) :: determinant, points(2, 2)
        integer :: interface_status

        positive_area = 0.0_dp
        positive_centroid = 0.0_dp
        negative_area = 0.0_dp
        negative_centroid = 0.0_dp
        interface_length = 0.0_dp
        normal = 0.0_dp
        status = 1
        if (any(.not. ieee_is_finite(vertices)) .or. &
            any(.not. ieee_is_finite(level_values))) return
        if (all(abs(level_values) <= topology_tolerance)) return
        determinant = (vertices(2, 2) - vertices(2, 1))* &
            (vertices(1, 3) - vertices(1, 1)) - &
            (vertices(1, 2) - vertices(1, 1))* &
            (vertices(2, 3) - vertices(2, 1))
        if (abs(determinant) <= topology_tolerance) return

        call level_set_side_moments( &
            vertices, level_values, .true., positive_area, positive_centroid)
        call level_set_side_moments( &
            vertices, level_values, .false., negative_area, negative_centroid)
        call evaluate_level_set_triangle_interface_2d( &
            vertices, level_values, points, interface_length, normal, &
            interface_status)
        if (interface_status /= 0) then
            if (.not. all(level_values >= -topology_tolerance) .and. &
                .not. all(level_values <= topology_tolerance)) then
                positive_area = 0.0_dp
                positive_centroid = 0.0_dp
                negative_area = 0.0_dp
                negative_centroid = 0.0_dp
                return
            end if
        end if
        status = 0
    end subroutine evaluate_level_set_triangle_cut_quadrature_2d

    pure function level_set_side_area(vertices, level_values, keep_positive) &
            result(area)
        real(dp), intent(in) :: vertices(2, 3), level_values(3)
        logical, intent(in) :: keep_positive
        real(dp) :: area
        real(dp) :: centroid(2)

        call level_set_side_moments( &
            vertices, level_values, keep_positive, area, centroid)
    end function level_set_side_area

    pure subroutine level_set_side_moments( &
            vertices, level_values, keep_positive, area, centroid)
        real(dp), intent(in) :: vertices(2, 3), level_values(3)
        logical, intent(in) :: keep_positive
        real(dp), intent(out) :: area, centroid(2)

        real(dp) :: polygon(2, 6), clipped(2, 6)
        real(dp) :: signed_area, cross, first_moment(2)
        integer :: point_count, point, next_point

        area = 0.0_dp
        centroid = 0.0_dp
        polygon(:, 1:3) = vertices
        point_count = 3
        call clip_level_set_polygon( &
            polygon, point_count, level_values, keep_positive, clipped)
        if (point_count < 3) return
        signed_area = 0.0_dp
        first_moment = 0.0_dp
        do point = 1, point_count
            next_point = 1 + mod(point, point_count)
            cross = polygon(1, point)*polygon(2, next_point) - &
                polygon(2, point)*polygon(1, next_point)
            signed_area = signed_area + 0.5_dp*cross
            first_moment = first_moment + 0.5_dp*cross*( &
                polygon(:, point) + polygon(:, next_point))
        end do
        if (abs(signed_area) <= topology_tolerance) return
        area = abs(signed_area)
        centroid = first_moment/(3.0_dp*signed_area)
    end subroutine level_set_side_moments

    pure subroutine clip_level_set_polygon( &
            polygon, point_count, level_values, keep_positive, clipped)
        real(dp), intent(inout) :: polygon(2, 6)
        integer, intent(inout) :: point_count
        real(dp), intent(in) :: level_values(3)
        logical, intent(in) :: keep_positive
        real(dp), intent(out) :: clipped(2, 6)

        integer :: current, previous, clipped_count
        logical :: current_inside, previous_inside
        real(dp) :: current_value, previous_value, fraction
        real(dp) :: intersection(2)

        clipped = 0.0_dp
        clipped_count = 0
        previous = 3
        previous_value = level_values(previous)
        previous_inside = side_inside(previous_value, keep_positive)
        do current = 1, 3
            current_value = level_values(current)
            current_inside = side_inside(current_value, keep_positive)
            if (current_inside .and. .not. previous_inside) then
                fraction = previous_value/(previous_value - current_value)
                intersection = polygon(:, previous) + fraction* &
                    (polygon(:, current) - polygon(:, previous))
                clipped_count = clipped_count + 1
                clipped(:, clipped_count) = intersection
            else if (.not. current_inside .and. previous_inside) then
                fraction = previous_value/(previous_value - current_value)
                intersection = polygon(:, previous) + fraction* &
                    (polygon(:, current) - polygon(:, previous))
                clipped_count = clipped_count + 1
                clipped(:, clipped_count) = intersection
            end if
            if (current_inside) then
                clipped_count = clipped_count + 1
                clipped(:, clipped_count) = polygon(:, current)
            end if
            previous = current
            previous_value = current_value
            previous_inside = current_inside
        end do
        polygon = clipped
        point_count = clipped_count
    end subroutine clip_level_set_polygon

    pure logical function side_inside(value, keep_positive) result(inside)
        real(dp), intent(in) :: value
        logical, intent(in) :: keep_positive

        if (keep_positive) then
            inside = value >= -topology_tolerance
        else
            inside = value <= topology_tolerance
        end if
    end function side_inside

    pure subroutine add_unique_point(candidate, points, point_count)
        real(dp), intent(in) :: candidate(2)
        real(dp), intent(inout) :: points(2, 2)
        integer, intent(inout) :: point_count

        integer :: point

        do point = 1, min(point_count, 2)
            if (maxval(abs(candidate - points(:, point))) <= &
                topology_tolerance) return
        end do
        if (point_count < 2) then
            point_count = point_count + 1
            points(:, point_count) = candidate
        end if
    end subroutine add_unique_point

end module fortfem_level_set_triangle_interface_2d
