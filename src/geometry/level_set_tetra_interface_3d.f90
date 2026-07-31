module fortfem_level_set_tetra_interface_3d
    !! Linear level-set intersection of a physical tetrahedron.
    !!
    !! A proper cut returns an ordered triangular or quadrilateral polygon on
    !! the zero level set, its physical area, and the unit normal in the
    !! direction of the affine physical level-set gradient.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_linalg, only: det3, inv3
    use fortfem_kinds, only: dp
    implicit none

    private

    real(dp), parameter :: topology_tolerance = 64.0_dp*epsilon(1.0_dp)

    public :: evaluate_level_set_tetra_interface_3d
    public :: evaluate_level_set_tetra_cut_quadrature_3d

contains

    subroutine evaluate_level_set_tetra_interface_3d( &
            vertices, level_values, points, point_count, area, normal, status)
        !! Return an ordered polygon for a proper linear tetrahedron cut.
        real(dp), intent(in) :: vertices(3, 4), level_values(4)
        real(dp), intent(out) :: points(3, 4), area, normal(3)
        integer, intent(out) :: point_count, status

        integer, parameter :: edge_start(6) = [1, 1, 1, 2, 2, 3]
        integer, parameter :: edge_end(6) = [2, 3, 4, 3, 4, 4]
        integer :: edge, first, second, inverse_status
        real(dp) :: jacobian(3, 3), inverse(3, 3), determinant
        real(dp) :: gradient_reference(3), gradient(3), gradient_norm
        real(dp) :: first_value, second_value, fraction, edge_vector(3)
        real(dp) :: edge_point(3), centroid(3), basis_u(3), basis_v(3)
        real(dp) :: reference_axis(3), area_vector(3), centered(3)
        real(dp) :: angles(4), swap_angle, swap_point(3)
        real(dp) :: next_centered(3)
        integer :: point, next_point, order

        points = 0.0_dp
        point_count = 0
        area = 0.0_dp
        normal = 0.0_dp
        status = 1
        if (any(.not. ieee_is_finite(vertices)) .or. &
            any(.not. ieee_is_finite(level_values))) return
        if (any(abs(level_values) <= topology_tolerance)) return

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
        determinant = det3(jacobian)
        if (abs(determinant) <= topology_tolerance* &
            max(1.0_dp, maxval(abs(jacobian))**3)) return
        call inv3(jacobian, inverse, inverse_status)
        if (inverse_status /= 0) return

        gradient_reference = level_values(2:4) - level_values(1)
        gradient = matmul(transpose(inverse), gradient_reference)
        gradient_norm = sqrt(dot_product(gradient, gradient))
        if (gradient_norm <= topology_tolerance) return
        normal = gradient/gradient_norm

        do edge = 1, 6
            first = edge_start(edge)
            second = edge_end(edge)
            first_value = level_values(first)
            second_value = level_values(second)
            if (first_value*second_value < 0.0_dp) then
                fraction = first_value/(first_value - second_value)
                edge_vector = vertices(:, second) - vertices(:, first)
                edge_point = vertices(:, first) + fraction*edge_vector
                if (point_count >= 4) return
                point_count = point_count + 1
                points(:, point_count) = edge_point
            end if
        end do
        if (point_count < 3) then
            point_count = 0
            points = 0.0_dp
            return
        end if

        centroid = sum(points(:, :point_count), dim=2)/real(point_count, dp)
        if (abs(normal(1)) < 0.9_dp) then
            reference_axis = [1.0_dp, 0.0_dp, 0.0_dp]
        else
            reference_axis = [0.0_dp, 1.0_dp, 0.0_dp]
        end if
        basis_u = cross_product(normal, reference_axis)
        basis_u = basis_u/sqrt(dot_product(basis_u, basis_u))
        basis_v = cross_product(normal, basis_u)
        do point = 1, point_count
            centered = points(:, point) - centroid
            angles(point) = atan2( &
                dot_product(centered, basis_v), dot_product(centered, basis_u))
        end do
        do point = 1, point_count - 1
            order = point
            do next_point = point + 1, point_count
                if (angles(next_point) < angles(order)) order = next_point
            end do
            if (order /= point) then
                swap_angle = angles(point)
                angles(point) = angles(order)
                angles(order) = swap_angle
                swap_point = points(:, point)
                points(:, point) = points(:, order)
                points(:, order) = swap_point
            end if
        end do

        area_vector = 0.0_dp
        do point = 1, point_count
            next_point = 1 + mod(point, point_count)
            centered = points(:, point) - centroid
            next_centered = points(:, next_point) - centroid
            area_vector = area_vector + cross_product(centered, next_centered)
        end do
        area = 0.5_dp*sqrt(dot_product(area_vector, area_vector))
        if (area <= topology_tolerance) then
            points = 0.0_dp
            point_count = 0
            area = 0.0_dp
            normal = 0.0_dp
            return
        end if
        status = 0
    end subroutine evaluate_level_set_tetra_interface_3d

    subroutine evaluate_level_set_tetra_cut_quadrature_3d( &
            vertices, level_values, positive_volume, positive_centroid, &
            negative_volume, negative_centroid, interface_area, normal, status)
        !! Return exact degree-one quadrature data for a linear tetrahedral cut.
        !!
        !! The two side regions are clipped convex polyhedra.  Their volumes and
        !! first moments are accumulated from the oriented parent faces and the
        !! interface polygon, so constants and affine fields are integrated
        !! exactly without a volumetric numerical quadrature rule.
        real(dp), intent(in) :: vertices(3, 4), level_values(4)
        real(dp), intent(out) :: positive_volume, positive_centroid(3)
        real(dp), intent(out) :: negative_volume, negative_centroid(3)
        real(dp), intent(out) :: interface_area, normal(3)
        integer, intent(out) :: status

        real(dp) :: jacobian(3, 3), determinant
        real(dp) :: points(3, 4)
        logical :: has_positive, has_negative
        integer :: interface_status, point_count

        positive_volume = 0.0_dp
        positive_centroid = 0.0_dp
        negative_volume = 0.0_dp
        negative_centroid = 0.0_dp
        interface_area = 0.0_dp
        normal = 0.0_dp
        status = 1
        if (any(.not. ieee_is_finite(vertices)) .or. &
            any(.not. ieee_is_finite(level_values))) return
        if (any(abs(level_values) <= topology_tolerance)) return
        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
        determinant = det3(jacobian)
        if (abs(determinant) <= topology_tolerance* &
            max(1.0_dp, maxval(abs(jacobian))**3)) return

        has_positive = any(level_values > 0.0_dp)
        has_negative = any(level_values < 0.0_dp)
        point_count = 0
        if (has_positive .and. has_negative) then
            call evaluate_level_set_tetra_interface_3d( &
                vertices, level_values, points, point_count, interface_area, &
                normal, interface_status)
            if (interface_status /= 0) return
        end if
        call accumulate_tetra_side( &
            vertices, level_values, .true., points, point_count, positive_volume, &
            positive_centroid, status)
        if (status /= 0) return
        call accumulate_tetra_side( &
            vertices, level_values, .false., points, point_count, negative_volume, &
            negative_centroid, status)
        if (status /= 0) return
        status = 0
    end subroutine evaluate_level_set_tetra_cut_quadrature_3d

    subroutine accumulate_tetra_side( &
            vertices, level_values, keep_positive, interface_points, &
            interface_count, volume, centroid, status)
        real(dp), intent(in) :: vertices(3, 4), level_values(4)
        logical, intent(in) :: keep_positive
        real(dp), intent(in) :: interface_points(3, 4)
        integer, intent(in) :: interface_count
        real(dp), intent(out) :: volume, centroid(3)
        integer, intent(out) :: status

        integer, parameter :: face_opposite(4) = [1, 2, 3, 4]
        integer :: face_vertices(3, 4), face, polygon_count
        integer :: indices(3), point
        real(dp) :: polygon(3, 6)
        real(dp) :: volume_accum, first_moment(3)

        face_vertices(:, 1) = [2, 3, 4]
        face_vertices(:, 2) = [1, 4, 3]
        face_vertices(:, 3) = [1, 2, 4]
        face_vertices(:, 4) = [1, 3, 2]
        volume = 0.0_dp
        centroid = 0.0_dp
        status = 1
        volume_accum = 0.0_dp
        first_moment = 0.0_dp
        do face = 1, 4
            indices = face_vertices(:, face)
            call orient_tetra_face(vertices, face_opposite(face), indices)
            call clip_tetra_face( &
                vertices, level_values, indices, keep_positive, polygon, &
                polygon_count)
            if (polygon_count >= 3) then
                call accumulate_oriented_polygon( &
                    polygon, polygon_count, volume_accum, first_moment)
            end if
        end do

        if (interface_count > 0) then
            if (keep_positive) then
                do point = 1, interface_count
                    polygon(:, point) = &
                        interface_points(:, interface_count - point + 1)
                end do
            else
                polygon(:, :interface_count) = &
                    interface_points(:, :interface_count)
            end if
            call accumulate_oriented_polygon( &
                polygon, interface_count, volume_accum, first_moment)
        end if

        if (volume_accum <= topology_tolerance) then
            if (volume_accum < -topology_tolerance) return
            volume = 0.0_dp
            centroid = 0.0_dp
            status = 0
            return
        end if
        volume = volume_accum
        centroid = first_moment/volume
        status = 0
    end subroutine accumulate_tetra_side

    subroutine orient_tetra_face(vertices, opposite, indices)
        real(dp), intent(in) :: vertices(3, 4)
        integer, intent(in) :: opposite
        integer, intent(inout) :: indices(3)

        real(dp) :: first_edge(3), second_edge(3), face_normal(3)
        integer :: swap_index

        first_edge = vertices(:, indices(2)) - vertices(:, indices(1))
        second_edge = vertices(:, indices(3)) - vertices(:, indices(1))
        face_normal = cross_product(first_edge, second_edge)
        if (dot_product(face_normal, vertices(:, opposite) - &
            vertices(:, indices(1))) > 0.0_dp) then
            swap_index = indices(2)
            indices(2) = indices(3)
            indices(3) = swap_index
        end if
    end subroutine orient_tetra_face

    subroutine clip_tetra_face( &
            vertices, level_values, indices, keep_positive, polygon, point_count)
        real(dp), intent(in) :: vertices(3, 4), level_values(4)
        integer, intent(in) :: indices(3)
        logical, intent(in) :: keep_positive
        real(dp), intent(out) :: polygon(3, 6)
        integer, intent(out) :: point_count

        integer :: current, previous, vertex
        real(dp) :: current_value, previous_value, fraction
        logical :: current_inside, previous_inside

        polygon = 0.0_dp
        point_count = 0
        previous = indices(3)
        previous_value = level_values(previous)
        previous_inside = side_contains(previous_value, keep_positive)
        do vertex = 1, 3
            current = indices(vertex)
            current_value = level_values(current)
            current_inside = side_contains(current_value, keep_positive)
            if (current_inside .neqv. previous_inside) then
                fraction = previous_value/(previous_value - current_value)
                point_count = point_count + 1
                polygon(:, point_count) = vertices(:, previous) + fraction* &
                    (vertices(:, current) - vertices(:, previous))
            end if
            if (current_inside) then
                point_count = point_count + 1
                polygon(:, point_count) = vertices(:, current)
            end if
            previous = current
            previous_value = current_value
            previous_inside = current_inside
        end do
    end subroutine clip_tetra_face

    pure logical function side_contains(value, keep_positive)
        real(dp), intent(in) :: value
        logical, intent(in) :: keep_positive

        side_contains = merge(value > 0.0_dp, value < 0.0_dp, keep_positive)
    end function side_contains

    subroutine accumulate_oriented_polygon( &
            polygon, point_count, volume, first_moment)
        real(dp), intent(in) :: polygon(3, 6)
        integer, intent(in) :: point_count
        real(dp), intent(inout) :: volume, first_moment(3)

        integer :: point
        real(dp) :: signed_tetra_volume, first_point(3)
        real(dp) :: second_point(3), third_point(3)

        first_point = polygon(:, 1)
        do point = 2, point_count - 1
            second_point = polygon(:, point)
            third_point = polygon(:, point + 1)
            signed_tetra_volume = dot_product(first_point, cross_product( &
                second_point, third_point))/6.0_dp
            volume = volume + signed_tetra_volume
            first_moment = first_moment + signed_tetra_volume*( &
                first_point + second_point + third_point)/4.0_dp
        end do
    end subroutine accumulate_oriented_polygon

    pure function cross_product(first, second) result(cross)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: cross(3)

        cross = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_level_set_tetra_interface_3d
