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

    pure function cross_product(first, second) result(cross)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: cross(3)

        cross = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_level_set_tetra_interface_3d
