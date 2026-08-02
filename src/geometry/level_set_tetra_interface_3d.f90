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
    public :: evaluate_level_set_tetra_interface_3d_jvp
    public :: evaluate_level_set_tetra_cut_quadrature_3d
    public :: evaluate_level_set_tetra_cut_quadrature_3d_jvp
    public :: evaluate_level_set_tetra_cut_moments_3d
    public :: evaluate_level_set_tetra_cut_moments_3d_jvp
    public :: evaluate_level_set_tetra_cut_third_moments_3d
    public :: evaluate_level_set_tetra_cut_third_moments_3d_jvp
    public :: evaluate_level_set_tetra_cut_fourth_moments_3d
    public :: evaluate_level_set_tetra_cut_fourth_moments_3d_jvp

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

    subroutine evaluate_level_set_tetra_interface_3d_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, points_dot, &
            area_dot, normal_dot, status)
        !! Apply the fixed-topology JVP of a tetrahedral level-set interface.
        real(dp), intent(in) :: vertices(3, 4), level_values(4)
        real(dp), intent(in) :: vertices_dot(3, 4), level_values_dot(4)
        real(dp), intent(out) :: points_dot(3, 4), area_dot, normal_dot(3)
        integer, intent(out) :: status

        real(dp) :: points(3, 4), area, normal(3)
        integer :: point_count

        points_dot = 0.0_dp
        area_dot = 0.0_dp
        normal_dot = 0.0_dp
        call evaluate_tetra_interface_geometry_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, points, &
            points_dot, point_count, area, area_dot, normal, normal_dot, status)
    end subroutine evaluate_level_set_tetra_interface_3d_jvp

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

    subroutine evaluate_level_set_tetra_cut_quadrature_3d_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, &
            positive_volume_dot, positive_centroid_dot, negative_volume_dot, &
            negative_centroid_dot, interface_area_dot, normal_dot, status)
        !! Apply the fixed-topology JVP of tetrahedral cut quadrature data.
        real(dp), intent(in) :: vertices(3, 4), level_values(4)
        real(dp), intent(in) :: vertices_dot(3, 4), level_values_dot(4)
        real(dp), intent(out) :: positive_volume_dot
        real(dp), intent(out) :: positive_centroid_dot(3)
        real(dp), intent(out) :: negative_volume_dot
        real(dp), intent(out) :: negative_centroid_dot(3)
        real(dp), intent(out) :: interface_area_dot, normal_dot(3)
        integer, intent(out) :: status

        real(dp) :: positive_volume, positive_centroid(3)
        real(dp) :: negative_volume, negative_centroid(3)
        real(dp) :: interface_area, normal(3)
        real(dp) :: interface_points(3, 4), interface_points_dot(3, 4)
        logical :: has_positive, has_negative
        integer :: interface_count, interface_status

        positive_volume_dot = 0.0_dp
        positive_centroid_dot = 0.0_dp
        negative_volume_dot = 0.0_dp
        negative_centroid_dot = 0.0_dp
        interface_area_dot = 0.0_dp
        normal_dot = 0.0_dp
        status = 1
        if (any(.not. ieee_is_finite(vertices)) .or. &
            any(.not. ieee_is_finite(level_values)) .or. &
            any(.not. ieee_is_finite(vertices_dot)) .or. &
            any(.not. ieee_is_finite(level_values_dot))) return
        if (any(abs(level_values) <= topology_tolerance)) return
        call evaluate_level_set_tetra_cut_quadrature_3d( &
            vertices, level_values, positive_volume, positive_centroid, &
            negative_volume, negative_centroid, interface_area, normal, status)
        if (status /= 0) return

        has_positive = any(level_values > 0.0_dp)
        has_negative = any(level_values < 0.0_dp)
        interface_count = 0
        interface_points = 0.0_dp
        interface_points_dot = 0.0_dp
        if (has_positive .and. has_negative) then
            call evaluate_tetra_interface_geometry_jvp( &
                vertices, level_values, vertices_dot, level_values_dot, &
                interface_points, interface_points_dot, interface_count, &
                interface_area, interface_area_dot, normal, normal_dot, &
                interface_status)
            if (interface_status /= 0) then
                status = interface_status
                return
            end if
        end if

        call accumulate_tetra_side_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, .true., &
            interface_points, interface_points_dot, interface_count, &
            positive_volume, positive_centroid, positive_volume_dot, &
            positive_centroid_dot, status)
        if (status /= 0) return
        call accumulate_tetra_side_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, .false., &
            interface_points, interface_points_dot, interface_count, &
            negative_volume, negative_centroid, negative_volume_dot, &
            negative_centroid_dot, status)
    end subroutine evaluate_level_set_tetra_cut_quadrature_3d_jvp

    subroutine evaluate_level_set_tetra_cut_moments_3d( &
            vertices, level_values, positive_volume, positive_centroid, &
            positive_second_moment, negative_volume, negative_centroid, &
            negative_second_moment, interface_area, normal, status)
        !! Return exact degree-two raw moments for a linear tetrahedral cut.
        !!
        !! The tensor entries are physical integrals
        !! \\(M_{ab}=\\int_{\\Omega^\\pm}x_a x_b\\,dV\\).  The same oriented
        !! face/interface tetrahedral fans as the degree-one primitive are
        !! used, so positive and negative moments close to the parent moments.
        real(dp), intent(in) :: vertices(3, 4), level_values(4)
        real(dp), intent(out) :: positive_volume, positive_centroid(3)
        real(dp), intent(out) :: positive_second_moment(3, 3)
        real(dp), intent(out) :: negative_volume, negative_centroid(3)
        real(dp), intent(out) :: negative_second_moment(3, 3)
        real(dp), intent(out) :: interface_area, normal(3)
        integer, intent(out) :: status

        real(dp) :: interface_points(3, 4)
        logical :: has_positive, has_negative
        integer :: interface_count, interface_status

        positive_second_moment = 0.0_dp
        negative_second_moment = 0.0_dp
        call evaluate_level_set_tetra_cut_quadrature_3d( &
            vertices, level_values, positive_volume, positive_centroid, &
            negative_volume, negative_centroid, interface_area, normal, status)
        if (status /= 0) return

        interface_points = 0.0_dp
        interface_count = 0
        has_positive = any(level_values > 0.0_dp)
        has_negative = any(level_values < 0.0_dp)
        if (has_positive .and. has_negative) then
            call evaluate_level_set_tetra_interface_3d( &
                vertices, level_values, interface_points, interface_count, &
                interface_area, normal, interface_status)
            if (interface_status /= 0) then
                status = interface_status
                return
            end if
        end if
        call accumulate_tetra_side( &
            vertices, level_values, .true., interface_points, interface_count, &
            positive_volume, positive_centroid, status, positive_second_moment)
        if (status /= 0) return
        call accumulate_tetra_side( &
            vertices, level_values, .false., interface_points, interface_count, &
            negative_volume, negative_centroid, status, negative_second_moment)
    end subroutine evaluate_level_set_tetra_cut_moments_3d

    subroutine evaluate_level_set_tetra_cut_moments_3d_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, &
            positive_volume_dot, positive_centroid_dot, positive_second_moment_dot, &
            negative_volume_dot, negative_centroid_dot, negative_second_moment_dot, &
            interface_area_dot, normal_dot, status)
        !! Apply the fixed-topology JVP of degree-two tetrahedral moments.
        real(dp), intent(in) :: vertices(3, 4), level_values(4)
        real(dp), intent(in) :: vertices_dot(3, 4), level_values_dot(4)
        real(dp), intent(out) :: positive_volume_dot, positive_centroid_dot(3)
        real(dp), intent(out) :: positive_second_moment_dot(3, 3)
        real(dp), intent(out) :: negative_volume_dot, negative_centroid_dot(3)
        real(dp), intent(out) :: negative_second_moment_dot(3, 3)
        real(dp), intent(out) :: interface_area_dot, normal_dot(3)
        integer, intent(out) :: status

        real(dp) :: positive_volume, positive_centroid(3)
        real(dp) :: negative_volume, negative_centroid(3)
        real(dp) :: interface_area, normal(3)
        real(dp) :: interface_points(3, 4), interface_points_dot(3, 4)
        logical :: has_positive, has_negative
        integer :: interface_count, interface_status

        positive_second_moment_dot = 0.0_dp
        negative_second_moment_dot = 0.0_dp
        call evaluate_level_set_tetra_cut_quadrature_3d_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, &
            positive_volume_dot, positive_centroid_dot, negative_volume_dot, &
            negative_centroid_dot, interface_area_dot, normal_dot, status)
        if (status /= 0) return
        call evaluate_level_set_tetra_cut_quadrature_3d( &
            vertices, level_values, positive_volume, positive_centroid, &
            negative_volume, negative_centroid, interface_area, normal, status)
        if (status /= 0) return

        interface_points = 0.0_dp
        interface_points_dot = 0.0_dp
        interface_count = 0
        has_positive = any(level_values > 0.0_dp)
        has_negative = any(level_values < 0.0_dp)
        if (has_positive .and. has_negative) then
            call evaluate_tetra_interface_geometry_jvp( &
                vertices, level_values, vertices_dot, level_values_dot, &
                interface_points, interface_points_dot, interface_count, &
                interface_area, interface_area_dot, normal, normal_dot, &
                interface_status)
            if (interface_status /= 0) then
                status = interface_status
                return
            end if
        end if
        call accumulate_tetra_side_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, .true., &
            interface_points, interface_points_dot, interface_count, &
            positive_volume, positive_centroid, positive_volume_dot, &
            positive_centroid_dot, status, positive_second_moment_dot)
        if (status /= 0) return
        call accumulate_tetra_side_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, .false., &
            interface_points, interface_points_dot, interface_count, &
            negative_volume, negative_centroid, negative_volume_dot, &
            negative_centroid_dot, status, negative_second_moment_dot)
    end subroutine evaluate_level_set_tetra_cut_moments_3d_jvp

    subroutine evaluate_level_set_tetra_cut_third_moments_3d( &
            vertices, level_values, positive_third_moment, &
            negative_third_moment, status)
        !! Return exact degree-three raw moments for a linear tetrahedral cut.
        !!
        !! The symmetric tensor contains the physical integrals
        !! \(M_{abc}=\int_{\Omega_\pm}x_a x_b x_c\,dV\).  It is a separate
        !! entry point so existing degree-two callers retain their ABI.
        real(dp), intent(in) :: vertices(3, 4), level_values(4)
        real(dp), intent(out) :: positive_third_moment(3, 3, 3)
        real(dp), intent(out) :: negative_third_moment(3, 3, 3)
        integer, intent(out) :: status

        real(dp) :: positive_volume, positive_centroid(3)
        real(dp) :: negative_volume, negative_centroid(3)
        real(dp) :: interface_area, normal(3)
        real(dp) :: interface_points(3, 4)
        logical :: has_positive, has_negative
        integer :: interface_count, interface_status

        positive_third_moment = 0.0_dp
        negative_third_moment = 0.0_dp
        status = 1
        call evaluate_level_set_tetra_cut_quadrature_3d( &
            vertices, level_values, positive_volume, positive_centroid, &
            negative_volume, negative_centroid, interface_area, normal, status)
        if (status /= 0) return

        interface_points = 0.0_dp
        interface_count = 0
        has_positive = any(level_values > 0.0_dp)
        has_negative = any(level_values < 0.0_dp)
        if (has_positive .and. has_negative) then
            call evaluate_level_set_tetra_interface_3d( &
                vertices, level_values, interface_points, interface_count, &
                interface_area, normal, interface_status)
            if (interface_status /= 0) then
                status = interface_status
                return
            end if
        end if
        call accumulate_tetra_side( &
            vertices, level_values, .true., interface_points, interface_count, &
            positive_volume, positive_centroid, status, &
            third_moment=positive_third_moment)
        if (status /= 0) return
        call accumulate_tetra_side( &
            vertices, level_values, .false., interface_points, interface_count, &
            negative_volume, negative_centroid, status, &
            third_moment=negative_third_moment)
    end subroutine evaluate_level_set_tetra_cut_third_moments_3d

    subroutine evaluate_level_set_tetra_cut_third_moments_3d_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, &
            positive_third_moment_dot, negative_third_moment_dot, status)
        !! Apply the fixed-topology JVP of degree-three cut moments.
        real(dp), intent(in) :: vertices(3, 4), level_values(4)
        real(dp), intent(in) :: vertices_dot(3, 4), level_values_dot(4)
        real(dp), intent(out) :: positive_third_moment_dot(3, 3, 3)
        real(dp), intent(out) :: negative_third_moment_dot(3, 3, 3)
        integer, intent(out) :: status

        real(dp) :: positive_volume, positive_centroid(3)
        real(dp) :: negative_volume, negative_centroid(3)
        real(dp) :: positive_volume_dot, positive_centroid_dot(3)
        real(dp) :: negative_volume_dot, negative_centroid_dot(3)
        real(dp) :: interface_area, interface_area_dot, normal(3), normal_dot(3)
        real(dp) :: interface_points(3, 4), interface_points_dot(3, 4)
        logical :: has_positive, has_negative
        integer :: interface_count, interface_status

        positive_third_moment_dot = 0.0_dp
        negative_third_moment_dot = 0.0_dp
        status = 1
        call evaluate_level_set_tetra_cut_quadrature_3d_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, &
            positive_volume_dot, positive_centroid_dot, negative_volume_dot, &
            negative_centroid_dot, interface_area_dot, normal_dot, status)
        if (status /= 0) return
        call evaluate_level_set_tetra_cut_quadrature_3d( &
            vertices, level_values, positive_volume, positive_centroid, &
            negative_volume, negative_centroid, interface_area, normal, status)
        if (status /= 0) return

        interface_points = 0.0_dp
        interface_points_dot = 0.0_dp
        interface_count = 0
        has_positive = any(level_values > 0.0_dp)
        has_negative = any(level_values < 0.0_dp)
        if (has_positive .and. has_negative) then
            call evaluate_tetra_interface_geometry_jvp( &
                vertices, level_values, vertices_dot, level_values_dot, &
                interface_points, interface_points_dot, interface_count, &
                interface_area, interface_area_dot, normal, normal_dot, &
                interface_status)
            if (interface_status /= 0) then
                status = interface_status
                return
            end if
        end if
        call accumulate_tetra_side_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, .true., &
            interface_points, interface_points_dot, interface_count, &
            positive_volume, positive_centroid, positive_volume_dot, &
            positive_centroid_dot, status, &
            third_moment_dot=positive_third_moment_dot)
        if (status /= 0) return
        call accumulate_tetra_side_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, .false., &
            interface_points, interface_points_dot, interface_count, &
            negative_volume, negative_centroid, negative_volume_dot, &
            negative_centroid_dot, status, &
            third_moment_dot=negative_third_moment_dot)
    end subroutine evaluate_level_set_tetra_cut_third_moments_3d_jvp

    subroutine evaluate_level_set_tetra_cut_fourth_moments_3d( &
            vertices, level_values, positive_fourth_moment, &
            negative_fourth_moment, status)
        !! Return exact degree-four raw moments for a linear tetrahedral cut.
        real(dp), intent(in) :: vertices(3, 4), level_values(4)
        real(dp), intent(out) :: positive_fourth_moment(3, 3, 3, 3)
        real(dp), intent(out) :: negative_fourth_moment(3, 3, 3, 3)
        integer, intent(out) :: status

        real(dp) :: positive_volume, positive_centroid(3)
        real(dp) :: negative_volume, negative_centroid(3)
        real(dp) :: interface_area, normal(3)
        real(dp) :: interface_points(3, 4)
        logical :: has_positive, has_negative
        integer :: interface_count, interface_status

        positive_fourth_moment = 0.0_dp
        negative_fourth_moment = 0.0_dp
        status = 1
        call evaluate_level_set_tetra_cut_quadrature_3d( &
            vertices, level_values, positive_volume, positive_centroid, &
            negative_volume, negative_centroid, interface_area, normal, status)
        if (status /= 0) return

        interface_points = 0.0_dp
        interface_count = 0
        has_positive = any(level_values > 0.0_dp)
        has_negative = any(level_values < 0.0_dp)
        if (has_positive .and. has_negative) then
            call evaluate_level_set_tetra_interface_3d( &
                vertices, level_values, interface_points, interface_count, &
                interface_area, normal, interface_status)
            if (interface_status /= 0) then
                status = interface_status
                return
            end if
        end if
        call accumulate_tetra_side( &
            vertices, level_values, .true., interface_points, interface_count, &
            positive_volume, positive_centroid, status, &
            fourth_moment=positive_fourth_moment)
        if (status /= 0) return
        call accumulate_tetra_side( &
            vertices, level_values, .false., interface_points, interface_count, &
            negative_volume, negative_centroid, status, &
            fourth_moment=negative_fourth_moment)
    end subroutine evaluate_level_set_tetra_cut_fourth_moments_3d

    subroutine evaluate_level_set_tetra_cut_fourth_moments_3d_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, &
            positive_fourth_moment_dot, negative_fourth_moment_dot, status)
        !! Apply the fixed-topology JVP of degree-four tetrahedral cut moments.
        real(dp), intent(in) :: vertices(3, 4), level_values(4)
        real(dp), intent(in) :: vertices_dot(3, 4), level_values_dot(4)
        real(dp), intent(out) :: positive_fourth_moment_dot(3, 3, 3, 3)
        real(dp), intent(out) :: negative_fourth_moment_dot(3, 3, 3, 3)
        integer, intent(out) :: status

        real(dp) :: positive_volume, positive_centroid(3)
        real(dp) :: negative_volume, negative_centroid(3)
        real(dp) :: positive_volume_dot, positive_centroid_dot(3)
        real(dp) :: negative_volume_dot, negative_centroid_dot(3)
        real(dp) :: interface_area, interface_area_dot, normal(3), normal_dot(3)
        real(dp) :: interface_points(3, 4), interface_points_dot(3, 4)
        logical :: has_positive, has_negative
        integer :: interface_count, interface_status

        positive_fourth_moment_dot = 0.0_dp
        negative_fourth_moment_dot = 0.0_dp
        status = 1
        call evaluate_level_set_tetra_cut_quadrature_3d_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, &
            positive_volume_dot, positive_centroid_dot, negative_volume_dot, &
            negative_centroid_dot, interface_area_dot, normal_dot, status)
        if (status /= 0) return
        call evaluate_level_set_tetra_cut_quadrature_3d( &
            vertices, level_values, positive_volume, positive_centroid, &
            negative_volume, negative_centroid, interface_area, normal, status)
        if (status /= 0) return

        interface_points = 0.0_dp
        interface_points_dot = 0.0_dp
        interface_count = 0
        has_positive = any(level_values > 0.0_dp)
        has_negative = any(level_values < 0.0_dp)
        if (has_positive .and. has_negative) then
            call evaluate_tetra_interface_geometry_jvp( &
                vertices, level_values, vertices_dot, level_values_dot, &
                interface_points, interface_points_dot, interface_count, &
                interface_area, interface_area_dot, normal, normal_dot, &
                interface_status)
            if (interface_status /= 0) then
                status = interface_status
                return
            end if
        end if
        call accumulate_tetra_side_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, .true., &
            interface_points, interface_points_dot, interface_count, &
            positive_volume, positive_centroid, positive_volume_dot, &
            positive_centroid_dot, status, &
            fourth_moment_dot=positive_fourth_moment_dot)
        if (status /= 0) return
        call accumulate_tetra_side_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, .false., &
            interface_points, interface_points_dot, interface_count, &
            negative_volume, negative_centroid, negative_volume_dot, &
            negative_centroid_dot, status, &
            fourth_moment_dot=negative_fourth_moment_dot)
    end subroutine evaluate_level_set_tetra_cut_fourth_moments_3d_jvp

    subroutine accumulate_tetra_side( &
            vertices, level_values, keep_positive, interface_points, &
            interface_count, volume, centroid, status, second_moment, third_moment, &
            fourth_moment)
        real(dp), intent(in) :: vertices(3, 4), level_values(4)
        logical, intent(in) :: keep_positive
        real(dp), intent(in) :: interface_points(3, 4)
        integer, intent(in) :: interface_count
        real(dp), intent(out) :: volume, centroid(3)
        integer, intent(out) :: status
        real(dp), intent(out), optional :: second_moment(3, 3)
        real(dp), intent(out), optional :: third_moment(3, 3, 3)
        real(dp), intent(out), optional :: fourth_moment(3, 3, 3, 3)

        integer, parameter :: face_opposite(4) = [1, 2, 3, 4]
        integer :: face_vertices(3, 4), face, polygon_count
        integer :: indices(3), point
        real(dp) :: polygon(3, 6)
        real(dp) :: volume_accum, first_moment(3)
        real(dp) :: third_moment_accum(3, 3, 3)
        real(dp) :: fourth_moment_accum(3, 3, 3, 3)

        face_vertices(:, 1) = [2, 3, 4]
        face_vertices(:, 2) = [1, 4, 3]
        face_vertices(:, 3) = [1, 2, 4]
        face_vertices(:, 4) = [1, 3, 2]
        volume = 0.0_dp
        centroid = 0.0_dp
        if (present(second_moment)) second_moment = 0.0_dp
        if (present(third_moment)) third_moment = 0.0_dp
        if (present(fourth_moment)) fourth_moment = 0.0_dp
        status = 1
        volume_accum = 0.0_dp
        first_moment = 0.0_dp
        third_moment_accum = 0.0_dp
        fourth_moment_accum = 0.0_dp
        do face = 1, 4
            indices = face_vertices(:, face)
            call orient_tetra_face(vertices, face_opposite(face), indices)
            call clip_tetra_face( &
                vertices, level_values, indices, keep_positive, polygon, &
                polygon_count)
            if (polygon_count >= 3) then
                if (present(fourth_moment)) then
                    call accumulate_oriented_polygon( &
                        polygon, polygon_count, volume_accum, first_moment, &
                        second_moment, third_moment_accum, fourth_moment_accum)
                else
                    call accumulate_oriented_polygon( &
                        polygon, polygon_count, volume_accum, first_moment, &
                        second_moment, third_moment_accum)
                end if
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
            if (present(fourth_moment)) then
                call accumulate_oriented_polygon( &
                    polygon, interface_count, volume_accum, first_moment, &
                    second_moment, third_moment_accum, fourth_moment_accum)
            else
                call accumulate_oriented_polygon( &
                    polygon, interface_count, volume_accum, first_moment, &
                    second_moment, third_moment_accum)
            end if
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
        if (present(third_moment)) third_moment = third_moment_accum
        if (present(fourth_moment)) fourth_moment = fourth_moment_accum
        status = 0
    end subroutine accumulate_tetra_side

    subroutine evaluate_tetra_interface_geometry_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, points, &
            points_dot, point_count, area, area_dot, normal, normal_dot, status)
        real(dp), intent(in) :: vertices(3, 4), level_values(4)
        real(dp), intent(in) :: vertices_dot(3, 4), level_values_dot(4)
        real(dp), intent(out) :: points(3, 4), points_dot(3, 4)
        integer, intent(out) :: point_count
        real(dp), intent(out) :: area, area_dot, normal(3), normal_dot(3)
        integer, intent(out) :: status

        integer, parameter :: edge_start(6) = [1, 1, 1, 2, 2, 3]
        integer, parameter :: edge_end(6) = [2, 3, 4, 3, 4, 4]
        integer :: edge, first, second, candidate_count, candidate
        integer :: point, next_point, best_candidate
        logical :: used(4)
        real(dp) :: candidate_points(3, 4), candidate_points_dot(3, 4)
        real(dp) :: first_value, second_value, fraction, fraction_dot
        real(dp) :: edge_vector(3), edge_vector_dot(3)
        real(dp) :: edge_point(3), edge_point_dot(3), distance_squared
        real(dp) :: best_distance, area_vector(3), area_vector_dot(3)
        real(dp) :: vector_norm, jacobian(3, 3), jacobian_dot(3, 3)
        real(dp) :: inverse(3, 3), inverse_dot(3, 3), determinant
        real(dp) :: gradient_reference(3), gradient_reference_dot(3)
        real(dp) :: gradient(3), gradient_dot(3), gradient_norm
        integer :: interface_status, inverse_status

        points = 0.0_dp
        points_dot = 0.0_dp
        point_count = 0
        area = 0.0_dp
        area_dot = 0.0_dp
        normal = 0.0_dp
        normal_dot = 0.0_dp
        status = 1
        if (any(.not. ieee_is_finite(vertices)) .or. &
            any(.not. ieee_is_finite(level_values)) .or. &
            any(.not. ieee_is_finite(vertices_dot)) .or. &
            any(.not. ieee_is_finite(level_values_dot))) return
        if (any(abs(level_values) <= topology_tolerance)) return
        call evaluate_level_set_tetra_interface_3d( &
            vertices, level_values, points, point_count, area, normal, &
            interface_status)
        if (interface_status /= 0) return

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
        jacobian_dot(:, 1) = vertices_dot(:, 2) - vertices_dot(:, 1)
        jacobian_dot(:, 2) = vertices_dot(:, 3) - vertices_dot(:, 1)
        jacobian_dot(:, 3) = vertices_dot(:, 4) - vertices_dot(:, 1)
        determinant = det3(jacobian)
        if (abs(determinant) <= topology_tolerance* &
            max(1.0_dp, maxval(abs(jacobian))**3)) return
        call inv3(jacobian, inverse, inverse_status)
        if (inverse_status /= 0) return
        inverse_dot = -matmul(inverse, matmul(jacobian_dot, inverse))
        gradient_reference = level_values(2:4) - level_values(1)
        gradient_reference_dot = level_values_dot(2:4) - level_values_dot(1)
        gradient = matmul(transpose(inverse), gradient_reference)
        gradient_dot = matmul(transpose(inverse_dot), gradient_reference) + &
            matmul(transpose(inverse), gradient_reference_dot)
        gradient_norm = sqrt(dot_product(gradient, gradient))
        if (gradient_norm <= topology_tolerance) return
        normal_dot = (gradient_dot - normal*dot_product(normal, gradient_dot))/ &
            gradient_norm

        candidate_count = 0
        do edge = 1, 6
            first = edge_start(edge)
            second = edge_end(edge)
            first_value = level_values(first)
            second_value = level_values(second)
            if (first_value*second_value < 0.0_dp) then
                fraction = first_value/(first_value - second_value)
                fraction_dot = (level_values_dot(first)*(first_value - second_value) - &
                    first_value*(level_values_dot(first) - level_values_dot(second)))/ &
                    (first_value - second_value)**2
                edge_vector = vertices(:, second) - vertices(:, first)
                edge_vector_dot = vertices_dot(:, second) - vertices_dot(:, first)
                edge_point = vertices(:, first) + fraction*edge_vector
                edge_point_dot = vertices_dot(:, first) + fraction_dot*edge_vector + &
                    fraction*edge_vector_dot
                candidate_count = candidate_count + 1
                candidate_points(:, candidate_count) = edge_point
                candidate_points_dot(:, candidate_count) = edge_point_dot
            end if
        end do
        if (candidate_count /= point_count) then
            points_dot = 0.0_dp
            area_dot = 0.0_dp
            normal_dot = 0.0_dp
            return
        end if
        used = .false.
        do point = 1, point_count
            best_candidate = 0
            best_distance = huge(1.0_dp)
            do candidate = 1, candidate_count
                if (.not. used(candidate)) then
                    distance_squared = sum((points(:, point) - &
                        candidate_points(:, candidate))**2)
                    if (distance_squared < best_distance) then
                        best_distance = distance_squared
                        best_candidate = candidate
                    end if
                end if
            end do
            if (best_candidate == 0) then
                points_dot = 0.0_dp
                area_dot = 0.0_dp
                normal_dot = 0.0_dp
                return
            end if
            used(best_candidate) = .true.
            points_dot(:, point) = candidate_points_dot(:, best_candidate)
        end do

        area_vector = 0.0_dp
        area_vector_dot = 0.0_dp
        do point = 1, point_count
            next_point = 1 + mod(point, point_count)
            area_vector = area_vector + cross_product( &
                points(:, point), points(:, next_point))
            area_vector_dot = area_vector_dot + cross_product( &
                points_dot(:, point), points(:, next_point)) + cross_product( &
                points(:, point), points_dot(:, next_point))
        end do
        vector_norm = sqrt(dot_product(area_vector, area_vector))
        if (vector_norm <= topology_tolerance) return
        area_dot = 0.5_dp*dot_product(area_vector, area_vector_dot)/vector_norm
        status = 0
    end subroutine evaluate_tetra_interface_geometry_jvp

    subroutine accumulate_tetra_side_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, &
            keep_positive, interface_points, interface_points_dot, &
            interface_count, volume, centroid, volume_dot, centroid_dot, status, &
            second_moment_dot, third_moment_dot, fourth_moment_dot)
        real(dp), intent(in) :: vertices(3, 4), level_values(4)
        real(dp), intent(in) :: vertices_dot(3, 4), level_values_dot(4)
        logical, intent(in) :: keep_positive
        real(dp), intent(in) :: interface_points(3, 4), interface_points_dot(3, 4)
        integer, intent(in) :: interface_count
        real(dp), intent(in) :: volume, centroid(3)
        real(dp), intent(out) :: volume_dot, centroid_dot(3)
        integer, intent(out) :: status
        real(dp), intent(out), optional :: second_moment_dot(3, 3)
        real(dp), intent(out), optional :: third_moment_dot(3, 3, 3)
        real(dp), intent(out), optional :: fourth_moment_dot(3, 3, 3, 3)

        integer, parameter :: face_opposite(4) = [1, 2, 3, 4]
        integer :: face_vertices(3, 4), face, polygon_count
        integer :: indices(3), point
        real(dp) :: polygon(3, 6), polygon_dot(3, 6)
        real(dp) :: volume_accum, volume_accum_dot
        real(dp) :: first_moment_dot(3)
        real(dp) :: third_moment_dot_accum(3, 3, 3)
        real(dp) :: fourth_moment_dot_accum(3, 3, 3, 3)

        face_vertices(:, 1) = [2, 3, 4]
        face_vertices(:, 2) = [1, 4, 3]
        face_vertices(:, 3) = [1, 2, 4]
        face_vertices(:, 4) = [1, 3, 2]
        volume_dot = 0.0_dp
        centroid_dot = 0.0_dp
        if (present(second_moment_dot)) second_moment_dot = 0.0_dp
        if (present(third_moment_dot)) third_moment_dot = 0.0_dp
        if (present(fourth_moment_dot)) fourth_moment_dot = 0.0_dp
        status = 1
        volume_accum = 0.0_dp
        volume_accum_dot = 0.0_dp
        first_moment_dot = 0.0_dp
        third_moment_dot_accum = 0.0_dp
        fourth_moment_dot_accum = 0.0_dp
        do face = 1, 4
            indices = face_vertices(:, face)
            call orient_tetra_face(vertices, face_opposite(face), indices)
            call clip_tetra_face_jvp( &
                vertices, level_values, vertices_dot, level_values_dot, indices, &
                keep_positive, polygon, polygon_dot, polygon_count)
            if (polygon_count >= 3) then
                if (present(fourth_moment_dot)) then
                    call accumulate_oriented_polygon_jvp( &
                        polygon, polygon_dot, polygon_count, volume_accum, &
                        volume_accum_dot, first_moment_dot, second_moment_dot, &
                        third_moment_dot_accum, fourth_moment_dot_accum)
                else
                    call accumulate_oriented_polygon_jvp( &
                        polygon, polygon_dot, polygon_count, volume_accum, &
                        volume_accum_dot, first_moment_dot, second_moment_dot, &
                        third_moment_dot_accum)
                end if
            end if
        end do

        if (interface_count > 0) then
            if (keep_positive) then
                do point = 1, interface_count
                    polygon(:, point) = &
                        interface_points(:, interface_count - point + 1)
                    polygon_dot(:, point) = &
                        interface_points_dot(:, interface_count - point + 1)
                end do
            else
                polygon(:, :interface_count) = &
                    interface_points(:, :interface_count)
                polygon_dot(:, :interface_count) = &
                    interface_points_dot(:, :interface_count)
            end if
            if (present(fourth_moment_dot)) then
                call accumulate_oriented_polygon_jvp( &
                    polygon, polygon_dot, interface_count, volume_accum, &
                    volume_accum_dot, first_moment_dot, second_moment_dot, &
                    third_moment_dot_accum, fourth_moment_dot_accum)
            else
                call accumulate_oriented_polygon_jvp( &
                    polygon, polygon_dot, interface_count, volume_accum, &
                    volume_accum_dot, first_moment_dot, second_moment_dot, &
                    third_moment_dot_accum)
            end if
        end if

        if (volume <= topology_tolerance) then
            if (abs(volume_accum) > 8.0_dp*topology_tolerance) return
            volume_dot = 0.0_dp
            centroid_dot = 0.0_dp
            status = 0
            return
        end if
        volume_dot = volume_accum_dot
        centroid_dot = (first_moment_dot - volume_dot*centroid)/volume
        if (present(third_moment_dot)) third_moment_dot = third_moment_dot_accum
        if (present(fourth_moment_dot)) fourth_moment_dot = fourth_moment_dot_accum
        status = 0
    end subroutine accumulate_tetra_side_jvp

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

    subroutine clip_tetra_face_jvp( &
            vertices, level_values, vertices_dot, level_values_dot, indices, &
            keep_positive, polygon, polygon_dot, point_count)
        real(dp), intent(in) :: vertices(3, 4), level_values(4)
        real(dp), intent(in) :: vertices_dot(3, 4), level_values_dot(4)
        integer, intent(in) :: indices(3)
        logical, intent(in) :: keep_positive
        real(dp), intent(out) :: polygon(3, 6), polygon_dot(3, 6)
        integer, intent(out) :: point_count

        integer :: current, previous, vertex
        real(dp) :: current_value, previous_value, fraction, fraction_dot
        real(dp) :: current_value_dot, previous_value_dot
        real(dp) :: edge_vector(3), edge_vector_dot(3)
        logical :: current_inside, previous_inside

        polygon = 0.0_dp
        polygon_dot = 0.0_dp
        point_count = 0
        previous = indices(3)
        previous_value = level_values(previous)
        previous_value_dot = level_values_dot(previous)
        previous_inside = side_contains(previous_value, keep_positive)
        do vertex = 1, 3
            current = indices(vertex)
            current_value = level_values(current)
            current_value_dot = level_values_dot(current)
            current_inside = side_contains(current_value, keep_positive)
            if (current_inside .neqv. previous_inside) then
                fraction = previous_value/(previous_value - current_value)
                fraction_dot = (previous_value_dot*(previous_value - current_value) - &
                    previous_value*(previous_value_dot - current_value_dot))/ &
                    (previous_value - current_value)**2
                edge_vector = vertices(:, current) - vertices(:, previous)
                edge_vector_dot = vertices_dot(:, current) - vertices_dot(:, previous)
                point_count = point_count + 1
                polygon(:, point_count) = vertices(:, previous) + fraction*edge_vector
                polygon_dot(:, point_count) = vertices_dot(:, previous) + &
                    fraction_dot*edge_vector + fraction*edge_vector_dot
            end if
            if (current_inside) then
                point_count = point_count + 1
                polygon(:, point_count) = vertices(:, current)
                polygon_dot(:, point_count) = vertices_dot(:, current)
            end if
            previous = current
            previous_value = current_value
            previous_value_dot = current_value_dot
            previous_inside = current_inside
        end do
    end subroutine clip_tetra_face_jvp

    pure logical function side_contains(value, keep_positive)
        real(dp), intent(in) :: value
        logical, intent(in) :: keep_positive

        side_contains = merge(value > 0.0_dp, value < 0.0_dp, keep_positive)
    end function side_contains

    subroutine accumulate_oriented_polygon( &
            polygon, point_count, volume, first_moment, second_moment, third_moment, &
            fourth_moment)
        real(dp), intent(in) :: polygon(3, 6)
        integer, intent(in) :: point_count
        real(dp), intent(inout) :: volume, first_moment(3)
        real(dp), intent(inout), optional :: second_moment(3, 3)
        real(dp), intent(inout), optional :: third_moment(3, 3, 3)
        real(dp), intent(inout), optional :: fourth_moment(3, 3, 3, 3)

        integer :: point
        real(dp) :: signed_tetra_volume, first_point(3)
        real(dp) :: second_point(3), third_point(3)
        real(dp) :: pair_moment(3, 3)
        real(dp) :: triple_moment(3, 3, 3)
        real(dp) :: fourth_moment_local(3, 3, 3, 3)

        if (present(second_moment)) pair_moment = 0.0_dp
        if (present(third_moment)) triple_moment = 0.0_dp
        if (present(fourth_moment)) fourth_moment_local = 0.0_dp

        first_point = polygon(:, 1)
        do point = 2, point_count - 1
            second_point = polygon(:, point)
            third_point = polygon(:, point + 1)
            signed_tetra_volume = dot_product(first_point, cross_product( &
                second_point, third_point))/6.0_dp
            volume = volume + signed_tetra_volume
            first_moment = first_moment + signed_tetra_volume*( &
                first_point + second_point + third_point)/4.0_dp
            if (present(second_moment)) then
                pair_moment = (outer_product(first_point, first_point) + &
                    outer_product(second_point, second_point) + &
                    outer_product(third_point, third_point))/10.0_dp + &
                    (outer_product(first_point, second_point) + &
                    outer_product(second_point, first_point) + &
                    outer_product(first_point, third_point) + &
                    outer_product(third_point, first_point) + &
                    outer_product(second_point, third_point) + &
                    outer_product(third_point, second_point))/20.0_dp
                second_moment = second_moment + signed_tetra_volume*pair_moment
            end if
            if (present(third_moment)) then
                triple_moment = tetra_third_moment( &
                    first_point, second_point, third_point)
                third_moment = third_moment + signed_tetra_volume*triple_moment
            end if
            if (present(fourth_moment)) then
                fourth_moment_local = tetra_fourth_moment( &
                    first_point, second_point, third_point)
                fourth_moment = fourth_moment + &
                    signed_tetra_volume*fourth_moment_local
            end if
        end do
    end subroutine accumulate_oriented_polygon

    subroutine accumulate_oriented_polygon_jvp( &
            polygon, polygon_dot, point_count, volume, volume_dot, first_moment_dot, &
            second_moment_dot, third_moment_dot, fourth_moment_dot)
        real(dp), intent(in) :: polygon(3, 6), polygon_dot(3, 6)
        integer, intent(in) :: point_count
        real(dp), intent(inout) :: volume, volume_dot, first_moment_dot(3)
        real(dp), intent(inout), optional :: second_moment_dot(3, 3)
        real(dp), intent(inout), optional :: third_moment_dot(3, 3, 3)
        real(dp), intent(inout), optional :: fourth_moment_dot(3, 3, 3, 3)

        integer :: point
        real(dp) :: first_point(3), first_point_dot(3)
        real(dp) :: second_point(3), second_point_dot(3)
        real(dp) :: third_point(3), third_point_dot(3)
        real(dp) :: signed_tetra_volume, signed_tetra_volume_dot
        real(dp) :: pair_moment(3, 3), pair_moment_dot(3, 3)
        real(dp) :: triple_moment(3, 3, 3), triple_moment_dot(3, 3, 3)
        real(dp) :: fourth_moment_local(3, 3, 3, 3)
        real(dp) :: fourth_moment_local_dot(3, 3, 3, 3)

        if (present(second_moment_dot)) then
            pair_moment = 0.0_dp
            pair_moment_dot = 0.0_dp
        end if
        if (present(third_moment_dot)) then
            triple_moment = 0.0_dp
            triple_moment_dot = 0.0_dp
        end if
        if (present(fourth_moment_dot)) then
            fourth_moment_local = 0.0_dp
            fourth_moment_local_dot = 0.0_dp
        end if

        first_point = polygon(:, 1)
        first_point_dot = polygon_dot(:, 1)
        do point = 2, point_count - 1
            second_point = polygon(:, point)
            second_point_dot = polygon_dot(:, point)
            third_point = polygon(:, point + 1)
            third_point_dot = polygon_dot(:, point + 1)
            signed_tetra_volume = dot_product(first_point, cross_product( &
                second_point, third_point))/6.0_dp
            signed_tetra_volume_dot = ( &
                dot_product(first_point_dot, cross_product(second_point, third_point)) + &
                dot_product(first_point, cross_product(second_point_dot, third_point)) + &
                dot_product(first_point, cross_product(second_point, third_point_dot)))/6.0_dp
            volume = volume + signed_tetra_volume
            volume_dot = volume_dot + signed_tetra_volume_dot
            first_moment_dot = first_moment_dot + &
                (signed_tetra_volume_dot*(first_point + second_point + third_point) + &
                signed_tetra_volume*(first_point_dot + second_point_dot + &
                third_point_dot))/4.0_dp
            if (present(second_moment_dot)) then
                pair_moment = (outer_product(first_point, first_point) + &
                    outer_product(second_point, second_point) + &
                    outer_product(third_point, third_point))/10.0_dp + &
                    (outer_product(first_point, second_point) + &
                    outer_product(second_point, first_point) + &
                    outer_product(first_point, third_point) + &
                    outer_product(third_point, first_point) + &
                    outer_product(second_point, third_point) + &
                    outer_product(third_point, second_point))/20.0_dp
                pair_moment_dot = (outer_product(first_point_dot, first_point) + &
                    outer_product(first_point, first_point_dot) + &
                    outer_product(second_point_dot, second_point) + &
                    outer_product(second_point, second_point_dot) + &
                    outer_product(third_point_dot, third_point) + &
                    outer_product(third_point, third_point_dot))/10.0_dp + &
                    (outer_product(first_point_dot, second_point) + &
                    outer_product(first_point, second_point_dot) + &
                    outer_product(second_point_dot, first_point) + &
                    outer_product(second_point, first_point_dot) + &
                    outer_product(first_point_dot, third_point) + &
                    outer_product(first_point, third_point_dot) + &
                    outer_product(third_point_dot, first_point) + &
                    outer_product(third_point, first_point_dot) + &
                    outer_product(second_point_dot, third_point) + &
                    outer_product(second_point, third_point_dot) + &
                    outer_product(third_point_dot, second_point) + &
                    outer_product(third_point, second_point_dot))/20.0_dp
                second_moment_dot = second_moment_dot + &
                    signed_tetra_volume_dot*pair_moment + &
                    signed_tetra_volume*pair_moment_dot
            end if
            if (present(third_moment_dot)) then
                triple_moment = tetra_third_moment( &
                    first_point, second_point, third_point)
                triple_moment_dot = tetra_third_moment_jvp( &
                    first_point, second_point, third_point, first_point_dot, &
                    second_point_dot, third_point_dot)
                third_moment_dot = third_moment_dot + &
                    signed_tetra_volume_dot*triple_moment + &
                    signed_tetra_volume*triple_moment_dot
            end if
            if (present(fourth_moment_dot)) then
                fourth_moment_local = tetra_fourth_moment( &
                    first_point, second_point, third_point)
                fourth_moment_local_dot = tetra_fourth_moment_jvp( &
                    first_point, second_point, third_point, first_point_dot, &
                    second_point_dot, third_point_dot)
                fourth_moment_dot = fourth_moment_dot + &
                    signed_tetra_volume_dot*fourth_moment_local + &
                    signed_tetra_volume*fourth_moment_local_dot
            end if
        end do
    end subroutine accumulate_oriented_polygon_jvp

    pure function tetra_fourth_moment(first, second, third) result(moment)
        !! Exact rank-four moments of the tetrahedron (0, first, second, third).
        real(dp), intent(in) :: first(3), second(3), third(3)
        real(dp) :: moment(3, 3, 3, 3)

        real(dp) :: vectors(3, 3), product, coefficient
        integer :: first_index, second_index, third_index, fourth_index
        integer :: assignment_first, assignment_second, assignment_third
        integer :: assignment_fourth, counts(3)

        vectors(:, 1) = first
        vectors(:, 2) = second
        vectors(:, 3) = third
        moment = 0.0_dp
        do first_index = 1, 3
            do second_index = 1, 3
                do third_index = 1, 3
                    do fourth_index = 1, 3
                        do assignment_first = 1, 3
                            do assignment_second = 1, 3
                                do assignment_third = 1, 3
                                    do assignment_fourth = 1, 3
                                        counts = 0
                                        counts(assignment_first) = &
                                            counts(assignment_first) + 1
                                        counts(assignment_second) = &
                                            counts(assignment_second) + 1
                                        counts(assignment_third) = &
                                            counts(assignment_third) + 1
                                        counts(assignment_fourth) = &
                                            counts(assignment_fourth) + 1
                                        product = vectors(first_index, assignment_first)* &
                                            vectors(second_index, assignment_second)* &
                                            vectors(third_index, assignment_third)* &
                                            vectors(fourth_index, assignment_fourth)
                                        coefficient = 6.0_dp*real( &
                                            factorial_integer(counts(1))* &
                                            factorial_integer(counts(2))* &
                                            factorial_integer(counts(3)), dp)/5040.0_dp
                                        moment(first_index, second_index, third_index, &
                                            fourth_index) = &
                                            moment(first_index, second_index, third_index, &
                                            fourth_index) + coefficient*product
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end function tetra_fourth_moment

    pure function tetra_fourth_moment_jvp( &
            first, second, third, first_dot, second_dot, third_dot) &
            result(moment_dot)
        real(dp), intent(in) :: first(3), second(3), third(3)
        real(dp), intent(in) :: first_dot(3), second_dot(3), third_dot(3)
        real(dp) :: moment_dot(3, 3, 3, 3)

        real(dp) :: vectors(3, 3), vectors_dot(3, 3)
        real(dp) :: product, product_dot, coefficient
        integer :: first_index, second_index, third_index, fourth_index
        integer :: assignment_first, assignment_second, assignment_third
        integer :: assignment_fourth, counts(3)

        vectors(:, 1) = first
        vectors(:, 2) = second
        vectors(:, 3) = third
        vectors_dot(:, 1) = first_dot
        vectors_dot(:, 2) = second_dot
        vectors_dot(:, 3) = third_dot
        moment_dot = 0.0_dp
        do first_index = 1, 3
            do second_index = 1, 3
                do third_index = 1, 3
                    do fourth_index = 1, 3
                        do assignment_first = 1, 3
                            do assignment_second = 1, 3
                                do assignment_third = 1, 3
                                    do assignment_fourth = 1, 3
                                        counts = 0
                                        counts(assignment_first) = &
                                            counts(assignment_first) + 1
                                        counts(assignment_second) = &
                                            counts(assignment_second) + 1
                                        counts(assignment_third) = &
                                            counts(assignment_third) + 1
                                        counts(assignment_fourth) = &
                                            counts(assignment_fourth) + 1
                                        product = vectors(first_index, assignment_first)* &
                                            vectors(second_index, assignment_second)* &
                                            vectors(third_index, assignment_third)* &
                                            vectors(fourth_index, assignment_fourth)
                                        product_dot = &
                                            vectors_dot(first_index, assignment_first)* &
                                            vectors(second_index, assignment_second)* &
                                            vectors(third_index, assignment_third)* &
                                            vectors(fourth_index, assignment_fourth) + &
                                            vectors(first_index, assignment_first)* &
                                            vectors_dot(second_index, assignment_second)* &
                                            vectors(third_index, assignment_third)* &
                                            vectors(fourth_index, assignment_fourth) + &
                                            vectors(first_index, assignment_first)* &
                                            vectors(second_index, assignment_second)* &
                                            vectors_dot(third_index, assignment_third)* &
                                            vectors(fourth_index, assignment_fourth) + &
                                            vectors(first_index, assignment_first)* &
                                            vectors(second_index, assignment_second)* &
                                            vectors(third_index, assignment_third)* &
                                            vectors_dot(fourth_index, assignment_fourth)
                                        coefficient = 6.0_dp*real( &
                                            factorial_integer(counts(1))* &
                                            factorial_integer(counts(2))* &
                                            factorial_integer(counts(3)), dp)/5040.0_dp
                                        moment_dot(first_index, second_index, third_index, &
                                            fourth_index) = &
                                            moment_dot(first_index, second_index, third_index, &
                                            fourth_index) + coefficient*product_dot
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end function tetra_fourth_moment_jvp

    pure integer function factorial_integer(number) result(value)
        integer, intent(in) :: number
        integer :: factor

        value = 1
        do factor = 2, number
            value = value*factor
        end do
    end function factorial_integer

    pure function tetra_third_moment(first, second, third) result(moment)
        real(dp), intent(in) :: first(3), second(3), third(3)
        real(dp) :: moment(3, 3, 3)

        moment = ( &
            symmetric_triple_product(first, first, first) + &
            symmetric_triple_product(second, second, second) + &
            symmetric_triple_product(third, third, third) + &
            symmetric_triple_product(first, first, second) + &
            symmetric_triple_product(first, first, third) + &
            symmetric_triple_product(second, second, first) + &
            symmetric_triple_product(second, second, third) + &
            symmetric_triple_product(third, third, first) + &
            symmetric_triple_product(third, third, second) + &
            symmetric_triple_product(first, second, third))/120.0_dp
    end function tetra_third_moment

    pure function tetra_third_moment_jvp( &
            first, second, third, first_dot, second_dot, third_dot) result(moment_dot)
        real(dp), intent(in) :: first(3), second(3), third(3)
        real(dp), intent(in) :: first_dot(3), second_dot(3), third_dot(3)
        real(dp) :: moment_dot(3, 3, 3)

        moment_dot = ( &
            symmetric_triple_product_jvp( &
            first, first, first, first_dot, first_dot, first_dot) + &
            symmetric_triple_product_jvp( &
            second, second, second, second_dot, second_dot, second_dot) + &
            symmetric_triple_product_jvp( &
            third, third, third, third_dot, third_dot, third_dot) + &
            symmetric_triple_product_jvp( &
            first, first, second, first_dot, first_dot, second_dot) + &
            symmetric_triple_product_jvp( &
            first, first, third, first_dot, first_dot, third_dot) + &
            symmetric_triple_product_jvp( &
            second, second, first, second_dot, second_dot, first_dot) + &
            symmetric_triple_product_jvp( &
            second, second, third, second_dot, second_dot, third_dot) + &
            symmetric_triple_product_jvp( &
            third, third, first, third_dot, third_dot, first_dot) + &
            symmetric_triple_product_jvp( &
            third, third, second, third_dot, third_dot, second_dot) + &
            symmetric_triple_product_jvp( &
            first, second, third, first_dot, second_dot, third_dot))/120.0_dp
    end function tetra_third_moment_jvp

    pure function symmetric_triple_product(first, second, third) result(product)
        real(dp), intent(in) :: first(3), second(3), third(3)
        real(dp) :: product(3, 3, 3)

        product = outer_product3(first, second, third) + &
            outer_product3(first, third, second) + &
            outer_product3(second, first, third) + &
            outer_product3(second, third, first) + &
            outer_product3(third, first, second) + &
            outer_product3(third, second, first)
    end function symmetric_triple_product

    pure function symmetric_triple_product_jvp( &
            first, second, third, first_dot, second_dot, third_dot) result(product_dot)
        real(dp), intent(in) :: first(3), second(3), third(3)
        real(dp), intent(in) :: first_dot(3), second_dot(3), third_dot(3)
        real(dp) :: product_dot(3, 3, 3)

        product_dot = outer_product3_jvp( &
            first, second, third, first_dot, second_dot, third_dot) + &
            outer_product3_jvp( &
            first, third, second, first_dot, third_dot, second_dot) + &
            outer_product3_jvp( &
            second, first, third, second_dot, first_dot, third_dot) + &
            outer_product3_jvp( &
            second, third, first, second_dot, third_dot, first_dot) + &
            outer_product3_jvp( &
            third, first, second, third_dot, first_dot, second_dot) + &
            outer_product3_jvp( &
            third, second, first, third_dot, second_dot, first_dot)
    end function symmetric_triple_product_jvp

    pure function outer_product3(first, second, third) result(product)
        real(dp), intent(in) :: first(3), second(3), third(3)
        real(dp) :: product(3, 3, 3)
        integer :: first_index, second_index, third_index

        do first_index = 1, 3
            do second_index = 1, 3
                do third_index = 1, 3
                    product(first_index, second_index, third_index) = &
                        first(first_index)*second(second_index)*third(third_index)
                end do
            end do
        end do
    end function outer_product3

    pure function outer_product3_jvp( &
            first, second, third, first_dot, second_dot, third_dot) result(product_dot)
        real(dp), intent(in) :: first(3), second(3), third(3)
        real(dp), intent(in) :: first_dot(3), second_dot(3), third_dot(3)
        real(dp) :: product_dot(3, 3, 3)

        product_dot = outer_product3(first_dot, second, third) + &
            outer_product3(first, second_dot, third) + &
            outer_product3(first, second, third_dot)
    end function outer_product3_jvp

    pure function outer_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3, 3)

        product(1, :) = first(1)*second
        product(2, :) = first(2)*second
        product(3, :) = first(3)*second
    end function outer_product

    pure function cross_product(first, second) result(cross)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: cross(3)

        cross = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_level_set_tetra_interface_3d
