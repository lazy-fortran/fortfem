module fortfem_fci_terminal_polyline_3d
    !! First material-surface event for a traced 3D FCI polyline.
    !!
    !! The polyline may contain several field-line integration segments.  The
    !! first transverse triangle hit is returned with its segment index,
    !! oriented normal, and physical connection length.  A valid no-hit path
    !! returns its final endpoint and total polyline length.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_fci_terminal_triangle_3d, only: &
        find_fci_first_hit_triangle_3d, find_fci_first_hit_triangle_3d_jvp
    use fortfem_kinds, only: dp
    implicit none

    private

    real(dp), parameter :: topology_tolerance = 64.0_dp*epsilon(1.0_dp)

    public :: find_fci_first_hit_polyline_3d
    public :: find_fci_first_hit_polyline_3d_jvp

contains

    subroutine find_fci_first_hit_polyline_3d( &
            polyline, surface_vertices, surface_triangles, hit_point, &
            connection_length, hit_segment, hit_triangle, normal, status)
        !! Find the first material triangle hit along a traced polyline.
        real(dp), intent(in) :: polyline(:, :), surface_vertices(:, :)
        integer, intent(in) :: surface_triangles(:, :)
        real(dp), intent(out) :: hit_point(3), connection_length, normal(3)
        integer, intent(out) :: hit_segment, hit_triangle, status

        integer :: segment, point_count, segment_status
        real(dp) :: segment_vector(3), segment_length, cumulative_length
        real(dp) :: local_point(3), local_parameter, local_normal(3)
        integer :: local_triangle

        hit_point = 0.0_dp
        connection_length = 0.0_dp
        hit_segment = 0
        hit_triangle = 0
        normal = 0.0_dp
        status = 1
        if (.not. valid_polyline( &
            polyline, surface_vertices, surface_triangles)) return
        point_count = size(polyline, 2)
        cumulative_length = 0.0_dp
        do segment = 1, point_count - 1
            segment_vector = polyline(:, segment + 1) - polyline(:, segment)
            segment_length = vector_norm(segment_vector)
            call find_fci_first_hit_triangle_3d( &
                polyline(:, segment), polyline(:, segment + 1), surface_vertices, &
                surface_triangles, local_point, local_parameter, local_triangle, &
                local_normal, segment_status)
            if (segment_status /= 0) return
            if (local_triangle /= 0) then
                hit_point = local_point
                connection_length = cumulative_length + local_parameter*segment_length
                hit_segment = segment
                hit_triangle = local_triangle
                normal = local_normal
                status = 0
                return
            end if
            cumulative_length = cumulative_length + segment_length
        end do
        hit_point = polyline(:, point_count)
        connection_length = cumulative_length
        status = 0
    end subroutine find_fci_first_hit_polyline_3d

    subroutine find_fci_first_hit_polyline_3d_jvp( &
            polyline, surface_vertices, surface_triangles, polyline_dot, &
            surface_vertices_dot, hit_point_dot, connection_length_dot, normal_dot, &
            status)
        !! Apply the fixed-event JVP of the first-hit polyline geometry.
        real(dp), intent(in) :: polyline(:, :), surface_vertices(:, :)
        integer, intent(in) :: surface_triangles(:, :)
        real(dp), intent(in) :: polyline_dot(:, :), surface_vertices_dot(:, :)
        real(dp), intent(out) :: hit_point_dot(3), connection_length_dot
        real(dp), intent(out) :: normal_dot(3)
        integer, intent(out) :: status

        real(dp) :: hit_point(3), connection_length, normal(3)
        integer :: hit_segment, hit_triangle, primal_status
        integer :: segment, segment_status
        real(dp) :: segment_vector(3), segment_vector_dot(3), segment_length
        real(dp) :: segment_length_dot, cumulative_length_dot
        real(dp) :: local_point(3), local_parameter, local_normal(3)
        real(dp) :: local_point_dot(3), local_parameter_dot, local_normal_dot(3)
        integer :: local_triangle

        hit_point_dot = 0.0_dp
        connection_length_dot = 0.0_dp
        normal_dot = 0.0_dp
        status = 1
        if (.not. valid_polyline( &
            polyline, surface_vertices, surface_triangles)) return
        if (any(shape(polyline_dot) /= shape(polyline)) .or. &
            any(shape(surface_vertices_dot) /= shape(surface_vertices))) return
        if (any(.not. ieee_is_finite(polyline_dot)) .or. &
            any(.not. ieee_is_finite(surface_vertices_dot))) return
        call find_fci_first_hit_polyline_3d( &
            polyline, surface_vertices, surface_triangles, hit_point, &
            connection_length, hit_segment, hit_triangle, normal, primal_status)
        if (primal_status /= 0) return

        if (hit_segment == 0) then
            hit_point_dot = polyline_dot(:, size(polyline, 2))
            do segment = 1, size(polyline, 2) - 1
                segment_vector = polyline(:, segment + 1) - polyline(:, segment)
                segment_vector_dot = polyline_dot(:, segment + 1) - &
                    polyline_dot(:, segment)
                segment_length = vector_norm(segment_vector)
                segment_length_dot = dot_product(segment_vector, segment_vector_dot)/ &
                    segment_length
                connection_length_dot = connection_length_dot + segment_length_dot
            end do
            status = 0
            return
        end if

        cumulative_length_dot = 0.0_dp
        do segment = 1, hit_segment
            segment_vector = polyline(:, segment + 1) - polyline(:, segment)
            segment_vector_dot = polyline_dot(:, segment + 1) - &
                polyline_dot(:, segment)
            segment_length = vector_norm(segment_vector)
            segment_length_dot = dot_product(segment_vector, segment_vector_dot)/ &
                segment_length
            if (segment < hit_segment) then
                cumulative_length_dot = cumulative_length_dot + segment_length_dot
            else
                call find_fci_first_hit_triangle_3d( &
                    polyline(:, segment), polyline(:, segment + 1), surface_vertices, &
                    surface_triangles, local_point, local_parameter, local_triangle, &
                    local_normal, segment_status)
                if (segment_status /= 0 .or. local_triangle /= hit_triangle) return
                call find_fci_first_hit_triangle_3d_jvp( &
                    polyline(:, segment), polyline(:, segment + 1), surface_vertices, &
                    surface_triangles, polyline_dot(:, segment), &
                    polyline_dot(:, segment + 1), surface_vertices_dot, &
                    local_point_dot, &
                    local_parameter_dot, local_normal_dot, segment_status)
                if (segment_status /= 0) return
                connection_length_dot = cumulative_length_dot + &
                    local_parameter_dot*segment_length + &
                    local_parameter*segment_length_dot
                hit_point_dot = local_point_dot
                normal_dot = local_normal_dot
            end if
        end do
        status = 0
    end subroutine find_fci_first_hit_polyline_3d_jvp

    pure logical function valid_polyline( &
            polyline, surface_vertices, surface_triangles) result(valid)
        real(dp), intent(in) :: polyline(:, :), surface_vertices(:, :)
        integer, intent(in) :: surface_triangles(:, :)

        integer :: segment
        real(dp) :: scale, segment_vector(3)

        valid = .false.
        if (size(polyline, 1) /= 3 .or. size(polyline, 2) < 2 .or. &
            size(surface_vertices, 1) /= 3 .or. size(surface_vertices, 2) < 3 .or. &
            size(surface_triangles, 1) /= 3 .or. size(surface_triangles, 2) < 1) return
        if (any(.not. ieee_is_finite(polyline)) .or. &
            any(.not. ieee_is_finite(surface_vertices))) return
        if (any(surface_triangles < 1) .or. &
            any(surface_triangles > size(surface_vertices, 2))) return
        scale = max(1.0_dp, maxval(abs(polyline)), maxval(abs(surface_vertices)))
        do segment = 1, size(polyline, 2) - 1
            segment_vector = polyline(:, segment + 1) - polyline(:, segment)
            if (vector_norm(segment_vector) <= topology_tolerance*scale) return
        end do
        valid = .true.
    end function valid_polyline

    pure function vector_norm(vector) result(norm)
        real(dp), intent(in) :: vector(3)
        real(dp) :: norm

        norm = sqrt(dot_product(vector, vector))
    end function vector_norm

end module fortfem_fci_terminal_polyline_3d
