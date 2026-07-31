module fortfem_fci_terminal_segment_2d
    !! Exact first-hit search for a straight FCI trace in a 2D facet set.
    !!
    !! The trace is parameterized as `x(t) = start + t*(finish-start)` with
    !! `0 < t <= 1`.  Facets are oriented segments; their returned right-hand
    !! normal is `[tangent(2), -tangent(1)]/norm(tangent)`.  A valid search with
    !! no intersection returns status zero and hit_segment zero.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none

    private

    real(dp), parameter :: topology_tolerance = 64.0_dp*epsilon(1.0_dp)

    public :: find_fci_first_hit_segment_2d
    public :: find_fci_first_hit_segment_2d_jvp

contains

    subroutine find_fci_first_hit_segment_2d( &
            start_point, finish_point, wall_vertices, wall_segments, hit_point, &
            hit_parameter, hit_segment, normal, status)
        !! Find the nearest transverse wall-facet intersection on a trace.
        real(dp), intent(in) :: start_point(2), finish_point(2)
        real(dp), intent(in) :: wall_vertices(:, :)
        integer, intent(in) :: wall_segments(:, :)
        real(dp), intent(out) :: hit_point(2), hit_parameter, normal(2)
        integer, intent(out) :: hit_segment, status

        integer :: facet, first, second, vertex_count, facet_count
        real(dp) :: trace(2), edge(2), offset(2), denominator
        real(dp) :: trace_norm, edge_norm, parameter, edge_parameter
        real(dp) :: best_parameter, scale

        hit_point = 0.0_dp
        hit_parameter = 0.0_dp
        hit_segment = 0
        normal = 0.0_dp
        status = 1
        vertex_count = size(wall_vertices, 2)
        facet_count = size(wall_segments, 2)
        if (vertex_count < 2 .or. facet_count < 1) return
        if (size(wall_vertices, 1) /= 2 .or. size(wall_segments, 1) /= 2) return
        if (any(.not. ieee_is_finite(start_point)) .or. &
            any(.not. ieee_is_finite(finish_point)) .or. &
            any(.not. ieee_is_finite(wall_vertices))) return
        if (any(wall_segments < 1) .or. any(wall_segments > vertex_count)) return
        trace = finish_point - start_point
        trace_norm = sqrt(dot_product(trace, trace))
        scale = max(1.0_dp, maxval(abs(wall_vertices)), maxval(abs(start_point)), &
            maxval(abs(finish_point)))
        if (trace_norm <= topology_tolerance*scale) return

        do facet = 1, facet_count
            first = wall_segments(1, facet)
            second = wall_segments(2, facet)
            edge = wall_vertices(:, second) - wall_vertices(:, first)
            edge_norm = sqrt(dot_product(edge, edge))
            if (edge_norm <= topology_tolerance*scale) return
        end do

        best_parameter = huge(1.0_dp)
        do facet = 1, facet_count
            first = wall_segments(1, facet)
            second = wall_segments(2, facet)
            edge = wall_vertices(:, second) - wall_vertices(:, first)
            offset = wall_vertices(:, first) - start_point
            denominator = cross_2d(trace, edge)
            if (abs(denominator) <= topology_tolerance*trace_norm* &
                sqrt(dot_product(edge, edge))) cycle
            parameter = cross_2d(offset, edge)/denominator
            edge_parameter = cross_2d(offset, trace)/denominator
            if (parameter <= topology_tolerance .or. &
                parameter > 1.0_dp + topology_tolerance .or. &
                edge_parameter < -topology_tolerance .or. &
                edge_parameter > 1.0_dp + topology_tolerance) cycle
            if (parameter < best_parameter - topology_tolerance) then
                best_parameter = parameter
                hit_segment = facet
                hit_parameter = min(1.0_dp, max(0.0_dp, parameter))
                hit_point = start_point + hit_parameter*trace
                normal = [edge(2), -edge(1)]/sqrt(dot_product(edge, edge))
            end if
        end do
        status = 0
    end subroutine find_fci_first_hit_segment_2d

    subroutine find_fci_first_hit_segment_2d_jvp( &
            start_point, finish_point, wall_vertices, wall_segments, &
            start_point_dot, finish_point_dot, wall_vertices_dot, hit_point_dot, &
            hit_parameter_dot, normal_dot, status)
        !! Apply the fixed-topology JVP of the first-hit geometry.
        real(dp), intent(in) :: start_point(2), finish_point(2)
        real(dp), intent(in) :: wall_vertices(:, :)
        integer, intent(in) :: wall_segments(:, :)
        real(dp), intent(in) :: start_point_dot(2), finish_point_dot(2)
        real(dp), intent(in) :: wall_vertices_dot(:, :)
        real(dp), intent(out) :: hit_point_dot(2), hit_parameter_dot, normal_dot(2)
        integer, intent(out) :: status

        real(dp) :: hit_point(2), hit_parameter, normal(2)
        real(dp) :: trace(2), trace_dot(2), edge(2), edge_dot(2)
        real(dp) :: offset(2), offset_dot(2), denominator, denominator_dot
        real(dp) :: numerator, numerator_dot, edge_norm, edge_norm_dot
        integer :: hit_segment, first, second, primal_status

        hit_point_dot = 0.0_dp
        hit_parameter_dot = 0.0_dp
        normal_dot = 0.0_dp
        status = 1
        if (size(wall_vertices, 1) /= 2 .or. &
            size(wall_vertices_dot, 1) /= 2 .or. &
            any(shape(wall_vertices_dot) /= shape(wall_vertices)) .or. &
            size(wall_segments, 1) /= 2) return
        if (any(.not. ieee_is_finite(start_point_dot)) .or. &
            any(.not. ieee_is_finite(finish_point_dot)) .or. &
            any(.not. ieee_is_finite(wall_vertices_dot))) return
        call find_fci_first_hit_segment_2d( &
            start_point, finish_point, wall_vertices, wall_segments, hit_point, &
            hit_parameter, hit_segment, normal, primal_status)
        if (primal_status /= 0) return
        if (hit_segment == 0) then
            status = 0
            return
        end if

        first = wall_segments(1, hit_segment)
        second = wall_segments(2, hit_segment)
        trace = finish_point - start_point
        trace_dot = finish_point_dot - start_point_dot
        edge = wall_vertices(:, second) - wall_vertices(:, first)
        edge_dot = wall_vertices_dot(:, second) - wall_vertices_dot(:, first)
        offset = wall_vertices(:, first) - start_point
        offset_dot = wall_vertices_dot(:, first) - start_point_dot
        denominator = cross_2d(trace, edge)
        denominator_dot = cross_2d(trace_dot, edge) + cross_2d(trace, edge_dot)
        if (abs(denominator) <= topology_tolerance* &
            sqrt(dot_product(trace, trace))*sqrt(dot_product(edge, edge))) return
        numerator = cross_2d(offset, edge)
        numerator_dot = cross_2d(offset_dot, edge) + cross_2d(offset, edge_dot)
        hit_parameter_dot = (numerator_dot*denominator - numerator*denominator_dot)/ &
            denominator**2
        hit_point_dot = start_point_dot + hit_parameter_dot*trace + &
            hit_parameter*trace_dot
        edge_norm = sqrt(dot_product(edge, edge))
        edge_norm_dot = dot_product(edge, edge_dot)/edge_norm
        normal_dot = [edge_dot(2), -edge_dot(1)]/edge_norm - &
            normal*edge_norm_dot/edge_norm
        status = 0
    end subroutine find_fci_first_hit_segment_2d_jvp

    pure real(dp) function cross_2d(first, second)
        real(dp), intent(in) :: first(2), second(2)

        cross_2d = first(1)*second(2) - first(2)*second(1)
    end function cross_2d

end module fortfem_fci_terminal_segment_2d
