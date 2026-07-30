module fortfem_maxwell_localized_rwg_surface
    !! Element-local RWG functions used by barycentric Maxwell trace spaces.
    !!
    !! The local edges follow the Buffa--Christiansen/Bempp convention
    !! (1,2), (3,1), (2,3). Unlike a conforming RWG space, every panel owns
    !! three independent functions.
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: evaluate_maxwell_localized_rwg_basis

contains

    subroutine evaluate_maxwell_localized_rwg_basis( &
            vertices, triangles, triangle, local_edge, point, value, &
            surface_divergence, status)
        real(dp), intent(in) :: vertices(:, :), point(3)
        integer, intent(in) :: triangles(:, :), triangle, local_edge
        real(dp), intent(out) :: value(3), surface_divergence
        integer, intent(out) :: status

        real(dp) :: area, edge_length, panel(3, 3)
        integer :: edge_vertices(2), opposite

        value = 0.0_dp
        surface_divergence = 0.0_dp
        status = 1
        if (size(vertices, 1) /= 3 .or. size(triangles, 1) /= 3) return
        if (triangle < 1 .or. triangle > size(triangles, 2)) return
        if (local_edge < 1 .or. local_edge > 3) return
        if (any(triangles(:, triangle) < 1)) return
        if (any(triangles(:, triangle) > size(vertices, 2))) return
        panel = vertices(:, triangles(:, triangle))
        if (.not. point_in_triangle(panel, point)) return
        call local_edge_data(local_edge, edge_vertices, opposite)
        area = triangle_area(panel)
        if (area <= tiny(1.0_dp)) return
        edge_length = norm2( &
            panel(:, edge_vertices(2)) - panel(:, edge_vertices(1)))
        value = edge_length/(2.0_dp*area)*(point - panel(:, opposite))
        surface_divergence = edge_length/area
        status = 0
    end subroutine evaluate_maxwell_localized_rwg_basis

    pure subroutine local_edge_data(local_edge, edge_vertices, opposite)
        integer, intent(in) :: local_edge
        integer, intent(out) :: edge_vertices(2), opposite

        select case (local_edge)
        case (1)
            edge_vertices = [1, 2]
            opposite = 3
        case (2)
            edge_vertices = [3, 1]
            opposite = 2
        case default
            edge_vertices = [2, 3]
            opposite = 1
        end select
    end subroutine local_edge_data

    pure logical function point_in_triangle(vertices, point) result(inside)
        real(dp), intent(in) :: vertices(3, 3), point(3)

        real(dp) :: displacement(3), dot00, dot01, dot02, dot11, dot12
        real(dp) :: denominator, normal(3), u, v

        normal = cross_product( &
            vertices(:, 2) - vertices(:, 1), &
            vertices(:, 3) - vertices(:, 1))
        displacement = point - vertices(:, 1)
        if (abs(dot_product(displacement, normal)) > &
            128.0_dp*epsilon(1.0_dp)*max(1.0_dp, norm2(normal))) then
            inside = .false.
            return
        end if
        dot00 = dot_product( &
            vertices(:, 2) - vertices(:, 1), &
            vertices(:, 2) - vertices(:, 1))
        dot01 = dot_product( &
            vertices(:, 2) - vertices(:, 1), &
            vertices(:, 3) - vertices(:, 1))
        dot02 = dot_product(vertices(:, 2) - vertices(:, 1), displacement)
        dot11 = dot_product( &
            vertices(:, 3) - vertices(:, 1), &
            vertices(:, 3) - vertices(:, 1))
        dot12 = dot_product(vertices(:, 3) - vertices(:, 1), displacement)
        denominator = dot00*dot11 - dot01*dot01
        if (denominator <= tiny(1.0_dp)) then
            inside = .false.
            return
        end if
        u = (dot11*dot02 - dot01*dot12)/denominator
        v = (dot00*dot12 - dot01*dot02)/denominator
        inside = u >= -128.0_dp*epsilon(1.0_dp) .and. &
            v >= -128.0_dp*epsilon(1.0_dp) .and. &
            u + v <= 1.0_dp + 128.0_dp*epsilon(1.0_dp)
    end function point_in_triangle

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

end module fortfem_maxwell_localized_rwg_surface
