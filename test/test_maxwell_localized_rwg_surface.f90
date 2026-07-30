program test_maxwell_localized_rwg_surface
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_maxwell_localized_rwg_basis
    use fortfem_kinds, only: dp
    implicit none

    real(dp) :: divergence, edge_midpoint(3), value(3), vertices(3, 3)
    integer :: local_edge, status, triangles(3, 1)
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 2) = [2.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 3) = [0.0_dp, 3.0_dp, 0.0_dp]
    triangles(:, 1) = [1, 2, 3]

    do local_edge = 1, 3
        edge_midpoint = local_edge_midpoint(vertices, triangles(:, 1), local_edge)
        call evaluate_maxwell_localized_rwg_basis( &
            vertices, triangles, 1, local_edge, edge_midpoint, value, &
            divergence, status)
        call record_condition(status == 0, &
            "localized RWG evaluates on its owning panel")
        call record_condition(abs( &
            dot_product(value, outward_conormal( &
            vertices, triangles(:, 1), local_edge)) - 1.0_dp) < 2.0e-14_dp, &
            "localized RWG has unit conormal trace on its reference edge")
        call record_condition(abs(divergence - &
            edge_length(vertices, triangles(:, 1), local_edge)/3.0_dp) < &
            2.0e-14_dp, &
            "localized RWG divergence follows the affine Piola scaling")
    end do

    call evaluate_maxwell_localized_rwg_basis( &
        vertices, triangles, 1, 1, [4.0_dp, 0.0_dp, 0.0_dp], value, &
        divergence, status)
    call record_condition(status /= 0, &
        "localized RWG rejects points outside its owning panel")

    call check_summary("Localized Maxwell RWG surface space")
    if (.not. all_passed) error stop 1

contains

    pure function local_edge_vertices(local_edge) result(local_vertices)
        integer, intent(in) :: local_edge
        integer :: local_vertices(2)

        select case (local_edge)
        case (1)
            local_vertices = [1, 2]
        case (2)
            local_vertices = [3, 1]
        case default
            local_vertices = [2, 3]
        end select
    end function local_edge_vertices

    pure function local_edge_midpoint(points, element, local_edge) result(point)
        real(dp), intent(in) :: points(:, :)
        integer, intent(in) :: element(3), local_edge
        real(dp) :: point(3)
        integer :: local_vertices(2)

        local_vertices = local_edge_vertices(local_edge)
        point = 0.5_dp*( &
            points(:, element(local_vertices(1))) + &
            points(:, element(local_vertices(2))))
    end function local_edge_midpoint

    pure function edge_length(points, element, local_edge) result(length)
        real(dp), intent(in) :: points(:, :)
        integer, intent(in) :: element(3), local_edge
        real(dp) :: length
        integer :: local_vertices(2)

        local_vertices = local_edge_vertices(local_edge)
        length = norm2( &
            points(:, element(local_vertices(2))) - &
            points(:, element(local_vertices(1))))
    end function edge_length

    pure function outward_conormal(points, element, local_edge) result(conormal)
        real(dp), intent(in) :: points(:, :)
        integer, intent(in) :: element(3), local_edge
        real(dp) :: conormal(3), normal(3), tangent(3)
        integer :: local_vertices(2)

        local_vertices = local_edge_vertices(local_edge)
        tangent = points(:, element(local_vertices(2))) - &
            points(:, element(local_vertices(1)))
        normal = cross_product( &
            points(:, element(2)) - points(:, element(1)), &
            points(:, element(3)) - points(:, element(1)))
        normal = normal/norm2(normal)
        conormal = cross_product(tangent, normal)
        conormal = conormal/norm2(conormal)
    end function outward_conormal

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_localized_rwg_surface
