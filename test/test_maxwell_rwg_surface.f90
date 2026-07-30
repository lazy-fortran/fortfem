program test_maxwell_rwg_surface
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_maxwell_rwg_surface_space, evaluate_maxwell_rwg_basis
    use fortfem_kinds, only: dp
    implicit none

    real(dp), allocatable :: divergence(:, :), values(:, :, :)
    real(dp) :: conormal(3), edge_midpoint(3), normal(3), vertices(3, 4)
    integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
    integer :: basis, status, triangle, triangles(3, 2)
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 3) = [1.0_dp, 1.0_dp, 1.0_dp]
    vertices(:, 4) = [0.0_dp, 1.0_dp, 1.0_dp]
    triangles(:, 1) = [1, 2, 3]
    triangles(:, 2) = [1, 3, 4]

    call build_maxwell_rwg_surface_space( &
        vertices, triangles, edge_vertices, edge_triangles, status)
    call record_condition(status == 0 .and. size(edge_vertices, 2) == 1, &
        "RWG topology creates one basis function on a shared surface edge")
    basis = 1
    edge_midpoint = 0.5_dp*( &
        vertices(:, edge_vertices(1, basis)) + &
        vertices(:, edge_vertices(2, basis)))
    allocate(values(3, 2, 1), divergence(2, 1))
    do triangle = 1, 2
        call evaluate_maxwell_rwg_basis( &
            vertices, triangles, edge_vertices, edge_triangles, basis, &
            triangle, edge_midpoint, values(:, triangle, basis), &
            divergence(triangle, basis), status)
        normal = triangle_normal(vertices(:, triangles(:, triangle)))
        call record_condition(status == 0 .and. &
            abs(dot_product(values(:, triangle, basis), normal)) < 2.0e-14_dp, &
            "RWG basis is tangential to each physical panel")
    end do

    conormal = outward_conormal( &
        vertices(:, triangles(:, 1)), triangles(:, 1), edge_vertices(:, basis))
    call record_condition(abs( &
        dot_product(values(:, 1, basis), conormal) + &
        dot_product(values(:, 2, basis), outward_conormal( &
        vertices(:, triangles(:, 2)), triangles(:, 2), &
        edge_vertices(:, basis)))) < &
        2.0e-14_dp, &
        "RWG conormal flux is continuous across the shared edge")
    call record_condition(abs(sum(divergence(:, basis)*[ &
        triangle_area(vertices(:, triangles(:, 1))), &
        triangle_area(vertices(:, triangles(:, 2))) ])) < 2.0e-14_dp, &
        "RWG surface divergence has zero integral over its support")

    call evaluate_maxwell_rwg_basis( &
        vertices, triangles, edge_vertices, edge_triangles, basis, 1, &
        edge_midpoint + [2.0_dp, 0.0_dp, 0.0_dp], values(:, 1, basis), &
        divergence(1, basis), status)
    call record_condition(status /= 0, &
        "RWG evaluation rejects a point outside the selected panel")
    call check_summary("Maxwell RWG surface trace space")
    if (.not. all_passed) error stop 1

contains

    pure function triangle_normal(points) result(normal)
        real(dp), intent(in) :: points(3, 3)
        real(dp) :: normal(3)

        normal = cross_product( &
            points(:, 2) - points(:, 1), points(:, 3) - points(:, 1))
        normal = normal/norm2(normal)
    end function triangle_normal

    pure function triangle_area(points) result(area)
        real(dp), intent(in) :: points(3, 3)
        real(dp) :: area

        area = 0.5_dp*norm2(cross_product( &
            points(:, 2) - points(:, 1), points(:, 3) - points(:, 1)))
    end function triangle_area

    pure function outward_conormal(points, vertex_ids, edge) result(conormal)
        real(dp), intent(in) :: points(3, 3)
        integer, intent(in) :: vertex_ids(3), edge(2)
        real(dp) :: conormal(3), tangent(3)
        integer :: local, next

        do local = 1, 3
            next = modulo(local, 3) + 1
            if ((vertex_ids(local) == edge(1) .and. &
                vertex_ids(next) == edge(2)) .or. &
                (vertex_ids(local) == edge(2) .and. &
                vertex_ids(next) == edge(1))) exit
        end do
        tangent = points(:, next) - points(:, local)
        conormal = cross_product(tangent, triangle_normal(points))
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

end program test_maxwell_rwg_surface
