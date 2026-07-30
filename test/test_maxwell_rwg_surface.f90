program test_maxwell_rwg_surface
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_maxwell_rwg_surface_space, evaluate_maxwell_rwg_basis, &
        map_maxwell_rwg_to_tetra_nedelec_edges
    use fortfem_kinds, only: dp
    implicit none

    real(dp), allocatable :: divergence(:, :), values(:, :, :)
    real(dp) :: conormal(3), edge_midpoint(3), normal(3), vertices(3, 4)
    integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
    integer, allocatable :: trace_dofs(:)
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

    call test_tetra_nedelec_trace()
    call check_summary("Maxwell RWG surface trace space")
    if (.not. all_passed) error stop 1

contains

    subroutine test_tetra_nedelec_trace()
        real(dp), allocatable :: basis_value(:, :), basis_divergence(:)
        real(dp), allocatable :: trace_coefficients(:), trace_scales(:)
        real(dp) :: centroid(3), electric_field(3), expected(3), normal(3)
        real(dp) :: tetra_vertices(3, 4)
        integer, allocatable :: rwg_triangles(:, :), rwg_vertices(:, :)
        integer :: boundary_triangles(3, 4), edge, face, local_status
        integer :: tetrahedra(4, 1)

        tetra_vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
        tetra_vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
        tetra_vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
        tetra_vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
        tetrahedra(:, 1) = [1, 2, 3, 4]
        boundary_triangles(:, 1) = [1, 3, 2]
        boundary_triangles(:, 2) = [1, 2, 4]
        boundary_triangles(:, 3) = [1, 4, 3]
        boundary_triangles(:, 4) = [2, 3, 4]
        call build_maxwell_rwg_surface_space( &
            tetra_vertices, boundary_triangles, rwg_vertices, rwg_triangles, &
            local_status)
        call map_maxwell_rwg_to_tetra_nedelec_edges( &
            tetra_vertices, tetrahedra, rwg_vertices, trace_dofs, &
            trace_scales, local_status)
        call record_condition(local_status == 0 .and. &
            all(trace_dofs >= 1) .and. all(trace_dofs <= 6), &
            "RWG trace unknowns map to global tetrahedral Nedelec edges")

        electric_field = [0.7_dp, -0.4_dp, 0.2_dp]
        allocate( &
            trace_coefficients(6), basis_value(3, 6), basis_divergence(6))
        do edge = 1, 6
            trace_coefficients(edge) = trace_scales(edge)*dot_product( &
                electric_field, &
                tetra_vertices(:, rwg_vertices(2, edge)) - &
                tetra_vertices(:, rwg_vertices(1, edge)))
        end do
        do face = 1, 4
            centroid = sum( &
                tetra_vertices(:, boundary_triangles(:, face)), dim=2)/3.0_dp
            basis_value = 0.0_dp
            do edge = 1, 6
                if (.not. any(rwg_triangles(:, edge) == face)) cycle
                call evaluate_maxwell_rwg_basis( &
                    tetra_vertices, boundary_triangles, rwg_vertices, &
                    rwg_triangles, edge, face, centroid, &
                    basis_value(:, edge), basis_divergence(edge), local_status)
            end do
            normal = triangle_normal( &
                tetra_vertices(:, boundary_triangles(:, face)))
            expected = cross_product(electric_field, normal)
            call record_condition(norm2( &
                matmul(basis_value, trace_coefficients) - expected) < &
                2.0e-14_dp, &
                "RWG functions reproduce the rotated Nedelec tangential trace")
        end do
    end subroutine test_tetra_nedelec_trace

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
