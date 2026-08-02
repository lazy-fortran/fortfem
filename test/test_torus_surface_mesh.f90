program test_torus_surface_mesh
    use check, only: check_condition, check_summary
    use fortfem_core, only: generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), allocatable :: vertices(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp) :: centroid(3), normal(3), radial(3), residual
    real(dp) :: first_edge(3), second_edge(3)
    integer :: element
    logical :: all_outward, all_passed

    all_passed = .true.
    call generate_torus_surface_mesh(2.0_dp, 0.6_dp, 12, 9, &
        vertices, triangles)
    call record_condition(size(vertices, 2) == 108, &
        "Torus mesh has one vertex per periodic parameter node")
    call record_condition(size(triangles, 2) == 216, &
        "Torus mesh has two triangles per parameter cell")

    residual = 0.0_dp
    do element = 1, size(vertices, 2)
        residual = max(residual, abs( &
            (sqrt(sum(vertices(1:2, element)**2)) - 2.0_dp)**2 + &
            vertices(3, element)**2 - 0.6_dp**2))
    end do
    call record_condition(residual < 2.0e-14_dp, &
        "Every torus vertex satisfies the implicit torus equation")

    all_outward = .true.
    do element = 1, size(triangles, 2)
        centroid = sum(vertices(:, triangles(:, element)), dim=2)/3.0_dp
        first_edge = vertices(:, triangles(2, element)) - &
            vertices(:, triangles(1, element))
        second_edge = vertices(:, triangles(3, element)) - &
            vertices(:, triangles(1, element))
        normal = cross_product(first_edge, second_edge)
        radial = [centroid(1), centroid(2), 0.0_dp]
        radial = radial/max(norm2(radial), tiny(1.0_dp))
        radial = centroid - 2.0_dp*radial
        all_outward = all_outward .and. dot_product(normal, radial) > 0.0_dp
    end do
    call record_condition(all_outward, &
        "Torus triangles have outward orientation")

    call check_summary("Toroidal surface mesh")
    if (.not. all_passed) error stop 1

contains

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

        call check_condition(condition, description)
        all_passed = all_passed .and. condition
    end subroutine record_condition
end program test_torus_surface_mesh
