program test_bem_sphere_surface_solution
    !! Independent oracle for the exterior monopole used by the 3-D gallery.
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_laplace_representation_triangles_3d, &
        generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), allocatable :: vertices(:, :), dirichlet(:), neumann(:)
    integer, allocatable :: triangles(:, :)
    real(dp) :: targets(3, 4), value, radius
    real(dp) :: first_edge(3), second_edge(3), normal(3), centroid(3)
    integer :: element, point, status
    logical :: all_passed

    all_passed = .true.
    call generate_sphere_surface_mesh(1.0_dp, 2, vertices, triangles)
    allocate(dirichlet(size(vertices, 2)), neumann(size(triangles, 2)))
    dirichlet = 1.0_dp
    do element = 1, size(triangles, 2)
        first_edge = vertices(:, triangles(2, element)) - &
            vertices(:, triangles(1, element))
        second_edge = vertices(:, triangles(3, element)) - &
            vertices(:, triangles(1, element))
        normal = cross_product(first_edge, second_edge)
        normal = normal/norm2(normal)
        centroid = sum(vertices(:, triangles(:, element)), dim=2)/3.0_dp
        neumann(element) = -dot_product(normal, centroid)/norm2(centroid)**3
    end do

    targets = reshape([ &
        0.0_dp, 0.0_dp, 1.25_dp, &
        1.40_dp, 0.0_dp, 0.15_dp, &
        -0.8_dp, 0.7_dp, 1.0_dp, &
        0.0_dp, -1.6_dp, 0.0_dp], [3, 4])
    do point = 1, size(targets, 2)
        radius = norm2(targets(:, point))
        call evaluate_laplace_representation_triangles_3d( &
            vertices, triangles, dirichlet, neumann, targets(:, point), &
            8, value, status)
        call record_condition(status == 0 .and. &
            abs(value - 1.0_dp/radius) < 3.5e-2_dp, &
            "sphere exterior representation reproduces 1/r")
    end do

    call check_summary("3-D exterior BEM surface solution")
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

end program test_bem_sphere_surface_solution
