program test_sphere_surface_mesh
    use check, only: check_condition, check_summary
    use fortfem_core, only: generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: vertices(:, :)
    real(dp) :: normal(3)
    integer :: triangle
    logical :: all_passed

    all_passed = .true.
    call generate_sphere_surface_mesh(2.5_dp, 2, vertices, triangles)
    call record_condition(size(vertices, 2) == 66 .and. &
        size(triangles, 2) == 128, &
        "Refined octahedral sphere has the analytical mesh counts")
    call record_condition(maxval(abs( &
        sqrt(sum(vertices**2, dim=1)) - 2.5_dp)) < 2.0e-14_dp, &
        "Every generated sphere vertex lies on the requested radius")
    do triangle = 1, size(triangles, 2)
        normal = cross_product( &
            vertices(:, triangles(2, triangle)) - &
            vertices(:, triangles(1, triangle)), &
            vertices(:, triangles(3, triangle)) - &
            vertices(:, triangles(1, triangle)))
        call record_condition(dot_product(normal, sum( &
            vertices(:, triangles(:, triangle)), dim=2)) > 0.0_dp, &
            "Sphere surface panels have outward orientation")
    end do
    call check_summary("Sphere surface mesh")
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

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_sphere_surface_mesh
