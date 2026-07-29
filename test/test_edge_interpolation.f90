program test_edge_interpolation
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_edge_interpolation_2d, only: &
        interpolate_nedelec_edge_dofs, interpolate_rt_edge_dofs
    use check, only: check_condition, check_summary
    implicit none

    type(mesh_2d_t) :: mesh
    real(dp), allocatable :: nedelec_dofs(:), rt_dofs(:)
    real(dp) :: expected, point_a(2), point_b(2), edge_vector(2)
    integer :: edge_id, dof_id
    logical :: all_passed

    all_passed = .true.
    mesh%n_vertices = 4
    mesh%n_triangles = 2
    mesh%has_triangles = .true.
    allocate(mesh%vertices(2, 4), mesh%triangles(3, 2))
    mesh%vertices(:, 1) = [0.0_dp, 0.0_dp]
    mesh%vertices(:, 2) = [1.0_dp, 0.0_dp]
    mesh%vertices(:, 3) = [1.0_dp, 1.0_dp]
    mesh%vertices(:, 4) = [0.0_dp, 1.0_dp]
    mesh%triangles(:, 1) = [1, 2, 3]
    mesh%triangles(:, 2) = [1, 3, 4]
    call mesh%build_edge_connectivity()

    allocate(nedelec_dofs(mesh%n_edges), rt_dofs(mesh%n_edges))
    call interpolate_nedelec_edge_dofs(mesh, gradient_field, 2, nedelec_dofs)
    call interpolate_rt_edge_dofs(mesh, flux_field, 2, rt_dofs)

    do edge_id = 1, mesh%n_edges
        point_a = mesh%vertices(:, mesh%edges(1, edge_id))
        point_b = mesh%vertices(:, mesh%edges(2, edge_id))
        edge_vector = point_b - point_a
        dof_id = mesh%edge_to_dof(edge_id) + 1

        expected = potential(point_b) - potential(point_a)
        call record_condition( &
            abs(nedelec_dofs(dof_id) - expected) < 1.0e-13_dp, &
            "Nedelec interpolation satisfies the gradient line theorem")

        expected = edge_vector(2) * quadratic_average( &
            point_a(1), point_b(1)) - edge_vector(1) * quadratic_average( &
            point_a(2), point_b(2))
        call record_condition(abs(rt_dofs(dof_id) - expected) < 1.0e-13_dp, &
            "RT interpolation matches the exact polynomial normal flux")
    end do

    call check_summary("Oriented edge interpolation")
    if (.not. all_passed) error stop 1

contains

    pure subroutine gradient_field(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        value = [x**2 + y, x - y**2]
    end subroutine gradient_field

    pure subroutine flux_field(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        value = [x**2, y**2]
    end subroutine flux_field

    pure real(dp) function potential(point) result(value)
        real(dp), intent(in) :: point(2)

        value = point(1)**3 / 3.0_dp + point(1) * point(2) - &
            point(2)**3 / 3.0_dp
    end function potential

    pure real(dp) function quadratic_average(a, b) result(value)
        real(dp), intent(in) :: a, b

        value = (a**2 + a * b + b**2) / 3.0_dp
    end function quadratic_average

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_edge_interpolation
