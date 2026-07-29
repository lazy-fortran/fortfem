program test_nedelec_reconstruction
    use check, only: check_condition, check_summary
    use fortfem_edge_interpolation_2d, only: interpolate_nedelec_edge_dofs
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_nedelec_field_2d, only: evaluate_nedelec_field_2d
    implicit none

    type(mesh_2d_t) :: mesh
    real(dp), allocatable :: real_dofs(:)
    complex(dp), allocatable :: dofs(:)
    complex(dp) :: curl, value(2)
    complex(dp), parameter :: scale = (1.0_dp, 2.0_dp)
    real(dp) :: point(2)
    integer :: triangle
    logical :: all_passed

    all_passed = .true.
    call create_unit_square(mesh)
    allocate(real_dofs(mesh%n_edges), dofs(mesh%n_edges))
    call interpolate_nedelec_edge_dofs( &
        mesh, nedelec_field, 2, real_dofs)
    dofs = scale * real_dofs

    do triangle = 1, mesh%n_triangles
        if (triangle == 1) then
            point = [0.75_dp, 0.25_dp]
        else
            point = [0.25_dp, 0.75_dp]
        end if
        call evaluate_nedelec_field_2d( &
            mesh, triangle, point(1), point(2), dofs, value, curl)
        call record_condition(maxval(abs(value - scale * [ &
            1.0_dp + 2.0_dp * point(2), &
            3.0_dp - 2.0_dp * point(1)])) < 1.0e-13_dp, &
            "Nedelec reconstruction reproduces an affine edge field")
        call record_condition(abs(curl + 4.0_dp * scale) < 1.0e-13_dp, &
            "Nedelec reconstruction returns the exact scalar curl")
    end do

    call check_summary("Nedelec reconstruction")
    if (.not. all_passed) error stop 1

contains

    subroutine create_unit_square(unit_mesh)
        type(mesh_2d_t), intent(out) :: unit_mesh

        unit_mesh%n_vertices = 4
        unit_mesh%n_triangles = 2
        unit_mesh%has_triangles = .true.
        allocate(unit_mesh%vertices(2, 4), unit_mesh%triangles(3, 2))
        unit_mesh%vertices(:, 1) = [0.0_dp, 0.0_dp]
        unit_mesh%vertices(:, 2) = [1.0_dp, 0.0_dp]
        unit_mesh%vertices(:, 3) = [1.0_dp, 1.0_dp]
        unit_mesh%vertices(:, 4) = [0.0_dp, 1.0_dp]
        unit_mesh%triangles(:, 1) = [1, 2, 3]
        unit_mesh%triangles(:, 2) = [1, 3, 4]
        call unit_mesh%build_edge_connectivity()
    end subroutine create_unit_square

    pure subroutine nedelec_field(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        value = [1.0_dp + 2.0_dp * y, 3.0_dp - 2.0_dp * x]
    end subroutine nedelec_field

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_nedelec_reconstruction
