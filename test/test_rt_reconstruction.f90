program test_rt_reconstruction
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_edge_interpolation_2d, only: interpolate_rt_edge_dofs
    use fortfem_rt_field_2d, only: evaluate_rt_field_2d, &
        reconstruct_axisymmetric_fourier_toroidal, rt_l2_norm
    use check, only: check_condition, check_summary
    implicit none

    type(mesh_2d_t) :: mesh
    real(dp), allocatable :: real_dofs(:)
    complex(dp), allocatable :: dofs(:), toroidal(:)
    complex(dp) :: value(2), divergence
    complex(dp), parameter :: scale = (1.0_dp, 2.0_dp)
    real(dp) :: point(2), expected_norm
    integer :: triangle_id
    logical :: all_passed

    all_passed = .true.
    call create_unit_square(mesh)
    allocate(real_dofs(mesh%n_edges), dofs(mesh%n_edges))
    allocate(toroidal(mesh%n_triangles))
    call interpolate_rt_edge_dofs(mesh, rt0_field, 1, real_dofs)
    dofs = scale * real_dofs

    do triangle_id = 1, mesh%n_triangles
        if (triangle_id == 1) then
            point = [0.75_dp, 0.25_dp]
        else
            point = [0.25_dp, 0.75_dp]
        end if
        call evaluate_rt_field_2d( &
            mesh, triangle_id, point(1), point(2), dofs, value, divergence)
        call record_condition(maxval(abs(value - scale * [ &
            1.0_dp + 2.0_dp * point(1), &
            3.0_dp + 2.0_dp * point(2)])) < 1.0e-13_dp, &
            "RT0 reconstruction reproduces an affine RT field")
        call record_condition(abs(divergence - 4.0_dp * scale) < 1.0e-13_dp, &
            "RT0 reconstruction has the exact element divergence")
    end do

    expected_norm = sqrt(5.0_dp * 62.0_dp / 3.0_dp)
    call record_condition(abs(rt_l2_norm(mesh, dofs) - expected_norm) < &
        1.0e-13_dp, "RT0 L2 norm matches an analytical polynomial integral")

    call reconstruct_axisymmetric_fourier_toroidal( &
        mesh, 2, dofs, toroidal)
    call record_condition(maxval(abs(toroidal - cmplx( &
        -4.0_dp, 2.0_dp, dp))) < 1.0e-13_dp, &
        "Fourier toroidal reconstruction enforces zero divergence")

    call check_summary("RT0 reconstruction")
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

    pure subroutine rt0_field(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        value = [1.0_dp + 2.0_dp * x, 3.0_dp + 2.0_dp * y]
    end subroutine rt0_field

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_rt_reconstruction
