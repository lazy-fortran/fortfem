program test_vector_plot_interpolation
    use fortfem_core, only: mesh_t, unit_square_mesh, vector_function_space_t
    use fortfem_feec, only: vector_function, vector_function_space, &
        vector_function_t
    use fortfem_edge_interpolation_2d, only: interpolate_rt_edge_dofs
    use fortfem_plot, only: compute_vector_plot_grid
    use fortfem_kinds, only: dp
    use check, only: check_condition, check_summary
    implicit none

    type(mesh_t) :: mesh
    type(vector_function_space_t) :: space
    type(vector_function_t) :: field
    real(dp) :: x_grid(5), y_grid(5), u_grid(5, 5), v_grid(5, 5)
    real(dp), allocatable :: rt_dofs(:)
    real(dp) :: edge_vector(2)
    integer :: edge, dof
    logical :: all_passed

    all_passed = .true.
    mesh = unit_square_mesh(4)
    allocate(rt_dofs(mesh%data%n_edges))
    space = vector_function_space(mesh, "Nedelec", 1)
    field = vector_function(space)

    ! The lowest-order Nedelec coefficient on an oriented edge is the
    ! line integral of the field.  These coefficients therefore represent
    ! the exact constant field (1, 0), independently of the mesh resolution.
    do edge = 1, mesh%data%n_edges
        dof = mesh%data%edge_to_dof(edge) + 1
        edge_vector = mesh%data%vertices(:, mesh%data%edges(2, edge)) - &
            mesh%data%vertices(:, mesh%data%edges(1, edge))
        field%values(dof, 1) = edge_vector(1)
    end do

    call compute_vector_plot_grid(field, 5, 5, x_grid, y_grid, u_grid, &
        v_grid)
    call record_condition(maxval(abs(u_grid - 1.0_dp)) < 1.0e-12_dp, &
        "Nedelec plot reconstruction preserves a constant x field")
    call record_condition(maxval(abs(v_grid)) < 1.0e-12_dp, &
        "Nedelec plot reconstruction preserves a zero y component")

    space = vector_function_space(mesh, "RT", 0)
    field = vector_function(space)
    call interpolate_rt_edge_dofs(mesh%data, constant_field, 2, rt_dofs)
    field%values(:, 1) = rt_dofs
    field%values(:, 2) = 0.0_dp
    call compute_vector_plot_grid(field, 5, 5, x_grid, y_grid, u_grid, &
        v_grid)
    call record_condition(maxval(abs(u_grid - 1.0_dp)) < 1.0e-12_dp, &
        "RT plot reconstruction preserves a constant x field")
    call record_condition(maxval(abs(v_grid)) < 1.0e-12_dp, &
        "RT plot reconstruction preserves a zero y component")
    call check_summary("Vector plot interpolation")

    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

    pure subroutine constant_field(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        associate(unused => x + y)
        end associate
        value = [1.0_dp, 0.0_dp]
    end subroutine constant_field

end program test_vector_plot_interpolation
