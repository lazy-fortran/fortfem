program test_cell_vector_source
    use check, only: check_condition, check_summary
    use fortfem_core, only: mesh_t, rectangle_mesh, vector_function_space_t
    use fortfem_feec, only: cell_vector_source, cell_vector_source_t, &
        compile_vector_form_rhs, dx, form_expr_t, init_measures, inner, &
        operator(*), vector_function_space, vector_test_function, &
        vector_test_function_t
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    type(cell_vector_source_t) :: source
    type(form_expr_t) :: form
    type(fortsparse_status_t) :: status
    type(mesh_t) :: mesh
    type(vector_function_space_t) :: space
    type(vector_test_function_t) :: test_field
    real(dp), allocatable :: field_dofs_x(:), field_dofs_y(:), load(:)
    real(dp), allocatable :: source_values(:, :)
    real(dp) :: edge_vector(2), expected_x, expected_y
    integer :: degree_of_freedom, edge
    logical :: all_passed

    all_passed = .true.
    call init_measures()
    mesh = rectangle_mesh(2, 2, [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp])
    space = vector_function_space(mesh, "Nedelec", 1)
    test_field = vector_test_function(space)

    allocate(source_values(2, mesh%data%n_triangles))
    source_values(:, 1) = [2.0_dp, -1.0_dp]
    source_values(:, 2) = [-3.0_dp, 4.0_dp]
    source = cell_vector_source(source_values)
    source_values = -100.0_dp
    form = 2.5_dp * inner(source, test_field) * dx
    source%values = -200.0_dp

    allocate(load(space%ndof))
    call compile_vector_form_rhs( &
        form, mesh%data, "Nedelec", 1, 4, load, status)
    allocate(field_dofs_x(space%ndof), field_dofs_y(space%ndof))
    field_dofs_x = 0.0_dp
    field_dofs_y = 0.0_dp
    do edge = 1, mesh%data%n_edges
        degree_of_freedom = mesh%data%edge_to_dof(edge) + 1
        edge_vector = mesh%data%vertices(:, mesh%data%edges(2, edge)) - &
            mesh%data%vertices(:, mesh%data%edges(1, edge))
        field_dofs_x(degree_of_freedom) = edge_vector(1)
        field_dofs_y(degree_of_freedom) = edge_vector(2)
    end do
    expected_x = 2.5_dp * 0.5_dp * (2.0_dp - 3.0_dp)
    expected_y = 2.5_dp * 0.5_dp * (-1.0_dp + 4.0_dp)
    call record_condition(status%code == 0, &
        "Cell vector source compiles for an order-one Nedelec space")
    call record_condition(abs(dot_product(field_dofs_x, load) - expected_x) < &
        3.0e-14_dp, "Cell source has the exact constant-x test integral")
    call record_condition(abs(dot_product(field_dofs_y, load) - expected_y) < &
        3.0e-14_dp, "Cell source has the exact constant-y test integral")

    source = cell_vector_source(reshape([1.0_dp, 2.0_dp], [2, 1]))
    form = inner(source, test_field) * dx
    call compile_vector_form_rhs( &
        form, mesh%data, "Nedelec", 1, 4, load, status)
    call record_condition(status%code /= 0, &
        "Cell source rejects a value count that differs from the mesh")

    call check_summary("Cell vector source")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_cell_vector_source
