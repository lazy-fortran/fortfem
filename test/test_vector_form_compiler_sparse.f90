program test_vector_form_compiler_sparse
    use check, only: check_condition, check_summary
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    use fortfem_api, only: compile_vector_form_csc, curl, div, dx, &
        form_expr_t, init_measures, inner, interpolate_nedelec_edge_dofs, &
        interpolate_rt_edge_dofs, operator(*), operator(+), &
        vector_test_function_t, vector_trial_function_t
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    implicit none

    type(form_expr_t) :: form
    type(csc_t) :: matrix
    type(fortsparse_status_t) :: sparse_status
    type(mesh_2d_t) :: mesh
    type(vector_test_function_t) :: test_field
    type(vector_trial_function_t) :: trial_field
    real(dp), allocatable :: dofs(:), matrix_times_dofs(:)
    real(dp) :: energy
    logical :: all_passed

    all_passed = .true.
    call init_measures()
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

    form = (2.0_dp * inner(curl(trial_field), curl(test_field)) + &
        3.0_dp * inner(trial_field, test_field)) * dx
    call compile_vector_form_csc( &
        form, mesh, "Nedelec", 1, 4, matrix, sparse_status)
    allocate(dofs(mesh%n_edges), matrix_times_dofs(mesh%n_edges))
    call interpolate_nedelec_edge_dofs(mesh, constant_field, 2, dofs)
    matrix_times_dofs = csc_matvec(matrix, dofs)
    energy = dot_product(dofs, matrix_times_dofs)
    call record_condition(sparse_status%code == 0 .and. &
        matrix%nrow == mesh%n_edges .and. abs(energy - 3.0_dp) < 2.0e-12_dp, &
        "Sparse Nedelec compiler preserves exact shared-edge mass energy")

    form = (2.0_dp * inner(div(trial_field), div(test_field)) + &
        3.0_dp * inner(trial_field, test_field)) * dx
    call compile_vector_form_csc( &
        form, mesh, "RT", 0, 4, matrix, sparse_status)
    call interpolate_rt_edge_dofs(mesh, constant_field, 2, dofs)
    matrix_times_dofs = csc_matvec(matrix, dofs)
    energy = dot_product(dofs, matrix_times_dofs)
    call record_condition(sparse_status%code == 0 .and. &
        matrix%nrow == mesh%n_edges .and. abs(energy - 3.0_dp) < 2.0e-12_dp, &
        "Sparse RT compiler preserves exact shared-edge mass energy")

    call compile_vector_form_csc( &
        form, mesh, "unknown", 0, 4, matrix, sparse_status)
    call record_condition(sparse_status%code /= 0, &
        "Sparse vector compiler rejects an unknown family")

    call check_summary("Sparse vector form compiler")
    if (.not. all_passed) error stop 1

contains

    pure subroutine constant_field(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        value = [1.0_dp, 0.0_dp]
        associate (unused_coordinates => [x, y])
            if (size(unused_coordinates) /= 2) error stop
        end associate
    end subroutine constant_field

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_vector_form_compiler_sparse
