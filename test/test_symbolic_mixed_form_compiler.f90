program test_symbolic_mixed_form_compiler
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        compile_mixed_form_csc, curl, div, dx, init_measures, inner, &
        interpolate_nedelec_edge_dofs, interpolate_rt_edge_dofs, operator(*), &
        test_function_t, trial_function_t, vector_test_function_t, &
        vector_trial_function_t
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    type(csc_t) :: matrix
    type(fortsparse_status_t) :: status
    type(mesh_2d_t) :: mesh
    type(test_function_t) :: scalar_test
    type(trial_function_t) :: scalar_trial
    type(vector_test_function_t) :: vector_test
    type(vector_trial_function_t) :: vector_trial
    real(dp), allocatable :: scalar_dofs(:), vector_dofs(:)
    real(dp) :: pairing
    logical :: all_passed

    all_passed = .true.
    call init_measures()
    call mesh%create_rectangular(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)

    call compile_mixed_form_csc( &
        2.5_dp*inner(div(vector_trial), scalar_test)*dx, mesh, "RT", 0, 4, &
        matrix, status)
    allocate(scalar_dofs(matrix%nrow), vector_dofs(matrix%ncol))
    scalar_dofs = 1.0_dp
    call interpolate_rt_edge_dofs(mesh, linear_rt_field, 4, vector_dofs)
    pairing = dot_product(scalar_dofs, csc_matvec(matrix, vector_dofs))
    call record_condition( &
        status%code == 0 .and. abs(pairing - 2.5_dp) < 2.0e-12_dp, &
        "symbolic RT-DG divergence block reproduces its continuum pairing")
    deallocate(scalar_dofs, vector_dofs)

    call compile_mixed_form_csc( &
        (-3.0_dp)*inner(scalar_trial, div(vector_test))*dx, mesh, "RT", 0, 4, &
        matrix, status)
    allocate(scalar_dofs(matrix%ncol), vector_dofs(matrix%nrow))
    scalar_dofs = 1.0_dp
    call interpolate_rt_edge_dofs(mesh, linear_rt_field, 4, vector_dofs)
    pairing = dot_product(vector_dofs, csc_matvec(matrix, scalar_dofs))
    call record_condition( &
        status%code == 0 .and. abs(pairing + 3.0_dp) < 2.0e-12_dp, &
        "symbolic scalar-RT block preserves transpose direction and sign")
    deallocate(scalar_dofs, vector_dofs)

    call compile_mixed_form_csc( &
        1.75_dp*inner(curl(vector_trial), scalar_test)*dx, mesh, &
        "Nedelec", 1, 4, matrix, status)
    allocate(scalar_dofs(matrix%nrow), vector_dofs(matrix%ncol))
    scalar_dofs = 1.0_dp
    call interpolate_nedelec_edge_dofs( &
        mesh, unit_curl_field, 4, vector_dofs)
    pairing = dot_product(scalar_dofs, csc_matvec(matrix, vector_dofs))
    call record_condition( &
        status%code == 0 .and. abs(pairing - 1.75_dp) < 2.0e-12_dp, &
        "symbolic Nedelec-DG curl block reproduces its continuum pairing")
    deallocate(scalar_dofs, vector_dofs)

    call check_summary("Symbolic rectangular mixed-form compiler")
    if (.not. all_passed) error stop 1

contains

    pure subroutine linear_rt_field(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        value = [x, 0.0_dp]
        associate(unused => y)
        end associate
    end subroutine linear_rt_field

    pure subroutine unit_curl_field(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        value = [-0.5_dp*y, 0.5_dp*x]
    end subroutine unit_curl_field

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_symbolic_mixed_form_compiler
