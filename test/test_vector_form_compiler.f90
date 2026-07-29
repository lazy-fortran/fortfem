program test_vector_form_compiler
    use check, only: check_condition, check_summary
    use fortfem_api, only: compile_vector_form_element, curl, div, dx, &
        form_expr_t, init_measures, inner, &
        interpolate_triangle_bdm, &
        interpolate_triangle_nedelec_second_kind, operator(*), operator(+), &
        vector_test_function_t, vector_trial_function_t
    use fortfem_kinds, only: dp
    implicit none

    type(form_expr_t) :: form
    type(vector_test_function_t) :: test_field
    type(vector_trial_function_t) :: trial_field
    real(dp), allocatable :: dofs(:), matrix(:, :)
    real(dp) :: energy
    real(dp) :: vertices(2, 3)
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call init_measures()
    vertices(:, 1) = [0.0_dp, 0.0_dp]
    vertices(:, 2) = [1.0_dp, 0.0_dp]
    vertices(:, 3) = [0.0_dp, 1.0_dp]

    form = (2.0_dp * inner(curl(trial_field), curl(test_field)) + &
        3.0_dp * inner(trial_field, test_field)) * dx
    call compile_vector_form_element( &
        form, "Nedelec", 1, vertices, 4, matrix, status)
    call record_condition(status == 0 .and. size(matrix, 1) == 3 .and. &
        abs(matrix(2, 2) - 4.5_dp) < 2.0e-13_dp, &
        "Compiled Nedelec form reproduces exact curl and mass energy")

    form = (2.0_dp * inner(div(trial_field), div(test_field)) + &
        3.0_dp * inner(trial_field, test_field)) * dx
    call compile_vector_form_element( &
        form, "RT", 0, vertices, 4, matrix, status)
    call record_condition(status == 0 .and. size(matrix, 1) == 3 .and. &
        abs(matrix(2, 2) - 4.5_dp) < 2.0e-13_dp, &
        "Compiled RT form reproduces exact divergence and mass energy")

    form = (2.0_dp * inner(curl(trial_field), curl(test_field)) + &
        3.0_dp * inner(trial_field, test_field)) * dx
    call compile_vector_form_element( &
        form, "Nedelec2", 1, vertices, 4, matrix, status)
    energy = huge(1.0_dp)
    if (status == 0) then
        call interpolate_triangle_nedelec_second_kind( &
            vertices, 1, 4, linear_field, dofs, status)
        if (status == 0) energy = dot_product(dofs, matmul(matrix, dofs))
    end if
    call record_condition(status == 0 .and. size(matrix, 1) == 6 .and. &
        abs(energy - 0.5_dp) < 3.0e-12_dp, &
        "Compiled second-kind Nedelec form has exact polynomial energy")

    form = (2.0_dp * inner(div(trial_field), div(test_field)) + &
        3.0_dp * inner(trial_field, test_field)) * dx
    call compile_vector_form_element( &
        form, "BDM", 1, vertices, 4, matrix, status)
    energy = huge(1.0_dp)
    if (status == 0) then
        call interpolate_triangle_bdm( &
            vertices, 1, 4, linear_field, dofs, status)
        if (status == 0) energy = dot_product(dofs, matmul(matrix, dofs))
    end if
    call record_condition(status == 0 .and. size(matrix, 1) == 6 .and. &
        abs(energy - 4.5_dp) < 3.0e-12_dp, &
        "Compiled BDM form has exact polynomial energy")

    call compile_vector_form_element( &
        form, "unknown", 1, vertices, 4, matrix, status)
    call record_condition(status /= 0, &
        "Vector compiler rejects an unknown element family")

    call check_summary("Vector form compiler")
    if (.not. all_passed) error stop 1

contains

    pure subroutine linear_field(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        value = [x, y]
    end subroutine linear_field

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_vector_form_compiler
