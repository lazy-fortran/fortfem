program test_scalar_form_compiler
    use check, only: check_condition, check_summary
    use fortfem_api, only: compile_form_matrix, compile_form_vector, &
        constant, dx, form_expr_t, function_t, grad, init_measures, inner, &
        operator(*), operator(+), test_function_t, trial_function_t
    use fortfem_kinds, only: dp
    implicit none

    type(form_expr_t) :: form
    type(test_function_t) :: v
    type(trial_function_t) :: u
    type(function_t) :: source
    real(dp) :: load(3)
    real(dp) :: expected(3, 3), mass(3, 3), matrix(3, 3)
    real(dp) :: stiffness(3, 3), vertices(2, 3)
    integer :: status, triangles(3, 1)
    logical :: all_passed

    all_passed = .true.
    call init_measures()
    vertices(:, 1) = [0.0_dp, 0.0_dp]
    vertices(:, 2) = [2.0_dp, 0.0_dp]
    vertices(:, 3) = [0.0_dp, 1.0_dp]
    triangles(:, 1) = [1, 2, 3]
    stiffness = reshape([ &
        1.25_dp, -0.25_dp, -1.0_dp, &
        -0.25_dp, 0.25_dp, 0.0_dp, &
        -1.0_dp, 0.0_dp, 1.0_dp], [3, 3])
    mass = reshape([ &
        2.0_dp, 1.0_dp, 1.0_dp, &
        1.0_dp, 2.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 2.0_dp], [3, 3]) / 12.0_dp

    form = (2.0_dp * inner(grad(u), grad(v)) + &
        3.0_dp * inner(u, v)) * dx
    call compile_form_matrix(form, vertices, triangles, matrix, status)
    expected = 2.0_dp * stiffness + 3.0_dp * mass
    call record_condition(status == 0 .and. &
        maxval(abs(matrix - expected)) < 2.0e-14_dp, &
        "Compiled scalar form matches exact stiffness and mass matrices")

    form = dx * (inner(u, v) + inner(grad(u), grad(v)))
    call compile_form_matrix(form, vertices, triangles, matrix, status)
    call record_condition(status == 0 .and. &
        maxval(abs(matrix - stiffness - mass)) < 2.0e-14_dp, &
        "Measure placement does not change compiled bilinear form")

    form = inner(grad(u), v) * dx
    call compile_form_matrix(form, vertices, triangles, matrix, status)
    call record_condition(status /= 0, &
        "Compiler rejects an unsupported mixed-rank integrand")

    source = constant(4.0_dp)
    form = source * v * dx
    call compile_form_vector(form, vertices, triangles, load, status)
    call record_condition(status == 0 .and. &
        maxval(abs(load - 4.0_dp / 3.0_dp)) < 2.0e-14_dp, &
        "Compiled constant load matches the exact P1 element integral")

    call check_summary("Scalar form compiler")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_scalar_form_compiler
