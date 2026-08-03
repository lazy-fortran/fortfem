program test_laplace_bem_single_layer
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: assemble_laplace_single_layer_constant
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: panel_start(2, 3), panel_end(2, 3)
    real(dp) :: matrix(3, 3), expected
    integer :: status
    logical :: all_passed

    all_passed = .true.

    panel_start(:, 1) = [0.0_dp, 0.0_dp]
    panel_end(:, 1) = [2.0_dp, 0.0_dp]
    call assemble_laplace_single_layer_constant( &
        panel_start(:, 1:1), panel_end(:, 1:1), 24, matrix(1:1, 1:1), status)
    expected = 4.0_dp * (1.5_dp - log(2.0_dp)) / (2.0_dp * pi)
    call record_condition(status == 0 .and. &
        abs(matrix(1, 1) - expected) < 2.0e-14_dp, &
        "Single-layer self panel uses the exact logarithmic integral")

    panel_start(:, 1) = [1.0_dp, 0.0_dp]
    panel_end(:, 1) = [0.0_dp, 0.0_dp]
    panel_start(:, 2) = [0.0_dp, 0.0_dp]
    panel_end(:, 2) = [0.0_dp, 1.0_dp]
    call assemble_laplace_single_layer_constant( &
        panel_start(:, 1:2), panel_end(:, 1:2), 24, matrix(1:2, 1:2), status)
    expected = -(0.5_dp * log(2.0_dp) - 1.5_dp + 0.25_dp * pi) / &
        (2.0_dp * pi)
    call record_condition(status == 0 .and. &
        abs(matrix(1, 2) - expected) < 2.0e-14_dp, &
        "Duffy quadrature resolves a right-angle endpoint singularity")
    call record_condition(abs(matrix(2, 1) - matrix(1, 2)) < 1.0e-15_dp, &
        "Single-layer Galerkin matrix is symmetric")

    panel_start(:, 1) = [0.0_dp, 0.0_dp]
    panel_end(:, 1) = [1.0_dp, 0.0_dp]
    panel_start(:, 2) = [2.0_dp, 0.0_dp]
    panel_end(:, 2) = [3.0_dp, 0.0_dp]
    call assemble_laplace_single_layer_constant( &
        panel_start(:, 1:2), panel_end(:, 1:2), 24, matrix(1:2, 1:2), status)
    expected = -(4.5_dp * log(3.0_dp) - 4.0_dp * log(2.0_dp) - 1.5_dp) / &
        (2.0_dp * pi)
    call record_condition(status == 0 .and. &
        abs(matrix(1, 2) - expected) < 2.0e-14_dp, &
        "Regular panel interaction matches an exact logarithmic integral")

    panel_start(:, 1) = [0.0_dp, 0.0_dp]
    panel_end(:, 1) = [0.0_dp, 0.0_dp]
    call assemble_laplace_single_layer_constant( &
        panel_start(:, 1:1), panel_end(:, 1:1), 24, matrix(1:1, 1:1), status)
    call record_condition(status /= 0, &
        "Single-layer assembly rejects a zero-length panel")

    call check_summary("Laplace BEM single-layer")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_laplace_bem_single_layer
