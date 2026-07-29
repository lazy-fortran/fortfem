program test_laplace_bem_adjoint_double_layer
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_laplace_adjoint_double_layer_constant
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: matrix(4, 4), panel_end(2, 4), panel_start(2, 4)
    real(dp) :: expected
    integer :: panel, status
    logical :: all_passed

    all_passed = .true.

    panel_start(:, 1) = [0.0_dp, 0.0_dp]
    panel_end(:, 1) = [1.0_dp, 0.0_dp]
    panel_start(:, 2) = [0.0_dp, 0.0_dp]
    panel_end(:, 2) = [0.0_dp, 1.0_dp]
    call assemble_laplace_adjoint_double_layer_constant( &
        panel_start(:, 1:2), panel_end(:, 1:2), 24, matrix(1:2, 1:2), status)
    expected = -(0.25_dp * pi + 0.5_dp * log(2.0_dp)) / (2.0_dp * pi)
    call record_condition(status == 0 .and. &
        abs(matrix(1, 2) - expected) < 2.0e-14_dp, &
        "Adjoint double layer matches an exact target-normal integral")

    panel_start(:, 1) = [0.0_dp, 0.0_dp]
    panel_end(:, 1) = [1.0_dp, 0.0_dp]
    panel_start(:, 2) = [1.0_dp, 0.0_dp]
    panel_end(:, 2) = [1.0_dp, 1.0_dp]
    panel_start(:, 3) = [1.0_dp, 1.0_dp]
    panel_end(:, 3) = [0.0_dp, 1.0_dp]
    panel_start(:, 4) = [0.0_dp, 1.0_dp]
    panel_end(:, 4) = [0.0_dp, 0.0_dp]
    call assemble_laplace_adjoint_double_layer_constant( &
        panel_start, panel_end, 24, matrix, status)
    do panel = 1, 4
        call record_condition( &
            abs(sum(matrix(:, panel)) + 0.5_dp) < 3.0e-14_dp, &
            "Adjoint double layer satisfies the dual closed-boundary identity")
    end do

    call check_summary("Laplace BEM adjoint double-layer")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_laplace_bem_adjoint_double_layer
