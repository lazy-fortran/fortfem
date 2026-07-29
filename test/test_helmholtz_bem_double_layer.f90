program test_helmholtz_bem_double_layer
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_helmholtz_double_layer_constant
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), parameter :: adjacent_reference = &
        (0.19708074187392954_dp, 0.11070749377973130_dp)
    complex(dp), parameter :: parallel_reference = &
        (0.01546902324733324_dp, 0.26101062907468830_dp)
    complex(dp) :: matrix(2, 2), reversed_matrix(2, 2)
    real(dp) :: panel_end(2, 2), panel_start(2, 2)
    integer :: status
    logical :: all_passed

    all_passed = .true.

    panel_start(:, 1) = [0.0_dp, 0.0_dp]
    panel_end(:, 1) = [1.0_dp, 0.0_dp]
    panel_start(:, 2) = [0.0_dp, 0.0_dp]
    panel_end(:, 2) = [0.0_dp, 1.0_dp]
    call assemble_helmholtz_double_layer_constant( &
        panel_start, panel_end, 1.5_dp, 48, matrix, status)
    call record_condition(status == 0 .and. &
        abs(matrix(1, 2) - adjacent_reference) < 2.0e-9_dp, &
        "Endpoint-singular Helmholtz double layer matches SciPy 1.18")
    call record_condition(abs(matrix(1, 1)) < 1.0e-15_dp .and. &
        abs(matrix(2, 2)) < 1.0e-15_dp, &
        "Straight-panel Helmholtz double-layer self terms vanish")

    panel_start(:, 2) = [0.0_dp, 1.0_dp]
    panel_end(:, 2) = [0.0_dp, 0.0_dp]
    call assemble_helmholtz_double_layer_constant( &
        panel_start, panel_end, 1.5_dp, 48, reversed_matrix, status)
    call record_condition(abs(reversed_matrix(1, 2) + matrix(1, 2)) < &
        2.0e-9_dp, "Reversing a source panel reverses the Helmholtz double layer")

    panel_start(:, 1) = [0.0_dp, 0.0_dp]
    panel_end(:, 1) = [1.0_dp, 0.0_dp]
    panel_start(:, 2) = [0.0_dp, 1.0_dp]
    panel_end(:, 2) = [1.0_dp, 1.0_dp]
    call assemble_helmholtz_double_layer_constant( &
        panel_start, panel_end, 2.0_dp, 48, matrix, status)
    call record_condition(status == 0 .and. &
        abs(matrix(1, 2) - parallel_reference) < 5.0e-11_dp, &
        "Regular Helmholtz double layer matches SciPy 1.18 quadrature")

    call check_summary("Helmholtz BEM double-layer")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_helmholtz_bem_double_layer
