program test_helmholtz_near_boundary_evaluation
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        evaluate_helmholtz_combined_potential_adaptive_constant
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), parameter :: reference = &
        (-0.18622628721349554_dp, -0.30653456272199375_dp)
    complex(dp) :: density(1), field(1)
    real(dp) :: error_estimate(1), panel_end(2, 1), panel_start(2, 1)
    real(dp) :: point(2, 1)
    integer :: status
    logical :: all_passed

    all_passed = .true.
    panel_start(:, 1) = [0.0_dp, 0.0_dp]
    panel_end(:, 1) = [1.0_dp, 0.0_dp]
    density(1) = (1.0_dp, 0.0_dp)
    point(:, 1) = [0.5_dp, 1.0e-4_dp]

    call evaluate_helmholtz_combined_potential_adaptive_constant( &
        panel_start, panel_end, 1.3_dp, 1.3_dp, density, point, 8, &
        1.0e-11_dp, 20, field, error_estimate, status)
    call record_condition(status == 0 .and. &
        abs(field(1) - reference) < 2.0e-10_dp, &
        "Adaptive combined potential matches a near-boundary SciPy integral")
    call record_condition(error_estimate(1) < 2.0e-10_dp, &
        "Adaptive combined potential reports a converged quadrature estimate")

    point(:, 1) = [0.5_dp, 0.0_dp]
    call evaluate_helmholtz_combined_potential_adaptive_constant( &
        panel_start, panel_end, 1.3_dp, 1.3_dp, density, point, 8, &
        1.0e-11_dp, 20, field, error_estimate, status)
    call record_condition(status /= 0, &
        "Off-surface evaluator rejects a target on the source panel")

    call check_summary("Helmholtz near-boundary field evaluation")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_helmholtz_near_boundary_evaluation
