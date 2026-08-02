program test_equation_objective_merit
    !! Independent weighted objective/constraint merit oracle.
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        OBJECTIVE_METADATA_KIND_OBJECTIVE, &
        equation_objective_metadata_t, &
        evaluate_equation_objective_merit, &
        evaluate_equation_objective_merit_jvp, &
        evaluate_equation_objective_merit_vjp, &
        initialize_equation_objective_metadata
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t
    implicit none

    integer, parameter :: dp = real64
    real(dp), parameter :: finite_difference_step = 1.0e-7_dp
    type(equation_objective_metadata_t) :: metadata, inactive_metadata
    type(fortsparse_status_t) :: status
    real(dp) :: values(2), values_dot(2), values_plus(2), values_minus(2)
    real(dp) :: values_bar(2), merit, merit_dot, merit_plus, merit_minus
    real(dp) :: expected_merit, expected_merit_dot, fd_dot, adjoint_left
    real(dp) :: adjoint_right

    call initialize_equation_objective_metadata( &
        metadata, "force-balance", OBJECTIVE_METADATA_KIND_OBJECTIVE, &
        target=[1.0_dp, -1.0_dp], lower_bound=[-2.0_dp, -2.0_dp], &
        upper_bound=[2.0_dp, 2.0_dp], weight=[2.0_dp, 0.5_dp], &
        scale=[2.0_dp, 0.5_dp], active=.true., fixed=.false., kkt_id=3, &
        nullspace_id=0, status=status, units="force", provenance="mms")
    call check_condition(status%code == FORTSPARSE_OK, &
        "objective metadata initializes for merit evaluation")

    values = [1.4_dp, -0.8_dp]
    values_dot = [0.3_dp, -0.2_dp]
    call evaluate_equation_objective_merit(metadata, values, merit, status)
    call check_condition(status%code == FORTSPARSE_OK, &
        "objective merit accepts finite values")
    expected_merit = 0.5_dp*(2.0_dp*(0.4_dp/2.0_dp)**2 + &
        0.5_dp*(0.2_dp/0.5_dp)**2)
    call check_condition(abs(merit - expected_merit) < 1.0e-14_dp, &
        "objective merit matches independent weighted residual oracle")

    call evaluate_equation_objective_merit_jvp( &
        metadata, values, values_dot, merit_dot, status)
    expected_merit_dot = 2.0_dp*(0.4_dp/2.0_dp)*(0.3_dp/2.0_dp) + &
        0.5_dp*(0.2_dp/0.5_dp)*(-0.2_dp/0.5_dp)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        abs(merit_dot - expected_merit_dot) < 1.0e-14_dp, &
        "objective merit JVP matches independent product rule")

    values_plus = values + finite_difference_step*values_dot
    values_minus = values - finite_difference_step*values_dot
    call evaluate_equation_objective_merit( &
        metadata, values_plus, merit_plus, status)
    call evaluate_equation_objective_merit( &
        metadata, values_minus, merit_minus, status)
    fd_dot = (merit_plus - merit_minus)/(2.0_dp*finite_difference_step)
    call check_condition(abs(merit_dot - fd_dot) < 2.0e-9_dp, &
        "objective merit JVP matches centered difference")

    call evaluate_equation_objective_merit_vjp( &
        metadata, values, 1.0_dp, values_bar, status)
    adjoint_left = dot_product(values_bar, values_dot)
    adjoint_right = merit_dot
    call check_condition(status%code == FORTSPARSE_OK .and. &
        abs(adjoint_left - adjoint_right) < 1.0e-14_dp, &
        "objective merit VJP satisfies real adjoint identity")

    metadata%active = .false.
    call evaluate_equation_objective_merit(metadata, values, merit, status)
    call check_condition(status%code == FORTSPARSE_OK .and. merit == 0.0_dp, &
        "inactive objective block contributes no merit")

    call initialize_equation_objective_metadata( &
        inactive_metadata, "bad", OBJECTIVE_METADATA_KIND_OBJECTIVE, &
        target=[0.0_dp], lower_bound=[-1.0_dp], upper_bound=[1.0_dp], &
        weight=[1.0_dp], scale=[1.0_dp], active=.true., fixed=.false., &
        kkt_id=0, nullspace_id=0, status=status)
    call evaluate_equation_objective_merit( &
        inactive_metadata, [0.0_dp, 1.0_dp], merit, status)
    call check_condition(status%code == FORTSPARSE_INVALID_MATRIX, &
        "objective merit rejects incompatible value shape")

    call check_summary("equation objective merit")
end program test_equation_objective_merit
