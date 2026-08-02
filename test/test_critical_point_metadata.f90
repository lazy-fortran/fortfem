program test_critical_point_metadata
    use check, only: check_condition, check_summary
    use fortfem_critical_point_metadata, only: &
        CRITICAL_CONTOUR_LIMITER, CRITICAL_CONTOUR_NONE, &
        CRITICAL_CONTOUR_SEPARATRIX, CRITICAL_POINT_DEGENERATE, &
        CRITICAL_POINT_O_POINT, CRITICAL_POINT_REGULAR, CRITICAL_POINT_X_POINT, &
        critical_point_metadata_t, evaluate_critical_point_metadata, &
        evaluate_critical_point_metadata_jvp, evaluate_critical_point_metadata_vjp, &
        validate_critical_point_metadata
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: candidate_count = 3
    real(dp), parameter :: gradient_tolerance = 1.0e-10_dp
    real(dp), parameter :: hessian_tolerance = 1.0e-10_dp
    real(dp), parameter :: step = 1.0e-6_dp
    real(dp) :: points(2, candidate_count), scalar_values(candidate_count)
    real(dp) :: gradients(2, candidate_count), hessians(2, 2, candidate_count)
    real(dp) :: gradients_dot(2, candidate_count), hessians_dot(2, 2, candidate_count)
    real(dp) :: event_margin_dot(candidate_count), event_margin_plus(candidate_count)
    real(dp) :: event_margin_minus(candidate_count), event_margin_bar(candidate_count)
    real(dp) :: points_bar(2, candidate_count), scalar_values_bar(candidate_count)
    real(dp) :: gradients_bar(2, candidate_count), hessians_bar(2, 2, candidate_count)
    integer :: contour_kind(candidate_count), bad_contour(1)
    type(critical_point_metadata_t) :: metadata
    type(fortsparse_status_t) :: status
    real(dp) :: lhs, rhs
    logical :: all_passed
    integer :: i

    all_passed = .true.
points = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 0.0_dp], [2, candidate_count])
    scalar_values = [0.2_dp, -0.3_dp, 0.7_dp]
    gradients = reshape([0.3_dp, 0.4_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
                        [2, candidate_count])
    hessians(:, :, 1) = reshape([1.0_dp, 0.1_dp, 0.2_dp, 1.5_dp], [2, 2])
    hessians(:, :, 2) = reshape([2.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
    hessians(:, :, 3) = reshape([1.0_dp, 0.0_dp, 0.0_dp, -2.0_dp], [2, 2])
    contour_kind = [CRITICAL_CONTOUR_NONE, CRITICAL_CONTOUR_LIMITER, &
                    CRITICAL_CONTOUR_SEPARATRIX]
    call evaluate_critical_point_metadata( &
        points, scalar_values, gradients, hessians, gradient_tolerance, &
        hessian_tolerance, contour_kind, metadata, status)
    call record(status%code == 0, "critical-point metadata evaluates")
    call record(all(metadata%classification == [CRITICAL_POINT_REGULAR, &
                                     CRITICAL_POINT_O_POINT, CRITICAL_POINT_X_POINT]), &
                "gradient/Hessian oracle classifies regular, O, and X candidates")
    call record(maxval(abs(metadata%event_margin - [ &
                       sqrt(0.25_dp) - gradient_tolerance, 2.0_dp - hessian_tolerance, &
                           2.0_dp - hessian_tolerance])) < 2.0e-12_dp, &
                "event margins match independent scalar oracle")
    call validate_critical_point_metadata(metadata, status)
    call record(status%code == 0, "metadata validates after construction")

    gradients_dot = reshape([0.02_dp, -0.01_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
                            [2, candidate_count])
    hessians_dot = reshape([(0.01_dp*real(i, dp), i=1, 12)], [2, 2, candidate_count])
    call evaluate_critical_point_metadata_jvp( &
        points, scalar_values, gradients, hessians, gradient_tolerance, &
        hessian_tolerance, contour_kind, gradients_dot, hessians_dot, &
        event_margin_dot, status)
    call record(status%code == 0, "critical-point margin JVP evaluates")
    call evaluate_critical_point_metadata( &
        points, scalar_values, gradients + step*gradients_dot, &
        hessians + step*hessians_dot, gradient_tolerance, hessian_tolerance, &
        contour_kind, metadata, status)
    event_margin_plus = metadata%event_margin
    call evaluate_critical_point_metadata( &
        points, scalar_values, gradients - step*gradients_dot, &
        hessians - step*hessians_dot, gradient_tolerance, hessian_tolerance, &
        contour_kind, metadata, status)
    event_margin_minus = metadata%event_margin
   call record(maxval(abs(event_margin_dot - (event_margin_plus - event_margin_minus)/ &
                  (2.0_dp*step))) < 2.0e-8_dp, "margin JVP matches central differences")

    event_margin_bar = [0.7_dp, -0.4_dp, 0.2_dp]
    call evaluate_critical_point_metadata_vjp( &
        points, scalar_values, gradients, hessians, gradient_tolerance, &
        hessian_tolerance, contour_kind, event_margin_bar, points_bar, &
        scalar_values_bar, gradients_bar, hessians_bar, status)
    lhs = sum(event_margin_bar*event_margin_dot)
    rhs = sum(gradients_bar*gradients_dot) + sum(hessians_bar*hessians_dot)
    call record(status%code == 0 .and. abs(lhs - rhs) < 2.0e-10_dp, &
                "margin VJP satisfies real dot-product oracle")

    bad_contour = [99]
    call evaluate_critical_point_metadata( &
        points(:, :1), scalar_values(:1), gradients(:, :1), hessians(:, :, :1), &
        gradient_tolerance, hessian_tolerance, bad_contour, metadata, status)
    call record(status%code /= 0, "invalid contour labels are rejected")
    call check_summary("critical point metadata")
    if (.not. all_passed) error stop 1

contains

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description
        call check_condition(condition, description)
        all_passed = all_passed .and. condition
    end subroutine record

end program test_critical_point_metadata
