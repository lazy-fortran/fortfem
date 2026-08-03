program test_weighted_boundary_response_diagnostics
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: evaluate_weighted_boundary_response_diagnostics
    implicit none

    integer, parameter :: dp = real64, n = 3
    complex(dp) :: response(n, n), weighted_response(n, n), hermitian(n, n)
    real(dp) :: weights(n), expected_reciprocity, expected_passivity
    real(dp) :: reciprocity_error, passivity_lower_bound, scale, radius
    integer :: row, status

    response = reshape([ &
        cmplx(2.0_dp, 0.0_dp, dp), cmplx(0.2_dp, 0.3_dp, dp), &
        cmplx(-0.1_dp, 0.0_dp, dp), cmplx(0.2_dp, -0.3_dp, dp), &
        cmplx(1.5_dp, 0.0_dp, dp), cmplx(0.25_dp, 0.1_dp, dp), &
        cmplx(-0.1_dp, 0.0_dp, dp), cmplx(0.25_dp, -0.1_dp, dp), &
        cmplx(1.0_dp, 0.0_dp, dp)], shape(response))
    weights = [2.0_dp, 1.0_dp, 3.0_dp]
    weighted_response = spread(weights, 2, n)*response
    hermitian = 0.5_dp*(weighted_response + conjg(transpose(weighted_response)))
    scale = max(1.0_dp, maxval(abs(weighted_response)))
    expected_reciprocity = maxval(abs( &
        weighted_response - transpose(weighted_response)))/scale
    expected_passivity = huge(1.0_dp)
    do row = 1, n
        radius = sum(abs(hermitian(row, :))) - abs(hermitian(row, row))
        expected_passivity = min(expected_passivity, &
            real(hermitian(row, row), dp) - radius)
    end do

    call evaluate_weighted_boundary_response_diagnostics( &
        response, weights, reciprocity_error, passivity_lower_bound, status)
    call check_condition(status == 0 .and. &
        abs(reciprocity_error - expected_reciprocity) < 1.0e-14_dp, &
        "weighted reciprocity certificate matches direct work-pairing oracle")
    call check_condition(status == 0 .and. &
        abs(passivity_lower_bound - expected_passivity) < 1.0e-14_dp, &
        "weighted passivity certificate matches Hermitian Gershgorin oracle")

    call evaluate_weighted_boundary_response_diagnostics( &
        response, [2.0_dp, -1.0_dp, 3.0_dp], reciprocity_error, &
        passivity_lower_bound, status)
    call check_condition(status /= 0, &
        "weighted response diagnostics reject nonpositive work weights")

    call evaluate_weighted_boundary_response_diagnostics( &
        response(:, 1:2), weights, reciprocity_error, passivity_lower_bound, status)
    call check_condition(status /= 0, &
        "weighted response diagnostics reject nonsquare responses")

    call check_summary("weighted boundary response diagnostics")
end program test_weighted_boundary_response_diagnostics
