program test_geometry_mortar_trace_coupling
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_geometry_mortar_trace_coupling, &
        assemble_geometry_mortar_trace_coupling_jvp, &
        assemble_geometry_mortar_trace_coupling_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: tolerance = 2.0e-11_dp
    real(dp), parameter :: finite_difference_tolerance = 5.0e-10_dp
    real(dp), parameter :: step = 1.0e-6_dp
    real(dp) :: test_trace(3, 2), trial_trace(3, 2)
    real(dp) :: reference_weights(3), surface_jacobian(3)
    real(dp) :: matrix(2, 2), physical_weights(3)
    real(dp) :: expected_matrix(2, 2), expected_weights(3)
    real(dp) :: test_trace_dot(3, 2), trial_trace_dot(3, 2)
    real(dp) :: reference_weights_dot(3), surface_jacobian_dot(3)
    real(dp) :: matrix_dot(2, 2), physical_weights_dot(3)
    real(dp) :: matrix_plus(2, 2), matrix_minus(2, 2)
    real(dp) :: weights_plus(3), weights_minus(3)
    real(dp) :: matrix_bar(2, 2), physical_weights_bar(3)
    real(dp) :: test_trace_bar(3, 2), trial_trace_bar(3, 2)
    real(dp) :: reference_weights_bar(3), surface_jacobian_bar(3)
    real(dp) :: forward_pairing, reverse_pairing
    real(dp) :: bad_jacobian(3)
    real(dp) :: bad_test_trace_dot(2, 2)
    integer :: quadrature, test_dof, trial_dof
    logical :: all_passed
    type(fortsparse_status_t) :: status

    all_passed = .true.
    test_trace = reshape([ &
        0.8_dp, 0.2_dp, 0.1_dp, &
        0.2_dp, 0.7_dp, 0.3_dp], [3, 2])
    trial_trace = reshape([ &
        0.5_dp, 0.4_dp, 0.2_dp, &
        0.1_dp, 0.3_dp, 0.6_dp], [3, 2])
    reference_weights = [2.0_dp, 1.0_dp, 3.0_dp]
    surface_jacobian = [1.5_dp, 2.0_dp, 0.5_dp]
    expected_weights = reference_weights*surface_jacobian
    expected_matrix = 0.0_dp
    do quadrature = 1, 3
        do test_dof = 1, 2
            do trial_dof = 1, 2
                expected_matrix(test_dof, trial_dof) = &
                    expected_matrix(test_dof, trial_dof) + &
                    expected_weights(quadrature)* &
                    test_trace(quadrature, test_dof)* &
                    trial_trace(quadrature, trial_dof)
            end do
        end do
    end do

    call assemble_geometry_mortar_trace_coupling( &
        test_trace, trial_trace, reference_weights, surface_jacobian, matrix, &
        physical_weights, status)
    call record_condition(status%code == 0, &
        "geometry mortar accepts positive physical surface weights")
    call record_condition(maxval(abs(physical_weights - expected_weights)) < &
        tolerance, "geometry mortar returns reference-weight times metric")
    call record_condition(maxval(abs(matrix - expected_matrix)) < tolerance, &
        "geometry mortar matches an independent physical cross-mass oracle")

    test_trace_dot = reshape([ &
        0.03_dp, -0.02_dp, 0.05_dp, &
        -0.04_dp, 0.01_dp, 0.02_dp], [3, 2])
    trial_trace_dot = reshape([ &
        -0.01_dp, 0.06_dp, 0.02_dp, &
        0.04_dp, -0.03_dp, 0.05_dp], [3, 2])
    reference_weights_dot = [0.2_dp, -0.1_dp, 0.3_dp]
    surface_jacobian_dot = [0.04_dp, 0.05_dp, -0.02_dp]
    call assemble_geometry_mortar_trace_coupling_jvp( &
        test_trace, trial_trace, reference_weights, surface_jacobian, &
        test_trace_dot, trial_trace_dot, reference_weights_dot, &
        surface_jacobian_dot, matrix_dot, physical_weights_dot, status)
    call record_condition(status%code == 0, &
        "geometry mortar fixed-topology JVP succeeds")
    call assemble_geometry_mortar_trace_coupling( &
        test_trace + step*test_trace_dot, trial_trace + step*trial_trace_dot, &
        reference_weights + step*reference_weights_dot, &
        surface_jacobian + step*surface_jacobian_dot, matrix_plus, weights_plus, &
        status)
    call assemble_geometry_mortar_trace_coupling( &
        test_trace - step*test_trace_dot, trial_trace - step*trial_trace_dot, &
        reference_weights - step*reference_weights_dot, &
        surface_jacobian - step*surface_jacobian_dot, matrix_minus, weights_minus, &
        status)
    call record_condition(maxval(abs(matrix_dot - &
        (matrix_plus - matrix_minus)/(2.0_dp*step))) < &
        finite_difference_tolerance, &
        "geometry mortar matrix JVP matches central differences")
    call record_condition(maxval(abs(physical_weights_dot - &
        (weights_plus - weights_minus)/(2.0_dp*step))) < &
        finite_difference_tolerance, &
        "geometry mortar metric JVP matches central differences")

    bad_test_trace_dot = 0.0_dp
    call assemble_geometry_mortar_trace_coupling_jvp( &
        test_trace, trial_trace, reference_weights, surface_jacobian, &
        bad_test_trace_dot, trial_trace_dot, reference_weights_dot, &
        surface_jacobian_dot, matrix_plus, weights_plus, status)
    call record_condition(status%code /= 0, &
        "geometry mortar JVP rejects incompatible increment shapes")

    matrix_bar = reshape([0.4_dp, -0.3_dp, 0.2_dp, 0.5_dp], [2, 2])
    physical_weights_bar = [0.7_dp, -0.2_dp, 0.6_dp]
    call assemble_geometry_mortar_trace_coupling_vjp( &
        test_trace, trial_trace, reference_weights, surface_jacobian, &
        matrix_bar, physical_weights_bar, test_trace_bar, trial_trace_bar, &
        reference_weights_bar, surface_jacobian_bar, status)
    call record_condition(status%code == 0, &
        "geometry mortar fixed-topology VJP succeeds")
    forward_pairing = sum(matrix_bar*matrix_dot) + &
        sum(physical_weights_bar*physical_weights_dot)
    reverse_pairing = sum(test_trace_bar*test_trace_dot) + &
        sum(trial_trace_bar*trial_trace_dot) + &
        sum(reference_weights_bar*reference_weights_dot) + &
        sum(surface_jacobian_bar*surface_jacobian_dot)
    call record_condition(abs(forward_pairing - reverse_pairing) < tolerance, &
        "geometry mortar VJP satisfies the real dot-product identity")

    bad_jacobian = [1.0_dp, 0.0_dp, 1.0_dp]
    call assemble_geometry_mortar_trace_coupling( &
        test_trace, trial_trace, reference_weights, bad_jacobian, matrix, &
        physical_weights, status)
    call record_condition(status%code /= 0, &
        "geometry mortar rejects nonpositive surface metrics")
    call check_summary("geometry mortar trace coupling")

    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_geometry_mortar_trace_coupling
