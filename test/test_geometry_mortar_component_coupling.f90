program test_geometry_mortar_component_coupling
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        assemble_geometry_mortar_component_coupling, &
        assemble_geometry_mortar_component_coupling_jvp, &
        assemble_geometry_mortar_component_coupling_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: tolerance = 3.0e-11_dp
    real(dp), parameter :: finite_difference_tolerance = 8.0e-10_dp
    real(dp), parameter :: step = 1.0e-6_dp
    real(dp) :: test_trace(2, 3, 2), trial_trace(2, 3, 2)
    real(dp) :: component_metric(2, 2, 3)
    real(dp) :: reference_weights(3), surface_jacobian(3)
    real(dp) :: matrix(2, 2), physical_weights(3), expected_matrix(2, 2)
    real(dp) :: test_trace_dot(2, 3, 2), trial_trace_dot(2, 3, 2)
    real(dp) :: component_metric_dot(2, 2, 3)
    real(dp) :: reference_weights_dot(3), surface_jacobian_dot(3)
    real(dp) :: matrix_dot(2, 2), physical_weights_dot(3)
    real(dp) :: matrix_plus(2, 2), matrix_minus(2, 2)
    real(dp) :: weights_plus(3), weights_minus(3)
    real(dp) :: matrix_bar(2, 2), physical_weights_bar(3)
    real(dp) :: test_trace_bar(2, 3, 2), trial_trace_bar(2, 3, 2)
    real(dp) :: component_metric_bar(2, 2, 3)
    real(dp) :: reference_weights_bar(3), surface_jacobian_bar(3)
    real(dp) :: tensor_test_trace(4, 2, 1), tensor_trial_trace(4, 2, 1)
    real(dp) :: tensor_metric(4, 4, 2), tensor_matrix(1, 1), tensor_weights(2)
    real(dp) :: tensor_reference_weights(2), tensor_surface_jacobian(2)
    real(dp) :: tensor_expected_matrix(1, 1)
    real(dp) :: forward_pairing, reverse_pairing
    real(dp) :: bad_surface_jacobian(3)
    logical :: all_passed
    type(fortsparse_status_t) :: status

    all_passed = .true.
    test_trace = reshape([ &
        0.8_dp, -0.2_dp, 0.3_dp, 0.1_dp, 0.4_dp, -0.7_dp, &
        -0.5_dp, 0.6_dp, 0.2_dp, 0.9_dp, -0.1_dp, 0.5_dp], [2, 3, 2])
    trial_trace = reshape([ &
        0.5_dp, 0.4_dp, -0.2_dp, 0.6_dp, -0.3_dp, 0.1_dp, &
        0.7_dp, -0.8_dp, 0.2_dp, 0.3_dp, 0.9_dp, -0.4_dp], [2, 3, 2])
    component_metric = reshape([ &
        1.0_dp, 0.2_dp, -0.1_dp, 1.3_dp, &
        0.7_dp, -0.4_dp, 0.5_dp, 1.1_dp, &
        1.2_dp, 0.3_dp, -0.6_dp, 0.8_dp], [2, 2, 3])
    reference_weights = [2.0_dp, 1.0_dp, 3.0_dp]
    surface_jacobian = [1.5_dp, 2.0_dp, 0.5_dp]
    expected_matrix = 0.0_dp
    call independent_cross_mass( &
        test_trace, trial_trace, component_metric, reference_weights, &
        surface_jacobian, expected_matrix)

    call assemble_geometry_mortar_component_coupling( &
        test_trace, trial_trace, reference_weights, surface_jacobian, &
        component_metric, matrix, physical_weights, status)
    call record_condition(status%code == 0, &
        "component mortar accepts vector-valued physical traces")
    call record_condition(maxval(abs(matrix - expected_matrix)) < tolerance, &
        "component mortar matches an independent tensor cross-mass oracle")
    call record_condition(maxval(abs(physical_weights - &
        reference_weights*surface_jacobian)) < tolerance, &
        "component mortar returns the physical quadrature measure")

    tensor_test_trace = reshape([ &
        0.2_dp, -0.3_dp, 0.4_dp, 0.5_dp, 0.6_dp, 0.1_dp, &
        -0.7_dp, 0.8_dp], [4, 2, 1])
    tensor_trial_trace = reshape([ &
        -0.4_dp, 0.9_dp, 0.3_dp, -0.2_dp, 0.5_dp, -0.6_dp, &
        0.7_dp, 0.1_dp], [4, 2, 1])
    tensor_metric = 0.0_dp
    tensor_metric(1, 1, :) = [1.0_dp, 0.8_dp]
    tensor_metric(2, 2, :) = [0.9_dp, 1.1_dp]
    tensor_metric(3, 3, :) = [1.2_dp, 0.7_dp]
    tensor_metric(4, 4, :) = [0.6_dp, 1.3_dp]
    tensor_metric(1, 2, :) = [0.1_dp, -0.2_dp]
    tensor_metric(2, 1, :) = [-0.3_dp, 0.4_dp]
    tensor_reference_weights = [1.0_dp, 0.8_dp]
    tensor_surface_jacobian = [1.2_dp, 0.9_dp]
    call independent_cross_mass( &
        tensor_test_trace, tensor_trial_trace, tensor_metric, &
        tensor_reference_weights, tensor_surface_jacobian, &
        tensor_expected_matrix)
    call assemble_geometry_mortar_component_coupling( &
        tensor_test_trace, tensor_trial_trace, tensor_reference_weights, &
        tensor_surface_jacobian, tensor_metric, tensor_matrix, tensor_weights, &
        status)
    call record_condition(status%code == 0 .and. &
        abs(tensor_matrix(1, 1) - tensor_expected_matrix(1, 1)) < tolerance, &
        "component mortar also contracts flattened tensor-valued traces")

    test_trace_dot = reshape([ &
        0.03_dp, -0.02_dp, 0.05_dp, -0.04_dp, 0.01_dp, 0.02_dp, &
        -0.01_dp, 0.06_dp, 0.02_dp, 0.04_dp, -0.03_dp, 0.05_dp], [2, 3, 2])
    trial_trace_dot = reshape([ &
        -0.01_dp, 0.06_dp, 0.02_dp, 0.04_dp, -0.03_dp, 0.05_dp, &
        0.02_dp, -0.04_dp, 0.03_dp, 0.01_dp, 0.05_dp, -0.02_dp], [2, 3, 2])
    component_metric_dot = reshape([ &
        0.02_dp, -0.03_dp, 0.04_dp, 0.01_dp, &
        -0.05_dp, 0.06_dp, 0.01_dp, -0.02_dp, &
        0.03_dp, 0.05_dp, -0.04_dp, 0.02_dp], [2, 2, 3])
    reference_weights_dot = [0.2_dp, -0.1_dp, 0.3_dp]
    surface_jacobian_dot = [0.04_dp, 0.05_dp, -0.02_dp]
    call assemble_geometry_mortar_component_coupling_jvp( &
        test_trace, trial_trace, reference_weights, surface_jacobian, &
        component_metric, test_trace_dot, trial_trace_dot, &
        reference_weights_dot, surface_jacobian_dot, component_metric_dot, &
        matrix_dot, physical_weights_dot, status)
    call record_condition(status%code == 0, &
        "component mortar fixed-topology JVP succeeds")
    call assemble_geometry_mortar_component_coupling( &
        test_trace + step*test_trace_dot, trial_trace + step*trial_trace_dot, &
        reference_weights + step*reference_weights_dot, &
        surface_jacobian + step*surface_jacobian_dot, &
        component_metric + step*component_metric_dot, matrix_plus, weights_plus, &
        status)
    call assemble_geometry_mortar_component_coupling( &
        test_trace - step*test_trace_dot, trial_trace - step*trial_trace_dot, &
        reference_weights - step*reference_weights_dot, &
        surface_jacobian - step*surface_jacobian_dot, &
        component_metric - step*component_metric_dot, matrix_minus, weights_minus, &
        status)
    call record_condition(maxval(abs(matrix_dot - &
        (matrix_plus - matrix_minus)/(2.0_dp*step))) < &
        finite_difference_tolerance, &
        "component mortar matrix JVP matches central differences")
    call record_condition(maxval(abs(physical_weights_dot - &
        (weights_plus - weights_minus)/(2.0_dp*step))) < &
        finite_difference_tolerance, &
        "component mortar metric JVP matches central differences")

    matrix_bar = reshape([0.4_dp, -0.3_dp, 0.2_dp, 0.5_dp], [2, 2])
    physical_weights_bar = [0.7_dp, -0.2_dp, 0.6_dp]
    call assemble_geometry_mortar_component_coupling_vjp( &
        test_trace, trial_trace, reference_weights, surface_jacobian, &
        component_metric, matrix_bar, physical_weights_bar, test_trace_bar, &
        trial_trace_bar, reference_weights_bar, surface_jacobian_bar, &
        component_metric_bar, status)
    call record_condition(status%code == 0, &
        "component mortar fixed-topology VJP succeeds")
    forward_pairing = sum(matrix_bar*matrix_dot) + &
        sum(physical_weights_bar*physical_weights_dot)
    reverse_pairing = sum(test_trace_bar*test_trace_dot) + &
        sum(trial_trace_bar*trial_trace_dot) + &
        sum(reference_weights_bar*reference_weights_dot) + &
        sum(surface_jacobian_bar*surface_jacobian_dot) + &
        sum(component_metric_bar*component_metric_dot)
    call record_condition(abs(forward_pairing - reverse_pairing) < tolerance, &
        "component mortar VJP satisfies the real dot-product identity")

    bad_surface_jacobian = [1.0_dp, 0.0_dp, 1.0_dp]
    call assemble_geometry_mortar_component_coupling( &
        test_trace, trial_trace, reference_weights, bad_surface_jacobian, &
        component_metric, matrix, physical_weights, status)
    call record_condition(status%code /= 0, &
        "component mortar rejects nonpositive surface metrics")

    call check_summary("geometry mortar component coupling")
    if (.not. all_passed) error stop 1

contains

    subroutine independent_cross_mass( &
            test_values, trial_values, metric, weights, jacobian, result)
        real(dp), intent(in) :: test_values(:, :, :), trial_values(:, :, :)
        real(dp), intent(in) :: metric(:, :, :), weights(:), jacobian(:)
        real(dp), intent(out) :: result(:, :)
        integer :: q, i, j, first_component, second_component

        result = 0.0_dp
        do q = 1, size(weights)
            do i = 1, size(test_values, 3)
                do j = 1, size(trial_values, 3)
                    do first_component = 1, size(test_values, 1)
                        do second_component = 1, size(trial_values, 1)
                            result(i, j) = result(i, j) + &
                                weights(q)*jacobian(q)* &
                                test_values(first_component, q, i)* &
                                metric(first_component, second_component, q)* &
                                trial_values(second_component, q, j)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine independent_cross_mass

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_geometry_mortar_component_coupling
