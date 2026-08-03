program test_enrichment_support_vector_diagnostics
    use check, only: check_condition, check_summary
    use fortfem_feec, only: evaluate_enrichment_support_vector_gram, &
        evaluate_enrichment_support_vector_gram_jvp, &
        evaluate_enrichment_support_vector_gram_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, FORTSPARSE_OK
    implicit none

    integer, parameter :: component_count = 2
    integer, parameter :: point_count = 3
    integer, parameter :: basis_count = 3
    real(dp), parameter :: epsilon_fd = 1.0e-7_dp
    real(dp), parameter :: sample_weights(point_count) = [ &
        0.5_dp, 1.25_dp, 0.75_dp]
    real(dp), parameter :: sample_weights_dot(point_count) = [ &
        0.04_dp, -0.03_dp, 0.02_dp]
    logical, parameter :: active_mask(basis_count) = [.true., .true., .false.]
    real(dp), parameter :: gram_bar(basis_count, basis_count) = reshape([ &
        0.2_dp, -0.1_dp, 0.3_dp, &
        -0.4_dp, 0.5_dp, 0.2_dp, &
        0.1_dp, -0.2_dp, 0.6_dp], [basis_count, basis_count])
    real(dp) :: reference_values(component_count, point_count, basis_count)
    real(dp) :: jacobians(component_count, component_count, point_count)
    real(dp) :: basis_values(component_count, point_count, basis_count)
    real(dp) :: basis_dot(component_count, point_count, basis_count)
    real(dp) :: metric(component_count, component_count, point_count)
    real(dp) :: metric_dot(component_count, component_count, point_count)
    real(dp) :: metric_plus(component_count, component_count, point_count)
    real(dp) :: metric_minus(component_count, component_count, point_count)
    real(dp) :: basis_plus(component_count, point_count, basis_count)
    real(dp) :: basis_minus(component_count, point_count, basis_count)
    real(dp) :: weights_plus(point_count), weights_minus(point_count)
    real(dp) :: gram(basis_count, basis_count), expected_gram(basis_count, basis_count)
    real(dp) :: gram_dot(basis_count, basis_count)
    real(dp) :: gram_plus(basis_count, basis_count)
    real(dp) :: gram_minus(basis_count, basis_count)
    real(dp) :: basis_bar(component_count, point_count, basis_count)
    real(dp) :: metric_bar(component_count, component_count, point_count)
    real(dp) :: weights_bar(point_count)
    real(dp) :: lhs, rhs
    real(dp) :: invalid_metric(component_count, component_count, point_count)
    integer :: point, first_basis, second_basis
    type(fortsparse_status_t) :: status
    logical :: all_passed

    all_passed = .true.
    call initialize_reference_data(reference_values, jacobians, metric, metric_dot)
    call map_covariant_values(reference_values, jacobians, basis_values)
    basis_dot = reshape([ &
        0.11_dp, -0.07_dp, 0.05_dp, 0.13_dp, -0.09_dp, 0.02_dp, &
        -0.03_dp, 0.08_dp, 0.12_dp, -0.04_dp, 0.06_dp, -0.1_dp, &
        0.02_dp, 0.03_dp, -0.05_dp, 0.04_dp, -0.08_dp, 0.01_dp], &
        shape(basis_values))

    expected_gram = 0.0_dp
    do first_basis = 1, basis_count
        if (.not. active_mask(first_basis)) cycle
        do second_basis = 1, basis_count
            if (.not. active_mask(second_basis)) cycle
            do point = 1, point_count
                expected_gram(first_basis, second_basis) = &
                    expected_gram(first_basis, second_basis) + &
                    sample_weights(point)*dot_product( &
                    basis_values(:, point, first_basis), &
                    matmul(metric(:, :, point), &
                    basis_values(:, point, second_basis)))
            end do
        end do
    end do

    call evaluate_enrichment_support_vector_gram( &
        basis_values, metric, sample_weights, active_mask, gram, status)
    call record(status%code == FORTSPARSE_OK .and. &
        maxval(abs(gram - expected_gram)) < 1.0e-14_dp .and. &
        maxval(abs(gram - transpose(gram))) < 1.0e-14_dp, &
        "vector support Gram matches the independent Piola/metric oracle")

    call evaluate_enrichment_support_vector_gram_jvp( &
        basis_values, metric, sample_weights, active_mask, basis_dot, &
        metric_dot, sample_weights_dot, gram_dot, status)
    basis_plus = basis_values + epsilon_fd*basis_dot
    basis_minus = basis_values - epsilon_fd*basis_dot
    metric_plus = metric + epsilon_fd*metric_dot
    metric_minus = metric - epsilon_fd*metric_dot
    weights_plus = sample_weights + epsilon_fd*sample_weights_dot
    weights_minus = sample_weights - epsilon_fd*sample_weights_dot
    call evaluate_enrichment_support_vector_gram( &
        basis_plus, metric_plus, weights_plus, active_mask, gram_plus, status)
    call evaluate_enrichment_support_vector_gram( &
        basis_minus, metric_minus, weights_minus, active_mask, gram_minus, status)
    call record(status%code == FORTSPARSE_OK .and. &
        maxval(abs(gram_dot - (gram_plus - gram_minus)/ &
        (2.0_dp*epsilon_fd))) < 2.0e-8_dp, &
        "vector support Gram JVP matches fixed-topology finite differences")

    call evaluate_enrichment_support_vector_gram_vjp( &
        basis_values, metric, sample_weights, active_mask, gram_bar, basis_bar, &
        metric_bar, weights_bar, status)
    lhs = sum(gram_bar*gram_dot)
    rhs = sum(basis_bar*basis_dot) + sum(metric_bar*metric_dot) + &
        sum(weights_bar*sample_weights_dot)
    call record(status%code == FORTSPARSE_OK .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "vector support Gram VJP satisfies the real dot-product identity")

    invalid_metric = metric
    invalid_metric(1, 2, 1) = invalid_metric(1, 2, 1) + 0.25_dp
    call evaluate_enrichment_support_vector_gram( &
        basis_values, invalid_metric, sample_weights, active_mask, gram, status)
    call record(status%code /= FORTSPARSE_OK, &
        "vector support Gram rejects a nonsymmetric metric")

    call check_summary("XFEM vector support diagnostics")
    if (.not. all_passed) error stop 1

contains

    subroutine initialize_reference_data( &
            values, jacobian, metric_value, metric_dot_value)
        real(dp), intent(out) :: values(:, :, :), jacobian(:, :, :)
        real(dp), intent(out) :: metric_value(:, :, :), metric_dot_value(:, :, :)

        values(:, :, 1) = reshape([ &
            1.0_dp, 2.0_dp, -0.5_dp, 1.5_dp, 2.0_dp, -1.0_dp], [2, 3])
        values(:, :, 2) = reshape([ &
            -1.0_dp, 0.5_dp, 0.75_dp, -2.0_dp, 1.25_dp, 0.25_dp], [2, 3])
        values(:, :, 3) = reshape([ &
            0.4_dp, -0.8_dp, 1.1_dp, 0.3_dp, -1.2_dp, 1.7_dp], [2, 3])
        jacobian(:, :, 1) = reshape([2.0_dp, 0.0_dp, 0.25_dp, 1.5_dp], [2, 2])
        jacobian(:, :, 2) = reshape([1.4_dp, 0.2_dp, -0.1_dp, 1.8_dp], [2, 2])
        jacobian(:, :, 3) = reshape([1.7_dp, -0.15_dp, 0.3_dp, 1.2_dp], [2, 2])
        metric_value(:, :, 1) = reshape([2.0_dp, 0.2_dp, 0.2_dp, 1.5_dp], [2, 2])
        metric_value(:, :, 2) = reshape([1.5_dp, -0.1_dp, -0.1_dp, 2.0_dp], [2, 2])
        metric_value(:, :, 3) = reshape([1.2_dp, 0.15_dp, 0.15_dp, 1.8_dp], [2, 2])
        metric_dot_value(:, :, 1) = reshape([ &
            0.03_dp, 0.01_dp, 0.01_dp, -0.02_dp], [2, 2])
        metric_dot_value(:, :, 2) = reshape([ &
            -0.02_dp, 0.02_dp, 0.02_dp, 0.04_dp], [2, 2])
        metric_dot_value(:, :, 3) = reshape([ &
            0.01_dp, -0.015_dp, -0.015_dp, 0.02_dp], [2, 2])
    end subroutine initialize_reference_data

    subroutine map_covariant_values(reference, jacobian, physical)
        real(dp), intent(in) :: reference(:, :, :), jacobian(:, :, :)
        real(dp), intent(out) :: physical(:, :, :)
        real(dp) :: determinant
        integer :: point, basis

        physical = 0.0_dp
        do point = 1, size(reference, 2)
            determinant = jacobian(1, 1, point)*jacobian(2, 2, point) - &
                jacobian(1, 2, point)*jacobian(2, 1, point)
            do basis = 1, size(reference, 3)
                physical(1, point, basis) = ( &
                    jacobian(2, 2, point)*reference(1, point, basis) - &
                    jacobian(2, 1, point)*reference(2, point, basis))/determinant
                physical(2, point, basis) = ( &
                    -jacobian(1, 2, point)*reference(1, point, basis) + &
                    jacobian(1, 1, point)*reference(2, point, basis))/determinant
            end do
        end do
    end subroutine map_covariant_values

    subroutine record(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        call check_condition(condition, message)
        all_passed = all_passed .and. condition
    end subroutine record

end program test_enrichment_support_vector_diagnostics
