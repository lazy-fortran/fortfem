program test_enrichment_support_tensor_diagnostics
    use check, only: check_condition, check_summary
    use fortfem_enrichment_support_tensor_diagnostics, only: &
        evaluate_enrichment_support_tensor_gram, &
        evaluate_enrichment_support_tensor_gram_jvp, &
        evaluate_enrichment_support_tensor_gram_vjp, &
        evaluate_enrichment_support_tensor_rank_condition
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp) :: basis(3, 3, 2), basis_dot(3, 3, 2), metric(2, 2, 3), metric_dot(2, 2, 3)
    real(dp) :: weights(3), weights_dot(3), gram(3, 3), gram_dot(3, 3)
    real(dp) :: gram_plus(3, 3), gram_minus(3, 3), gram_bar(3, 3)
    real(dp) :: basis_bar(3, 3, 2), metric_bar(2, 2, 3), weights_bar(3)
    real(dp) :: lhs, rhs, step, energy_plus, energy_minus
    logical :: active(3)
    integer :: rank
    real(dp) :: condition_estimate, eigenvalues(2)
    type(fortsparse_status_t) :: status

    step = 1.0e-7_dp
    weights = [0.5_dp, 1.0_dp, 1.5_dp]
    weights_dot = [0.1_dp, -0.05_dp, 0.2_dp]
    active = [.true., .false., .true.]
    basis = 0.0_dp
    basis(1, 1, :) = [1.0_dp, 0.0_dp]
    basis(2, 1, :) = [0.0_dp, 1.0_dp]
    basis(3, 1, :) = [1.0_dp, 1.0_dp]
    basis(:, 3, :) = 2.0_dp*basis(:, 1, :)
    basis(:, 2, :) = reshape([0.3_dp, -0.2_dp, 0.4_dp, 0.1_dp, -0.5_dp, 0.6_dp], [3, 2])
    basis_dot = reshape([0.05_dp, -0.1_dp, 0.2_dp, 0.1_dp, -0.04_dp, 0.03_dp, &
        0.07_dp, -0.02_dp, 0.06_dp, 0.08_dp, -0.03_dp, 0.09_dp, 0.02_dp, 0.05_dp, &
        -0.06_dp, 0.01_dp, 0.04_dp, -0.02_dp], [3, 3, 2])
    metric(:, :, 1) = reshape([2.0_dp, 0.2_dp, 0.2_dp, 1.5_dp], [2, 2])
    metric(:, :, 2) = reshape([1.5_dp, -0.1_dp, -0.1_dp, 2.0_dp], [2, 2])
    metric(:, :, 3) = reshape([1.2_dp, 0.15_dp, 0.15_dp, 1.8_dp], [2, 2])
    metric_dot(:, :, 1) = reshape([0.1_dp, 0.02_dp, 0.02_dp, -0.05_dp], [2, 2])
    metric_dot(:, :, 2) = reshape([-0.03_dp, 0.01_dp, 0.01_dp, 0.04_dp], [2, 2])
    metric_dot(:, :, 3) = reshape([0.02_dp, -0.01_dp, -0.01_dp, 0.05_dp], [2, 2])

    call evaluate_enrichment_support_tensor_gram(basis, metric, weights, active, gram, status)
    call check_condition(status%code == 0, "point-first tensor Gram accepts active enrichment data")
    gram_bar = 0.0_dp
    call independent_gram(basis, metric, weights, active, gram_bar)
    call check_condition(maxval(abs(gram - gram_bar)) < 1.0e-14_dp, &
        "tensor Gram matches an independent nested-loop oracle")

    call evaluate_enrichment_support_tensor_gram_jvp(basis, metric, weights, active, basis_dot, &
        metric_dot, weights_dot, gram_dot, status)
    call evaluate_enrichment_support_tensor_gram(basis + step*basis_dot, metric + step*metric_dot, &
        weights + step*weights_dot, active, gram_plus, status)
    call evaluate_enrichment_support_tensor_gram(basis - step*basis_dot, metric - step*metric_dot, &
        weights - step*weights_dot, active, gram_minus, status)
    call check_condition(maxval(abs(gram_dot - (gram_plus - gram_minus)/(2.0_dp*step))) < 1.0e-7_dp, &
        "tensor Gram JVP matches central differences")

    gram_bar = reshape([0.4_dp, -0.2_dp, 0.1_dp, 0.3_dp, 0.5_dp, -0.6_dp, &
        -0.7_dp, 0.8_dp, 0.2_dp], [3, 3])
    call evaluate_enrichment_support_tensor_gram_vjp(basis, metric, weights, active, gram_bar, basis_bar, &
        metric_bar, weights_bar, status)
    lhs = sum(gram_bar*gram_dot)
    rhs = sum(basis_bar*basis_dot) + sum(metric_bar*metric_dot) + dot_product(weights_bar, weights_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "tensor Gram VJP satisfies the real transpose oracle")

    call evaluate_enrichment_support_tensor_rank_condition(gram, active, 1.0e-10_dp, rank, &
        condition_estimate, eigenvalues, status)
    call check_condition(status%code == 0 .and. rank == 1 .and. condition_estimate > 1.0e100_dp, &
        "rank diagnostic detects a linearly dependent active enrichment")

    call evaluate_enrichment_support_tensor_gram(basis, metric, [0.5_dp, 0.0_dp, 1.5_dp], active, gram, status)
    call check_condition(status%code /= 0, "tensor Gram rejects nonpositive sample weights")
    call check_summary("enrichment support tensor diagnostics")

contains

    subroutine independent_gram(values, metric_values, sample_weights, mask, result)
        real(dp), intent(in) :: values(:, :, :), metric_values(:, :, :), sample_weights(:)
        logical, intent(in) :: mask(:)
        real(dp), intent(out) :: result(:, :)
        integer :: point, first, second

        result = 0.0_dp
        do first = 1, size(mask)
            if (.not. mask(first)) cycle
            do second = 1, size(mask)
                if (.not. mask(second)) cycle
                do point = 1, size(sample_weights)
                    result(first, second) = result(first, second) + sample_weights(point)* &
                        sum(values(point, first, :)*matmul(metric_values(:, :, point), &
                        values(point, second, :)))
                end do
            end do
        end do
    end subroutine independent_gram

end program test_enrichment_support_tensor_diagnostics
