program test_tensor_volume_work
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_tensor_volume_work, assemble_tensor_volume_work_jvp, &
        assemble_tensor_volume_work_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: quadrature_count = 2, test_count = 3
    real(dp), parameter :: eps = 1.0e-7_dp
    real(dp) :: test_gradient(quadrature_count, test_count, 3, 3)
    real(dp) :: test_gradient_dot(quadrature_count, test_count, 3, 3)
    real(dp) :: stress(quadrature_count, 3, 3), stress_dot(quadrature_count, 3, 3)
    real(dp) :: weights(quadrature_count), weights_dot(quadrature_count)
    real(dp) :: residual(test_count), residual_dot(test_count)
    real(dp) :: residual_plus(test_count), residual_minus(test_count)
    real(dp) :: residual_bar(test_count)
    real(dp) :: test_gradient_bar(quadrature_count, test_count, 3, 3)
    real(dp) :: stress_bar(quadrature_count, 3, 3), weights_bar(quadrature_count)
    real(dp) :: oracle(test_count), oracle_dot(test_count)
    real(dp) :: lhs, rhs
    real(dp) :: bad_residual(1)
    type(fortsparse_status_t) :: status
    integer :: quadrature, test_function, row, column

    do quadrature = 1, quadrature_count
        weights(quadrature) = 1.1_dp + 0.2_dp*real(quadrature, dp)
        weights_dot(quadrature) = -0.04_dp*real(quadrature, dp)
        do row = 1, 3
            do column = 1, 3
                stress(quadrature, row, column) = &
                    0.07_dp*real(row + 2*column + quadrature, dp)
                stress_dot(quadrature, row, column) = &
                    -0.02_dp*real(2*row + column + quadrature, dp)
                do test_function = 1, test_count
                    test_gradient(quadrature, test_function, row, column) = &
                        0.03_dp*real(test_function + row + column + quadrature, dp)
                    test_gradient_dot(quadrature, test_function, row, column) = &
                        -0.01_dp*real(test_function + 2*row + column, dp)
                end do
            end do
        end do
    end do

    oracle = 0.0_dp
    oracle_dot = 0.0_dp
    do quadrature = 1, quadrature_count
        do test_function = 1, test_count
            do row = 1, 3
                do column = 1, 3
                    oracle(test_function) = oracle(test_function) + &
                        weights(quadrature)*test_gradient( &
                        quadrature, test_function, row, column)* &
                        stress(quadrature, row, column)
                    oracle_dot(test_function) = oracle_dot(test_function) + &
                        weights_dot(quadrature)*test_gradient( &
                        quadrature, test_function, row, column)* &
                        stress(quadrature, row, column) + &
                        weights(quadrature)*test_gradient_dot( &
                        quadrature, test_function, row, column)* &
                        stress(quadrature, row, column) + &
                        weights(quadrature)*test_gradient( &
                        quadrature, test_function, row, column)* &
                        stress_dot(quadrature, row, column)
                end do
            end do
        end do
    end do

    call assemble_tensor_volume_work( &
        test_gradient, stress, weights, residual, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(residual - oracle)) < 1.0e-13_dp, &
        "tensor volume work matches the independent contraction oracle")
    call assemble_tensor_volume_work_jvp( &
        test_gradient, stress, weights, test_gradient_dot, stress_dot, &
        weights_dot, residual_dot, status)
    call assemble_tensor_volume_work( &
        test_gradient + eps*test_gradient_dot, stress + eps*stress_dot, &
        weights + eps*weights_dot, residual_plus, status)
    call assemble_tensor_volume_work( &
        test_gradient - eps*test_gradient_dot, stress - eps*stress_dot, &
        weights - eps*weights_dot, residual_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(residual_dot - oracle_dot)) < 1.0e-13_dp .and. &
        maxval(abs((residual_plus - residual_minus)/(2.0_dp*eps) - residual_dot)) &
        < 1.0e-8_dp, "tensor volume-work JVP matches product rule")

    residual_bar = [0.4_dp, -0.3_dp, 0.2_dp]
    call assemble_tensor_volume_work_vjp( &
        test_gradient, stress, weights, residual_bar, test_gradient_bar, &
        stress_bar, weights_bar, status)
    lhs = dot_product(residual_bar, residual_dot)
    rhs = sum(test_gradient_bar*test_gradient_dot) + sum(stress_bar*stress_dot) + &
        dot_product(weights_bar, weights_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "tensor volume-work VJP satisfies the real dot-product identity")

    call assemble_tensor_volume_work( &
        test_gradient, stress, weights, bad_residual, status)
    call check_condition(status%code /= 0, &
        "tensor volume work rejects an incompatible residual shape")
    call check_summary("tensor volume work")
end program test_tensor_volume_work
