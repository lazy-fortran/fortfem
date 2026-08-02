program test_tensor_diffusion_matrix
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_tensor_diffusion_matrix, &
        assemble_tensor_diffusion_matrix_jvp, assemble_tensor_diffusion_matrix_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: quadrature_count = 2, basis_count = 2
    real(dp), parameter :: finite_difference_step = 1.0e-7_dp
    real(dp) :: basis_gradients(quadrature_count, basis_count, 2)
    real(dp) :: tensor(quadrature_count, 2, 2), weights(quadrature_count)
    real(dp) :: matrix(basis_count, basis_count)
    real(dp) :: gradients_dot(quadrature_count, basis_count, 2)
    real(dp) :: tensor_dot(quadrature_count, 2, 2), weights_dot(quadrature_count)
    real(dp) :: matrix_dot(basis_count, basis_count)
    real(dp) :: gradients_plus(quadrature_count, basis_count, 2)
    real(dp) :: gradients_minus(quadrature_count, basis_count, 2)
    real(dp) :: tensor_plus(quadrature_count, 2, 2), tensor_minus(quadrature_count, 2, 2)
    real(dp) :: weights_plus(quadrature_count), weights_minus(quadrature_count)
    real(dp) :: matrix_plus(basis_count, basis_count), matrix_minus(basis_count, basis_count)
    real(dp) :: gradients_bar(quadrature_count, basis_count, 2)
    real(dp) :: tensor_bar(quadrature_count, 2, 2), weights_bar(quadrature_count)
    real(dp) :: matrix_bar(basis_count, basis_count), lhs, rhs
    real(dp), parameter :: expected_matrix(2, 2) = reshape([ &
        14.25_dp, 1.875_dp, 1.875_dp, 20.0_dp], [2, 2])
    type(fortsparse_status_t) :: status

    basis_gradients = 0.0_dp
    basis_gradients(1, 1, :) = [1.0_dp, 0.0_dp]
    basis_gradients(1, 2, :) = [0.0_dp, 1.0_dp]
    basis_gradients(2, 1, :) = [2.0_dp, 1.0_dp]
    basis_gradients(2, 2, :) = [-1.0_dp, 2.0_dp]
    tensor = 0.0_dp
    tensor(1, :, :) = reshape([4.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
    tensor(2, :, :) = reshape([2.0_dp, -0.3_dp, -0.3_dp, 3.0_dp], [2, 2])
    weights = [0.5_dp, 1.25_dp]
    call assemble_tensor_diffusion_matrix( &
        basis_gradients, tensor, weights, matrix, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(matrix - expected_matrix)) < 1.0e-14_dp, &
        "tensor diffusion matrix matches the independent quadrature oracle")

    gradients_dot = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp, -0.5_dp, 0.6_dp, &
        -0.7_dp, 0.8_dp], shape(gradients_dot))
    tensor_dot = reshape([ &
        0.2_dp, -0.1_dp, 0.05_dp, 0.3_dp, -0.4_dp, 0.15_dp, &
        0.25_dp, -0.35_dp], shape(tensor_dot))
    weights_dot = [0.07_dp, -0.09_dp]
    call assemble_tensor_diffusion_matrix_jvp( &
        basis_gradients, tensor, weights, gradients_dot, tensor_dot, &
        weights_dot, matrix_dot, status)
    gradients_plus = basis_gradients + finite_difference_step*gradients_dot
    gradients_minus = basis_gradients - finite_difference_step*gradients_dot
    tensor_plus = tensor + finite_difference_step*tensor_dot
    tensor_minus = tensor - finite_difference_step*tensor_dot
    weights_plus = weights + finite_difference_step*weights_dot
    weights_minus = weights - finite_difference_step*weights_dot
    call assemble_tensor_diffusion_matrix( &
        gradients_plus, tensor_plus, weights_plus, matrix_plus, status)
    call assemble_tensor_diffusion_matrix( &
        gradients_minus, tensor_minus, weights_minus, matrix_minus, status)
    call check_condition(maxval(abs(matrix_dot - &
        (matrix_plus - matrix_minus)/(2.0_dp*finite_difference_step))) < 1.0e-8_dp, &
        "tensor diffusion matrix JVP matches a central difference")

    matrix_bar = reshape([0.6_dp, -0.4_dp, 0.3_dp, 0.8_dp], [2, 2])
    call assemble_tensor_diffusion_matrix_vjp( &
        basis_gradients, tensor, weights, matrix_bar, gradients_bar, tensor_bar, &
        weights_bar, status)
    lhs = sum(matrix_bar*matrix_dot)
    rhs = sum(gradients_bar*gradients_dot) + sum(tensor_bar*tensor_dot) + &
        dot_product(weights_bar, weights_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "tensor diffusion matrix VJP satisfies the real dot-product identity")

    call assemble_tensor_diffusion_matrix( &
        basis_gradients, tensor, [-0.5_dp, 1.25_dp], matrix, status)
    call check_condition(status%code /= 0, &
        "tensor diffusion matrix rejects a non-positive quadrature weight")
    call check_summary("tensor diffusion matrix")
end program test_tensor_diffusion_matrix
