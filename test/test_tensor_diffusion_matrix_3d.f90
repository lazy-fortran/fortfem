program test_tensor_diffusion_matrix_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_tensor_diffusion_matrix_3d, &
        assemble_tensor_diffusion_matrix_3d_jvp, &
        assemble_tensor_diffusion_matrix_3d_vjp, &
        assemble_tensor_diffusion_matrix_nd
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: quadrature_count = 2, basis_count = 2
    real(dp), parameter :: finite_difference_step = 1.0e-7_dp
    real(dp) :: basis_gradients(quadrature_count, basis_count, 3)
    real(dp) :: tensor(quadrature_count, 3, 3), weights(quadrature_count)
    real(dp) :: matrix(basis_count, basis_count), matrix_nd(basis_count, basis_count)
    real(dp) :: gradients_dot(quadrature_count, basis_count, 3)
    real(dp) :: tensor_dot(quadrature_count, 3, 3), weights_dot(quadrature_count)
    real(dp) :: matrix_dot(basis_count, basis_count)
    real(dp) :: gradients_plus(quadrature_count, basis_count, 3)
    real(dp) :: gradients_minus(quadrature_count, basis_count, 3)
    real(dp) :: tensor_plus(quadrature_count, 3, 3)
    real(dp) :: tensor_minus(quadrature_count, 3, 3)
    real(dp) :: weights_plus(quadrature_count), weights_minus(quadrature_count)
    real(dp) :: matrix_plus(basis_count, basis_count)
    real(dp) :: matrix_minus(basis_count, basis_count)
    real(dp) :: gradients_bar(quadrature_count, basis_count, 3)
    real(dp) :: tensor_bar(quadrature_count, 3, 3), weights_bar(quadrature_count)
    real(dp) :: matrix_bar(basis_count, basis_count), lhs, rhs
    real(dp) :: basis_gradients_2d(quadrature_count, basis_count, 2)
    real(dp) :: tensor_2d(quadrature_count, 2, 2)
    real(dp), parameter :: expected_matrix(2, 2) = reshape([ &
        15.89875_dp, 1.16875_dp, 1.16875_dp, 12.506875_dp], [2, 2])
    type(fortsparse_status_t) :: status

    basis_gradients = 0.0_dp
    basis_gradients(1, 1, :) = [1.0_dp, 2.0_dp, 0.5_dp]
    basis_gradients(1, 2, :) = [-1.0_dp, 0.0_dp, 2.0_dp]
    basis_gradients(2, 1, :) = [0.5_dp, -1.0_dp, 1.5_dp]
    basis_gradients(2, 2, :) = [2.0_dp, 0.25_dp, -0.5_dp]
    tensor(1, :, :) = reshape([ &
        3.0_dp, 0.2_dp, 0.1_dp, 0.2_dp, 2.0_dp, 0.3_dp, &
        0.1_dp, 0.3_dp, 1.5_dp], [3, 3])
    tensor(2, :, :) = reshape([ &
        1.5_dp, -0.1_dp, 0.4_dp, -0.1_dp, 2.5_dp, 0.2_dp, &
        0.4_dp, 0.2_dp, 1.2_dp], [3, 3])
    weights = [0.75_dp, 1.1_dp]

    call assemble_tensor_diffusion_matrix_3d( &
        basis_gradients, tensor, weights, matrix, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(matrix - expected_matrix)) < 1.0e-13_dp, &
        "3D tensor diffusion matrix matches an independent oracle")
    call assemble_tensor_diffusion_matrix_nd( &
        basis_gradients, tensor, weights, matrix_nd, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(matrix_nd - expected_matrix)) < 1.0e-13_dp, &
        "general-dimension tensor diffusion matches the 3D wrapper")

    gradients_dot = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp, -0.5_dp, 0.6_dp, &
        -0.7_dp, 0.8_dp, -0.9_dp, 1.0_dp, -1.1_dp, 1.2_dp], &
        shape(gradients_dot))
    tensor_dot = reshape([ &
        0.2_dp, -0.1_dp, 0.05_dp, 0.3_dp, -0.4_dp, 0.15_dp, &
        0.25_dp, -0.35_dp, 0.45_dp, -0.55_dp, 0.65_dp, -0.75_dp, &
        0.85_dp, -0.95_dp, 1.05_dp, 1.15_dp, -1.25_dp, 1.35_dp], &
        shape(tensor_dot))
    weights_dot = [0.07_dp, -0.09_dp]
    call assemble_tensor_diffusion_matrix_3d_jvp( &
        basis_gradients, tensor, weights, gradients_dot, tensor_dot, &
        weights_dot, matrix_dot, status)
    gradients_plus = basis_gradients + finite_difference_step*gradients_dot
    gradients_minus = basis_gradients - finite_difference_step*gradients_dot
    tensor_plus = tensor + finite_difference_step*tensor_dot
    tensor_minus = tensor - finite_difference_step*tensor_dot
    weights_plus = weights + finite_difference_step*weights_dot
    weights_minus = weights - finite_difference_step*weights_dot
    call assemble_tensor_diffusion_matrix_3d( &
        gradients_plus, tensor_plus, weights_plus, matrix_plus, status)
    call assemble_tensor_diffusion_matrix_3d( &
        gradients_minus, tensor_minus, weights_minus, matrix_minus, status)
    call check_condition(maxval(abs(matrix_dot - &
        (matrix_plus - matrix_minus)/(2.0_dp*finite_difference_step))) < 1.0e-8_dp, &
        "3D tensor diffusion JVP matches a central difference")

    matrix_bar = reshape([0.6_dp, -0.4_dp, 0.3_dp, 0.8_dp], [2, 2])
    call assemble_tensor_diffusion_matrix_3d_vjp( &
        basis_gradients, tensor, weights, matrix_bar, gradients_bar, tensor_bar, &
        weights_bar, status)
    lhs = sum(matrix_bar*matrix_dot)
    rhs = sum(gradients_bar*gradients_dot) + sum(tensor_bar*tensor_dot) + &
        dot_product(weights_bar, weights_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "3D tensor diffusion VJP satisfies the real dot-product identity")

    basis_gradients_2d = 0.0_dp
    tensor_2d = 0.0_dp
    call assemble_tensor_diffusion_matrix_3d( &
        basis_gradients_2d, tensor_2d, weights, matrix, status)
    call check_condition(status%code /= 0, &
        "3D tensor diffusion rejects a two-dimensional tensor")
    call check_summary("3D tensor diffusion matrix")
end program test_tensor_diffusion_matrix_3d
