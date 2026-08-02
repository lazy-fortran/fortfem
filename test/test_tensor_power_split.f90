program test_tensor_power_split
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_tensor_power_split, &
        evaluate_tensor_power_split_jvp, evaluate_tensor_power_split_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: tensor(3, 3) = reshape([ &
        4.0_dp, 1.5_dp, -0.4_dp, 0.5_dp, 2.0_dp, 0.8_dp, &
        0.4_dp, -0.2_dp, 3.0_dp], [3, 3])
    real(dp), parameter :: vector(3) = [1.2_dp, -0.7_dp, 0.9_dp]
    real(dp), parameter :: tensor_dot(3, 3) = reshape([ &
        -0.2_dp, 0.3_dp, 0.1_dp, -0.4_dp, 0.5_dp, -0.6_dp, &
        0.2_dp, 0.7_dp, -0.1_dp], [3, 3])
    real(dp), parameter :: vector_dot(3) = [0.25_dp, -0.1_dp, 0.3_dp]
    real(dp), parameter :: power_bar(3) = [0.7_dp, -0.2_dp, 0.4_dp]
    real(dp), parameter :: step = 1.0e-7_dp
    real(dp) :: symmetric_power, skew_power, total_power
    real(dp) :: symmetric_power_dot, skew_power_dot, total_power_dot
    real(dp) :: symmetric_plus, skew_plus, total_plus
    real(dp) :: symmetric_minus, skew_minus, total_minus
    real(dp) :: tensor_bar(3, 3), vector_bar(3), left, right
    real(dp) :: symmetric_tensor(3, 3), skew_tensor(3, 3)
    real(dp) :: bad_tensor(3, 3)
    type(fortsparse_status_t) :: status

    symmetric_tensor = 0.5_dp*(tensor + transpose(tensor))
    skew_tensor = 0.5_dp*(tensor - transpose(tensor))
    call evaluate_tensor_power_split( &
        tensor, vector, symmetric_power, skew_power, total_power, status)
    call check_condition(status%code == 0, &
        "tensor power split accepts finite three-dimensional data")
    call check_condition(abs(symmetric_power - dot_product(vector, &
        matmul(symmetric_tensor, vector))) < 1.0e-14_dp .and. &
        abs(skew_power - dot_product(vector, matmul(skew_tensor, vector))) < &
        1.0e-14_dp .and. abs(total_power - dot_product(vector, &
        matmul(tensor, vector))) < 1.0e-14_dp, &
        "symmetric, skew, and total powers match an independent oracle")
    call check_condition(abs(skew_power) < 1.0e-14_dp .and. &
        abs(total_power - symmetric_power) < 1.0e-14_dp, &
        "a skew constitutive part contributes no instantaneous power")

    call evaluate_tensor_power_split_jvp( &
        tensor, vector, tensor_dot, vector_dot, symmetric_power_dot, &
        skew_power_dot, total_power_dot, status)
    call evaluate_tensor_power_split( &
        tensor + step*tensor_dot, vector + step*vector_dot, symmetric_plus, &
        skew_plus, total_plus, status)
    call evaluate_tensor_power_split( &
        tensor - step*tensor_dot, vector - step*vector_dot, symmetric_minus, &
        skew_minus, total_minus, status)
    call check_condition(maxval(abs([symmetric_power_dot, skew_power_dot, &
        total_power_dot] - [ (symmetric_plus - symmetric_minus)/(2.0_dp*step), &
        (skew_plus - skew_minus)/(2.0_dp*step), &
        (total_plus - total_minus)/(2.0_dp*step) ])) < 3.0e-8_dp, &
        "tensor power split JVP matches an independent central difference")

    call evaluate_tensor_power_split_vjp( &
        tensor, vector, power_bar(1), power_bar(2), power_bar(3), tensor_bar, &
        vector_bar, status)
    left = power_bar(1)*symmetric_power_dot + power_bar(2)*skew_power_dot + &
        power_bar(3)*total_power_dot
    right = sum(tensor_bar*tensor_dot) + dot_product(vector_bar, vector_dot)
    call check_condition(status%code == 0 .and. abs(left - right) < 2.0e-13_dp, &
        "tensor power split VJP satisfies the real transpose oracle")

    bad_tensor = ieee_value(0.0_dp, ieee_quiet_nan)
    call evaluate_tensor_power_split( &
        bad_tensor, vector, symmetric_power, skew_power, total_power, status)
    call check_condition(status%code /= 0, &
        "tensor power split rejects non-finite constitutive data")
    call check_summary("tensor power split")
end program test_tensor_power_split
