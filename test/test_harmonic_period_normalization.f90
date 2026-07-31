program test_harmonic_period_normalization
    use check, only: check_condition, check_summary
    use fortfem_api, only: normalize_harmonic_one_forms, &
        normalize_harmonic_one_forms_jvp, normalize_harmonic_one_forms_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: harmonic(3, 2) = reshape([ &
        2.0_dp, 0.25_dp, 1.0_dp, 0.5_dp, 1.5_dp, -0.2_dp], [3, 2])
    real(dp), parameter :: cycles(3, 2) = reshape([ &
        1.0_dp, 0.0_dp, 0.3_dp, 0.0_dp, 1.0_dp, -0.2_dp], [3, 2])
    real(dp), parameter :: target(2, 2) = reshape([ &
        1.0_dp, -0.1_dp, 0.2_dp, 1.1_dp], [2, 2])
    real(dp), parameter :: harmonic_dot(3, 2) = reshape([ &
        0.03_dp, -0.02_dp, 0.04_dp, -0.01_dp, 0.02_dp, 0.05_dp], [3, 2])
    real(dp), parameter :: cycles_dot(3, 2) = reshape([ &
        -0.01_dp, 0.02_dp, 0.03_dp, 0.04_dp, -0.02_dp, 0.01_dp], [3, 2])
    real(dp), parameter :: target_dot(2, 2) = reshape([ &
        0.02_dp, -0.03_dp, 0.01_dp, 0.04_dp], [2, 2])
    real(dp), parameter :: normalized_bar(3, 2) = reshape([ &
        0.4_dp, -0.3_dp, 0.2_dp, -0.1_dp, 0.5_dp, 0.6_dp], [3, 2])
    real(dp), parameter :: normalization_bar(2, 2) = reshape([ &
        -0.2_dp, 0.7_dp, 0.3_dp, -0.4_dp], [2, 2])
    real(dp), parameter :: eps = 1.0e-7_dp
    real(dp) :: normalized(3, 2), normalization(2, 2)
    real(dp) :: normalized_dot(3, 2), normalization_dot(2, 2)
    real(dp) :: normalized_plus(3, 2), normalized_minus(3, 2)
    real(dp) :: normalization_plus(2, 2), normalization_minus(2, 2)
    real(dp) :: harmonic_bar(3, 2), cycles_bar(3, 2), target_bar(2, 2)
    real(dp) :: period_matrix(2, 2), expected_normalization(2, 2)
    real(dp) :: determinant, lhs, rhs
    real(dp) :: bad_normalized(1, 2)
    integer :: status

    call normalize_harmonic_one_forms( &
        harmonic, cycles, target, normalized, normalization, status)
    period_matrix = matmul(transpose(cycles), harmonic)
    determinant = period_matrix(1, 1)*period_matrix(2, 2) - &
        period_matrix(1, 2)*period_matrix(2, 1)
    expected_normalization = reshape([ &
        period_matrix(2, 2), -period_matrix(2, 1), -period_matrix(1, 2), &
        period_matrix(1, 1)], [2, 2])/determinant
    call check_condition(status == 0 .and. &
        maxval(abs(normalization - matmul(expected_normalization, target))) < &
        1.0e-13_dp .and. maxval(abs(matmul(transpose(cycles), normalized) - &
        target)) < 1.0e-13_dp, &
        "harmonic period normalization matches the independent inverse oracle")

    call normalize_harmonic_one_forms_jvp( &
        harmonic, cycles, target, harmonic_dot, cycles_dot, target_dot, &
        normalized_dot, normalization_dot, status)
    call normalize_harmonic_one_forms( &
        harmonic + eps*harmonic_dot, cycles + eps*cycles_dot, &
        target + eps*target_dot, normalized_plus, normalization_plus, status)
    call normalize_harmonic_one_forms( &
        harmonic - eps*harmonic_dot, cycles - eps*cycles_dot, &
        target - eps*target_dot, normalized_minus, normalization_minus, status)
    call check_condition(status == 0 .and. &
        maxval(abs(normalized_dot - (normalized_plus - normalized_minus)/ &
        (2.0_dp*eps))) &
        < 1.0e-8_dp .and. maxval(abs(normalization_dot - &
        (normalization_plus - normalization_minus)/(2.0_dp*eps))) < 1.0e-8_dp, &
        "harmonic period normalization JVP matches central differences")

    call normalize_harmonic_one_forms_vjp( &
        harmonic, cycles, target, normalized, normalization, normalized_bar, &
        normalization_bar, harmonic_bar, cycles_bar, target_bar, status)
    lhs = sum(normalized_bar*normalized_dot) + &
        sum(normalization_bar*normalization_dot)
    rhs = sum(harmonic_bar*harmonic_dot) + sum(cycles_bar*cycles_dot) + &
        sum(target_bar*target_dot)
    call check_condition(status == 0 .and. abs(lhs - rhs) < 1.0e-11_dp, &
        "harmonic period normalization VJP satisfies the real dot identity")

    call normalize_harmonic_one_forms( &
        harmonic, cycles, target, bad_normalized, normalization, status)
    call check_condition(status /= 0, &
        "harmonic period normalization rejects an incompatible output shape")
    call check_summary("Harmonic period normalization")
end program test_harmonic_period_normalization
