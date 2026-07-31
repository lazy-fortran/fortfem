program test_surface_current
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_interface_surface_current, &
        assemble_interface_surface_current_jvp, &
        assemble_interface_surface_current_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: plus_field(2, 3) = reshape([ &
        2.0_dp, 1.0_dp, 3.0_dp, 0.0_dp, 4.0_dp, 2.0_dp], [2, 3])
    real(dp), parameter :: minus_field(2, 3) = reshape([ &
        1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 4.0_dp, 0.0_dp], [2, 3])
    real(dp), parameter :: normals(2, 3) = reshape([ &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], [2, 3])
    real(dp), parameter :: weights(2) = [2.0_dp, 3.0_dp]
    real(dp), parameter :: expected_current(2, 3) = reshape([ &
        -2.0_dp, 0.0_dp, 1.0_dp, -2.0_dp, 0.0_dp, -2.0_dp], [2, 3])
    real(dp), parameter :: expected_integrated(3) = [-4.0_dp, -4.0_dp, -6.0_dp]
    real(dp), parameter :: plus_dot(2, 3) = reshape([ &
        0.2_dp, -0.4_dp, -0.1_dp, 0.5_dp, 0.3_dp, 0.2_dp], [2, 3])
    real(dp), parameter :: minus_dot(2, 3) = reshape([ &
        -0.3_dp, 0.1_dp, 0.4_dp, -0.2_dp, -0.2_dp, 0.3_dp], [2, 3])
    real(dp), parameter :: normals_dot(2, 3) = reshape([ &
        0.1_dp, 0.0_dp, -0.2_dp, 0.2_dp, 0.0_dp, -0.1_dp], [2, 3])
    real(dp), parameter :: weights_dot(2) = [0.2_dp, -0.3_dp]
    real(dp), parameter :: current_bar(2, 3) = reshape([ &
        0.7_dp, -0.3_dp, -0.4_dp, 0.5_dp, 0.2_dp, 0.9_dp], [2, 3])
    real(dp), parameter :: integrated_bar(3) = [0.6_dp, -0.8_dp, 0.4_dp]
    real(dp), parameter :: eps = 1.0e-6_dp
    real(dp) :: current(2, 3), integrated(3)
    real(dp) :: current_dot(2, 3), integrated_dot(3)
    real(dp) :: current_plus(2, 3), current_minus(2, 3)
    real(dp) :: integrated_plus(3), integrated_minus(3)
    real(dp) :: plus_bar(2, 3), minus_bar(2, 3), normals_bar(2, 3)
    real(dp) :: weights_bar(2), lhs, rhs
    real(dp) :: reversed_current(2, 3), reversed_integrated(3)
    real(dp), parameter :: bad_normals(2, 3) = reshape([ &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp], [2, 3])
    real(dp), parameter :: bad_weights(2) = [2.0_dp, -1.0_dp]
    type(fortsparse_status_t) :: status

    call assemble_interface_surface_current( &
        plus_field, minus_field, normals, weights, current, integrated, status)
    call check_condition(status%code == 0, &
        "surface current accepts oriented unit-normal traces")
    call check_condition(maxval(abs(current - expected_current)) < 1.0e-14_dp .and. &
        maxval(abs(integrated - expected_integrated)) < 1.0e-14_dp, &
        "Ampere rotated jump and integrated ledger match the independent oracle")

    call assemble_interface_surface_current( &
        minus_field, plus_field, -normals, weights, reversed_current, &
        reversed_integrated, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(reversed_current - current)) < 1.0e-14_dp .and. &
        maxval(abs(reversed_integrated - integrated)) < 1.0e-14_dp, &
        "reversing interface orientation preserves the physical sheet current")

    call assemble_interface_surface_current_jvp( &
        plus_field, minus_field, normals, weights, plus_dot, minus_dot, &
        normals_dot, weights_dot, current_dot, integrated_dot, status)
    call assemble_interface_surface_current( &
        plus_field + eps*plus_dot, minus_field + eps*minus_dot, &
        normals + eps*normals_dot, weights + eps*weights_dot, current_plus, &
        integrated_plus, status)
    call assemble_interface_surface_current( &
        plus_field - eps*plus_dot, minus_field - eps*minus_dot, &
        normals - eps*normals_dot, weights - eps*weights_dot, current_minus, &
        integrated_minus, status)
    call check_condition(maxval(abs((current_plus - current_minus)/(2.0_dp*eps) - &
        current_dot)) < 1.0e-8_dp .and. &
        maxval(abs((integrated_plus - integrated_minus)/(2.0_dp*eps) - &
        integrated_dot)) < 1.0e-8_dp, &
        "surface-current JVP matches independent central differences")

    call assemble_interface_surface_current_vjp( &
        plus_field, minus_field, normals, weights, current_bar, integrated_bar, &
        plus_bar, minus_bar, normals_bar, weights_bar, status)
    lhs = sum(current_bar*current_dot) + dot_product(integrated_bar, integrated_dot)
    rhs = sum(plus_bar*plus_dot) + sum(minus_bar*minus_dot) + &
        sum(normals_bar*normals_dot) + dot_product(weights_bar, weights_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "surface-current VJP satisfies the real dot-product identity")

    call assemble_interface_surface_current( &
        plus_field, minus_field, bad_normals, weights, current, integrated, status)
    call check_condition(status%code /= 0, &
        "surface current rejects non-unit normals")
    call assemble_interface_surface_current( &
        plus_field, minus_field, normals, bad_weights, current, integrated, status)
    call check_condition(status%code /= 0, &
        "surface current rejects non-positive surface weights")

    call check_summary("surface current")
end program test_surface_current
