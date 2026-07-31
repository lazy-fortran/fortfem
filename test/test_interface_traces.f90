program test_interface_traces
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        compute_interface_scalar_jump_average, &
        compute_interface_vector_traces
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp) :: scalar_average, scalar_jump
    real(dp), parameter :: plus_vector(3) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: minus_vector(3) = [-1.0_dp, 4.0_dp, 5.0_dp]
    real(dp), parameter :: normal(3) = [0.0_dp, 0.0_dp, 1.0_dp]
    real(dp), parameter :: bad_normal(3) = [0.0_dp, 0.0_dp, 2.0_dp]
    real(dp) :: plus_normal, minus_normal, normal_average, normal_jump
    real(dp) :: plus_tangent(3), minus_tangent(3)
    real(dp) :: tangent_average(3), tangent_jump(3), rotated_jump(3)
    type(fortsparse_status_t) :: status

    call compute_interface_scalar_jump_average( &
        3.0_dp, -1.0_dp, scalar_average, scalar_jump, status)
    call check_condition(status%code == 0 .and. &
        abs(scalar_average - 1.0_dp) < 1.0e-14_dp .and. &
        abs(scalar_jump - 4.0_dp) < 1.0e-14_dp, &
        "interface scalar average and jump use plus-minus orientation")

    call compute_interface_vector_traces( &
        plus_vector, minus_vector, normal, plus_normal, minus_normal, &
        normal_average, normal_jump, plus_tangent, minus_tangent, &
        tangent_average, tangent_jump, rotated_jump, status)
    call check_condition(status%code == 0, &
        "interface vector traces accept a unit oriented normal")
    call check_condition(abs(plus_normal - 3.0_dp) < 1.0e-14_dp .and. &
        abs(minus_normal - 5.0_dp) < 1.0e-14_dp .and. &
        abs(normal_average - 4.0_dp) < 1.0e-14_dp .and. &
        abs(normal_jump + 2.0_dp) < 1.0e-14_dp, &
        "interface normal traces and jump match the independent oracle")
    call check_condition(maxval(abs(plus_tangent - [1.0_dp, 2.0_dp, 0.0_dp])) < &
        1.0e-14_dp .and. maxval(abs(minus_tangent - [-1.0_dp, 4.0_dp, 0.0_dp])) < &
        1.0e-14_dp, "interface tangential projections remove normal traces")
    call check_condition(maxval(abs(tangent_average - [0.0_dp, 3.0_dp, 0.0_dp])) < &
        1.0e-14_dp .and. maxval(abs(tangent_jump - [2.0_dp, -2.0_dp, 0.0_dp])) < &
        1.0e-14_dp, "interface tangential average and jump match the oracle")
    call check_condition(maxval(abs(rotated_jump - [2.0_dp, 2.0_dp, 0.0_dp])) < &
        1.0e-14_dp .and. abs(dot_product(normal, rotated_jump)) < 1.0e-14_dp, &
        "rotated tangential jump gives the oriented Ampere surface-current sign")

    call compute_interface_vector_traces( &
        plus_vector, minus_vector, bad_normal, plus_normal, minus_normal, &
        normal_average, normal_jump, plus_tangent, minus_tangent, &
        tangent_average, tangent_jump, rotated_jump, status)
    call check_condition(status%code /= 0, &
        "interface vector traces reject a non-unit normal")
    call check_summary("interface traces")
end program test_interface_traces
