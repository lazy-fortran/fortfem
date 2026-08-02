program test_surface_vector_trace
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_surface_vector_trace, only: &
        evaluate_surface_vector_trace, &
        evaluate_surface_vector_trace_jvp, &
        evaluate_surface_vector_trace_vjp
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    implicit none

    integer, parameter :: sample_count = 3
    real(dp), parameter :: epsilon = 1.0e-7_dp, tolerance = 3.0e-8_dp
    real(dp) :: field(3, sample_count), normals(3, sample_count)
    real(dp) :: field_dot(3, sample_count), normals_dot(3, sample_count)
    real(dp) :: field_plus(3, sample_count), field_minus(3, sample_count)
    real(dp) :: normals_plus(3, sample_count), normals_minus(3, sample_count)
    real(dp) :: normal_component(sample_count), tangential_field(3, sample_count)
    real(dp) :: normal_component_dot(sample_count)
    real(dp) :: tangential_field_dot(3, sample_count)
    real(dp) :: normal_component_plus(sample_count)
    real(dp) :: normal_component_minus(sample_count)
    real(dp) :: tangential_field_plus(3, sample_count)
    real(dp) :: tangential_field_minus(3, sample_count)
    real(dp) :: normal_component_bar(sample_count)
    real(dp) :: tangential_field_bar(3, sample_count)
    real(dp) :: field_bar(3, sample_count), normals_bar(3, sample_count)
    real(dp) :: lhs, rhs, expected_normal, normal_length
    type(fortsparse_status_t) :: status
    logical :: all_passed
    integer :: sample

    all_passed = .true.
    field = reshape([ &
        1.0_dp, 2.0_dp, 3.0_dp, &
        -2.0_dp, 0.5_dp, 4.0_dp, &
        0.25_dp, -1.0_dp, 2.5_dp], shape(field))
    normals = reshape([ &
        3.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, -4.0_dp, 0.0_dp, &
        1.0_dp, 2.0_dp, 2.0_dp], shape(normals))
    field_dot = reshape([ &
        0.2_dp, -0.3_dp, 0.1_dp, &
        0.5_dp, 0.4_dp, -0.2_dp, &
        -0.1_dp, 0.7_dp, 0.6_dp], shape(field_dot))
    normals_dot = reshape([ &
        0.1_dp, 0.2_dp, -0.3_dp, &
        -0.2_dp, 0.3_dp, 0.4_dp, &
        0.5_dp, -0.1_dp, 0.2_dp], shape(normals_dot))

    call evaluate_surface_vector_trace( &
        field, normals, normal_component, tangential_field, status)
    call record(status%code == FORTSPARSE_OK, &
        "surface trace accepts non-unit caller-owned normals")
    do sample = 1, sample_count
        normal_length = sqrt(dot_product(normals(:, sample), normals(:, sample)))
        expected_normal = dot_product(normals(:, sample), field(:, sample))/ &
            normal_length
        call record(abs(normal_component(sample) - expected_normal) < tolerance .and. &
            abs(dot_product(normals(:, sample), tangential_field(:, sample))) < &
            tolerance .and. &
            maxval(abs(field(:, sample) - (tangential_field(:, sample) + &
            normal_component(sample)*normals(:, sample)/normal_length))) < tolerance, &
            "normal and tangential traces reconstruct each physical field sample")
    end do

    call evaluate_surface_vector_trace_jvp( &
        field, normals, field_dot, normals_dot, normal_component_dot, &
        tangential_field_dot, status)
    field_plus = field + epsilon*field_dot
    field_minus = field - epsilon*field_dot
    normals_plus = normals + epsilon*normals_dot
    normals_minus = normals - epsilon*normals_dot
    call evaluate_surface_vector_trace( &
        field_plus, normals_plus, normal_component_plus, tangential_field_plus, status)
    call evaluate_surface_vector_trace( &
        field_minus, normals_minus, normal_component_minus, &
        tangential_field_minus, status)
    call record(maxval(abs(normal_component_dot - (normal_component_plus - &
        normal_component_minus)/(2.0_dp*epsilon))) < tolerance .and. &
        maxval(abs(tangential_field_dot - (tangential_field_plus - &
        tangential_field_minus)/(2.0_dp*epsilon))) < tolerance, &
        "surface trace JVP matches a complete central difference")

    normal_component_bar = [0.4_dp, -0.7_dp, 0.2_dp]
    tangential_field_bar = reshape([ &
        0.1_dp, 0.3_dp, -0.2_dp, &
        -0.4_dp, 0.2_dp, 0.6_dp, &
        0.5_dp, -0.1_dp, 0.7_dp], shape(tangential_field_bar))
    call evaluate_surface_vector_trace_vjp( &
        field, normals, normal_component_bar, tangential_field_bar, &
        field_bar, normals_bar, status)
    lhs = dot_product(normal_component_bar, normal_component_dot) + &
        sum(tangential_field_bar*tangential_field_dot)
    rhs = sum(field_bar*field_dot) + sum(normals_bar*normals_dot)
    call record(status%code == FORTSPARSE_OK .and. abs(lhs - rhs) < tolerance, &
        "surface trace VJP satisfies the real dot-product oracle")

    call evaluate_surface_vector_trace( &
        field, reshape([3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 2.0_dp, 2.0_dp], [3, sample_count]), normal_component, &
        tangential_field, status)
    call record(status%code /= FORTSPARSE_OK, "zero surface normal is rejected")

    call evaluate_surface_vector_trace( &
        field(:2, :), normals, normal_component, tangential_field, status)
    call record(status%code /= FORTSPARSE_OK, &
        "non-three-dimensional vectors are rejected")

    call check_summary("Surface vector trace")
    if (.not. all_passed) error stop 1

contains

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_surface_vector_trace
