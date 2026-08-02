program test_surface_shape_objective_properties
    use, intrinsic :: iso_fortran_env, only: int32
    use check, only: check_condition, check_property, check_summary, &
        property_random_integer, property_random_unit, property_rng_initialize, &
        property_rng_t
    use fortfem_api, only: &
        evaluate_surface_shape_objective, &
        evaluate_surface_shape_objective_jvp, &
        evaluate_surface_shape_objective_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    logical :: all_passed
    integer(int32) :: first_failed_seed, shrunk_seed

    call check_property( &
        "random surface-shape value/JVP/VJP and validation identities", &
        20260802_int32, 24, surface_shape_case, all_passed, first_failed_seed, &
        shrunk_seed)
    call check_condition(all_passed .and. first_failed_seed == 0_int32 .and. &
        shrunk_seed == 0_int32, &
        "random surface-shape property reports no failure seed")
    call check_summary("random surface-shape objective properties")
    if (.not. all_passed) error stop 1

contains

    logical function surface_shape_case(case_seed)
        integer(int32), intent(in) :: case_seed
        type(property_rng_t) :: rng
        type(fortsparse_status_t) :: status
        integer :: dimension, point_count, component, point
        real(dp), allocatable :: candidate(:, :), target(:, :), candidate_dot(:, :)
        real(dp), allocatable :: target_dot(:, :), weights(:), weights_dot(:)
        real(dp), allocatable :: candidate_bar(:, :), target_bar(:, :), weights_bar(:)
        real(dp), allocatable :: bad_target(:, :)
        real(dp) :: objective, objective_dot, objective_plus, objective_minus
        real(dp) :: objective_bar, expected, expected_dot, finite_difference_dot
        real(dp) :: difference, difference_dot, adjoint_left, adjoint_right
        real(dp), parameter :: step = 1.0e-7_dp

        surface_shape_case = .false.
        call property_rng_initialize(rng, case_seed)
        dimension = property_random_integer(rng, 1, 3)
        point_count = property_random_integer(rng, 2, 6)
        allocate(candidate(dimension, point_count), target(dimension, point_count))
        allocate(candidate_dot(dimension, point_count))
        allocate(target_dot(dimension, point_count))
        allocate(weights(point_count), weights_dot(point_count))
        allocate(candidate_bar(dimension, point_count))
        allocate(target_bar(dimension, point_count), weights_bar(point_count))
        allocate(bad_target(dimension, point_count - 1))

        do point = 1, point_count
            weights(point) = 0.2_dp + 1.8_dp*property_random_unit(rng)
            weights_dot(point) = 0.6_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)
            do component = 1, dimension
                candidate(component, point) = &
                    2.0_dp*property_random_unit(rng) - 1.0_dp
                target(component, point) = &
                    2.0_dp*property_random_unit(rng) - 1.0_dp
                candidate_dot(component, point) = &
                    0.6_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)
                target_dot(component, point) = &
                    0.6_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)
            end do
        end do

        call evaluate_surface_shape_objective( &
            candidate, target, weights, objective, status)
        if (status%code /= 0) return
        expected = 0.0_dp
        do point = 1, point_count
            do component = 1, dimension
                difference = candidate(component, point) - target(component, point)
                expected = expected + 0.5_dp*weights(point)*difference**2
            end do
        end do
        if (abs(objective - expected) > 5.0e-13_dp) return

        call evaluate_surface_shape_objective_jvp( &
            candidate, target, weights, candidate_dot, target_dot, weights_dot, &
            objective_dot, status)
        if (status%code /= 0) return
        expected_dot = 0.0_dp
        do point = 1, point_count
            do component = 1, dimension
                difference = candidate(component, point) - target(component, point)
                difference_dot = candidate_dot(component, point) - &
                    target_dot(component, point)
                expected_dot = expected_dot + weights(point)*difference* &
                    difference_dot + 0.5_dp*weights_dot(point)*difference**2
            end do
        end do
        if (abs(objective_dot - expected_dot) > 5.0e-12_dp) return

        call evaluate_surface_shape_objective( &
            candidate + step*candidate_dot, target + step*target_dot, &
            weights + step*weights_dot, objective_plus, status)
        if (status%code /= 0) return
        call evaluate_surface_shape_objective( &
            candidate - step*candidate_dot, target - step*target_dot, &
            weights - step*weights_dot, objective_minus, status)
        if (status%code /= 0) return
        finite_difference_dot = (objective_plus - objective_minus)/(2.0_dp*step)
        if (abs(objective_dot - finite_difference_dot) > 3.0e-8_dp) return

        objective_bar = 0.2_dp + 1.6_dp*property_random_unit(rng)
        call evaluate_surface_shape_objective_vjp( &
            candidate, target, weights, objective_bar, candidate_bar, target_bar, &
            weights_bar, status)
        if (status%code /= 0) return
        adjoint_left = sum(candidate_bar*candidate_dot) + &
            sum(target_bar*target_dot) + dot_product(weights_bar, weights_dot)
        adjoint_right = objective_bar*objective_dot
        if (abs(adjoint_left - adjoint_right) > &
            5.0e-12_dp*max(1.0_dp, abs(adjoint_left), abs(adjoint_right))) return

        weights(1) = 0.0_dp
        call evaluate_surface_shape_objective( &
            candidate, target, weights, objective, status)
        if (status%code == 0) return
        weights(1) = 0.2_dp
        bad_target = target(:, 1:point_count - 1)
        call evaluate_surface_shape_objective( &
            candidate, bad_target, weights, objective, status)
        if (status%code == 0) return
        surface_shape_case = .true.
    end function surface_shape_case

end program test_surface_shape_objective_properties
