program test_field_aligned_constitutive_properties
    use, intrinsic :: iso_fortran_env, only: int32
    use check, only: check_condition, check_property, check_summary, &
        property_random_unit, property_rng_initialize, property_rng_t
    use fortfem_feec, only: &
        evaluate_field_aligned_constitutive_tensor, evaluate_tensor_power_split
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    logical :: all_passed
    integer(int32) :: first_failed_seed, shrunk_seed

    call check_property( &
        "seeded field-aligned projector/Hall and power identities", &
        20260802_int32, 16, constitutive_case, all_passed, first_failed_seed, &
        shrunk_seed)
    call check_condition(all_passed .and. first_failed_seed == 0_int32 .and. &
        shrunk_seed == 0_int32, &
        "random constitutive property reports no failure seed")
    call check_summary("random field-aligned constitutive properties")
    if (.not. all_passed) error stop 1

contains

    logical function constitutive_case(case_seed)
        integer(int32), intent(in) :: case_seed
        integer, parameter :: trial_count = 2
        type(property_rng_t) :: rng
        type(fortsparse_status_t) :: status
        real(dp) :: direction(3), vector(3), tensor(3, 3), expected(3, 3)
        real(dp) :: symmetric_tensor(3, 3), skew_tensor(3, 3)
        real(dp) :: parallel, perpendicular, hall
        real(dp) :: symmetric_power, skew_power, total_power
        real(dp) :: expected_symmetric_power, expected_skew_power
        real(dp) :: expected_total_power, norm_squared
        real(dp) :: bad_direction(3), near_direction(3)
        integer :: trial
        logical :: direction_found

        constitutive_case = .false.
        call property_rng_initialize(rng, case_seed)
        do trial = 1, trial_count
            call random_unit_direction(rng, direction, direction_found)
            if (.not. direction_found) return
            vector = [ &
                2.0_dp*property_random_unit(rng) - 1.0_dp, &
                2.0_dp*property_random_unit(rng) - 1.0_dp, &
                2.0_dp*property_random_unit(rng) - 1.0_dp]
            parallel = 0.25_dp + 4.75_dp*property_random_unit(rng)
            perpendicular = 0.10_dp + 2.90_dp*property_random_unit(rng)
            hall = -2.0_dp + 4.0_dp*property_random_unit(rng)

            call evaluate_field_aligned_constitutive_tensor( &
                parallel, perpendicular, direction, tensor, status, hall)
            if (status%code /= 0) return
            call independent_projector_hall( &
                parallel, perpendicular, hall, direction, expected)
            if (maxval(abs(tensor - expected)) > 2.0e-13_dp) return

            symmetric_tensor = 0.5_dp*(tensor + transpose(tensor))
            skew_tensor = 0.5_dp*(tensor - transpose(tensor))
            expected_symmetric_power = dot_product(vector, matmul( &
                symmetric_tensor, vector))
            expected_skew_power = dot_product(vector, matmul(skew_tensor, vector))
            expected_total_power = dot_product(vector, matmul(tensor, vector))
            call evaluate_tensor_power_split( &
                tensor, vector, symmetric_power, skew_power, total_power, status)
            if (status%code /= 0) return
            if (abs(symmetric_power - expected_symmetric_power) > 2.0e-13_dp .or. &
                abs(skew_power - expected_skew_power) > 2.0e-13_dp .or. &
                abs(total_power - expected_total_power) > 2.0e-13_dp .or. &
                abs(skew_power) > 2.0e-13_dp .or. &
                abs(total_power - symmetric_power) > 2.0e-13_dp) return

            near_direction = direction*(1.0_dp + 1.0e-13_dp)
            call evaluate_field_aligned_constitutive_tensor( &
                parallel, perpendicular, near_direction, tensor, status, hall)
            if (status%code /= 0) return
            bad_direction = direction*(1.0_dp + 1.0e-6_dp)
            call evaluate_field_aligned_constitutive_tensor( &
                parallel, perpendicular, bad_direction, tensor, status, hall)
            if (status%code == 0) return

            norm_squared = dot_product(direction, direction)
            if (abs(norm_squared - 1.0_dp) > 2.0e-13_dp) return
        end do
        constitutive_case = .true.
    end function constitutive_case

    subroutine random_unit_direction(rng, direction, found)
        type(property_rng_t), intent(inout) :: rng
        real(dp), intent(out) :: direction(3)
        logical, intent(out) :: found
        real(dp) :: norm_squared
        integer :: attempt

        direction = 0.0_dp
        found = .false.
        do attempt = 1, 12
            direction = [ &
                2.0_dp*property_random_unit(rng) - 1.0_dp, &
                2.0_dp*property_random_unit(rng) - 1.0_dp, &
                2.0_dp*property_random_unit(rng) - 1.0_dp]
            norm_squared = dot_product(direction, direction)
            if (norm_squared > 0.05_dp) then
                direction = direction/sqrt(norm_squared)
                found = .true.
                return
            end if
        end do
    end subroutine random_unit_direction

    subroutine independent_projector_hall( &
            parallel, perpendicular, hall, direction, tensor)
        real(dp), intent(in) :: parallel, perpendicular, hall, direction(3)
        real(dp), intent(out) :: tensor(3, 3)
        real(dp) :: cross_matrix(3, 3), projector(3, 3)
        integer :: row, column

        projector = 0.0_dp
        cross_matrix = 0.0_dp
        do row = 1, 3
            projector(row, row) = perpendicular
        end do
        do row = 1, 3
            do column = 1, 3
                projector(row, column) = projector(row, column) + &
                    (parallel - perpendicular)*direction(row)*direction(column)
            end do
        end do
        cross_matrix(1, 2) = -direction(3)
        cross_matrix(1, 3) = direction(2)
        cross_matrix(2, 1) = direction(3)
        cross_matrix(2, 3) = -direction(1)
        cross_matrix(3, 1) = -direction(2)
        cross_matrix(3, 2) = direction(1)
        tensor = projector + hall*cross_matrix
    end subroutine independent_projector_hall

end program test_field_aligned_constitutive_properties
