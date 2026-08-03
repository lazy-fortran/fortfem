program test_piola_enriched_vector
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        evaluate_piola_enriched_vector_values, &
        evaluate_piola_enriched_vector_values_jvp, &
        evaluate_piola_enriched_vector_values_vjp, &
        PIOLA_COVARIANT, PIOLA_CONTRAVARIANT
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, FORTSPARSE_OK
    implicit none

    real(dp), parameter :: epsilon_fd = 1.0e-7_dp
    logical :: all_passed

    all_passed = .true.
    call check_2d(PIOLA_COVARIANT, "2D covariant Piola enrichment")
    call check_2d(PIOLA_CONTRAVARIANT, "2D contravariant Piola enrichment")
    call check_3d(PIOLA_COVARIANT, "3D covariant Piola enrichment")
    call check_3d(PIOLA_CONTRAVARIANT, "3D contravariant Piola enrichment")
    call check_summary("Piola-enriched vector values")
    if (.not. all_passed) error stop 1

contains

    subroutine check_2d(map_kind, label)
        integer, intent(in) :: map_kind
        character(len=*), intent(in) :: label
        real(dp) :: jacobian(2, 2, 2), jacobian_dot(2, 2, 2)
        real(dp) :: reference(2, 2, 2), reference_dot(2, 2, 2)
        real(dp) :: activation(2, 2), activation_dot(2, 2)
        real(dp) :: physical(2, 2, 2), expected(2, 2, 2)
        real(dp) :: physical_dot(2, 2, 2)
        real(dp) :: physical_plus(2, 2, 2), physical_minus(2, 2, 2)
        real(dp) :: jacobian_plus(2, 2, 2), jacobian_minus(2, 2, 2)
        real(dp) :: reference_plus(2, 2, 2), reference_minus(2, 2, 2)
        real(dp) :: activation_plus(2, 2), activation_minus(2, 2)
        real(dp) :: physical_bar(2, 2, 2), jacobian_bar(2, 2, 2)
        real(dp) :: reference_bar(2, 2, 2), activation_bar(2, 2)
        real(dp) :: jacobian_bad(2, 2, 2)
        real(dp) :: lhs, rhs
        type(fortsparse_status_t) :: status

        jacobian(:, :, 1) = reshape([1.7_dp, 0.2_dp, -0.1_dp, 1.3_dp], [2, 2])
        jacobian(:, :, 2) = reshape([1.2_dp, -0.3_dp, 0.15_dp, 1.8_dp], [2, 2])
        jacobian_dot(:, :, 1) = reshape([0.03_dp, -0.02_dp, 0.01_dp, 0.04_dp], [2, 2])
        jacobian_dot(:, :, 2) = reshape([-0.02_dp, 0.01_dp, 0.03_dp, -0.01_dp], [2, 2])
        reference(:, :, 1) = reshape([1.0_dp, -0.5_dp, 0.4_dp, 1.3_dp], [2, 2])
        reference(:, :, 2) = reshape([-0.7_dp, 0.8_dp, 1.1_dp, -0.2_dp], [2, 2])
        reference_dot(:, :, 1) = reshape([0.1_dp, -0.04_dp, 0.02_dp, 0.06_dp], [2, 2])
        reference_dot(:, :, 2) = reshape([-0.03_dp, 0.05_dp, 0.07_dp, -0.02_dp], [2, 2])
        activation = reshape([1.1_dp, -0.4_dp, 0.7_dp, 1.3_dp], [2, 2])
        activation_dot = reshape([0.02_dp, -0.03_dp, 0.04_dp, 0.01_dp], [2, 2])
        physical_bar = reshape([0.2_dp, -0.1_dp, 0.4_dp, 0.3_dp, &
            -0.5_dp, 0.6_dp, -0.2_dp, 0.7_dp], [2, 2, 2])

        call expected_2d(map_kind, jacobian, reference, activation, expected)
        call evaluate_piola_enriched_vector_values( &
            map_kind, jacobian, reference, activation, physical, status)
        call record(status%code == FORTSPARSE_OK .and. &
            maxval(abs(physical - expected)) < 1.0e-14_dp, label//" value oracle")

        call evaluate_piola_enriched_vector_values_jvp( &
            map_kind, jacobian, reference, activation, jacobian_dot, &
            reference_dot, activation_dot, physical_dot, status)
        jacobian_plus = jacobian + epsilon_fd*jacobian_dot
        jacobian_minus = jacobian - epsilon_fd*jacobian_dot
        reference_plus = reference + epsilon_fd*reference_dot
        reference_minus = reference - epsilon_fd*reference_dot
        activation_plus = activation + epsilon_fd*activation_dot
        activation_minus = activation - epsilon_fd*activation_dot
        call evaluate_piola_enriched_vector_values( &
            map_kind, jacobian_plus, reference_plus, activation_plus, &
            physical_plus, status)
        call evaluate_piola_enriched_vector_values( &
            map_kind, jacobian_minus, reference_minus, activation_minus, &
            physical_minus, status)
        call record(status%code == FORTSPARSE_OK .and. &
            maxval(abs(physical_dot - (physical_plus - physical_minus)/ &
            (2.0_dp*epsilon_fd))) < 2.0e-8_dp, label//" JVP oracle")

        call evaluate_piola_enriched_vector_values_vjp( &
            map_kind, jacobian, reference, activation, physical_bar, &
            jacobian_bar, reference_bar, activation_bar, status)
        lhs = sum(physical_bar*physical_dot)
        rhs = sum(jacobian_bar*jacobian_dot) + sum(reference_bar*reference_dot) + &
            sum(activation_bar*activation_dot)
        call record(status%code == FORTSPARSE_OK .and. abs(lhs - rhs) < 1.0e-13_dp, &
            label//" VJP dot-product oracle")

        jacobian_bad = jacobian
        jacobian_bad(:, :, 1) = reshape([1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp], [2, 2])
        call evaluate_piola_enriched_vector_values( &
            map_kind, jacobian_bad, reference, activation, physical, status)
        call record(status%code /= FORTSPARSE_OK, label//" rejects singular Jacobian")
        call evaluate_piola_enriched_vector_values( &
            99, jacobian, reference, activation, physical, status)
        call record(status%code /= FORTSPARSE_OK, label//" rejects invalid map kind")
    end subroutine check_2d

    subroutine check_3d(map_kind, label)
        integer, intent(in) :: map_kind
        character(len=*), intent(in) :: label
        real(dp) :: jacobian(3, 3, 2), jacobian_dot(3, 3, 2)
        real(dp) :: reference(3, 2, 2), reference_dot(3, 2, 2)
        real(dp) :: activation(2, 2), activation_dot(2, 2)
        real(dp) :: physical(3, 2, 2), expected(3, 2, 2)
        real(dp) :: physical_dot(3, 2, 2)
        real(dp) :: physical_plus(3, 2, 2), physical_minus(3, 2, 2)
        real(dp) :: jacobian_plus(3, 3, 2), jacobian_minus(3, 3, 2)
        real(dp) :: reference_plus(3, 2, 2), reference_minus(3, 2, 2)
        real(dp) :: activation_plus(2, 2), activation_minus(2, 2)
        real(dp) :: physical_bar(3, 2, 2), jacobian_bar(3, 3, 2)
        real(dp) :: reference_bar(3, 2, 2), activation_bar(2, 2)
        real(dp) :: jacobian_bad(3, 3, 2)
        real(dp) :: lhs, rhs
        type(fortsparse_status_t) :: status

        jacobian(:, :, 1) = reshape([1.4_dp, 0.1_dp, 0.2_dp, &
            0.0_dp, 1.6_dp, 0.15_dp, 0.1_dp, -0.05_dp, 1.3_dp], [3, 3])
        jacobian(:, :, 2) = reshape([1.2_dp, -0.1_dp, 0.05_dp, &
            0.2_dp, 1.5_dp, 0.1_dp, -0.08_dp, 0.04_dp, 1.7_dp], [3, 3])
        jacobian_dot(:, :, 1) = reshape([0.02_dp, -0.01_dp, 0.03_dp, &
            0.01_dp, 0.04_dp, -0.02_dp, -0.03_dp, 0.02_dp, 0.01_dp], [3, 3])
        jacobian_dot(:, :, 2) = reshape([-0.01_dp, 0.03_dp, 0.02_dp, &
            0.02_dp, -0.02_dp, 0.01_dp, 0.04_dp, -0.01_dp, 0.03_dp], [3, 3])
        reference(:, :, 1) = reshape([1.0_dp, -0.5_dp, 0.8_dp, &
            0.4_dp, 1.3_dp, -0.7_dp], [3, 2])
        reference(:, :, 2) = reshape([-0.3_dp, 0.9_dp, 0.2_dp, &
            1.1_dp, -0.6_dp, 0.5_dp], [3, 2])
        reference_dot(:, :, 1) = reshape([0.04_dp, -0.03_dp, 0.01_dp, &
            0.02_dp, 0.05_dp, -0.02_dp], [3, 2])
        reference_dot(:, :, 2) = reshape([-0.02_dp, 0.01_dp, 0.03_dp, &
            0.04_dp, -0.01_dp, 0.02_dp], [3, 2])
        activation = reshape([0.9_dp, -0.6_dp, 1.2_dp, 0.5_dp], [2, 2])
        activation_dot = reshape([0.03_dp, 0.02_dp, -0.01_dp, 0.04_dp], [2, 2])
        physical_bar = reshape([0.2_dp, -0.1_dp, 0.3_dp, 0.4_dp, -0.2_dp, 0.5_dp, &
            -0.3_dp, 0.6_dp, 0.1_dp, 0.7_dp, -0.4_dp, 0.2_dp], [3, 2, 2])

        call expected_3d(map_kind, jacobian, reference, activation, expected)
        call evaluate_piola_enriched_vector_values( &
            map_kind, jacobian, reference, activation, physical, status)
        call record(status%code == FORTSPARSE_OK .and. &
            maxval(abs(physical - expected)) < 1.0e-14_dp, label//" value oracle")

        call evaluate_piola_enriched_vector_values_jvp( &
            map_kind, jacobian, reference, activation, jacobian_dot, &
            reference_dot, activation_dot, physical_dot, status)
        jacobian_plus = jacobian + epsilon_fd*jacobian_dot
        jacobian_minus = jacobian - epsilon_fd*jacobian_dot
        reference_plus = reference + epsilon_fd*reference_dot
        reference_minus = reference - epsilon_fd*reference_dot
        activation_plus = activation + epsilon_fd*activation_dot
        activation_minus = activation - epsilon_fd*activation_dot
        call evaluate_piola_enriched_vector_values( &
            map_kind, jacobian_plus, reference_plus, activation_plus, &
            physical_plus, status)
        call evaluate_piola_enriched_vector_values( &
            map_kind, jacobian_minus, reference_minus, activation_minus, &
            physical_minus, status)
        call record(status%code == FORTSPARSE_OK .and. &
            maxval(abs(physical_dot - (physical_plus - physical_minus)/ &
            (2.0_dp*epsilon_fd))) < 2.0e-8_dp, label//" JVP oracle")

        call evaluate_piola_enriched_vector_values_vjp( &
            map_kind, jacobian, reference, activation, physical_bar, &
            jacobian_bar, reference_bar, activation_bar, status)
        lhs = sum(physical_bar*physical_dot)
        rhs = sum(jacobian_bar*jacobian_dot) + sum(reference_bar*reference_dot) + &
            sum(activation_bar*activation_dot)
        call record(status%code == FORTSPARSE_OK .and. abs(lhs - rhs) < 1.0e-13_dp, &
            label//" VJP dot-product oracle")

        jacobian_bad = jacobian
        jacobian_bad(:, :, 1) = 0.0_dp
        call evaluate_piola_enriched_vector_values( &
            map_kind, jacobian_bad, reference, activation, physical, status)
        call record(status%code /= FORTSPARSE_OK, label//" rejects singular Jacobian")
    end subroutine check_3d

    subroutine expected_2d(map_kind, jacobian, reference, activation, physical)
        integer, intent(in) :: map_kind
        real(dp), intent(in) :: jacobian(:, :, :), reference(:, :, :)
        real(dp), intent(in) :: activation(:, :)
        real(dp), intent(out) :: physical(:, :, :)
        real(dp) :: determinant, inverse(2, 2), mapped(2)
        integer :: point, basis

        do point = 1, size(reference, 2)
            determinant = jacobian(1, 1, point)*jacobian(2, 2, point) - &
                jacobian(1, 2, point)*jacobian(2, 1, point)
            inverse(1, 1) = jacobian(2, 2, point)/determinant
            inverse(1, 2) = -jacobian(1, 2, point)/determinant
            inverse(2, 1) = -jacobian(2, 1, point)/determinant
            inverse(2, 2) = jacobian(1, 1, point)/determinant
            do basis = 1, size(reference, 3)
                if (map_kind == PIOLA_COVARIANT) then
                    mapped = matmul(transpose(inverse), reference(:, point, basis))
                else
                    mapped = matmul(jacobian(:, :, point), &
                        reference(:, point, basis))/determinant
                end if
                physical(:, point, basis) = activation(point, basis)*mapped
            end do
        end do
    end subroutine expected_2d

    subroutine expected_3d(map_kind, jacobian, reference, activation, physical)
        integer, intent(in) :: map_kind
        real(dp), intent(in) :: jacobian(:, :, :), reference(:, :, :)
        real(dp), intent(in) :: activation(:, :)
        real(dp), intent(out) :: physical(:, :, :)
        real(dp) :: determinant, inverse(3, 3), mapped(3)
        integer :: point, basis

        do point = 1, size(reference, 2)
            call inverse_3x3(jacobian(:, :, point), inverse, determinant)
            do basis = 1, size(reference, 3)
                if (map_kind == PIOLA_COVARIANT) then
                    mapped = matmul(transpose(inverse), reference(:, point, basis))
                else
                    mapped = matmul(jacobian(:, :, point), &
                        reference(:, point, basis))/determinant
                end if
                physical(:, point, basis) = activation(point, basis)*mapped
            end do
        end do
    end subroutine expected_3d

    subroutine inverse_3x3(matrix, inverse, determinant)
        real(dp), intent(in) :: matrix(3, 3)
        real(dp), intent(out) :: inverse(3, 3), determinant

        determinant = matrix(1, 1)*(matrix(2, 2)*matrix(3, 3) - &
            matrix(2, 3)*matrix(3, 2)) - matrix(1, 2)*( &
            matrix(2, 1)*matrix(3, 3) - matrix(2, 3)*matrix(3, 1)) + &
            matrix(1, 3)*(matrix(2, 1)*matrix(3, 2) - &
            matrix(2, 2)*matrix(3, 1))
        inverse(1, 1) = (matrix(2, 2)*matrix(3, 3) - &
            matrix(2, 3)*matrix(3, 2))/determinant
        inverse(1, 2) = (matrix(1, 3)*matrix(3, 2) - &
            matrix(1, 2)*matrix(3, 3))/determinant
        inverse(1, 3) = (matrix(1, 2)*matrix(2, 3) - &
            matrix(1, 3)*matrix(2, 2))/determinant
        inverse(2, 1) = (matrix(2, 3)*matrix(3, 1) - &
            matrix(2, 1)*matrix(3, 3))/determinant
        inverse(2, 2) = (matrix(1, 1)*matrix(3, 3) - &
            matrix(1, 3)*matrix(3, 1))/determinant
        inverse(2, 3) = (matrix(1, 3)*matrix(2, 1) - &
            matrix(1, 1)*matrix(2, 3))/determinant
        inverse(3, 1) = (matrix(2, 1)*matrix(3, 2) - &
            matrix(2, 2)*matrix(3, 1))/determinant
        inverse(3, 2) = (matrix(1, 2)*matrix(3, 1) - &
            matrix(1, 1)*matrix(3, 2))/determinant
        inverse(3, 3) = (matrix(1, 1)*matrix(2, 2) - &
            matrix(1, 2)*matrix(2, 1))/determinant
    end subroutine inverse_3x3

    subroutine record(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        call check_condition(condition, message)
        all_passed = all_passed .and. condition
    end subroutine record

end program test_piola_enriched_vector
