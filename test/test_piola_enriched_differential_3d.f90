program test_piola_enriched_differential_3d
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        evaluate_piola_enriched_vector_differential_3d, &
        evaluate_piola_enriched_vector_differential_3d_jvp, &
        evaluate_piola_enriched_vector_differential_3d_vjp, &
        PIOLA_COVARIANT, PIOLA_CONTRAVARIANT
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    implicit none

    real(dp), parameter :: epsilon_fd = 1.0e-7_dp
    logical :: all_passed

    all_passed = .true.
    call check_map(PIOLA_COVARIANT, "3D covariant Piola H(curl)")
    call check_map(PIOLA_CONTRAVARIANT, "3D contravariant Piola H(div)")
    call check_summary("Piola-enriched 3D differential")
    if (.not. all_passed) error stop 1

contains

    subroutine check_map(map_kind, label)
        integer, intent(in) :: map_kind
        character(len=*), intent(in) :: label
        integer, parameter :: points = 2, basis = 2
        real(dp) :: jacobian(3, 3, points), jacobian_dot(3, 3, points)
        real(dp) :: reference(3, points, basis), reference_dot(3, points, basis)
        real(dp) :: reference_curl(3, points, basis), reference_curl_dot(3, points, basis)
        real(dp) :: reference_divergence(points, basis), reference_divergence_dot(points, basis)
        real(dp) :: activation(points, basis), activation_dot(points, basis)
        real(dp) :: activation_gradient(3, points, basis), activation_gradient_dot(3, points, basis)
        real(dp) :: values(3, points, basis), curl(3, points, basis), divergence(points, basis)
        real(dp) :: values_dot(3, points, basis), curl_dot(3, points, basis)
        real(dp) :: divergence_dot(points, basis)
        real(dp) :: values_plus(3, points, basis), curl_plus(3, points, basis)
        real(dp) :: divergence_plus(points, basis)
        real(dp) :: values_minus(3, points, basis), curl_minus(3, points, basis)
        real(dp) :: divergence_minus(points, basis)
        real(dp) :: expected_values(3, points, basis), expected_curl(3, points, basis)
        real(dp) :: expected_divergence(points, basis)
        real(dp) :: values_bar(3, points, basis), curl_bar(3, points, basis)
        real(dp) :: divergence_bar(points, basis)
        real(dp) :: jacobian_bar(3, 3, points), reference_bar(3, points, basis)
        real(dp) :: reference_curl_bar(3, points, basis)
        real(dp) :: reference_divergence_bar(points, basis)
        real(dp) :: activation_bar(points, basis), activation_gradient_bar(3, points, basis)
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
        reference_curl(:, :, 1) = reshape([0.7_dp, -0.4_dp, 0.9_dp, 0.2_dp, 0.6_dp, -0.8_dp], [3, 2])
        reference_curl(:, :, 2) = reshape([-0.5_dp, 0.3_dp, 0.1_dp, 0.8_dp, -0.2_dp, 0.4_dp], [3, 2])
        reference_curl_dot(:, :, 1) = reshape([0.03_dp, 0.04_dp, -0.02_dp, 0.01_dp, -0.05_dp, 0.02_dp], [3, 2])
        reference_curl_dot(:, :, 2) = reshape([-0.01_dp, 0.02_dp, 0.04_dp, -0.03_dp, 0.05_dp, 0.01_dp], [3, 2])
        reference_divergence = reshape([0.6_dp, -0.3_dp, 0.8_dp, 0.2_dp], [points, basis])
        reference_divergence_dot = reshape([0.02_dp, -0.01_dp, 0.03_dp, -0.04_dp], [points, basis])
        activation = reshape([0.9_dp, -0.6_dp, 1.2_dp, 0.5_dp], [points, basis])
        activation_dot = reshape([0.03_dp, 0.02_dp, -0.01_dp, 0.04_dp], [points, basis])
        activation_gradient(:, :, 1) = reshape([0.2_dp, -0.1_dp, 0.3_dp, &
            0.4_dp, -0.2_dp, 0.5_dp], [3, 2])
        activation_gradient(:, :, 2) = reshape([-0.3_dp, 0.6_dp, 0.1_dp, &
            0.7_dp, -0.4_dp, 0.2_dp], [3, 2])
        activation_gradient_dot(:, :, 1) = reshape([0.01_dp, 0.02_dp, -0.03_dp, &
            -0.02_dp, 0.04_dp, 0.01_dp], [3, 2])
        activation_gradient_dot(:, :, 2) = reshape([-0.04_dp, 0.03_dp, 0.02_dp, &
            0.01_dp, -0.02_dp, 0.05_dp], [3, 2])
        values_bar = reshape([0.2_dp, -0.1_dp, 0.3_dp, 0.4_dp, -0.2_dp, 0.5_dp, &
            -0.3_dp, 0.6_dp, 0.1_dp, 0.7_dp, -0.4_dp, 0.2_dp], [3, points, basis])
        curl_bar = reshape([-0.1_dp, 0.2_dp, -0.3_dp, 0.4_dp, 0.5_dp, -0.6_dp, &
            0.7_dp, -0.8_dp, 0.9_dp, -1.0_dp, 1.1_dp, -1.2_dp], [3, points, basis])
        divergence_bar = reshape([0.3_dp, -0.4_dp, 0.5_dp, -0.6_dp], [points, basis])

        call expected(map_kind, jacobian, reference, reference_curl, reference_divergence, &
            activation, activation_gradient, expected_values, expected_curl, expected_divergence)
        call evaluate_piola_enriched_vector_differential_3d( &
            map_kind, jacobian, reference, reference_curl, reference_divergence, activation, &
            activation_gradient, values, curl, divergence, status)
        call record(status%code == FORTSPARSE_OK .and. &
            maxval(abs(values - expected_values)) < 1.0e-14_dp .and. &
            maxval(abs(curl - expected_curl)) < 1.0e-14_dp .and. &
            maxval(abs(divergence - expected_divergence)) < 1.0e-14_dp, label//" value oracle")

        call evaluate_piola_enriched_vector_differential_3d_jvp( &
            map_kind, jacobian, reference, reference_curl, reference_divergence, activation, &
            activation_gradient, jacobian_dot, reference_dot, reference_curl_dot, &
            reference_divergence_dot, activation_dot, activation_gradient_dot, values_dot, &
            curl_dot, divergence_dot, status)
        call evaluate_piola_enriched_vector_differential_3d( &
            map_kind, jacobian + epsilon_fd*jacobian_dot, reference + epsilon_fd*reference_dot, &
            reference_curl + epsilon_fd*reference_curl_dot, &
            reference_divergence + epsilon_fd*reference_divergence_dot, &
            activation + epsilon_fd*activation_dot, &
            activation_gradient + epsilon_fd*activation_gradient_dot, values_plus, curl_plus, &
            divergence_plus, status)
        call evaluate_piola_enriched_vector_differential_3d( &
            map_kind, jacobian - epsilon_fd*jacobian_dot, reference - epsilon_fd*reference_dot, &
            reference_curl - epsilon_fd*reference_curl_dot, &
            reference_divergence - epsilon_fd*reference_divergence_dot, &
            activation - epsilon_fd*activation_dot, &
            activation_gradient - epsilon_fd*activation_gradient_dot, values_minus, curl_minus, &
            divergence_minus, status)
        call record(status%code == FORTSPARSE_OK .and. &
            maxval(abs(values_dot - (values_plus - values_minus)/(2.0_dp*epsilon_fd))) < 3.0e-8_dp .and. &
            maxval(abs(curl_dot - (curl_plus - curl_minus)/(2.0_dp*epsilon_fd))) < 3.0e-8_dp .and. &
            maxval(abs(divergence_dot - (divergence_plus - divergence_minus)/(2.0_dp*epsilon_fd))) < 3.0e-8_dp, &
            label//" JVP central difference")

        call evaluate_piola_enriched_vector_differential_3d_vjp( &
            map_kind, jacobian, reference, reference_curl, reference_divergence, activation, &
            activation_gradient, values_bar, curl_bar, divergence_bar, jacobian_bar, &
            reference_bar, reference_curl_bar, reference_divergence_bar, activation_bar, &
            activation_gradient_bar, status)
        lhs = sum(values_bar*values_dot) + sum(curl_bar*curl_dot) + sum(divergence_bar*divergence_dot)
        rhs = sum(jacobian_bar*jacobian_dot) + sum(reference_bar*reference_dot) + &
            sum(reference_curl_bar*reference_curl_dot) + sum(reference_divergence_bar*reference_divergence_dot) + &
            sum(activation_bar*activation_dot) + sum(activation_gradient_bar*activation_gradient_dot)
        call record(status%code == FORTSPARSE_OK .and. abs(lhs - rhs) < 3.0e-12_dp, &
            label//" VJP dot-product oracle")
    end subroutine check_map

    subroutine expected(map_kind, jacobian, reference, reference_curl, reference_divergence, &
            activation, activation_gradient, values, curl, divergence)
        integer, intent(in) :: map_kind
        real(dp), intent(in) :: jacobian(:, :, :), reference(:, :, :), reference_curl(:, :, :)
        real(dp), intent(in) :: reference_divergence(:, :), activation(:, :), activation_gradient(:, :, :)
        real(dp), intent(out) :: values(:, :, :), curl(:, :, :), divergence(:, :)
        real(dp) :: inverse(3, 3), determinant, mapped(3), mapped_curl(3), mapped_divergence
        integer :: point, basis

        do point = 1, size(reference, 2)
            call inverse_3x3(jacobian(:, :, point), inverse, determinant)
            do basis = 1, size(reference, 3)
                if (map_kind == PIOLA_COVARIANT) then
                    mapped = matmul(transpose(inverse), reference(:, point, basis))
                    mapped_curl = matmul(jacobian(:, :, point), reference_curl(:, point, basis))/determinant
                    values(:, point, basis) = activation(point, basis)*mapped
                    curl(:, point, basis) = activation(point, basis)*mapped_curl + &
                        cross3(activation_gradient(:, point, basis), mapped)
                    divergence(point, basis) = 0.0_dp
                else
                    mapped = matmul(jacobian(:, :, point), reference(:, point, basis))/determinant
                    mapped_divergence = reference_divergence(point, basis)/determinant
                    values(:, point, basis) = activation(point, basis)*mapped
                    curl(:, point, basis) = 0.0_dp
                    divergence(point, basis) = activation(point, basis)*mapped_divergence + &
                        dot_product(activation_gradient(:, point, basis), mapped)
                end if
            end do
        end do
    end subroutine expected

    subroutine inverse_3x3(matrix, inverse, determinant)
        real(dp), intent(in) :: matrix(3, 3)
        real(dp), intent(out) :: inverse(3, 3), determinant

        determinant = matrix(1, 1)*(matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2)) - &
            matrix(1, 2)*(matrix(2, 1)*matrix(3, 3) - matrix(2, 3)*matrix(3, 1)) + &
            matrix(1, 3)*(matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1))
        inverse(1, 1) = (matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2))/determinant
        inverse(1, 2) = (matrix(1, 3)*matrix(3, 2) - matrix(1, 2)*matrix(3, 3))/determinant
        inverse(1, 3) = (matrix(1, 2)*matrix(2, 3) - matrix(1, 3)*matrix(2, 2))/determinant
        inverse(2, 1) = (matrix(2, 3)*matrix(3, 1) - matrix(2, 1)*matrix(3, 3))/determinant
        inverse(2, 2) = (matrix(1, 1)*matrix(3, 3) - matrix(1, 3)*matrix(3, 1))/determinant
        inverse(2, 3) = (matrix(1, 3)*matrix(2, 1) - matrix(1, 1)*matrix(2, 3))/determinant
        inverse(3, 1) = (matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1))/determinant
        inverse(3, 2) = (matrix(1, 2)*matrix(3, 1) - matrix(1, 1)*matrix(3, 2))/determinant
        inverse(3, 3) = (matrix(1, 1)*matrix(2, 2) - matrix(1, 2)*matrix(2, 1))/determinant
    end subroutine inverse_3x3

    pure function cross3(left, right) result(product)
        real(dp), intent(in) :: left(3), right(3)
        real(dp) :: product(3)

        product = [left(2)*right(3) - left(3)*right(2), &
            left(3)*right(1) - left(1)*right(3), &
            left(1)*right(2) - left(2)*right(1)]
    end function cross3

    subroutine record(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        call check_condition(condition, message)
        all_passed = all_passed .and. condition
    end subroutine record

end program test_piola_enriched_differential_3d
