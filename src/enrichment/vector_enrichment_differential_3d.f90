module fortfem_vector_enrichment_differential_3d
    !! Product-rule curl/divergence actions for a 3D vector enrichment.
    !!
    !! For an activation psi and vector field b, this neutral layer exposes
    !!
    !!   curl(psi b) = psi curl(b) + grad(psi) x b,
    !!   div (psi b) = psi div (b) + grad(psi) . b.
    !!
    !! The caller supplies the base value and gradient and may obtain psi
    !! from a shifted Heaviside, kink, singular, or IGA enrichment.  No
    !! Piola map or physical constitutive law is selected here; the explicit
    !! product terms make the de Rham effect of an enrichment inspectable.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_vector_enrichment_differential_3d
    public :: evaluate_vector_enrichment_differential_3d_jvp
    public :: evaluate_vector_enrichment_differential_3d_vjp

contains

    subroutine evaluate_vector_enrichment_differential_3d( &
            base_values, base_gradient, activation, activation_gradient, &
            enriched_values, curl_values, divergence, status)
        real(dp), intent(in) :: base_values(:, :), base_gradient(:, :, :)
        real(dp), intent(in) :: activation(:), activation_gradient(:, :)
        real(dp), intent(out) :: enriched_values(:, :), curl_values(:, :)
        real(dp), intent(out) :: divergence(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: point
        real(dp) :: base_curl(3), base_divergence

        enriched_values = 0.0_dp
        curl_values = 0.0_dp
        divergence = 0.0_dp
        if (.not. validate_inputs( &
            base_values, base_gradient, activation, activation_gradient, &
            enriched_values, curl_values, divergence, status)) return
        do point = 1, size(activation)
            call base_differential( &
                base_gradient(:, :, point), base_curl, base_divergence)
            enriched_values(:, point) = activation(point)*base_values(:, point)
            curl_values(:, point) = activation(point)*base_curl + &
                cross3(activation_gradient(:, point), base_values(:, point))
            divergence(point) = activation(point)*base_divergence + &
                dot_product(activation_gradient(:, point), base_values(:, point))
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_vector_enrichment_differential_3d

    subroutine evaluate_vector_enrichment_differential_3d_jvp( &
            base_values, base_gradient, activation, activation_gradient, &
            base_values_dot, base_gradient_dot, activation_dot, &
            activation_gradient_dot, enriched_values_dot, curl_values_dot, &
            divergence_dot, status)
        real(dp), intent(in) :: base_values(:, :), base_gradient(:, :, :)
        real(dp), intent(in) :: activation(:), activation_gradient(:, :)
        real(dp), intent(in) :: base_values_dot(:, :), base_gradient_dot(:, :, :)
        real(dp), intent(in) :: activation_dot(:), activation_gradient_dot(:, :)
        real(dp), intent(out) :: enriched_values_dot(:, :), curl_values_dot(:, :)
        real(dp), intent(out) :: divergence_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: point
        real(dp) :: base_curl(3), base_curl_dot(3)
        real(dp) :: base_divergence, base_divergence_dot

        enriched_values_dot = 0.0_dp
        curl_values_dot = 0.0_dp
        divergence_dot = 0.0_dp
        if (.not. validate_inputs( &
            base_values, base_gradient, activation, activation_gradient, &
            enriched_values_dot, curl_values_dot, divergence_dot, status)) return
        if (.not. validate_direction( &
            base_values_dot, base_gradient_dot, activation_dot, &
            activation_gradient_dot, size(activation))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "vector enrichment JVP has incompatible increments")
            return
        end if
        do point = 1, size(activation)
            call base_differential( &
                base_gradient(:, :, point), base_curl, base_divergence)
            call base_differential( &
                base_gradient_dot(:, :, point), base_curl_dot, &
                base_divergence_dot)
            enriched_values_dot(:, point) = activation_dot(point)* &
                base_values(:, point) + activation(point)*base_values_dot(:, point)
            curl_values_dot(:, point) = activation_dot(point)*base_curl + &
                activation(point)*base_curl_dot + &
                cross3(activation_gradient_dot(:, point), base_values(:, point)) + &
                cross3(activation_gradient(:, point), base_values_dot(:, point))
            divergence_dot(point) = activation_dot(point)*base_divergence + &
                activation(point)*base_divergence_dot + &
                dot_product(activation_gradient_dot(:, point), &
                base_values(:, point)) + &
                dot_product(activation_gradient(:, point), base_values_dot(:, point))
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_vector_enrichment_differential_3d_jvp

    subroutine evaluate_vector_enrichment_differential_3d_vjp( &
            base_values, base_gradient, activation, activation_gradient, &
            enriched_values_bar, curl_values_bar, divergence_bar, base_values_bar, &
            base_gradient_bar, activation_bar, activation_gradient_bar, status)
        real(dp), intent(in) :: base_values(:, :), base_gradient(:, :, :)
        real(dp), intent(in) :: activation(:), activation_gradient(:, :)
        real(dp), intent(in) :: enriched_values_bar(:, :), curl_values_bar(:, :)
        real(dp), intent(in) :: divergence_bar(:)
        real(dp), intent(out) :: base_values_bar(:, :), base_gradient_bar(:, :, :)
        real(dp), intent(out) :: activation_bar(:), activation_gradient_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: point
        real(dp) :: base_curl(3), base_divergence
        real(dp) :: curl_base_bar(3), divergence_base_bar

        base_values_bar = 0.0_dp
        base_gradient_bar = 0.0_dp
        activation_bar = 0.0_dp
        activation_gradient_bar = 0.0_dp
        if (.not. validate_inputs( &
            base_values, base_gradient, activation, activation_gradient, &
            enriched_values_bar, curl_values_bar, divergence_bar, status)) return
        if (.not. validate_adjoint( &
            base_values_bar, base_gradient_bar, activation_bar, &
            activation_gradient_bar, size(activation))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "vector enrichment VJP has incompatible cotangents")
            return
        end if
        do point = 1, size(activation)
            call base_differential( &
                base_gradient(:, :, point), base_curl, base_divergence)
            base_values_bar(:, point) = base_values_bar(:, point) + &
                activation(point)*enriched_values_bar(:, point)
            activation_bar(point) = activation_bar(point) + &
                dot_product(enriched_values_bar(:, point), base_values(:, point))

            curl_base_bar = activation(point)*curl_values_bar(:, point)
            activation_bar(point) = activation_bar(point) + &
                dot_product(curl_values_bar(:, point), base_curl)
            base_values_bar(:, point) = base_values_bar(:, point) + &
                cross3(curl_values_bar(:, point), activation_gradient(:, point))
            activation_gradient_bar(:, point) = &
                activation_gradient_bar(:, point) + &
                cross3(base_values(:, point), curl_values_bar(:, point))

            divergence_base_bar = activation(point)*divergence_bar(point)
            activation_bar(point) = activation_bar(point) + &
                divergence_bar(point)*base_divergence
            base_values_bar(:, point) = base_values_bar(:, point) + &
                divergence_bar(point)*activation_gradient(:, point)
            activation_gradient_bar(:, point) = &
                activation_gradient_bar(:, point) + &
                divergence_bar(point)*base_values(:, point)
            call base_differential_vjp( &
                curl_base_bar, divergence_base_bar, base_gradient_bar(:, :, point))
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_vector_enrichment_differential_3d_vjp

    subroutine base_differential(gradient, curl_value, divergence_value)
        real(dp), intent(in) :: gradient(3, 3)
        real(dp), intent(out) :: curl_value(3), divergence_value

        curl_value = [gradient(3, 2) - gradient(2, 3), &
            gradient(1, 3) - gradient(3, 1), &
            gradient(2, 1) - gradient(1, 2)]
        divergence_value = gradient(1, 1) + gradient(2, 2) + gradient(3, 3)
    end subroutine base_differential

    subroutine base_differential_vjp(curl_bar, divergence_bar, gradient_bar)
        real(dp), intent(in) :: curl_bar(3), divergence_bar
        real(dp), intent(inout) :: gradient_bar(3, 3)

        gradient_bar(3, 2) = gradient_bar(3, 2) + curl_bar(1)
        gradient_bar(2, 3) = gradient_bar(2, 3) - curl_bar(1)
        gradient_bar(1, 3) = gradient_bar(1, 3) + curl_bar(2)
        gradient_bar(3, 1) = gradient_bar(3, 1) - curl_bar(2)
        gradient_bar(2, 1) = gradient_bar(2, 1) + curl_bar(3)
        gradient_bar(1, 2) = gradient_bar(1, 2) - curl_bar(3)
        gradient_bar(1, 1) = gradient_bar(1, 1) + divergence_bar
        gradient_bar(2, 2) = gradient_bar(2, 2) + divergence_bar
        gradient_bar(3, 3) = gradient_bar(3, 3) + divergence_bar
    end subroutine base_differential_vjp

    logical function validate_inputs( &
            base_values, base_gradient, activation, activation_gradient, &
            enriched_values, curl_values, divergence, status) result(valid)
        real(dp), intent(in) :: base_values(:, :), base_gradient(:, :, :)
        real(dp), intent(in) :: activation(:), activation_gradient(:, :)
        real(dp), intent(in) :: enriched_values(:, :), curl_values(:, :)
        real(dp), intent(in) :: divergence(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: point_count

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "vector enrichment has incompatible arrays")
        point_count = size(base_values, 2)
        if (size(base_values, 1) /= 3 .or. point_count < 1 .or. &
            size(base_gradient, 1) /= 3 .or. size(base_gradient, 2) /= 3 .or. &
            size(base_gradient, 3) /= point_count .or. &
            size(activation) /= point_count .or. &
            size(activation_gradient, 1) /= 3 .or. &
            size(activation_gradient, 2) /= point_count .or. &
            size(enriched_values, 1) /= 3 .or. &
            size(enriched_values, 2) /= point_count .or. &
            size(curl_values, 1) /= 3 .or. size(curl_values, 2) /= point_count .or. &
            size(divergence) /= point_count) return
        if (.not. finite_real_2d(base_values) .or. &
            .not. finite_real_3d(base_gradient) .or. &
            .not. finite_real_1d(activation) .or. &
            .not. finite_real_2d(activation_gradient) .or. &
            .not. finite_real_2d(enriched_values) .or. &
            .not. finite_real_2d(curl_values) .or. &
            .not. finite_real_1d(divergence)) return
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_inputs

    logical function validate_direction( &
            base_values_dot, base_gradient_dot, activation_dot, &
            activation_gradient_dot, point_count) result(valid)
        real(dp), intent(in) :: base_values_dot(:, :), base_gradient_dot(:, :, :)
        real(dp), intent(in) :: activation_dot(:), activation_gradient_dot(:, :)
        integer, intent(in) :: point_count

        valid = size(base_values_dot, 1) == 3 .and. &
            size(base_values_dot, 2) == point_count .and. &
            size(base_gradient_dot, 1) == 3 .and. &
            size(base_gradient_dot, 2) == 3 .and. &
            size(base_gradient_dot, 3) == point_count .and. &
            size(activation_dot) == point_count .and. &
            size(activation_gradient_dot, 1) == 3 .and. &
            size(activation_gradient_dot, 2) == point_count .and. &
            finite_real_2d(base_values_dot) .and. &
            finite_real_3d(base_gradient_dot) .and. &
            finite_real_1d(activation_dot) .and. &
            finite_real_2d(activation_gradient_dot)
    end function validate_direction

    logical function validate_adjoint( &
            base_values_bar, base_gradient_bar, activation_bar, &
            activation_gradient_bar, point_count) result(valid)
        real(dp), intent(in) :: base_values_bar(:, :), base_gradient_bar(:, :, :)
        real(dp), intent(in) :: activation_bar(:), activation_gradient_bar(:, :)
        integer, intent(in) :: point_count

        valid = size(base_values_bar, 1) == 3 .and. &
            size(base_values_bar, 2) == point_count .and. &
            size(base_gradient_bar, 1) == 3 .and. &
            size(base_gradient_bar, 2) == 3 .and. &
            size(base_gradient_bar, 3) == point_count .and. &
            size(activation_bar) == point_count .and. &
            size(activation_gradient_bar, 1) == 3 .and. &
            size(activation_gradient_bar, 2) == point_count
    end function validate_adjoint

    pure function cross3(left, right) result(product)
        real(dp), intent(in) :: left(3), right(3)
        real(dp) :: product(3)

        product = [left(2)*right(3) - left(3)*right(2), &
            left(3)*right(1) - left(1)*right(3), &
            left(1)*right(2) - left(2)*right(1)]
    end function cross3

    pure logical function finite_real_1d(values) result(valid)
        real(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(values))
    end function finite_real_1d

    pure logical function finite_real_2d(values) result(valid)
        real(dp), intent(in) :: values(:, :)

        valid = all(ieee_is_finite(values))
    end function finite_real_2d

    pure logical function finite_real_3d(values) result(valid)
        real(dp), intent(in) :: values(:, :, :)

        valid = all(ieee_is_finite(values))
    end function finite_real_3d

end module fortfem_vector_enrichment_differential_3d
