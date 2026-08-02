module fortfem_batched_vector_enrichment_differential_3d
    !! Batched product-rule differential for vector XFEM/XIGA basis values.
    !!
    !! For each basis function and quadrature point this evaluates
    !!   b_e = psi b,
    !!   curl(b_e) = psi curl(b) + grad(psi) x b,
    !!   div(b_e) = psi div(b) + grad(psi) . b.
    !! The operator is a space-construction primitive: Piola maps, cut
    !! quadrature, continuity, and exact-sequence policy remain caller-owned.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: evaluate_batched_vector_enrichment_differential_3d
    public :: evaluate_batched_vector_enrichment_differential_3d_jvp
    public :: evaluate_batched_vector_enrichment_differential_3d_vjp

    interface finite
        module procedure finite_vector
        module procedure finite_matrix
        module procedure finite_tensor
    end interface finite

contains

    subroutine evaluate_batched_vector_enrichment_differential_3d( &
            base_values, base_curls, base_divergences, activation, activation_gradient, &
            enriched_values, enriched_curls, enriched_divergences, status)
        real(dp), intent(in) :: base_values(:, :, :), base_curls(:, :, :)
        real(dp), intent(in) :: base_divergences(:, :), activation(:, :)
        real(dp), intent(in) :: activation_gradient(:, :, :)
        real(dp), intent(out) :: enriched_values(:, :, :), enriched_curls(:, :, :)
        real(dp), intent(out) :: enriched_divergences(:, :)
        integer, intent(out) :: status
        integer :: basis, point

        enriched_values = 0.0_dp
        enriched_curls = 0.0_dp
        enriched_divergences = 0.0_dp
        status = 1
        if (.not. valid_inputs(base_values, base_curls, base_divergences, activation, &
            activation_gradient, enriched_values, enriched_curls, enriched_divergences)) return
        do point = 1, size(base_values, 3)
            do basis = 1, size(base_values, 2)
                enriched_values(:, basis, point) = activation(basis, point)*base_values(:, basis, point)
                enriched_curls(:, basis, point) = activation(basis, point)*base_curls(:, basis, point) + &
                    cross_product(activation_gradient(:, basis, point), base_values(:, basis, point))
                enriched_divergences(basis, point) = activation(basis, point)*base_divergences(basis, point) + &
                    dot_product(activation_gradient(:, basis, point), base_values(:, basis, point))
            end do
        end do
        status = 0
    end subroutine evaluate_batched_vector_enrichment_differential_3d

    subroutine evaluate_batched_vector_enrichment_differential_3d_jvp( &
            base_values, base_curls, base_divergences, activation, activation_gradient, &
            base_values_dot, base_curls_dot, base_divergences_dot, activation_dot, &
            activation_gradient_dot, enriched_values_dot, enriched_curls_dot, &
            enriched_divergences_dot, status)
        real(dp), intent(in) :: base_values(:, :, :), base_curls(:, :, :)
        real(dp), intent(in) :: base_divergences(:, :), activation(:, :)
        real(dp), intent(in) :: activation_gradient(:, :, :)
        real(dp), intent(in) :: base_values_dot(:, :, :), base_curls_dot(:, :, :)
        real(dp), intent(in) :: base_divergences_dot(:, :), activation_dot(:, :)
        real(dp), intent(in) :: activation_gradient_dot(:, :, :)
        real(dp), intent(out) :: enriched_values_dot(:, :, :), enriched_curls_dot(:, :, :)
        real(dp), intent(out) :: enriched_divergences_dot(:, :)
        integer, intent(out) :: status
        integer :: basis, point

        enriched_values_dot = 0.0_dp
        enriched_curls_dot = 0.0_dp
        enriched_divergences_dot = 0.0_dp
        status = 1
        if (.not. valid_inputs(base_values, base_curls, base_divergences, activation, &
            activation_gradient, enriched_values_dot, enriched_curls_dot, enriched_divergences_dot)) return
        if (.not. same_shapes(base_values_dot, base_values, base_curls_dot, base_curls, &
            base_divergences_dot, base_divergences, activation_dot, activation, &
            activation_gradient_dot, activation_gradient)) return
        do point = 1, size(base_values, 3)
            do basis = 1, size(base_values, 2)
                enriched_values_dot(:, basis, point) = activation_dot(basis, point)*base_values(:, basis, point) + &
                    activation(basis, point)*base_values_dot(:, basis, point)
                enriched_curls_dot(:, basis, point) = activation_dot(basis, point)*base_curls(:, basis, point) + &
                    activation(basis, point)*base_curls_dot(:, basis, point) + &
                    cross_product(activation_gradient_dot(:, basis, point), base_values(:, basis, point)) + &
                    cross_product(activation_gradient(:, basis, point), base_values_dot(:, basis, point))
                enriched_divergences_dot(basis, point) = activation_dot(basis, point)*base_divergences(basis, point) + &
                    activation(basis, point)*base_divergences_dot(basis, point) + &
                    dot_product(activation_gradient_dot(:, basis, point), base_values(:, basis, point)) + &
                    dot_product(activation_gradient(:, basis, point), base_values_dot(:, basis, point))
            end do
        end do
        status = 0
    end subroutine evaluate_batched_vector_enrichment_differential_3d_jvp

    subroutine evaluate_batched_vector_enrichment_differential_3d_vjp( &
            base_values, base_curls, base_divergences, activation, activation_gradient, &
            enriched_values_bar, enriched_curls_bar, enriched_divergences_bar, &
            base_values_bar, base_curls_bar, base_divergences_bar, activation_bar, &
            activation_gradient_bar, status)
        real(dp), intent(in) :: base_values(:, :, :), base_curls(:, :, :)
        real(dp), intent(in) :: base_divergences(:, :), activation(:, :)
        real(dp), intent(in) :: activation_gradient(:, :, :)
        real(dp), intent(in) :: enriched_values_bar(:, :, :), enriched_curls_bar(:, :, :)
        real(dp), intent(in) :: enriched_divergences_bar(:, :)
        real(dp), intent(out) :: base_values_bar(:, :, :), base_curls_bar(:, :, :)
        real(dp), intent(out) :: base_divergences_bar(:, :), activation_bar(:, :)
        real(dp), intent(out) :: activation_gradient_bar(:, :, :)
        integer, intent(out) :: status
        integer :: basis, point
        real(dp) :: scalar_bar

        base_values_bar = 0.0_dp
        base_curls_bar = 0.0_dp
        base_divergences_bar = 0.0_dp
        activation_bar = 0.0_dp
        activation_gradient_bar = 0.0_dp
        status = 1
        if (.not. valid_inputs(base_values, base_curls, base_divergences, activation, &
            activation_gradient, enriched_values_bar, enriched_curls_bar, enriched_divergences_bar)) return
        if (.not. same_shapes(base_values_bar, base_values, base_curls_bar, base_curls, &
            base_divergences_bar, base_divergences, activation_bar, activation, &
            activation_gradient_bar, activation_gradient)) return
        if (any(shape(enriched_values_bar) /= shape(base_values)) .or. &
            any(shape(enriched_curls_bar) /= shape(base_curls)) .or. &
            any(shape(enriched_divergences_bar) /= shape(base_divergences)) .or. &
            .not. finite(enriched_values_bar) .or. .not. finite(enriched_curls_bar) .or. &
            .not. finite(enriched_divergences_bar)) return
        do point = 1, size(base_values, 3)
            do basis = 1, size(base_values, 2)
                scalar_bar = dot_product(enriched_values_bar(:, basis, point), base_values(:, basis, point))
                activation_bar(basis, point) = activation_bar(basis, point) + scalar_bar
                base_values_bar(:, basis, point) = base_values_bar(:, basis, point) + &
                    activation(basis, point)*enriched_values_bar(:, basis, point)
                scalar_bar = dot_product(enriched_curls_bar(:, basis, point), base_curls(:, basis, point))
                activation_bar(basis, point) = activation_bar(basis, point) + scalar_bar
                base_curls_bar(:, basis, point) = base_curls_bar(:, basis, point) + &
                    activation(basis, point)*enriched_curls_bar(:, basis, point)
                base_values_bar(:, basis, point) = base_values_bar(:, basis, point) + &
                    cross_product(enriched_curls_bar(:, basis, point), activation_gradient(:, basis, point))
                activation_gradient_bar(:, basis, point) = activation_gradient_bar(:, basis, point) + &
                    cross_product(base_values(:, basis, point), enriched_curls_bar(:, basis, point))
                scalar_bar = enriched_divergences_bar(basis, point)
                activation_bar(basis, point) = activation_bar(basis, point) + &
                    scalar_bar*base_divergences(basis, point)
                base_divergences_bar(basis, point) = base_divergences_bar(basis, point) + &
                    scalar_bar*activation(basis, point)
                base_values_bar(:, basis, point) = base_values_bar(:, basis, point) + &
                    scalar_bar*activation_gradient(:, basis, point)
                activation_gradient_bar(:, basis, point) = activation_gradient_bar(:, basis, point) + &
                    scalar_bar*base_values(:, basis, point)
            end do
        end do
        status = 0
    end subroutine evaluate_batched_vector_enrichment_differential_3d_vjp

    logical function valid_inputs( &
            base_values, base_curls, base_divergences, activation, activation_gradient, &
            enriched_values, enriched_curls, enriched_divergences) result(valid)
        real(dp), intent(in) :: base_values(:, :, :), base_curls(:, :, :)
        real(dp), intent(in) :: base_divergences(:, :), activation(:, :)
        real(dp), intent(in) :: activation_gradient(:, :, :)
        real(dp), intent(in) :: enriched_values(:, :, :), enriched_curls(:, :, :)
        real(dp), intent(in) :: enriched_divergences(:, :)
        integer :: basis_count, point_count

        basis_count = size(base_values, 2)
        point_count = size(base_values, 3)
        valid = size(base_values, 1) == 3 .and. basis_count > 0 .and. point_count > 0 .and. &
            all(shape(base_curls) == [3, basis_count, point_count]) .and. &
            all(shape(base_divergences) == [basis_count, point_count]) .and. &
            all(shape(activation) == [basis_count, point_count]) .and. &
            all(shape(activation_gradient) == [3, basis_count, point_count]) .and. &
            all(shape(enriched_values) == [3, basis_count, point_count]) .and. &
            all(shape(enriched_curls) == [3, basis_count, point_count]) .and. &
            all(shape(enriched_divergences) == [basis_count, point_count]) .and. &
            finite(base_values) .and. finite(base_curls) .and. finite(base_divergences) .and. &
            finite(activation) .and. finite(activation_gradient) .and. finite(enriched_values) .and. &
            finite(enriched_curls) .and. finite(enriched_divergences)
    end function valid_inputs

    logical function same_shapes(first, first_reference, second, second_reference, &
            third, third_reference, fourth, fourth_reference, fifth, fifth_reference) result(same)
        real(dp), intent(in) :: first(:, :, :), first_reference(:, :, :)
        real(dp), intent(in) :: second(:, :, :), second_reference(:, :, :)
        real(dp), intent(in) :: third(:, :), third_reference(:, :)
        real(dp), intent(in) :: fourth(:, :), fourth_reference(:, :)
        real(dp), intent(in) :: fifth(:, :, :), fifth_reference(:, :, :)

        same = all(shape(first) == shape(first_reference)) .and. all(shape(second) == shape(second_reference)) .and. &
            all(shape(third) == shape(third_reference)) .and. all(shape(fourth) == shape(fourth_reference)) .and. &
            all(shape(fifth) == shape(fifth_reference)) .and. finite(first) .and. finite(second) .and. &
            finite(third) .and. finite(fourth) .and. finite(fifth)
    end function same_shapes

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

    logical function finite_vector(value) result(valid)
        real(dp), intent(in) :: value(:)
        valid = all(ieee_is_finite(value))
    end function finite_vector

    logical function finite_matrix(value) result(valid)
        real(dp), intent(in) :: value(:, :)
        valid = all(ieee_is_finite(value))
    end function finite_matrix

    logical function finite_tensor(value) result(valid)
        real(dp), intent(in) :: value(:, :, :)
        valid = all(ieee_is_finite(value))
    end function finite_tensor

end module fortfem_batched_vector_enrichment_differential_3d
