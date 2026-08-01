module fortfem_enrichment_support_vector_diagnostics
    !! Metric-aware support Gram matrices for Piola-mapped vector enrichments.
    !!
    !! The caller supplies physical vector values.  They may be covariant or
    !! contravariant Piola values, mapped IGA fields, or any other compatible
    !! vector representation.  The metric and quadrature weights stay explicit
    !! so conditioning is independent of the upstream geometry map.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_enrichment_support_vector_gram
    public :: evaluate_enrichment_support_vector_gram_jvp
    public :: evaluate_enrichment_support_vector_gram_vjp

contains

    subroutine evaluate_enrichment_support_vector_gram( &
            basis_values, metric, sample_weights, active_mask, gram, status)
        real(dp), intent(in) :: basis_values(:, :, :), metric(:, :, :)
        real(dp), intent(in) :: sample_weights(:)
        logical, intent(in) :: active_mask(:)
        real(dp), intent(out) :: gram(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: point, first_basis, second_basis

        gram = 0.0_dp
        call validate_vector_gram_inputs( &
            basis_values, metric, sample_weights, active_mask, gram, status)
        if (status%code /= FORTSPARSE_OK) return

        do first_basis = 1, size(active_mask)
            if (.not. active_mask(first_basis)) cycle
            do second_basis = 1, size(active_mask)
                if (.not. active_mask(second_basis)) cycle
                do point = 1, size(sample_weights)
                    gram(first_basis, second_basis) = &
                        gram(first_basis, second_basis) + &
                        sample_weights(point)*dot_product( &
                        basis_values(:, point, first_basis), &
                        matmul(metric(:, :, point), &
                        basis_values(:, point, second_basis)))
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_enrichment_support_vector_gram

    subroutine evaluate_enrichment_support_vector_gram_jvp( &
            basis_values, metric, sample_weights, active_mask, basis_dot, &
            metric_dot, weights_dot, gram_dot, status)
        real(dp), intent(in) :: basis_values(:, :, :), metric(:, :, :)
        real(dp), intent(in) :: sample_weights(:)
        logical, intent(in) :: active_mask(:)
        real(dp), intent(in) :: basis_dot(:, :, :), metric_dot(:, :, :)
        real(dp), intent(in) :: weights_dot(:)
        real(dp), intent(out) :: gram_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: point, first_basis, second_basis
        real(dp) :: first_value, second_value
        real(dp) :: first_dot, second_dot, metric_value

        gram_dot = 0.0_dp
        call validate_vector_gram_inputs( &
            basis_values, metric, sample_weights, active_mask, gram_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        call validate_vector_gram_direction( &
            basis_values, metric, sample_weights, active_mask, basis_dot, &
            metric_dot, weights_dot, gram_dot, status)
        if (status%code /= FORTSPARSE_OK) return

        do first_basis = 1, size(active_mask)
            if (.not. active_mask(first_basis)) cycle
            do second_basis = 1, size(active_mask)
                if (.not. active_mask(second_basis)) cycle
                do point = 1, size(sample_weights)
                    first_value = dot_product( &
                        basis_values(:, point, first_basis), &
                        matmul(metric(:, :, point), &
                        basis_values(:, point, second_basis)))
                    first_dot = dot_product( &
                        basis_dot(:, point, first_basis), &
                        matmul(metric(:, :, point), &
                        basis_values(:, point, second_basis)))
                    second_dot = dot_product( &
                        basis_values(:, point, first_basis), &
                        matmul(metric_dot(:, :, point), &
                        basis_values(:, point, second_basis)))
                    metric_value = dot_product( &
                        basis_values(:, point, first_basis), &
                        matmul(metric(:, :, point), &
                        basis_dot(:, point, second_basis)))
                    gram_dot(first_basis, second_basis) = &
                        gram_dot(first_basis, second_basis) + &
                        weights_dot(point)*first_value + sample_weights(point)*&
                        (first_dot + second_dot + metric_value)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_enrichment_support_vector_gram_jvp

    subroutine evaluate_enrichment_support_vector_gram_vjp( &
            basis_values, metric, sample_weights, active_mask, gram_bar, &
            basis_bar, metric_bar, weights_bar, status)
        real(dp), intent(in) :: basis_values(:, :, :), metric(:, :, :)
        real(dp), intent(in) :: sample_weights(:)
        logical, intent(in) :: active_mask(:)
        real(dp), intent(in) :: gram_bar(:, :)
        real(dp), intent(out) :: basis_bar(:, :, :), metric_bar(:, :, :)
        real(dp), intent(out) :: weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: point, first_basis, second_basis
        integer :: first_component, second_component
        real(dp) :: metric_times_second(size(basis_values, 1))
        real(dp) :: transpose_metric_times_first(size(basis_values, 1))
        real(dp) :: pairing

        basis_bar = 0.0_dp
        metric_bar = 0.0_dp
        weights_bar = 0.0_dp
        call validate_vector_gram_inputs( &
            basis_values, metric, sample_weights, active_mask, gram_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        call validate_vector_gram_adjoint( &
            basis_values, metric, sample_weights, active_mask, gram_bar, &
            basis_bar, metric_bar, weights_bar, status)
        if (status%code /= FORTSPARSE_OK) return

        do point = 1, size(sample_weights)
            do first_basis = 1, size(active_mask)
                if (.not. active_mask(first_basis)) cycle
                do second_basis = 1, size(active_mask)
                    if (.not. active_mask(second_basis)) cycle
                    metric_times_second = matmul(metric(:, :, point), &
                        basis_values(:, point, second_basis))
                    transpose_metric_times_first = matmul( &
                        transpose(metric(:, :, point)), &
                        basis_values(:, point, first_basis))
                    basis_bar(:, point, first_basis) = &
                        basis_bar(:, point, first_basis) + &
                        sample_weights(point)*gram_bar(first_basis, second_basis)*&
                        metric_times_second
                    basis_bar(:, point, second_basis) = &
                        basis_bar(:, point, second_basis) + &
                        sample_weights(point)*gram_bar(first_basis, second_basis)*&
                        transpose_metric_times_first
                    pairing = dot_product( &
                        basis_values(:, point, first_basis), metric_times_second)
                    weights_bar(point) = weights_bar(point) + &
                        gram_bar(first_basis, second_basis)*pairing
                    do first_component = 1, size(basis_values, 1)
                        do second_component = 1, size(basis_values, 1)
                            metric_bar(first_component, second_component, point) = &
                                metric_bar(first_component, second_component, point) + &
                                sample_weights(point)*gram_bar(first_basis, &
                                second_basis)*basis_values(first_component, point, &
                                first_basis)*basis_values(second_component, point, &
                                second_basis)
                        end do
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_enrichment_support_vector_gram_vjp

    subroutine validate_vector_gram_inputs( &
            basis_values, metric, sample_weights, active_mask, gram, status)
        real(dp), intent(in) :: basis_values(:, :, :), metric(:, :, :)
        real(dp), intent(in) :: sample_weights(:)
        logical, intent(in) :: active_mask(:)
        real(dp), intent(in) :: gram(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: point

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "vector support Gram has incompatible arrays")
        if (size(basis_values, 1) < 2) return
        if (size(basis_values, 1) > 3) return
        if (size(basis_values, 2) < 1) return
        if (size(basis_values, 3) < 1) return
        if (size(metric, 1) /= size(basis_values, 1)) return
        if (size(metric, 2) /= size(basis_values, 1)) return
        if (size(metric, 3) /= size(basis_values, 2)) return
        if (size(sample_weights) /= size(basis_values, 2)) return
        if (size(active_mask) /= size(basis_values, 3)) return
        if (size(gram, 1) /= size(active_mask)) return
        if (size(gram, 2) /= size(active_mask)) return
        if (any(.not. ieee_is_finite(basis_values))) return
        if (any(.not. ieee_is_finite(metric))) return
        if (any(.not. ieee_is_finite(sample_weights))) return
        if (any(sample_weights < 0.0_dp)) return
        if (any(.not. ieee_is_finite(gram))) return
        do point = 1, size(sample_weights)
            if (.not. metric_is_spd(metric(:, :, point))) return
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_vector_gram_inputs

    subroutine validate_vector_gram_direction( &
            basis_values, metric, sample_weights, active_mask, basis_dot, &
            metric_dot, weights_dot, gram_dot, status)
        real(dp), intent(in) :: basis_values(:, :, :), metric(:, :, :)
        real(dp), intent(in) :: sample_weights(:)
        logical, intent(in) :: active_mask(:)
        real(dp), intent(in) :: basis_dot(:, :, :), metric_dot(:, :, :)
        real(dp), intent(in) :: weights_dot(:), gram_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "vector support Gram JVP has incompatible increments")
        if (size(basis_dot, 1) /= size(basis_values, 1)) return
        if (size(basis_dot, 2) /= size(basis_values, 2)) return
        if (size(basis_dot, 3) /= size(basis_values, 3)) return
        if (size(metric_dot, 1) /= size(metric, 1)) return
        if (size(metric_dot, 2) /= size(metric, 2)) return
        if (size(metric_dot, 3) /= size(metric, 3)) return
        if (size(weights_dot) /= size(sample_weights)) return
        if (size(gram_dot, 1) /= size(active_mask)) return
        if (size(gram_dot, 2) /= size(active_mask)) return
        if (any(.not. ieee_is_finite(basis_dot))) return
        if (any(.not. ieee_is_finite(metric_dot))) return
        if (any(.not. ieee_is_finite(weights_dot))) return
        if (.not. symmetric_metric_direction(metric_dot)) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_vector_gram_direction

    subroutine validate_vector_gram_adjoint( &
            basis_values, metric, sample_weights, active_mask, gram_bar, &
            basis_bar, metric_bar, weights_bar, status)
        real(dp), intent(in) :: basis_values(:, :, :), metric(:, :, :)
        real(dp), intent(in) :: sample_weights(:)
        logical, intent(in) :: active_mask(:)
        real(dp), intent(in) :: gram_bar(:, :)
        real(dp), intent(in) :: basis_bar(:, :, :), metric_bar(:, :, :)
        real(dp), intent(in) :: weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "vector support Gram VJP has incompatible cotangents")
        if (size(basis_bar, 1) /= size(basis_values, 1)) return
        if (size(basis_bar, 2) /= size(basis_values, 2)) return
        if (size(basis_bar, 3) /= size(basis_values, 3)) return
        if (size(metric_bar, 1) /= size(metric, 1)) return
        if (size(metric_bar, 2) /= size(metric, 2)) return
        if (size(metric_bar, 3) /= size(metric, 3)) return
        if (size(weights_bar) /= size(sample_weights)) return
        if (any(.not. ieee_is_finite(gram_bar))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_vector_gram_adjoint

    pure logical function metric_is_spd(matrix) result(valid)
        real(dp), intent(in) :: matrix(:, :)
        real(dp) :: scale, tolerance, determinant

        scale = max(1.0_dp, maxval(abs(matrix)))
        tolerance = 128.0_dp*epsilon(1.0_dp)*scale
        valid = symmetric_matrix(matrix)
        if (.not. valid) return
        if (size(matrix, 1) == 2) then
            determinant = matrix(1, 1)*matrix(2, 2) - &
                matrix(1, 2)*matrix(2, 1)
            valid = matrix(1, 1) > tolerance .and. determinant > tolerance
        else
            determinant = matrix(1, 1)*matrix(2, 2)*matrix(3, 3) + &
                2.0_dp*matrix(1, 2)*matrix(1, 3)*matrix(2, 3) - &
                matrix(1, 1)*matrix(2, 3)**2 - &
                matrix(2, 2)*matrix(1, 3)**2 - &
                matrix(3, 3)*matrix(1, 2)**2
            valid = matrix(1, 1) > tolerance .and. &
                matrix(1, 1)*matrix(2, 2) - matrix(1, 2)**2 > tolerance .and. &
                determinant > tolerance
        end if
    end function metric_is_spd

    pure logical function symmetric_metric_direction(matrix) result(valid)
        real(dp), intent(in) :: matrix(:, :, :)
        integer :: point

        valid = .true.
        do point = 1, size(matrix, 3)
            if (.not. symmetric_matrix(matrix(:, :, point))) then
                valid = .false.
                return
            end if
        end do
    end function symmetric_metric_direction

    pure logical function symmetric_matrix(matrix) result(valid)
        real(dp), intent(in) :: matrix(:, :)
        real(dp) :: scale, tolerance

        scale = max(1.0_dp, maxval(abs(matrix)))
        tolerance = 128.0_dp*epsilon(1.0_dp)*scale
        valid = maxval(abs(matrix - transpose(matrix))) <= tolerance
    end function symmetric_matrix

end module fortfem_enrichment_support_vector_diagnostics
