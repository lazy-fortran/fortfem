module fortfem_enrichment_support_tensor_diagnostics
    !! Point-first tensor Gram contraction for vector/tensor enrichments.
    !!
    !! ``basis_values(point,basis,component)`` is caller-owned and may contain
    !! Piola-mapped or XIGA-enriched vector/tensor samples.  A fixed active mask
    !! selects the current enrichment space; activation is not differentiated.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_enrichment_support_tensor_gram
    public :: evaluate_enrichment_support_tensor_gram_jvp
    public :: evaluate_enrichment_support_tensor_gram_vjp
    public :: evaluate_enrichment_support_tensor_rank_condition

contains

    subroutine evaluate_enrichment_support_tensor_gram( &
            basis_values, component_metric, sample_weights, active_mask, gram, status)
        real(dp), intent(in) :: basis_values(:, :, :), component_metric(:, :, :)
        real(dp), intent(in) :: sample_weights(:)
        logical, intent(in) :: active_mask(:)
        real(dp), intent(out) :: gram(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: point, first_basis, second_basis

        gram = 0.0_dp
        call validate_inputs(basis_values, component_metric, sample_weights, active_mask, gram, status)
        if (status%code /= FORTSPARSE_OK) return
        do first_basis = 1, size(active_mask)
            if (.not. active_mask(first_basis)) cycle
            do second_basis = 1, size(active_mask)
                if (.not. active_mask(second_basis)) cycle
                do point = 1, size(sample_weights)
                    gram(first_basis, second_basis) = gram(first_basis, second_basis) + &
                        sample_weights(point)*dot_product(basis_values(point, first_basis, :), &
                        matmul(component_metric(:, :, point), basis_values(point, second_basis, :)))
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_enrichment_support_tensor_gram

    subroutine evaluate_enrichment_support_tensor_gram_jvp( &
            basis_values, component_metric, sample_weights, active_mask, basis_dot, &
            metric_dot, weights_dot, gram_dot, status)
        real(dp), intent(in) :: basis_values(:, :, :), component_metric(:, :, :)
        real(dp), intent(in) :: sample_weights(:), basis_dot(:, :, :), metric_dot(:, :, :)
        logical, intent(in) :: active_mask(:)
        real(dp), intent(in) :: weights_dot(:)
        real(dp), intent(out) :: gram_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: point, first_basis, second_basis
        real(dp) :: base_pairing, first_dot, metric_pairing, second_dot

        gram_dot = 0.0_dp
        call validate_inputs(basis_values, component_metric, sample_weights, active_mask, gram_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (any(shape(basis_dot) /= shape(basis_values)) .or. &
            any(shape(metric_dot) /= shape(component_metric)) .or. &
            size(weights_dot) /= size(sample_weights) .or. any(.not. ieee_is_finite(basis_dot)) .or. &
            any(.not. ieee_is_finite(metric_dot)) .or. any(.not. ieee_is_finite(weights_dot)) .or. &
            .not. symmetric_metric_direction(metric_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "tensor support Gram JVP has incompatible increments")
            return
        end if
        do first_basis = 1, size(active_mask)
            if (.not. active_mask(first_basis)) cycle
            do second_basis = 1, size(active_mask)
                if (.not. active_mask(second_basis)) cycle
                do point = 1, size(sample_weights)
                    base_pairing = dot_product(basis_values(point, first_basis, :), &
                        matmul(component_metric(:, :, point), basis_values(point, second_basis, :)))
                    first_dot = dot_product(basis_dot(point, first_basis, :), &
                        matmul(component_metric(:, :, point), basis_values(point, second_basis, :)))
                    metric_pairing = dot_product(basis_values(point, first_basis, :), &
                        matmul(metric_dot(:, :, point), basis_values(point, second_basis, :)))
                    second_dot = dot_product(basis_values(point, first_basis, :), &
                        matmul(component_metric(:, :, point), basis_dot(point, second_basis, :)))
                    gram_dot(first_basis, second_basis) = gram_dot(first_basis, second_basis) + &
                        weights_dot(point)*base_pairing + sample_weights(point)*(first_dot + &
                        metric_pairing + second_dot)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_enrichment_support_tensor_gram_jvp

    subroutine evaluate_enrichment_support_tensor_gram_vjp( &
            basis_values, component_metric, sample_weights, active_mask, gram_bar, basis_bar, &
            metric_bar, weights_bar, status)
        real(dp), intent(in) :: basis_values(:, :, :), component_metric(:, :, :)
        real(dp), intent(in) :: sample_weights(:), gram_bar(:, :)
        logical, intent(in) :: active_mask(:)
        real(dp), intent(out) :: basis_bar(:, :, :), metric_bar(:, :, :), weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: point, first_basis, second_basis, first_component, second_component
        real(dp) :: metric_second(size(basis_values, 3)), metric_transpose_first(size(basis_values, 3))
        real(dp) :: pairing

        basis_bar = 0.0_dp
        metric_bar = 0.0_dp
        weights_bar = 0.0_dp
        call validate_inputs(basis_values, component_metric, sample_weights, active_mask, gram_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        if (any(shape(basis_bar) /= shape(basis_values)) .or. &
            any(shape(metric_bar) /= shape(component_metric)) .or. size(weights_bar) /= size(sample_weights) .or. &
            any(.not. ieee_is_finite(gram_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "tensor support Gram VJP has incompatible cotangents")
            return
        end if
        do point = 1, size(sample_weights)
            do first_basis = 1, size(active_mask)
                if (.not. active_mask(first_basis)) cycle
                do second_basis = 1, size(active_mask)
                    if (.not. active_mask(second_basis)) cycle
                    metric_second = matmul(component_metric(:, :, point), basis_values(point, second_basis, :))
                    metric_transpose_first = matmul(transpose(component_metric(:, :, point)), &
                        basis_values(point, first_basis, :))
                    basis_bar(point, first_basis, :) = basis_bar(point, first_basis, :) + &
                        sample_weights(point)*gram_bar(first_basis, second_basis)*metric_second
                    basis_bar(point, second_basis, :) = basis_bar(point, second_basis, :) + &
                        sample_weights(point)*gram_bar(first_basis, second_basis)*metric_transpose_first
                    pairing = dot_product(basis_values(point, first_basis, :), metric_second)
                    weights_bar(point) = weights_bar(point) + gram_bar(first_basis, second_basis)*pairing
                    do first_component = 1, size(basis_values, 3)
                        do second_component = 1, size(basis_values, 3)
                            metric_bar(first_component, second_component, point) = &
                                metric_bar(first_component, second_component, point) + sample_weights(point)* &
                                gram_bar(first_basis, second_basis)*basis_values(point, first_basis, first_component)* &
                                basis_values(point, second_basis, second_component)
                        end do
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_enrichment_support_tensor_gram_vjp

    subroutine evaluate_enrichment_support_tensor_rank_condition( &
            gram, active_mask, tolerance, rank, condition_estimate, eigenvalues, status)
        real(dp), intent(in) :: gram(:, :), tolerance
        logical, intent(in) :: active_mask(:)
        integer, intent(out) :: rank
        real(dp), intent(out) :: condition_estimate, eigenvalues(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: active_gram(:, :), values(:)
        integer, allocatable :: active_indices(:)
        integer :: active_count, first, second, index
        real(dp) :: scale

        rank = 0
        condition_estimate = 0.0_dp
        eigenvalues = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "tensor support rank diagnostic has incompatible arrays")
        if (size(gram, 1) < 1 .or. size(gram, 2) /= size(gram, 1) .or. &
            size(active_mask) /= size(gram, 1) .or. size(eigenvalues) /= count(active_mask) .or. &
            .not. ieee_is_finite(tolerance) .or. tolerance <= 0.0_dp .or. &
            any(.not. ieee_is_finite(gram))) return
        active_count = count(active_mask)
        if (active_count == 0) then
            call status_set(status, FORTSPARSE_OK, "")
            return
        end if
        allocate(active_indices(active_count), active_gram(active_count, active_count), values(active_count))
        index = 0
        do first = 1, size(active_mask)
            if (active_mask(first)) then
                index = index + 1
                active_indices(index) = first
            end if
        end do
        do first = 1, active_count
            do second = 1, active_count
                active_gram(first, second) = 0.5_dp*(gram(active_indices(first), active_indices(second)) + &
                    gram(active_indices(second), active_indices(first)))
            end do
        end do
        call jacobi_eigenvalues(active_gram, values)
        call sort_descending(values)
        eigenvalues = values
        scale = maxval(abs(values))
        if (scale <= tiny(1.0_dp)) then
            condition_estimate = huge(1.0_dp)
            call status_set(status, FORTSPARSE_OK, "")
            return
        end if
        rank = count(values > tolerance*scale)
        if (rank == active_count .and. values(active_count) > 0.0_dp) then
            condition_estimate = sqrt(values(1)/values(active_count))
        else
            condition_estimate = huge(1.0_dp)
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_enrichment_support_tensor_rank_condition

    subroutine validate_inputs(basis_values, metric, weights, active_mask, gram, status)
        real(dp), intent(in) :: basis_values(:, :, :), metric(:, :, :), weights(:), gram(:, :)
        logical, intent(in) :: active_mask(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: point

        call status_set(status, FORTSPARSE_INVALID_MATRIX, "tensor support Gram has incompatible arrays")
        if (size(basis_values, 1) < 1 .or. size(basis_values, 2) < 1 .or. size(basis_values, 3) < 1 .or. &
            size(metric, 1) /= size(basis_values, 3) .or. size(metric, 2) /= size(basis_values, 3) .or. &
            size(metric, 3) /= size(basis_values, 1) .or. size(weights) /= size(basis_values, 1) .or. &
            size(active_mask) /= size(basis_values, 2) .or. size(gram, 1) /= size(active_mask) .or. &
            size(gram, 2) /= size(active_mask) .or. any(.not. ieee_is_finite(basis_values)) .or. &
            any(.not. ieee_is_finite(metric)) .or. any(.not. ieee_is_finite(weights)) .or. &
            any(weights <= 0.0_dp) .or. any(.not. ieee_is_finite(gram))) return
        do point = 1, size(weights)
            if (.not. metric_is_spd(metric(:, :, point))) return
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_inputs

    logical function symmetric_metric_direction(metric_dot) result(valid)
        real(dp), intent(in) :: metric_dot(:, :, :)
        integer :: point
        real(dp) :: scale

        valid = .true.
        do point = 1, size(metric_dot, 3)
            scale = max(1.0_dp, maxval(abs(metric_dot(:, :, point))))
            if (maxval(abs(metric_dot(:, :, point) - transpose(metric_dot(:, :, point)))) > &
                64.0_dp*epsilon(1.0_dp)*scale) then
                valid = .false.
                return
            end if
        end do
    end function symmetric_metric_direction

    logical function metric_is_spd(matrix) result(valid)
        real(dp), intent(in) :: matrix(:, :)
        real(dp) :: leading, determinant, scale

        valid = .false.
        if (size(matrix, 1) /= size(matrix, 2) .or. size(matrix, 1) < 1) return
        scale = max(1.0_dp, maxval(abs(matrix)))
        if (maxval(abs(matrix - transpose(matrix))) > 64.0_dp*epsilon(1.0_dp)*scale) return
        if (size(matrix, 1) == 1) then
            valid = matrix(1, 1) > 0.0_dp
        else if (size(matrix, 1) == 2) then
            valid = matrix(1, 1) > 0.0_dp .and. matrix(1, 1)*matrix(2, 2) - matrix(1, 2)**2 > 0.0_dp
        else
            leading = matrix(1, 1)
            determinant = matrix(1, 1)*(matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2)) - &
                matrix(1, 2)*(matrix(2, 1)*matrix(3, 3) - matrix(2, 3)*matrix(3, 1)) + &
                matrix(1, 3)*(matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1))
            valid = leading > 0.0_dp .and. matrix(1, 1)*matrix(2, 2) - &
                matrix(1, 2)**2 > 0.0_dp .and. determinant > 0.0_dp
        end if
    end function metric_is_spd

    subroutine jacobi_eigenvalues(matrix, values)
        real(dp), intent(in) :: matrix(:, :)
        real(dp), intent(out) :: values(:)
        real(dp), allocatable :: work(:, :)
        real(dp) :: off, scale, tau, tangent, cosine, sine, app, aqq, apq, akp, akq
        integer :: n, p, q, i, j, sweep

        n = size(matrix, 1)
        allocate(work(n, n))
        work = 0.5_dp*(matrix + transpose(matrix))
        do sweep = 1, max(1, 50*n*n)
            off = 0.0_dp; scale = 0.0_dp; p = 1; q = min(2, n)
            do i = 1, n
                scale = max(scale, abs(work(i, i)))
                do j = i + 1, n
                    if (abs(work(i, j)) > off) then
                        off = abs(work(i, j)); p = i; q = j
                    end if
                end do
            end do
            if (off <= 10.0_dp*epsilon(1.0_dp)*max(1.0_dp, scale)) exit
            apq = work(p, q); app = work(p, p); aqq = work(q, q)
            tau = (aqq - app)/(2.0_dp*apq)
            tangent = sign(1.0_dp, tau)/(abs(tau) + sqrt(1.0_dp + tau*tau))
            cosine = 1.0_dp/sqrt(1.0_dp + tangent*tangent); sine = tangent*cosine
            do i = 1, n
                if (i == p .or. i == q) cycle
                akp = work(i, p); akq = work(i, q)
                work(i, p) = cosine*akp - sine*akq; work(p, i) = work(i, p)
                work(i, q) = sine*akp + cosine*akq; work(q, i) = work(i, q)
            end do
            work(p, p) = cosine*cosine*app - 2.0_dp*sine*cosine*apq + sine*sine*aqq
            work(q, q) = sine*sine*app + 2.0_dp*sine*cosine*apq + cosine*cosine*aqq
            work(p, q) = 0.0_dp; work(q, p) = 0.0_dp
        end do
        values = [(work(i, i), i=1, n)]
    end subroutine jacobi_eigenvalues

    pure subroutine sort_descending(values)
        real(dp), intent(inout) :: values(:)
        real(dp) :: candidate
        integer :: first, second, largest

        do first = 1, size(values) - 1
            largest = first
            do second = first + 1, size(values)
                if (values(second) > values(largest)) largest = second
            end do
            if (largest /= first) then
                candidate = values(first); values(first) = values(largest); values(largest) = candidate
            end if
        end do
    end subroutine sort_descending

end module fortfem_enrichment_support_tensor_diagnostics
