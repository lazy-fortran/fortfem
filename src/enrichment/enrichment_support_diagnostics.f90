module fortfem_enrichment_support_diagnostics
    !! Weighted support Gram matrices and fixed-activation rank diagnostics.
    !!
    !! The active mask is a discrete topology/space decision owned by the
    !! caller.  The Gram value, JVP, and VJP treat that mask as fixed, so the
    !! differentiable contract applies between enrichment-space rebuilds.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_enrichment_support_gram
    public :: evaluate_enrichment_support_gram_jvp
    public :: evaluate_enrichment_support_gram_vjp
    public :: evaluate_enrichment_support_rank_condition

contains

    subroutine evaluate_enrichment_support_gram( &
            basis_values, sample_weights, active_mask, gram, status)
        real(dp), intent(in) :: basis_values(:, :), sample_weights(:)
        logical, intent(in) :: active_mask(:)
        real(dp), intent(out) :: gram(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: point, first_basis, second_basis

        gram = 0.0_dp
        call validate_gram_inputs( &
            basis_values, sample_weights, active_mask, gram, status)
        if (status%code /= FORTSPARSE_OK) return

        do first_basis = 1, size(basis_values, 2)
            if (.not. active_mask(first_basis)) cycle
            do second_basis = 1, size(basis_values, 2)
                if (.not. active_mask(second_basis)) cycle
                do point = 1, size(basis_values, 1)
                    gram(first_basis, second_basis) = &
                        gram(first_basis, second_basis) + sample_weights(point)* &
                        basis_values(point, first_basis)* &
                        basis_values(point, second_basis)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_enrichment_support_gram

    subroutine evaluate_enrichment_support_gram_jvp( &
            basis_values, sample_weights, active_mask, basis_dot, weights_dot, &
            gram_dot, status)
        real(dp), intent(in) :: basis_values(:, :), sample_weights(:)
        logical, intent(in) :: active_mask(:)
        real(dp), intent(in) :: basis_dot(:, :), weights_dot(:)
        real(dp), intent(out) :: gram_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: point, first_basis, second_basis

        gram_dot = 0.0_dp
        call validate_gram_inputs( &
            basis_values, sample_weights, active_mask, gram_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(basis_dot, 1) /= size(basis_values, 1) .or. &
            size(basis_dot, 2) /= size(basis_values, 2) .or. &
            size(weights_dot) /= size(sample_weights) .or. &
            any(.not. ieee_is_finite(basis_dot)) .or. &
            any(.not. ieee_is_finite(weights_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support Gram JVP has incompatible increments")
            return
        end if

        do first_basis = 1, size(basis_values, 2)
            if (.not. active_mask(first_basis)) cycle
            do second_basis = 1, size(basis_values, 2)
                if (.not. active_mask(second_basis)) cycle
                do point = 1, size(basis_values, 1)
                    gram_dot(first_basis, second_basis) = &
                        gram_dot(first_basis, second_basis) + weights_dot(point)* &
                        basis_values(point, first_basis)* &
                        basis_values(point, second_basis) + sample_weights(point)* &
                        (basis_dot(point, first_basis)*basis_values(point, second_basis) + &
                        basis_values(point, first_basis)*basis_dot(point, second_basis))
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_enrichment_support_gram_jvp

    subroutine evaluate_enrichment_support_gram_vjp( &
            basis_values, sample_weights, active_mask, gram_bar, basis_bar, &
            weights_bar, status)
        real(dp), intent(in) :: basis_values(:, :), sample_weights(:)
        logical, intent(in) :: active_mask(:)
        real(dp), intent(in) :: gram_bar(:, :)
        real(dp), intent(out) :: basis_bar(:, :), weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: point, first_basis, second_basis

        basis_bar = 0.0_dp
        weights_bar = 0.0_dp
        call validate_gram_inputs( &
            basis_values, sample_weights, active_mask, gram_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(basis_bar, 1) /= size(basis_values, 1) .or. &
            size(basis_bar, 2) /= size(basis_values, 2) .or. &
            size(weights_bar) /= size(sample_weights) .or. &
            any(.not. ieee_is_finite(gram_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support Gram VJP has incompatible cotangents")
            return
        end if

        do point = 1, size(basis_values, 1)
            do first_basis = 1, size(basis_values, 2)
                if (.not. active_mask(first_basis)) cycle
                do second_basis = 1, size(basis_values, 2)
                    if (.not. active_mask(second_basis)) cycle
                    weights_bar(point) = weights_bar(point) + &
                        gram_bar(first_basis, second_basis)* &
                        basis_values(point, first_basis)* &
                        basis_values(point, second_basis)
                    basis_bar(point, first_basis) = &
                        basis_bar(point, first_basis) + sample_weights(point)* &
                        (gram_bar(first_basis, second_basis)* &
                        basis_values(point, second_basis) + &
                        gram_bar(second_basis, first_basis)* &
                        basis_values(point, second_basis))
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_enrichment_support_gram_vjp

    subroutine evaluate_enrichment_support_rank_condition( &
            gram, active_mask, tolerance, rank, condition_estimate, &
            eigenvalues, status)
        real(dp), intent(in) :: gram(:, :), tolerance
        logical, intent(in) :: active_mask(:)
        integer, intent(out) :: rank
        real(dp), intent(out) :: condition_estimate, eigenvalues(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: active_gram(:, :), work_eigenvalues(:)
        real(dp) :: scale, threshold, largest, smallest
        integer, allocatable :: active_indices(:)
        integer :: active_count, first_basis, second_basis, index

        rank = 0
        condition_estimate = 0.0_dp
        eigenvalues = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "support rank diagnostic has incompatible arrays")
        if (size(gram, 1) < 1 .or. size(gram, 2) /= size(gram, 1) .or. &
            size(active_mask) /= size(gram, 1) .or. tolerance <= 0.0_dp .or. &
            .not. ieee_is_finite(tolerance) .or. &
            any(.not. ieee_is_finite(gram))) return

        active_count = count(active_mask)
        if (size(eigenvalues) /= active_count) return
        if (active_count == 0) then
            call status_set(status, FORTSPARSE_OK, "")
            return
        end if

        allocate(active_indices(active_count), active_gram(active_count, active_count), &
            work_eigenvalues(active_count))
        index = 0
        do first_basis = 1, size(active_mask)
            if (.not. active_mask(first_basis)) cycle
            index = index + 1
            active_indices(index) = first_basis
        end do
        do first_basis = 1, active_count
            do second_basis = 1, active_count
                active_gram(first_basis, second_basis) = 0.5_dp*( &
                    gram(active_indices(first_basis), active_indices(second_basis)) + &
                    gram(active_indices(second_basis), active_indices(first_basis)))
            end do
        end do

        call symmetric_jacobi_eigenvalues(active_gram, work_eigenvalues)
        call sort_descending(work_eigenvalues)
        eigenvalues = work_eigenvalues
        scale = maxval(abs(work_eigenvalues))
        if (scale <= tiny(1.0_dp)) then
            call status_set(status, FORTSPARSE_OK, "")
            return
        end if
        threshold = tolerance*scale
        rank = count(work_eigenvalues > threshold)
        if (rank == active_count) then
            largest = work_eigenvalues(1)
            smallest = work_eigenvalues(active_count)
            condition_estimate = sqrt(max(largest, 0.0_dp)/smallest)
        else
            condition_estimate = huge(1.0_dp)
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_enrichment_support_rank_condition

    subroutine validate_gram_inputs( &
            basis_values, sample_weights, active_mask, gram, status)
        real(dp), intent(in) :: basis_values(:, :), sample_weights(:)
        logical, intent(in) :: active_mask(:)
        real(dp), intent(in) :: gram(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "support Gram matrix has incompatible arrays")
        if (size(basis_values, 1) < 1 .or. size(basis_values, 2) < 1 .or. &
            size(sample_weights) /= size(basis_values, 1) .or. &
            size(active_mask) /= size(basis_values, 2) .or. &
            size(gram, 1) /= size(basis_values, 2) .or. &
            size(gram, 2) /= size(basis_values, 2)) return
        if (any(.not. ieee_is_finite(basis_values)) .or. &
            any(.not. ieee_is_finite(sample_weights)) .or. &
            any(sample_weights < 0.0_dp) .or. &
            any(.not. ieee_is_finite(gram))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support Gram matrix received invalid data")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_gram_inputs

    subroutine symmetric_jacobi_eigenvalues(matrix, values)
        real(dp), intent(in) :: matrix(:, :)
        real(dp), intent(out) :: values(:)
        real(dp), allocatable :: work(:, :)
        real(dp) :: diagonal_scale, off_diagonal, tau, tangent
        real(dp) :: cosine, sine, app, aqq, apq, akp, akq
        integer :: n, p, q, j, k, sweep, max_sweeps

        n = size(matrix, 1)
        allocate(work(n, n))
        work = 0.5_dp*(matrix + transpose(matrix))
        max_sweeps = max(1, 50*n*n)
        do sweep = 1, max_sweeps
            off_diagonal = 0.0_dp
            diagonal_scale = 0.0_dp
            p = 1
            q = min(2, n)
            do k = 1, n
                diagonal_scale = max(diagonal_scale, abs(work(k, k)))
                do j = k + 1, n
                    if (abs(work(k, j)) > off_diagonal) then
                        off_diagonal = abs(work(k, j))
                        p = k
                        q = j
                    end if
                end do
            end do
            diagonal_scale = max(1.0_dp, diagonal_scale)
            if (off_diagonal <= 10.0_dp*epsilon(1.0_dp)*diagonal_scale) exit
            apq = work(p, q)
            app = work(p, p)
            aqq = work(q, q)
            tau = (aqq - app)/(2.0_dp*apq)
            tangent = sign(1.0_dp, tau)/(abs(tau) + sqrt(1.0_dp + tau*tau))
            cosine = 1.0_dp/sqrt(1.0_dp + tangent*tangent)
            sine = tangent*cosine
            do k = 1, n
                if (k == p .or. k == q) cycle
                akp = work(k, p)
                akq = work(k, q)
                work(k, p) = cosine*akp - sine*akq
                work(p, k) = work(k, p)
                work(k, q) = sine*akp + cosine*akq
                work(q, k) = work(k, q)
            end do
            work(p, p) = cosine*cosine*app - 2.0_dp*sine*cosine*apq + &
                sine*sine*aqq
            work(q, q) = sine*sine*app + 2.0_dp*sine*cosine*apq + &
                cosine*cosine*aqq
            work(p, q) = 0.0_dp
            work(q, p) = 0.0_dp
        end do
        values = [(work(k, k), k=1, n)]
    end subroutine symmetric_jacobi_eigenvalues

    pure subroutine sort_descending(values)
        real(dp), intent(inout) :: values(:)
        real(dp) :: candidate
        integer :: first, second, largest_index

        do first = 1, size(values) - 1
            largest_index = first
            do second = first + 1, size(values)
                if (values(second) > values(largest_index)) largest_index = second
            end do
            if (largest_index == first) cycle
            candidate = values(first)
            values(first) = values(largest_index)
            values(largest_index) = candidate
        end do
    end subroutine sort_descending

end module fortfem_enrichment_support_diagnostics
