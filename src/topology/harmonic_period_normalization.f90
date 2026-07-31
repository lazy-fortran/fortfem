module fortfem_harmonic_period_normalization
    !! Fixed-topology period normalization for harmonic one-form bases.
    !!
    !! If H stores harmonic one-forms on edges and C stores oriented cycles,
    !! the period matrix is P = C^T H.  This primitive returns A from
    !! P A = T and the normalized basis H A, where T is caller-owned target
    !! period data.  Cycle topology is frozen; geometry, flux units, and
    !! physical interpretation of T remain outside this module.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve
    implicit none
    private

    public :: normalize_harmonic_one_forms
    public :: normalize_harmonic_one_forms_jvp
    public :: normalize_harmonic_one_forms_vjp

contains

    subroutine normalize_harmonic_one_forms( &
            harmonic_forms, period_cycles, target_periods, normalized_forms, &
            normalization_matrix, status)
        real(dp), intent(in) :: harmonic_forms(:, :), period_cycles(:, :)
        real(dp), intent(in) :: target_periods(:, :)
        real(dp), intent(out) :: normalized_forms(:, :), normalization_matrix(:, :)
        integer, intent(out) :: status

        integer :: form_count, info
        real(dp), allocatable :: period_matrix(:, :)

        normalized_forms = 0.0_dp
        normalization_matrix = 0.0_dp
        if (.not. validate_inputs( &
            harmonic_forms, period_cycles, target_periods, normalized_forms, &
            normalization_matrix)) then
            status = 1
            return
        end if
        form_count = size(harmonic_forms, 2)
        allocate(period_matrix(form_count, form_count))
        period_matrix = matmul(transpose(period_cycles), harmonic_forms)
        call dense_solve(period_matrix, target_periods, normalization_matrix, info)
        if (info /= 0) then
            status = 2
            return
        end if
        normalized_forms = matmul(harmonic_forms, normalization_matrix)
        status = 0
    end subroutine normalize_harmonic_one_forms

    subroutine normalize_harmonic_one_forms_jvp( &
            harmonic_forms, period_cycles, target_periods, harmonic_forms_dot, &
            period_cycles_dot, target_periods_dot, normalized_forms_dot, &
            normalization_matrix_dot, status)
        real(dp), intent(in) :: harmonic_forms(:, :), period_cycles(:, :)
        real(dp), intent(in) :: target_periods(:, :)
        real(dp), intent(in) :: harmonic_forms_dot(:, :), period_cycles_dot(:, :)
        real(dp), intent(in) :: target_periods_dot(:, :)
        real(dp), intent(out) :: normalized_forms_dot(:, :)
        real(dp), intent(out) :: normalization_matrix_dot(:, :)
        integer, intent(out) :: status

        integer :: form_count, info
        real(dp), allocatable :: period_matrix(:, :), period_matrix_dot(:, :)
        real(dp), allocatable :: normalization_matrix(:, :), rhs_dot(:, :)

        normalized_forms_dot = 0.0_dp
        normalization_matrix_dot = 0.0_dp
        if (.not. validate_inputs( &
            harmonic_forms, period_cycles, target_periods, normalized_forms_dot, &
            normalization_matrix_dot)) then
            status = 1
            return
        end if
        form_count = size(harmonic_forms, 2)
        if (.not. validate_direction( &
            harmonic_forms_dot, period_cycles_dot, target_periods_dot, &
            size(harmonic_forms, 1), form_count)) then
            status = 1
            return
        end if
        allocate(period_matrix(form_count, form_count), &
            period_matrix_dot(form_count, form_count), &
            normalization_matrix(form_count, form_count), &
            rhs_dot(form_count, form_count))
        period_matrix = matmul(transpose(period_cycles), harmonic_forms)
        period_matrix_dot = matmul(transpose(period_cycles_dot), harmonic_forms) + &
            matmul(transpose(period_cycles), harmonic_forms_dot)
        call dense_solve(period_matrix, target_periods, normalization_matrix, info)
        if (info /= 0) then
            status = 2
            return
        end if
        rhs_dot = target_periods_dot - matmul( &
            period_matrix_dot, normalization_matrix)
        call dense_solve(period_matrix, rhs_dot, normalization_matrix_dot, info)
        if (info /= 0) then
            status = 2
            return
        end if
        normalized_forms_dot = matmul(harmonic_forms_dot, normalization_matrix) + &
            matmul(harmonic_forms, normalization_matrix_dot)
        status = 0
    end subroutine normalize_harmonic_one_forms_jvp

    subroutine normalize_harmonic_one_forms_vjp( &
            harmonic_forms, period_cycles, target_periods, normalized_forms, &
            normalization_matrix, normalized_forms_bar, normalization_matrix_bar, &
            harmonic_forms_bar, period_cycles_bar, target_periods_bar, status)
        real(dp), intent(in) :: harmonic_forms(:, :), period_cycles(:, :)
        real(dp), intent(in) :: target_periods(:, :), normalized_forms(:, :)
        real(dp), intent(in) :: normalization_matrix(:, :)
        real(dp), intent(in) :: normalized_forms_bar(:, :)
        real(dp), intent(in) :: normalization_matrix_bar(:, :)
        real(dp), intent(out) :: harmonic_forms_bar(:, :), period_cycles_bar(:, :)
        real(dp), intent(out) :: target_periods_bar(:, :)
        integer, intent(out) :: status

        integer :: form_count, info
        real(dp), allocatable :: period_matrix(:, :), normalization_bar(:, :)
        real(dp), allocatable :: period_matrix_bar(:, :), period_adjoint(:, :)

        harmonic_forms_bar = 0.0_dp
        period_cycles_bar = 0.0_dp
        target_periods_bar = 0.0_dp
        if (.not. validate_inputs( &
            harmonic_forms, period_cycles, target_periods, normalized_forms, &
            normalization_matrix)) then
            status = 1
            return
        end if
        form_count = size(harmonic_forms, 2)
        if (size(normalized_forms_bar, 1) /= size(harmonic_forms, 1) .or. &
            size(normalized_forms_bar, 2) /= form_count .or. &
            size(normalization_matrix_bar, 1) /= form_count .or. &
            size(normalization_matrix_bar, 2) /= form_count .or. &
            .not. finite_matrix(normalized_forms_bar) .or. &
            .not. finite_matrix(normalization_matrix_bar)) then
            status = 1
            return
        end if
        allocate(period_matrix(form_count, form_count), &
            normalization_bar(form_count, form_count), &
            period_matrix_bar(form_count, form_count), &
            period_adjoint(form_count, form_count))
        period_matrix = matmul(transpose(period_cycles), harmonic_forms)
        normalization_bar = normalization_matrix_bar + matmul( &
            transpose(harmonic_forms), normalized_forms_bar)
        call dense_solve(transpose(period_matrix), normalization_bar, &
            period_adjoint, info)
        if (info /= 0) then
            status = 2
            return
        end if
        target_periods_bar = period_adjoint
        period_matrix_bar = -matmul(period_adjoint, transpose(normalization_matrix))
        harmonic_forms_bar = matmul(normalized_forms_bar, &
            transpose(normalization_matrix)) + matmul(period_cycles, period_matrix_bar)
        period_cycles_bar = matmul(harmonic_forms, transpose(period_matrix_bar))
        status = 0
    end subroutine normalize_harmonic_one_forms_vjp

    logical function validate_inputs( &
            harmonic_forms, period_cycles, target_periods, normalized_forms, &
            normalization_matrix) result(valid)
        real(dp), intent(in) :: harmonic_forms(:, :), period_cycles(:, :)
        real(dp), intent(in) :: target_periods(:, :)
        real(dp), intent(in) :: normalized_forms(:, :), normalization_matrix(:, :)
        integer :: edge_count, form_count

        valid = .false.
        edge_count = size(harmonic_forms, 1)
        form_count = size(harmonic_forms, 2)
        if (edge_count < 1 .or. form_count < 1 .or. &
            size(period_cycles, 1) /= edge_count .or. &
            size(period_cycles, 2) /= form_count .or. &
            size(target_periods, 1) /= form_count .or. &
            size(target_periods, 2) /= form_count .or. &
            size(normalized_forms, 1) /= edge_count .or. &
            size(normalized_forms, 2) /= form_count .or. &
            size(normalization_matrix, 1) /= form_count .or. &
            size(normalization_matrix, 2) /= form_count) return
        if (.not. finite_matrix(harmonic_forms) .or. &
            .not. finite_matrix(period_cycles) .or. &
            .not. finite_matrix(target_periods) .or. &
            .not. finite_matrix(normalized_forms) .or. &
            .not. finite_matrix(normalization_matrix)) return
        valid = .true.
    end function validate_inputs

    logical function validate_direction( &
            harmonic_forms_dot, period_cycles_dot, target_periods_dot, &
            edge_count, form_count) result(valid)
        real(dp), intent(in) :: harmonic_forms_dot(:, :), period_cycles_dot(:, :)
        real(dp), intent(in) :: target_periods_dot(:, :)
        integer, intent(in) :: edge_count, form_count

        valid = size(harmonic_forms_dot, 1) == edge_count .and. &
            size(harmonic_forms_dot, 2) == form_count .and. &
            size(period_cycles_dot, 1) == edge_count .and. &
            size(period_cycles_dot, 2) == form_count .and. &
            size(target_periods_dot, 1) == form_count .and. &
            size(target_periods_dot, 2) == form_count .and. &
            finite_matrix(harmonic_forms_dot) .and. &
            finite_matrix(period_cycles_dot) .and. &
            finite_matrix(target_periods_dot)
    end function validate_direction

    pure logical function finite_matrix(values) result(valid)
        real(dp), intent(in) :: values(:, :)

        valid = all(ieee_is_finite(values))
    end function finite_matrix

end module fortfem_harmonic_period_normalization
