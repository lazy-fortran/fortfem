program test_triangle_feec_exact_sequence
    use check, only: check_condition, check_summary
    use fortfem_feec, only: build_triangle_discrete_gradient, &
        evaluate_triangle_nedelec_first_kind, &
        initialize_triangle_nedelec_first_kind, triangle_duffy_quadrature, &
        triangle_nedelec_first_kind_t
    use fortfem_kinds, only: dp
    implicit none

    type(triangle_nedelec_first_kind_t) :: nedelec_basis
    real(dp), allocatable :: curls(:), eta(:), gradient_matrix(:, :)
    real(dp), allocatable :: nedelec_values(:, :), weights(:), xi(:)
    real(dp), allocatable :: zero_curls(:)
    integer :: expected_columns, expected_rows, node, order, status
    logical :: all_passed

    all_passed = .true.
    do order = 1, 4
        call build_triangle_discrete_gradient(order, gradient_matrix, status)
        expected_columns = (order + 1) * (order + 2) / 2
        expected_rows = order * (order + 2)
        call record_condition(status == 0 .and. &
            size(gradient_matrix, 1) == expected_rows .and. &
            size(gradient_matrix, 2) == expected_columns, &
            "Discrete triangle gradient has exact FEEC dimensions")
        call record_condition(maxval(abs(matmul( &
            gradient_matrix, [(1.0_dp, node=1, expected_columns)]))) < &
            2.0e-11_dp, "Discrete gradient annihilates constants")
        call record_condition(matrix_rank(gradient_matrix, 2.0e-10_dp) == &
            expected_columns - 1, &
            "Discrete gradient has only the constant nullspace")

        call initialize_triangle_nedelec_first_kind( &
            order, nedelec_basis, status)
        allocate(curls(expected_rows), nedelec_values(2, expected_rows))
        allocate(zero_curls(expected_columns))
        call triangle_duffy_quadrature( &
            2 * order, xi, eta, weights, status)
        do node = 1, size(xi)
            call evaluate_triangle_nedelec_first_kind( &
                nedelec_basis, xi(node), eta(node), nedelec_values, &
                curls, status)
            zero_curls = matmul(curls, gradient_matrix)
            call record_condition(maxval(abs(zero_curls)) < 3.0e-9_dp, &
                "Discrete triangle complex satisfies curl grad equals zero")
        end do
        deallocate(curls, nedelec_values, zero_curls, xi, eta, weights)
        deallocate(gradient_matrix)
    end do

    call build_triangle_discrete_gradient(0, gradient_matrix, status)
    call record_condition(status /= 0, &
        "Discrete gradient rejects order zero")

    call check_summary("Arbitrary-order triangle FEEC exact sequence")
    if (.not. all_passed) error stop 1

contains

    function matrix_rank(matrix, tolerance) result(rank)
        real(dp), intent(in) :: matrix(:, :)
        real(dp), intent(in) :: tolerance
        integer :: rank

        real(dp), allocatable :: work(:, :), temporary_row(:)
        real(dp) :: pivot_value
        integer :: column, pivot, row

        allocate(work(size(matrix, 1), size(matrix, 2)))
        allocate(temporary_row(size(matrix, 2)))
        work = matrix
        rank = 0
        do column = 1, size(work, 2)
            pivot = 0
            pivot_value = tolerance
            do row = rank + 1, size(work, 1)
                if (abs(work(row, column)) > pivot_value) then
                    pivot = row
                    pivot_value = abs(work(row, column))
                end if
            end do
            if (pivot == 0) cycle
            rank = rank + 1
            if (pivot /= rank) then
                temporary_row = work(rank, :)
                work(rank, :) = work(pivot, :)
                work(pivot, :) = temporary_row
            end if
            work(rank, :) = work(rank, :) / work(rank, column)
            do row = 1, size(work, 1)
                if (row == rank) cycle
                work(row, :) = work(row, :) - &
                    work(row, column) * work(rank, :)
            end do
        end do
    end function matrix_rank

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_triangle_feec_exact_sequence
