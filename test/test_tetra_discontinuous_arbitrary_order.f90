program test_tetra_discontinuous_arbitrary_order
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_tetra_discontinuous, &
        initialize_tetra_discontinuous, &
        tetra_discontinuous_dof_count, tetra_discontinuous_t, &
        tetra_duffy_quadrature
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: expected_counts(0:5) = [1, 4, 10, 20, 35, 56]
    type(tetra_discontinuous_t) :: basis, copied_basis
    real(dp), allocatable :: values(:), x(:), y(:), z(:), weights(:)
    real(dp) :: exact_integral, numerical_integral
    integer :: basis_id, degree, exponent_sum, status
    logical :: all_passed

    all_passed = .true.
    do degree = 0, 5
        call initialize_tetra_discontinuous(degree, basis, status)
        call record_condition(status == 0 .and. &
            tetra_discontinuous_dof_count(basis) == expected_counts(degree), &
            "Tetrahedral discontinuous dimension matches P_k")
        copied_basis = basis
        call record_condition( &
            tetra_discontinuous_dof_count(copied_basis) == &
            expected_counts(degree), &
            "Tetrahedral discontinuous basis supports deep copy")

        call tetra_duffy_quadrature(2*degree + 2, x, y, z, weights, status)
        call record_condition(status == 0, &
            "Tetrahedral discontinuous test quadrature succeeds")
        allocate(values(expected_counts(degree)))
        numerical_integral = 0.0_dp
        do basis_id = 1, expected_counts(degree)
            call integrate_basis( &
                basis, basis_id, x, y, z, weights, numerical_integral)
            exponent_sum = sum(basis%exponents(:, basis_id))
            exact_integral = real( &
                factorial(basis%exponents(1, basis_id))* &
                factorial(basis%exponents(2, basis_id))* &
                factorial(basis%exponents(3, basis_id)), dp)/ &
                real(factorial(exponent_sum + 3), dp)
            call record_condition(abs(numerical_integral - exact_integral) < &
                2.0e-13_dp, &
                "Tetrahedral discontinuous monomial has exact integral")
        end do
        deallocate(values)
    end do

    call initialize_tetra_discontinuous(-1, basis, status)
    call record_condition(status /= 0, &
        "Tetrahedral discontinuous basis rejects negative degree")
    call initialize_tetra_discontinuous(6, basis, status)
    call record_condition(status /= 0, &
        "Tetrahedral discontinuous basis rejects unsupported degree")

    call check_summary("Tetrahedral discontinuous arbitrary order")
    if (.not. all_passed) error stop 1

contains

    subroutine integrate_basis( &
            selected_basis, selected_id, x, y, z, weights, integral)
        type(tetra_discontinuous_t), intent(in) :: selected_basis
        integer, intent(in) :: selected_id
        real(dp), intent(in) :: x(:), y(:), z(:), weights(:)
        real(dp), intent(out) :: integral

        integer :: point, local_status

        integral = 0.0_dp
        do point = 1, size(weights)
            call evaluate_tetra_discontinuous( &
                selected_basis, x(point), y(point), z(point), values, &
                local_status)
            if (local_status /= 0) error stop "DG basis evaluation failed"
            integral = integral + weights(point)*values(selected_id)
        end do
    end subroutine integrate_basis

    pure integer function factorial(value) result(product)
        integer, intent(in) :: value
        integer :: factor

        product = 1
        do factor = 2, value
            product = product*factor
        end do
    end function factorial

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition
end program test_tetra_discontinuous_arbitrary_order
