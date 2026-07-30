program test_tetra_lagrange_arbitrary_order
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_tetra_lagrange, initialize_tetra_lagrange, &
        tetra_lagrange_dof_count, tetra_lagrange_nodes, tetra_lagrange_t
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: expected_counts(0:4) = [1, 4, 10, 20, 35]
    type(tetra_lagrange_t) :: basis
    real(dp), allocatable :: gradients(:, :), nodes(:, :), values(:)
    real(dp), allocatable :: minus_values(:), shifted_gradients(:, :)
    real(dp), allocatable :: shifted_values(:)
    real(dp) :: point(3), reproduced, target
    integer :: component, degree, node, status
    logical :: all_passed

    all_passed = .true.
    do degree = 0, 4
        call initialize_tetra_lagrange(degree, basis, status)
        call tetra_lagrange_nodes(basis, nodes)
        call record_condition(status == 0 .and. &
            tetra_lagrange_dof_count(basis) == expected_counts(degree), &
            "Tetrahedral Lagrange dimension matches P_k")
        allocate(values(size(nodes, 2)), gradients(3, size(nodes, 2)))
        allocate( &
            shifted_values(size(nodes, 2)), &
            minus_values(size(nodes, 2)), &
            shifted_gradients(3, size(nodes, 2)))
        do node = 1, size(nodes, 2)
            call evaluate_tetra_lagrange( &
                basis, nodes(:, node), values, gradients, status)
            values(node) = values(node) - 1.0_dp
            call record_condition(status == 0 .and. &
                maxval(abs(values)) < 3.0e-13_dp, &
                "Tetrahedral Lagrange basis is cardinal at every node")
        end do
        point = [0.17_dp, 0.21_dp, 0.13_dp]
        call evaluate_tetra_lagrange( &
            basis, point, values, gradients, status)
        call record_condition(abs(sum(values) - 1.0_dp) < 3.0e-13_dp .and. &
            maxval(abs(sum(gradients, dim=2))) < 2.0e-12_dp, &
            "Tetrahedral Lagrange basis forms a partition of unity")
        if (degree >= 1) then
            reproduced = sum(values*( &
                2.0_dp*nodes(1, :) - 3.0_dp*nodes(2, :) + &
                0.5_dp*nodes(3, :) + 1.25_dp))
            target = 2.0_dp*point(1) - 3.0_dp*point(2) + &
                0.5_dp*point(3) + 1.25_dp
            call record_condition(abs(reproduced - target) < 3.0e-13_dp, &
                "Tetrahedral Lagrange basis reproduces linear polynomials")
        end if
        do component = 1, 3
            point(component) = point(component) + 1.0e-6_dp
            call evaluate_tetra_lagrange( &
                basis, point, shifted_values, shifted_gradients, status)
            point(component) = point(component) - 2.0e-6_dp
            call evaluate_tetra_lagrange( &
                basis, point, minus_values, shifted_gradients, status)
            point(component) = point(component) + 1.0e-6_dp
            call record_condition(maxval(abs( &
                (shifted_values - minus_values)/(2.0e-6_dp) - &
                gradients(component, :))) < &
                2.0e-7_dp, &
                "Tetrahedral Lagrange gradients match finite differences")
        end do
        deallocate( &
            values, gradients, minus_values, shifted_values, shifted_gradients)
        deallocate(nodes)
    end do

    call initialize_tetra_lagrange(-1, basis, status)
    call record_condition(status /= 0, &
        "Tetrahedral Lagrange basis rejects a negative degree")
    call check_summary("Tetrahedral Lagrange arbitrary order")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_lagrange_arbitrary_order
