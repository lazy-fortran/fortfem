program test_tetra_feec_exact_sequence
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_tetra_discrete_gradient, evaluate_tetra_lagrange, &
        evaluate_tetra_nedelec_first_kind, initialize_tetra_lagrange, &
        initialize_tetra_nedelec_first_kind, tetra_lagrange_dof_count, &
        tetra_lagrange_t, tetra_nedelec_dof_count, &
        tetra_nedelec_first_kind_t
    use fortfem_kinds, only: dp
    implicit none

    type(tetra_lagrange_t) :: lagrange
    type(tetra_nedelec_first_kind_t) :: nedelec
    real(dp), allocatable :: curls(:, :), gradient_matrix(:, :)
    real(dp), allocatable :: gradients(:, :), nedelec_values(:, :)
    real(dp), allocatable :: values(:)
    real(dp) :: point(3)
    integer :: lagrange_dof, order, status
    logical :: all_passed

    all_passed = .true.
    point = [0.19_dp, 0.23_dp, 0.17_dp]
    do order = 1, 4
        call initialize_tetra_lagrange(order, lagrange, status)
        call initialize_tetra_nedelec_first_kind(order, nedelec, status)
        call build_tetra_discrete_gradient(order, gradient_matrix, status)
        allocate( &
            values(tetra_lagrange_dof_count(lagrange)), &
            gradients(3, tetra_lagrange_dof_count(lagrange)), &
            nedelec_values(3, tetra_nedelec_dof_count(nedelec)), &
            curls(3, tetra_nedelec_dof_count(nedelec)))
        call evaluate_tetra_lagrange( &
            lagrange, point, values, gradients, status)
        call evaluate_tetra_nedelec_first_kind( &
            nedelec, point, nedelec_values, curls, status)
        do lagrange_dof = 1, size(gradient_matrix, 2)
            call record_condition(maxval(abs( &
                matmul(nedelec_values, gradient_matrix(:, lagrange_dof)) - &
                gradients(:, lagrange_dof))) < 3.0e-10_dp, &
                "Discrete tetrahedral gradient reproduces H1 gradients")
        end do
        call record_condition(maxval(abs( &
            matmul(curls, gradient_matrix))) < 3.0e-10_dp, &
            "Tetrahedral discrete curl annihilates discrete gradients")
        deallocate(values, gradients, nedelec_values, curls, gradient_matrix)
    end do

    call build_tetra_discrete_gradient(0, gradient_matrix, status)
    call record_condition(status /= 0, &
        "Tetrahedral discrete gradient rejects order zero")
    call check_summary("Tetrahedral FEEC exact sequence")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_feec_exact_sequence
