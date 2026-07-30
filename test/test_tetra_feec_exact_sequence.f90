program test_tetra_feec_exact_sequence
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_tetra_discrete_curl, build_tetra_discrete_gradient, &
        evaluate_tetra_lagrange, evaluate_tetra_rt, &
        evaluate_tetra_nedelec_first_kind, initialize_tetra_lagrange, &
        initialize_tetra_nedelec_first_kind, initialize_tetra_rt, &
        tetra_lagrange_dof_count, tetra_lagrange_t, &
        tetra_nedelec_dof_count, tetra_nedelec_first_kind_t, &
        tetra_rt_dof_count, tetra_rt_t
    use fortfem_kinds, only: dp
    implicit none

    type(tetra_lagrange_t) :: lagrange
    type(tetra_nedelec_first_kind_t) :: nedelec
    type(tetra_rt_t) :: rt
    real(dp), allocatable :: curl_matrix(:, :), curls(:, :)
    real(dp), allocatable :: gradient_matrix(:, :)
    real(dp), allocatable :: gradients(:, :), nedelec_values(:, :)
    real(dp), allocatable :: rt_divergences(:), rt_values(:, :)
    real(dp), allocatable :: values(:)
    real(dp) :: point(3)
    integer :: lagrange_dof, order, status
    logical :: all_passed

    all_passed = .true.
    point = [0.19_dp, 0.23_dp, 0.17_dp]
    do order = 1, 4
        call initialize_tetra_lagrange(order, lagrange, status)
        call initialize_tetra_nedelec_first_kind(order, nedelec, status)
        call initialize_tetra_rt(order - 1, rt, status)
        call build_tetra_discrete_gradient(order, gradient_matrix, status)
        call build_tetra_discrete_curl(order, curl_matrix, status)
        allocate( &
            values(tetra_lagrange_dof_count(lagrange)), &
            gradients(3, tetra_lagrange_dof_count(lagrange)), &
            nedelec_values(3, tetra_nedelec_dof_count(nedelec)), &
            curls(3, tetra_nedelec_dof_count(nedelec)), &
            rt_values(3, tetra_rt_dof_count(rt)), &
            rt_divergences(tetra_rt_dof_count(rt)))
        call evaluate_tetra_lagrange( &
            lagrange, point, values, gradients, status)
        call evaluate_tetra_nedelec_first_kind( &
            nedelec, point, nedelec_values, curls, status)
        call evaluate_tetra_rt( &
            rt, point, rt_values, rt_divergences, status)
        do lagrange_dof = 1, size(gradient_matrix, 2)
            call record_condition(maxval(abs( &
                matmul(nedelec_values, gradient_matrix(:, lagrange_dof)) - &
                gradients(:, lagrange_dof))) < 3.0e-10_dp, &
                "Discrete tetrahedral gradient reproduces H1 gradients")
        end do
        call record_condition(maxval(abs( &
            matmul(curls, gradient_matrix))) < 3.0e-10_dp, &
            "Tetrahedral discrete curl annihilates discrete gradients")
        call record_condition(maxval(abs( &
            matmul(rt_values, curl_matrix) - curls)) < 4.0e-9_dp, &
            "Discrete tetrahedral curl reproduces Nedelec curls")
        call record_condition(maxval(abs( &
            matmul(reshape(rt_divergences, [1, size(rt_divergences)]), &
            curl_matrix))) < 4.0e-9_dp, &
            "Tetrahedral discrete divergence annihilates discrete curls")
        deallocate( &
            values, gradients, nedelec_values, curls, rt_values, &
            rt_divergences, gradient_matrix, curl_matrix)
    end do

    call build_tetra_discrete_gradient(0, gradient_matrix, status)
    call record_condition(status /= 0, &
        "Tetrahedral discrete gradient rejects order zero")
    call build_tetra_discrete_curl(0, curl_matrix, status)
    call record_condition(status /= 0, &
        "Tetrahedral discrete curl rejects order zero")
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
