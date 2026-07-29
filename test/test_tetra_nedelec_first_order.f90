program test_tetra_nedelec_first_order
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_tetra_nedelec_first_order
    use fortfem_kinds, only: dp
    implicit none

    real(dp) :: curls(3, 6), expected, midpoint(3), tangent(3)
    real(dp) :: values(3, 6), vertices(3, 4)
    integer :: edge, edge_vertices(2, 6), local_status, other_edge
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
    vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
    edge_vertices(:, 1) = [1, 2]
    edge_vertices(:, 2) = [1, 3]
    edge_vertices(:, 3) = [1, 4]
    edge_vertices(:, 4) = [2, 3]
    edge_vertices(:, 5) = [2, 4]
    edge_vertices(:, 6) = [3, 4]

    do edge = 1, 6
        midpoint = 0.5_dp * ( &
            vertices(:, edge_vertices(1, edge)) + &
            vertices(:, edge_vertices(2, edge)))
        tangent = vertices(:, edge_vertices(2, edge)) - &
            vertices(:, edge_vertices(1, edge))
        call evaluate_tetra_nedelec_first_order( &
            midpoint, values, curls, local_status)
        do other_edge = 1, 6
            expected = 0.0_dp
            if (other_edge == edge) expected = 1.0_dp
            call record_condition(abs( &
                dot_product(values(:, other_edge), tangent) - expected) < &
                2.0e-14_dp, &
                "Tetrahedral basis has Kronecker oriented edge moments")
        end do
    end do

    call check_curls_from_finite_differences()

    call check_summary("First-order tetrahedral Nedelec basis")
    if (.not. all_passed) error stop 1

contains

    subroutine check_curls_from_finite_differences()
        real(dp), parameter :: step = 1.0e-6_dp
        real(dp) :: center(3), curls_at_center(3, 6)
        real(dp) :: numerical_curl(3), shifted(3)
        real(dp) :: values_minus(3, 6), values_plus(3, 6)
        integer :: basis, direction, status

        center = [0.2_dp, 0.2_dp, 0.2_dp]
        call evaluate_tetra_nedelec_first_order( &
            center, values, curls_at_center, status)
        do basis = 1, 6
            numerical_curl = 0.0_dp
            do direction = 1, 3
                shifted = center
                shifted(direction) = shifted(direction) + step
                call evaluate_tetra_nedelec_first_order( &
                    shifted, values_plus, curls, status)
                shifted(direction) = center(direction) - step
                call evaluate_tetra_nedelec_first_order( &
                    shifted, values_minus, curls, status)
                select case (direction)
                case (1)
                    numerical_curl(2) = numerical_curl(2) - &
                        (values_plus(3, basis) - values_minus(3, basis)) / &
                        (2.0_dp * step)
                    numerical_curl(3) = numerical_curl(3) + &
                        (values_plus(2, basis) - values_minus(2, basis)) / &
                        (2.0_dp * step)
                case (2)
                    numerical_curl(1) = numerical_curl(1) + &
                        (values_plus(3, basis) - values_minus(3, basis)) / &
                        (2.0_dp * step)
                    numerical_curl(3) = numerical_curl(3) - &
                        (values_plus(1, basis) - values_minus(1, basis)) / &
                        (2.0_dp * step)
                case (3)
                    numerical_curl(1) = numerical_curl(1) - &
                        (values_plus(2, basis) - values_minus(2, basis)) / &
                        (2.0_dp * step)
                    numerical_curl(2) = numerical_curl(2) + &
                        (values_plus(1, basis) - values_minus(1, basis)) / &
                        (2.0_dp * step)
                end select
            end do
            call record_condition(maxval(abs( &
                numerical_curl - curls_at_center(:, basis))) < 2.0e-9_dp, &
                "Tetrahedral curls agree with independent differentiation")
        end do
    end subroutine check_curls_from_finite_differences

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_nedelec_first_order
