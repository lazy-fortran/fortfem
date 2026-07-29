program test_generated_tetra_nedelec_candidates
    use check, only: check_condition, check_summary
    use fortfem_generated_tetra_nedelec_candidates_order_1, only: &
        evaluate_candidates_order_1
    use fortfem_generated_tetra_nedelec_candidates_order_2, only: &
        evaluate_candidates_order_2
    use fortfem_generated_tetra_nedelec_candidates_order_3, only: &
        evaluate_candidates_order_3
    use fortfem_generated_tetra_nedelec_candidates_order_4, only: &
        evaluate_candidates_order_4
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: point(3) = [0.23_dp, 0.19_dp, 0.17_dp]
    real(dp), allocatable :: curls(:, :), expected(:, :), values(:, :)
    integer :: dof_count, order
    logical :: all_passed

    all_passed = .true.
    do order = 1, 4
        dof_count = order * (order + 2) * (order + 3) / 2
        allocate( &
            values(3, dof_count), curls(3, dof_count), &
            expected(3, dof_count))
        call evaluate_generated(order, point, values, curls)
        call expected_candidate_values(order, point, expected)
        call record_condition(maxval(abs(values - expected)) < 2.0e-14_dp, &
            "Generated tetrahedral candidates match their polynomial definition")
        call check_curls_by_finite_difference(order, dof_count, curls)
        deallocate(values, curls, expected)
    end do

    call check_summary("Generated tetrahedral Nedelec candidates")
    if (.not. all_passed) error stop 1

contains

    subroutine evaluate_generated(order, location, values, curls)
        integer, intent(in) :: order
        real(dp), intent(in) :: location(3)
        real(dp), intent(out) :: values(:, :), curls(:, :)

        select case (order)
        case (1)
            call evaluate_candidates_order_1( &
                location(1), location(2), location(3), values, curls)
        case (2)
            call evaluate_candidates_order_2( &
                location(1), location(2), location(3), values, curls)
        case (3)
            call evaluate_candidates_order_3( &
                location(1), location(2), location(3), values, curls)
        case (4)
            call evaluate_candidates_order_4( &
                location(1), location(2), location(3), values, curls)
        end select
    end subroutine evaluate_generated

    subroutine expected_candidate_values(order, location, values)
        integer, intent(in) :: order
        real(dp), intent(in) :: location(3)
        real(dp), intent(out) :: values(:, :)

        real(dp) :: monomial
        integer :: candidate, component, total_degree
        integer :: x_degree, y_degree, z_degree

        values = 0.0_dp
        candidate = 0
        do component = 1, 3
            do total_degree = 0, order - 1
                do x_degree = 0, total_degree
                    do y_degree = 0, total_degree - x_degree
                        z_degree = total_degree - x_degree - y_degree
                        candidate = candidate + 1
                        values(component, candidate) = &
                            monomial_value( &
                            location, [x_degree, y_degree, z_degree])
                    end do
                end do
            end do
        end do
        total_degree = order - 1
        do component = 4, 5
            do x_degree = 0, total_degree
                do y_degree = 0, total_degree - x_degree
                    z_degree = total_degree - x_degree - y_degree
                    candidate = candidate + 1
                    monomial = monomial_value( &
                        location, [x_degree, y_degree, z_degree])
                    if (component == 4) then
                        values(:, candidate) = [ &
                            -location(2) * monomial, &
                            location(1) * monomial, 0.0_dp]
                    else
                        values(:, candidate) = [ &
                            -location(3) * monomial, 0.0_dp, &
                            location(1) * monomial]
                    end if
                end do
            end do
        end do
        do y_degree = 0, total_degree
            z_degree = total_degree - y_degree
            candidate = candidate + 1
            monomial = monomial_value( &
                location, [0, y_degree, z_degree])
            values(:, candidate) = [ &
                0.0_dp, -location(3) * monomial, &
                location(2) * monomial]
        end do
    end subroutine expected_candidate_values

    subroutine check_curls_by_finite_difference(order, dof_count, curls)
        integer, intent(in) :: order, dof_count
        real(dp), intent(in) :: curls(:, :)

        real(dp), parameter :: step = 2.0e-6_dp
        real(dp), allocatable :: scratch_curls(:, :)
        real(dp), allocatable :: values_minus(:, :), values_plus(:, :)
        real(dp) :: numerical(3), shifted(3)
        integer :: candidate, direction

        allocate( &
            scratch_curls(3, dof_count), values_minus(3, dof_count), &
            values_plus(3, dof_count))
        do candidate = 1, dof_count
            numerical = 0.0_dp
            do direction = 1, 3
                shifted = point
                shifted(direction) = shifted(direction) + step
                call evaluate_generated( &
                    order, shifted, values_plus, scratch_curls)
                shifted(direction) = shifted(direction) - 2.0_dp * step
                call evaluate_generated( &
                    order, shifted, values_minus, scratch_curls)
                select case (direction)
                case (1)
                    numerical(2) = numerical(2) - &
                        centered_difference( &
                        values_plus(3, candidate), &
                        values_minus(3, candidate), step)
                    numerical(3) = numerical(3) + &
                        centered_difference( &
                        values_plus(2, candidate), &
                        values_minus(2, candidate), step)
                case (2)
                    numerical(1) = numerical(1) + &
                        centered_difference( &
                        values_plus(3, candidate), &
                        values_minus(3, candidate), step)
                    numerical(3) = numerical(3) - &
                        centered_difference( &
                        values_plus(1, candidate), &
                        values_minus(1, candidate), step)
                case (3)
                    numerical(1) = numerical(1) - &
                        centered_difference( &
                        values_plus(2, candidate), &
                        values_minus(2, candidate), step)
                    numerical(2) = numerical(2) + &
                        centered_difference( &
                        values_plus(1, candidate), &
                        values_minus(1, candidate), step)
                end select
            end do
            call record_condition(maxval(abs( &
                numerical - curls(:, candidate))) < 3.0e-10_dp, &
                "Generated tetrahedral candidate curl matches finite differences")
        end do
    end subroutine check_curls_by_finite_difference

    pure function centered_difference(plus, minus, step) result(value)
        real(dp), intent(in) :: plus, minus, step
        real(dp) :: value

        value = (plus - minus) / (2.0_dp * step)
    end function centered_difference

    pure function monomial_value(location, powers) result(value)
        real(dp), intent(in) :: location(3)
        integer, intent(in) :: powers(3)
        real(dp) :: value

        value = location(1)**powers(1) * location(2)**powers(2) * &
            location(3)**powers(3)
    end function monomial_value

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_generated_tetra_nedelec_candidates
