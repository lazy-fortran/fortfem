program test_level_set_triangle_cut_fourth_moments_2d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_level_set_triangle_cut_fourth_moments_2d, &
        evaluate_level_set_triangle_cut_fourth_moments_2d_jvp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: vertices(2, 3) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
    real(dp), parameter :: level_values(3) = [1.0_dp, -1.0_dp, -1.0_dp]
    real(dp), parameter :: vertices_dot(2, 3) = reshape([ &
        0.03_dp, -0.02_dp, -0.01_dp, 0.04_dp, 0.02_dp, -0.01_dp], [2, 3])
    real(dp), parameter :: level_values_dot(3) = [0.07_dp, -0.03_dp, 0.05_dp]
    real(dp) :: positive_fourth(2, 2, 2, 2), negative_fourth(2, 2, 2, 2)
    real(dp) :: positive_fourth_dot(2, 2, 2, 2), negative_fourth_dot(2, 2, 2, 2)
    real(dp) :: plus_fourth(2, 2, 2, 2), minus_fourth(2, 2, 2, 2)
    real(dp) :: plus_negative_fourth(2, 2, 2, 2), minus_negative_fourth(2, 2, 2, 2)
    real(dp) :: fd_positive_dot(2, 2, 2, 2), fd_negative_dot(2, 2, 2, 2)
    real(dp) :: expected_positive(2, 2, 2, 2), expected_parent(2, 2, 2, 2)
    real(dp) :: step
    integer :: status, plus_status, minus_status
    logical :: all_passed

    all_passed = .true.
    expected_positive = 0.0_dp
    expected_parent = 0.0_dp
    call set_simplex_moments(expected_positive, 0.5_dp)
    call set_simplex_moments(expected_parent, 1.0_dp)

    call evaluate_level_set_triangle_cut_fourth_moments_2d( &
        vertices, level_values, positive_fourth, negative_fourth, status)
    call record_condition(status == 0 .and. maxval(abs(positive_fourth - &
        expected_positive)) < 1.0e-14_dp, &
        "positive cut integrates quartic monomials exactly")
    call record_condition(maxval(abs(positive_fourth + negative_fourth - &
        expected_parent)) < 1.0e-14_dp, &
        "positive and negative quartic moments conserve the parent triangle")

    call evaluate_level_set_triangle_cut_fourth_moments_2d_jvp( &
        vertices, level_values, vertices_dot, level_values_dot, &
        positive_fourth_dot, negative_fourth_dot, status)
    step = 1.0e-6_dp
    call evaluate_level_set_triangle_cut_fourth_moments_2d( &
        vertices + step*vertices_dot, level_values + step*level_values_dot, &
        plus_fourth, plus_negative_fourth, plus_status)
    call evaluate_level_set_triangle_cut_fourth_moments_2d( &
        vertices - step*vertices_dot, level_values - step*level_values_dot, &
        minus_fourth, minus_negative_fourth, minus_status)
    fd_positive_dot = (plus_fourth - minus_fourth)/(2.0_dp*step)
    fd_negative_dot = (plus_negative_fourth - minus_negative_fourth)/(2.0_dp*step)
    call record_condition(status == 0 .and. plus_status == 0 .and. &
        minus_status == 0 .and. maxval(abs(positive_fourth_dot - &
        fd_positive_dot)) < 1.0e-7_dp .and. maxval(abs(negative_fourth_dot - &
        fd_negative_dot)) < 1.0e-7_dp, &
        "quartic cut-moment JVP matches a fixed-topology difference oracle")

    call check_summary("2D level-set quartic cut moments")
    if (.not. all_passed) error stop 1

contains

    subroutine set_simplex_moments(moment, scale)
        real(dp), intent(out) :: moment(2, 2, 2, 2)
        real(dp), intent(in) :: scale
        integer :: first, second, third, fourth, x_degree

        moment = 0.0_dp
        do first = 1, 2
            do second = 1, 2
                do third = 1, 2
                    do fourth = 1, 2
                        x_degree = count([first, second, third, fourth] == 1)
                        moment(first, second, third, fourth) = simplex_moment( &
                            x_degree, 4 - x_degree, scale)
                    end do
                end do
            end do
        end do
    end subroutine set_simplex_moments

    pure real(dp) function simplex_moment(x_degree, y_degree, scale) result(value)
        integer, intent(in) :: x_degree, y_degree
        real(dp), intent(in) :: scale

        value = scale**(x_degree + y_degree + 2) * &
            factorial(x_degree)*factorial(y_degree)/ &
            factorial(x_degree + y_degree + 2)
    end function simplex_moment

    pure real(dp) function factorial(number) result(value)
        integer, intent(in) :: number
        integer :: factor

        value = 1.0_dp
        do factor = 2, number
            value = value*real(factor, dp)
        end do
    end function factorial

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_level_set_triangle_cut_fourth_moments_2d
