program test_level_set_tetra_cut_fourth_moments_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_level_set_tetra_cut_fourth_moments_3d, &
        evaluate_level_set_tetra_cut_fourth_moments_3d_jvp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: vertices(3, 4) = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 4])
    real(dp), parameter :: level_values(4) = [1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp]
    real(dp), parameter :: vertices_dot(3, 4) = reshape([ &
        0.02_dp, -0.01_dp, 0.03_dp, -0.02_dp, 0.04_dp, -0.01_dp, &
        0.01_dp, -0.03_dp, 0.02_dp, 0.03_dp, -0.02_dp, 0.01_dp], [3, 4])
    real(dp), parameter :: level_values_dot(4) = [0.06_dp, -0.02_dp, 0.04_dp, -0.03_dp]
    real(dp) :: positive_fourth(3, 3, 3, 3), negative_fourth(3, 3, 3, 3)
    real(dp) :: positive_fourth_dot(3, 3, 3, 3), negative_fourth_dot(3, 3, 3, 3)
    real(dp) :: plus_fourth(3, 3, 3, 3), minus_fourth(3, 3, 3, 3)
    real(dp) :: plus_negative_fourth(3, 3, 3, 3), minus_negative_fourth(3, 3, 3, 3)
    real(dp) :: expected_positive(3, 3, 3, 3), expected_parent(3, 3, 3, 3)
    real(dp) :: step
    integer :: status, plus_status, minus_status
    logical :: all_passed

    all_passed = .true.
    expected_positive = simplex_tensor(0.5_dp)
    expected_parent = simplex_tensor(1.0_dp)

    call evaluate_level_set_tetra_cut_fourth_moments_3d( &
        vertices, level_values, positive_fourth, negative_fourth, status)
    call record_condition(status == 0 .and. maxval(abs(positive_fourth - &
        expected_positive)) < 1.0e-14_dp, &
        "positive tetra cut integrates quartic monomials exactly")
    call record_condition(maxval(abs(positive_fourth + negative_fourth - &
        expected_parent)) < 1.0e-14_dp, &
        "positive and negative quartic moments conserve the parent tetrahedron")

    call evaluate_level_set_tetra_cut_fourth_moments_3d_jvp( &
        vertices, level_values, vertices_dot, level_values_dot, &
        positive_fourth_dot, negative_fourth_dot, status)
    step = 1.0e-6_dp
    call evaluate_level_set_tetra_cut_fourth_moments_3d( &
        vertices + step*vertices_dot, level_values + step*level_values_dot, &
        plus_fourth, plus_negative_fourth, plus_status)
    call evaluate_level_set_tetra_cut_fourth_moments_3d( &
        vertices - step*vertices_dot, level_values - step*level_values_dot, &
        minus_fourth, minus_negative_fourth, minus_status)
    call record_condition(status == 0 .and. plus_status == 0 .and. &
        minus_status == 0 .and. maxval(abs(positive_fourth_dot - &
        (plus_fourth - minus_fourth)/(2.0_dp*step))) < 1.0e-7_dp .and. &
        maxval(abs(negative_fourth_dot - (plus_negative_fourth - &
        minus_negative_fourth)/(2.0_dp*step))) < 1.0e-7_dp, &
        "tetrahedral quartic cut-moment JVP matches a fixed-topology difference oracle")

    call check_summary("3D level-set tetrahedral quartic cut moments")
    if (.not. all_passed) error stop 1

contains

    pure function simplex_tensor(scale) result(moment)
        real(dp), intent(in) :: scale
        real(dp) :: moment(3, 3, 3, 3)
        integer :: first, second, third, fourth
        integer :: x_degree, y_degree, z_degree

        do first = 1, 3
            do second = 1, 3
                do third = 1, 3
                    do fourth = 1, 3
                        x_degree = count([first, second, third, fourth] == 1)
                        y_degree = count([first, second, third, fourth] == 2)
                        z_degree = 4 - x_degree - y_degree
                        moment(first, second, third, fourth) = simplex_moment( &
                            x_degree, y_degree, z_degree, scale)
                    end do
                end do
            end do
        end do
    end function simplex_tensor

    pure real(dp) function simplex_moment(x_degree, y_degree, z_degree, scale) &
            result(value)
        integer, intent(in) :: x_degree, y_degree, z_degree
        real(dp), intent(in) :: scale

        value = scale**(x_degree + y_degree + z_degree + 3) * &
            factorial(x_degree)*factorial(y_degree)*factorial(z_degree)/ &
            factorial(x_degree + y_degree + z_degree + 3)
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

end program test_level_set_tetra_cut_fourth_moments_3d
