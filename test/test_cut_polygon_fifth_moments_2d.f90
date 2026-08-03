program test_cut_polygon_fifth_moments_2d
    use check, only: check_condition, check_summary
    use fortfem_cut_polygon_fifth_moments_2d, only: &
        evaluate_cut_polygon_fifth_moments_2d, &
        evaluate_cut_polygon_fifth_moments_2d_jvp, &
        evaluate_cut_polygon_fifth_moments_2d_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp) :: polygon(2, 4), polygon_dot(2, 4), polygon_plus(2, 4), polygon_minus(2, 4)
    real(dp) :: moment(2, 2, 2, 2, 2), moment_dot(2, 2, 2, 2, 2)
    real(dp) :: moment_plus(2, 2, 2, 2, 2), moment_minus(2, 2, 2, 2, 2)
    real(dp) :: moment_bar(2, 2, 2, 2, 2), polygon_bar(2, 4), lhs, rhs
    real(dp) :: expected, step
    integer :: first, second, third, fourth, fifth, x_degree
    type(fortsparse_status_t) :: status

    step = 1.0e-7_dp
    polygon = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], [2, 4])
    polygon_dot = reshape([0.02_dp, -0.01_dp, -0.03_dp, 0.04_dp, 0.05_dp, 0.01_dp, &
        -0.02_dp, 0.03_dp], [2, 4])
    call evaluate_cut_polygon_fifth_moments_2d(polygon, moment, status)
    call check_condition(status%code == 0, "degree-five polygon moment accepts a nondegenerate cut cell")
    do first = 1, 2
        do second = 1, 2
            do third = 1, 2
                do fourth = 1, 2
                    do fifth = 1, 2
                        x_degree = count([first, second, third, fourth, fifth] == 1)
                        expected = 1.0_dp/(real(x_degree + 1, dp)*real(6 - x_degree, dp))
                        if (abs(moment(first, second, third, fourth, fifth) - expected) > 1.0e-13_dp) then
                            call check_condition(.false., "unit-square degree-five polynomial oracle")
                        end if
                    end do
                end do
            end do
        end do
    end do
    call check_condition(maxval(abs(moment(1, 1, 1, 1, :) - [1.0_dp/6.0_dp, 1.0_dp/10.0_dp])) < 1.0e-13_dp, &
        "unit-square fifth moments preserve tensor symmetry")

    call evaluate_cut_polygon_fifth_moments_2d_jvp(polygon, polygon_dot, moment_dot, status)
    polygon_plus = polygon + step*polygon_dot
    polygon_minus = polygon - step*polygon_dot
    call evaluate_cut_polygon_fifth_moments_2d(polygon_plus, moment_plus, status)
    call evaluate_cut_polygon_fifth_moments_2d(polygon_minus, moment_minus, status)
    call check_condition(maxval(abs(moment_dot - (moment_plus - moment_minus)/(2.0_dp*step))) < 1.0e-8_dp, &
        "degree-five moment JVP matches central differences")

    moment_bar = reshape([(0.01_dp*real(first, dp), first=1, 32)], [2, 2, 2, 2, 2])
    call evaluate_cut_polygon_fifth_moments_2d_vjp(polygon, moment_bar, polygon_bar, status)
    lhs = sum(moment_bar*moment_dot)
    rhs = sum(polygon_bar*polygon_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "degree-five moment VJP satisfies the real transpose oracle")
    call evaluate_cut_polygon_fifth_moments_2d(polygon(:, [1, 2, 2, 1]), moment, status)
    call check_condition(status%code /= 0, "degree-five moments reject degenerate topology")
    call check_summary("cut polygon fifth moments")
end program test_cut_polygon_fifth_moments_2d
