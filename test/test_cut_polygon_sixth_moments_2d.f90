program test_cut_polygon_sixth_moments_2d
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        evaluate_cut_polygon_sixth_moments_2d, &
        evaluate_cut_polygon_sixth_moments_2d_jvp, &
        evaluate_cut_polygon_sixth_moments_2d_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp) :: polygon(2, 4), polygon_dot(2, 4), polygon_plus(2, 4), polygon_minus(2, 4)
    real(dp) :: moment(2, 2, 2, 2, 2, 2), moment_dot(2, 2, 2, 2, 2, 2)
    real(dp) :: moment_plus(2, 2, 2, 2, 2, 2), moment_minus(2, 2, 2, 2, 2, 2)
    real(dp) :: moment_bar(2, 2, 2, 2, 2, 2), polygon_bar(2, 4), lhs, rhs
    real(dp) :: flat(64), flat_bar(64), expected, conserved, conserved_moments, step
    integer :: index, bit, x_degree, y_degree, work
    type(fortsparse_status_t) :: status

    step = 1.0e-7_dp
    polygon = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, &
        0.0_dp, 1.0_dp], [2, 4])
    polygon_dot = reshape([0.02_dp, -0.01_dp, -0.03_dp, 0.04_dp, 0.05_dp, 0.01_dp, &
        -0.02_dp, 0.03_dp], [2, 4])

    call evaluate_cut_polygon_sixth_moments_2d(polygon, moment, status)
    call check_condition(status%code == 0, &
        "degree-six polygon moment accepts a nondegenerate cut cell")
    flat = reshape(moment, [64])
    do index = 0, 63
        work = index
        x_degree = 0
        do bit = 1, 6
            if (mod(work, 2) == 1) x_degree = x_degree + 1
            work = work/2
        end do
        y_degree = 6 - x_degree
        expected = 1.0_dp/(real(x_degree + 1, dp)*real(y_degree + 1, dp))
        call check_condition(abs(flat(index + 1) - expected) < 1.0e-13_dp, &
            "unit-square degree-six polynomial Gauss-Green oracle")
    end do

    conserved = 0.0_dp
    conserved_moments = 0.0_dp
    do index = 0, 6
        conserved = conserved + binomial(6, index)/(real(index + 1, dp)*real(7 - index, dp))
        conserved_moments = conserved_moments + binomial(6, index)*flat(2**index)
    end do
    call check_condition(abs(conserved - conserved_moments) &
        < 1.0e-13_dp, "degree-six moments preserve polynomial conservation")
    call check_condition(abs(flat(1) - flat(64)) < 1.0e-13_dp, &
        "degree-six moments preserve tensor symmetry")

    call evaluate_cut_polygon_sixth_moments_2d_jvp(polygon, polygon_dot, moment_dot, status)
    polygon_plus = polygon + step*polygon_dot
    polygon_minus = polygon - step*polygon_dot
    call evaluate_cut_polygon_sixth_moments_2d(polygon_plus, moment_plus, status)
    call evaluate_cut_polygon_sixth_moments_2d(polygon_minus, moment_minus, status)
    call check_condition(maxval(abs(moment_dot - (moment_plus - moment_minus)/(2.0_dp*step))) < 1.0e-8_dp, &
        "degree-six moment JVP matches central differences")

    moment_bar = reshape([(0.01_dp*real(index, dp), index=1, 64)], [2, 2, 2, 2, 2, 2])
    call evaluate_cut_polygon_sixth_moments_2d_vjp(polygon, moment_bar, polygon_bar, status)
    flat_bar = reshape(moment_bar, [64])
    lhs = sum(flat_bar*reshape(moment_dot, [64]))
    rhs = sum(polygon_bar*polygon_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "degree-six moment VJP satisfies the real transpose oracle")

    call evaluate_cut_polygon_sixth_moments_2d(polygon(:, [1, 2, 2, 1]), moment, status)
    call check_condition(status%code /= 0, "degree-six moments reject degenerate topology")
    call check_summary("cut polygon sixth moments")

contains

    pure function binomial(number, choose) result(value)
        integer, intent(in) :: number, choose
        real(dp) :: value
        integer :: factor

        value = 0.0_dp
        if (choose < 0 .or. choose > number) return
        value = 1.0_dp
        do factor = 1, choose
            value = value*real(number - choose + factor, dp)/real(factor, dp)
        end do
    end function binomial

end program test_cut_polygon_sixth_moments_2d
