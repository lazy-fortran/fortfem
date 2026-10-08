program test_triangle_level_quadrature
    use fortfem_kinds, only: dp
    use check, only: check_condition, check_summary
    use fortfem_triangle_level_quadrature
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    type(triangle_level_workspace_t) :: workspace, limited
    real(dp) :: coefficients(6), expected(6), nodal(6), integral(6), error(6)
    real(dp) :: replay(6), values(6), absolute(6), levels(2), gauge
    real(dp), parameter :: nodes(2, 6) = reshape([0._dp, 0._dp, 1._dp, 0._dp, &
        0._dp, 1._dp, .5_dp, 0._dp, .5_dp, .5_dp, 0._dp, .5_dp], [2, 6])
    real(dp), parameter :: pi = acos(-1._dp), radius = .125_dp
    integer :: status, i, trial, mode
    character(180) :: message
    mode = 0
    call initialize_triangle_level_workspace(workspace, 6, 4, status)
    call check_condition(status == 0, 'Reusable vector workspace initializes')
    do trial = 1, 20
        expected = real([trial, 2 - trial, trial + 3, 2*trial, -trial, 1 - trial], dp)/8
        do i = 1, 6
            nodal(i) = polynomial(expected, nodes(:, i))
        end do
        call triangle_level_coefficients(2, nodal, coefficients, status)
        call check_condition(status == 0 .and. maxval(abs(coefficients - expected)) &
            < 1e-14_dp, 'P2 nodal reconstruction equals an independent quadratic')
        gauge = 2._dp**40
        call triangle_level_coefficients(2, nodal + gauge, coefficients, status)
        coefficients(1) = coefficients(1) - gauge
        call check_condition(status == 0 .and. maxval(abs(coefficients - expected)) &
            < 1e-14_dp, 'Dyadic P2 coefficients survive a large common gauge')
        expected(4:6) = 0
        do i = 1, 3
            nodal(i) = polynomial(expected, nodes(:, i))
        end do
        call triangle_level_coefficients(1, nodal(:3), coefficients, status)
        call check_condition(status == 0 .and. maxval(abs(coefficients - expected)) &
            < 1e-14_dp, 'P1 nodal reconstruction equals an independent affine scalar')
    end do
    call triangle_level_coefficients(2, nodal(:3), coefficients, status)
    call check_condition(status /= 0 .and. maxval(abs(coefficients)) <= 0, &
        'Nodal helper rejects an incompatible element size')
    coefficients = 0
    absolute = 1e-12_dp
    call integrate_triangle_levels(coefficients, levels(:0), callback, absolute, &
        1e-11_dp, workspace, integral, error, status, message)
    expected = [.5_dp, 1._dp/6, 1._dp/6, 1._dp/24, 1._dp/12, 0._dp]
    call check_condition(status == 0, 'Polynomial vector integration converges')
    call check_condition(maxval(abs(integral - expected)) < 2e-13_dp, &
        'Constant, affine, quadratic and cancelling moments have exact simplex values')
    call verify_trace(expected)
    mode = 1
    coefficients = [0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 1._dp]
    levels(1) = .25_dp
    call integrate_triangle_levels(coefficients, levels(:1), callback, absolute, &
        1e-11_dp, workspace, integral, error, status, message)
    expected = [1._dp/8, 5._dp/192, 1._dp/8, -1._dp/8, 0._dp, .5_dp]
    call check_condition(status == 0, 'Quadratic horizontal level integration converges')
    call check_condition(maxval(abs(integral - expected)) < 2e-12_dp, &
        'Eta squared step and continuous ramp match exact area and moment')
    call verify_trace(expected)
    mode = 2
    coefficients = [.25_dp**2*2 - radius**2, -.5_dp, -.5_dp, 1._dp, 0._dp, 1._dp]
    levels(1) = 0
    call integrate_triangle_levels(coefficients, levels(:1), callback, absolute, &
        1e-11_dp, workspace, integral, error, status, message)
    expected = [pi*radius**2, pi*radius**4/2, .25_dp*pi*radius**2, &
        0._dp, .5_dp, .5_dp]
    print *, 'circle status/work/panels/error: ', status, workspace%nevaluations, &
        workspace%npanels, maxval(abs(integral - expected)), trim(message)
    call check_condition(status == 0, 'Interior circular cut with tangency events converges')
    call check_condition(maxval(abs(integral - expected)) < 3e-12_dp, &
        'Interior disk step, ramp and first moment match independent polar integrals')
    if (status == 0) call verify_trace(expected)
    mode = 3
    coefficients = [0._dp, 1._dp, 0._dp, 0._dp, 0._dp, 0._dp]
    levels = [.25_dp, .75_dp]
    call integrate_triangle_levels(coefficients, levels, callback, absolute, &
        1e-11_dp, workspace, integral, error, status, message)
    expected = [9._dp/32, 1._dp/32, 9._dp/128, 1._dp/384, .5_dp, 0._dp]
    call check_condition(status == 0 .and. maxval(abs(integral - expected)) < 2e-12_dp, &
        'Two affine level cuts match independent triangular wedge moments')
    call verify_trace(expected)
    mode = 9
    coefficients = [0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 1._dp]
    levels(1) = .25_dp
    call integrate_triangle_levels(coefficients, levels(:1), callback, absolute, &
        1e-11_dp, workspace, integral, error, status, message)
    expected = [1._dp/8, -.25_dp, 0._dp, 1._dp/32, .5_dp, 0._dp]
    call check_condition(status == 0 .and. maxval(abs(integral - expected)) < 3e-12_dp, &
        'Continuous kink and discontinuous slope moments agree on the same quadratic cut')
    if (status == 0) call verify_trace(expected)
    mode = 6
    coefficients = 0
    call integrate_triangle_levels(coefficients, levels(:0), callback, absolute, &
        1e-11_dp, workspace, integral, error, status, message)
    expected = [2*log(2._dp) - 1, 1.5_dp - 2*log(2._dp), &
        2*log(2._dp) - 1.25_dp, 1._dp/306, 1.05_dp*log(21._dp) - 1, 0._dp]
    print *, 'rational status/work/depth/error: ', status, workspace%nevaluations, &
        workspace%inner_depth, maxval(abs(integral - expected)), trim(message)
    call check_condition(status == 0, 'Rational and degree sixteen vector integrals converge')
    call check_condition(workspace%inner_depth > 0, &
        'Inner error refinement is exercised independently of outer refinement')
    call check_condition(maxval(abs(integral - expected)) < 3e-12_dp, &
        'Rational logarithmic integrals and high-degree simplex moments match closed forms')
    if (status == 0) call verify_trace(expected)
    mode = 7
    coefficients = [0._dp, 0._dp, 0._dp, 1._dp, 0._dp, 1._dp]
    levels(1) = radius**2
    call integrate_triangle_levels(coefficients, levels(:1), callback, absolute, &
        1e-11_dp, workspace, integral, error, status, message)
    expected = [pi*radius**2/4, pi*radius**4/8, radius**3/3, &
        0._dp, .5_dp, .5_dp]
    call check_condition(status == 0 .and. maxval(abs(integral - expected)) < 3e-12_dp, &
        'A circular cut tangent to both coordinate edges matches polar quarter-disk moments')
    if (status == 0) call verify_trace(expected)
    mode = 8
    coefficients = [.25_dp**2*2, -.5_dp, -.5_dp, 1._dp, 0._dp, 1._dp]
    levels = [(.0625_dp)**2, radius**2]
    call integrate_triangle_levels(coefficients, levels, callback, absolute, &
        1e-11_dp, workspace, integral, error, status, message)
    expected = [pi*(radius**2 - .0625_dp**2), pi*.0625_dp**2, &
        pi*radius**2, 0._dp, .5_dp, .5_dp]
    call check_condition(status == 0 .and. maxval(abs(integral - expected)) < 3e-12_dp, &
        'Two quadratic knots partition an annulus without dropping either topology event')
    if (status == 0) call verify_trace(expected)
    mode = 0
    coefficients = [0._dp, 1._dp, 0._dp, 0._dp, 0._dp, 0._dp]
    levels = [.25_dp, .75_dp]
    call initialize_triangle_level_workspace(limited, 6, 4, status, max_panels=1)
    call integrate_triangle_levels(coefficients, levels, callback, absolute, &
        1e-11_dp, limited, integral, error, status, message)
    call check_condition(status == TRIANGLE_LEVEL_BUDGET .and. limited%npoints == 0, &
        'Initial topology rejects an insufficient panel budget explicitly')
    mode = 0
    call initialize_triangle_level_workspace(limited, 6, 4, status, max_evaluations=1)
    call integrate_triangle_levels(coefficients, levels, callback, absolute, &
        1e-11_dp, limited, integral, error, status, message)
    call check_condition(status == TRIANGLE_LEVEL_BUDGET .and. limited%npoints == 0 &
        .and. maxval(abs(integral)) <= 0, 'Evaluation exhaustion fails with no usable trace or load')
    call initialize_triangle_level_workspace(limited, 6, 4, status, max_points=1)
    call integrate_triangle_levels(coefficients, levels(:0), callback, absolute, &
        1e-11_dp, limited, integral, error, status, message)
    call check_condition(status == TRIANGLE_LEVEL_BUDGET .and. limited%npoints == 0 &
        .and. maxval(abs(integral)) <= 0, 'Trace capacity exhaustion fails without a partial trace')
    mode = 4
    call integrate_triangle_levels(coefficients, levels(:0), callback, absolute, &
        1e-11_dp, workspace, integral, error, status, message)
    call check_condition(status == TRIANGLE_LEVEL_CALLBACK .and. workspace%npoints == 0, &
        'A failing callback propagates explicitly')
    mode = 5
    call integrate_triangle_levels(coefficients, levels(:0), callback, absolute, &
        1e-11_dp, workspace, integral, error, status, message)
    call check_condition(status == TRIANGLE_LEVEL_CALLBACK .and. workspace%npoints == 0, &
        'Nonfinite callback values are rejected')
    mode = 0
    call integrate_triangle_levels(coefficients, levels(2:1:-1), callback, absolute, &
        1e-11_dp, workspace, integral, error, status, message)
    call check_condition(status == TRIANGLE_LEVEL_INVALID, 'Unsorted levels are rejected')
    call integrate_triangle_levels(coefficients, levels(:0), callback, -absolute, &
        1e-11_dp, workspace, integral, error, status, message)
    call check_condition(status == TRIANGLE_LEVEL_INVALID, 'Negative component tolerances are rejected')
    call integrate_triangle_levels(coefficients, levels(:0), callback, 0*absolute, &
        0._dp, workspace, integral, error, status, message)
    call check_condition(status == TRIANGLE_LEVEL_INVALID, 'An empty total error budget is rejected')
    coefficients(2) = ieee_value(0._dp, ieee_quiet_nan)
    call integrate_triangle_levels(coefficients, levels(:0), callback, absolute, &
        1e-11_dp, workspace, integral, error, status, message)
    call check_condition(status == TRIANGLE_LEVEL_INVALID, 'Nonfinite geometry is rejected before root evaluation')
    nodal(1) = ieee_value(0._dp, ieee_quiet_nan)
    call triangle_level_coefficients(2, nodal, coefficients, status)
    call check_condition(status == TRIANGLE_LEVEL_INVALID, 'Nonfinite nodal values are rejected before algebra')
    call initialize_triangle_level_workspace(limited, 6, huge(1), status)
    call check_condition(status == TRIANGLE_LEVEL_INVALID .and. .not. limited%initialized, &
        'Impossible topology capacities are rejected without integer overflow')
    call check_summary('Source-conforming quadratic-level triangle quadrature')
contains
    pure real(dp) function polynomial(c, point) result(value)
        real(dp), intent(in) :: c(6), point(2)
        value = c(1) + c(2)*point(1) + c(3)*point(2) + c(4)*point(1)**2 + &
            c(5)*point(1)*point(2) + c(6)*point(2)**2
    end function polynomial
    subroutine callback(point, values, status)
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: values(:)
        integer, intent(out) :: status
        real(dp) :: x, y, distance, indicator
        status = 0
        x = point(1)
        y = point(2)
        select case (mode)
        case (0)
            values = [1._dp, x, y, x*y, x*x, x - y]
        case (1)
            indicator = 0
            if (y*y > .25_dp) indicator = 1
            values = [indicator, max(y*y - .25_dp, 0._dp), indicator, &
                -indicator, 0._dp, 1._dp]
        case (2)
            distance = (x - .25_dp)**2 + (y - .25_dp)**2
            indicator = 0
            if (distance < radius**2) indicator = 1
            values = [indicator, max(radius**2 - distance, 0._dp), x*indicator, &
                (x - y)*indicator, 1._dp, 1._dp]
        case (3)
            values = [merge(1._dp, 0._dp, x > .25_dp), &
                merge(1._dp, 0._dp, x > .75_dp), max(x - .25_dp, 0._dp), &
                max(x - .75_dp, 0._dp), 1._dp, 0._dp]
        case (9)
            indicator = -1
            if (y > .5_dp) indicator = 1
            values = [abs(y - .5_dp), indicator, y*indicator, y*y*indicator, 1._dp, 0._dp]
        case (6)
            values = [1/(1 + x), x/(1 + x), y/(1 + x), y**16, &
                1/(.05_dp + y), (x - y)/(1 + x + y)]
        case (7)
            distance = x*x + y*y
            indicator = 0
            if (distance < radius**2) indicator = 1
            values = [indicator, max(radius**2 - distance, 0._dp), x*indicator, &
                (x - y)*indicator, 1._dp, 1._dp]
        case (8)
            distance = (x - .25_dp)**2 + (y - .25_dp)**2
            values = [merge(1._dp, 0._dp, distance < radius**2 &
                .and. distance > .0625_dp**2), &
                merge(1._dp, 0._dp, distance < .0625_dp**2), &
                merge(1._dp, 0._dp, distance < radius**2), &
                0._dp, 1._dp, 1._dp]
        case (4)
            status = 42
            values = 0
        case (5)
            values = ieee_value(0._dp, ieee_quiet_nan)
        end select
    end subroutine callback
    subroutine verify_trace(expected)
        real(dp), intent(in) :: expected(6)
        integer :: index, callback_status
        call check_condition(workspace%npoints > 0, 'Accepted integration retains reference nodes')
        call check_condition(all(workspace%weights(:workspace%npoints) > 0), &
            'Every retained reference weight is positive')
        call check_condition(abs(sum(workspace%weights(:workspace%npoints)) - .5_dp) &
            < 2e-13_dp, 'Retained weights cover the full triangle area')
        replay = 0
        do index = 1, workspace%npoints
            call callback(workspace%points(:, index), values, callback_status)
            replay = replay + workspace%weights(index)*values
        end do
        call check_condition(maxval(abs(replay - expected)) < 3e-12_dp, &
            'Independent replay on accepted positive trace preserves the integral')
    end subroutine verify_trace
end program test_triangle_level_quadrature
