program test_triangle_level_conditioning
    use fortfem_kinds, only: dp
    use fortfem_triangle_level_quadrature
    use check, only: check_condition, check_summary
    implicit none
    type(triangle_level_workspace_t) :: workspace
    real(dp) :: coefficients(6), levels(1), integral(4), error(4), exact(4), replay(4)
    real(dp) :: point_values(4), absolute(4), scale, gauge, radius, center(2), tolerance(4)
    real(dp), parameter :: pi = acos(-1._dp)
    real(dp), parameter :: scales(5) = [1e-20_dp, 1e-8_dp, 1._dp, 1e8_dp, 1e20_dp]
    integer :: status, trial, point_id, mode, evaluations(7), panels(7), points(7)
    character(180) :: message
    call initialize_triangle_level_workspace(workspace, 4, 1, status)
    absolute = 1e-13_dp
    mode = 1
    do trial = 1, size(scales)
        scale = scales(trial)
        coefficients = [0._dp, 0._dp, 0._dp, 0._dp, 0._dp, scale]
        levels = .25_dp*scale
        call integrate_triangle_levels(coefficients, levels, callback, absolute, &
            1e-11_dp, workspace, integral, error, status, message)
        exact = [1._dp/8, 5._dp/192, 1._dp/12, .5_dp]
        call qualify('Cut topology is invariant under twenty-decade coefficient scaling')
    end do
    gauge = 2._dp**40
    coefficients = [gauge, 0._dp, 0._dp, 0._dp, 0._dp, 1._dp]
    levels = gauge + .25_dp
    call integrate_triangle_levels(coefficients, levels, callback, absolute, &
        1e-11_dp, workspace, integral, error, status, message)
    call qualify('A representable common large gauge preserves cut moments')
    mode = 2
    center = [.25_dp, .25_dp]
    do trial = 1, 7
        radius = 2._dp**(-trial - 2)
        coefficients = [sum(center**2) - radius**2, -2*center(1), -2*center(2), &
            1._dp, 0._dp, 1._dp]
        levels = 0
        call integrate_triangle_levels(coefficients, levels, callback, absolute, &
            1e-11_dp, workspace, integral, error, status, message)
        exact = [pi*radius**2, pi*radius**4/2, center(2)*pi*radius**2, .5_dp]
        call qualify('Shrinking closed circular cuts retain both tangency events')
        evaluations(trial) = workspace%nevaluations
        panels(trial) = workspace%npanels
        points(trial) = workspace%npoints
        call check_condition(workspace%nevaluations <= 10000, &
            'Generated tangent regularization bounds work on circular exact oracles')
    end do
    print *, 'circle conditioning callback counts: ', evaluations
    print *, 'circle conditioning panel counts: ', panels
    print *, 'circle conditioning retained points: ', points
    mode = 3
    coefficients = 0
    levels = 0
    absolute = huge(1._dp)/2
    call integrate_triangle_levels(coefficients, levels(:0), callback, absolute, &
        huge(1._dp)/2, workspace, integral, error, status, message)
    call check_condition(status == 0 .and. all(integral > 4.9e299_dp) &
        .and. all(integral < 5.1e299_dp), &
        'Finite large callback and loose budgets do not overflow tolerance arithmetic')
    call check_summary('Quadratic-level triangle conditioning and empirical error estimates')
contains
    subroutine callback(point, values, status)
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: values(:)
        integer, intent(out) :: status
        real(dp) :: indicator, distance
        status = 0
        indicator = 0
        if (mode == 1) then
            if (point(2) > .5_dp) indicator = 1
            values = [indicator, max(point(2)**2 - .25_dp, 0._dp), &
                point(2)*indicator, 1._dp]
        else if (mode == 3) then
            values = 1e300_dp
        else
            distance = sum((point - center)**2)
            if (distance < radius**2) indicator = 1
            values = [indicator, max(radius**2 - distance, 0._dp), &
                point(2)*indicator, 1._dp]
        end if
    end subroutine callback
    subroutine qualify(label)
        character(*), intent(in) :: label
        call check_condition(status == 0, label)
        if (status /= 0) then
            print *, trim(message)
            return
        end if
        tolerance = absolute + 1e-11_dp*abs(integral)
        call check_condition(all(error <= tolerance), &
            'Returned component estimates meet the requested absolute-plus-relative budget')
        call check_condition(all(abs(integral - exact) <= tolerance), &
            'Independent analytical integrals meet the same requested budget')
        ! This empirical check allows explicitly stated floating-point accumulation
        ! error. Embedded-rule estimates are not a rigorous enclosure for all callbacks.
        call check_condition(all(abs(integral - exact) <= &
            error + 64*epsilon(1._dp)*max(1._dp, abs(exact))), &
            'Analytical error is consistent with the estimate and floating-point floor')
        call check_condition(all(workspace%weights(:workspace%npoints) > 0), &
            'Conditioned cuts retain positive weights')
        replay = 0
        do point_id = 1, workspace%npoints
            call callback(workspace%points(:, point_id), point_values, status)
            replay = replay + workspace%weights(point_id)*point_values
        end do
        call check_condition(maxval(abs(replay - integral)) < 1e-13_dp, &
            'Reused accepted trace agrees with the qualified vector integral')
    end subroutine qualify
end program test_triangle_level_conditioning
