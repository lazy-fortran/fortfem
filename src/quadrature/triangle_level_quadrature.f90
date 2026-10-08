module fortfem_triangle_level_quadrature
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_roots, only: root_brent
    use fortnum_status, only: fortnum_status_t
    use fortfem_generated_triangle_level_events, only: generated_triangle_level_events
    use fortfem_generated_triangle_level_slice, only: generated_triangle_level_slice
    use fortfem_generated_triangle_level_polynomial, only: &
        generated_triangle_level_polynomial
    use fortfem_generated_triangle_level_tangent_interval, only: &
        generated_triangle_level_tangent_interval
    use fortfem_generated_triangle_level_interval, only: generated_triangle_level_interval
    use fortfem_generated_reference_p1_scalar_coefficients, only: &
        generated_reference_p1_scalar_coefficients
    use fortfem_generated_reference_p2_scalar_coefficients, only: &
        generated_reference_p2_scalar_coefficients
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: triangle_level_workspace_t, triangle_level_integrand_t
    public :: initialize_triangle_level_workspace, triangle_level_coefficients
    public :: integrate_triangle_levels
    integer, parameter, public :: TRIANGLE_LEVEL_INVALID = 1
    integer, parameter, public :: TRIANGLE_LEVEL_CALLBACK = 2
    integer, parameter, public :: TRIANGLE_LEVEL_BUDGET = 3

    abstract interface
        subroutine triangle_level_integrand_t(point, values, status)
            import :: dp
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: values(:)
            integer, intent(out) :: status
        end subroutine triangle_level_integrand_t
    end interface

    type :: triangle_level_workspace_t
        integer :: nvalues = 0, max_levels = 0, max_panels = 0
        integer :: max_evaluations = 0, max_points = 0
        integer :: npoints = 0, nevaluations = 0, npanels = 0, inner_depth = 0
        logical :: initialized = .false.
        ! Successful calls expose a positive weighted reference trace in 1:npoints.
        real(dp), allocatable :: points(:, :), weights(:)
        real(dp), allocatable :: panel_lower(:), panel_upper(:)
        real(dp), allocatable :: panel_value(:, :), panel_error(:, :), panel_inner_error(:, :)
        real(dp), allocatable :: cuts(:, :), events(:), inner_events(:), tangencies(:)
        logical, allocatable :: panel_tangent(:)
        real(dp), allocatable :: scratch(:, :)
        real(dp) :: nodes4(4), weights4(4), nodes8(8), weights8(8)
        real(dp) :: nodes16(16), weights16(16), nodes32(32), weights32(32)
    end type triangle_level_workspace_t
contains
    subroutine initialize_triangle_level_workspace(workspace, nvalues, max_levels, &
            status, max_panels, max_evaluations, max_points)
        type(triangle_level_workspace_t), intent(out) :: workspace
        integer, intent(in) :: nvalues, max_levels
        integer, intent(out) :: status
        integer, intent(in), optional :: max_panels, max_evaluations, max_points
        integer :: panels, evaluations, points, ierr
        status = TRIANGLE_LEVEL_INVALID
        if (nvalues < 1 .or. max_levels < 0) return
        if (max_levels > (huge(1) - 2)/6) return
        panels = 256
        evaluations = 200000
        points = 65536
        if (present(max_panels)) panels = max_panels
        if (present(max_evaluations)) evaluations = max_evaluations
        if (present(max_points)) points = max_points
        if (panels < 1 .or. evaluations < 1 .or. points < 1) return
        allocate(workspace%points(2, points), workspace%weights(points), &
            workspace%panel_lower(panels), workspace%panel_upper(panels), &
            workspace%panel_value(nvalues, panels), &
            workspace%panel_error(nvalues, panels), &
            workspace%panel_inner_error(nvalues, panels), workspace%cuts(6, max_levels), &
            workspace%events(6*max_levels + 2), workspace%tangencies(2*max_levels), &
            workspace%panel_tangent(panels), &
            workspace%inner_events(2*max_levels + 2), &
            workspace%scratch(nvalues, 8), stat=ierr)
        if (ierr /= 0) return
        call gauss_legendre_ab(4, 0._dp, 1._dp, workspace%nodes4, workspace%weights4)
        call gauss_legendre_ab(8, 0._dp, 1._dp, workspace%nodes8, workspace%weights8)
        call gauss_legendre_ab(16, 0._dp, 1._dp, workspace%nodes16, workspace%weights16)
        call gauss_legendre_ab(32, 0._dp, 1._dp, workspace%nodes32, workspace%weights32)
        workspace%nvalues = nvalues
        workspace%max_levels = max_levels
        workspace%max_panels = panels
        workspace%max_evaluations = evaluations
        workspace%max_points = points
        workspace%initialized = .true.
        status = 0
    end subroutine initialize_triangle_level_workspace

    subroutine triangle_level_coefficients(order, nodal_values, coefficients, status)
        integer, intent(in) :: order
        real(dp), intent(in) :: nodal_values(:)
        real(dp), intent(out) :: coefficients(6)
        integer, intent(out) :: status
        real(dp) :: centered(6)
        integer :: count
        coefficients = 0
        status = TRIANGLE_LEVEL_INVALID
        select case (order)
        case (1)
            count = 3
        case (2)
            count = 6
        case default
            return
        end select
        if (size(nodal_values) /= count) return
        if (.not. all(ieee_is_finite(nodal_values))) return
        ! Reject unrepresentable coefficient ranges before centered subtraction.
        if (maxval(abs(nodal_values)) > huge(1._dp)/64) return
        centered = 0
        centered(:count) = nodal_values - nodal_values(1)
        if (order == 1) then
            call generated_reference_p1_scalar_coefficients(centered(:3), coefficients)
        else
            call generated_reference_p2_scalar_coefficients(centered, coefficients)
        end if
        coefficients(1) = nodal_values(1)
        status = 0
    end subroutine triangle_level_coefficients

    subroutine interval_roots(coefficients, lower, upper, roots, count, status)
        real(dp), intent(in) :: coefficients(3), lower, upper
        real(dp), intent(out) :: roots(2)
        integer, intent(out) :: count, status
        type(fortnum_status_t) :: root_status
        real(dp) :: c(3), knots(3), values(3), slopes(2), scale, turning, root
        integer :: nknots, i
        roots = 0
        count = 0
        status = 0
        scale = maxval(abs(coefficients))
        if (scale <= 0) return
        c = coefficients/scale
        call generated_triangle_level_polynomial(c, lower, values(1), slopes(1))
        call generated_triangle_level_polynomial(c, upper, values(2), slopes(2))
        nknots = 2
        knots(:2) = [lower, upper]
        if ((slopes(1) < 0 .and. slopes(2) > 0) .or. &
            (slopes(1) > 0 .and. slopes(2) < 0)) then
            call root_brent(slope_at, lower, upper, turning, root_status, &
                xtol=4*epsilon(1._dp), max_iter=80)
            if (root_status%code /= 0) then
                status = TRIANGLE_LEVEL_BUDGET
                return
            end if
            nknots = 3
            knots = [lower, turning, upper]
            call generated_triangle_level_polynomial(c, turning, values(2), root)
            call generated_triangle_level_polynomial(c, upper, values(3), root)
        end if
        do i = 1, nknots
            if (abs(values(i)) <= 0) call append_root(knots(i))
        end do
        do i = 1, nknots - 1
            if (.not. ((values(i) < 0 .and. values(i + 1) > 0) .or. &
                (values(i) > 0 .and. values(i + 1) < 0))) cycle
            call root_brent(value_at, knots(i), knots(i + 1), root, root_status, &
                xtol=4*epsilon(1._dp), max_iter=80)
            if (root_status%code /= 0) then
                status = TRIANGLE_LEVEL_BUDGET
                return
            end if
            call append_root(root)
        end do
    contains
        pure real(dp) function value_at(x) result(value)
            real(dp), intent(in) :: x
            real(dp) :: slope
            call generated_triangle_level_polynomial(c, x, value, slope)
        end function value_at
        pure real(dp) function slope_at(x) result(slope)
            real(dp), intent(in) :: x
            real(dp) :: value
            call generated_triangle_level_polynomial(c, x, value, slope)
        end function slope_at
        subroutine append_root(x)
            real(dp), intent(in) :: x
            if (count > 0) then
                if (any(abs(roots(:count) - x) <= 0)) return
            end if
            if (count == 2) return
            count = count + 1
            roots(count) = x
        end subroutine append_root
    end subroutine interval_roots

    subroutine sort_unique(points, count)
        real(dp), intent(inout) :: points(:)
        integer, intent(inout) :: count
        real(dp) :: value
        integer :: i, j, retained
        do i = 2, count
            value = points(i)
            j = i - 1
            do while (j >= 1)
                if (points(j) <= value) exit
                points(j + 1) = points(j)
                j = j - 1
            end do
            points(j + 1) = value
        end do
        retained = 0
        do i = 1, count
            if (retained > 0) then
                if (abs(points(i) - points(retained)) <= 0) cycle
            end if
            retained = retained + 1
            points(retained) = points(i)
        end do
        count = retained
    end subroutine sort_unique

    subroutine inner_partition(workspace, xi, nlevels, count, status)
        type(triangle_level_workspace_t), intent(inout) :: workspace
        real(dp), intent(in) :: xi
        integer, intent(in) :: nlevels
        integer, intent(out) :: count, status
        real(dp) :: inner(3), upper, roots(2), empty(6)
        integer :: level, nroots
        empty = 0
        call generated_triangle_level_slice(empty, xi, inner, upper)
        workspace%inner_events(:2) = [0._dp, upper]
        count = 2
        status = 0
        do level = 1, nlevels
            call generated_triangle_level_slice(workspace%cuts(:, level), xi, inner, upper)
            call interval_roots(inner, 0._dp, upper, roots, nroots, status)
            if (status /= 0) return
            workspace%inner_events(count + 1:count + nroots) = roots(:nroots)
            count = count + nroots
        end do
        call sort_unique(workspace%inner_events, count)
    end subroutine inner_partition

    subroutine evaluate_inner(workspace, xi, nlevels, integrand, status, &
            record_trace, outer_weight)
        type(triangle_level_workspace_t), intent(inout) :: workspace
        real(dp), intent(in) :: xi, outer_weight
        integer, intent(in) :: nlevels
        procedure(triangle_level_integrand_t) :: integrand
        integer, intent(out) :: status
        logical, intent(in) :: record_trace
        real(dp) :: lower, upper, eta, weight, point(2), base_lower, base_upper, ignored
        integer :: ncuts, cut, node, callback_status, subdivision, nsub, low_order, high_order
        real(dp) :: unit_node, unit_weight
        call inner_partition(workspace, xi, nlevels, ncuts, status)
        if (status /= 0) return
        workspace%scratch(:, 4:5) = 0
        nsub = 1
        low_order = 4
        high_order = 8
        if (workspace%inner_depth > 0) then
            ! Upgrade the embedded pair before uniform subdivision; smooth rational
            ! branches converge much faster with degree fifteen low-rule exactness.
            nsub = 2**(workspace%inner_depth - 1)
            low_order = 8
            high_order = 16
        end if
        do cut = 1, ncuts - 1
            base_lower = workspace%inner_events(cut)
            base_upper = workspace%inner_events(cut + 1)
            do subdivision = 1, nsub
                call generated_triangle_level_interval(base_lower, base_upper, &
                    real(subdivision - 1, dp)/real(nsub, dp), 1._dp, lower, ignored)
                call generated_triangle_level_interval(base_lower, base_upper, &
                    real(subdivision, dp)/real(nsub, dp), 1._dp, upper, ignored)
                if (upper <= lower) cycle
                workspace%scratch(:, 1:2) = 0
                if (.not. record_trace) then
                    do node = 1, low_order
                        if (low_order == 4) then
                            unit_node = workspace%nodes4(node)
                            unit_weight = workspace%weights4(node)
                        else
                            unit_node = workspace%nodes8(node)
                            unit_weight = workspace%weights8(node)
                        end if
                        call generated_triangle_level_interval(lower, upper, &
                            unit_node, unit_weight, eta, weight)
                        point = [xi, eta]
                        call callback(point, callback_status)
                        if (callback_status /= 0) then
                            status = callback_status
                            return
                        end if
                        workspace%scratch(:, 1) = workspace%scratch(:, 1) + &
                            weight*workspace%scratch(:, 8)
                    end do
                end if
                do node = 1, high_order
                    if (high_order == 8) then
                        unit_node = workspace%nodes8(node)
                        unit_weight = workspace%weights8(node)
                    else
                        unit_node = workspace%nodes16(node)
                        unit_weight = workspace%weights16(node)
                    end if
                    call generated_triangle_level_interval(lower, upper, &
                        unit_node, unit_weight, eta, weight)
                    point = [xi, eta]
                    if (record_trace) then
                        if (workspace%npoints == workspace%max_points) then
                            status = TRIANGLE_LEVEL_BUDGET
                            return
                        end if
                        workspace%npoints = workspace%npoints + 1
                        workspace%points(:, workspace%npoints) = point
                        workspace%weights(workspace%npoints) = outer_weight*weight
                    else
                        call callback(point, callback_status)
                        if (callback_status /= 0) then
                            status = callback_status
                            return
                        end if
                        workspace%scratch(:, 2) = workspace%scratch(:, 2) + &
                            weight*workspace%scratch(:, 8)
                    end if
                end do
                workspace%scratch(:, 4) = workspace%scratch(:, 4) + &
                    workspace%scratch(:, 2)
                workspace%scratch(:, 5) = workspace%scratch(:, 5) + &
                    abs(workspace%scratch(:, 2) - workspace%scratch(:, 1))
            end do
        end do
    contains
        subroutine callback(point, status)
            real(dp), intent(in) :: point(2)
            integer, intent(out) :: status
            status = TRIANGLE_LEVEL_BUDGET
            if (workspace%nevaluations >= workspace%max_evaluations) return
            workspace%nevaluations = workspace%nevaluations + 1
            call integrand(point, workspace%scratch(:, 8), status)
            if (status /= 0) then
                status = TRIANGLE_LEVEL_CALLBACK
                return
            end if
            if (.not. all(ieee_is_finite(workspace%scratch(:, 8)))) then
                status = TRIANGLE_LEVEL_CALLBACK
                return
            end if
            if (maxval(abs(workspace%scratch(:, 8))) > huge(1._dp)/64) then
                status = TRIANGLE_LEVEL_CALLBACK
            end if
        end subroutine callback
    end subroutine evaluate_inner

    subroutine evaluate_panel(workspace, lower, upper, nlevels, integrand, panel, status)
        type(triangle_level_workspace_t), intent(inout) :: workspace
        real(dp), intent(in) :: lower, upper
        integer, intent(in) :: nlevels, panel
        procedure(triangle_level_integrand_t) :: integrand
        integer, intent(out) :: status
        real(dp) :: xi, weight, unit_node, unit_weight
        integer :: node, low_order, high_order
        workspace%scratch(:, 6:7) = 0
        workspace%panel_inner_error(:, panel) = 0
        low_order = 8
        high_order = 16
        if (workspace%inner_depth > 0 .or. workspace%panel_tangent(panel)) then
            low_order = 16
            high_order = 32
        end if
        do node = 1, low_order
            if (low_order == 8) then
                unit_node = workspace%nodes8(node)
                unit_weight = workspace%weights8(node)
            else
                unit_node = workspace%nodes16(node)
                unit_weight = workspace%weights16(node)
            end if
            call outer_node(workspace, panel, lower, upper, &
                unit_node, unit_weight, xi, weight)
            call evaluate_inner(workspace, xi, nlevels, integrand, status, .false., 0._dp)
            if (status /= 0) return
            workspace%scratch(:, 7) = workspace%scratch(:, 7) + &
                weight*workspace%scratch(:, 4)
        end do
        do node = 1, high_order
            if (high_order == 16) then
                unit_node = workspace%nodes16(node)
                unit_weight = workspace%weights16(node)
            else
                unit_node = workspace%nodes32(node)
                unit_weight = workspace%weights32(node)
            end if
            call outer_node(workspace, panel, lower, upper, &
                unit_node, unit_weight, xi, weight)
            call evaluate_inner(workspace, xi, nlevels, integrand, status, .false., 0._dp)
            if (status /= 0) return
            workspace%scratch(:, 6) = workspace%scratch(:, 6) + &
                weight*workspace%scratch(:, 4)
            workspace%panel_inner_error(:, panel) = &
                workspace%panel_inner_error(:, panel) + weight*workspace%scratch(:, 5)
        end do
        workspace%panel_value(:, panel) = workspace%scratch(:, 6)
        workspace%panel_error(:, panel) = workspace%panel_inner_error(:, panel) + &
            abs(workspace%scratch(:, 6) - workspace%scratch(:, 7))
        workspace%panel_lower(panel) = lower
        workspace%panel_upper(panel) = upper
        status = 0
    end subroutine evaluate_panel

    subroutine outer_node(workspace, panel, lower, upper, unit_node, unit_weight, &
            node, weight)
        type(triangle_level_workspace_t), intent(in) :: workspace
        integer, intent(in) :: panel
        real(dp), intent(in) :: lower, upper, unit_node, unit_weight
        real(dp), intent(out) :: node, weight
        if (workspace%panel_tangent(panel)) then
            call generated_triangle_level_tangent_interval(lower, upper, &
                unit_node, unit_weight, node, weight)
        else
            call generated_triangle_level_interval(lower, upper, &
                unit_node, unit_weight, node, weight)
        end if
    end subroutine outer_node

    pure real(dp) function component_budget(absolute, relative, integral) result(budget)
        real(dp), intent(in) :: absolute, relative, integral
        real(dp) :: magnitude
        ! Saturation preserves a valid loose tolerance without infinite arithmetic.
        magnitude = abs(integral)
        budget = huge(1._dp)
        if (magnitude > 1) then
            if (relative > huge(1._dp)/magnitude) return
        end if
        budget = relative*magnitude
        if (absolute > huge(1._dp) - budget) then
            budget = huge(1._dp)
        else
            budget = budget + absolute
        end if
    end function component_budget

    subroutine integrate_triangle_levels(coefficients, levels, integrand, epsabs, &
            epsrel, workspace, integral, error, status, message)
        real(dp), intent(in) :: coefficients(6), levels(:), epsabs(:), epsrel
        procedure(triangle_level_integrand_t) :: integrand
        type(triangle_level_workspace_t), intent(inout) :: workspace
        real(dp), intent(out) :: integral(:), error(:)
        integer, intent(out) :: status
        character(*), intent(out) :: message
        real(dp) :: edge(3, 2), discriminant(3), roots(2), scale, target
        real(dp) :: lower, upper, midpoint, xi, weight, worst, score, ratio
        integer :: nlevels, count, nroots, level, edge_id, panel, component, chosen, node
        integer :: trace_order, ntangencies
        real(dp) :: unit_node, unit_weight
        logical :: converged, refine_inner
        integral = 0
        error = 0
        message = 'invalid triangle level quadrature input or workspace'
        status = TRIANGLE_LEVEL_INVALID
        workspace%npoints = 0
        workspace%nevaluations = 0
        workspace%npanels = 0
        workspace%inner_depth = 0
        if (.not. workspace%initialized) return
        if (size(integral) /= workspace%nvalues) return
        if (size(error) /= workspace%nvalues) return
        if (size(epsabs) /= workspace%nvalues) return
        nlevels = size(levels)
        if (nlevels > workspace%max_levels) return
        if (.not. all(ieee_is_finite(coefficients))) return
        if (.not. all(ieee_is_finite(levels))) return
        if (.not. all(ieee_is_finite(epsabs))) return
        if (.not. ieee_is_finite(epsrel)) return
        if (any(epsabs < 0) .or. epsrel < 0) return
        if (all(epsabs <= 0) .and. epsrel <= 0) return
        if (maxval(abs(coefficients)) > huge(1._dp)/64) return
        if (nlevels > 0) then
            if (maxval(abs(levels)) > huge(1._dp)/64) return
        end if
        do level = 2, nlevels
            if (levels(level) < levels(level - 1)) return
        end do
        ! Shift the level before normalization, preserving common large gauges.
        workspace%events(:2) = [0._dp, 1._dp]
        count = 2
        ntangencies = 0
        do level = 1, nlevels
            workspace%cuts(:, level) = coefficients
            workspace%cuts(1, level) = coefficients(1) - levels(level)
            scale = maxval(abs(workspace%cuts(:, level)))
            if (scale > 0) workspace%cuts(:, level) = workspace%cuts(:, level)/scale
            call generated_triangle_level_events(workspace%cuts(:, level), &
                edge, discriminant)
            do edge_id = 1, 2
                call interval_roots(edge(:, edge_id), 0._dp, 1._dp, &
                    roots, nroots, status)
                if (status /= 0) goto 900
                workspace%events(count + 1:count + nroots) = roots(:nroots)
                count = count + nroots
            end do
            call interval_roots(discriminant, 0._dp, 1._dp, roots, nroots, status)
            if (status /= 0) goto 900
            workspace%tangencies(ntangencies + 1:ntangencies + nroots) = roots(:nroots)
            ntangencies = ntangencies + nroots
            workspace%events(count + 1:count + nroots) = roots(:nroots)
            count = count + nroots
        end do
        call sort_unique(workspace%events, count)
        if (count - 1 > workspace%max_panels) then
            status = TRIANGLE_LEVEL_BUDGET
            message = 'initial level topology exceeds the panel budget'
            goto 910
        end if
        workspace%npanels = count - 1
        do panel = 1, workspace%npanels
            workspace%panel_tangent(panel) = &
                any(abs(workspace%tangencies(:ntangencies) - workspace%events(panel)) <= 0) &
                .or. any(abs(workspace%tangencies(:ntangencies) - workspace%events(panel + 1)) <= 0)
            call evaluate_panel(workspace, workspace%events(panel), &
                workspace%events(panel + 1), nlevels, integrand, panel, status)
            if (status /= 0) goto 900
        end do
        do
            integral = 0
            error = 0
            workspace%scratch(:, 3) = 0
            do panel = 1, workspace%npanels
                integral = integral + workspace%panel_value(:, panel)
                error = error + workspace%panel_error(:, panel)
                workspace%scratch(:, 3) = workspace%scratch(:, 3) + &
                    workspace%panel_inner_error(:, panel)
            end do
            converged = .true.
            refine_inner = .false.
            do component = 1, workspace%nvalues
                target = component_budget(epsabs(component), epsrel, integral(component))
                if (error(component) > target) converged = .false.
                if (workspace%scratch(component, 3) > target/2) refine_inner = .true.
            end do
            if (converged) exit
            if (refine_inner) then
                if (workspace%inner_depth >= 15) then
                    status = TRIANGLE_LEVEL_BUDGET
                    message = 'inner quadrature resolution limit reached'
                    goto 910
                end if
                workspace%inner_depth = workspace%inner_depth + 1
                do panel = 1, workspace%npanels
                    lower = workspace%panel_lower(panel)
                    upper = workspace%panel_upper(panel)
                    call evaluate_panel(workspace, lower, upper, nlevels, &
                        integrand, panel, status)
                    if (status /= 0) goto 900
                end do
                cycle
            end if
            if (workspace%npanels == workspace%max_panels) then
                status = TRIANGLE_LEVEL_BUDGET
                message = 'adaptive outer panel budget exhausted'
                goto 910
            end if
            worst = -1
            chosen = 1
            do panel = 1, workspace%npanels
                score = 0
                do component = 1, workspace%nvalues
                    target = component_budget(epsabs(component), epsrel, integral(component))
                    ratio = 0
                    if (target > 0) then
                        ratio = workspace%panel_error(component, panel)
                        if (target < 1) then
                            if (ratio > huge(1._dp)*target) then
                                ratio = huge(1._dp)
                            else
                                ratio = ratio/target
                            end if
                        else
                            ratio = ratio/target
                        end if
                    else if (workspace%panel_error(component, panel) > 0) then
                        ratio = huge(1._dp)
                    end if
                    score = max(score, ratio)
                end do
                if (score > worst) then
                    worst = score
                    chosen = panel
                end if
            end do
            lower = workspace%panel_lower(chosen)
            upper = workspace%panel_upper(chosen)
            call generated_triangle_level_interval(lower, upper, .5_dp, 1._dp, &
                midpoint, weight)
            if (midpoint <= lower .or. midpoint >= upper) then
                status = TRIANGLE_LEVEL_BUDGET
                message = 'outer quadrature floating-point resolution limit reached'
                goto 910
            end if
            workspace%npanels = workspace%npanels + 1
            workspace%panel_tangent(workspace%npanels) = workspace%panel_tangent(chosen)
            call evaluate_panel(workspace, lower, midpoint, nlevels, &
                integrand, chosen, status)
            if (status /= 0) goto 900
            call evaluate_panel(workspace, midpoint, upper, nlevels, &
                integrand, workspace%npanels, status)
            if (status /= 0) goto 900
        end do
        ! Materialize accepted high-rule nodes without reevaluating the callback.
        do panel = 1, workspace%npanels
            trace_order = 16
            if (workspace%inner_depth > 0 .or. workspace%panel_tangent(panel)) trace_order = 32
            do node = 1, trace_order
                if (trace_order == 16) then
                    unit_node = workspace%nodes16(node)
                    unit_weight = workspace%weights16(node)
                else
                    unit_node = workspace%nodes32(node)
                    unit_weight = workspace%weights32(node)
                end if
                call outer_node(workspace, panel, workspace%panel_lower(panel), &
                    workspace%panel_upper(panel), unit_node, unit_weight, xi, weight)
                call evaluate_inner(workspace, xi, nlevels, integrand, status, &
                    .true., weight)
                if (status /= 0) goto 900
            end do
        end do
        status = 0
        message = ''
        return
        900     continue
        select case (status)
        case (TRIANGLE_LEVEL_CALLBACK)
            message = 'integrand callback failed or returned an unrepresentable value'
        case default
            message = 'quadrature evaluation, trace capacity or root resolution budget exhausted'
        end select
        910     continue
        integral = 0
        workspace%npoints = 0
    end subroutine integrate_triangle_levels
end module fortfem_triangle_level_quadrature
