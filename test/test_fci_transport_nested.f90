program test_fci_transport_nested
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_fci_field_line_tracer, only: trace_fci_field_line_rk4
    use fortfem_fci_interpolation_map, only: &
        build_fci_bilinear_interpolation_map_2d
    use fortfem_fci_parallel_operator, only: assemble_fci_parallel_gradient_csc
    use fortfem_fci_transport_reference, only: fci_transport_reference_t, &
        setup_fci_transport_reference, solve_fci_transport_reference
    use fortsparse, only: csc_t, csc_from_triplet, csc_matvec, fortsparse_status_t
    use omp_lib, only: omp_get_wtime
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp), dz = 0.5_dp
    character(16) :: mode
    integer :: refinement

    call get_command_argument(1, mode)
    write (*, '(a)') 'cells,unknowns,epsilon,relative_mass_error,'// &
        'scaled_residual,weak_residual,parallel_defect,baseline_error,'// &
        'trace_flux_defect,trace_step_difference,setup_s,solve_s'
    do refinement = 1, 3
        call nested_ladder(4*2**(refinement - 1) + 1, trim(mode) == 'gate')
    end do
    call check_summary('FCI nested-surface transport')
contains
    pure real(dp) function flux(x, y)
        real(dp), intent(in) :: x, y
        flux = sin(pi*x)*sin(pi*y)
    end function flux

    pure subroutine field_line_rhs(z, point, derivative)
        real(dp), intent(in) :: z, point(:)
        real(dp), intent(out) :: derivative(:)
        real(dp) :: speed, bx, by

        ! Independent analytic fixture, not a production transport leaf.
        associate (unused => z)
            speed = 0.15_dp*(1.0_dp + flux(point(1), point(2)))
            bx = speed*pi*sin(pi*point(1))*cos(pi*point(2))
            by = -speed*pi*cos(pi*point(1))*sin(pi*point(2))
            derivative = [bx, by, sqrt(1.0_dp + bx*bx + by*by)]
        end associate
    end subroutine field_line_rhs

    subroutine nested_ladder(cells, enforce_gate)
        integer, intent(in) :: cells
        logical, intent(in) :: enforce_gate
        type(csc_t) :: gradient, unfolded, perpendicular
        type(fci_transport_reference_t) :: solver
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: grid(:), tx(:), ty(:), sx(:), sy(:)
        real(dp), allocatable :: forward(:, :, :), backward(:, :, :), ell(:, :)
        real(dp), allocatable :: mass(:), weights(:), exact(:), rhs(:), solution(:)
        real(dp), allocatable :: values(:), boundary(:), image(:), weak_image(:)
        integer, allocatable :: rows(:), columns(:), anchors(:)
        logical, allocatable :: reaches_boundary(:)
        real(dp) :: point(3), endpoint(3), check_endpoint(3)
        real(dp) :: eps, error, residual, defect, baseline, t, setup_s, solve_s
        real(dp) :: weak_residual
        real(dp) :: h, trace_defect, trace_difference, eigenvalue, factor
        integer :: side, plane_size, sample_count, n, i, j, row, k, col
        integer :: plane, exponent, count, anchor_count, global, neighbor, sweep
        logical :: trace_ok

        side = cells + 1
        plane_size = side*side
        sample_count = (cells - 1)**2
        n = 2*plane_size
        h = 1.0_dp/cells
        allocate (grid(side), sx(sample_count), sy(sample_count))
        allocate (tx(sample_count), ty(sample_count))
        allocate (forward(sample_count, plane_size, 2), source=0.0_dp)
        allocate (backward(sample_count, plane_size, 2), source=0.0_dp)
        allocate (ell(sample_count, 2))
        grid = [(real(i - 1, dp)*h, i=1, side)]
        row = 0
        trace_defect = 0.0_dp
        trace_difference = 0.0_dp
        trace_ok = .true.
        do j = 2, cells
            do i = 2, cells
                row = row + 1
                sx(row) = grid(i)
                sy(row) = grid(j)
                point = [sx(row), sy(row), 0.0_dp]
                call trace_fci_field_line_rk4(point, dz, 64, &
                    field_line_rhs, endpoint, status)
                trace_ok = trace_ok .and. status%code == 0
                call trace_fci_field_line_rk4(point, dz, 128, &
                    field_line_rhs, check_endpoint, status)
                trace_ok = trace_ok .and. status%code == 0
                trace_difference = max(trace_difference, &
                    maxval(abs(endpoint - check_endpoint)))
                trace_defect = max(trace_defect, &
                    abs(flux(endpoint(1), endpoint(2)) - flux(sx(row), sy(row))))
                tx(row) = endpoint(1)
                ty(row) = endpoint(2)
                ell(row, :) = endpoint(3)
            end do
        end do
        call check_condition(trace_ok, 'trace nested field')
        call check_condition(trace_defect < 1.0e-9_dp, &
            'traced endpoints preserve analytic flux')
        call check_condition(trace_difference < 1.0e-9_dp, &
            'trace step doubling controls endpoint error')
        call build_fci_bilinear_interpolation_map_2d(grid, grid, tx, ty, &
            forward(:, :, 1), status)
        call check_condition(status%code == 0, 'build traced endpoint map')
        call build_fci_bilinear_interpolation_map_2d(grid, grid, sx, sy, &
            backward(:, :, 1), status)
        call check_condition(status%code == 0, 'build departure map')
        forward(:, :, 2) = forward(:, :, 1)
        backward(:, :, 2) = backward(:, :, 1)
        allocate (reaches_boundary(plane_size), source=.false.)
        do j = 1, side
            do i = 1, side
                global = (j - 1)*side + i
                reaches_boundary(global) = &
                    i == 1 .or. i == side .or. j == 1 .or. j == side
            end do
        end do
        ! Positive interpolation leaks every sampled closed orbit to the wall.
        do sweep = 1, plane_size
            row = 0
            do j = 2, cells
                do i = 2, cells
                    row = row + 1
                    global = (j - 1)*side + i
                    if (any(forward(row, :, 1) > 0.0_dp .and. &
                        reaches_boundary)) reaches_boundary(global) = .true.
                end do
            end do
            if (all(reaches_boundary)) exit
        end do
        call check_condition(all(reaches_boundary), &
            'interpolation graph connects every interior node to the wall')
        call assemble_fci_parallel_gradient_csc(forward, backward, ell, &
            unfolded, status)
        call check_condition(status%code == 0, 'assemble native FCI gradient')
        allocate (rows(unfolded%nnz), columns(unfolded%nnz))
        do col = 1, unfolded%ncol
            do k = unfolded%col_ptr(col), unfolded%col_ptr(col + 1) - 1
                rows(k) = unfolded%row_idx(k)
                ! Identify the two copies of the periodic end plane.
                columns(k) = 1 + modulo(col - 1, n)
            end do
        end do
        call csc_from_triplet(unfolded%nrow, n, rows, columns, &
            unfolded%val, gradient, status)
        call check_condition(status%code == 0, 'periodic identification')
        deallocate (rows, columns)
        allocate (rows(6*n), columns(6*n), values(6*n))
        allocate (mass(n), exact(n), rhs(n), solution(n), weights(2*sample_count))
        allocate (anchors(n - 2*sample_count))
        mass = h*h*dz
        weights = h*h*dz
        count = 0
        anchor_count = 0
        do plane = 0, 1
            do j = 1, side
                do i = 1, side
                    global = plane*plane_size + (j - 1)*side + i
                    exact(global) = flux(grid(i), grid(j))
                    if (i == 1 .or. i == side .or. j == 1 .or. j == side) then
                        anchor_count = anchor_count + 1
                        anchors(anchor_count) = global
                        exact(global) = 0.0_dp
                    end if
                    count = count + 1
                    rows(count) = global
                    columns(count) = global
                    values(count) = 4.0_dp*dz + 2.0_dp*h*h/dz
                    count = count + 1
                    rows(count) = global
                    columns(count) = modulo(global - 1 + plane_size, n) + 1
                    values(count) = -2.0_dp*h*h/dz
                    do k = 1, 4
                        neighbor = 0
                        if (k == 1 .and. i > 1) neighbor = global - 1
                        if (k == 2 .and. i < side) neighbor = global + 1
                        if (k == 3 .and. j > 1) neighbor = global - side
                        if (k == 4 .and. j < side) neighbor = global + side
                        if (neighbor == 0) cycle
                        count = count + 1
                        rows(count) = global
                        columns(count) = neighbor
                        values(count) = -dz
                    end do
                end do
            end do
        end do
        call csc_from_triplet(n, n, rows(:count), columns(:count), &
            values(:count), perpendicular, status)
        call check_condition(status%code == 0, 'plane Laplacian assembly')
        boundary = 0.0_dp*real(anchors, dp)
        image = csc_matvec(gradient, exact)
        defect = sqrt(sum(weights*image**2)/sum(mass*exact**2))
        ! The sampled sine is an independent five-point eigenvector.
        eigenvalue = 8.0_dp*sin(pi*h/2.0_dp)**2/h**2
        image = csc_matvec(perpendicular, exact) - eigenvalue*mass*exact
        image(anchors) = 0.0_dp
        call check_condition(maxval(abs(image)) < 1.0e-13_dp, &
            'plane operator matches analytic difference eigenvalue')
        rhs = 2.0_dp*pi*pi*mass*exact
        factor = 2.0_dp*pi*pi/eigenvalue
        baseline = abs(factor - 1.0_dp)
        t = omp_get_wtime()
        call setup_fci_transport_reference(gradient, weights, mass, &
            perpendicular, anchors, 0, solver, status)
        setup_s = omp_get_wtime() - t
        call check_condition(status%code == 0, 'nested discrete nullspace setup')
        if (status%code /= 0) then
            write (*, '(a)') status%msg
            return
        end if
        ! The oracle requires positive excess conductivity, even for this limit.
        call solve_fci_transport_reference(solver, 0.0_dp, 1.0e-20_dp, 1.0_dp, &
            rhs, boundary, solution, residual, status)
        call check_condition(status%code == 0, 'isotropic baseline solve')
        call check_condition(maxval(abs(solution - factor*exact)) < 1.0e-11_dp, &
            'baseline solve matches independent sine eigenvalue')
        do exponent = 6, 13
            eps = 10.0_dp**(-exponent)
            t = omp_get_wtime()
            call solve_fci_transport_reference(solver, 0.0_dp, &
                1.0_dp/eps - 1.0_dp, 1.0_dp, rhs, boundary, &
                solution, residual, status)
            solve_s = omp_get_wtime() - t
            call check_condition(status%code == 0, 'nested oracle solve')
            call check_condition(residual < 1.0e-10_dp, 'nested scaled residual')
            image = weights*csc_matvec(gradient, solution)
            weak_image = csc_matvec(perpendicular, solution) - rhs
            do col = 1, n
                do k = gradient%col_ptr(col), gradient%col_ptr(col + 1) - 1
                    weak_image(col) = weak_image(col) + &
                        (1.0_dp/eps - 1.0_dp)*gradient%val(k)* &
                        image(gradient%row_idx(k))
                end do
            end do
            weak_image(anchors) = 0.0_dp
            weak_residual = sqrt(sum(weak_image**2)/sum(rhs**2))
            call check_condition(weak_residual < 1.0e-9_dp, &
                'oracle satisfies original weak equations on free nodes')
            error = sqrt(sum(mass*(solution - exact)**2)/sum(mass*exact**2))
            write (*, '(i0,",",i0,10(",",es16.8))') cells, n, eps, error, &
                residual, weak_residual, defect, baseline, trace_defect, &
                trace_difference, setup_s, solve_s
            if (enforce_gate) call check_condition(error < 0.1_dp, &
                'gate 1: nested relative mass error below ten percent')
        end do
    end subroutine nested_ladder
end program test_fci_transport_nested
