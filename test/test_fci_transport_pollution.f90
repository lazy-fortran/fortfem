program test_fci_transport_pollution
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use check, only: check_condition, check_summary
    use fortfem_gll_line, only: build_gll_line, gll_rule, lagrange_jet
    use fortfem_tensor_line_tables, only: tensor_line_tables_t
    use fortfem_generated_pollution_sovinec, only: pollution_sovinec
    use fortfem_generated_pollution_circle, only: pollution_circle
    use fortfem_generated_pollution_basis, only: pollution_basis
    use fortfem_generated_pollution_integrand, only: pollution_integrand
    use fortfem_generated_pollution_error, only: pollution_error
    use fortfem_fci_field_line_tracer, only: trace_fci_field_line_rk4
    use fortfem_fci_interpolation_map, only: build_fci_bilinear_interpolation_map_2d
    use fortfem_generated_fci_parallel_gradient, only: generated_fci_parallel_gradient
    use fortfem_sparse_direct, only: sparse_direct_factor_t, &
        sparse_direct_factor_csc, sparse_direct_solve_factored, sparse_direct_free
    use fortsparse, only: sparse_options_t, FORTSPARSE_SYM_SPD, fortsparse_status_t
    implicit none

    real(dp), parameter :: residual_limit = 1.0e3_dp*epsilon(1.0_dp)

    type :: fixture_t
        real(dp), allocatable :: slow(:, :), parallel(:, :), load(:)
        real(dp), allocatable :: boundary_slow(:), boundary_parallel(:)
        integer, allocatable :: col_ptr(:), row_idx(:)
        real(dp), allocatable :: sparse_slow(:), sparse_parallel(:)
        integer :: centre = 0, cells = 0, degree = 0
        real(dp) :: setup_seconds = 0.0_dp
    end type
    interface
        subroutine dpbsv(uplo, n, kd, nrhs, ab, ldab, b, ldb, info)
            import dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, kd, nrhs, ldab, ldb
            real(dp), intent(inout) :: ab(ldab, *), b(ldb, *)
            integer, intent(out) :: info
        end subroutine
    end interface
    type(fixture_t) :: fixture
    character(4096) :: mode, output
    character(16) :: trace_benchmark
    real(dp) :: centre, pollution, residual, seconds, expected, pi, h
    integer :: cells, unit, exponent, degree

    pi = acos(-1.0_dp)
    call get_command_argument(1, mode)
    if (trim(mode) == 'sparse') then
        call sparse_oracle()
    else if (trim(mode) == 'cost') then
        call cost_ladder()
    else if (trim(mode) == 'parallel') then
        call get_command_argument(2, output)
        call parallel_gate(trim(output))
    else if (trim(mode) == 'isotropic') then
        call get_command_argument(2, output)
        call isotropic_gate(trim(output))
    else if (trim(mode) == 'scaling') then
        call get_command_argument(2, output)
        read (output, *) degree
        call get_command_argument(3, output)
        call conforming_report(degree, trim(output))
    else if (trim(mode) == 'paper') then
        call get_command_argument(2, output)
        if (len_trim(output) == 0) error stop 'paper mode needs a CSV path'
        open (newunit=unit, file=trim(output), status='replace')
        write (unit, '(a)') 'benchmark,method,p,cells,unknowns,chi_ratio,'// &
            'centre,effective_perp_error,residual,setup_s,solve_s,solve_s_per_unknown'
        call assemble_fixture('sovinec', 'h1', 1, 8, fixture)
        call solve_fixture(fixture, 1.0e3_dp, centre, pollution, residual, seconds)
        call write_row(unit, 'sovinec', 'h1', fixture, 1.0e3_dp, &
            centre, pollution, residual, seconds)
        close (unit)
        call check_condition(residual <= residual_limit, &
            'normwise original weak residual')
    else
        call residual_oracle()
        call get_command_argument(2, output)
        unit = 0
        if (len_trim(output) > 0) then
            open (newunit=unit, file=trim(output), status='replace')
            write (unit, '(a)') 'cells,chi_ratio,centre,'// &
                'normwise_residual,residual_limit'
        end if
        do cells = 8, 16, 8
            call assemble_fixture('sovinec', 'symmetric', 1, cells, fixture)
            h = 1.0_dp/real(cells, dp)
            expected = (pi*h/(2.0_dp*sin(pi*h/2.0_dp)))**2
            do exponent = 0, 6, 3
                call solve_fixture(fixture, 10.0_dp**exponent, &
                    centre, pollution, residual, seconds)
                write (*, '(a,i0,a,i0,2(a,es16.8))') 'sine cells=', cells, &
                    ' exponent=', exponent, ' centre=', centre, ' residual=', residual
                if (unit /= 0) write (unit, '(i0,4(",",es24.16))') cells, &
                    10.0_dp**exponent, centre, residual, residual_limit
                call check_condition(abs(centre - expected) < &
                    2.0e-9_dp*expected, 'symmetric sine eigenmode oracle')
                call check_condition(residual <= residual_limit, &
                    'symmetric original residual')
            end do
        end do
        if (unit /= 0) close (unit)
        call assemble_fixture('sovinec', 'h1', 1, 8, fixture)
        call solve_fixture(fixture, 1.0_dp, centre, pollution, residual, seconds)
        call check_condition(abs(centre - 1.0_dp) < 0.02_dp, &
            'H1 isotropic central value against analytic solution')
        call check_condition(residual <= residual_limit, 'H1 original residual')
    end if
    call check_summary('anisotropic pollution')
contains
    subroutine residual_oracle()
        real(dp) :: band(2, 2), x(2), rhs(2), r
        ! The unused band corner must not enter the matrix norm.
        band(:, 1) = [3.0_dp, 2.0_dp]
        band(:, 2) = [6.0_dp, huge(1.0_dp)]
        x = [2.0_dp, -1.0_dp]
        rhs = [1.0_dp, 4.0_dp]
        r = backward_residual(band, x, rhs)
        call check_condition(abs(r - 3.0_dp/10.0_dp) < epsilon(r), &
            'normwise residual against exact 2x2 oracle')
        band(2, 2) = 0.0_dp
        r = backward_residual(1.0e6_dp*band, x, 1.0e6_dp*rhs)
        call check_condition(abs(r - 3.0_dp/10.0_dp) < epsilon(r), &
            'normwise residual invariant under equation scaling')
    end subroutine

    subroutine gate(passed, message)
        logical, intent(in) :: passed
        character(*), intent(in) :: message
        call check_condition(passed, message)
        if (.not. passed) call check_summary('anisotropic pollution')
    end subroutine

    subroutine isotropic_gate(path)
        character(*), intent(in) :: path
        type(fixture_t) :: data
        real(dp), allocatable :: solution(:)
        real(dp) :: c, e, r, s, errors(4), rate
        integer :: p, level, ne, out
        if (len_trim(path) == 0) error stop 'isotropic mode needs a CSV path'
        open (newunit=out, file=path, status='replace')
        write (out, '(a)') 'p,cells,unknowns,l2_error,rate,normwise_residual'
        do p = 1, 4
            do level = 1, 4
                ne = 2**level
                call assemble_fixture('sovinec', 'h1', p, ne, data)
                call solve_fixture(data, 1.0_dp, c, e, r, s, solution)
                errors(level) = l2_error('sovinec', p, ne, solution)
                rate = 0.0_dp
                if (level > 1) rate = log(errors(level - 1)/errors(level))/log(2.0_dp)
                write (out, '(3(i0,","),3(es24.16,:,","))') p, ne, &
                    size(solution), errors(level), rate, r
                flush (out)
                write (*, '(a,i0,a,i0,3(a,es14.6))') 'isotropic p=', p, &
                    ' cells=', ne, ' L2=', errors(level), ' order=', rate, &
                    ' residual=', r
                call gate(r <= residual_limit, 'isotropic normwise residual')
                if (level >= 3) call gate(abs(rate - real(p + 1, dp)) <= 0.3_dp, &
                    'isotropic L2 design order p+1 within 0.3')
            end do
        end do
        close (out)
    end subroutine

    subroutine parallel_gate(path)
        character(*), intent(in) :: path
        type(tensor_line_tables_t) :: line
        real(dp), allocatable :: nodes(:), weights(:), nodal(:, :)
        real(dp) :: errors(4), rate, h, x, y, value, source, bx, by
        real(dp) :: v, gx, gy, dx, dy, slow, parallel, load, exact
        integer :: out, p, level, ne, status, ex, ey, ix, iy, qx, qy
        open (newunit=out, file=path, status='replace')
        write (out, '(a)') 'p,cells,l2_parallel_error,rate'
        do p = 1, 4
            call gll_rule(p, nodes, weights, status)
            if (status /= 0) error stop 'parallel interpolation nodes failed'
            allocate (nodal(p + 1, p + 1))
            do level = 1, 4
                ne = 2**level
                h = 1.0_dp/real(ne, dp)
                call build_gll_line(p, ne, line, status, nq=p + 4, &
                    x0=-0.5_dp, x1=0.5_dp)
                if (status /= 0) error stop 'parallel quadrature failed'
                errors(level) = 0.0_dp
                do ey = 1, ne
                    do ex = 1, ne
                        do iy = 1, p + 1
                            y = -0.5_dp + (real(ey - 1, dp) + &
                                0.5_dp*(nodes(iy) + 1.0_dp))*h
                            do ix = 1, p + 1
                                x = -0.5_dp + (real(ex - 1, dp) + &
                                    0.5_dp*(nodes(ix) + 1.0_dp))*h
                                nodal(ix, iy) = exp(x + 2.0_dp*y)
                            end do
                        end do
                        do qy = 1, line%nq
                            do qx = 1, line%nq
                                x = line%xq(qx, ex)
                                y = line%xq(qy, ey)
                                call benchmark('sovinec', x, y, value, source, bx, by)
                                dx = 0.0_dp
                                dy = 0.0_dp
                                do iy = 1, p + 1
                                    do ix = 1, p + 1
                                        call pollution_basis(line%b(qx, ix, ex), &
                                            line%b(qy, iy, ey), line%d(qx, ix, ex), &
                                            line%d(qy, iy, ey), v, gx, gy)
                                        dx = dx + nodal(ix, iy)*gx
                                        dy = dy + nodal(ix, iy)*gy
                                    end do
                                end do
                                ! Analytic derivatives of exp(x+2y) are independent.
                                exact = exp(x + 2.0_dp*y)
                                dx = dx - exact
                                dy = dy - 2.0_dp*exact
                                call pollution_integrand(dx, dy, dx, dy, bx, by, &
                                    line%wq(qx, ex)*line%wq(qy, ey), 0.0_dp, &
                                    0.0_dp, slow, parallel, load)
                                errors(level) = errors(level) + parallel
                            end do
                        end do
                    end do
                end do
                errors(level) = sqrt(errors(level))
                rate = 0.0_dp
                if (level > 1) rate = log(errors(level - 1)/errors(level))/log(2.0_dp)
                write (out, '(2(i0,","),2(es24.16,:,","))') &
                    p, ne, errors(level), rate
                flush (out)
                write (*, '(a,i0,a,i0,2(a,es14.6))') 'parallel p=', p, &
                    ' cells=', ne, ' L2=', errors(level), ' order=', rate
                if (level >= 3) call gate(abs(rate - real(p, dp)) <= 0.3_dp, &
                    'parallel interpolant L2 design order p within 0.3')
            end do
            deallocate (nodes, weights, nodal)
        end do
        close (out)
    end subroutine

    function l2_error(name, p, ne, solution) result(error)
        character(*), intent(in) :: name
        integer, intent(in) :: p, ne
        real(dp), intent(in) :: solution(:)
        real(dp) :: error, numerical, exact, source, bx, by, term, coefficient
        type(tensor_line_tables_t) :: line
        integer :: status, n, ex, ey, ix, iy, qx, qy, gi, gj
        ! Independent overintegration measures between-node error as well.
        call build_gll_line(p, ne, line, status, nq=p + 4, x0=-0.5_dp, x1=0.5_dp)
        if (status /= 0) error stop 'error quadrature failed'
        n = p*ne - 1
        error = 0.0_dp
        do ey = 1, ne
            do ex = 1, ne
                do qy = 1, line%nq
                    do qx = 1, line%nq
                        numerical = 0.0_dp
                        do iy = 1, p + 1
                            gj = (ey - 1)*p + iy - 1
                            do ix = 1, p + 1
                                gi = (ex - 1)*p + ix - 1
                                coefficient = 0.0_dp
                                if (gi > 0 .and. gi <= n .and. &
                                    gj > 0 .and. gj <= n) &
                                    coefficient = solution((gj - 1)*n + gi)
                                numerical = numerical + coefficient* &
                                    line%b(qx, ix, ex)*line%b(qy, iy, ey)
                            end do
                        end do
                        call benchmark(name, line%xq(qx, ex), line%xq(qy, ey), &
                            exact, source, bx, by)
                        call pollution_error(numerical, exact, &
                            line%wq(qx, ex)*line%wq(qy, ey), term)
                        error = error + term
                    end do
                end do
            end do
        end do
        error = sqrt(error)
    end function

    subroutine conforming_report(p, path)
        integer, intent(in) :: p
        character(*), intent(in) :: path
        type(fixture_t) :: data
        real(dp) :: c, e, r, s, baseline, chi, signed_error, h(4), values(4, 4)
        real(dp) :: aa(4), bb(4), a, b, constant, intercept, log_chi(4)
        integer :: out, level, power, ne, first
        logical :: fit_available
        fit_available = .true.
        if (p < 1 .or. p > 4) error stop 'scaling degree must be 1..4'
        open (newunit=out, file=path, status='replace')
        write (out, '(a)') 'benchmark,method,p,cells,unknowns,chi_ratio,'// &
            'centre,effective_perp_error,residual,setup_s,solve_s,'// &
            'solve_s_per_unknown,baseline_effective_perp,signed_pollution'
        first = 2
        if (p == 1) first = 3
        do level = 1, 4
            ne = 2**(first + level - 1)
            h(level) = log(1.0_dp/real(ne, dp))
            call assemble_fixture('sovinec', 'h1', p, ne, data)
            call solve_fixture(data, 1.0_dp, c, e, r, s)
            baseline = 1.0_dp/c
            call gate(r <= residual_limit, 'isotropic baseline normwise residual')
            do power = 1, 4
                chi = 10.0_dp**(6 + power)
                log_chi(power) = log(chi)
                call solve_fixture(data, chi, c, e, r, s)
                signed_error = 1.0_dp/c - baseline
                call write_row(out, 'sovinec', 'h1', data, chi, c, e, r, s, &
                    baseline, signed_error)
                flush (out)
                call gate(r <= residual_limit, 'scaling normwise original residual')
                if (signed_error > 0.0_dp) then
                    values(level, power) = log(signed_error)
                else
                    fit_available = .false.
                end if
            end do
        end do
        close (out)
        if (.not. fit_available) then
            write (*, '(a)') 'No logarithmic fit: signed pollution is nonpositive.'
            return
        end if
        do level = 1, 4
            bb(level) = slope(log_chi, values(level, :))
            aa(level) = slope(h, values(:, level))
        end do
        a = sum(aa)/4.0_dp
        b = sum(bb)/4.0_dp
        intercept = sum(values)/16.0_dp - a*sum(h)/4.0_dp - &
            b*sum(log_chi)/4.0_dp
        constant = exp(intercept)
        write (*, '(a,i0,3(a,es16.8))') 'conforming p=', p, &
            ' a=', a, ' b=', b, ' C=', constant
        write (*, '(a,4es16.8)') 'a by anisotropy: ', aa
        write (*, '(a,4es16.8)') 'b by mesh: ', bb
        ! The window crosses response regimes; the fit is descriptive only.
    end subroutine

    pure function slope(x, y) result(value)
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: value, mx, my
        mx = sum(x)/real(size(x), dp)
        my = sum(y)/real(size(y), dp)
        value = sum((x - mx)*(y - my))/sum((x - mx)**2)
    end function

    pure subroutine benchmark(name, x, y, value, source, bx, by)
        character(*), intent(in) :: name
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value, source, bx, by
        real(dp) :: vx, vy, norm
        if (name == 'sovinec') then
            call pollution_sovinec(x, y, pi, value, source, vx, vy)
        else
            call pollution_circle(x, y, value, source, vx, vy)
        end if
        norm = sqrt(vx*vx + vy*vy)
        bx = 0.0_dp
        by = 0.0_dp
        if (norm > tiny(norm)) then
            bx = vx/norm
            by = vy/norm
        end if
    end subroutine

    subroutine assemble_fixture(name, method, p, ne, data, sparse)
        character(*), intent(in) :: name, method
        integer, intent(in) :: p, ne
        type(fixture_t), intent(out) :: data
        logical, intent(in), optional :: sparse
        type(tensor_line_tables_t) :: line
        real(dp), allocatable :: v(:), gx(:), gy(:), boundary(:), nodes(:), weights(:)
        integer, allocatable :: ids(:)
        integer :: n, nd, bw, ex, ey, ix, iy, qx, qy, a, b, gi, gj, status
        integer(int64) :: start, finish, rate
        real(dp) :: value, source, bx, by, slow, parallel, load, weight, h, x, y
        real(dp) :: endpoints(2, 2), points(2), maps(2, 4), gridx(2), gridy(2), q(4)
        type(fortsparse_status_t) :: trace_status
        call system_clock(start, rate)
        n = p*ne - 1
        nd = n*n
        bw = min(p*n + p, nd - 1)
        data%cells = ne
        data%degree = p
        data%centre = (n/2)*n + n/2 + 1
        if (present(sparse)) then
            if (sparse) call sparse_pattern(p, ne, data)
        end if
        if (.not. allocated(data%col_ptr)) then
            allocate (data%slow(bw + 1, nd), data%parallel(bw + 1, nd))
            data%slow = 0.0_dp
            data%parallel = 0.0_dp
        end if
        allocate (data%load(nd), data%boundary_slow(nd), &
            data%boundary_parallel(nd))
        data%load = 0.0_dp
        data%boundary_slow = 0.0_dp
        data%boundary_parallel = 0.0_dp
        allocate (v((p + 1)**2), gx((p + 1)**2), gy((p + 1)**2))
        allocate (ids((p + 1)**2), boundary((p + 1)**2))
        call build_gll_line(p, ne, line, status, x0=-0.5_dp, x1=0.5_dp)
        if (status /= 0) error stop 'line construction failed'
        call gll_rule(p, nodes, weights, status)
        if (status /= 0) error stop 'boundary nodes failed'
        h = 1.0_dp/real(ne, dp)
        trace_benchmark = name
        do ey = 1, ne
            do ex = 1, ne
                a = 0
                do iy = 1, p + 1
                    do ix = 1, p + 1
                        a = a + 1
                        gi = (ex - 1)*p + ix - 1
                        gj = (ey - 1)*p + iy - 1
                        ids(a) = 0
                        boundary(a) = 0.0_dp
                        if (gi > 0 .and. gi <= n .and. &
                            gj > 0 .and. gj <= n) then
                            ids(a) = (gj - 1)*n + gi
                        else
                            x = -0.5_dp + (real(ex - 1, dp) + &
                                0.5_dp*(nodes(ix) + 1.0_dp))*h
                            y = -0.5_dp + (real(ey - 1, dp) + &
                                0.5_dp*(nodes(iy) + 1.0_dp))*h
                            call benchmark(name, x, y, boundary(a), &
                                source, bx, by)
                        end if
                    end do
                end do
                do qy = 1, line%nq
                    do qx = 1, line%nq
                        if (method /= 'h1') then
                            if (qx /= 1 .or. qy /= 1) cycle
                            x = -0.5_dp + (real(ex, dp) - 0.5_dp)*h
                            y = -0.5_dp + (real(ey, dp) - 0.5_dp)*h
                            weight = h*h
                        else
                            x = line%xq(qx, ex)
                            y = line%xq(qy, ey)
                            weight = line%wq(qx, ex)*line%wq(qy, ey)
                        end if
                        call benchmark(name, x, y, value, source, bx, by)
                        if (method == 'fci') then
                            points = [x, y]
                            call trace_fci_field_line_rk4(points, 0.5_dp*h, 64, &
                                pollution_rhs, endpoints(:, 1), trace_status)
                            if (trace_status%code /= 0) &
                                call gate(.false., 'FCI forward trace')
                            call trace_fci_field_line_rk4(points, -0.5_dp*h, 64, &
                                pollution_rhs, endpoints(:, 2), trace_status)
                            if (trace_status%code /= 0) &
                                call gate(.false., 'FCI backward trace')
                            gridx = [x - 0.5_dp*h, x + 0.5_dp*h]
                            gridy = [y - 0.5_dp*h, y + 0.5_dp*h]
                            call build_fci_bilinear_interpolation_map_2d(gridx, gridy, &
                                endpoints(1, :), endpoints(2, :), maps, trace_status)
                            if (trace_status%code /= 0) &
                                call gate(.false., 'FCI endpoint map')
                            do a = 1, 4
                                call generated_fci_parallel_gradient(maps(1, a), &
                                    1.0_dp, maps(2, a), 1.0_dp, h, q(a))
                            end do
                        end if
                        a = 0
                        do iy = 1, p + 1
                            do ix = 1, p + 1
                                a = a + 1
                                if (method /= 'h1') then
                                    call pollution_basis(0.5_dp, 0.5_dp, &
                                        real(2*ix - 3, dp)/h, &
                                        real(2*iy - 3, dp)/h, v(a), gx(a), gy(a))
                                else
                                    call pollution_basis(line%b(qx, ix, ex), &
                                        line%b(qy, iy, ey), line%d(qx, ix, ex), &
                                        line%d(qy, iy, ey), v(a), gx(a), gy(a))
                                end if
                            end do
                        end do
                        do a = 1, size(ids)
                            if (ids(a) == 0) cycle
                            do b = 1, size(ids)
                                call pollution_integrand(gx(a), gy(a), gx(b), &
                                    gy(b), bx, by, weight, source, v(a), &
                                    slow, parallel, load)
                                if (method /= 'h1') slow = 0.0_dp
                                if (method == 'fci') then
                                    call pollution_integrand(q(a), 0.0_dp, q(b), &
                                        0.0_dp, 1.0_dp, 0.0_dp, weight, 0.0_dp, &
                                        0.0_dp, slow, parallel, load)
                                    slow = 0.0_dp
                                end if
                                if (ids(b) == 0) then
                                    data%boundary_slow(ids(a)) = &
                                        data%boundary_slow(ids(a)) - slow*boundary(b)
                                    data%boundary_parallel(ids(a)) = &
                                        data%boundary_parallel(ids(a)) - &
                                        parallel*boundary(b)
                                else if (ids(a) >= ids(b)) then
                                    call add_entry(data, ids(a), ids(b), slow, parallel)
                                end if
                            end do
                            if (method == 'h1') &
                                data%load(ids(a)) = data%load(ids(a)) + load
                        end do
                    end do
                end do
            end do
        end do
        if (method /= 'h1') then
            do gj = 1, n
                do gi = 1, n
                    a = (gj - 1)*n + gi
                    x = -0.5_dp + real(gi, dp)*h
                    y = -0.5_dp + real(gj, dp)*h
                    call benchmark(name, x, y, value, source, bx, by)
                    data%load(a) = h*h*source
                    call add_entry(data, a, a, 4.0_dp, 0.0_dp)
                    if (gi < n) call add_entry(data, a + 1, a, -1.0_dp, 0.0_dp)
                    if (gj < n) call add_entry(data, a + n, a, -1.0_dp, 0.0_dp)
                    do b = 1, 4
                        ix = gi
                        iy = gj
                        if (b == 1) ix = ix - 1
                        if (b == 2) ix = ix + 1
                        if (b == 3) iy = iy - 1
                        if (b == 4) iy = iy + 1
                        if (ix > 0 .and. ix <= n .and. iy > 0 .and. iy <= n) cycle
                        call benchmark(name, -0.5_dp + real(ix, dp)*h, &
                            -0.5_dp + real(iy, dp)*h, value, source, bx, by)
                        data%boundary_slow(a) = data%boundary_slow(a) + value
                    end do
                end do
            end do
        end if
        call system_clock(finish)
        data%setup_seconds = real(finish - start, dp)/real(rate, dp)
    end subroutine

    pure subroutine pollution_rhs(parameter, point, derivative)
        real(dp), intent(in) :: parameter, point(:)
        real(dp), intent(out) :: derivative(:)
        real(dp) :: value, source, bx, by
        associate (unused => parameter)
            call benchmark(trim(trace_benchmark), point(1), point(2), &
                value, source, bx, by)
            derivative = [bx, by]
        end associate
    end subroutine

    subroutine sparse_pattern(p, ne, data)
        integer, intent(in) :: p, ne
        type(fixture_t), intent(inout) :: data
        integer :: n, j, x, y, lo, hi, bottom, top, k, ix, iy
        n = p*ne - 1
        allocate (data%col_ptr(n*n + 1))
        data%col_ptr(1) = 1
        do j = 1, n*n
            x = mod(j - 1, n) + 1
            y = (j - 1)/n + 1
            lo = max(1, ((x - 1)/p)*p)
            hi = min(n, (x/p + 1)*p)
            bottom = max(1, ((y - 1)/p)*p)
            top = min(n, (y/p + 1)*p)
            data%col_ptr(j + 1) = data%col_ptr(j) + (hi - lo + 1)*(top - bottom + 1)
        end do
        k = data%col_ptr(n*n + 1) - 1
        allocate (data%row_idx(k), data%sparse_slow(k), data%sparse_parallel(k))
        data%sparse_slow = 0.0_dp
        data%sparse_parallel = 0.0_dp
        do j = 1, n*n
            x = mod(j - 1, n) + 1
            y = (j - 1)/n + 1
            k = data%col_ptr(j)
            do iy = max(1, ((y - 1)/p)*p), min(n, (y/p + 1)*p)
                do ix = max(1, ((x - 1)/p)*p), min(n, (x/p + 1)*p)
                    data%row_idx(k) = (iy - 1)*n + ix
                    k = k + 1
                end do
            end do
        end do
    end subroutine

    integer function sparse_slot(data, row, col) result(k)
        type(fixture_t), intent(in) :: data
        integer, intent(in) :: row, col
        integer :: n, p, x, y, lo, hi, bottom
        p = data%degree
        n = p*data%cells - 1
        x = mod(col - 1, n) + 1
        y = (col - 1)/n + 1
        lo = max(1, ((x - 1)/p)*p)
        hi = min(n, (x/p + 1)*p)
        bottom = max(1, ((y - 1)/p)*p)
        k = data%col_ptr(col) + mod(row - 1, n) + 1 - lo + &
            ((row - 1)/n + 1 - bottom)*(hi - lo + 1)
    end function

    subroutine add_entry(data, row, col, slow, parallel)
        type(fixture_t), intent(inout) :: data
        integer, intent(in) :: row, col
        real(dp), intent(in) :: slow, parallel
        integer :: k
        if (allocated(data%col_ptr)) then
            k = sparse_slot(data, row, col)
            data%sparse_slow(k) = data%sparse_slow(k) + slow
            data%sparse_parallel(k) = data%sparse_parallel(k) + parallel
            if (row == col) return
            k = sparse_slot(data, col, row)
            data%sparse_slow(k) = data%sparse_slow(k) + slow
            data%sparse_parallel(k) = data%sparse_parallel(k) + parallel
        else
            k = row - col + 1
            data%slow(k, col) = data%slow(k, col) + slow
            data%parallel(k, col) = data%parallel(k, col) + parallel
        end if
    end subroutine

    subroutine solve_sparse(data, chi, centre, pollution, residual, seconds, solution)
        type(fixture_t), intent(in) :: data
        real(dp), intent(in) :: chi
        real(dp), intent(out) :: centre, pollution, residual, seconds
        real(dp), allocatable, intent(out), optional :: solution(:)
        type(sparse_direct_factor_t) :: factor
        type(sparse_options_t) :: options
        real(dp), allocatable :: values(:), rhs(:), x(:), r(:), rowsum(:)
        integer :: n, status, col, k, row
        integer(int64) :: start, finish, rate
        n = size(data%load)
        values = data%sparse_slow + (chi - 1.0_dp)*data%sparse_parallel
        rhs = data%load + data%boundary_slow + (chi - 1.0_dp)*data%boundary_parallel
        allocate (x(n), rowsum(n))
        options%symmetry = FORTSPARSE_SYM_SPD
        call system_clock(start, rate)
        call sparse_direct_factor_csc(factor, n, data%col_ptr, data%row_idx, &
            values, status, options)
        if (status /= 0) write (*, '(a,i0)') 'native sparse factor status=', status
        call gate(status == 0, 'native sparse SPD factorization')
        call sparse_direct_solve_factored(factor, rhs, x, status)
        call gate(status == 0, 'native sparse SPD solve')
        call system_clock(finish)
        seconds = real(finish - start, dp)/real(rate, dp)
        call sparse_direct_free(factor)
        r = -rhs
        rowsum = 0.0_dp
        do col = 1, n
            do k = data%col_ptr(col), data%col_ptr(col + 1) - 1
                row = data%row_idx(k)
                r(row) = r(row) + values(k)*x(col)
                rowsum(row) = rowsum(row) + abs(values(k))
            end do
        end do
        residual = maxval(abs(r))/(maxval(rowsum)*maxval(abs(x)) + maxval(abs(rhs)))
        centre = probe_centre(data, x)
        pollution = abs(1.0_dp/centre - 1.0_dp)
        if (present(solution)) solution = x
    end subroutine

    function probe_centre(data, solution) result(value)
        type(fixture_t), intent(in) :: data
        real(dp), intent(in) :: solution(:)
        real(dp) :: value, basis, gx, gy
        real(dp), allocatable :: nodes(:), weights(:), b(:), d(:)
        integer :: p, ne, n, ix, iy, gi, gj, status
        p = data%degree
        ne = data%cells
        if (mod(p*ne, 2) == 0) then
            value = solution(data%centre)
            return
        end if
        call gll_rule(p, nodes, weights, status)
        if (status /= 0) error stop 'central interpolation nodes failed'
        allocate (b(p + 1), d(p + 1))
        call lagrange_jet(nodes, 0.0_dp, b, d)
        n = p*ne - 1
        value = 0.0_dp
        do iy = 1, p + 1
            gj = (ne/2)*p + iy - 1
            do ix = 1, p + 1
                gi = (ne/2)*p + ix - 1
                call pollution_basis(b(ix), b(iy), d(ix), d(iy), basis, gx, gy)
                value = value + basis*solution((gj - 1)*n + gi)
            end do
        end do
    end function

    subroutine sparse_oracle()
        type(fixture_t) :: band, sparse
        real(dp), allocatable :: xb(:), xs(:)
        real(dp) :: c, e, r, s, chi
        integer :: p, power
        do p = 1, 4
            call assemble_fixture('sovinec', 'h1', p, 4, band)
            call assemble_fixture('sovinec', 'h1', p, 4, sparse, .true.)
            do power = 0, 10, 5
                chi = 10.0_dp**power
                call solve_fixture(band, chi, c, e, r, s, xb)
                call gate(r <= residual_limit, 'banded oracle residual')
                call solve_fixture(sparse, chi, c, e, r, s, xs)
                call gate(r <= residual_limit, 'sparse oracle residual')
                call gate(maxval(abs(xs - xb)) <= 1.0e-6_dp*maxval(abs(xb)), &
                    'sparse solution against independent LAPACK banded solve')
            end do
        end do
        call assemble_fixture('sovinec', 'h1', 1, 31, sparse, .true.)
        call solve_fixture(sparse, 1.0_dp, c, e, r, s)
        call gate(abs(c - 1.0_dp) < 0.002_dp, &
            'between-node central value against analytic cosine solution')
        call gate(r <= residual_limit, 'between-node central-value residual')
    end subroutine

    subroutine cost_ladder()
        type(fixture_t) :: data
        character(4096) :: name, method, argument, path
        real(dp) :: c, e, r, s, baseline, baseline_s, chi, signed_error
        integer :: p, ne, out, power
        call get_command_argument(2, name)
        call get_command_argument(3, method)
        call get_command_argument(4, argument)
        read (argument, *) p
        call get_command_argument(5, argument)
        read (argument, *) power
        chi = 10.0_dp**power
        call get_command_argument(6, path)
        open (newunit=out, file=trim(path), status='replace')
        write (out, '(a)') 'benchmark,method,p,cells,unknowns,chi_ratio,'// &
            'centre,effective_perp_error,residual,setup_s,solve_s,'// &
            'solve_s_per_unknown,baseline_effective_perp,signed_pollution'
        ne = 4
        do
            call assemble_fixture(trim(name), trim(method), p, ne, data, .true.)
            call solve_fixture(data, 1.0_dp, c, e, r, baseline_s)
            baseline = 1.0_dp/c
            call gate(r <= residual_limit, 'cost baseline normwise residual')
            call solve_fixture(data, chi, c, e, r, s)
            signed_error = 1.0_dp/c - baseline
            call write_row(out, trim(name), trim(method), data, chi, c, e, r, s, &
                baseline, signed_error)
            flush (out)
            call gate(r <= residual_limit, 'cost anisotropic normwise residual')
            if (signed_error <= 0.1_dp) exit
            if (ne == 1001/p) exit
            ne = min(2*ne, 1001/p)
        end do
        close (out)
    end subroutine

    subroutine solve_fixture(data, chi, centre, pollution, residual, seconds, solution)
        type(fixture_t), intent(in) :: data
        real(dp), intent(in) :: chi
        real(dp), intent(out) :: centre, pollution, residual, seconds
        real(dp), allocatable, intent(out), optional :: solution(:)
        real(dp), allocatable :: band(:, :), factor(:, :), rhs(:, :), original(:)
        integer :: n, bw, info
        integer(int64) :: start, finish, rate
        if (allocated(data%col_ptr)) then
            call solve_sparse(data, chi, centre, pollution, residual, seconds, solution)
            return
        end if
        n = size(data%load)
        bw = size(data%slow, 1) - 1
        band = data%slow + (chi - 1.0_dp)*data%parallel
        factor = band
        allocate (rhs(n, 1))
        original = data%load + data%boundary_slow + &
            (chi - 1.0_dp)*data%boundary_parallel
        rhs(:, 1) = original
        call system_clock(start, rate)
        call dpbsv('L', n, bw, 1, factor, bw + 1, rhs, n, info)
        call system_clock(finish)
        if (info /= 0) error stop 'banded reference factorization failed'
        seconds = real(finish - start, dp)/real(rate, dp)
        residual = backward_residual(band, rhs(:, 1), original)
        centre = probe_centre(data, rhs(:, 1))
        pollution = abs(1.0_dp/centre - 1.0_dp)
        if (present(solution)) solution = rhs(:, 1)
    end subroutine

    function backward_residual(band, x, rhs) result(value)
        real(dp), intent(in) :: band(:, :), x(:), rhs(:)
        real(dp) :: value
        real(dp), allocatable :: r(:), row_norm(:)
        integer :: n, bw, j, k, i
        n = size(rhs)
        bw = size(band, 1) - 1
        r = -rhs
        row_norm = abs(band(1, :))
        do j = 1, n
            r(j) = r(j) + band(1, j)*x(j)
            do k = 1, min(bw, n - j)
                i = j + k
                r(i) = r(i) + band(k + 1, j)*x(j)
                r(j) = r(j) + band(k + 1, j)*x(i)
                row_norm(i) = row_norm(i) + abs(band(k + 1, j))
                row_norm(j) = row_norm(j) + abs(band(k + 1, j))
            end do
        end do
        value = maxval(abs(r))/(maxval(row_norm)*maxval(abs(x)) + maxval(abs(rhs)))
    end function

    subroutine write_row(unit, name, method, data, chi, centre, pollution, residual, &
            s, baseline, signed_error)
        integer, intent(in) :: unit
        character(*), intent(in) :: name, method
        type(fixture_t), intent(in) :: data
        real(dp), intent(in) :: chi, centre, pollution, residual, s
        real(dp), intent(in), optional :: baseline, signed_error
        write (unit, '(a,",",a,3(",",i0),7(",",es24.16))', advance='no') &
            name, method, &
            data%degree, data%cells, size(data%load), chi, centre, pollution, &
            residual, data%setup_seconds, s, s/real(size(data%load), dp)
        if (present(baseline) .and. present(signed_error)) then
            write (unit, '(2(",",es24.16))') baseline, signed_error
        else
            write (unit, *)
        end if
        write (*, '(a,1x,a,2(1x,i0),3(1x,es13.5))') name, method, &
            data%degree, data%cells, chi, pollution, residual
    end subroutine
end program
