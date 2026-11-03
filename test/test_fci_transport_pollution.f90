program test_fci_transport_pollution
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use check, only: check_condition, check_summary
    use fortfem_gll_line, only: build_gll_line
    use fortfem_tensor_line_tables, only: tensor_line_tables_t
    use fortfem_generated_pollution_sovinec, only: pollution_sovinec
    use fortfem_generated_pollution_circle, only: pollution_circle
    use fortfem_generated_pollution_basis, only: pollution_basis
    use fortfem_generated_pollution_integrand, only: pollution_integrand
    implicit none

    type :: fixture_t
        real(dp), allocatable :: slow(:, :), parallel(:, :), load(:)
        real(dp), allocatable :: boundary_slow(:), boundary_parallel(:)
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
    real(dp) :: centre, pollution, residual, seconds, expected, pi, h
    integer :: cells, unit, exponent

    pi = acos(-1.0_dp)
    call get_command_argument(1, mode)
    if (trim(mode) == 'paper') then
        ! Figure-reading intervals are fixed before running the harness.
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
        call check_condition(residual < 1.0e-10_dp, 'original weak residual')
        call check_condition(pollution >= 5.0_dp .and. pollution <= 7.0_dp, &
            'Sovinec Fig. 9: p=1, h=1/8, chi=1e3 graphical interval [5,7]')
    else
        do cells = 8, 16, 8
            call assemble_fixture('sovinec', 'symmetric', 1, cells, fixture)
            h = 1.0_dp/real(cells, dp)
            expected = (pi*h/(2.0_dp*sin(pi*h/2.0_dp)))**2
            do exponent = 0, 6, 3
                call solve_fixture(fixture, 10.0_dp**exponent, &
                    centre, pollution, residual, seconds)
                write (*, '(a,i0,a,i0,2(a,es16.8))') 'sine cells=', cells, &
                    ' exponent=', exponent, ' centre=', centre, ' residual=', residual
                call check_condition(abs(centre - expected) < &
                    2.0e-9_dp*expected, 'symmetric sine eigenmode oracle')
                call check_condition(residual < 1.0e-9_dp, &
                    'symmetric original residual')
            end do
        end do
        call assemble_fixture('sovinec', 'h1', 1, 8, fixture)
        call solve_fixture(fixture, 1.0_dp, centre, pollution, residual, seconds)
        call check_condition(abs(centre - 1.0_dp) < 0.02_dp, &
            'H1 isotropic central value against analytic solution')
        call check_condition(residual < 1.0e-12_dp, 'H1 original residual')
    end if
    call check_summary('anisotropic pollution')
contains
    subroutine benchmark(name, x, y, value, source, bx, by)
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

    subroutine assemble_fixture(name, method, p, ne, data)
        character(*), intent(in) :: name, method
        integer, intent(in) :: p, ne
        type(fixture_t), intent(out) :: data
        type(tensor_line_tables_t) :: line
        real(dp), allocatable :: v(:), gx(:), gy(:), boundary(:)
        integer, allocatable :: ids(:)
        integer :: n, nd, bw, ex, ey, ix, iy, qx, qy, a, b, gi, gj, status
        integer(int64) :: start, finish, rate
        real(dp) :: value, source, bx, by, slow, parallel, load, weight, h, x, y
        call system_clock(start, rate)
        n = p*ne - 1
        nd = n*n
        bw = min(p*n + p, nd - 1)
        data%cells = ne
        data%degree = p
        data%centre = (n/2)*n + n/2 + 1
        allocate (data%slow(bw + 1, nd), data%parallel(bw + 1, nd))
        allocate (data%load(nd), data%boundary_slow(nd), &
            data%boundary_parallel(nd))
        data%slow = 0.0_dp
        data%parallel = 0.0_dp
        data%load = 0.0_dp
        data%boundary_slow = 0.0_dp
        data%boundary_parallel = 0.0_dp
        allocate (v((p + 1)**2), gx((p + 1)**2), gy((p + 1)**2))
        allocate (ids((p + 1)**2), boundary((p + 1)**2))
        call build_gll_line(p, ne, line, status, x0=-0.5_dp, x1=0.5_dp)
        if (status /= 0) error stop 'line construction failed'
        h = 1.0_dp/real(ne, dp)
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
                            ! Boundary values use physical nodes; p=1 for circle.
                            x = -0.5_dp + real(gi, dp)/real(p*ne, dp)
                            y = -0.5_dp + real(gj, dp)/real(p*ne, dp)
                            call benchmark(name, x, y, boundary(a), &
                                source, bx, by)
                        end if
                    end do
                end do
                do qy = 1, line%nq
                    do qx = 1, line%nq
                        if (method == 'symmetric') then
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
                        a = 0
                        do iy = 1, p + 1
                            do ix = 1, p + 1
                                a = a + 1
                                if (method == 'symmetric') then
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
                                if (method == 'symmetric') slow = 0.0_dp
                                if (ids(b) == 0) then
                                    data%boundary_slow(ids(a)) = &
                                        data%boundary_slow(ids(a)) - slow*boundary(b)
                                    data%boundary_parallel(ids(a)) = &
                                        data%boundary_parallel(ids(a)) - &
                                        parallel*boundary(b)
                                else if (ids(a) >= ids(b)) then
                                    gi = ids(a) - ids(b) + 1
                                    data%slow(gi, ids(b)) = &
                                        data%slow(gi, ids(b)) + slow
                                    data%parallel(gi, ids(b)) = &
                                        data%parallel(gi, ids(b)) + parallel
                                end if
                            end do
                            if (method /= 'symmetric') &
                                data%load(ids(a)) = data%load(ids(a)) + load
                        end do
                    end do
                end do
            end do
        end do
        if (method == 'symmetric') then
            do gj = 1, n
                do gi = 1, n
                    a = (gj - 1)*n + gi
                    x = -0.5_dp + real(gi, dp)*h
                    y = -0.5_dp + real(gj, dp)*h
                    call benchmark(name, x, y, value, source, bx, by)
                    data%load(a) = h*h*source
                    data%slow(1, a) = 4.0_dp
                    if (gi < n) data%slow(2, a) = -1.0_dp
                    if (gj < n) data%slow(n + 1, a) = -1.0_dp
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

    subroutine solve_fixture(data, chi, centre, pollution, residual, seconds)
        type(fixture_t), intent(in) :: data
        real(dp), intent(in) :: chi
        real(dp), intent(out) :: centre, pollution, residual, seconds
        real(dp), allocatable :: band(:, :), factor(:, :), rhs(:, :), original(:), r(:)
        integer :: n, bw, info, j, k, i
        integer(int64) :: start, finish, rate
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
        r = -original
        do j = 1, n
            r(j) = r(j) + band(1, j)*rhs(j, 1)
            do k = 1, min(bw, n - j)
                i = j + k
                r(i) = r(i) + band(k + 1, j)*rhs(j, 1)
                r(j) = r(j) + band(k + 1, j)*rhs(i, 1)
            end do
        end do
        residual = norm2(r)/norm2(original)
        centre = rhs(data%centre, 1)
        pollution = abs(1.0_dp/centre - 1.0_dp)
    end subroutine

    subroutine write_row(unit, name, method, data, chi, centre, pollution, residual, s)
        integer, intent(in) :: unit
        character(*), intent(in) :: name, method
        type(fixture_t), intent(in) :: data
        real(dp), intent(in) :: chi, centre, pollution, residual, s
        write (unit, '(a,",",a,3(",",i0),7(",",es24.16))') name, method, &
            data%degree, data%cells, size(data%load), chi, centre, pollution, &
            residual, data%setup_seconds, s, s/real(size(data%load), dp)
        write (*, '(a,1x,a,2(1x,i0),3(1x,es13.5))') name, method, &
            data%degree, data%cells, chi, pollution, residual
    end subroutine
end program
