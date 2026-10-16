module meridian_edge_remainder_fixture
    use fortfem_kinds, only: dp
    implicit none
    real(dp), allocatable :: left(:), width(:)
    logical :: polynomial = .false.
contains
    subroutine load(cell, edge, tau, g)
        integer, intent(in) :: cell, edge
        real(dp), intent(in) :: tau
        complex(dp), intent(out) :: g
        real(dp) :: x
        x = left(cell) + width(cell)*tau
        if (polynomial) then
            g = cmplx(1 + x + x*x + x**3, 0.0_dp, dp)
        else
            g = cmplx(cos(acos(-1.0_dp)*x), 0.0_dp, dp)
        end if
        if (edge < 1) error stop 'invalid edge'
    end subroutine load

    subroutine load_derivative(cell, edge, tau, dg)
        integer, intent(in) :: cell, edge
        real(dp), intent(in) :: tau
        complex(dp), intent(out) :: dg
        real(dp) :: x, pi
        pi = acos(-1.0_dp)
        x = left(cell) + width(cell)*tau
        if (polynomial) then
            dg = cmplx(width(cell)*(1 + 2*x + 3*x*x), 0.0_dp, dp)
        else
            dg = cmplx(-pi*width(cell)*sin(pi*x), 0.0_dp, dp)
        end if
        if (edge < 1) error stop 'invalid edge'
    end subroutine load_derivative
end module meridian_edge_remainder_fixture

program test_meridian_edge_remainder
    !! Independent continuous oracle: a sheet in the periodic-cylinder
    !! rectangle [0,1] x [-1/2,1/2], zero tangential A on the outer boundary.
    !! For a cosine sheet coefficient r_j, the exact correction energy is
    !! R0/n^2 * alpha_j*tanh(alpha_j/2)*|r_j|^2/4 (j > 0).
    !! The Fourier coefficients of G-I_3 G are integrated analytically below;
    !! doubling the series length checks the independent reference tail.
    use meridian_edge_remainder_fixture
    use, intrinsic :: iso_fortran_env, only: int64
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build
    use fortfem_equilibrated_hcurl_meridian_n, only: &
        meridian_hcurl_edge_remainder, HCURL_OK
    use fortfem_generated_meridian_edge_lift, only: generated_meridian_edge_lift
    use fortfem_generated_meridian_edge_lift_jvp, only: &
        generated_meridian_edge_lift_jvp
    use fortfem_generated_meridian_edge_lift_vjp, only: &
        generated_meridian_edge_lift_vjp
    !$ use omp_lib, only: omp_set_num_threads
    implicit none
    type(mesh_2d_t) :: mesh
    type(meridian_mesh_t) :: mm
    type(meridian_feec_t) :: fs
    integer, allocatable :: cells(:), edges(:)
    real(dp) :: eta(3), truth, tail_check, ratio, rzero, rnegative, rdouble
    real(dp) :: xa, xb
    integer :: lev, nx, t, e, i, count_, st

    call check_leaf_derivatives()
    do lev = 1, 3
        nx = 2**lev
        call mesh%create_rectangular(nx + 1, nx + 1, &
            0.0_dp, 1.0_dp, -0.5_dp, 0.5_dp)
        call meridian_mesh_straight(mesh, 1, mm, st)
        call require(st == 0, 'mesh')
        call meridian_feec_build(mm, 3, fs, st, quad_degree=24, cylinder_r0=2.0_dp)
        call require(st == 0, 'FEEC')
        allocate(cells(nx), edges(nx), left(mesh%n_triangles), &
            width(mesh%n_triangles))
        left = 0.0_dp
        width = 0.0_dp
        count_ = 0
        do t = 1, mesh%n_triangles
            do e = 1, 3
                i = mod(e, 3) + 1
                if (abs(mesh%vertices(2, mesh%triangles(e, t))) > 1.0e-14_dp) cycle
                if (abs(mesh%vertices(2, mesh%triangles(i, t))) > 1.0e-14_dp) cycle
                xa = mesh%vertices(1, mesh%triangles(e, t))
                xb = mesh%vertices(1, mesh%triangles(i, t))
                if (xb <= xa) cycle
                count_ = count_ + 1
                cells(count_) = t
                edges(count_) = e
                left(t) = xa
                width(t) = xb - xa
            end do
        end do
        call require(count_ == nx, 'sheet edge count')
        polynomial = .true.
        call meridian_hcurl_edge_remainder(fs, 2, cells, edges, load, &
            load_derivative, rzero, st)
        call require(st == HCURL_OK, 'polynomial status')
        call require(rzero < 2.0e-12_dp, 'cubic recovery to rounding')
        print '(a,i3,es12.4)', 'cubic recovery ', nx, rzero
        polynomial = .false.
        call meridian_hcurl_edge_remainder(fs, 2, cells, edges, load, &
            load_derivative, eta(lev), st)
        call require(st == HCURL_OK, 'remainder status')
        call meridian_hcurl_edge_remainder(fs, -2, cells, edges, load, &
            load_derivative, rnegative, st)
        call require(st == HCURL_OK, 'negative mode status')
        call require(abs(rnegative - eta(lev)) < 1.0e-13_dp, 'mode parity')
        call meridian_hcurl_edge_remainder(fs, 4, cells, edges, load, &
            load_derivative, rdouble, st)
        call require(st == HCURL_OK, 'larger mode status')
        call require(rdouble < eta(lev), 'mode-dependent constant')
        truth = exact_energy(nx, 8192)
        tail_check = exact_energy(nx, 16384)
        print '(a,i3,3es16.8)', 'Fourier reference ', nx, truth, tail_check, &
            truth/tail_check - 1
        call require(abs(truth/tail_check - 1) < 1.0e-5_dp, 'Fourier tail')
        ratio = eta(lev)/tail_check
        print '(a,i3,3(a,es12.4))', 'nx=', nx, ' eta=', eta(lev), &
            ' error=', tail_check, ' index=', ratio
        call require(eta(lev) >= tail_check, 'continuous error containment')
        call require(ratio < 3.0_dp, 'efficiency cap')
        call meridian_hcurl_edge_remainder(fs, 0, cells, edges, load, &
            load_derivative, rzero, st)
        call require(st /= HCURL_OK, 'zero mode rejected')
        call meridian_hcurl_edge_remainder(fs, 2, [cells(1), cells(1)], &
            [edges(1), edges(1)], load, load_derivative, rzero, st)
        call require(st /= HCURL_OK, 'duplicate oriented edge rejected')
        if (lev == 3) call scaling()
        deallocate(cells, edges, left, width)
    end do
    print '(a,2f8.3)', 'orders=', log(eta(:2)/eta(2:))/log(2.0_dp)
    call require(all(eta(2:) < eta(:2)/8), 'at least cubic convergence')
    print '(a)', 'PASS'
contains
    subroutine check_leaf_derivatives()
        real(dp) :: v(8), d(8), vp(8), vm(8), fp(3), fm(3), fd(3), adj(8)
        real(dp), parameter :: step = 1.0e-6_dp, bar(3) = [0.7_dp,-0.3_dp,0.2_dp]
        integer :: j
        v = [0.3_dp, 0.5_dp, -1.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, 0.2_dp, -0.4_dp]
        d = [(sin(real(j, dp)), j=1, 8)]
        vp = v + step*d
        vm = v - step*d
        call generated_meridian_edge_lift(vp(1),vp(2),vp(3),vp(4),vp(5),vp(6), &
            vp(7),vp(8),fp(1),fp(2),fp(3))
        call generated_meridian_edge_lift(vm(1),vm(2),vm(3),vm(4),vm(5),vm(6), &
            vm(7),vm(8),fm(1),fm(2),fm(3))
        call generated_meridian_edge_lift_jvp(v(1),v(2),v(3),v(4),v(5),v(6), &
            v(7),v(8),d(1),d(2),d(3),d(4),d(5),d(6),d(7),d(8),fd(1),fd(2),fd(3))
        call require(maxval(abs(fd - (fp - fm)/(2*step))) < 1.0e-9_dp, 'lift JVP')
        call generated_meridian_edge_lift_vjp(v(1),v(2),v(3),v(4),v(5),v(6), &
            v(7),v(8),bar(1),bar(2),bar(3),adj(1),adj(2),adj(3),adj(4),adj(5), &
            adj(6),adj(7),adj(8))
        call require(abs(dot_product(bar,fd) - dot_product(adj,d)) < 1.0e-13_dp, &
            'lift JVP/VJP pairing')
    end subroutine check_leaf_derivatives

    subroutine scaling()
        integer :: threads, rep, iteration, st
        integer(int64) :: tick, tock, ticks_per_second
        real(dp) :: seconds(3), answer
        do threads = 1, 4
            if (threads == 3) cycle
            !$ call omp_set_num_threads(threads)
            do rep = 1, 3
                call system_clock(tick, ticks_per_second)
                do iteration = 1, 20
                    call meridian_hcurl_edge_remainder(fs, 2, cells, edges, load, &
                        load_derivative, answer, st)
                    call require(st == HCURL_OK, 'parallel remainder status')
                    call require(abs(answer/eta(3) - 1) < 1.0e-12_dp, &
                        'thread count preserves remainder')
                end do
                call system_clock(tock)
                seconds(rep) = real(tock - tick, dp)/ticks_per_second/20
            end do
            print '(a,i2,3es12.4)', 'seconds/call at threads=', threads, seconds
        end do
    end subroutine scaling

    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: message
        if (.not. ok) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine require

    real(dp) function exact_energy(nx, cutoff) result(norm_)
        integer, intent(in) :: nx, cutoff
        real(dp) :: pi, h, a, alpha, coeff, sum_, g(0:3), c(0:3)
        complex(dp) :: iw, ee, moment(0:3), integral_
        integer :: j, k, panel
        pi = acos(-1.0_dp)
        h = 1.0_dp/nx
        sum_ = 0.0_dp
        do j = 1, cutoff
            iw = cmplx(0.0_dp, pi*j*h, dp)
            ee = exp(iw)
            moment(0) = (ee - 1)/iw
            do k = 1, 3
                moment(k) = (ee - k*moment(k - 1))/iw
            end do
            integral_ = (0.0_dp, 0.0_dp)
            do panel = 0, nx - 1
                a = panel*h
                g = [(cos(pi*(a + h*k/3)), k=0, 3)]
                c(0) = g(0)
                c(1) = (-11*g(0) + 18*g(1) - 9*g(2) + 2*g(3))/2
                c(2) = (18*g(0) - 45*g(1) + 36*g(2) - 9*g(3))/2
                c(3) = (-9*g(0) + 27*g(1) - 27*g(2) + 9*g(3))/2
                integral_ = integral_ + h*exp(cmplx(0.0_dp, pi*j*a, dp)) &
                    *sum(c*moment)
            end do
            coeff = -2*real(integral_, dp)
            if (j == 1) coeff = coeff + 1
            alpha = sqrt((pi*j)**2 + 1.0_dp)
            sum_ = sum_ + coeff**2*alpha*tanh(alpha/2)/8
        end do
        norm_ = sqrt(sum_)
    end function exact_energy
end program test_meridian_edge_remainder
