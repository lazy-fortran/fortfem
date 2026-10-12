program test_axisymmetric_gs_boundary_operators_2d
    !! Independent oracles for the axisymmetric Grad-Shafranov boundary
    !! operators on curved P2 panels:
    !!  - the constant-density single layer against a ring of filaments
    !!    summed with AGM elliptic integrals,
    !!  - K 1 = -1/2, K' = K^T, symmetry/positivity of V, W and S,
    !!  - the Galerkin Calderon identity V W = 1/4 - K^2 against the
    !!    Maue-type hypersingular form,
    !!  - the exterior Dirichlet problem for the field of an interior
    !!    filament on a circle and a D-shaped curve: Neumann data, energy
    !!    and observed convergence orders,
    !!  - consistency of the Johnson-Nedelec form, the Neumann-to-Dirichlet
    !!    map.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_axisymmetric_gs_kernel, only: gs_green, gs_green_gradient
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: &
        gs_bem_curve_t, gs_bem_operators_t, gs_bem_curve_from_nodes, &
        gs_bem_assemble, gs_bem_steklov_symmetric, gs_bem_steklov_jn, &
        gs_bem_neumann_to_dirichlet, gs_bem_dual_norm_squared, &
        gs_bem_exterior_neumann, gs_bem_panel_point, gs_bem_trace_value, &
        gs_bem_flux_value, gs_bem_trace_node_index
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: circle_center(2) = [3.0_dp, 0.5_dp]
    real(dp), parameter :: circle_source(2) = [3.25_dp, 0.62_dp]
    real(dp), parameter :: d_source(2) = [3.3_dp, 0.3_dp]
    logical :: all_passed

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character(len=1), intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine dsyev
    end interface

    all_passed = .true.

    call test_uniform_sheet()
    call test_structure(1)
    call test_structure(2)
    call test_exterior_dirichlet(1, 0)
    call test_exterior_dirichlet(1, 1)
    call test_exterior_dirichlet(2, 1)
    call test_exterior_dirichlet(2, 2)
    call test_calderon(1)
    call test_calderon(2)
    call test_jn_and_ntd()

    call check_summary("Axisymmetric Grad-Shafranov boundary operators")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

    ! ------------------------------------------------------------------
    ! Curves
    ! ------------------------------------------------------------------

    pure function curve_point(shape, theta) result(x)
        !! 1: circle of radius 1 about (3, 0.5); 2: D shape
        !! R = 3 + cos(theta + 0.4 sin theta), Z = 0.5 + 1.4 sin(theta).
        integer, intent(in) :: shape
        real(dp), intent(in) :: theta
        real(dp) :: x(2)

        if (shape == 1) then
            x = circle_center + [cos(theta), sin(theta)]
        else
            x = [3.0_dp + cos(theta + 0.4_dp*sin(theta)), &
                0.5_dp + 1.4_dp*sin(theta)]
        end if
    end function curve_point

    subroutine make_curve(shape, np, curve)
        integer, intent(in) :: shape, np
        type(gs_bem_curve_t), intent(out) :: curve
        real(dp), allocatable :: nodes(:, :)
        integer :: j, status

        allocate(nodes(2, 2*np))
        do j = 1, 2*np
            nodes(:, j) = curve_point(shape, pi*real(j - 1, dp)/np)
        end do
        call gs_bem_curve_from_nodes(nodes, curve, status)
        if (status /= 0) error stop "curve construction failed"
    end subroutine make_curve

    function source_point(shape) result(y)
        integer, intent(in) :: shape
        real(dp) :: y(2)

        if (shape == 1) then
            y = circle_source
        else
            y = d_source
        end if
    end function source_point

    ! ------------------------------------------------------------------
    ! Independent kernel: AGM elliptic integrals and exact differences
    ! ------------------------------------------------------------------

    function green_agm(r, rp, dd2, dz) result(g)
        !! G = (1/(2 pi)) sqrt(R R') ((2 - m) K - 2 E)/sqrt(m) from the
        !! arithmetic-geometric mean, with q = |x - y|^2/s^2 given exactly.
        real(dp), intent(in) :: r, rp, dd2, dz
        real(dp) :: g, s2, q, m, a, b, c, kk, ee, sumc, pw
        integer :: i

        s2 = (r + rp)**2 + dz*dz
        q = dd2/s2
        m = 4.0_dp*r*rp/s2
        a = 1.0_dp
        b = sqrt(q)
        c = sqrt(m)
        sumc = 0.5_dp*m
        pw = 0.5_dp
        do i = 1, 60
            c = 0.5_dp*(a - b)
            pw = pw*2.0_dp
            sumc = sumc + pw*c*c
            call agm_step(a, b)
            if (abs(c) < 1.0e-17_dp) exit
        end do
        kk = 0.5_dp*pi/a
        ee = kk*(1.0_dp - sumc)
        g = sqrt(r*rp)*((2.0_dp - m)*kk - 2.0_dp*ee)/(2.0_dp*pi*sqrt(m))
    end function green_agm

    pure subroutine agm_step(a, b)
        real(dp), intent(inout) :: a, b
        real(dp) :: an

        an = 0.5_dp*(a + b)
        b = sqrt(a*b)
        a = an
    end subroutine agm_step

    function sheet_potential(phi_x) result(val)
        !! Potential at the circle point of angle phi_x of a uniform unit
        !! line density on the unit circle about circle_center, summed as
        !! a ring of filaments: int G ds over the circle.
        real(dp), intent(in) :: phi_x
        real(dp) :: val
        real(dp), parameter :: sigma = 0.2_dp
        real(dp) :: gx(10), gw(10), lo, hi, t, f, x(2), y(2), dd2
        real(dp) :: ctr
        integer :: layer, k, sgn

        call gauss_unit(10, gx, gw)
        x = circle_center + [cos(phi_x), sin(phi_x)]
        val = 0.0_dp
        hi = pi
        do layer = 1, 24
            lo = hi*sigma
            if (layer == 24) lo = 0.0_dp
            do k = 1, 10
                t = lo + (hi - lo)*gx(k)
                f = 0.0_dp
                do sgn = -1, 1, 2
                    ctr = phi_x + sgn*t
                    y = circle_center + [cos(ctr), sin(ctr)]
                    dd2 = 4.0_dp*sin(0.5_dp*t)**2
                    f = f + green_agm(x(1), y(1), dd2, x(2) - y(2))
                end do
                val = val + (hi - lo)*gw(k)*f
            end do
            hi = lo
        end do
    end function sheet_potential

    subroutine gauss_unit(n, x, w)
        use fortnum_quadrature, only: gauss_legendre
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)

        call gauss_legendre(n, x, w)
        x = 0.5_dp*(x + 1.0_dp)
        w = 0.5_dp*w
    end subroutine gauss_unit

    ! ------------------------------------------------------------------
    ! Exact exterior data of the filament
    ! ------------------------------------------------------------------

    subroutine filament_data(src, x, tau, u, t)
        !! Trace u = G(x, src) and t = (1/R) n . grad G, n = (tau_Z, -tau_R).
        real(dp), intent(in) :: src(2), x(2), tau(2)
        real(dp), intent(out) :: u, t
        real(dp) :: grad(4)

        call gs_green_gradient(x(1), x(2), src(1), src(2), u, grad)
        t = (tau(2)*grad(1) - tau(1)*grad(2))/x(1)
    end subroutine filament_data

    subroutine nodal_trace(curve, src, u)
        type(gs_bem_curve_t), intent(in) :: curve
        real(dp), intent(in) :: src(2)
        real(dp), allocatable, intent(out) :: u(:)
        real(dp) :: x(2), tau(2), sp, t, xi
        integer :: p, k

        allocate(u(2*curve%n_panels))
        do p = 1, curve%n_panels
            do k = 1, 3
                xi = 0.5_dp*(k - 1)
                call gs_bem_panel_point(curve, p, xi, x, tau, sp)
                call filament_data(src, x, tau, u(gs_bem_trace_node_index( &
                    curve, p, k)), t)
            end do
        end do
    end subroutine nodal_trace

    ! ------------------------------------------------------------------
    ! Tests
    ! ------------------------------------------------------------------

    subroutine test_uniform_sheet()
        !! V applied to the constant density against the circle sheet summed
        !! from filaments.  The P2 curve approximates the circle with an
        !! O(h^3) position error, so the comparison converges with the mesh.
        integer, parameter :: nq = 8, levels(3) = [12, 24, 48]
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ops
        real(dp), allocatable :: ones_y(:), vones(:), oracle(:)
        real(dp) :: gx(nq), gw(nq), x(2), tau(2), sp, phi(2), err(3, 0:1), scale
        integer :: np, p, k, deg, status, l, ii, i

        call gauss_unit(nq, gx, gw)
        do i = 1, 3
            np = levels(i)
            call make_curve(1, np, curve)
            do deg = 0, 1
                call gs_bem_assemble(curve, deg, ops, status, &
                    hypersingular=.false.)
                allocate(ones_y(ops%ny), oracle(ops%ny))
                ones_y = 1.0_dp
                vones = matmul(ops%v, ones_y)
                oracle = 0.0_dp
                do p = 1, np
                    do k = 1, nq
                        call gs_bem_panel_point(curve, p, gx(k), x, tau, sp)
                        scale = gw(k)*sp*sheet_potential(atan2(x(2) - &
                            circle_center(2), x(1) - circle_center(1)))
                        if (deg == 0) then
                            phi(1) = 1.0_dp
                        else
                            phi = [1.0_dp - gx(k), gx(k)]
                        end if
                        do l = 1, deg + 1
                            ii = (p - 1)*(deg + 1) + l
                            oracle(ii) = oracle(ii) + scale*phi(l)
                        end do
                    end do
                end do
                err(i, deg) = maxval(abs(vones - oracle))/maxval(abs(oracle))
                deallocate(ones_y, oracle)
            end do
        end do
        do deg = 0, 1
            print '(a,i0,3es11.3,a,f6.2)', "  sheet potential, flux degree ", &
                deg, err(:, deg), " rate", log(err(2, deg)/err(3, deg))/log(2.0_dp)
            call record_condition(err(3, deg) < 1.0e-6_dp .and. &
                log(err(2, deg)/err(3, deg))/log(2.0_dp) > 3.5_dp, &
                "V 1 converges to the filament-ring sheet potential")
        end do
    end subroutine test_uniform_sheet

    subroutine test_structure(shape)
        integer, intent(in) :: shape
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ops
        real(dp), allocatable :: ones_x(:), r(:), s(:, :), ev(:)
        real(dp) :: wnorm, scale
        integer :: status
        character(len=24) :: tag

        write (tag, '(a,i0)') "shape ", shape
        call make_curve(shape, 32, curve)
        call gs_bem_assemble(curve, 1, ops, status)
        call record_condition(status == 0, "assembly status, "//trim(tag))
        allocate(ones_x(ops%nx))
        ones_x = 1.0_dp
        r = matmul(ops%k, ones_x) + 0.5_dp*matmul(ops%m, ones_x)
        call record_condition(maxval(abs(r)) < 1.0e-8_dp, &
            "K 1 = -1/2 (interior jump), "//trim(tag))
        call record_condition(maxval(abs(ops%kt - transpose(ops%k))) &
            < 1.0e-12_dp, "K' assembled from d_nx G equals K^T, "//trim(tag))
        scale = maxval(abs(ops%v))
        call record_condition(maxval(abs(ops%v - transpose(ops%v))) &
            < 1.0e-12_dp*scale, "V symmetric, "//trim(tag))
        wnorm = maxval(abs(ops%w))
        call record_condition(maxval(abs(ops%w - transpose(ops%w))) &
            < 1.0e-6_dp*wnorm, "W symmetric (Maue-type form), "//trim(tag))
        call record_condition(maxval(abs(matmul(ops%w, ones_x))) &
            < 1.0e-10_dp*wnorm, "W 1 = 0, "//trim(tag))
        ev = eigenvalues(ops%v)
        call record_condition(ev(1) > 0.0_dp, "V positive definite, "//trim(tag))
        ev = eigenvalues(0.5_dp*(ops%w + transpose(ops%w)))
        call record_condition(ev(1) > -1.0e-7_dp*ev(size(ev)), &
            "W positive semidefinite, "//trim(tag))
        call gs_bem_steklov_symmetric(ops, s, status)
        ev = eigenvalues(s)
        call record_condition(status == 0 .and. ev(1) > 0.0_dp, &
            "S_h symmetric positive definite, "//trim(tag))
    end subroutine test_structure

    function eigenvalues(a) result(w)
        real(dp), intent(in) :: a(:, :)
        real(dp), allocatable :: w(:)
        real(dp), allocatable :: b(:, :), work(:)
        integer :: n, info

        n = size(a, 1)
        b = a
        allocate(w(n), work(max(1, 3*n)))
        call dsyev('N', 'U', n, b, n, w, work, size(work), info)
        if (info /= 0) w = -huge(1.0_dp)
    end function eigenvalues

    subroutine neumann_l2_error(curve, ops, src, tcoef, err)
        !! L2(Gamma) error of -t_h(flux) against t of the filament field.
        type(gs_bem_curve_t), intent(in) :: curve
        type(gs_bem_operators_t), intent(in) :: ops
        real(dp), intent(in) :: src(2), tcoef(:)
        real(dp), intent(out) :: err
        integer, parameter :: nq = 8
        real(dp) :: gx(nq), gw(nq), x(2), tau(2), sp, u, t, th
        integer :: p, k

        call gauss_unit(nq, gx, gw)
        err = 0.0_dp
        do p = 1, curve%n_panels
            do k = 1, nq
                call gs_bem_panel_point(curve, p, gx(k), x, tau, sp)
                call filament_data(src, x, tau, u, t)
                th = gs_bem_flux_value(ops%flux_degree, p, gx(k), tcoef)
                err = err + gw(k)*sp*(th - t)**2
            end do
        end do
        err = sqrt(err)
    end subroutine neumann_l2_error

    subroutine exact_energy(curve, src, energy)
        !! int (1/R)|grad u|^2 over the exterior = -int t u ds.
        type(gs_bem_curve_t), intent(in) :: curve
        real(dp), intent(in) :: src(2)
        real(dp), intent(out) :: energy
        integer, parameter :: nq = 10
        real(dp) :: gx(nq), gw(nq), x(2), tau(2), sp, u, t
        integer :: p, k

        call gauss_unit(nq, gx, gw)
        energy = 0.0_dp
        do p = 1, curve%n_panels
            do k = 1, nq
                call gs_bem_panel_point(curve, p, gx(k), x, tau, sp)
                call filament_data(src, x, tau, u, t)
                energy = energy - gw(k)*sp*t*u
            end do
        end do
    end subroutine exact_energy

    subroutine test_exterior_dirichlet(shape, deg)
        integer, intent(in) :: shape, deg
        integer, parameter :: levels(3) = [12, 24, 48]
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ops
        real(dp), allocatable :: u(:), t(:), s(:, :)
        real(dp) :: err(3), eerr(3), energy, src(2), rate_t, rate_e
        integer :: i, status
        character(len=40) :: tag

        write (tag, '(a,i0,a,i0)') "shape ", shape, ", flux degree ", deg
        src = source_point(shape)
        do i = 1, 3
            call make_curve(shape, levels(i), curve)
            call gs_bem_assemble(curve, deg, ops, status)
            call nodal_trace(curve, src, u)
            call gs_bem_exterior_neumann(ops, u, t, status)
            call neumann_l2_error(curve, ops, src, t, err(i))
            call gs_bem_steklov_symmetric(ops, s, status)
            call exact_energy(curve, src, energy)
            eerr(i) = abs(dot_product(u, matmul(s, u)) - energy)/energy
        end do
        rate_t = log(err(2)/err(3))/log(2.0_dp)
        rate_e = log(eerr(2)/eerr(3))/log(2.0_dp)
        print '(a,a,3es11.3,a,f6.2,a,3es11.3,a,f6.2)', "  Neumann ", trim(tag), &
            err, " rate", rate_t, "  energy", eerr, " rate", rate_e
        call record_condition(err(3) < err(2) .and. err(2) < err(1), &
            "Neumann data converge monotonically, "//trim(tag))
        ! The trace is only P2, so the flux order saturates at 2.
        call record_condition(rate_t > real(min(deg + 1, 2), dp) - 0.3_dp, &
            "Neumann L2 order min(flux degree + 1, 2), "//trim(tag))
        call record_condition(eerr(3) < eerr(2) .and. rate_e > 2.0_dp, &
            "exterior energy <S u, u> converges, "//trim(tag))
    end subroutine test_exterior_dirichlet

    subroutine test_calderon(shape)
        !! <W u, v> against <M^T V^{-1}(1/2 M - K) P (1/2 + K) u, v>,
        !! P the L2 projection onto the trace space (V W = 1/4 - K^2).
        integer, intent(in) :: shape
        integer, parameter :: levels(3) = [12, 24, 48]
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ops
        real(dp), allocatable :: u(:), v(:), g(:), c(:, :), mx(:, :)
        real(dp) :: err(3), lhs, rhs, src(2), srcv(2)
        integer :: i, status, info, ipiv(512)
        character(len=16) :: tag

        write (tag, '(a,i0)') "shape ", shape
        src = source_point(shape)
        srcv = src + [-0.15_dp, 0.1_dp]
        do i = 1, 3
            call make_curve(shape, levels(i), curve)
            call gs_bem_assemble(curve, 2, ops, status)
            call nodal_trace(curve, src, u)
            call nodal_trace(curve, srcv, v)
            ! g = Mxx^{-1} (1/2 Mxx + Kxx) u
            mx = ops%mxx
            g = matmul(0.5_dp*ops%mxx + ops%kxx, u)
            call solve_dense(mx, g, info, ipiv)
            ! c = V^{-1}(1/2 M - K) g
            c = reshape(matmul(0.5_dp*ops%m - ops%k, g), [ops%ny, 1])
            call solve_spd(ops%v, c)
            rhs = dot_product(matmul(transpose(ops%m), c(:, 1)), v)
            lhs = dot_product(matmul(ops%w, u), v)
            err(i) = abs(lhs - rhs)/abs(lhs)
        end do
        print '(a,a,3es11.3,a,f6.2)', "  Calderon ", tag, err, " rate", &
            log(err(2)/err(3))/log(2.0_dp)
        call record_condition(err(2) < 0.05_dp*err(1), &
            "Galerkin Calderon V W = 1/4 - K^2 converges, "//trim(tag))
        call record_condition(err(3) < 1.0e-6_dp, &
            "Calderon identity small on the fine mesh, "//trim(tag))
    end subroutine test_calderon

    subroutine solve_dense(a, b, info, ipiv)
        real(dp), intent(inout) :: a(:, :), b(:)
        integer, intent(out) :: info, ipiv(:)
        interface
            subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                import :: dp
                integer, intent(in) :: n, nrhs, lda, ldb
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
                integer, intent(out) :: ipiv(*), info
            end subroutine dgesv
        end interface

        call dgesv(size(a, 1), 1, a, size(a, 1), ipiv, b, size(a, 1), info)
    end subroutine solve_dense

    subroutine solve_spd(a, b)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(inout) :: b(:, :)
        real(dp), allocatable :: f(:, :)
        integer :: info
        interface
            subroutine dposv(uplo, n, nrhs, a, lda, b, ldb, info)
                import :: dp
                character(len=1), intent(in) :: uplo
                integer, intent(in) :: n, nrhs, lda, ldb
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
                integer, intent(out) :: info
            end subroutine dposv
        end interface

        f = a
        call dposv('L', size(a, 1), size(b, 2), f, size(a, 1), b, size(a, 1), &
            info)
    end subroutine solve_spd

    subroutine test_jn_and_ntd()
        integer, parameter :: levels(3) = [12, 24, 48]
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ops
        real(dp), allocatable :: u(:), s(:, :), sj(:, :), n(:, :), x(:)
        real(dp) :: diff(3), esym, ejn, res, dual
        integer :: i, status

        do i = 1, 3
            call make_curve(2, levels(i), curve)
            call gs_bem_assemble(curve, 1, ops, status)
            call nodal_trace(curve, d_source, u)
            call gs_bem_steklov_symmetric(ops, s, status)
            call gs_bem_steklov_jn(ops, sj, status)
            esym = dot_product(u, matmul(s, u))
            ejn = dot_product(u, matmul(sj, u))
            diff(i) = abs(esym - ejn)/esym
        end do
        print '(a,3es11.3)', "  JN vs symmetric energy difference ", diff
        call record_condition(diff(3) < diff(2) .and. diff(2) < diff(1), &
            "Johnson-Nedelec energy converges to the symmetric form")

        call gs_bem_neumann_to_dirichlet(s, n, status)
        call record_condition(status == 0, "NtD map from Cholesky inverse")
        x = matmul(n, matmul(s, u))
        res = maxval(abs(x - u))/maxval(abs(u))
        call record_condition(res < 1.0e-8_dp, "N_h S_h u = u")
        dual = gs_bem_dual_norm_squared(n, matmul(s, u))
        call record_condition(abs(dual - esym)/esym < 1.0e-8_dp, &
            "dual norm of S u equals the energy <S u, u>")
    end subroutine test_jn_and_ntd

end program test_axisymmetric_gs_boundary_operators_2d
