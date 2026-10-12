program test_axisymmetric_gs_kernel
    !! Independent oracles for the axisymmetric Delta* Green's function:
    !! AGM elliptic integrals, the Neumann integral of the mutual
    !! inductance, on-axis Biot-Savart, the dipole far field, Green's
    !! identity on closed curves, the single-layer jump relation, high-order
    !! finite differences for Delta* G, all derivatives and the coil average.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_axisymmetric_gs_kernel, only: gs_green, gs_green_gradient, &
        gs_green_hessian, gs_delta_star, gs_green_conormal, gs_green_log_split, &
        gs_green_log_split_gradient, gs_green_jvp, gs_green_vjp, &
        gs_green_gradient_jvp, gs_green_gradient_vjp, gs_rect_coil_t, &
        gs_rect_coil_flux, gs_rect_coil_flux_gradient, &
        gs_rect_coil_flux_converged, gs_rect_coil_flux_jvp, &
        gs_rect_coil_flux_vjp, gs_mu0
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    logical :: all_passed

    all_passed = .true.

    call test_symmetry()
    call test_mutual_inductance()
    call test_delta_star_fd()
    call test_helmholtz_pair()
    call test_dipole_far_field()
    call test_green_identity()
    call test_single_layer_jump()
    call test_derivatives_fd()
    call test_log_split()
    call test_products()
    call test_rect_coil()

    call check_summary("Axisymmetric Grad-Shafranov kernel")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

    function configs(i) result(x)
        !! Test configurations (R, Z, R', Z'): moderate, near, far, tiny m.
        integer, intent(in) :: i
        real(dp) :: x(4)

        select case (i)
        case (1)
            x = [1.0_dp, 0.3_dp, 1.4_dp, -0.2_dp]
        case (2)
            x = [2.0_dp, 0.5_dp, 2.01_dp, 0.5_dp]
        case (3)
            x = [0.7_dp, 0.1_dp, 3.0_dp, 2.0_dp]
        case (4)
            x = [1.0_dp, 0.0_dp, 1.0_dp, 200.0_dp]
        case (5)
            x = [1.0_dp, 0.0_dp, 1.0_dp, 1.808_dp]
        case (6)
            x = [1.0_dp, 0.0_dp, 1.0_dp, 2.211_dp]
        case (7)
            x = [1.0_dp, 0.0_dp, 1.001_dp, 0.002_dp]
        case default
            x = [2.0_dp, 1.0_dp, 0.5_dp, -1.0_dp]
        end select
    end function configs

    function gval(x) result(g)
        real(dp), intent(in) :: x(4)
        real(dp) :: g

        g = gs_green(x(1), x(2), x(3), x(4))
    end function gval

    function ggrad(x) result(grad)
        real(dp), intent(in) :: x(4)
        real(dp) :: grad(4), g

        call gs_green_gradient(x(1), x(2), x(3), x(4), g, grad)
    end function ggrad

    function dist(x) result(d)
        real(dp), intent(in) :: x(4)
        real(dp) :: d

        d = hypot(x(1) - x(3), x(2) - x(4))
    end function dist

    function lcg_unit(state) result(u)
        integer(8), intent(inout) :: state
        real(dp) :: u

        state = mod(1103515245_8*state + 12345_8, 2147483648_8)
        u = real(state, dp)/2147483648.0_dp
    end function lcg_unit

    subroutine test_symmetry()
        real(dp) :: x(4), g1, g2, a(4), b(4), worst, worst_grad
        integer :: i
        integer(8) :: state

        state = 7_8
        worst = 0.0_dp
        worst_grad = 0.0_dp
        do i = 1, 200
            x = [0.2_dp + 3.0_dp*lcg_unit(state), 2.0_dp*lcg_unit(state), &
                0.2_dp + 3.0_dp*lcg_unit(state), 2.0_dp*lcg_unit(state)]
            g1 = gs_green(x(1), x(2), x(3), x(4))
            g2 = gs_green(x(3), x(4), x(1), x(2))
            worst = worse(worst, abs(g1 - g2)/abs(g1))
            call gs_green_gradient(x(1), x(2), x(3), x(4), g1, a)
            call gs_green_gradient(x(3), x(4), x(1), x(2), g2, b)
            worst_grad = worse(worst_grad, maxval(abs(a - [b(3), b(4), b(1), b(2)])) &
                /maxval(abs(a)))
        end do
        call record_condition(worst < 1.0e-14_dp, "G(x,x') = G(x',x)")
        call record_condition(worst_grad < 1.0e-13_dp, &
            "gradient respects the point exchange symmetry")
        call record_condition(gs_green(0.0_dp, 0.3_dp, 1.0_dp, 0.0_dp) == 0.0_dp, &
            "G vanishes on the symmetry axis")
    end subroutine test_symmetry

    subroutine agm_ke(q, kk, ee)
        !! K and E at m = 1 - q by the arithmetic-geometric mean.
        real(dp), intent(in) :: q
        real(dp), intent(out) :: kk, ee
        real(dp) :: a, b, c, summ, pw
        integer :: n

        a = 1.0_dp
        b = sqrt(q)
        c = sqrt(1.0_dp - q)
        pw = 0.5_dp
        summ = pw*c*c
        do n = 1, 40
            c = 0.5_dp*(a - b)
            pw = 2.0_dp*pw
            summ = summ + pw*c*c
            call agm_step(a, b)
            if (abs(c) < 1.0e-17_dp) exit
        end do
        kk = pi/(2.0_dp*a)
        ee = kk*(1.0_dp - summ)
    end subroutine agm_ke

    subroutine agm_step(a, b)
        real(dp), intent(inout) :: a, b
        real(dp) :: an

        an = 0.5_dp*(a + b)
        b = sqrt(a*b)
        a = an
    end subroutine agm_step

    function neumann_inductance(x, npts) result(mu)
        !! Mutual inductance / mu0 of two coaxial loops from the periodic
        !! trapezoid rule applied to the Neumann integral.
        real(dp), intent(in) :: x(4)
        integer, intent(in) :: npts
        real(dp) :: mu, phi, dd, rr, summ
        integer :: i

        dd = (x(1) - x(3))**2 + (x(2) - x(4))**2
        rr = x(1)*x(3)
        summ = 0.0_dp
        do i = 0, npts - 1
            phi = 2.0_dp*pi*(real(i, dp) + 0.5_dp)/real(npts, dp)
            summ = summ + cos(phi)/sqrt(dd + 4.0_dp*rr*sin(0.5_dp*phi)**2)
        end do
        mu = 0.5_dp*rr*summ*2.0_dp*pi/real(npts, dp)
    end function neumann_inductance

    subroutine test_mutual_inductance()
        real(dp) :: x(4), m_green, m_agm, m_neu, q, kk, ee, k, s2
        real(dp) :: worst_agm, worst_neu
        integer :: i

        worst_agm = 0.0_dp
        worst_neu = 0.0_dp
        do i = 1, 8
            x = configs(i)
            if (i == 4) cycle
            s2 = (x(1) + x(3))**2 + (x(2) - x(4))**2
            q = ((x(1) - x(3))**2 + (x(2) - x(4))**2)/s2
            call agm_ke(q, kk, ee)
            k = sqrt(4.0_dp*x(1)*x(3)/s2)
            m_agm = sqrt(x(1)*x(3))*((2.0_dp/k - k)*kk - 2.0_dp/k*ee)
            m_green = 2.0_dp*pi*gs_green(x(1), x(2), x(3), x(4))
            m_neu = neumann_inductance(x, 400000)
            worst_agm = worse(worst_agm, abs(m_green - m_agm)/abs(m_agm))
            worst_neu = worse(worst_neu, abs(m_green - m_neu)/abs(m_neu))
        end do
        call record_condition(worst_agm < 2.0e-13_dp, &
            "2 pi G equals mu0-normalised M = sqrt(RR')((2/k-k)K-(2/k)E) (AGM)")
        call record_condition(worst_neu < 1.0e-12_dp, &
            "2 pi G equals the Neumann double line integral (to q = 2e-7)")
    end subroutine test_mutual_inductance

    subroutine test_delta_star_fd()
        !! Delta* G from 8th-order differences of G values only.
        real(dp), parameter :: c1(4) = [4.0_dp/5.0_dp, -1.0_dp/5.0_dp, &
            4.0_dp/105.0_dp, -1.0_dp/280.0_dp]
        real(dp), parameter :: c2(4) = [8.0_dp/5.0_dp, -1.0_dp/5.0_dp, &
            8.0_dp/315.0_dp, -1.0_dp/560.0_dp]
        real(dp), parameter :: c20 = -205.0_dp/72.0_dp
        real(dp) :: x(4), h, grr, gr, gzz, gc, fd, scale, an, an_s, an_f
        real(dp) :: worst_fd, worst_an
        integer :: i, j, which

        worst_fd = 0.0_dp
        worst_an = 0.0_dp
        do which = 1, 2
            do i = 1, 3
                x = configs(i)
                if (i == 3) x = [1.0_dp, 0.3_dp, 1.4_dp, -0.2_dp]
                if (which == 2) x = [x(3), x(4), x(1), x(2)]
                h = 0.04_dp*dist(x)
                gc = gval(x)
                grr = c20*gc/h**2
                gr = 0.0_dp
                gzz = c20*gc/h**2
                do j = 1, 4
                    gr = gr + c1(j)*(gval(shift(x, 1, j*h)) - gval(shift(x, 1, -j*h)))/h
                    grr = grr + c2(j)*(gval(shift(x, 1, j*h)) &
                        + gval(shift(x, 1, -j*h)))/h**2
                    gzz = gzz + c2(j)*(gval(shift(x, 2, j*h)) &
                        + gval(shift(x, 2, -j*h)))/h**2
                end do
                fd = grr - gr/x(1) + gzz
                scale = abs(grr) + abs(gzz)
                worst_fd = worse(worst_fd, abs(fd)/scale)
                call gs_delta_star(x(1), x(2), x(3), x(4), an_f, an_s)
                an = an_f
                scale = 1.0_dp/dist(x)**2
                worst_an = worse(worst_an, abs(an_f)/scale)
            end do
        end do
        call record_condition(worst_fd < 1.0e-9_dp, &
            "Delta* G = 0 away from the source (8th-order finite differences)")
        call record_condition(worst_an < 1.0e-9_dp, &
            "analytic Hessian satisfies Delta* G = 0")
        x = configs(1)
        call gs_delta_star(x(1), x(2), x(3), x(4), an_f, an_s)
        call record_condition(abs(an_s) < 1.0e-9_dp, &
            "Delta* G = 0 also in the source variables")
    end subroutine test_delta_star_fd

    function worse(a, b) result(w)
        !! max(a, b) that propagates NaN as a failure.
        real(dp), intent(in) :: a, b
        real(dp) :: w

        if (a /= a .or. b /= b) then
            w = huge(1.0_dp)
        else
            w = max(a, b)
        end if
    end function worse

    function step(x) result(h)
        !! Finite-difference step that keeps R, R' positive.
        real(dp), intent(in) :: x(4)
        real(dp) :: h

        h = 0.02_dp*min(dist(x), 0.5_dp*min(x(1), x(3)))
    end function step

    function shift(x, idx, delta) result(y)
        real(dp), intent(in) :: x(4), delta
        integer, intent(in) :: idx
        real(dp) :: y(4)

        y = x
        y(idx) = y(idx) + delta
    end function shift

    subroutine test_helmholtz_pair()
        real(dp), parameter :: radius = 1.0_dp, rax = 1.0e-4_dp
        real(dp) :: zs(4), bz, br, bz_ref, br_ref, g1(4), g2(4), g
        real(dp) :: worst_bz, worst_br, zc
        integer :: i, k

        zs = [0.0_dp, 0.3_dp, -0.7_dp, 2.0_dp]
        worst_bz = 0.0_dp
        worst_br = 0.0_dp
        do i = 1, 4
            call gs_green_gradient(rax, zs(i), radius, 0.5_dp*radius, g, g1)
            call gs_green_gradient(rax, zs(i), radius, -0.5_dp*radius, g, g2)
            bz = gs_mu0*(g1(1) + g2(1))/rax
            br = -gs_mu0*(g1(2) + g2(2))/rax
            bz_ref = 0.0_dp
            br_ref = 0.0_dp
            do k = -1, 1, 2
                zc = 0.5_dp*real(k, dp)*radius
                bz_ref = bz_ref + gs_mu0*radius**2 &
                    /(2.0_dp*(radius**2 + (zs(i) - zc)**2)**1.5_dp)
                br_ref = br_ref + 0.5_dp*rax*gs_mu0*3.0_dp*radius**2*(zs(i) - zc)/ &
                    (2.0_dp*(radius**2 + (zs(i) - zc)**2)**2.5_dp)
            end do
            worst_bz = worse(worst_bz, abs(bz - bz_ref)/abs(bz_ref))
            worst_br = worse(worst_br, abs(br - br_ref)/max(abs(br_ref), 1.0e-30_dp))
        end do
        call record_condition(worst_bz < 1.0e-6_dp, &
            "Helmholtz pair: axial B_z from psi matches Biot-Savart")
        call record_condition(worst_br < 1.0e-6_dp, &
            "Helmholtz pair: near-axis B_R = -(R/2) dB_z/dz")
        call gs_green_gradient(1.0e-6_dp, 0.0_dp, radius, 0.5_dp*radius, g, g1)
        call gs_green_gradient(1.0e-6_dp, 0.0_dp, radius, -0.5_dp*radius, g, g2)
        bz = gs_mu0*(g1(1) + g2(1))/1.0e-6_dp
        call record_condition(abs(bz/(gs_mu0*(0.8_dp**1.5_dp/radius)) - 1.0_dp) &
            < 1.0e-9_dp, &
            "Helmholtz centre field (4/5)^(3/2) mu0 I/a")
    end subroutine test_helmholtz_pair

    subroutine test_dipole_far_field()
        real(dp) :: radius, rr, zz, g, dip, theta
        real(dp) :: worst_1e3, worst_1e5
        integer :: i

        worst_1e3 = 0.0_dp
        worst_1e5 = 0.0_dp
        radius = 1.0_dp
        do i = 1, 3
            theta = (0.3_dp + 0.4_dp*real(i, dp))
            rr = 1.0e3_dp*sin(theta)
            zz = 1.0e3_dp*cos(theta)
            g = gs_green(rr, zz, radius, 0.0_dp)
            dip = radius**2*rr**2/(4.0_dp*(rr**2 + zz**2)**1.5_dp)
            worst_1e3 = worse(worst_1e3, abs(g/dip - 1.0_dp))
            rr = 1.0e5_dp*sin(theta)
            zz = 1.0e5_dp*cos(theta)
            g = gs_green(rr, zz, radius, 0.0_dp)
            dip = radius**2*rr**2/(4.0_dp*(rr**2 + zz**2)**1.5_dp)
            worst_1e5 = worse(worst_1e5, abs(g/dip - 1.0_dp))
        end do
        call record_condition(worst_1e3 < 5.0e-6_dp, &
            "far field: dipole psi = mu0 m R^2/(4 pi r^3) at r/a = 1e3")
        call record_condition(worst_1e5 < 5.0e-10_dp, &
            "far field keeps full relative accuracy at r/a = 1e5 (m ~ 1e-10)")
    end subroutine test_dipole_far_field

    subroutine curve_point(theta, r0, z0, r, z, n, wdl)
        !! Smooth closed CCW curve, outward unit normal n and |t| dtheta.
        real(dp), intent(in) :: theta, r0, z0
        real(dp), intent(out) :: r, z, n(2), wdl
        real(dp) :: tr, tz

        r = r0 + 0.5_dp*cos(theta) + 0.1_dp*cos(2.0_dp*theta)
        z = z0 + 0.7_dp*sin(theta)
        tr = -0.5_dp*sin(theta) - 0.2_dp*sin(2.0_dp*theta)
        tz = 0.7_dp*cos(theta)
        wdl = hypot(tr, tz)
        n = [tz, -tr]/wdl
    end subroutine curve_point

    subroutine test_green_identity()
        integer, parameter :: npts = 800
        real(dp), parameter :: r0 = 2.0_dp, z0 = 0.1_dp
        real(dp) :: theta, r, z, n(2), wdl, dn_f, dn_s, dn_b
        real(dp) :: src(2), total_f, total_s, nrm(2)
        integer :: i, k
        real(dp) :: expected(3), got_f(3), got_s(3)

        do k = 1, 3
            select case (k)
            case (1)
                src = [r0 + 0.25_dp, z0 + 0.15_dp]
                expected(k) = -1.0_dp
            case (2)
                src = [r0 + 1.3_dp, z0 + 0.4_dp]
                expected(k) = 0.0_dp
            case default
                src = [r0 + 0.58_dp, z0]
                expected(k) = -1.0_dp
            end select
            total_f = 0.0_dp
            total_s = 0.0_dp
            do i = 0, npts - 1
                theta = 2.0_dp*pi*(real(i, dp) + 0.5_dp)/real(npts, dp)
                call curve_point(theta, r0, z0, r, z, n, wdl)
                nrm = n
                ! curve point as field point, fixed source
                call gs_green_conormal(r, z, src(1), src(2), nrm, nrm, dn_f, dn_s, dn_b)
                total_f = total_f + dn_f*wdl
                ! curve point as source point, fixed field point
                call gs_green_conormal(src(1), src(2), r, z, nrm, nrm, dn_f, dn_s, dn_b)
                total_s = total_s + dn_s*wdl
            end do
            got_f(k) = total_f*2.0_dp*pi/real(npts, dp)
            got_s(k) = total_s*2.0_dp*pi/real(npts, dp)
        end do
        call record_condition(maxval(abs(got_f - expected)) < 1.0e-10_dp, &
            "Green identity: closed-curve integral of (1/R) dG/dn = -1 or 0")
        call record_condition(maxval(abs(got_s - expected)) < 1.0e-10_dp, &
            "Green identity with the source-point conormal derivative")
    end subroutine test_green_identity

    subroutine test_single_layer_jump()
        !! Unit density on a circle: [ (1/R) d_n u ] = -sigma across the sheet.
        integer, parameter :: npts = 40000
        real(dp), parameter :: r0 = 2.0_dp, a = 0.5_dp
        real(dp) :: jump(2), deltas(2), theta, y(2), g, grad(4), side(2)
        real(dp) :: dl, estimate, xr
        integer :: i, kd, ks

        deltas = [0.004_dp, 0.002_dp]
        dl = 2.0_dp*pi*a/real(npts, dp)
        do kd = 1, 2
            do ks = 1, 2
                xr = r0 + a + merge(1.0_dp, -1.0_dp, ks == 1)*deltas(kd)
                side(ks) = 0.0_dp
                do i = 0, npts - 1
                    theta = 2.0_dp*pi*(real(i, dp) + 0.5_dp)/real(npts, dp)
                    y = [r0 + a*cos(theta), a*sin(theta)]
                    call gs_green_gradient(xr, 0.0_dp, y(1), y(2), g, grad)
                    side(ks) = side(ks) + grad(1)*dl
                end do
                side(ks) = side(ks)/xr
            end do
            jump(kd) = side(1) - side(2)
        end do
        estimate = 2.0_dp*jump(2) - jump(1)
        call record_condition(abs(estimate + 1.0_dp) < 1.0e-4_dp, &
            "single-layer sheet: jump of (1/R) dn(u) is -sigma")
    end subroutine test_single_layer_jump

    subroutine test_derivatives_fd()
        real(dp), parameter :: c1(4) = [4.0_dp/5.0_dp, -1.0_dp/5.0_dp, &
            4.0_dp/105.0_dp, -1.0_dp/280.0_dp]
        real(dp) :: x(4), h, fdg(4), fdh(4, 4), g, grad(4), hess(4, 4)
        real(dp) :: worst_g, worst_h, worst_sym
        integer :: i, a, j

        worst_g = 0.0_dp
        worst_h = 0.0_dp
        worst_sym = 0.0_dp
        do i = 1, 8
            x = configs(i)
            h = step(x)
            call gs_green_hessian(x(1), x(2), x(3), x(4), g, grad, hess)
            do a = 1, 4
                fdg(a) = 0.0_dp
                fdh(:, a) = 0.0_dp
                do j = 1, 4
                    fdg(a) = fdg(a) + c1(j)*(gval(shift(x, a, j*h)) &
                        - gval(shift(x, a, -j*h)))/h
                    fdh(:, a) = fdh(:, a) + c1(j)*(ggrad(shift(x, a, j*h)) &
                        - ggrad(shift(x, a, -j*h)))/h
                end do
            end do
            worst_g = worse(worst_g, maxval(abs(fdg - grad))/maxval(abs(grad)))
            worst_h = worse(worst_h, maxval(abs(fdh - hess))/maxval(abs(hess)))
            worst_sym = worse(worst_sym, maxval(abs(hess - transpose(hess))) &
                /maxval(abs(hess)))
        end do
        call record_condition(worst_g < 1.0e-9_dp, &
            "gradient matches finite differences (q from 1e-6 to m = 1e-10)")
        call record_condition(worst_h < 1.0e-8_dp, &
            "Hessian matches finite differences of the gradient")
        call record_condition(worst_sym < 1.0e-14_dp, "Hessian is symmetric")
    end subroutine test_derivatives_fd

    subroutine test_log_split()
        real(dp), parameter :: c1(4) = [4.0_dp/5.0_dp, -1.0_dp/5.0_dp, &
            4.0_dp/105.0_dp, -1.0_dp/280.0_dp]
        real(dp) :: x(4), coeff, smooth, dcoeff(4), dsmooth(4), g, grad(4)
        real(dp) :: cp, cm, sp, sm, h, fdc(4), fds(4), lnd, dgrad(4), dd
        real(dp) :: xp(4), xm(4)
        real(dp) :: worst_val, worst_grad, worst_fd
        integer :: i, a, j

        worst_val = 0.0_dp
        worst_grad = 0.0_dp
        worst_fd = 0.0_dp
        do i = 1, 7
            if (i == 3 .or. i == 4) cycle
            x = configs(i)
            call gs_green_log_split_gradient(x(1), x(2), x(3), x(4), coeff, &
                dcoeff, smooth, dsmooth)
            call gs_green_gradient(x(1), x(2), x(3), x(4), g, grad)
            dd = dist(x)
            lnd = log(dd)
            worst_val = worse(worst_val, abs(g - (-lnd*coeff + smooth))/abs(g))
            dgrad = [x(1) - x(3), x(2) - x(4), x(3) - x(1), x(4) - x(2)]/dd**2
            worst_grad = worse(worst_grad, &
                maxval(abs(grad - (-dgrad*coeff - lnd*dcoeff &
                + dsmooth)))/maxval(abs(grad)))
            h = 0.01_dp*(x(1) + x(3))
            do a = 1, 4
                fdc(a) = 0.0_dp
                fds(a) = 0.0_dp
                do j = 1, 4
                    xp = shift(x, a, j*h)
                    xm = shift(x, a, -j*h)
                    call gs_green_log_split(xp(1), xp(2), xp(3), xp(4), cp, sp)
                    call gs_green_log_split(xm(1), xm(2), xm(3), xm(4), cm, sm)
                    fdc(a) = fdc(a) + c1(j)*(cp - cm)/h
                    fds(a) = fds(a) + c1(j)*(sp - sm)/h
                end do
            end do
            worst_fd = worse(worst_fd, maxval(abs(fdc - dcoeff))/maxval(abs(dcoeff)))
            worst_fd = worse(worst_fd, maxval(abs(fds - dsmooth))/maxval(abs(dsmooth)))
        end do
        call record_condition(worst_val < 1.0e-13_dp, &
            "log split: G = -ln(dist) coeff + smooth")
        call record_condition(worst_grad < 1.0e-9_dp, &
            "log split gradients reproduce grad G")
        call record_condition(worst_fd < 1.0e-9_dp, &
            "log split gradients match finite differences (both q branches)")
        call gs_green_log_split_gradient(1.3_dp, 0.2_dp, 1.3_dp, 0.2_dp, coeff, &
            dcoeff, smooth, dsmooth)
        call record_condition(abs(coeff - 1.3_dp/(2.0_dp*pi)) < 1.0e-15_dp .and. &
            smooth == smooth, &
            "log split at coincident points: coeff = R/(2 pi), smooth finite")
    end subroutine test_log_split

    subroutine test_products()
        real(dp) :: x(4), t(4), w(4), g, dg, bar(4), grad(4), dgrad(4), bar2(4)
        real(dp) :: lhs, rhs, wbar(4)

        x = configs(1)
        t = [0.3_dp, -0.7_dp, 0.2_dp, 0.5_dp]
        w = [-0.4_dp, 0.9_dp, 0.1_dp, 0.6_dp]
        call gs_green_jvp(x(1), x(2), x(3), x(4), t, g, dg)
        call gs_green_vjp(x(1), x(2), x(3), x(4), 1.7_dp, bar)
        call record_condition(abs(dot_product(bar, t) - 1.7_dp*dg) < 1.0e-14_dp, &
            "G JVP/VJP dot-product identity")
        call gs_green_gradient_jvp(x(1), x(2), x(3), x(4), t, grad, dgrad)
        call gs_green_gradient_vjp(x(1), x(2), x(3), x(4), w, bar2)
        lhs = dot_product(w, dgrad)
        rhs = dot_product(bar2, t)
        wbar = grad
        call record_condition(abs(lhs - rhs) < 1.0e-12_dp*maxval(abs(wbar)), &
            "gradient JVP/VJP dot-product identity")
    end subroutine test_products

    subroutine test_rect_coil()
        type(gs_rect_coil_t) :: coil, line
        real(dp), parameter :: c1(4) = [4.0_dp/5.0_dp, -1.0_dp/5.0_dp, &
            4.0_dp/105.0_dp, -1.0_dp/280.0_dp]
        real(dp) :: ref, e(5), psi, grad(6), fd(6), p(6), h(6), tol, bar(6)
        real(dp) :: dpsi, psi2
        integer :: orders(5), k, a, j, used, status
        real(dp) :: pp, pm

        coil = gs_rect_coil_t(r=3.0_dp, z=1.0_dp, dr=0.2_dp, dz=0.4_dp, turns=5.0_dp)
        ref = gs_rect_coil_flux(coil, 2.0_dp, 0.0_dp, 48, 48)
        orders = [2, 3, 4, 6, 8]
        do k = 1, 5
            e(k) = abs(gs_rect_coil_flux(coil, 2.0_dp, 0.0_dp, orders(k), &
                orders(k)) - ref)/abs(ref)
        end do
        call record_condition(e(1) > e(2) .and. e(2) > e(3) .and. e(3) > e(4) &
            .and. e(4) > e(5) .and. e(5) < 1.0e-10_dp, &
            "rectangle average converges monotonically under Gauss refinement")
        call gs_rect_coil_flux_converged(coil, 2.0_dp, 0.0_dp, 1.0e-12_dp, psi, &
            used, status)
        call record_condition(status == 0 .and. abs(psi - ref) < 1.0e-12_dp*abs(ref), &
            "adaptive order reaches the stated tolerance")
        call gs_rect_coil_flux_converged(coil, 3.05_dp, 1.1_dp, 1.0e-12_dp, psi, &
            used, status)
        call record_condition(status == 2, "field point inside the rectangle flagged")

        line = gs_rect_coil_t(r=3.0_dp, z=1.0_dp, dr=0.0_dp, dz=0.0_dp, turns=5.0_dp)
        call record_condition(abs(gs_rect_coil_flux(line, 2.0_dp, 0.0_dp, 3, 3) &
            - 5.0_dp*gs_mu0*gs_green(2.0_dp, 0.0_dp, 3.0_dp, &
            1.0_dp)) < 1.0e-18_dp, &
            "zero-size rectangle reduces to the filament flux mu0 N G")

        call gs_rect_coil_flux_gradient(coil, 2.0_dp, 0.0_dp, 12, 12, psi, grad)
        p = [2.0_dp, 0.0_dp, coil%r, coil%z, coil%dr, coil%dz]
        h = 0.03_dp
        do a = 1, 6
            fd(a) = 0.0_dp
            do j = 1, 4
                pp = flux_at(shift6(p, a, j*h(a)))
                pm = flux_at(shift6(p, a, -j*h(a)))
                fd(a) = fd(a) + c1(j)*(pp - pm)/h(a)
            end do
        end do
        tol = 1.0e-8_dp*maxval(abs(grad))
        call record_condition(maxval(abs(fd - grad)) < tol, &
            "rectangle flux gradient matches finite differences")
        call gs_rect_coil_flux_jvp(coil, 2.0_dp, 0.0_dp, 12, 12, &
            [0.1_dp, 0.2_dp, -0.3_dp, 0.4_dp, 0.5_dp, -0.6_dp], &
            psi2, dpsi)
        call gs_rect_coil_flux_vjp(coil, 2.0_dp, 0.0_dp, 12, 12, 2.0_dp, bar)
        call record_condition(abs(dot_product(bar, [0.1_dp, 0.2_dp, -0.3_dp, 0.4_dp, &
            0.5_dp, -0.6_dp]) - 2.0_dp*dpsi) &
            < 1.0e-12_dp*maxval(abs(grad)), &
            "coil flux JVP/VJP dot-product identity")
    end subroutine test_rect_coil

    function flux_at(pv) result(f)
        real(dp), intent(in) :: pv(6)
        real(dp) :: f
        type(gs_rect_coil_t) :: c

        c = gs_rect_coil_t(r=pv(3), z=pv(4), dr=pv(5), dz=pv(6), turns=5.0_dp)
        f = gs_rect_coil_flux(c, pv(1), pv(2), 12, 12)
    end function flux_at

    function shift6(pv, idx, delta) result(y)
        real(dp), intent(in) :: pv(6), delta
        integer, intent(in) :: idx
        real(dp) :: y(6)

        y = pv
        y(idx) = y(idx) + delta
    end function shift6

end program test_axisymmetric_gs_kernel
