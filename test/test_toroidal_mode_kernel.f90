program test_toroidal_mode_kernel
    !! Independent oracles for the toroidal-mode Green's function G_n =
    !! Q_{|n|-1/2}(chi)/(2 pi sqrt(R R')): the phi'-quadrature that defines
    !! it, the closed form for n = 0, the identity G_GS = R R' G_1 with the
    !! Grad-Shafranov kernel, the mode sum back to 1/(4 pi |x - x'|),
    !! symmetry, the n-mode Laplacian by finite differences, Gauss's law on a
    !! small circle, tabulated mpmath values of the log-split series,
    !! derivatives, the log split and the coil flux.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortnum_special_elliptic, only: ellip_k_from_complement
    use fortfem_axisymmetric_gs_kernel, only: gs_green, gs_green_gradient, &
        gs_rect_coil_t, gs_rect_coil_flux, gs_mu0
    use fortfem_toroidal_mode_kernel, only: tm_green, tm_green_modes, &
        tm_green_gradient, tm_green_hessian, tm_green_operator, &
        tm_green_conormal, tm_green_log_split, tm_green_log_split_gradient, &
        tm_green_modes_log_split, tm_log_split_limit, tm_filament_flux, &
        tm_filament_flux_gradient, tm_rect_coil_t, tm_rect_coil_flux, &
        tm_rect_coil_flux_modes, tm_rect_coil_flux_gradient, &
        tm_rect_coil_flux_converged, tm_rect_coil_flux_table
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: orders(4) = [0, 1, 3, 6]
    logical :: all_passed

    all_passed = .true.

    call test_phi_quadrature()
    call test_n0_closed_form()
    call test_grad_shafranov_relation()
    call test_mode_sum()
    call test_symmetry()
    call test_laplacian_fd()
    call test_gauss_law()
    call test_series_reference()
    call test_derivatives_fd()
    call test_log_split()
    call test_coil_flux()

    call check_summary("Toroidal-mode Green's function")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

    function configs(i) result(x)
        !! Test configurations (R, Z, R', Z'): moderate, near, far, tall.
        integer, intent(in) :: i
        real(dp) :: x(4)

        select case (i)
        case (1)
            x = [1.0_dp, 0.3_dp, 1.4_dp, -0.2_dp]
        case (2)
            x = [2.0_dp, 0.5_dp, 2.01_dp, 0.52_dp]
        case (3)
            x = [0.7_dp, 0.1_dp, 3.0_dp, 2.0_dp]
        case (4)
            x = [1.0_dp, 0.0_dp, 1.0_dp, 20.0_dp]
        case (5)
            x = [3.0_dp, 0.4_dp, 3.3_dp, 0.1_dp]
        case default
            x = [2.0_dp, 1.0_dp, 0.5_dp, -1.0_dp]
        end select
    end function configs

    function rel_err(a, b) result(e)
        real(dp), intent(in) :: a, b
        real(dp) :: e

        e = abs(a - b)/max(abs(b), tiny(1.0_dp))
    end function rel_err

    subroutine test_phi_quadrature()
        !! G_n = int_0^{2 pi} cos(n D) dD/(4 pi r(D)): the periodic trapezoid
        !! rule converges geometrically.
        integer, parameter :: nq = 4096
        real(dp) :: x(4), s, d, r3, worst
        integer :: i, j, k

        worst = 0.0_dp
        do i = 1, 6
            x = configs(i)
            do k = 1, size(orders)
                s = 0.0_dp
                do j = 0, nq - 1
                    d = 2.0_dp*pi*real(j, dp)/nq
                    r3 = sqrt((x(1) - x(3))**2 + 4.0_dp*x(1)*x(3)*sin(0.5_dp*d)**2 &
                        + (x(2) - x(4))**2)
                    s = s + cos(orders(k)*d)/(4.0_dp*pi*r3)
                end do
                s = s*2.0_dp*pi/nq
                ! The oscillatory sum cancels for large n: scale by G_0.
                worst = max(worst, abs(tm_green(orders(k), x(1), x(2), x(3), &
                    x(4)) - s)/tm_green(0, x(1), x(2), x(3), x(4)))
            end do
        end do
        call record_condition(worst < 1.0e-13_dp, &
            "G_n equals the phi' integral of cos(n D)/(4 pi r) (n=0,1,3,6)")
    end subroutine test_phi_quadrature

    subroutine test_n0_closed_form()
        !! G_0 = K(m)/(pi s), m = 4 R R'/s^2 (independent: fortnum K).
        real(dp) :: x(4), s, q, worst
        integer :: i

        worst = 0.0_dp
        do i = 1, 6
            x = configs(i)
            s = sqrt((x(1) + x(3))**2 + (x(2) - x(4))**2)
            q = ((x(1) - x(3))**2 + (x(2) - x(4))**2)/s**2
            worst = max(worst, rel_err(tm_green(0, x(1), x(2), x(3), x(4)), &
                ellip_k_from_complement(q)/(pi*s)))
            worst = max(worst, rel_err(tm_green(-3, x(1), x(2), x(3), x(4)), &
                tm_green(3, x(1), x(2), x(3), x(4))))
        end do
        call record_condition(worst < 1.0e-13_dp, &
            "n = 0 equals the elliptic-K kernel K(m)/(pi s); G_-n = G_n")
    end subroutine test_n0_closed_form

    subroutine test_grad_shafranov_relation()
        !! The existing Grad-Shafranov kernel is R R' G_1 (value and gradient).
        real(dp) :: x(4), g, grad(4), gs, gsg(4), worst, worst_g
        integer :: i

        worst = 0.0_dp
        worst_g = 0.0_dp
        do i = 1, 6
            x = configs(i)
            call tm_green_gradient(1, x(1), x(2), x(3), x(4), g, grad)
            call gs_green_gradient(x(1), x(2), x(3), x(4), gs, gsg)
            worst = max(worst, rel_err(x(1)*x(3)*g, gs))
            worst_g = max(worst_g, maxval(abs(x(1)*x(3)*grad &
                + g*[x(3), 0.0_dp, x(1), 0.0_dp] - gsg))/maxval(abs(gsg)))
        end do
        call record_condition(worst < 1.0e-12_dp, &
            "G_GS = R R' G_1 (value)")
        call record_condition(worst_g < 1.0e-11_dp, &
            "G_GS = R R' G_1 (gradient)")
    end subroutine test_grad_shafranov_relation

    subroutine test_mode_sum()
        !! 1/(4 pi |x - x'|) = (1/(2 pi)) sum_n G_|n| e^{i n D} at off-axis
        !! points.
        integer, parameter :: nmodes = 120
        real(dp) :: x(4), g(0:nmodes), s, d, r3, worst
        integer :: i, j, k
        real(dp), parameter :: deltas(3) = [0.0_dp, 1.3_dp, pi]

        worst = 0.0_dp
        do i = 1, 6
            if (i == 2 .or. i == 5) cycle   ! near points converge too slowly
            x = configs(i)
            call tm_green_modes(nmodes, x(1), x(2), x(3), x(4), g)
            do j = 1, size(deltas)
                d = deltas(j)
                s = g(0)
                do k = 1, nmodes
                    s = s + 2.0_dp*g(k)*cos(k*d)
                end do
                s = s/(2.0_dp*pi)
                r3 = sqrt(x(1)**2 + x(3)**2 - 2.0_dp*x(1)*x(3)*cos(d) &
                    + (x(2) - x(4))**2)
                worst = max(worst, rel_err(s, 1.0_dp/(4.0_dp*pi*r3)))
            end do
        end do
        call record_condition(worst < 1.0e-9_dp, &
            "mode sum reconstructs the 3D point-source Green's function")
    end subroutine test_mode_sum

    subroutine test_symmetry()
        real(dp) :: x(4), g1, g2, a(4), b(4), worst, worst_g
        integer :: i, k

        worst = 0.0_dp
        worst_g = 0.0_dp
        do i = 1, 6
            x = configs(i)
            do k = 1, size(orders)
                call tm_green_gradient(orders(k), x(1), x(2), x(3), x(4), g1, a)
                call tm_green_gradient(orders(k), x(3), x(4), x(1), x(2), g2, b)
                worst = max(worst, rel_err(g1, g2))
                worst_g = max(worst_g, maxval(abs(a - [b(3), b(4), b(1), &
                    b(2)]))/maxval(abs(a)))
            end do
        end do
        call record_condition(worst < 1.0e-13_dp, "G_n(x, x') = G_n(x', x)")
        call record_condition(worst_g < 1.0e-12_dp, &
            "gradient exchanges field and source blocks")
    end subroutine test_symmetry

    subroutine test_laplacian_fd()
        !! The n-mode Laplacian vanishes away from the source, by a
        !! fourth-order finite difference in the field point, and by the
        !! analytic Hessian (field and source).  The same stencil without the
        !! n^2 term does not vanish.
        real(dp) :: x(4), h, lap, lap_no_n2, scale, worst, worst_h, worst_s
        real(dp) :: field, source, g0, f(-2:2, 2), nn
        integer :: i, k, a

        worst = 0.0_dp
        worst_h = 0.0_dp
        worst_s = 0.0_dp
        lap_no_n2 = 0.0_dp
        do i = 1, 6
            x = configs(i)
            h = 0.02_dp*min(x(1), 0.5_dp*sqrt((x(1) - x(3))**2 + (x(2) - x(4))**2))
            do k = 1, size(orders)
                nn = real(orders(k), dp)**2
                g0 = tm_green(orders(k), x(1), x(2), x(3), x(4))
                do a = -2, 2
                    f(a, 1) = tm_green(orders(k), x(1) + a*h, x(2), x(3), x(4))
                    f(a, 2) = tm_green(orders(k), x(1), x(2) + a*h, x(3), x(4))
                end do
                ! (1/R) d_R (R d_R u) = u_RR + u_R/R
                lap = d2(f(:, 1), h) + d1(f(:, 1), h)/x(1) + d2(f(:, 2), h)
                scale = abs(g0)/h**2
                worst = max(worst, abs(lap - nn*g0/x(1)**2)/scale)
                if (orders(k) == 6 .and. i == 1) &
                    lap_no_n2 = abs(lap)/(nn*abs(g0)/x(1)**2)
                call tm_green_operator(orders(k), x(1), x(2), x(3), x(4), &
                    field, source)
                worst_h = max(worst_h, abs(field)/scale)
                worst_s = max(worst_s, abs(source)/scale)
            end do
        end do
        call record_condition(worst < 1.0e-7_dp, &
            "n-mode Laplacian is zero away from the source (finite differences)")
        call record_condition(worst_h < 1.0e-8_dp .and. worst_s < 1.0e-8_dp, &
            "n-mode Laplacian is zero (analytic Hessian, field and source)")
        call record_condition(abs(lap_no_n2 - 1.0_dp) < 1.0e-4_dp, &
            "the stencil sees the n^2/R^2 term (n = 6 control)")
    end subroutine test_laplacian_fd

    pure function d1(f, h) result(v)
        real(dp), intent(in) :: f(-2:2), h
        real(dp) :: v

        v = (f(-2) - 8.0_dp*f(-1) + 8.0_dp*f(1) - f(2))/(12.0_dp*h)
    end function d1

    pure function d2(f, h) result(v)
        real(dp), intent(in) :: f(-2:2), h
        real(dp) :: v

        v = (-f(-2) + 16.0_dp*f(-1) - 30.0_dp*f(0) + 16.0_dp*f(1) &
            - f(2))/(12.0_dp*h*h)
    end function d2

    subroutine test_gauss_law()
        !! -div(R grad G_n) + n^2 G_n/R = delta: the flux of R grad G_n through
        !! a small circle around the source is -1 (an independent check of
        !! the 1/(2 pi) normalisation).
        integer, parameter :: nq = 256
        real(dp) :: rho, th, xs(2), xp(2), nrm(2), flux, df, ds, db, worst
        integer :: j, k

        worst = 0.0_dp
        xs = [2.0_dp, 0.5_dp]
        rho = 1.0e-4_dp
        do k = 1, size(orders)
            flux = 0.0_dp
            do j = 0, nq - 1
                th = 2.0_dp*pi*(real(j, dp) + 0.5_dp)/nq
                nrm = [cos(th), sin(th)]
                xp = xs + rho*nrm
                call tm_green_conormal(orders(k), xp(1), xp(2), xs(1), xs(2), &
                    nrm, nrm, df, ds, db)
                flux = flux + df*rho*2.0_dp*pi/nq
            end do
            worst = max(worst, abs(flux + 1.0_dp))
        end do
        call record_condition(worst < 1.0e-6_dp, &
            "flux of R grad G_n around the source is -1 (Gauss's law)")
    end subroutine test_gauss_law

    subroutine test_series_reference()
        !! Q_{n-1/2}(1 + 2 t) and dQ/dchi from mpmath (60 digits) at
        !! R = R' = 1, Z - Z' = 2 sqrt(t): G = Q/(2 pi), dG/dZ = Q' zeta/(2 pi).
        real(dp), parameter :: ts(4) = [1.0e-15_dp, 1.0e-8_dp, 1.0e-3_dp, &
            0.05_dp]
        real(dp), parameter :: q_ref(4, 4) = reshape([ &
            18.655682558575229_dp, 10.596634709104487_dp, &
            4.8392124738455304_dp, 2.8611928721988956_dp, &
            16.655682558575247_dp, 8.5966348100708336_dp, &
            2.8435512378841936_dp, 0.97876028288694116_dp, &
            15.589015891908712_dp, 7.5299688178016262_dp, &
            1.7983259875998871_dp, 0.26068388880330731_dp, &
            14.899261202154446_dp, 6.8402161845512524_dp, &
            1.1595447746429294_dp, 0.050525025814959992_dp], [4, 4])
        real(dp), parameter :: dq_ref(4, 4) = reshape([ &
            -250000000000002.14_dp, -25000001.137079329_dp, &
            -250.41704057305543_dp, -5.1632187536472477_dp, &
            -249999999999993.32_dp, -24999996.338761997_dp, &
            -248.49826841727878_dp, -4.2489441929125246_dp, &
            -249999999999929.36_dp, -24999964.618888372_dp, &
            -239.75461390048635_dp, -2.3218076962409974_dp, &
            -249999999999724.49_dp, -24999868.543661312_dp, &
            -220.71454166643393_dp, -0.78579047907753799_dp], [4, 4])
        real(dp) :: zeta, g, grad(4), worst_q, worst_d
        integer :: i, k

        worst_q = 0.0_dp
        worst_d = 0.0_dp
        do k = 1, 4
            do i = 1, 4
                zeta = 2.0_dp*sqrt(ts(i))
                call tm_green_gradient(orders(k), 1.0_dp, zeta, 1.0_dp, &
                    0.0_dp, g, grad)
                worst_q = max(worst_q, rel_err(2.0_dp*pi*g, q_ref(i, k)))
                worst_d = max(worst_d, rel_err(2.0_dp*pi*grad(2), &
                    dq_ref(i, k)*zeta))
            end do
        end do
        call record_condition(worst_q < 1.0e-13_dp, &
            "Q_{n-1/2} matches mpmath from t = 1e-15 to 0.05")
        call record_condition(worst_d < 1.0e-10_dp, &
            "dQ/dchi matches mpmath from t = 1e-15 to 0.05")
    end subroutine test_series_reference

    subroutine test_derivatives_fd()
        !! Gradient and Hessian against central differences of the lower
        !! order, including a Hessian in the series regime.
        real(dp) :: x(4), g, grad(4), hess(4, 4), gp, gm, h, e(4), worst_g, &
            worst_h, gradp(4), gradm(4), g1(4), scale
        integer :: i, j, k, a

        worst_g = 0.0_dp
        worst_h = 0.0_dp
        do i = 1, 6
            x = configs(i)
            do k = 1, size(orders)
                call tm_green_hessian(orders(k), x(1), x(2), x(3), x(4), g, &
                    grad, hess)
                do a = 1, 4
                    h = 5.0e-4_dp*min(x(1), x(3), sqrt((x(1) - x(3))**2 &
                        + (x(2) - x(4))**2))
                    e = 0.0_dp
                    e(a) = h
                    gp = tm_green(orders(k), x(1) + e(1), x(2) + e(2), &
                        x(3) + e(3), x(4) + e(4))
                    gm = tm_green(orders(k), x(1) - e(1), x(2) - e(2), &
                        x(3) - e(3), x(4) - e(4))
                    worst_g = max(worst_g, abs((gp - gm)/(2.0_dp*h) - grad(a)) &
                        /max(maxval(abs(grad)), 1.0e-300_dp))
                    call tm_green_gradient(orders(k), x(1) + e(1), x(2) + e(2), &
                        x(3) + e(3), x(4) + e(4), gp, gradp)
                    call tm_green_gradient(orders(k), x(1) - e(1), x(2) - e(2), &
                        x(3) - e(3), x(4) - e(4), gm, gradm)
                    g1 = (gradp - gradm)/(2.0_dp*h)
                    scale = maxval(abs(hess))
                    do j = 1, 4
                        worst_h = max(worst_h, abs(g1(j) - hess(j, a))/scale)
                    end do
                end do
            end do
        end do
        call record_condition(worst_g < 1.0e-5_dp, &
            "gradient matches central differences of G_n")
        call record_condition(worst_h < 1.0e-5_dp, &
            "Hessian matches central differences of the gradient")
    end subroutine test_derivatives_fd

    subroutine test_log_split()
        !! G_n = -ln(d) coeff + smooth (value and gradient) against the
        !! direct kernel where both are accurate; coeff -> 1/(2 pi R) at
        !! coincidence; mode-set split equals the single-mode split.
        real(dp), parameter :: dists(3) = [3.0e-2_dp, 1.0e-2_dp, 1.0e-3_dp]
        real(dp) :: x(4), g, grad(4), c, dc(4), s, ds(4), d, lnd, rec, recg(4)
        real(dp) :: worst, worst_g, worst_m, cm(0:7), sm(0:7), c0, s0
        integer :: j, k, m

        worst = 0.0_dp
        worst_g = 0.0_dp
        worst_m = 0.0_dp
        do j = 1, size(dists)
            do k = 1, size(orders)
                x = [2.0_dp, 0.5_dp, 2.0_dp + 0.6_dp*dists(j), &
                    0.5_dp - 0.8_dp*dists(j)]
                d = sqrt((x(1) - x(3))**2 + (x(2) - x(4))**2)
                lnd = log(d)
                call tm_green_gradient(orders(k), x(1), x(2), x(3), x(4), g, &
                    grad)
                call tm_green_log_split_gradient(orders(k), x(1), x(2), x(3), &
                    x(4), c, dc, s, ds)
                rec = -c*lnd + s
                recg(1) = -c*(x(1) - x(3))/d**2 - lnd*dc(1) + ds(1)
                recg(2) = -c*(x(2) - x(4))/d**2 - lnd*dc(2) + ds(2)
                recg(3) = c*(x(1) - x(3))/d**2 - lnd*dc(3) + ds(3)
                recg(4) = c*(x(2) - x(4))/d**2 - lnd*dc(4) + ds(4)
                worst = max(worst, rel_err(rec, g))
                worst_g = max(worst_g, maxval(abs(recg - grad))/ &
                    maxval(abs(grad)))
            end do
        end do
        call tm_green_modes_log_split(7, 2.0_dp, 0.5_dp, 2.01_dp, 0.49_dp, &
            cm, sm)
        do m = 0, 7
            call tm_green_log_split(m, 2.0_dp, 0.5_dp, 2.01_dp, 0.49_dp, c0, s0)
            worst_m = max(worst_m, rel_err(cm(m), c0), rel_err(sm(m), s0))
        end do
        call tm_green_log_split(3, 2.0_dp, 0.5_dp, 2.0_dp, 0.5_dp, c, s)
        call record_condition(worst < 1.0e-10_dp, &
            "G_n = -ln(d) coeff + smooth (value)")
        call record_condition(worst_g < 1.0e-9_dp, &
            "log split reproduces the gradient")
        call record_condition(worst_m == 0.0_dp, &
            "mode-set split equals the single-mode split")
        call record_condition(rel_err(c, 1.0_dp/(4.0_dp*pi)) < 1.0e-14_dp, &
            "coeff = 1/(2 pi R) at coincidence")
        call record_condition(tm_log_split_limit(0) >= tm_log_split_limit(7), &
            "the split limit shrinks with the mode number")
    end subroutine test_log_split

    subroutine test_coil_flux()
        type(tm_rect_coil_t) :: coil
        type(gs_rect_coil_t) :: gcoil
        real(dp) :: r, z, psi, grad(6), psi_m(0:6), fil, fg(4), fd, worst
        real(dp) :: p1, p2, e(6), tol, table(3, 2), vals(0:6)
        integer :: k, a, order, status, i

        coil = tm_rect_coil_t(r=2.4_dp, z=0.6_dp, dr=0.2_dp, dz=0.3_dp, &
            turns=3.0_dp)
        gcoil = gs_rect_coil_t(r=2.4_dp, z=0.6_dp, dr=0.2_dp, dz=0.3_dp, &
            turns=3.0_dp)
        r = 1.7_dp
        z = -0.2_dp
        ! n = 1 flux times mu0 R is the flux per radian of the GS kernel.
        psi = tm_rect_coil_flux(1, coil, r, z, 8, 8)
        call record_condition(rel_err(gs_mu0*r*psi, &
            gs_rect_coil_flux(gcoil, r, z, 8, 8)) < 1.0e-12_dp, &
            "n = 1 rectangle flux equals the Grad-Shafranov coil flux")

        ! All modes at once.
        call tm_rect_coil_flux_modes(6, coil, r, z, 6, 6, psi_m)
        worst = 0.0_dp
        do k = 0, 6
            worst = max(worst, rel_err(psi_m(k), tm_rect_coil_flux(k, coil, r, &
                z, 6, 6)))
        end do
        call record_condition(worst < 1.0e-13_dp, &
            "mode-set coil flux equals the per-n flux")

        ! The thin-rectangle limit is the filament times the turns.
        coil%dr = 1.0e-6_dp
        coil%dz = 1.0e-6_dp
        fil = tm_filament_flux(3, coil%r, coil%z, r, z)
        call record_condition(rel_err(tm_rect_coil_flux(3, coil, r, z, 2, 2), &
            coil%turns*fil) < 1.0e-9_dp, "thin rectangle equals the filament")
        coil%dr = 0.2_dp
        coil%dz = 0.3_dp

        ! Direct 3D check of the filament: u(R, Z, phi = 0) of the line
        ! density cos(n phi') is Rc int cos(n phi') dphi'/(4 pi r).
        fd = 0.0_dp
        do i = 0, 2047
            p1 = 2.0_dp*pi*real(i, dp)/2048
            p2 = sqrt(r**2 + coil%r**2 - 2.0_dp*r*coil%r*cos(p1) &
                + (z - coil%z)**2)
            fd = fd + cos(6.0_dp*p1)/(4.0_dp*pi*p2)
        end do
        fd = fd*2.0_dp*pi/2048*coil%r
        call record_condition(rel_err(tm_filament_flux(6, coil%r, coil%z, r, &
            z), fd) < 1.0e-12_dp, &
            "filament flux equals the 3D line-source potential projected on n")

        ! Gradients against central differences.
        call tm_rect_coil_flux_gradient(3, coil, r, z, 8, 8, psi, grad)
        worst = 0.0_dp
        do a = 1, 6
            e = 0.0_dp
            e(a) = 1.0e-5_dp
            worst = max(worst, abs((shifted(3, coil, r, z, e) &
                - shifted(3, coil, r, z, -e))/2.0e-5_dp - grad(a)) &
                /maxval(abs(grad)))
        end do
        call record_condition(worst < 1.0e-7_dp, &
            "rectangle flux gradient matches central differences")
        call tm_filament_flux_gradient(3, coil%r, coil%z, r, z, fil, fg)
        worst = 0.0_dp
        worst = max(worst, abs((tm_filament_flux(3, coil%r, coil%z, r + 1e-5_dp, z) &
            - tm_filament_flux(3, coil%r, coil%z, r - 1e-5_dp, z))/2e-5_dp &
            - fg(1))/abs(fg(1)))
        worst = max(worst, abs((tm_filament_flux(3, coil%r + 1e-5_dp, coil%z, r, z) &
            - tm_filament_flux(3, coil%r - 1e-5_dp, coil%z, r, z))/2e-5_dp &
            - fg(3))/abs(fg(3)))
        call record_condition(worst < 1.0e-7_dp, &
            "filament flux gradient (field and source radius) matches differences")

        tol = 1.0e-10_dp
        call tm_rect_coil_flux_converged(3, coil, r, z, tol, psi, order, status)
        call record_condition(status == 0 .and. rel_err(psi, &
            tm_rect_coil_flux(3, coil, r, z, 64, 64)) < 1.0e-9_dp, &
            "converged coil flux")
        call tm_rect_coil_flux_converged(3, coil, coil%r, coil%z, tol, psi, &
            order, status)
        call record_condition(status == 2, "field point inside the coil flagged")

        call tm_rect_coil_flux_table(3, [coil, coil], [1.5_dp, 1.7_dp, 1.9_dp], &
            [0.0_dp, -0.2_dp, 0.1_dp], 6, 6, table)
        vals(0) = tm_rect_coil_flux(3, coil, 1.7_dp, -0.2_dp, 6, 6)
        call record_condition(table(2, 2) == vals(0) .and. &
            table(2, 1) == table(2, 2), "coil table entries are independent")
    end subroutine test_coil_flux

    function shifted(n, coil, r, z, e) result(psi)
        !! Flux with (r, z, coil r, z, dr, dz) shifted by e.
        integer, intent(in) :: n
        type(tm_rect_coil_t), intent(in) :: coil
        real(dp), intent(in) :: r, z, e(6)
        real(dp) :: psi
        type(tm_rect_coil_t) :: c

        c = coil
        c%r = c%r + e(3)
        c%z = c%z + e(4)
        c%dr = c%dr + e(5)
        c%dz = c%dz + e(6)
        psi = tm_rect_coil_flux(n, c, r + e(1), z + e(2), 8, 8)
    end function shifted

end program test_toroidal_mode_kernel
