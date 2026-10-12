module fortfem_axisymmetric_gs_kernel
    !! Free-space Green's function of the Grad-Shafranov operator
    !!
    !!     Delta* psi = R d/dR ((1/R) dpsi/dR) + d^2 psi/dZ^2
    !!
    !! for a unit toroidal current filament at (R', Z').  With the weak form
    !! a(psi, v) = int grad psi . grad v / R the kernel solves
    !!
    !!     div ((1/R) grad G) = -delta(R - R') delta(Z - Z')      (G of x, x')
    !!
    !! and is
    !!
    !!     G = (1/(2 pi)) sqrt(R R') ((2 - k^2) K(k) - 2 E(k)) / k
    !!       = s P(m) / (4 pi),    P(m) = (2 - m) K(m) - 2 E(m),
    !!
    !! with m = k^2 = 4 R R'/s^2 and s^2 = (R + R')^2 + (Z - Z')^2.  G carries
    !! no mu0: the poloidal flux per radian (psi = R A_phi, Wb/rad) of a
    !! current I in the filament is psi = mu0 I G.  The mutual inductance of
    !! two coaxial loops is M = 2 pi mu0 G.  G vanishes on the axis R = 0 and
    !! G(x, x') = G(x', x).
    !!
    !! Accuracy.  Derivatives are evaluated in the variable that is regular
    !! at the evaluation point: the complement q = d^2/s^2 (d = distance of
    !! the two points, formed without cancellation) for q < 1/2, and m for
    !! q >= 1/2, where P is summed from its power series
    !! P = (pi/2) sum_{n>=2} 4 n (n - 1) a_n m^n/(2 n - 1)^2,
    !! a_n = ((1/2)_n/n!)^2, so that the far field (P ~ pi m^2/16) keeps full
    !! relative accuracy.  Derivatives follow from dK/dm = (E - qK)/(2 m q),
    !! dP/dm = (E - qK)/(2 q) and the Legendre ODE for K''.
    !!
    !! Point order.  Gradients and Hessians are with respect to the four
    !! coordinates (R, Z, R', Z') in this order.  The kernel is singular at
    !! coincident points (G -> +Inf); R and R' must be positive or zero.
    !!
    !! Log-singular split.  Near coincidence G = -ln(dist) * coeff + smooth,
    !! dist = sqrt((R - R')^2 + (Z - Z')^2); coeff -> R/(2 pi) and both
    !! coeff and smooth are analytic in the coordinates, so singular
    !! quadrature integrates ln(dist) analytically.  The split is
    !! intended for q <= 1/2; farther away use G itself.
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre
    use fortnum_special_elliptic, only: ellip_k_from_complement, &
        ellip_e_from_complement, ellip_k_log_a, ellip_k_log_b, ellip_e_log_a, &
        ellip_e_log_b, ellip_k_derivative, ellip_e_derivative
    implicit none
    private

    real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp
    real(dp), parameter :: ln4 = 1.38629436111989061883446424291635313_dp
    !> Vacuum permeability used by the coil flux (4 pi 1e-7, FreeGS value).
    real(dp), parameter, public :: gs_mu0 = 4.0e-7_dp*pi

    !> Rectangular coil cross-section: centre (r, z), full sizes dr, dz and
    !> the number of turns carrying the circuit current.
    type, public :: gs_rect_coil_t
        real(dp) :: r = 1.0_dp
        real(dp) :: z = 0.0_dp
        real(dp) :: dr = 0.0_dp
        real(dp) :: dz = 0.0_dp
        real(dp) :: turns = 1.0_dp
    end type gs_rect_coil_t

    public :: gs_green
    public :: gs_green_gradient
    public :: gs_green_hessian
    public :: gs_delta_star
    public :: gs_green_conormal
    public :: gs_green_log_split
    public :: gs_green_log_split_gradient
    public :: gs_green_jvp, gs_green_vjp
    public :: gs_green_gradient_jvp, gs_green_gradient_vjp
    public :: gs_rect_coil_flux
    public :: gs_rect_coil_flux_gradient
    public :: gs_rect_coil_flux_converged
    public :: gs_rect_coil_flux_jvp, gs_rect_coil_flux_vjp

contains

    pure function gs_green(r, z, rp, zp) result(g)
        !! G(R, Z; R', Z') (no mu0).
        real(dp), intent(in) :: r, z, rp, zp
        real(dp) :: g
        real(dp) :: h, h1(3), h2(3, 3)

        call kernel_jet(r, rp, z - zp, 0, h, h1, h2)
        g = h/(4.0_dp*pi)
    end function gs_green

    pure subroutine gs_green_gradient(r, z, rp, zp, g, grad)
        !! G and its gradient (dG/dR, dG/dZ, dG/dR', dG/dZ').
        real(dp), intent(in) :: r, z, rp, zp
        real(dp), intent(out) :: g, grad(4)
        real(dp) :: h, h1(3), h2(3, 3)

        call kernel_jet(r, rp, z - zp, 1, h, h1, h2)
        g = h/(4.0_dp*pi)
        grad = gradient_to_points(h1)/(4.0_dp*pi)
    end subroutine gs_green_gradient

    pure subroutine gs_green_hessian(r, z, rp, zp, g, grad, hess)
        !! G, gradient and the symmetric Hessian in (R, Z, R', Z').
        real(dp), intent(in) :: r, z, rp, zp
        real(dp), intent(out) :: g, grad(4), hess(4, 4)
        real(dp) :: h, h1(3), h2(3, 3)

        call kernel_jet(r, rp, z - zp, 2, h, h1, h2)
        g = h/(4.0_dp*pi)
        grad = gradient_to_points(h1)/(4.0_dp*pi)
        hess = hessian_to_points(h2)/(4.0_dp*pi)
    end subroutine gs_green_hessian

    pure subroutine gs_delta_star(r, z, rp, zp, field, source)
        !! Delta* G with respect to the field point (R, Z) and with respect
        !! to the source point (R', Z'); both vanish away from coincidence.
        real(dp), intent(in) :: r, z, rp, zp
        real(dp), intent(out) :: field, source
        real(dp) :: g, grad(4), hess(4, 4)

        call gs_green_hessian(r, z, rp, zp, g, grad, hess)
        field = hess(1, 1) - grad(1)/r + hess(2, 2)
        source = hess(3, 3) - grad(3)/rp + hess(4, 4)
    end subroutine gs_delta_star

    pure subroutine gs_green_conormal(r, z, rp, zp, n, np, d_field, d_source, &
            d_both)
        !! Conormal derivatives with the 1/R weighting of Delta*:
        !!   d_field  = (1/R)  n  . grad_x  G   (normal n at the field point),
        !!   d_source = (1/R') n' . grad_x' G   (normal np at the source),
        !!   d_both   = (1/(R R')) n . (d^2 G/dx dx') . n'.
        !! n and np are used as given (unit normals expected).
        real(dp), intent(in) :: r, z, rp, zp, n(2), np(2)
        real(dp), intent(out) :: d_field, d_source, d_both
        real(dp) :: g, grad(4), hess(4, 4)

        call gs_green_hessian(r, z, rp, zp, g, grad, hess)
        d_field = (n(1)*grad(1) + n(2)*grad(2))/r
        d_source = (np(1)*grad(3) + np(2)*grad(4))/rp
        d_both = (n(1)*(hess(1, 3)*np(1) + hess(1, 4)*np(2)) &
            + n(2)*(hess(2, 3)*np(1) + hess(2, 4)*np(2)))/(r*rp)
    end subroutine gs_green_conormal

    pure subroutine gs_green_log_split(r, z, rp, zp, coeff, smooth)
        !! G = -ln(dist) * coeff + smooth, see the module header.
        real(dp), intent(in) :: r, z, rp, zp
        real(dp), intent(out) :: coeff, smooth
        real(dp) :: dcoeff(4), dsmooth(4)

        call gs_green_log_split_gradient(r, z, rp, zp, coeff, dcoeff, smooth, &
            dsmooth)
    end subroutine gs_green_log_split

    pure subroutine gs_green_log_split_gradient(r, z, rp, zp, coeff, dcoeff, &
            smooth, dsmooth)
        !! Split pieces and their gradients in (R, Z, R', Z').  Valid at
        !! coincident points, where coeff = R/(2 pi).
        real(dp), intent(in) :: r, z, rp, zp
        real(dp), intent(out) :: coeff, dcoeff(4), smooth, dsmooth(4)
        real(dp) :: zeta, dd, uu, d_a(3), u_a(3), s, s_a(3), q, q_a(3)
        real(dp) :: pa, pa1, pb, pb1, ln_s, c_t(3), g_t(3)
        integer :: a

        zeta = z - zp
        dd = (r - rp)**2 + zeta**2
        uu = (r + rp)**2 + zeta**2
        d_a = [2.0_dp*(r - rp), -2.0_dp*(r - rp), 2.0_dp*zeta]
        u_a = [2.0_dp*(r + rp), 2.0_dp*(r + rp), 2.0_dp*zeta]
        s = sqrt(uu)
        s_a = u_a/(2.0_dp*s)
        q = dd/uu
        q_a = (d_a - q*u_a)/uu
        call split_pieces(q, pa, pa1, pb, pb1)
        ln_s = log(s)
        coeff = s*pa/(4.0_dp*pi)
        smooth = s*(ln_s*pa + pb)/(4.0_dp*pi)
        do a = 1, 3
            c_t(a) = (s_a(a)*pa + s*pa1*q_a(a))/(4.0_dp*pi)
            g_t(a) = (s_a(a)*(ln_s + 1.0_dp)*pa + s*ln_s*pa1*q_a(a) &
                + s_a(a)*pb + s*pb1*q_a(a))/(4.0_dp*pi)
        end do
        dcoeff = gradient_to_points(c_t)
        dsmooth = gradient_to_points(g_t)
    end subroutine gs_green_log_split_gradient

    pure subroutine gs_green_jvp(r, z, rp, zp, tangent, g, dg)
        !! Forward product dg = grad G . tangent, tangent in (R, Z, R', Z').
        real(dp), intent(in) :: r, z, rp, zp, tangent(4)
        real(dp), intent(out) :: g, dg
        real(dp) :: grad(4)

        call gs_green_gradient(r, z, rp, zp, g, grad)
        dg = dot_product(grad, tangent)
    end subroutine gs_green_jvp

    pure subroutine gs_green_vjp(r, z, rp, zp, gbar, bar)
        !! Reverse product bar = gbar * grad G for the cotangent gbar of G.
        real(dp), intent(in) :: r, z, rp, zp, gbar
        real(dp), intent(out) :: bar(4)
        real(dp) :: g, grad(4)

        call gs_green_gradient(r, z, rp, zp, g, grad)
        bar = gbar*grad
    end subroutine gs_green_vjp

    pure subroutine gs_green_gradient_jvp(r, z, rp, zp, tangent, grad, dgrad)
        !! Forward product of the gradient: dgrad = Hessian . tangent.
        real(dp), intent(in) :: r, z, rp, zp, tangent(4)
        real(dp), intent(out) :: grad(4), dgrad(4)
        real(dp) :: g, hess(4, 4)

        call gs_green_hessian(r, z, rp, zp, g, grad, hess)
        dgrad = matmul(hess, tangent)
    end subroutine gs_green_gradient_jvp

    pure subroutine gs_green_gradient_vjp(r, z, rp, zp, gradbar, bar)
        !! Reverse product of the gradient: bar = Hessian . gradbar (the
        !! Hessian is symmetric).
        real(dp), intent(in) :: r, z, rp, zp, gradbar(4)
        real(dp), intent(out) :: bar(4)
        real(dp) :: g, grad(4), hess(4, 4)

        call gs_green_hessian(r, z, rp, zp, g, grad, hess)
        bar = matmul(transpose(hess), gradbar)
    end subroutine gs_green_gradient_vjp

    function gs_rect_coil_flux(coil, r, z, nr, nz) result(psi)
        !! Flux per radian (Wb/rad) at (r, z) per ampere of the circuit
        !! current: mu0 * turns * (cross-section average of G), by the
        !! tensor Gauss-Legendre rule with nr x nz points.  The field point
        !! must be outside the rectangle.
        type(gs_rect_coil_t), intent(in) :: coil
        real(dp), intent(in) :: r, z
        integer, intent(in) :: nr, nz
        real(dp) :: psi
        real(dp) :: xr(nr), wr(nr), xz(nz), wz(nz)
        integer :: i, j

        call gauss_legendre(nr, xr, wr)
        call gauss_legendre(nz, xz, wz)
        psi = 0.0_dp
        do j = 1, nz
            do i = 1, nr
                psi = psi + wr(i)*wz(j)*gs_green(r, z, &
                    coil%r + 0.5_dp*coil%dr*xr(i), coil%z + 0.5_dp*coil%dz*xz(j))
            end do
        end do
        psi = 0.25_dp*gs_mu0*coil%turns*psi
    end function gs_rect_coil_flux

    subroutine gs_rect_coil_flux_gradient(coil, r, z, nr, nz, psi, grad)
        !! Flux and its gradient with respect to
        !! (r, z, coil%r, coil%z, coil%dr, coil%dz); the flux is linear in
        !! coil%turns.  Same quadrature as gs_rect_coil_flux.
        type(gs_rect_coil_t), intent(in) :: coil
        real(dp), intent(in) :: r, z
        integer, intent(in) :: nr, nz
        real(dp), intent(out) :: psi, grad(6)
        real(dp) :: xr(nr), wr(nr), xz(nz), wz(nz)
        real(dp) :: g, gg(4), w
        integer :: i, j

        call gauss_legendre(nr, xr, wr)
        call gauss_legendre(nz, xz, wz)
        psi = 0.0_dp
        grad = 0.0_dp
        do j = 1, nz
            do i = 1, nr
                w = wr(i)*wz(j)
                call gs_green_gradient(r, z, coil%r + 0.5_dp*coil%dr*xr(i), &
                    coil%z + 0.5_dp*coil%dz*xz(j), g, gg)
                psi = psi + w*g
                grad(1) = grad(1) + w*gg(1)
                grad(2) = grad(2) + w*gg(2)
                grad(3) = grad(3) + w*gg(3)
                grad(4) = grad(4) + w*gg(4)
                grad(5) = grad(5) + w*gg(3)*0.5_dp*xr(i)
                grad(6) = grad(6) + w*gg(4)*0.5_dp*xz(j)
            end do
        end do
        psi = 0.25_dp*gs_mu0*coil%turns*psi
        grad = 0.25_dp*gs_mu0*coil%turns*grad
    end subroutine gs_rect_coil_flux_gradient

    subroutine gs_rect_coil_flux_converged(coil, r, z, tol, psi, order, status)
        !! Coil flux with the Gauss order (same in both directions) raised
        !! through 2, 3, 4, 6, 8, 12, ..., 64 until two successive values
        !! agree to the relative tolerance tol.  status: 0 converged,
        !! 1 not converged at order 64 (the last value is returned),
        !! 2 the field point lies in or on the rectangle.
        type(gs_rect_coil_t), intent(in) :: coil
        real(dp), intent(in) :: r, z, tol
        real(dp), intent(out) :: psi
        integer, intent(out) :: order, status
        integer, parameter :: orders(10) = [2, 3, 4, 6, 8, 12, 16, 24, 32, 64]
        real(dp) :: previous
        integer :: k

        order = 0
        psi = 0.0_dp
        status = 2
        if (abs(r - coil%r) <= 0.5_dp*abs(coil%dr) .and. &
            abs(z - coil%z) <= 0.5_dp*abs(coil%dz)) return
        status = 1
        previous = gs_rect_coil_flux(coil, r, z, orders(1), orders(1))
        order = orders(1)
        psi = previous
        do k = 2, size(orders)
            psi = gs_rect_coil_flux(coil, r, z, orders(k), orders(k))
            order = orders(k)
            if (abs(psi - previous) <= tol*abs(psi)) then
                status = 0
                return
            end if
            previous = psi
        end do
    end subroutine gs_rect_coil_flux_converged

    subroutine gs_rect_coil_flux_jvp(coil, r, z, nr, nz, tangent, psi, dpsi)
        !! Forward product dpsi = grad psi . tangent, with tangent in
        !! (r, z, coil%r, coil%z, coil%dr, coil%dz).
        type(gs_rect_coil_t), intent(in) :: coil
        real(dp), intent(in) :: r, z, tangent(6)
        integer, intent(in) :: nr, nz
        real(dp), intent(out) :: psi, dpsi
        real(dp) :: grad(6)

        call gs_rect_coil_flux_gradient(coil, r, z, nr, nz, psi, grad)
        dpsi = dot_product(grad, tangent)
    end subroutine gs_rect_coil_flux_jvp

    subroutine gs_rect_coil_flux_vjp(coil, r, z, nr, nz, psibar, bar)
        !! Reverse product bar = psibar * grad psi in the same ordering.
        type(gs_rect_coil_t), intent(in) :: coil
        real(dp), intent(in) :: r, z, psibar
        integer, intent(in) :: nr, nz
        real(dp), intent(out) :: bar(6)
        real(dp) :: psi

        call gs_rect_coil_flux_gradient(coil, r, z, nr, nz, psi, bar)
        bar = psibar*bar
    end subroutine gs_rect_coil_flux_vjp

    ! ------------------------------------------------------------------
    ! Private kernel machinery
    ! ------------------------------------------------------------------

    pure function gradient_to_points(h1) result(grad)
        !! Map a gradient in t = (R, R', zeta = Z - Z') to (R, Z, R', Z').
        real(dp), intent(in) :: h1(3)
        real(dp) :: grad(4)

        grad = [h1(1), h1(3), h1(2), -h1(3)]
    end function gradient_to_points

    pure function hessian_to_points(h2) result(hess)
        !! Map a Hessian in t to (R, Z, R', Z') as J^T H J.
        real(dp), intent(in) :: h2(3, 3)
        real(dp) :: hess(4, 4)
        real(dp) :: jac(3, 4)

        jac = 0.0_dp
        jac(1, 1) = 1.0_dp
        jac(2, 3) = 1.0_dp
        jac(3, 2) = 1.0_dp
        jac(3, 4) = -1.0_dp
        hess = matmul(transpose(jac), matmul(h2, jac))
    end function hessian_to_points

    pure subroutine kernel_jet(r, rp, zeta, order, h, h1, h2)
        !! h = s P = 4 pi G and its derivatives (up to `order`) with
        !! respect to t = (R, R', zeta).  Entries beyond `order` are zero.
        real(dp), intent(in) :: r, rp, zeta
        integer, intent(in) :: order
        real(dp), intent(out) :: h, h1(3), h2(3, 3)
        real(dp) :: dd, uu, nn, s, q, m
        real(dp) :: d_a(3), d_ab(3, 3), u_a(3), u_ab(3, 3), n_a(3), n_ab(3, 3)
        real(dp) :: num_a(3), num_ab(3, 3), y, y_a(3), y_ab(3, 3)
        real(dp) :: s_a(3), s_ab(3, 3), f, f1, f2
        integer :: a, b

        h1 = 0.0_dp
        h2 = 0.0_dp
        if (.not. (r >= 0.0_dp .and. rp >= 0.0_dp)) then
            h = ieee_nan()
            return
        end if
        dd = (r - rp)**2 + zeta**2
        uu = (r + rp)**2 + zeta**2
        nn = 4.0_dp*r*rp
        s = sqrt(uu)
        q = dd/uu
        m = nn/uu

        d_a = [2.0_dp*(r - rp), -2.0_dp*(r - rp), 2.0_dp*zeta]
        u_a = [2.0_dp*(r + rp), 2.0_dp*(r + rp), 2.0_dp*zeta]
        n_a = [4.0_dp*rp, 4.0_dp*r, 0.0_dp]
        d_ab = 0.0_dp
        d_ab(1, 1) = 2.0_dp
        d_ab(1, 2) = -2.0_dp
        d_ab(2, 1) = -2.0_dp
        d_ab(2, 2) = 2.0_dp
        d_ab(3, 3) = 2.0_dp
        u_ab = 0.0_dp
        u_ab(1:2, 1:2) = 2.0_dp
        u_ab(3, 3) = 2.0_dp
        n_ab = 0.0_dp
        n_ab(1, 2) = 4.0_dp
        n_ab(2, 1) = 4.0_dp

        if (q < 0.5_dp) then
            y = q
            num_a = d_a
            num_ab = d_ab
            call p_near(q, m, f, f1, f2)
        else
            y = m
            num_a = n_a
            num_ab = n_ab
            call p_far(m, f, f1, f2)
        end if
        h = s*f
        if (order < 1) return

        y_a = (num_a - y*u_a)/uu
        s_a = u_a/(2.0_dp*s)
        h1 = s_a*f + s*f1*y_a
        if (order < 2) return

        do b = 1, 3
            do a = 1, 3
                y_ab(a, b) = (num_ab(a, b) - y_a(b)*u_a(a) - y_a(a)*u_a(b) &
                    - y*u_ab(a, b))/uu
                s_ab(a, b) = u_ab(a, b)/(2.0_dp*s) &
                    - u_a(a)*u_a(b)/(4.0_dp*s**3)
            end do
        end do
        do b = 1, 3
            do a = 1, 3
                h2(a, b) = s_ab(a, b)*f + (s_a(a)*y_a(b) + s_a(b)*y_a(a))*f1 &
                    + s*f2*y_a(a)*y_a(b) + s*f1*y_ab(a, b)
            end do
        end do
    end subroutine kernel_jet

    pure subroutine p_near(q, m, f, f1, f2)
        !! P = (2 - m) K - 2 E in the complement variable y = q (q < 1/2):
        !! f = P, f1 = dP/dq, f2 = d^2 P/dq^2.  m = 1 - q is passed in
        !! (formed from the coordinates) and enters only through factors
        !! that are bounded away from cancellation.
        real(dp), intent(in) :: q, m
        real(dp), intent(out) :: f, f1, f2
        real(dp) :: kk, ee, kp, kpp

        kk = ellip_k_from_complement(q)
        ee = ellip_e_from_complement(q)
        kp = (ee - q*kk)/(2.0_dp*m*q)
        kpp = (0.25_dp*kk + (2.0_dp*m - 1.0_dp)*kp)/(m*q)
        f = (1.0_dp + q)*kk - 2.0_dp*ee
        f1 = -m*kp
        f2 = kp + m*kpp
    end subroutine p_near

    pure subroutine p_far(m, f, f1, f2)
        !! P in the variable y = m (m <= 1/2) from its power series.
        real(dp), intent(in) :: m
        real(dp), intent(out) :: f, f1, f2
        real(dp) :: a, c, mp, t2
        integer :: n

        ! c_n = 4 n (n - 1) a_n/(2 n - 1)^2, P = (pi/2) sum c_n m^n.
        a = 0.25_dp
        mp = 1.0_dp
        f = 0.0_dp
        f1 = 0.0_dp
        f2 = 0.0_dp
        do n = 2, 400
            a = a*((2.0_dp*n - 1.0_dp)/(2.0_dp*n))**2
            c = 4.0_dp*n*(n - 1.0_dp)*a/(2.0_dp*n - 1.0_dp)**2
            f = f + c*mp*m*m
            f1 = f1 + n*c*mp*m
            t2 = n*(n - 1.0_dp)*c*mp
            f2 = f2 + t2
            mp = mp*m
            if (n > 3 .and. t2 <= 1.0e-18_dp*f2) exit
        end do
        f = 0.5_dp*pi*f
        f1 = 0.5_dp*pi*f1
        f2 = 0.5_dp*pi*f2
    end subroutine p_far

    pure subroutine split_pieces(q, pa, pa1, pb, pb1)
        !! P = -1/2 ln(q) P_A(q) + P_B(q) with P_A = (1 + q) A_K - 2 A_E and
        !! P_B = (1 + q) B_K - 2 B_E, and their q-derivatives pa1, pb1.
        real(dp), intent(in) :: q
        real(dp), intent(out) :: pa, pa1, pb, pb1
        real(dp) :: ak, bk, ae, be, dk, de, akp, aep, bkp, bep, lq

        ak = ellip_k_log_a(q)
        bk = ellip_k_log_b(q)
        ae = ellip_e_log_a(q)
        be = ellip_e_log_b(q)
        dk = ellip_k_derivative(q)
        de = ellip_e_derivative(q)
        akp = 2.0_dp/pi*dk
        aep = 2.0_dp/pi*(dk - de)
        if (q < 0.5_dp) then
            call b_series_derivative(q, bkp, bep)
        else
            lq = log(q)
            bkp = -ellip_k_derivative(1.0_dp - q) + ak/(2.0_dp*q) + 0.5_dp*lq*akp
            bep = -ellip_e_derivative(1.0_dp - q) + ae/(2.0_dp*q) + 0.5_dp*lq*aep
        end if
        pa = (1.0_dp + q)*ak - 2.0_dp*ae
        pb = (1.0_dp + q)*bk - 2.0_dp*be
        pa1 = ak + (1.0_dp + q)*akp - 2.0_dp*aep
        pb1 = bk + (1.0_dp + q)*bkp - 2.0_dp*bep
    end subroutine split_pieces

    pure subroutine b_series_derivative(q, bkp, bep)
        !! dB_K/dq and dB_E/dq for q < 1/2 from the DLMF 19.12.1-2 series:
        !! B_K = sum a_n d_n q^n, B_E = 1 + (1/2) sum e_n g_n q^(n+1) with
        !! a_n = ((1/2)_n/n!)^2, e_n = (1/2)_n (3/2)_n/((2)_n n!),
        !! d_n = psi(n+1) - psi(n+1/2), g_n = d_n - 1/((2n+1)(2n+2)).
        real(dp), intent(in) :: q
        real(dp), intent(out) :: bkp, bep
        real(dp) :: a, d, e, qn, tk, te, x
        integer :: n

        ! B_K' = sum_{n>=1} n a_n d_n q^(n-1)
        a = 1.0_dp
        d = ln4
        qn = 1.0_dp
        bkp = 0.0_dp
        do n = 1, 200
            a = a*((2.0_dp*n - 1.0_dp)/(2.0_dp*n))**2
            d = d - 1.0_dp/(n*(2.0_dp*n - 1.0_dp))
            tk = n*a*d*qn
            bkp = bkp + tk
            qn = qn*q
            if (abs(tk) <= 1.0e-18_dp*abs(bkp)) exit
        end do
        ! B_E' = (1/2) sum_{n>=0} (n + 1) e_n g_n q^n
        e = 1.0_dp
        d = ln4
        qn = 1.0_dp
        bep = 0.0_dp
        do n = 0, 200
            x = real(n, dp)
            te = 0.5_dp*(x + 1.0_dp)*e*(d - 1.0_dp/((2.0_dp*x + 1.0_dp)* &
                (2.0_dp*x + 2.0_dp)))*qn
            bep = bep + te
            e = e*(x + 0.5_dp)*(x + 1.5_dp)/((x + 2.0_dp)*(x + 1.0_dp))
            d = d - 1.0_dp/((x + 1.0_dp)*(2.0_dp*x + 1.0_dp))
            qn = qn*q
            if (n > 0 .and. abs(te) <= 1.0e-18_dp*abs(bep)) exit
        end do
    end subroutine b_series_derivative

    pure function ieee_nan() result(x)
        use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
        real(dp) :: x

        x = ieee_value(1.0_dp, ieee_quiet_nan)
    end function ieee_nan

end module fortfem_axisymmetric_gs_kernel
