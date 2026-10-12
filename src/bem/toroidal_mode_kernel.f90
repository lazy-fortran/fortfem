module fortfem_toroidal_mode_kernel
    !! Toroidal-mode Green's function of the Laplacian on the meridional half
    !! plane (R > 0) and the coil flux it gives per mode number n.
    !!
    !! The operator is the scalar form of the paper's per-mode reduction
    !! (Albert, Lainer and Biro, Comput. Phys. Commun. 277 (2022) 108401,
    !! doi:10.1016/j.cpc.2022.108401): grad U -> grad_t U + i n U e^3
    !! (Eq. (34)), div U -> div_t u + i n U^3 (Eq. (35)), with sqrt(g) g^{tt}
    !! = R and sqrt(g) g^{33} = 1/R in cylindrical coordinates (Eqs. (7), (44)),
    !! so -div(R grad u) + n^2 u / R is the e^{i n phi} component, with the
    !! Fourier convention of Eqs. (31)-(32).
    !!
    !! Normalisation (derived, and checked against a direct phi-quadrature in
    !! the tests).  The Laplace Green's function 1/(4 pi |x - x'|) has the
    !! toroidal expansion (Heine; Cohl and Tohline 1999)
    !!
    !!     1/|x - x'| = 1/(pi sqrt(R R')) sum_n Q_{|n|-1/2}(chi) e^{i n (phi - phi')},
    !!     chi = (R^2 + R'^2 + (Z - Z')^2)/(2 R R').
    !!
    !! A field e^{i n phi} u(R, Z) is therefore produced by a source
    !! e^{i n phi'} f(R', Z') through
    !!
    !!     u(R, Z) = int G_n(R, Z; R', Z') f(R', Z') R' dR' dZ',
    !!     G_n = int_0^{2 pi} e^{-i n D} dD/(4 pi r(D))
    !!         = Q_{|n|-1/2}(chi)/(2 pi sqrt(R R')),
    !!
    !! (D = phi - phi'; the factor 2 pi of the phi' integral turns the 1/pi of
    !! the expansion into 1/(2 pi)).  G_n is real, even in n, symmetric in
    !! (R, Z) <-> (R', Z'), and solves, with respect to the measure R' dR' dZ',
    !!
    !!     -[ (1/R) d_R (R d_R) + d_Z^2 - n^2/R^2 ] G_n = delta(R - R') delta(Z - Z')/R,
    !!
    !! equivalently -div(R grad G_n) + (n^2/R) G_n = delta_2D.  Near coincidence
    !! Q_{|n|-1/2}(chi) ~ -ln(d) + ... with d = |x - x'| in the meridian plane
    !! (chi - 1 = d^2/(2 R R')), so G_n ~ -ln(d)/(2 pi R), the 2D kernel of the
    !! weight R.
    !!
    !! Relation to the axisymmetric Grad-Shafranov kernel.  Q_{-1/2}(chi) =
    !! k K(k), k^2 = 2/(chi + 1) = 4 R R'/s^2, s^2 = (R + R')^2 + (Z - Z')^2,
    !! hence
    !!
    !!     G_0 = K(k)/(pi s)           (scalar potential of a ring source),
    !!
    !! and Q_{1/2}(chi) = ((2 - k^2) K - 2 E)/k, so the Green's function of
    !! fortfem_axisymmetric_gs_kernel (the poloidal flux psi = R A_phi of a ring
    !! current, whose A_phi is the n = 1 vector mode) is exactly
    !!
    !!     G_GS = R R' G_1.
    !!
    !! That kernel is the |n| = 1 member of this family, not the n = 0 one.
    !!
    !! Accuracy.  Q_{n-1/2} comes from fortnum (toroidal_q_sequence).  Close
    !! to coincidence t = (chi - 1)/2 = d^2/(4 R R') < t_s the logarithmic
    !! decomposition Q = L(t) p(t) + r(t), L = -ln(t)/2, with p = P_{n-1/2} and
    !! r analytic, is summed from its Taylor series in t.  Its coefficients
    !! follow from the Legendre equation t (1 + t) y'' + (1 + 2 t) y' = lam y,
    !! lam = n^2 - 1/4:
    !!     c_0 = 1,  c_{k+1} = c_k (lam - k (k + 1))/(k + 1)^2,
    !!     e_0 = -gamma - psi(n + 1/2) = 2 ln 2 - sum_{j<=n} 2/(2 j - 1),
    !!     (k + 1)^2 e_{k+1} = (k + 1) c_{k+1} + (k + 1/2) c_k
    !!                         + (lam - k (k + 1)) e_k,
    !! p = sum c_k t^k, r = sum e_k t^k.  This keeps Q and dQ/dchi accurate
    !! down to t = 1e-20 where the difference formulas of the recurrences lose
    !! all digits.  Second derivatives use the Legendre equation.
    !!
    !! Point order.  Gradients and Hessians are with respect to
    !! (R, Z, R', Z') in this order.  R and R' must be positive; the kernel is
    !! singular at coincidence (+Inf).  A negative n is treated as |n|.
    !!
    !! Log split.  G_n = -ln(d) coeff + smooth with
    !!     coeff = p(t)/(2 pi sqrt(R R')),
    !!     smooth = (ln(4 R R') p(t)/2 + r(t))/(2 pi sqrt(R R')),
    !! both analytic in the coordinates, coeff -> 1/(2 pi R) at coincidence.
    !!
    !! Coil flux.  A closed filament at (Rc, Zc) carrying the mode e^{i n phi'}
    !! with unit amplitude per unit current produces u = Rc G_n(R, Z; Rc, Zc);
    !! a rectangular winding is the cross-section average of R' G_n.  Neither
    !! carries mu0.  For n = 1, mu0 R times this flux is the flux per radian
    !! of fortfem_axisymmetric_gs_kernel (gs_rect_coil_flux).
    use, intrinsic :: ieee_arithmetic, only: ieee_positive_inf, ieee_quiet_nan, &
        ieee_value
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre
    use fortnum_special_toroidal, only: toroidal_q_sequence
    implicit none
    private

    real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp
    real(dp), parameter :: ln2 = 0.693147180559945309417232121458176568_dp
    !> Largest t for which the log split is evaluated at all.
    real(dp), parameter :: split_t_max = 0.5_dp

    !> Rectangular coil cross-section: centre (r, z), full sizes dr, dz and
    !> the number of turns carrying the circuit current.
    type, public :: tm_rect_coil_t
        real(dp) :: r = 1.0_dp
        real(dp) :: z = 0.0_dp
        real(dp) :: dr = 0.0_dp
        real(dp) :: dz = 0.0_dp
        real(dp) :: turns = 1.0_dp
    end type tm_rect_coil_t

    public :: tm_green
    public :: tm_green_modes
    public :: tm_green_gradient
    public :: tm_green_hessian
    public :: tm_green_operator
    public :: tm_green_conormal
    public :: tm_green_log_split
    public :: tm_green_log_split_gradient
    public :: tm_green_modes_log_split
    public :: tm_log_split_limit
    public :: tm_filament_flux
    public :: tm_filament_flux_gradient
    public :: tm_rect_coil_flux
    public :: tm_rect_coil_flux_modes
    public :: tm_rect_coil_flux_gradient
    public :: tm_rect_coil_flux_converged
    public :: tm_rect_coil_flux_table

contains

    ! ------------------------------------------------------------------
    ! Green's function
    ! ------------------------------------------------------------------

    pure function tm_green(n, r, z, rp, zp) result(g)
        !! G_n(R, Z; R', Z') for the measure R' dR' dZ'.
        integer, intent(in) :: n
        real(dp), intent(in) :: r, z, rp, zp
        real(dp) :: g
        real(dp) :: h, h1(3), h2(3, 3)

        call mode_jet(abs(n), 0, r, rp, z - zp, h, h1, h2)
        g = h
    end function tm_green

    pure subroutine tm_green_modes(nmax, r, z, rp, zp, g)
        !! G_m for m = 0..nmax at one point pair.
        integer, intent(in) :: nmax
        real(dp), intent(in) :: r, z, rp, zp
        real(dp), intent(out) :: g(0:nmax)
        real(dp) :: q(0:nmax), qx(0:nmax)

        if (.not. (r > 0.0_dp .and. rp > 0.0_dp)) then
            g = ieee_value(1.0_dp, ieee_quiet_nan)
            return
        end if
        call q_modes(nmax, r, rp, z - zp, q, qx)
        g = weight(r, rp)*q
    end subroutine tm_green_modes

    pure subroutine tm_green_gradient(n, r, z, rp, zp, g, grad)
        !! G_n and its gradient (dG/dR, dG/dZ, dG/dR', dG/dZ').
        integer, intent(in) :: n
        real(dp), intent(in) :: r, z, rp, zp
        real(dp), intent(out) :: g, grad(4)
        real(dp) :: h, h1(3), h2(3, 3)

        call mode_jet(abs(n), 1, r, rp, z - zp, h, h1, h2)
        g = h
        grad = gradient_to_points(h1)
    end subroutine tm_green_gradient

    pure subroutine tm_green_hessian(n, r, z, rp, zp, g, grad, hess)
        !! G_n, gradient and the symmetric Hessian in (R, Z, R', Z').
        integer, intent(in) :: n
        real(dp), intent(in) :: r, z, rp, zp
        real(dp), intent(out) :: g, grad(4), hess(4, 4)
        real(dp) :: h, h1(3), h2(3, 3)

        call mode_jet(abs(n), 2, r, rp, z - zp, h, h1, h2)
        g = h
        grad = gradient_to_points(h1)
        hess = hessian_to_points(h2)
    end subroutine tm_green_hessian

    pure subroutine tm_green_operator(n, r, z, rp, zp, field, source)
        !! The n-mode Laplacian
        !!     G_RR + G_R/R + G_ZZ - n^2 G/R^2
        !! with respect to the field point (R, Z) and, in the same form, with
        !! respect to the source point (R', Z'); both vanish away from
        !! coincidence.
        integer, intent(in) :: n
        real(dp), intent(in) :: r, z, rp, zp
        real(dp), intent(out) :: field, source
        real(dp) :: g, grad(4), hess(4, 4)
        real(dp) :: n2

        call tm_green_hessian(n, r, z, rp, zp, g, grad, hess)
        n2 = real(n, dp)**2
        field = hess(1, 1) + grad(1)/r + hess(2, 2) - n2*g/(r*r)
        source = hess(3, 3) + grad(3)/rp + hess(4, 4) - n2*g/(rp*rp)
    end subroutine tm_green_operator

    pure subroutine tm_green_conormal(n, r, z, rp, zp, nrm, nrmp, d_field, &
            d_source, d_both)
        !! Conormal derivatives with the weight R of the operator
        !! -div(R grad) + n^2/R:
        !!   d_field  = R  nrm  . grad_x  G_n   (normal nrm at the field point),
        !!   d_source = R' nrmp . grad_x' G_n   (normal nrmp at the source),
        !!   d_both   = R R' nrm . (d^2 G_n/dx dx') . nrmp.
        !! The normals are used as given (unit normals expected).
        integer, intent(in) :: n
        real(dp), intent(in) :: r, z, rp, zp, nrm(2), nrmp(2)
        real(dp), intent(out) :: d_field, d_source, d_both
        real(dp) :: g, grad(4), hess(4, 4)

        call tm_green_hessian(n, r, z, rp, zp, g, grad, hess)
        d_field = r*(nrm(1)*grad(1) + nrm(2)*grad(2))
        d_source = rp*(nrmp(1)*grad(3) + nrmp(2)*grad(4))
        d_both = r*rp*(nrm(1)*(hess(1, 3)*nrmp(1) + hess(1, 4)*nrmp(2)) &
            + nrm(2)*(hess(2, 3)*nrmp(1) + hess(2, 4)*nrmp(2)))
    end subroutine tm_green_conormal

    ! ------------------------------------------------------------------
    ! Logarithmic split
    ! ------------------------------------------------------------------

    pure function tm_log_split_limit(nmax) result(ts)
        !! Largest t = d^2/(4 R R') up to which the split pieces for modes up
        !! to nmax keep full accuracy (the series alternates and loses about
        !! exp(2 sqrt(lam t)) in digits).  Beyond it use G_n itself.
        integer, intent(in) :: nmax
        real(dp) :: ts

        ts = min(0.2_dp, 2.0_dp/real(nmax + 1, dp)**2)
    end function tm_log_split_limit

    pure subroutine tm_green_log_split(n, r, z, rp, zp, coeff, smooth)
        !! G_n = -ln(dist) coeff + smooth, see the module header.
        integer, intent(in) :: n
        real(dp), intent(in) :: r, z, rp, zp
        real(dp), intent(out) :: coeff, smooth
        real(dp) :: dcoeff(4), dsmooth(4)

        call tm_green_log_split_gradient(n, r, z, rp, zp, coeff, dcoeff, &
            smooth, dsmooth)
    end subroutine tm_green_log_split

    pure subroutine tm_green_log_split_gradient(n, r, z, rp, zp, coeff, &
            dcoeff, smooth, dsmooth)
        !! Split pieces and their gradients in (R, Z, R', Z').  Valid at
        !! coincident points, where coeff = 1/(2 pi R).  The pieces are
        !! accurate for t < tm_log_split_limit(|n|) and defined for t < 1/2.
        integer, intent(in) :: n
        real(dp), intent(in) :: r, z, rp, zp
        real(dp), intent(out) :: coeff, dcoeff(4), smooth, dsmooth(4)
        real(dp) :: zeta, t, p, pt, rr, rt, w, wa(3), lg, lga(3), ta(3)
        real(dp) :: c_t(3), s_t(3), inner
        integer :: a

        zeta = z - zp
        t = ((r - rp)**2 + zeta**2)/(4.0_dp*r*rp)
        call log_split_series(abs(n), t, p, pt, rr, rt)
        w = weight(r, rp)
        wa = [-0.5_dp*w/r, -0.5_dp*w/rp, 0.0_dp]
        lg = 0.5_dp*log(4.0_dp*r*rp)
        lga = [0.5_dp/r, 0.5_dp/rp, 0.0_dp]
        ta = [(r - rp)/(2.0_dp*r*rp) - t/r, -(r - rp)/(2.0_dp*r*rp) - t/rp, &
            zeta/(2.0_dp*r*rp)]
        coeff = w*p
        inner = lg*p + rr
        smooth = w*inner
        do a = 1, 3
            c_t(a) = wa(a)*p + w*pt*ta(a)
            s_t(a) = wa(a)*inner + w*(lga(a)*p + (lg*pt + rt)*ta(a))
        end do
        dcoeff = gradient_to_points(c_t)
        dsmooth = gradient_to_points(s_t)
    end subroutine tm_green_log_split_gradient

    pure subroutine tm_green_modes_log_split(nmax, r, z, rp, zp, coeff, smooth)
        !! Split pieces (values only) of G_m for m = 0..nmax.
        integer, intent(in) :: nmax
        real(dp), intent(in) :: r, z, rp, zp
        real(dp), intent(out) :: coeff(0:nmax), smooth(0:nmax)
        real(dp) :: t, p, pt, rr, rt, w, lg
        integer :: m

        t = ((r - rp)**2 + (z - zp)**2)/(4.0_dp*r*rp)
        w = weight(r, rp)
        lg = 0.5_dp*log(4.0_dp*r*rp)
        do m = 0, nmax
            call log_split_series(m, t, p, pt, rr, rt)
            coeff(m) = w*p
            smooth(m) = w*(lg*p + rr)
        end do
    end subroutine tm_green_modes_log_split

    ! ------------------------------------------------------------------
    ! Coil flux
    ! ------------------------------------------------------------------

    pure function tm_filament_flux(n, rc, zc, r, z) result(psi)
        !! Mode-n potential at (r, z) of a unit-amplitude filament at
        !! (rc, zc): rc G_n(r, z; rc, zc).
        integer, intent(in) :: n
        real(dp), intent(in) :: rc, zc, r, z
        real(dp) :: psi

        psi = rc*tm_green(n, r, z, rc, zc)
    end function tm_filament_flux

    pure subroutine tm_filament_flux_gradient(n, rc, zc, r, z, psi, grad)
        !! Filament flux and its gradient with respect to (r, z, rc, zc).
        integer, intent(in) :: n
        real(dp), intent(in) :: rc, zc, r, z
        real(dp), intent(out) :: psi, grad(4)
        real(dp) :: g, gg(4)

        call tm_green_gradient(n, r, z, rc, zc, g, gg)
        psi = rc*g
        grad = rc*gg
        grad(3) = grad(3) + g
    end subroutine tm_filament_flux_gradient

    function tm_rect_coil_flux(n, coil, r, z, nr, nz) result(psi)
        !! Mode-n potential at (r, z) per unit circuit current:
        !! turns times the cross-section average of R' G_n, by the tensor
        !! Gauss-Legendre rule with nr x nz points.  The field point must be
        !! outside the rectangle.
        integer, intent(in) :: n
        type(tm_rect_coil_t), intent(in) :: coil
        real(dp), intent(in) :: r, z
        integer, intent(in) :: nr, nz
        real(dp) :: psi
        real(dp) :: xr(nr), wr(nr), xz(nz), wz(nz), rq
        integer :: i, j

        call gauss_legendre(nr, xr, wr)
        call gauss_legendre(nz, xz, wz)
        psi = 0.0_dp
        do j = 1, nz
            do i = 1, nr
                rq = coil%r + 0.5_dp*coil%dr*xr(i)
                psi = psi + wr(i)*wz(j)*rq*tm_green(n, r, z, rq, &
                    coil%z + 0.5_dp*coil%dz*xz(j))
            end do
        end do
        psi = 0.25_dp*coil%turns*psi
    end function tm_rect_coil_flux

    subroutine tm_rect_coil_flux_modes(nmax, coil, r, z, nr, nz, psi)
        !! psi(m) = tm_rect_coil_flux(m, ...) for m = 0..nmax with one
        !! evaluation of the Q sequence per quadrature point.
        integer, intent(in) :: nmax
        type(tm_rect_coil_t), intent(in) :: coil
        real(dp), intent(in) :: r, z
        integer, intent(in) :: nr, nz
        real(dp), intent(out) :: psi(0:nmax)
        real(dp) :: xr(nr), wr(nr), xz(nz), wz(nz), rq, g(0:nmax)
        integer :: i, j

        call gauss_legendre(nr, xr, wr)
        call gauss_legendre(nz, xz, wz)
        psi = 0.0_dp
        do j = 1, nz
            do i = 1, nr
                rq = coil%r + 0.5_dp*coil%dr*xr(i)
                call tm_green_modes(nmax, r, z, rq, &
                    coil%z + 0.5_dp*coil%dz*xz(j), g)
                psi = psi + wr(i)*wz(j)*rq*g
            end do
        end do
        psi = 0.25_dp*coil%turns*psi
    end subroutine tm_rect_coil_flux_modes

    subroutine tm_rect_coil_flux_gradient(n, coil, r, z, nr, nz, psi, grad)
        !! Flux and its gradient with respect to
        !! (r, z, coil%r, coil%z, coil%dr, coil%dz); the flux is linear in
        !! coil%turns.  Same quadrature as tm_rect_coil_flux.
        integer, intent(in) :: n
        type(tm_rect_coil_t), intent(in) :: coil
        real(dp), intent(in) :: r, z
        integer, intent(in) :: nr, nz
        real(dp), intent(out) :: psi, grad(6)
        real(dp) :: xr(nr), wr(nr), xz(nz), wz(nz)
        real(dp) :: g, gg(4), w, rq, zq, f_r, f_z
        integer :: i, j

        call gauss_legendre(nr, xr, wr)
        call gauss_legendre(nz, xz, wz)
        psi = 0.0_dp
        grad = 0.0_dp
        do j = 1, nz
            do i = 1, nr
                w = wr(i)*wz(j)
                rq = coil%r + 0.5_dp*coil%dr*xr(i)
                zq = coil%z + 0.5_dp*coil%dz*xz(j)
                call tm_green_gradient(n, r, z, rq, zq, g, gg)
                ! f = R' G_n; df/dR' = G_n + R' dG/dR'
                f_r = g + rq*gg(3)
                f_z = rq*gg(4)
                psi = psi + w*rq*g
                grad(1) = grad(1) + w*rq*gg(1)
                grad(2) = grad(2) + w*rq*gg(2)
                grad(3) = grad(3) + w*f_r
                grad(4) = grad(4) + w*f_z
                grad(5) = grad(5) + w*f_r*0.5_dp*xr(i)
                grad(6) = grad(6) + w*f_z*0.5_dp*xz(j)
            end do
        end do
        psi = 0.25_dp*coil%turns*psi
        grad = 0.25_dp*coil%turns*grad
    end subroutine tm_rect_coil_flux_gradient

    subroutine tm_rect_coil_flux_converged(n, coil, r, z, tol, psi, order, &
            status)
        !! Coil flux with the Gauss order (same in both directions) raised
        !! through 2, 3, 4, 6, 8, 12, ..., 64 until two successive values
        !! agree to the relative tolerance tol.  status: 0 converged,
        !! 1 not converged at order 64 (the last value is returned),
        !! 2 the field point lies in or on the rectangle.
        integer, intent(in) :: n
        type(tm_rect_coil_t), intent(in) :: coil
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
        previous = tm_rect_coil_flux(n, coil, r, z, orders(1), orders(1))
        order = orders(1)
        psi = previous
        do k = 2, size(orders)
            psi = tm_rect_coil_flux(n, coil, r, z, orders(k), orders(k))
            order = orders(k)
            if (abs(psi - previous) <= tol*abs(psi)) then
                status = 0
                return
            end if
            previous = psi
        end do
    end subroutine tm_rect_coil_flux_converged

    subroutine tm_rect_coil_flux_table(n, coils, r, z, nr, nz, psi)
        !! psi(i, k) = tm_rect_coil_flux(n, coils(k), r(i), z(i), nr, nz) for
        !! all points and coils, evaluated in parallel over the entries.  Each
        !! entry is computed independently, so the table does not depend on
        !! the number of threads.
        integer, intent(in) :: n
        type(tm_rect_coil_t), intent(in) :: coils(:)
        real(dp), intent(in) :: r(:), z(:)
        integer, intent(in) :: nr, nz
        real(dp), intent(out) :: psi(size(r), size(coils))
        integer :: i, k

        !$omp parallel do collapse(2) schedule(static) default(shared) &
        !$omp& private(i, k)
        do k = 1, size(coils)
            do i = 1, size(r)
                psi(i, k) = tm_rect_coil_flux(n, coils(k), r(i), z(i), nr, nz)
            end do
        end do
        !$omp end parallel do
    end subroutine tm_rect_coil_flux_table

    ! ------------------------------------------------------------------
    ! Private machinery
    ! ------------------------------------------------------------------

    pure function weight(r, rp) result(w)
        !! 1/(2 pi sqrt(R R')).
        real(dp), intent(in) :: r, rp
        real(dp) :: w

        w = 0.5_dp/(pi*sqrt(r*rp))
    end function weight

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

    pure subroutine q_modes(nmax, r, rp, zeta, q, qx)
        !! Q_{m-1/2}(chi) and dQ/dchi for m = 0..nmax.
        integer, intent(in) :: nmax
        real(dp), intent(in) :: r, rp, zeta
        real(dp), intent(out) :: q(0:nmax), qx(0:nmax)
        real(dp) :: t, chi, p, pt, rr, rt
        integer :: m

        t = ((r - rp)**2 + zeta**2)/(4.0_dp*r*rp)
        if (.not. (t > 0.0_dp)) then
            q = ieee_value(1.0_dp, ieee_positive_inf)
            qx = -ieee_value(1.0_dp, ieee_positive_inf)
            return
        end if
        if (t < tm_log_split_limit(nmax)) then
            do m = 0, nmax
                call log_split_series(m, t, p, pt, rr, rt)
                q(m) = -0.5_dp*log(t)*p + rr
                ! dQ/dchi = (1/2) dQ/dt, dQ/dt = -p/(2 t) - ln(t) pt/2 + rt
                qx(m) = 0.5_dp*(-0.5_dp*p/t - 0.5_dp*log(t)*pt + rt)
            end do
        else
            chi = (r*r + rp*rp + zeta*zeta)/(2.0_dp*r*rp)
            call toroidal_q_sequence(nmax, chi, q, qx)
        end if
    end subroutine q_modes

    pure subroutine mode_jet(m, order, r, rp, zeta, h, h1, h2)
        !! h = G_m and its derivatives (up to `order`) with respect to
        !! t = (R, R', zeta).  Entries beyond `order` are zero.
        integer, intent(in) :: m, order
        real(dp), intent(in) :: r, rp, zeta
        real(dp), intent(out) :: h, h1(3), h2(3, 3)
        real(dp) :: q(0:m), qx(0:m), w, wa(3), wab(3, 3)
        real(dp) :: chi, ca(3), cab(3, 3), qxx, lam, t
        integer :: a, b

        h1 = 0.0_dp
        h2 = 0.0_dp
        if (.not. (r > 0.0_dp .and. rp > 0.0_dp)) then
            h = ieee_value(1.0_dp, ieee_quiet_nan)
            return
        end if
        call q_modes(m, r, rp, zeta, q, qx)
        w = weight(r, rp)
        h = w*q(m)
        if (order < 1) return

        chi = (r*r + rp*rp + zeta*zeta)/(2.0_dp*r*rp)
        ca = [((r - rp)*(r + rp) - zeta**2)/(2.0_dp*r*r*rp), &
            ((rp - r)*(rp + r) - zeta**2)/(2.0_dp*r*rp*rp), zeta/(r*rp)]
        wa = [-0.5_dp*w/r, -0.5_dp*w/rp, 0.0_dp]
        h1 = wa*q(m) + w*qx(m)*ca
        if (order < 2) return

        t = ((r - rp)**2 + zeta**2)/(4.0_dp*r*rp)
        lam = real(m, dp)**2 - 0.25_dp
        qxx = (lam*q(m) - 2.0_dp*chi*qx(m))/(2.0_dp*t*(chi + 1.0_dp))
        cab = 0.0_dp
        cab(1, 1) = 2.0_dp*chi/(r*r) - 1.0_dp/(r*rp)
        cab(2, 2) = 2.0_dp*chi/(rp*rp) - 1.0_dp/(r*rp)
        cab(1, 2) = chi/(r*rp) - 1.0_dp/(rp*rp) - 1.0_dp/(r*r)
        cab(2, 1) = cab(1, 2)
        cab(1, 3) = -zeta/(r*r*rp)
        cab(3, 1) = cab(1, 3)
        cab(2, 3) = -zeta/(r*rp*rp)
        cab(3, 2) = cab(2, 3)
        cab(3, 3) = 1.0_dp/(r*rp)
        wab = 0.0_dp
        wab(1, 1) = 0.75_dp*w/(r*r)
        wab(2, 2) = 0.75_dp*w/(rp*rp)
        wab(1, 2) = 0.25_dp*w/(r*rp)
        wab(2, 1) = wab(1, 2)
        do b = 1, 3
            do a = 1, 3
                h2(a, b) = wab(a, b)*q(m) + wa(a)*qx(m)*ca(b) &
                    + wa(b)*qx(m)*ca(a) + w*qxx*ca(a)*ca(b) &
                    + w*qx(m)*cab(a, b)
            end do
        end do
    end subroutine mode_jet

    pure subroutine log_split_series(m, t, p, pt, r, rt)
        !! Taylor sums of p_m(t) = P_{m-1/2}(1 + 2 t) and the analytic part
        !! r_m(t) of Q_{m-1/2}(1 + 2 t) = -ln(t) p_m/2 + r_m, and their
        !! t-derivatives pt, rt; see the module header for the recurrences.
        integer, intent(in) :: m
        real(dp), intent(in) :: t
        real(dp), intent(out) :: p, pt, r, rt
        real(dp) :: lam, c, cn, e, en, tk, kk, a, d1, d2, ratio
        integer :: k, j

        lam = real(m, dp)**2 - 0.25_dp
        c = 1.0_dp
        e = 2.0_dp*ln2
        do j = 1, m
            e = e - 2.0_dp/real(2*j - 1, dp)
        end do
        tk = 1.0_dp
        p = 0.0_dp
        pt = 0.0_dp
        r = 0.0_dp
        rt = 0.0_dp
        if (.not. (t < split_t_max)) then
            p = ieee_value(1.0_dp, ieee_quiet_nan)
            pt = p
            r = p
            rt = p
            return
        end if
        do k = 0, 300
            kk = real(k, dp)
            a = lam - kk*(kk + 1.0_dp)
            cn = c*a/(kk + 1.0_dp)**2
            en = ((kk + 1.0_dp)*cn + (kk + 0.5_dp)*c + a*e)/(kk + 1.0_dp)**2
            p = p + c*tk
            r = r + e*tk
            d1 = (kk + 1.0_dp)*cn*tk
            d2 = (kk + 1.0_dp)*en*tk
            pt = pt + d1
            rt = rt + d2
            ratio = abs(a)*t/(kk + 1.0_dp)**2
            tk = tk*t
            c = cn
            e = en
            if (ratio < 0.5_dp .and. k >= 2) then
                if (abs(d1) <= 1.0e-18_dp*abs(pt) .and. &
                    abs(d2) <= 1.0e-18_dp*abs(rt) .and. &
                    abs(c*tk) <= 1.0e-18_dp*abs(p)) exit
            end if
        end do
    end subroutine log_split_series

end module fortfem_toroidal_mode_kernel
