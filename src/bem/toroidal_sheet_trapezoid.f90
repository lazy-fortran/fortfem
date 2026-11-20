module fortfem_toroidal_sheet_trapezoid
    !! Magnetic energy of a surface current on a smooth toroidal surface
    !! sampled on a uniform grid, by the corrected trapezoidal rule:
    !!
    !!   E = nfp w^2/(8 pi) sum_{i in period 1} [ sum_{j /= i} k_i . k_j / |X_i - X_j|
    !!                                            - Z_i |k_i|^2 ],
    !!
    !! k = Phi_u X_v - Phi_v X_u, w = h_u h_v, j over the whole torus. Z_i is
    !! the Epstein zeta function at s = 1/2 of the lattice form Q(m, n) =
    !! |m h_u X_u + n h_v X_v|^2 at point i; subtracting it from the punctured
    !! sum gives the singular inner integral to O(h^3) (corrected trapezoidal
    !! rules: Marin, Runborg and Tornberg, Adv. Comput. Math. 40 (2014); Wu
    !! and Martinsson, SIAM J. Numer. Anal. 59 (2021)). For a smooth periodic
    !! surface the outer sum converges spectrally.
    !!
    !! Point inputs u(:, i, j) = (R, Z, R_u, Z_u, R_v, Z_v, Phi_u, Phi_v) at
    !! u = (i - 1) h_u, v = (j - 1) h_v of the first period (h_u = 1/nu, h_v =
    !! 2 pi/(nfp nv)); the other periods repeat them at v shifted by whole
    !! periods. The double sum is invariant under exchanging the two points
    !! with a shift by whole periods, so the gradient and the Hessian-vector
    !! product in the inputs of the first period take the target side twice.
    use fortfem_kinds, only: dp
    use fortfem_generated_ts_trap_pair_gradient_w, only: generated_ts_trap_pair_gradient_w
    use fortfem_generated_ts_trap_pair_hvp_w, only: generated_ts_trap_pair_hvp_w
    use fortfem_generated_ts_trap_pair_blocks_w, only: generated_ts_trap_pair_blocks_w
    use fortfem_generated_ts_trap_pair_gradient_pou, only: generated_ts_trap_pair_gradient_pou
    use fortfem_generated_ts_trap_pair_hvp_pou, only: generated_ts_trap_pair_hvp_pou
    use fortfem_generated_ts_trap_pair_blocks_pou, only: generated_ts_trap_pair_blocks_pou
    use fortfem_generated_ts_trap_zeta_term, only: generated_ts_trap_zeta_term
    use fortfem_generated_ts_trap_zeta_base, only: generated_ts_trap_zeta_base
    use fortfem_generated_ts_trap_diag, only: generated_ts_trap_diag
    implicit none
    private

    real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp

    !> r0 > 0: two-level partition of unity in the distance r. Pairs with
    !> r < r0 are summed on the grid with the weight psi(r/r0) (1 below r0/2,
    !> a C-infinity step to 0 at r0), the rest with 1 - psi on the grid of
    !> every second point in u and v (nu, nv even; weight 16 per pair): both
    !> integrands are smooth, so both sums stay spectral, and the far sum has
    !> 16 times fewer pairs. r0 = 0: every pair on the grid.
    !> jt /= 0: the energy minus that of the current Phi = jt v alone (the
    !> poloidal current of a toroidal-field coil, whose large self energy a
    !> caller adds exactly by other means); the subtraction is in the kernels.
    type, public :: tst_grid_t
        integer :: nu = 0, nv = 0, nfp = 1
        real(dp) :: r0 = 0, jt = 0
    end type

    public :: tst_energy, tst_hvp, tst_zeta, tst_current_modes
    public :: tst_hessian

contains

    subroutine tst_zeta(a, b, c, z)
        !! Epstein zeta at s = 1/2 of the form (a, b; b, c) with its gradient
        !! and Hessian (aa, ab, ac, bb, bc, cc) in (a, b, c): the Ewald split
        !! of the form scaled to unit determinant, cut where pi q > 36.
        real(dp), intent(in) :: a, b, c
        real(dp), intent(out) :: z(10)
        real(dp) :: t(10), l, tr, lmin
        integer :: m, n, mm
        l = sqrt(a*c - b*b)
        tr = (a + c)/l
        lmin = 0.5_dp*(tr - sqrt(max(tr*tr - 4, 0.0_dp)))
        mm = ceiling(sqrt(36/(pi*lmin))) + 1
        call generated_ts_trap_zeta_base(a, b, c, z(1), z(2), z(3), z(4), z(5), z(6), &
            z(7), z(8), z(9), z(10))
        do n = -mm, mm
            do m = -mm, mm
                if (m == 0 .and. n == 0) cycle
                call generated_ts_trap_zeta_term(a, b, c, real(m, dp), real(n, dp), &
                    t(1), t(2), t(3), t(4), t(5), t(6), t(7), t(8), t(9), t(10))
                z = z + t
            end do
        end do
    end subroutine tst_zeta

    pure real(dp) function inverse_distance(x, y, c1, c2) result(w)
        !! 1/|X_y - X_x| from the point inputs and (cos, sin) of both angles.
        real(dp), intent(in) :: x(8), y(8), c1(2), c2(2)
        w = 1/sqrt((y(1)*c2(1) - x(1)*c1(1))**2 + (y(1)*c2(2) - x(1)*c1(2))**2 &
            + (y(2) - x(2))**2)
    end function inverse_distance

    pure function metric(g, x) result(abc)
        type(tst_grid_t), intent(in) :: g
        real(dp), intent(in) :: x(8)
        real(dp) :: abc(3), hu, hv
        hu = 1.0_dp/g%nu
        hv = 2*pi/(g%nfp*g%nv)
        abc = [hu*hu*(x(3)**2 + x(4)**2), hu*hv*(x(3)*x(5) + x(4)*x(6)), &
            hv*hv*(x(5)**2 + x(1)**2 + x(6)**2)]
    end function metric

    subroutine tst_energy(g, u, energy, gradient)
        !! E and, when present, its gradient in u (same shape).
        type(tst_grid_t), intent(in) :: g
        real(dp), intent(in) :: u(:, :, :)
        real(dp), intent(out) :: energy
        real(dp), intent(out), optional :: gradient(:, :, :)
        real(dp) :: hu, hv, w, s, e, cs(2, g%nv*g%nfp), out(9), z(10), abc(3), d(17)
        real(dp) :: gx(8), zero8(8), wi, alpha, beta, scale
        integer :: i1, j1, i2, j2, jj, nvt, level, step, mode
        hu = 1.0_dp/g%nu
        hv = 2*pi/(g%nfp*g%nv)
        w = hu*hv
        nvt = g%nv*g%nfp
        do j2 = 1, nvt
            cs(:, j2) = [cos((j2 - 1)*hv), sin((j2 - 1)*hv)]
        end do
        zero8 = 0
        e = 0
        !$omp parallel do collapse(2) schedule(dynamic, 4) default(shared) &
        !$omp& private(i1, j1, i2, j2, jj, s, out, z, abc, d, gx, wi, alpha, beta, scale, &
        !$omp& level, step, mode) reduction(+:e)
        do j1 = 1, g%nv
            do i1 = 1, g%nu
                s = 0; gx = 0
                do level = 0, merge(1, 0, g%r0 > 0)
                    if (level == 1 .and. (mod(i1, 2) == 0 .or. mod(j1, 2) == 0)) exit
                    step = 1 + level
                    scale = merge(1.0_dp, 16.0_dp, level == 0)
                    do j2 = 1, nvt, step
                        jj = 1 + mod(j2 - 1, g%nv)
                        do i2 = 1, g%nu, step
                            if (i2 == i1 .and. j2 == j1) cycle
                            wi = inverse_distance(u(:, i1, j1), u(:, i2, jj), cs(:, j1), &
                                cs(:, j2))
                            call pair_class(g%r0, wi, level, mode, alpha, beta)
                            if (mode == 0) cycle
                            call pair_gradient(u(:, i1, j1), u(:, i2, jj), cs(:, j1), &
                                cs(:, j2), wi, mode, alpha, beta, g%r0, g%jt, out)
                            s = s + scale*out(1)
                            gx = gx + scale*out(2:9)
                        end do
                    end do
                end do
                abc = metric(g, u(:, i1, j1))
                call tst_zeta(abc(1), abc(2), abc(3), z)
                call diag(u(:, i1, j1), hu, hv, g%jt, z, zero8, d)
                e = e + s - d(1)
                if (present(gradient)) gradient(:, i1, j1) = &
                    g%nfp*w*w/(8*pi)*(2*gx - d(2:9))
            end do
        end do
        !$omp end parallel do
        energy = g%nfp*w*w/(8*pi)*e
    end subroutine tst_energy

    pure subroutine pair_class(r0, w, level, mode, alpha, beta)
        !! mode 0: no contribution, 1: the plain kernel, 2: times alpha + beta
        !! S(2 - 2 r/r0) (level 0: psi, level 1: 1 - psi).
        real(dp), intent(in) :: r0, w
        integer, intent(in) :: level
        integer, intent(out) :: mode
        real(dp), intent(out) :: alpha, beta
        real(dp) :: t
        alpha = 0; beta = 0
        if (r0 <= 0) then
            mode = 1
            return
        end if
        t = 1/(w*r0)
        if (level == 0) then
            mode = merge(0, merge(1, 2, t <= 0.5_dp), t >= 1)
            beta = 1
        else
            mode = merge(0, merge(1, 2, t >= 1), t <= 0.5_dp)
            alpha = 1; beta = -1
        end if
    end subroutine pair_class

    pure real(dp) function smooth_weight(r0, w, level) result(sig)
        !! The weight of pair_class as a number (for the bilinear current block).
        real(dp), intent(in) :: r0, w
        integer, intent(in) :: level
        real(dp) :: y, e1, e2
        integer :: mode
        real(dp) :: alpha, beta
        call pair_class(r0, w, level, mode, alpha, beta)
        sig = 0
        if (mode == 1) sig = 1
        if (mode /= 2) return
        y = 2 - 2/(w*r0)
        e1 = exp(-1/y); e2 = exp(-1/(1 - y))
        sig = alpha + beta*e1/(e1 + e2)
    end function smooth_weight

    subroutine pair_gradient(x, y, c1, c2, wi, mode, alpha, beta, r0, jt, o)
        real(dp), intent(in) :: x(8), y(8), c1(2), c2(2), wi, alpha, beta, r0, jt
        integer, intent(in) :: mode
        real(dp), intent(out) :: o(9)
        if (mode == 1) then
            call generated_ts_trap_pair_gradient_w(x(1), x(2), x(3), x(4), x(5), x(6), &
                x(7), x(8), y(1), y(2), y(3), y(4), y(5), y(6), y(7), y(8), 1.0_dp, &
                c1(1), c1(2), c2(1), c2(2), c2(1) - c1(1), c2(2) - c1(2), jt, wi, o(1), o(2), &
                o(3), o(4), o(5), o(6), o(7), o(8), o(9))
        else
            call generated_ts_trap_pair_gradient_pou(x(1), x(2), x(3), x(4), x(5), x(6), &
                x(7), x(8), y(1), y(2), y(3), y(4), y(5), y(6), y(7), y(8), 1.0_dp, &
                c1(1), c1(2), c2(1), c2(2), c2(1) - c1(1), c2(2) - c1(2), jt, wi, alpha, beta, &
                r0, o(1), o(2), o(3), o(4), o(5), o(6), o(7), o(8), o(9))
        end if
    end subroutine pair_gradient

    subroutine pair_hvp(x, y, c1, c2, wi, mode, alpha, beta, r0, jt, dx, dy, o)
        real(dp), intent(in) :: x(8), y(8), c1(2), c2(2), wi, alpha, beta, r0, jt, dx(8), &
            dy(8)
        integer, intent(in) :: mode
        real(dp), intent(out) :: o(8)
        if (mode == 1) then
            call generated_ts_trap_pair_hvp_w(x(1), x(2), x(3), x(4), x(5), x(6), x(7), &
                x(8), y(1), y(2), y(3), y(4), y(5), y(6), y(7), y(8), 1.0_dp, c1(1), c1(2), &
                c2(1), c2(2), c2(1) - c1(1), c2(2) - c1(2), jt, wi, dx(1), dx(2), dx(3), dx(4), &
                dx(5), dx(6), dx(7), dx(8), dy(1), dy(2), dy(3), dy(4), dy(5), dy(6), &
                dy(7), dy(8), o(1), o(2), o(3), o(4), o(5), o(6), o(7), o(8))
        else
            call generated_ts_trap_pair_hvp_pou(x(1), x(2), x(3), x(4), x(5), x(6), x(7), &
                x(8), y(1), y(2), y(3), y(4), y(5), y(6), y(7), y(8), 1.0_dp, c1(1), c1(2), &
                c2(1), c2(2), c2(1) - c1(1), c2(2) - c1(2), jt, wi, alpha, beta, r0, dx(1), &
                dx(2), dx(3), dx(4), dx(5), dx(6), dx(7), dx(8), dy(1), dy(2), dy(3), &
                dy(4), dy(5), dy(6), dy(7), dy(8), o(1), o(2), o(3), o(4), o(5), o(6), &
                o(7), o(8))
        end if
    end subroutine pair_hvp

    subroutine tst_hvp(g, u, du, hv_out)
        !! Hessian of E in u applied to du.
        type(tst_grid_t), intent(in) :: g
        real(dp), intent(in) :: u(:, :, :), du(:, :, :)
        real(dp), intent(out) :: hv_out(:, :, :)
        real(dp) :: hu, hv, w, cs(2, g%nv*g%nfp), out(8), z(10), abc(3), d(17), acc(8)
        real(dp) :: wi, alpha, beta, scale
        integer :: i1, j1, i2, j2, jj, nvt, level, step, mode
        hu = 1.0_dp/g%nu
        hv = 2*pi/(g%nfp*g%nv)
        w = hu*hv
        nvt = g%nv*g%nfp
        do j2 = 1, nvt
            cs(:, j2) = [cos((j2 - 1)*hv), sin((j2 - 1)*hv)]
        end do
        !$omp parallel do collapse(2) schedule(dynamic, 4) default(shared) &
        !$omp& private(i1, j1, i2, j2, jj, out, z, abc, d, acc, wi, alpha, beta, scale, &
        !$omp& level, step, mode)
        do j1 = 1, g%nv
            do i1 = 1, g%nu
                acc = 0
                do level = 0, merge(1, 0, g%r0 > 0)
                    if (level == 1 .and. (mod(i1, 2) == 0 .or. mod(j1, 2) == 0)) exit
                    step = 1 + level
                    scale = merge(1.0_dp, 16.0_dp, level == 0)
                    do j2 = 1, nvt, step
                        jj = 1 + mod(j2 - 1, g%nv)
                        do i2 = 1, g%nu, step
                            if (i2 == i1 .and. j2 == j1) cycle
                            wi = inverse_distance(u(:, i1, j1), u(:, i2, jj), cs(:, j1), &
                                cs(:, j2))
                            call pair_class(g%r0, wi, level, mode, alpha, beta)
                            if (mode == 0) cycle
                            call pair_hvp(u(:, i1, j1), u(:, i2, jj), cs(:, j1), cs(:, j2), &
                                wi, mode, alpha, beta, g%r0, g%jt, du(:, i1, j1), &
                                du(:, i2, jj), out)
                            acc = acc + scale*out
                        end do
                    end do
                end do
                abc = metric(g, u(:, i1, j1))
                call tst_zeta(abc(1), abc(2), abc(3), z)
                call diag(u(:, i1, j1), hu, hv, g%jt, z, du(:, i1, j1), d)
                hv_out(:, i1, j1) = g%nfp*w*w/(8*pi)*(2*acc - d(10:17))
            end do
        end do
        !$omp end parallel do
    end subroutine tst_hvp

    subroutine tst_current_modes(g, u, du, t, dt, mform)
        !! Hessian of E in the current potential phi(i, n) of the basis
        !! Phi = sum_n phi(i, n) T_n(v) at u_i (spectral in u: Phi_u = du
        !! (phi T)), T_n periodic over one field period, given by its values
        !! t(n, j) and v derivatives dt(n, j) at the grid of one period. mform
        !! (i + nu (n - 1), i' + nu (n' - 1)). The kernel is contracted over the
        !! source columns inside the pair loop, so the grid form is not stored.
        type(tst_grid_t), intent(in) :: g
        real(dp), intent(in) :: u(:, :, :), du(:, :), t(:, :), dt(:, :)
        real(dp), intent(out) :: mform(:, :)
        real(dp), allocatable :: fl(:, :, :, :, :), gg(:, :, :, :)
        real(dp) :: hu, hv, w, cs(2, g%nv*g%nfp), z(10), abc(3), c, blk(2, 2), sig
        real(dp) :: x1(3), xu1(3), xv1(3), x2(3), xu2(3), xv2(3), r
        real(dp), allocatable :: h1(:, :, :)
        integer :: i1, j1, i2, j2, jj, nvt, nm, nu, ia, n, ta, level
        hu = 1.0_dp/g%nu
        hv = 2*pi/(g%nfp*g%nv)
        w = hu*hv
        nvt = g%nv*g%nfp
        nu = g%nu
        nm = size(t, 1)
        c = g%nfp*w*w/(8*pi)
        do j2 = 1, nvt
            cs(:, j2) = [cos((j2 - 1)*hv), sin((j2 - 1)*hv)]
        end do
        ! fl(ta, i1, j1, :, :) = (F L)(target, type ta; basis (i', n'))
        allocate (fl(2, nu, g%nv, nu, nm))
        !$omp parallel do collapse(2) schedule(dynamic, 4) default(shared) &
        !$omp& private(i1, j1, i2, j2, jj, z, abc, x1, xu1, xv1, x2, xu2, xv2, r, gg, blk, n, ta, &
        !$omp& sig, level)
        do j1 = 1, g%nv
            do i1 = 1, nu
                if (.not. allocated(gg)) allocate (gg(2, 2, nu, nm))
                ! gg(ta, tb, i2, n') = sum over source columns of F(a, b) T_n'(b)
                gg = 0
                call cartesian(u(:, i1, j1), cs(:, j1), x1, xu1, xv1)
                do level = 0, merge(1, 0, g%r0 > 0)
                    if (level == 1 .and. (mod(i1, 2) == 0 .or. mod(j1, 2) == 0)) exit
                    do j2 = 1, nvt, 1 + level
                        jj = 1 + mod(j2 - 1, g%nv)
                        do i2 = 1, nu, 1 + level
                            if (i2 == i1 .and. j2 == j1) cycle
                            call cartesian(u(:, i2, jj), cs(:, j2), x2, xu2, xv2)
                            r = norm2(x2 - x1)
                            sig = smooth_weight(g%r0, 1/r, level)*merge(1, 16, level == 0)
                            if (sig == 0) cycle
                            blk(1, 1) = sig*dot_product(xv1, xv2)/r
                            blk(1, 2) = -sig*dot_product(xv1, xu2)/r
                            blk(2, 1) = -sig*dot_product(xu1, xv2)/r
                            blk(2, 2) = sig*dot_product(xu1, xu2)/r
                            do n = 1, nm
                                gg(:, 1, i2, n) = gg(:, 1, i2, n) + blk(:, 1)*t(n, jj)
                                gg(:, 2, i2, n) = gg(:, 2, i2, n) + blk(:, 2)*dt(n, jj)
                            end do
                        end do
                    end do
                end do
                gg = 2*c*gg
                abc = metric(g, u(:, i1, j1))
                call tst_zeta(abc(1), abc(2), abc(3), z)
                blk(1, 1) = -2*c*z(1)*dot_product(xv1, xv1)
                blk(1, 2) = 2*c*z(1)*dot_product(xv1, xu1)
                blk(2, 1) = blk(1, 2)
                blk(2, 2) = -2*c*z(1)*dot_product(xu1, xu1)
                do n = 1, nm
                    gg(:, 1, i1, n) = gg(:, 1, i1, n) + blk(:, 1)*t(n, j1)
                    gg(:, 2, i1, n) = gg(:, 2, i1, n) + blk(:, 2)*dt(n, j1)
                end do
                ! Phi_u of the source = du applied along i2 to the T-weighted sum
                do ta = 1, 2
                    fl(ta, i1, j1, :, :) = matmul(transpose(du), gg(ta, 1, :, :)) &
                        + gg(ta, 2, :, :)
                end do
            end do
        end do
        !$omp end parallel do
        ! M = L^T (F L): the target side the same way
        allocate (h1(nu, nm, nu*nm))
        h1 = 0
        do j1 = 1, g%nv
            do n = 1, nm
                do i1 = 1, nu
                    h1(i1, n, :) = h1(i1, n, :) + t(n, j1)*reshape(fl(1, i1, j1, :, :), [nu*nm])
                end do
            end do
        end do
        mform = 0
        do n = 1, nm
            do ia = 1, nu
                mform(ia + nu*(n - 1), :) = matmul(du(:, ia), h1(:, n, :))
            end do
        end do
        do j1 = 1, g%nv
            do n = 1, nm
                do i1 = 1, nu
                    mform(i1 + nu*(n - 1), :) = mform(i1 + nu*(n - 1), :) &
                        + dt(n, j1)*reshape(fl(2, i1, j1, :, :), [nu*nm])
                end do
            end do
        end do
    end subroutine tst_current_modes

    subroutine tst_hessian(g, u, h)
        !! Dense Hessian of E in the inputs of the first period, h(c, a, c', b)
        !! with a, b = i + nu (j - 1); for coarse grids (8 nu nv squared
        !! entries).
        type(tst_grid_t), intent(in) :: g
        real(dp), intent(in) :: u(:, :, :)
        real(dp), intent(out) :: h(:, :, :, :)
        real(dp) :: hu, hv, w, cs(2, g%nv*g%nfp), out(100), z(10), abc(3), d(17), c
        real(dp) :: e8(8), hxx(8, 8), hxy(8, 8), wi, alpha, beta
        integer :: i1, j1, i2, j2, jj, nvt, a, b, k, r, q, level, mode
        hu = 1.0_dp/g%nu
        hv = 2*pi/(g%nfp*g%nv)
        w = hu*hv
        nvt = g%nv*g%nfp
        c = g%nfp*w*w/(8*pi)
        do j2 = 1, nvt
            cs(:, j2) = [cos((j2 - 1)*hv), sin((j2 - 1)*hv)]
        end do
        !$omp parallel do collapse(2) schedule(dynamic, 4) default(shared) &
        !$omp& private(i1, j1, i2, j2, jj, a, b, k, r, q, out, z, abc, d, e8, hxx, hxy, wi, &
        !$omp& alpha, beta, level, mode)
        do j1 = 1, g%nv
            do i1 = 1, g%nu
                a = i1 + g%nu*(j1 - 1)
                h(:, a, :, :) = 0
                do level = 0, merge(1, 0, g%r0 > 0)
                    if (level == 1 .and. (mod(i1, 2) == 0 .or. mod(j1, 2) == 0)) exit
                    do j2 = 1, nvt, 1 + level
                        jj = 1 + mod(j2 - 1, g%nv)
                        do i2 = 1, g%nu, 1 + level
                            if (i2 == i1 .and. j2 == j1) cycle
                            b = i2 + g%nu*(jj - 1)
                            wi = inverse_distance(u(:, i1, j1), u(:, i2, jj), cs(:, j1), cs(:, j2))
                            call pair_class(g%r0, wi, level, mode, alpha, beta)
                            if (mode == 0) cycle
                            call blocks(u(:, i1, j1), u(:, i2, jj), cs(:, j1), cs(:, j2), wi, &
                                mode, alpha, beta, out)
                            out = merge(1, 16, level == 0)*out
                            k = 0
                            do r = 1, 8
                                do q = r, 8
                                    k = k + 1
                                    hxx(r, q) = out(k); hxx(q, r) = out(k)
                                end do
                            end do
                            do r = 1, 8
                                do q = 1, 8
                                    k = k + 1
                                    hxy(r, q) = out(k)
                                end do
                            end do
                            h(:, a, :, a) = h(:, a, :, a) + 2*c*hxx
                            h(:, a, :, b) = h(:, a, :, b) + 2*c*hxy
                        end do
                    end do
                end do
                abc = metric(g, u(:, i1, j1))
                call tst_zeta(abc(1), abc(2), abc(3), z)
                do q = 1, 8
                    e8 = 0; e8(q) = 1
                    call diag(u(:, i1, j1), hu, hv, g%jt, z, e8, d)
                    h(:, a, q, a) = h(:, a, q, a) - c*d(10:17)
                end do
            end do
        end do
        !$omp end parallel do
    contains
        subroutine blocks(x, y, c1, c2, wi, mode, alpha, beta, o)
            real(dp), intent(in) :: x(8), y(8), c1(2), c2(2), wi, alpha, beta
            integer, intent(in) :: mode
            real(dp), intent(out) :: o(100)
            if (mode == 2) then
                call generated_ts_trap_pair_blocks_pou(x(1), x(2), x(3), x(4), x(5), x(6), &
                    x(7), x(8), y(1), y(2), y(3), y(4), y(5), y(6), y(7), y(8), 1.0_dp, &
                    c1(1), c1(2), c2(1), c2(2), c2(1) - c1(1), c2(2) - c1(2), g%jt, wi, alpha, &
                    beta, g%r0, o(1), o(2), o(3), o(4), o(5), o(6), o(7), o(8), o(9), o(10), &
                    o(11), o(12), o(13), o(14), o(15), o(16), o(17), o(18), o(19), o(20), &
                    o(21), o(22), o(23), o(24), o(25), o(26), o(27), o(28), o(29), o(30), &
                    o(31), o(32), o(33), o(34), o(35), o(36), o(37), o(38), o(39), o(40), &
                    o(41), o(42), o(43), o(44), o(45), o(46), o(47), o(48), o(49), o(50), &
                    o(51), o(52), o(53), o(54), o(55), o(56), o(57), o(58), o(59), o(60), &
                    o(61), o(62), o(63), o(64), o(65), o(66), o(67), o(68), o(69), o(70), &
                    o(71), o(72), o(73), o(74), o(75), o(76), o(77), o(78), o(79), o(80), &
                    o(81), o(82), o(83), o(84), o(85), o(86), o(87), o(88), o(89), o(90), &
                    o(91), o(92), o(93), o(94), o(95), o(96), o(97), o(98), o(99), o(100))
                return
            end if
            call generated_ts_trap_pair_blocks_w(x(1), x(2), x(3), x(4), x(5), x(6), &
                x(7), x(8), y(1), y(2), y(3), y(4), y(5), y(6), y(7), y(8), 1.0_dp, &
                c1(1), c1(2), c2(1), c2(2), c2(1) - c1(1), c2(2) - c1(2), g%jt, &
                inverse_distance(x, y, c1, c2), o(1), o(2), o(3), o(4), o(5), o(6), &
                o(7), o(8), o(9), o(10), o(11), o(12), o(13), o(14), o(15), o(16), &
                o(17), o(18), o(19), o(20), o(21), o(22), o(23), o(24), o(25), o(26), &
                o(27), o(28), o(29), o(30), o(31), o(32), o(33), o(34), o(35), o(36), &
                o(37), o(38), o(39), o(40), o(41), o(42), o(43), o(44), o(45), o(46), &
                o(47), o(48), o(49), o(50), o(51), o(52), o(53), o(54), o(55), o(56), &
                o(57), o(58), o(59), o(60), o(61), o(62), o(63), o(64), o(65), o(66), &
                o(67), o(68), o(69), o(70), o(71), o(72), o(73), o(74), o(75), o(76), &
                o(77), o(78), o(79), o(80), o(81), o(82), o(83), o(84), o(85), o(86), &
                o(87), o(88), o(89), o(90), o(91), o(92), o(93), o(94), o(95), o(96), &
                o(97), o(98), o(99), o(100))
        end subroutine blocks
    end subroutine tst_hessian

    pure subroutine cartesian(x, cs, p, xu, xv)
        real(dp), intent(in) :: x(8), cs(2)
        real(dp), intent(out) :: p(3), xu(3), xv(3)
        p = [x(1)*cs(1), x(1)*cs(2), x(2)]
        xu = [x(3)*cs(1), x(3)*cs(2), x(4)]
        xv = [x(5)*cs(1) - x(1)*cs(2), x(5)*cs(2) + x(1)*cs(1), x(6)]
    end subroutine cartesian

    pure subroutine diag(x, hu, hv, jt, z, dx, d)
        !! |k|^2 Z of one point: value, gradient in x and its derivative along dx.
        real(dp), intent(in) :: x(8), hu, hv, jt, z(10), dx(8)
        real(dp), intent(out) :: d(17)
        call generated_ts_trap_diag(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), hu, &
            hv, jt, z(1), z(2), z(3), z(4), z(5), z(6), z(7), z(8), z(9), z(10), dx(1), &
            dx(2), dx(3), dx(4), dx(5), dx(6), dx(7), dx(8), d(1), d(2), d(3), d(4), &
            d(5), d(6), d(7), d(8), d(9), d(10), d(11), d(12), d(13), d(14), d(15), &
            d(16), d(17))
    end subroutine diag

end module fortfem_toroidal_sheet_trapezoid
