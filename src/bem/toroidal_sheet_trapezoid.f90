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
    use fortfem_generated_toroidal_sheet_pair_value, only: &
        generated_toroidal_sheet_pair_value
    use fortfem_generated_ts_trap_pair_gradient_w, only: generated_ts_trap_pair_gradient_w
    use fortfem_generated_ts_trap_pair_hvp_w, only: generated_ts_trap_pair_hvp_w
    use fortfem_generated_ts_trap_zeta_term, only: generated_ts_trap_zeta_term
    use fortfem_generated_ts_trap_zeta_base, only: generated_ts_trap_zeta_base
    use fortfem_generated_ts_trap_diag, only: generated_ts_trap_diag
    implicit none
    private

    real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp

    type, public :: tst_grid_t
        integer :: nu = 0, nv = 0, nfp = 1
    end type

    public :: tst_energy, tst_hvp, tst_zeta, tst_current_form, tst_current_modes

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
        real(dp) :: gx(8), zero8(8)
        integer :: i1, j1, i2, j2, jj, nvt
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
        !$omp& private(i1, j1, i2, j2, jj, s, out, z, abc, d, gx) reduction(+:e)
        do j1 = 1, g%nv
            do i1 = 1, g%nu
                s = 0; gx = 0
                do j2 = 1, nvt
                    jj = 1 + mod(j2 - 1, g%nv)
                    do i2 = 1, g%nu
                        if (i2 == i1 .and. j2 == j1) cycle
                        if (present(gradient)) then
                            call generated_ts_trap_pair_gradient_w(u(1, i1, j1), &
                                u(2, i1, j1), u(3, i1, j1), u(4, i1, j1), u(5, i1, j1), &
                                u(6, i1, j1), u(7, i1, j1), u(8, i1, j1), u(1, i2, jj), &
                                u(2, i2, jj), u(3, i2, jj), u(4, i2, jj), u(5, i2, jj), &
                                u(6, i2, jj), u(7, i2, jj), u(8, i2, jj), 1.0_dp, &
                                cs(1, j1), cs(2, j1), cs(1, j2), cs(2, j2), &
                                cs(1, j2) - cs(1, j1), cs(2, j2) - cs(2, j1), &
                                inverse_distance(u(:, i1, j1), u(:, i2, jj), cs(:, j1), &
                                cs(:, j2)), out(1), out(2), out(3), out(4), out(5), out(6), &
                                out(7), out(8), out(9))
                            gx = gx + out(2:9)
                        else
                            call generated_toroidal_sheet_pair_value(u(1, i1, j1), &
                                u(2, i1, j1), u(3, i1, j1), u(4, i1, j1), u(5, i1, j1), &
                                u(6, i1, j1), u(7, i1, j1), u(8, i1, j1), u(1, i2, jj), &
                                u(2, i2, jj), u(3, i2, jj), u(4, i2, jj), u(5, i2, jj), &
                                u(6, i2, jj), u(7, i2, jj), u(8, i2, jj), 1.0_dp, &
                                cs(1, j1), cs(2, j1), cs(1, j2), cs(2, j2), &
                                cs(1, j2) - cs(1, j1), cs(2, j2) - cs(2, j1), out(1))
                        end if
                        s = s + out(1)
                    end do
                end do
                abc = metric(g, u(:, i1, j1))
                call tst_zeta(abc(1), abc(2), abc(3), z)
                call diag(u(:, i1, j1), hu, hv, z, zero8, d)
                e = e + s - d(1)
                if (present(gradient)) gradient(:, i1, j1) = &
                    g%nfp*w*w/(8*pi)*(2*gx - d(2:9))
            end do
        end do
        !$omp end parallel do
        energy = g%nfp*w*w/(8*pi)*e
    end subroutine tst_energy

    subroutine tst_hvp(g, u, du, hv_out)
        !! Hessian of E in u applied to du.
        type(tst_grid_t), intent(in) :: g
        real(dp), intent(in) :: u(:, :, :), du(:, :, :)
        real(dp), intent(out) :: hv_out(:, :, :)
        real(dp) :: hu, hv, w, cs(2, g%nv*g%nfp), out(8), z(10), abc(3), d(17), acc(8)
        integer :: i1, j1, i2, j2, jj, nvt
        hu = 1.0_dp/g%nu
        hv = 2*pi/(g%nfp*g%nv)
        w = hu*hv
        nvt = g%nv*g%nfp
        do j2 = 1, nvt
            cs(:, j2) = [cos((j2 - 1)*hv), sin((j2 - 1)*hv)]
        end do
        !$omp parallel do collapse(2) schedule(dynamic, 4) default(shared) &
        !$omp& private(i1, j1, i2, j2, jj, out, z, abc, d, acc)
        do j1 = 1, g%nv
            do i1 = 1, g%nu
                acc = 0
                do j2 = 1, nvt
                    jj = 1 + mod(j2 - 1, g%nv)
                    do i2 = 1, g%nu
                        if (i2 == i1 .and. j2 == j1) cycle
                        call generated_ts_trap_pair_hvp_w(u(1, i1, j1), u(2, i1, j1), &
                            u(3, i1, j1), u(4, i1, j1), u(5, i1, j1), u(6, i1, j1), &
                            u(7, i1, j1), u(8, i1, j1), u(1, i2, jj), u(2, i2, jj), &
                            u(3, i2, jj), u(4, i2, jj), u(5, i2, jj), u(6, i2, jj), &
                            u(7, i2, jj), u(8, i2, jj), 1.0_dp, cs(1, j1), cs(2, j1), &
                            cs(1, j2), cs(2, j2), cs(1, j2) - cs(1, j1), &
                            cs(2, j2) - cs(2, j1), inverse_distance(u(:, i1, j1), &
                            u(:, i2, jj), cs(:, j1), cs(:, j2)), du(1, i1, j1), du(2, i1, j1), &
                            du(3, i1, j1), du(4, i1, j1), du(5, i1, j1), du(6, i1, j1), &
                            du(7, i1, j1), du(8, i1, j1), du(1, i2, jj), du(2, i2, jj), &
                            du(3, i2, jj), du(4, i2, jj), du(5, i2, jj), du(6, i2, jj), &
                            du(7, i2, jj), du(8, i2, jj), out(1), out(2), out(3), &
                            out(4), out(5), out(6), out(7), out(8))
                        acc = acc + out
                    end do
                end do
                abc = metric(g, u(:, i1, j1))
                call tst_zeta(abc(1), abc(2), abc(3), z)
                call diag(u(:, i1, j1), hu, hv, z, du(:, i1, j1), d)
                hv_out(:, i1, j1) = g%nfp*w*w/(8*pi)*(2*acc - d(10:17))
            end do
        end do
        !$omp end parallel do
    end subroutine tst_hvp

    subroutine tst_current_form(g, u, form)
        !! Hessian of E in (Phi_u, Phi_v) of the first-period points,
        !! form(a, i, b, j) with a, b = 1 (Phi_u), 2 (Phi_v) and i, j the point
        !! index i1 + nu (j1 - 1); E is bilinear in them, so form does not
        !! depend on Phi.
        type(tst_grid_t), intent(in) :: g
        real(dp), intent(in) :: u(:, :, :)
        real(dp), intent(out) :: form(:, :, :, :)
        real(dp) :: hu, hv, w, cs(2, g%nv*g%nfp), z(10), abc(3), c
        real(dp) :: x1(3), xu1(3), xv1(3), x2(3), xu2(3), xv2(3), r
        integer :: i1, j1, i2, j2, jj, nvt, a, b
        hu = 1.0_dp/g%nu
        hv = 2*pi/(g%nfp*g%nv)
        w = hu*hv
        nvt = g%nv*g%nfp
        c = g%nfp*w*w/(8*pi)
        do j2 = 1, nvt
            cs(:, j2) = [cos((j2 - 1)*hv), sin((j2 - 1)*hv)]
        end do
        !$omp parallel do collapse(2) schedule(dynamic, 4) default(shared) &
        !$omp& private(i1, j1, i2, j2, jj, a, b, z, abc, x1, xu1, xv1, x2, xu2, xv2, r)
        do j1 = 1, g%nv
            do i1 = 1, g%nu
                a = i1 + g%nu*(j1 - 1)
                form(:, a, :, :) = 0
                call cartesian(u(:, i1, j1), cs(:, j1), x1, xu1, xv1)
                do j2 = 1, nvt
                    jj = 1 + mod(j2 - 1, g%nv)
                    do i2 = 1, g%nu
                        if (i2 == i1 .and. j2 == j1) cycle
                        b = i2 + g%nu*(jj - 1)
                        call cartesian(u(:, i2, jj), cs(:, j2), x2, xu2, xv2)
                        r = norm2(x2 - x1)
                        form(1, a, 1, b) = form(1, a, 1, b) + 2*c*dot_product(xv1, xv2)/r
                        form(1, a, 2, b) = form(1, a, 2, b) - 2*c*dot_product(xv1, xu2)/r
                        form(2, a, 1, b) = form(2, a, 1, b) - 2*c*dot_product(xu1, xv2)/r
                        form(2, a, 2, b) = form(2, a, 2, b) + 2*c*dot_product(xu1, xu2)/r
                    end do
                end do
                abc = metric(g, u(:, i1, j1))
                call tst_zeta(abc(1), abc(2), abc(3), z)
                form(1, a, 1, a) = form(1, a, 1, a) - 2*c*z(1)*dot_product(xv1, xv1)
                form(1, a, 2, a) = form(1, a, 2, a) + 2*c*z(1)*dot_product(xv1, xu1)
                form(2, a, 1, a) = form(2, a, 1, a) + 2*c*z(1)*dot_product(xu1, xv1)
                form(2, a, 2, a) = form(2, a, 2, a) - 2*c*z(1)*dot_product(xu1, xu1)
            end do
        end do
        !$omp end parallel do
    end subroutine tst_current_form

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
        real(dp) :: hu, hv, w, cs(2, g%nv*g%nfp), z(10), abc(3), c, blk(2, 2)
        real(dp) :: x1(3), xu1(3), xv1(3), x2(3), xu2(3), xv2(3), r
        real(dp), allocatable :: h1(:, :, :)
        integer :: i1, j1, i2, j2, jj, nvt, nm, nu, ia, n, ta
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
        !$omp& private(i1, j1, i2, j2, jj, z, abc, x1, xu1, xv1, x2, xu2, xv2, r, gg, blk, n, ta)
        do j1 = 1, g%nv
            do i1 = 1, nu
                if (.not. allocated(gg)) allocate (gg(2, 2, nu, nm))
                ! gg(ta, tb, i2, n') = sum over source columns of F(a, b) T_n'(b)
                gg = 0
                call cartesian(u(:, i1, j1), cs(:, j1), x1, xu1, xv1)
                do j2 = 1, nvt
                    jj = 1 + mod(j2 - 1, g%nv)
                    do i2 = 1, nu
                        if (i2 == i1 .and. j2 == j1) cycle
                        call cartesian(u(:, i2, jj), cs(:, j2), x2, xu2, xv2)
                        r = norm2(x2 - x1)
                        blk(1, 1) = dot_product(xv1, xv2)/r
                        blk(1, 2) = -dot_product(xv1, xu2)/r
                        blk(2, 1) = -dot_product(xu1, xv2)/r
                        blk(2, 2) = dot_product(xu1, xu2)/r
                        do n = 1, nm
                            gg(:, 1, i2, n) = gg(:, 1, i2, n) + blk(:, 1)*t(n, jj)
                            gg(:, 2, i2, n) = gg(:, 2, i2, n) + blk(:, 2)*dt(n, jj)
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

    pure subroutine cartesian(x, cs, p, xu, xv)
        real(dp), intent(in) :: x(8), cs(2)
        real(dp), intent(out) :: p(3), xu(3), xv(3)
        p = [x(1)*cs(1), x(1)*cs(2), x(2)]
        xu = [x(3)*cs(1), x(3)*cs(2), x(4)]
        xv = [x(5)*cs(1) - x(1)*cs(2), x(5)*cs(2) + x(1)*cs(1), x(6)]
    end subroutine cartesian

    pure subroutine diag(x, hu, hv, z, dx, d)
        !! |k|^2 Z of one point: value, gradient in x and its derivative along dx.
        real(dp), intent(in) :: x(8), hu, hv, z(10), dx(8)
        real(dp), intent(out) :: d(17)
        call generated_ts_trap_diag(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), hu, &
            hv, z(1), z(2), z(3), z(4), z(5), z(6), z(7), z(8), z(9), z(10), dx(1), &
            dx(2), dx(3), dx(4), dx(5), dx(6), dx(7), dx(8), d(1), d(2), d(3), d(4), &
            d(5), d(6), d(7), d(8), d(9), d(10), d(11), d(12), d(13), d(14), d(15), &
            d(16), d(17))
    end subroutine diag

end module fortfem_toroidal_sheet_trapezoid
