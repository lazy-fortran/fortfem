module fortfem_toroidal_mode_boundary_operators_2d
    !! Galerkin boundary operators of the toroidal-mode Laplacian
    !!
    !!     L_n u = -div(R grad u) + (n^2/R) u,        (R, Z) in the half plane,
    !!
    !! (the e^{i n phi} component of the 3D Laplacian, weight R) on a closed
    !! meridional curve Gamma made of curved P2 panels, and the exterior
    !! Steklov-Poincare map S_n.  The curve type, the trace and flux spaces,
    !! the panel quadrature rules and the operator container are those of
    !! fortfem_axisymmetric_gs_boundary_operators_2d.
    !!
    !! Reference.  The weight R and the n^2/R term are those of the per-mode
    !! reduction of Albert, Lainer and Biro, Comput. Phys. Commun. 277 (2022)
    !! 108401, doi:10.1016/j.cpc.2022.108401 (Eqs. (34)-(36), (44), (47)); the
    !! Fourier convention f = sum_n f_n exp(i n phi) is its Eq. (31) and that
    !! of the MEPHIT paper, J. Phys.: Conf. Ser. 2397 (2022) 012014
    !! (doi:10.1088/1742-6596/2397/1/012014), Eq. (10).  The Poincare-Steklov
    !! operator S_n is the exterior closure that paper names as future work
    !! for the homogeneous Dirichlet box of its Eq. (11).
    !!
    !! Conventions.  Gamma is counter-clockwise, the unit normal n points out
    !! of the enclosed region (into the exterior), ds is arc length and the
    !! conormal derivative is R d_n.  With G_n of fortfem_toroidal_mode_kernel
    !! (L_n G_n = delta, G_n decays at infinity and is regular on the axis)
    !!
    !!   V t (x) = int G_n(x,y) t(y) ds_y
    !!   K u (x) = int R_y d_{n_y} G_n(x,y) u(y) ds_y
    !!   K' t(x) = R_x d_{n_x} int G_n(x,y) t(y) ds_y
    !!   W u (x) = -R_x d_{n_x} int R_y d_{n_y} G_n(x,y) u(y) ds_y.
    !!
    !! For an exterior field with L_n u = 0, trace u and t = R d_n u
    !! (n into the exterior)
    !!
    !!   V t = (K - 1/2) u,    W u = -(1/2 + K') t,    S_n u := -t = V^{-1}(1/2 - K) u,
    !!
    !! and <S_n u, u> is the exterior energy int (R |grad u|^2 + n^2 u^2/R).
    !! V is symmetric positive definite for every n.  K 1 = -1/2 holds for
    !! n = 0 only (n^2 u/R destroys harmonicity of the constant).  The
    !! Steklov form S_h = W + (M/2 - K)^T V^{-1} (M/2 - K) and the
    !! Johnson-Nedelec form are the ones of the GS module and are provided
    !! here under the tm_ prefix.
    !!
    !! Relation to the Grad-Shafranov operators.  The GS operators act on
    !! the poloidal flux psi = R A_phi and use the weight 1/R; their kernel is
    !! R R' G_1, so they are the n = 1 vector-potential operators rescaled by
    !! R, not the n = 0 scalar operators of this module.
    !!
    !! Hypersingular form (derivation).  For the surface of revolution S and
    !! the 3D Maue form <W3 u, v> = int_S int_S G3(x,y) curl_S u(y) . curl_S v(x)
    !! (G3 = 1/(4 pi |x-y|), curl_S = normal x grad_S) take u e^{i n phi} and
    !! v e^{-i n phi}.  With tau = (tau_R, tau_Z) the meridian tangent,
    !! n = (tau_Z, -tau_R) and (e_R, e_phi, e_Z) right handed,
    !!
    !!   curl_S (u e^{i n phi}) = e^{i n phi} [ (i n u/R) tau_R e_R - u' e_phi
    !!                                          + (i n u/R) tau_Z e_Z ],
    !!
    !! (' = d/ds).  The scalar products of the cylinder frames at phi and phi'
    !! are e_R.e_R' = e_phi.e_phi' = cos D, e_R.e_phi' = sin D,
    !! e_phi.e_R' = -sin D, e_Z.e_Z = 1 (D = phi - phi').  Integrating over
    !! phi' with G3 e^{-i n D} cos D -> C = (G_{n-1} + G_{n+1})/2,
    !! G3 e^{-i n D} sin D -> (G_{n-1} - G_{n+1})/(2 i), G3 e^{-i n D} -> G_n
    !! (G_m := G_|m|, each with the factor 2 pi of phi'), and dividing by the
    !! 2 pi of the phi integral, with dS = R ds dphi,
    !!
    !!   <W u, v> = int int { n^2 u(y) v(x) [ tau_R(y) tau_R(x) C
    !!                                        + tau_Z(y) tau_Z(x) G_n ]
    !!                      + R_x R_y u'(y) v'(x) C
    !!                      + n D [ R_x u(y) v'(x) tau_R(y)
    !!                              + R_y u'(y) v(x) tau_R(x) ] } ds ds,
    !!
    !! with D = (G_{n-1} - G_{n+1})/2.  All kernels are logarithmically
    !! singular (D is bounded), so no hypersingular integral remains and
    !! the form is symmetric.  For n = 0 it is int int R R' G_1 u' v' ds ds,
    !! the planar Maue form d_s V d_s with the kernel R R' G_1.
    !!
    !! Quadrature.  As in the GS module: composite Gauss rules for separated
    !! panels, Duffy rules with a geometrically graded radial rule on the same
    !! and on neighbouring panels, where the kernels are evaluated through
    !! the logarithmic split G_m = -ln(dist) coeff_m + smooth_m with the
    !! difference x - y formed without cancellation from the panel
    !! polynomials.
    !!
    !! Parallelism.  Panel pairs are coloured by parity (panels p and p + 2
    !! share no trace node; a third colour for an odd number of panels) and
    !! the pairs of one colour class run under OpenMP.  The accumulation
    !! order of each matrix entry is fixed, so the result is independent of
    !! the number of threads.
    use fortfem_kinds, only: dp
    use fortfem_toroidal_mode_kernel, only: tm_green_gradient, &
        tm_green_modes, tm_green_log_split_gradient, tm_green_modes_log_split, &
        tm_log_split_limit
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gs_bem_curve_t, &
        gs_bem_operators_t, gs_bem_panel_point, gs_bem_trace_node_index, &
        p2_basis, flux_basis, gauss_unit, graded_rule, build_self_rule, &
        build_corner_rule, gs_bem_steklov_symmetric, gs_bem_steklov_jn, &
        gs_bem_exterior_neumann, gs_bem_neumann_to_dirichlet
    implicit none
    private

    public :: tm_bem_assemble
    public :: tm_bem_steklov_symmetric
    public :: tm_bem_steklov_jn
    public :: tm_bem_exterior_neumann
    public :: tm_bem_neumann_to_dirichlet

contains

    ! ------------------------------------------------------------------
    ! Kernel values at one pair of points
    ! ------------------------------------------------------------------

    pure subroutine kernel_values(n, r, z, rp, zp, dr, dz, nx, ny, need_w, &
            g, kd, fd, wc, wd)
        !! G_n, K-kernel kd = R_y n_y . grad_y G_n, K'-kernel
        !! fd = R_x n_x . grad_x G_n and, with need_w, the W kernels
        !! wc = (G_{n-1} + G_{n+1})/2 and wd = (G_{n-1} - G_{n+1})/2.
        !! (dr, dz) is the accurately formed difference x - y.
        integer, intent(in) :: n
        real(dp), intent(in) :: r, z, rp, zp, dr, dz, nx(2), ny(2)
        logical, intent(in) :: need_w
        real(dp), intent(out) :: g, kd, fd, wc, wd
        integer :: na, nmax, lo, hi
        real(dp) :: d2, t, c, dc(4), sm, dsm(4), lnd, gx(2), gy(2), grad(4)
        real(dp) :: cm(0:abs(n) + 1), sm_m(0:abs(n) + 1), gm(0:abs(n) + 1)
        logical :: split

        na = abs(n)
        nmax = na + 1
        lo = abs(na - 1)
        hi = na + 1
        d2 = dr*dr + dz*dz
        t = d2/(4.0_dp*r*rp)
        split = t < tm_log_split_limit(nmax)
        wc = 0.0_dp
        wd = 0.0_dp
        if (split) then
            call tm_green_log_split_gradient(na, r, z, rp, zp, c, dc, sm, dsm)
            lnd = 0.5_dp*log(d2)
            g = -c*lnd + sm
            gx(1) = -c*dr/d2 - lnd*dc(1) + dsm(1)
            gx(2) = -c*dz/d2 - lnd*dc(2) + dsm(2)
            gy(1) = c*dr/d2 - lnd*dc(3) + dsm(3)
            gy(2) = c*dz/d2 - lnd*dc(4) + dsm(4)
            if (need_w) then
                call tm_green_modes_log_split(nmax, r, z, rp, zp, cm, sm_m)
                gm = -cm*lnd + sm_m
            end if
        else
            call tm_green_gradient(na, r, z, rp, zp, g, grad)
            gx = grad(1:2)
            gy = grad(3:4)
            if (need_w) call tm_green_modes(nmax, r, z, rp, zp, gm)
        end if
        fd = r*(nx(1)*gx(1) + nx(2)*gx(2))
        kd = rp*(ny(1)*gy(1) + ny(2)*gy(2))
        if (need_w) then
            wc = 0.5_dp*(gm(lo) + gm(hi))
            wd = 0.5_dp*(gm(lo) - gm(hi))
            ! For n = 0 the modes n - 1 and n + 1 coincide (D = 0).
            if (na == 0) wd = 0.0_dp
        end if
    end subroutine kernel_values

    ! ------------------------------------------------------------------
    ! Assembly
    ! ------------------------------------------------------------------

    subroutine tm_bem_assemble(curve, n, flux_degree, ops, status, hypersingular)
        !! Assemble V_n, K_n, K'_n, W_n and the mass matrices.  W (the most
        !! expensive kernel) is skipped when hypersingular = .false.
        type(gs_bem_curve_t), intent(in) :: curve
        integer, intent(in) :: n, flux_degree
        type(gs_bem_operators_t), intent(out) :: ops
        integer, intent(out) :: status
        logical, intent(in), optional :: hypersingular

        real(dp), allocatable :: rad_x(:), rad_w(:)
        real(dp), allocatable :: self_xi(:), self_eta(:), self_w(:), &
            self_oa(:)
        real(dp), allocatable :: nxt_xi(:), nxt_eta(:), nxt_w(:), &
            nxt_oa(:), nxt_ob(:)
        real(dp), allocatable :: prv_xi(:), prv_eta(:), prv_w(:), &
            prv_oa(:), prv_ob(:)
        real(dp) :: lens(curve%n_panels)
        real(dp) :: nn
        integer :: np, nx, ny, nl, a, b, kind, ncol, col_a, col_b, ia, ib
        integer, allocatable :: members(:, :), count_col(:)
        logical :: want_w

        status = 1
        want_w = .true.
        if (present(hypersingular)) want_w = hypersingular
        if (curve%n_panels < 4 .or. flux_degree < 0 .or. flux_degree > 4) return
        nn = real(n, dp)
        np = curve%n_panels
        nl = flux_degree + 1
        nx = 2*np
        ny = nl*np
        ops%nx = nx
        ops%ny = ny
        ops%flux_degree = flux_degree
        allocate(ops%v(ny, ny), ops%k(ny, nx), ops%kt(nx, ny), &
            ops%w(nx, nx), ops%m(ny, nx), ops%mxx(nx, nx), ops%kxx(nx, nx))
        ops%v = 0.0_dp
        ops%k = 0.0_dp
        ops%kt = 0.0_dp
        ops%w = 0.0_dp
        ops%m = 0.0_dp
        ops%mxx = 0.0_dp
        ops%kxx = 0.0_dp

        call graded_rule(rad_x, rad_w)
        call build_self_rule(rad_x, rad_w, self_xi, self_eta, self_w, self_oa)
        call build_corner_rule(rad_x, rad_w, .true., nxt_xi, nxt_eta, nxt_w, &
            nxt_oa, nxt_ob)
        call build_corner_rule(rad_x, rad_w, .false., prv_xi, prv_eta, prv_w, &
            prv_oa, prv_ob)

        call panel_lengths(lens)
        call assemble_mass()
        ! Panel-pair parity colouring: see the module header.
        ncol = 2
        if (mod(np, 2) == 1) ncol = 3
        allocate(members(np, ncol), count_col(ncol))
        count_col = 0
        do a = 1, np
            col_a = panel_colour(a)
            count_col(col_a) = count_col(col_a) + 1
            members(count_col(col_a), col_a) = a
        end do
        do col_a = 1, ncol
            do col_b = 1, ncol
                !$omp parallel do schedule(dynamic, 1) &
                !$omp& default(shared) private(ia, ib, a, b, kind)
                do ib = 1, count_col(col_b)
                    do ia = 1, count_col(col_a)
                        a = members(ia, col_a)
                        b = members(ib, col_b)
                        if (b == a) then
                            kind = 1
                        else if (b == mod(a, np) + 1) then
                            kind = 2
                        else if (a == mod(b, np) + 1) then
                            kind = 3
                        else
                            kind = 0
                        end if
                        call integrate_pair(a, b, kind)
                    end do
                end do
                !$omp end parallel do
            end do
        end do
        status = 0

    contains

        pure function panel_colour(p) result(c)
            integer, intent(in) :: p
            integer :: c

            c = 1 + mod(p - 1, 2)
            if (mod(np, 2) == 1 .and. p == np) c = 3
        end function panel_colour

        subroutine panel_lengths(len_p)
            real(dp), intent(out) :: len_p(:)
            real(dp) :: x0(2), x1(2), tau(2), sp
            integer :: p

            do p = 1, np
                call gs_bem_panel_point(curve, p, 0.0_dp, x0, tau, sp)
                call gs_bem_panel_point(curve, p, 1.0_dp, x1, tau, sp)
                len_p(p) = norm2(x1 - x0)
            end do
        end subroutine panel_lengths

        subroutine assemble_mass()
            real(dp) :: gx(10), gw(10), xi, x(2), tau(2), sp
            real(dp) :: psi(3), dpsi(3), phi(nl)
            integer :: p, q, i, j, ii, jj

            call gauss_unit(10, gx, gw)
            do p = 1, np
                do q = 1, 10
                    xi = gx(q)
                    call gs_bem_panel_point(curve, p, xi, x, tau, sp)
                    call p2_basis(xi, psi, dpsi)
                    call flux_basis(flux_degree, xi, phi)
                    do j = 1, 3
                        jj = gs_bem_trace_node_index(curve, p, j)
                        do i = 1, nl
                            ii = (p - 1)*nl + i
                            ops%m(ii, jj) = ops%m(ii, jj) &
                                + gw(q)*sp*phi(i)*psi(j)
                        end do
                        do i = 1, 3
                            ii = gs_bem_trace_node_index(curve, p, i)
                            ops%mxx(ii, jj) = ops%mxx(ii, jj) &
                                + gw(q)*sp*psi(i)*psi(j)
                        end do
                    end do
                end do
            end do
        end subroutine assemble_mass

        subroutine integrate_pair(pa, pb, pair_kind)
            integer, intent(in) :: pa, pb, pair_kind
            real(dp), allocatable :: fx(:), fw(:)
            real(dp) :: dmin, hmax
            integer :: ng, m, i, j, ia, ib, k

            select case (pair_kind)
            case (1)
                do k = 1, size(self_xi)
                    call add_point(pa, pb, self_xi(k), self_eta(k), &
                        self_w(k), 1, self_oa(k), 0.0_dp)
                end do
            case (2)
                do k = 1, size(nxt_xi)
                    call add_point(pa, pb, nxt_xi(k), nxt_eta(k), nxt_w(k), &
                        2, nxt_oa(k), nxt_ob(k))
                end do
            case (3)
                do k = 1, size(prv_xi)
                    call add_point(pa, pb, prv_xi(k), prv_eta(k), prv_w(k), &
                        3, prv_oa(k), prv_ob(k))
                end do
            case default
                hmax = max(lens(pa), lens(pb))
                dmin = pair_distance(pa, pb)
                m = max(1, min(12, ceiling(2.5_dp*hmax/max(dmin, 1.0e-14_dp))))
                ng = 8
                if (dmin > 8.0_dp*hmax) ng = 6
                allocate(fx(ng), fw(ng))
                call gauss_unit(ng, fx, fw)
                do ia = 1, m
                    do i = 1, ng
                        do ib = 1, m
                            do j = 1, ng
                                call add_point(pa, pb, &
                                    (real(ia - 1, dp) + fx(i))/m, &
                                    (real(ib - 1, dp) + fx(j))/m, &
                                    fw(i)*fw(j)/real(m*m, dp), 0, &
                                    0.0_dp, 0.0_dp)
                            end do
                        end do
                    end do
                end do
                deallocate(fx, fw)
            end select
        end subroutine integrate_pair

        function pair_distance(pa, pb) result(dmin)
            integer, intent(in) :: pa, pb
            real(dp) :: dmin, xa(2), xb(2), tau(2), sp
            integer :: i, j

            dmin = huge(1.0_dp)
            do i = 0, 4
                call gs_bem_panel_point(curve, pa, 0.25_dp*i, xa, tau, sp)
                do j = 0, 4
                    call gs_bem_panel_point(curve, pb, 0.25_dp*j, xb, tau, sp)
                    dmin = min(dmin, norm2(xa - xb))
                end do
            end do
        end function pair_distance

        subroutine add_point(pa, pb, xi, eta, wt, mode, oa, ob)
            !! Add the contribution of the point pair (xi on pa, eta on pb).
            !! mode 0: separated; 1: same panel (oa = xi - eta);
            !! 2: pb follows pa (oa = 1 - xi, ob = eta);
            !! 3: pb precedes pa (oa = xi, ob = 1 - eta).
            integer, intent(in) :: pa, pb, mode
            real(dp), intent(in) :: xi, eta, wt, oa, ob
            real(dp) :: x(2), y(2), ta(2), tb(2), spx, spy, nxv(2), nyv(2)
            real(dp) :: dr(2), g, kd, fd, wc, wd, dsarea, psi_x(3), dpsi_x(3)
            real(dp) :: psi_y(3), dpsi_y(3), phi_x(nl), phi_y(nl)
            real(dp) :: d_shape(2), dxs(3), dys(3), ent
            real(dp) :: ca(2), cb(2), c2a(2), c2b(2), rel_x(2), rel_y(2)
            integer :: i, j, ii, jj, dofx(3), dofy(3)

            call gs_bem_panel_point(curve, pa, xi, x, ta, spx)
            call gs_bem_panel_point(curve, pb, eta, y, tb, spy)
            nxv = [ta(2), -ta(1)]
            nyv = [tb(2), -tb(1)]
            select case (mode)
            case (1)
                d_shape = curve%geom(:, 2, pa) &
                    + curve%geom(:, 3, pa)*(xi + eta)
                dr = oa*d_shape
            case (2)
                ca = curve%geom(:, 2, pa)
                c2a = curve%geom(:, 3, pa)
                cb = curve%geom(:, 2, pb)
                c2b = curve%geom(:, 3, pb)
                rel_x = -(ca + 2.0_dp*c2a)*oa + c2a*oa*oa
                rel_y = cb*ob + c2b*ob*ob
                dr = rel_x - rel_y
            case (3)
                ca = curve%geom(:, 2, pa)
                c2a = curve%geom(:, 3, pa)
                cb = curve%geom(:, 2, pb)
                c2b = curve%geom(:, 3, pb)
                rel_x = ca*oa + c2a*oa*oa
                rel_y = -(cb + 2.0_dp*c2b)*ob + c2b*ob*ob
                dr = rel_x - rel_y
            case default
                dr = x - y
            end select
            call kernel_values(n, x(1), x(2), y(1), y(2), dr(1), dr(2), nxv, &
                nyv, want_w, g, kd, fd, wc, wd)
            dsarea = wt*spx*spy
            call p2_basis(xi, psi_x, dpsi_x)
            call p2_basis(eta, psi_y, dpsi_y)
            call flux_basis(flux_degree, xi, phi_x)
            call flux_basis(flux_degree, eta, phi_y)
            do i = 1, 3
                dofx(i) = gs_bem_trace_node_index(curve, pa, i)
                dofy(i) = gs_bem_trace_node_index(curve, pb, i)
            end do
            do j = 1, nl
                jj = (pb - 1)*nl + j
                do i = 1, nl
                    ii = (pa - 1)*nl + i
                    ops%v(ii, jj) = ops%v(ii, jj) + dsarea*g*phi_x(i)*phi_y(j)
                end do
                do i = 1, 3
                    ops%kt(dofx(i), jj) = ops%kt(dofx(i), jj) &
                        + dsarea*fd*psi_x(i)*phi_y(j)
                end do
            end do
            do j = 1, 3
                do i = 1, nl
                    ii = (pa - 1)*nl + i
                    ops%k(ii, dofy(j)) = ops%k(ii, dofy(j)) &
                        + dsarea*kd*phi_x(i)*psi_y(j)
                end do
                do i = 1, 3
                    ops%kxx(dofx(i), dofy(j)) = ops%kxx(dofx(i), dofy(j)) &
                        + dsarea*kd*psi_x(i)*psi_y(j)
                end do
            end do
            if (want_w) then
                dxs = dpsi_x/spx
                dys = dpsi_y/spy
                do j = 1, 3
                    do i = 1, 3
                        ! test function v = psi_x(i) at x, trial u = psi_y(j)
                        ! at y; see the module header for the form.
                        ent = nn*nn*psi_x(i)*psi_y(j)*(tb(1)*ta(1)*wc &
                            + tb(2)*ta(2)*g) &
                            + x(1)*y(1)*dys(j)*dxs(i)*wc &
                            + nn*wd*(x(1)*psi_y(j)*dxs(i)*tb(1) &
                            + y(1)*dys(j)*psi_x(i)*ta(1))
                        ops%w(dofx(i), dofy(j)) = ops%w(dofx(i), dofy(j)) &
                            + dsarea*ent
                    end do
                end do
            end if
        end subroutine add_point

    end subroutine tm_bem_assemble

    ! ------------------------------------------------------------------
    ! Steklov-Poincare operator and friends
    ! ------------------------------------------------------------------

    subroutine tm_bem_steklov_symmetric(ops, s, status)
        !! S_h = W + (M/2 - K)^T V^{-1} (M/2 - K), the symmetric (Costabel)
        !! form of the exterior Steklov-Poincare operator S_n.
        type(gs_bem_operators_t), intent(in) :: ops
        real(dp), allocatable, intent(out) :: s(:, :)
        integer, intent(out) :: status

        call gs_bem_steklov_symmetric(ops, s, status)
    end subroutine tm_bem_steklov_symmetric

    subroutine tm_bem_steklov_jn(ops, s, status)
        !! Johnson-Nedelec form S_h = M^T V^{-1} (M/2 - K); needs no W and
        !! is not symmetric.
        type(gs_bem_operators_t), intent(in) :: ops
        real(dp), allocatable, intent(out) :: s(:, :)
        integer, intent(out) :: status

        call gs_bem_steklov_jn(ops, s, status)
    end subroutine tm_bem_steklov_jn

    subroutine tm_bem_exterior_neumann(ops, u, t, status)
        !! Conormal flux t = R d_n u of the exterior-decaying extension of
        !! the trace u (n into the exterior), t = V^{-1} (K - M/2) u.
        type(gs_bem_operators_t), intent(in) :: ops
        real(dp), intent(in) :: u(:)
        real(dp), allocatable, intent(out) :: t(:)
        integer, intent(out) :: status

        call gs_bem_exterior_neumann(ops, u, t, status)
    end subroutine tm_bem_exterior_neumann

    subroutine tm_bem_neumann_to_dirichlet(s, nd, status)
        !! N_h = S_h^{-1} for the symmetric positive definite S_h.
        real(dp), intent(in) :: s(:, :)
        real(dp), allocatable, intent(out) :: nd(:, :)
        integer, intent(out) :: status

        call gs_bem_neumann_to_dirichlet(s, nd, status)
    end subroutine tm_bem_neumann_to_dirichlet

end module fortfem_toroidal_mode_boundary_operators_2d
