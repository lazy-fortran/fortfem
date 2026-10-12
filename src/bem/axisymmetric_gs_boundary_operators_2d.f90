module fortfem_axisymmetric_gs_boundary_operators_2d
    !! Galerkin boundary operators of the axisymmetric Grad-Shafranov operator
    !! L u = -div((1/R) grad u) on a closed meridional curve Gamma (R > 0)
    !! made of curved P2 panels, and the exterior Steklov-Poincare map.
    !!
    !! Conventions.  The curve is traversed counter-clockwise, the unit
    !! normal n = (t_Z, -t_R) points out of the enclosed region (into the
    !! exterior), ds is arc length and the conormal derivative is
    !! (1/R) d_n.  With the kernel G of fortfem_axisymmetric_gs_kernel
    !! (L G = delta, G decays at infinity and vanishes on the axis):
    !!
    !!   V t (x) = int G(x,y) t(y) ds_y
    !!   K u (x) = int (1/R_y) d_{n_y} G(x,y) u(y) ds_y
    !!   K' t(x) = int (1/R_x) d_{n_x} G(x,y) t(y) ds_y
    !!   W u (x) = -(1/R_x) d_{n_x} int (1/R_y) d_{n_y} G(x,y) u(y) ds_y
    !!
    !! so that K 1 = -1/2 and V W = 1/4 - K^2 hold (Calderon).  Because
    !! G decays, V is symmetric positive definite and the exterior Dirichlet
    !! problem is uniquely solvable.  For an exterior-harmonic field u with
    !! trace u and t = (1/R) d_n u (n into the exterior)
    !!
    !!   V t = (K - 1/2) u,   S u := -t = V^{-1} (1/2 - K) u,
    !!
    !! and <S u, u> is the exterior energy int (1/R)|grad u|^2.
    !!
    !! Spaces.  Traces u are continuous piecewise quadratic on the panel
    !! parameter (nodes: panel ends and mid point, 2 np dofs).  Fluxes t are
    !! discontinuous piecewise polynomials of degree flux_degree
    !! (flux_degree + 1 dofs per panel).
    !!
    !! Hypersingular form.  With a = 1/R and the conjugate (magnetic scalar)
    !! potential H of G, d_y H = a(y) rot grad_y G, one finds for the
    !! Galerkin form
    !!
    !!   <W u, v> = int int G(x,y) (a^2 v)'(x) u'(y) ds ds
    !!            + int int n_x . delta(x,y) v(x) u'(y) ds ds,
    !!
    !! where ' is the arc-length derivative and delta = V - a(x)^2 rot
    !! grad_x G with V_1 = (Z_x - Z_y) E(m)/(pi s |x-y|^2) and
    !! V_2 = -a(x) a(y) d_{R_y} G (E complete elliptic integral of the
    !! second kind, m = 4 R_x R_y/s^2, s^2 = (R_x+R_y)^2 + (Z_x-Z_y)^2).
    !! delta has only a logarithmic singularity, so both integrals are
    !! weakly singular.  The planar limit is the Maue form d_s V d_s.
    !!
    !! Quadrature.  Panel pairs that are well separated use composite Gauss
    !! rules refined by the distance.  On the same panel and on neighbouring
    !! panels the pair square is mapped by Duffy transformations onto the
    !! singular set, with a geometrically graded rule toward it.  On these
    !! pairs the kernel is evaluated through the logarithmic split of the
    !! kernel module, G = -ln(dist) coeff + smooth, with the geometric
    !! difference x - y formed without cancellation from the panel
    !! polynomials (x(s) - x(t) = (s - t) (c1 + c2 (s + t)) on one panel and
    !! offsets from the shared node on neighbouring panels).
    use fortfem_kinds, only: dp
    use fortfem_axisymmetric_gs_kernel, only: gs_green_gradient, &
        gs_green_log_split_gradient
    use fortnum_quadrature, only: gauss_legendre
    use fortnum_special_elliptic, only: ellip_e_from_complement
    implicit none
    private

    real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp
    !> Largest complement q = |x-y|^2/s^2 for which the log split is used.
    real(dp), parameter :: split_q_max = 0.4_dp
    integer, parameter :: radial_layers = 14, radial_points = 8
    integer, parameter :: angular_points = 12
    real(dp), parameter :: grading = 0.2_dp

    !> Closed curve of np curved P2 panels.  geom(:, k+1, p) are the
    !> polynomial coefficients of x_p(xi) = sum_k geom(:, k+1, p) xi^k.
    type, public :: gs_bem_curve_t
        integer :: n_panels = 0
        real(dp), allocatable :: nodes(:, :)       !! (2, 2 np) end, mid, ...
        real(dp), allocatable :: geom(:, :, :)     !! (2, 3, np)
    end type gs_bem_curve_t

    !> Dense Galerkin matrices.  X: P2 traces (nx = 2 np), Y: fluxes.
    type, public :: gs_bem_operators_t
        integer :: nx = 0, ny = 0, flux_degree = 1
        real(dp), allocatable :: v(:, :)      !! (ny, ny)  single layer
        real(dp), allocatable :: k(:, :)      !! (ny, nx)  <K u, phi>
        real(dp), allocatable :: kt(:, :)     !! (nx, ny)  <K' t, psi>
        real(dp), allocatable :: w(:, :)      !! (nx, nx)  hypersingular
        real(dp), allocatable :: m(:, :)      !! (ny, nx)  <psi, phi>
        real(dp), allocatable :: mxx(:, :)    !! (nx, nx)  <psi, psi>
        real(dp), allocatable :: kxx(:, :)    !! (nx, nx)  <K u, psi>
    end type gs_bem_operators_t

    public :: gs_bem_curve_from_nodes
    public :: gs_bem_assemble
    public :: gs_bem_steklov_symmetric
    public :: gs_bem_steklov_jn
    public :: gs_bem_neumann_to_dirichlet
    public :: gs_bem_dual_norm_squared
    public :: gs_bem_exterior_neumann
    public :: gs_bem_panel_point
    public :: gs_bem_trace_value
    public :: gs_bem_flux_value
    public :: gs_bem_trace_node_index

    interface
        subroutine dpotrf(uplo, n, a, lda, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine dpotrf
        subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dpotrs
        subroutine dpotri(uplo, n, a, lda, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine dpotri
    end interface

contains

    ! ------------------------------------------------------------------
    ! Geometry and basis functions
    ! ------------------------------------------------------------------

    subroutine gs_bem_curve_from_nodes(nodes, curve, status)
        !! Panel p has the nodes 2p-1 (start), 2p (middle) and 2p+1 (end,
        !! the start of panel p+1; the last panel closes on node 1).  The
        !! curve must be counter-clockwise with at least 4 panels.
        real(dp), intent(in) :: nodes(:, :)
        type(gs_bem_curve_t), intent(out) :: curve
        integer, intent(out) :: status
        real(dp) :: x0(2), xm(2), x1(2)
        integer :: np, p, p_end

        status = 1
        if (size(nodes, 1) /= 2) return
        if (mod(size(nodes, 2), 2) /= 0) return
        np = size(nodes, 2)/2
        if (np < 4) return
        if (any(nodes(1, :) <= 0.0_dp)) return
        curve%n_panels = np
        allocate(curve%nodes(2, 2*np), curve%geom(2, 3, np))
        curve%nodes = nodes
        do p = 1, np
            p_end = 2*p + 1
            if (p == np) p_end = 1
            x0 = nodes(:, 2*p - 1)
            xm = nodes(:, 2*p)
            x1 = nodes(:, p_end)
            curve%geom(:, 3, p) = 2.0_dp*(x0 + x1 - 2.0_dp*xm)
            curve%geom(:, 2, p) = x1 - x0 - curve%geom(:, 3, p)
            curve%geom(:, 1, p) = x0
        end do
        status = 0
    end subroutine gs_bem_curve_from_nodes

    pure function gs_bem_trace_node_index(curve, p, k) result(idx)
        !! Global trace dof of local node k = 1, 2, 3 of panel p.
        type(gs_bem_curve_t), intent(in) :: curve
        integer, intent(in) :: p, k
        integer :: idx

        idx = 2*p - 2 + k
        if (idx > 2*curve%n_panels) idx = 1
    end function gs_bem_trace_node_index

    pure subroutine gs_bem_panel_point(curve, p, xi, x, tau, speed)
        !! Position, unit tangent and |dx/dxi| on panel p at xi in [0, 1].
        !! The outward normal is (tau(2), -tau(1)).
        type(gs_bem_curve_t), intent(in) :: curve
        integer, intent(in) :: p
        real(dp), intent(in) :: xi
        real(dp), intent(out) :: x(2), tau(2), speed
        real(dp) :: t(2)

        x = curve%geom(:, 1, p) + xi*(curve%geom(:, 2, p) &
            + xi*curve%geom(:, 3, p))
        t = curve%geom(:, 2, p) + 2.0_dp*xi*curve%geom(:, 3, p)
        speed = norm2(t)
        tau = t/speed
    end subroutine gs_bem_panel_point

    pure subroutine p2_basis(xi, psi, dpsi)
        real(dp), intent(in) :: xi
        real(dp), intent(out) :: psi(3), dpsi(3)

        psi = [(1.0_dp - xi)*(1.0_dp - 2.0_dp*xi), 4.0_dp*xi*(1.0_dp - xi), &
            xi*(2.0_dp*xi - 1.0_dp)]
        dpsi = [4.0_dp*xi - 3.0_dp, 4.0_dp - 8.0_dp*xi, 4.0_dp*xi - 1.0_dp]
    end subroutine p2_basis

    pure subroutine flux_basis(deg, xi, phi)
        !! Lagrange basis of degree deg on equispaced nodes of [0, 1].
        integer, intent(in) :: deg
        real(dp), intent(in) :: xi
        real(dp), intent(out) :: phi(deg + 1)
        real(dp) :: xk, xj
        integer :: j, k

        do k = 0, deg
            phi(k + 1) = 1.0_dp
            xk = real(k, dp)/max(deg, 1)
            do j = 0, deg
                if (j == k) cycle
                xj = real(j, dp)/max(deg, 1)
                phi(k + 1) = phi(k + 1)*(xi - xj)/(xk - xj)
            end do
        end do
    end subroutine flux_basis

    pure function gs_bem_trace_value(curve, p, xi, u) result(val)
        type(gs_bem_curve_t), intent(in) :: curve
        integer, intent(in) :: p
        real(dp), intent(in) :: xi, u(:)
        real(dp) :: val
        real(dp) :: psi(3), dpsi(3)
        integer :: k

        call p2_basis(xi, psi, dpsi)
        val = 0.0_dp
        do k = 1, 3
            val = val + psi(k)*u(gs_bem_trace_node_index(curve, p, k))
        end do
    end function gs_bem_trace_value

    pure function gs_bem_flux_value(deg, p, xi, t) result(val)
        integer, intent(in) :: deg, p
        real(dp), intent(in) :: xi, t(:)
        real(dp) :: val
        real(dp) :: phi(deg + 1)

        call flux_basis(deg, xi, phi)
        val = dot_product(phi, t((p - 1)*(deg + 1) + 1:p*(deg + 1)))
    end function gs_bem_flux_value

    ! ------------------------------------------------------------------
    ! Kernel values at one pair of points
    ! ------------------------------------------------------------------

    pure subroutine kernel_values(r, z, rp, zp, dr, dz, nx, ny, need_w, &
            g, kd, fd, dn)
        !! G, K-kernel kd = (1/R_y) n_y . grad_y G, K'-kernel
        !! fd = (1/R_x) n_x . grad_x G and, with need_w, dn = n_x . delta.
        !! (dr, dz) is the accurately formed difference x - y.
        real(dp), intent(in) :: r, z, rp, zp, dr, dz, nx(2), ny(2)
        logical, intent(in) :: need_w
        real(dp), intent(out) :: g, kd, fd, dn
        real(dp) :: d2, s2, q, c, dc(4), sm, dsm(4), lnd, gx(2), gy(2)
        real(dp) :: grad(4), ax, ay, s, e, v1, v2, e1, e2

        d2 = dr*dr + dz*dz
        s2 = (r + rp)**2 + dz*dz
        q = d2/s2
        if (q < split_q_max) then
            call gs_green_log_split_gradient(r, z, rp, zp, c, dc, sm, dsm)
            lnd = 0.5_dp*log(d2)
            g = -c*lnd + sm
            gx(1) = -c*dr/d2 - lnd*dc(1) + dsm(1)
            gx(2) = -c*dz/d2 - lnd*dc(2) + dsm(2)
            gy(1) = c*dr/d2 - lnd*dc(3) + dsm(3)
            gy(2) = c*dz/d2 - lnd*dc(4) + dsm(4)
        else
            call gs_green_gradient(r, z, rp, zp, g, grad)
            gx = grad(1:2)
            gy = grad(3:4)
        end if
        ax = 1.0_dp/r
        ay = 1.0_dp/rp
        fd = ax*(nx(1)*gx(1) + nx(2)*gx(2))
        kd = ay*(ny(1)*gy(1) + ny(2)*gy(2))
        dn = 0.0_dp
        if (need_w) then
            s = sqrt(s2)
            e = ellip_e_from_complement(q)
            v1 = dz*e/(pi*s*d2)
            v2 = -ax*ay*gy(1)
            e1 = v1 + ax*ax*gx(2)
            e2 = v2 - ax*ax*gx(1)
            dn = nx(1)*e1 + nx(2)*e2
        end if
    end subroutine kernel_values

    ! ------------------------------------------------------------------
    ! Quadrature rules on the unit interval and pair rules
    ! ------------------------------------------------------------------

    subroutine gauss_unit(n, x, w)
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)

        call gauss_legendre(n, x, w)
        x = 0.5_dp*(x + 1.0_dp)
        w = 0.5_dp*w
    end subroutine gauss_unit

    subroutine graded_rule(x, w)
        !! Composite Gauss rule on [0, 1] graded geometrically toward 0.
        real(dp), allocatable, intent(out) :: x(:), w(:)
        real(dp) :: gx(radial_points), gw(radial_points), lo, hi
        integer :: layer, k

        allocate(x((radial_layers + 1)*radial_points), &
            w((radial_layers + 1)*radial_points))
        call gauss_unit(radial_points, gx, gw)
        hi = 1.0_dp
        k = 0
        do layer = 1, radial_layers + 1
            if (layer == radial_layers + 1) then
                lo = 0.0_dp
            else
                lo = hi*grading
            end if
            x(k + 1:k + radial_points) = lo + (hi - lo)*gx
            w(k + 1:k + radial_points) = (hi - lo)*gw
            k = k + radial_points
            hi = lo
        end do
    end subroutine graded_rule

    ! ------------------------------------------------------------------
    ! Assembly
    ! ------------------------------------------------------------------

    subroutine gs_bem_assemble(curve, flux_degree, ops, status, hypersingular)
        !! Assemble V, K, K', W and the mass matrices.  W (the most expensive
        !! kernel) is skipped when hypersingular = .false.
        type(gs_bem_curve_t), intent(in) :: curve
        integer, intent(in) :: flux_degree
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
        real(dp) :: cen(2, curve%n_panels)
        integer :: np, nx, ny, nl, a, b, kind
        logical :: want_w

        status = 1
        want_w = .true.
        if (present(hypersingular)) want_w = hypersingular
        if (curve%n_panels < 4 .or. flux_degree < 0 .or. flux_degree > 4) return
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

        call panel_summary(lens, cen)
        call assemble_mass()
        do a = 1, np
            do b = 1, np
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
        status = 0

    contains

        subroutine panel_summary(len_p, cen_p)
            real(dp), intent(out) :: len_p(:), cen_p(:, :)
            real(dp) :: x0(2), x1(2), tau(2), sp
            integer :: p

            do p = 1, np
                call gs_bem_panel_point(curve, p, 0.0_dp, x0, tau, sp)
                call gs_bem_panel_point(curve, p, 1.0_dp, x1, tau, sp)
                len_p(p) = norm2(x1 - x0)
                call gs_bem_panel_point(curve, p, 0.5_dp, cen_p(:, p), tau, sp)
            end do
        end subroutine panel_summary

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
            integer :: ng, m, i, j, ia, ib, n

            select case (pair_kind)
            case (1)
                do n = 1, size(self_xi)
                    call add_point(pa, pb, self_xi(n), self_eta(n), &
                        self_w(n), 1, self_oa(n), 0.0_dp)
                end do
            case (2)
                do n = 1, size(nxt_xi)
                    call add_point(pa, pb, nxt_xi(n), nxt_eta(n), nxt_w(n), &
                        2, nxt_oa(n), nxt_ob(n))
                end do
            case (3)
                do n = 1, size(prv_xi)
                    call add_point(pa, pb, prv_xi(n), prv_eta(n), prv_w(n), &
                        3, prv_oa(n), prv_ob(n))
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
            real(dp) :: dr(2), g, kd, fd, dn, dsarea, psi_x(3), dpsi_x(3)
            real(dp) :: psi_y(3), dpsi_y(3), phi_x(nl), phi_y(nl)
            real(dp) :: ax_s, awx(3), dpy(3), d_shape(2), dd(2)
            integer :: i, j, ii, jj, dofx(3), dofy(3)
            real(dp) :: ca(2), cb(2), c2a(2), c2b(2), rel_x(2), rel_y(2)

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
            dd = dr
            call kernel_values(x(1), x(2), y(1), y(2), dd(1), dd(2), nxv, nyv, &
                want_w, g, kd, fd, dn)
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
                ! (a^2 psi)'(x) = psi_s/R^2 - 2 tau_R psi/R^3
                ax_s = 1.0_dp/x(1)
                awx = ax_s*ax_s*dpsi_x/spx &
                    - 2.0_dp*ta(1)*ax_s**3*psi_x
                dpy = dpsi_y/spy
                do j = 1, 3
                    do i = 1, 3
                        ops%w(dofx(i), dofy(j)) = ops%w(dofx(i), dofy(j)) &
                            + dsarea*dpy(j)*(g*awx(i) + dn*psi_x(i))
                    end do
                end do
            end if
        end subroutine add_point

    end subroutine gs_bem_assemble

    subroutine build_self_rule(rx, rw, xi, eta, w, oa)
        !! Same-panel rule: the two triangles xi > eta and xi < eta, each
        !! parametrised by rho = |xi - eta| (graded) and the position v.
        real(dp), intent(in) :: rx(:), rw(:)
        real(dp), allocatable, intent(out) :: xi(:), eta(:), w(:), oa(:)
        real(dp) :: vx(angular_points), vw(angular_points), e0, rho
        integer :: n, k, j

        n = 2*size(rx)*angular_points
        allocate(xi(n), eta(n), w(n), oa(n))
        call gauss_unit(angular_points, vx, vw)
        n = 0
        do k = 1, size(rx)
            rho = rx(k)
            do j = 1, angular_points
                e0 = (1.0_dp - rho)*vx(j)
                n = n + 1
                xi(n) = e0 + rho
                eta(n) = e0
                w(n) = rw(k)*vw(j)*(1.0_dp - rho)
                oa(n) = rho
                n = n + 1
                xi(n) = e0
                eta(n) = e0 + rho
                w(n) = rw(k)*vw(j)*(1.0_dp - rho)
                oa(n) = -rho
            end do
        end do
    end subroutine build_self_rule

    subroutine build_corner_rule(rx, rw, follows, xi, eta, w, oa, ob)
        !! Neighbouring panels: Duffy rule around the shared corner.  With
        !! follows, the second panel starts where the first ends; otherwise
        !! it ends where the first starts.  oa and ob are the parameter
        !! distances of the two points from the corner.
        real(dp), intent(in) :: rx(:), rw(:)
        logical, intent(in) :: follows
        real(dp), allocatable, intent(out) :: xi(:), eta(:), w(:), oa(:), ob(:)
        real(dp) :: vx(angular_points), vw(angular_points), s, t
        integer :: n, k, j, region

        n = 2*size(rx)*angular_points
        allocate(xi(n), eta(n), w(n), oa(n), ob(n))
        call gauss_unit(angular_points, vx, vw)
        n = 0
        do k = 1, size(rx)
            do j = 1, angular_points
                do region = 1, 2
                    if (region == 1) then
                        t = rx(k)
                        s = rx(k)*vx(j)
                    else
                        s = rx(k)
                        t = rx(k)*vx(j)
                    end if
                    n = n + 1
                    w(n) = rw(k)*vw(j)*rx(k)
                    oa(n) = s
                    ob(n) = t
                    if (follows) then
                        xi(n) = 1.0_dp - s
                        eta(n) = t
                    else
                        xi(n) = s
                        eta(n) = 1.0_dp - t
                    end if
                end do
            end do
        end do
    end subroutine build_corner_rule

    ! ------------------------------------------------------------------
    ! Steklov-Poincare operator and friends
    ! ------------------------------------------------------------------

    subroutine v_solve(ops, rhs, status)
        !! rhs <- V^{-1} rhs (columns), V symmetric positive definite.
        type(gs_bem_operators_t), intent(in) :: ops
        real(dp), intent(inout) :: rhs(:, :)
        integer, intent(out) :: status
        real(dp), allocatable :: chol(:, :)
        integer :: info

        status = 1
        chol = 0.5_dp*(ops%v + transpose(ops%v))
        call dpotrf('L', ops%ny, chol, ops%ny, info)
        if (info /= 0) return
        call dpotrs('L', ops%ny, size(rhs, 2), chol, ops%ny, rhs, ops%ny, info)
        if (info /= 0) return
        status = 0
    end subroutine v_solve

    subroutine gs_bem_steklov_symmetric(ops, s, status)
        !! S_h = W + (M/2 - K)^T V^{-1} (M/2 - K), the symmetric (Costabel)
        !! form of the exterior Steklov-Poincare operator.
        type(gs_bem_operators_t), intent(in) :: ops
        real(dp), allocatable, intent(out) :: s(:, :)
        integer, intent(out) :: status
        real(dp), allocatable :: b(:, :), c(:, :)

        b = 0.5_dp*ops%m - ops%k
        c = b
        call v_solve(ops, c, status)
        if (status /= 0) return
        s = ops%w + matmul(transpose(b), c)
        s = 0.5_dp*(s + transpose(s))
    end subroutine gs_bem_steklov_symmetric

    subroutine gs_bem_steklov_jn(ops, s, status)
        !! Johnson-Nedelec form S_h = M^T V^{-1} (M/2 - K); it needs no W and
        !! is not symmetric.
        type(gs_bem_operators_t), intent(in) :: ops
        real(dp), allocatable, intent(out) :: s(:, :)
        integer, intent(out) :: status
        real(dp), allocatable :: c(:, :)

        c = 0.5_dp*ops%m - ops%k
        call v_solve(ops, c, status)
        if (status /= 0) return
        s = matmul(transpose(ops%m), c)
    end subroutine gs_bem_steklov_jn

    subroutine gs_bem_exterior_neumann(ops, u, t, status)
        !! Conormal flux t = (1/R) d_n u of the exterior-decaying extension
        !! of the trace u (n into the exterior), t = V^{-1} (K - M/2) u, as
        !! flux coefficients.
        type(gs_bem_operators_t), intent(in) :: ops
        real(dp), intent(in) :: u(:)
        real(dp), allocatable, intent(out) :: t(:)
        integer, intent(out) :: status
        real(dp), allocatable :: rhs(:, :)

        allocate(rhs(ops%ny, 1))
        rhs(:, 1) = matmul(ops%k - 0.5_dp*ops%m, u)
        call v_solve(ops, rhs, status)
        t = rhs(:, 1)
    end subroutine gs_bem_exterior_neumann

    subroutine gs_bem_neumann_to_dirichlet(s, n, status)
        !! N_h = S_h^{-1} for the symmetric positive definite S_h: the trace
        !! coefficients produced by the load vector mu (functional on the
        !! trace space).  mu^T N_h mu is the squared dual norm.
        real(dp), intent(in) :: s(:, :)
        real(dp), allocatable, intent(out) :: n(:, :)
        integer, intent(out) :: status
        integer :: info, i, j, m

        status = 1
        m = size(s, 1)
        n = 0.5_dp*(s + transpose(s))
        call dpotrf('L', m, n, m, info)
        if (info /= 0) return
        call dpotri('L', m, n, m, info)
        if (info /= 0) return
        do j = 1, m
            do i = j + 1, m
                n(j, i) = n(i, j)
            end do
        end do
        status = 0
    end subroutine gs_bem_neumann_to_dirichlet

    pure function gs_bem_dual_norm_squared(n, mu) result(val)
        real(dp), intent(in) :: n(:, :), mu(:)
        real(dp) :: val

        val = dot_product(mu, matmul(n, mu))
    end function gs_bem_dual_norm_squared

end module fortfem_axisymmetric_gs_boundary_operators_2d
