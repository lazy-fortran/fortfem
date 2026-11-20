module fortfem_gs_single_layer_shape_2d
    !! Galerkin single-layer operator of the axisymmetric Grad-Shafranov
    !! kernel G on a closed meridional curve of curved P_k panels, with its
    !! exact first and second derivatives in the panel node coordinates.
    !!
    !!   V_ij = int int G(x, y) phi_i(x) phi_j(y) ds_x ds_y
    !!   g_i  = int psi(x) phi_i(x) ds_x
    !!
    !! G and the flux basis phi (discontinuous Lagrange of degree
    !! flux_degree per panel) are those of
    !! fortfem_axisymmetric_gs_boundary_operators_2d; for k = 2 this V equals
    !! its ops%v.  Panel p has the nodes k(p-1)+1, ..., k(p-1)+k+1 at the
    !! equispaced parameters 0, 1/k, ..., 1 (the last node of the last panel
    !! is node 1).
    !!
    !! Exactness.  The quadrature is fixed in the panel parameters by a plan
    !! made once (gs_sl_plan_create): same-panel and neighbouring pairs use
    !! the Duffy rules of the operator module, separated pairs a composite
    !! Gauss rule whose subdivision the plan freezes.  At fixed parameters
    !! every integrand input (points, tangents and the cancellation-free
    !! difference of the two points) is linear in the node coordinates, so
    !! the generated value/gradient/Hessian jets of the integrand
    !! (tools/codegen gen_gs_single_layer_shape_products) give the exact
    !! derivatives of the discrete V and g.  Rebuild the plan when the
    !! geometry has moved far from the one it was made for.
    !!
    !! Coordinates.  X(2(n-1)+c) is component c (1: R, 2: Z) of node n.
    use fortfem_kinds, only: dp
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: flux_basis, &
        gauss_unit, graded_rule, build_self_rule, build_corner_rule
    use fortfem_generated_gs_single_layer_near_jet, only: &
        generated_gs_single_layer_near_jet
    use fortfem_generated_gs_single_layer_far_jet, only: &
        generated_gs_single_layer_far_jet
    use fortfem_generated_gs_boundary_load_jet, only: &
        generated_gs_boundary_load_jet
    implicit none
    private

    !> Largest q = |x-y|^2/((R+R')^2 + (Z-Z')^2) for the log-split form.
    real(dp), parameter :: split_q_max = 0.4_dp
    integer, parameter :: n_in = 10

    type, public :: gs_sl_curve_t
        integer :: degree = 0, n_panels = 0, n_nodes = 0
        real(dp), allocatable :: nodes(:, :)
        !> N_a(xi) = sum_j lagrange(a, j) xi**(j-1).
        real(dp), allocatable :: lagrange(:, :)
    end type gs_sl_curve_t

    type, public :: gs_sl_plan_t
        integer :: n_panels = 0, flux_degree = 0
        integer, allocatable :: kind(:, :), sub(:, :), order(:, :)
        real(dp), allocatable :: self_xi(:), self_eta(:), self_w(:), self_oa(:)
        real(dp), allocatable :: nxt_xi(:), nxt_eta(:), nxt_w(:), nxt_oa(:), &
            nxt_ob(:)
        real(dp), allocatable :: prv_xi(:), prv_eta(:), prv_w(:), prv_oa(:), &
            prv_ob(:)
    end type gs_sl_plan_t

    abstract interface
        subroutine gs_sl_field(x, value, gradient, hessian)
            !! A field and its jets at x = (R, Z); hessian = (RR, RZ, ZZ).
            import :: dp
            real(dp), intent(in) :: x(2)
            real(dp), intent(out) :: value, gradient(2), hessian(3)
        end subroutine gs_sl_field
    end interface

    public :: gs_sl_field
    public :: gs_sl_curve_create, gs_sl_plan_create
    public :: gs_sl_node, gs_sl_point
    public :: gs_sl_assemble, gs_sl_load

contains

    subroutine gs_sl_curve_create(degree, nodes, curve, status)
        integer, intent(in) :: degree
        real(dp), intent(in) :: nodes(:, :)
        type(gs_sl_curve_t), intent(out) :: curve
        integer, intent(out) :: status
        real(dp) :: vander(degree + 1, degree + 1), xa
        integer :: a, j

        status = 1
        if (degree < 1 .or. degree > 4) return
        if (size(nodes, 1) /= 2) return
        if (mod(size(nodes, 2), degree) /= 0) return
        if (size(nodes, 2)/degree < 4) return
        if (any(nodes(1, :) <= 0.0_dp)) return
        curve%degree = degree
        curve%n_nodes = size(nodes, 2)
        curve%n_panels = curve%n_nodes/degree
        curve%nodes = nodes
        do a = 1, degree + 1
            xa = real(a - 1, dp)/degree
            do j = 1, degree + 1
                vander(a, j) = xa**(j - 1)
            end do
        end do
        ! Columns of the inverse Vandermonde are the monomial coefficients
        ! of the Lagrange basis.
        call invert_small(vander, status)
        if (status /= 0) return
        curve%lagrange = transpose(vander)
        status = 0
    end subroutine gs_sl_curve_create

    pure function gs_sl_node(curve, p, a) result(n)
        !! Global node of local node a = 1, ..., k+1 of panel p.
        type(gs_sl_curve_t), intent(in) :: curve
        integer, intent(in) :: p, a
        integer :: n

        n = (p - 1)*curve%degree + a
        if (n > curve%n_nodes) n = n - curve%n_nodes
    end function gs_sl_node

    pure subroutine gs_sl_point(curve, p, xi, x, t)
        !! Position and (unnormalised) tangent dx/dxi of panel p at xi.
        type(gs_sl_curve_t), intent(in) :: curve
        integer, intent(in) :: p
        real(dp), intent(in) :: xi
        real(dp), intent(out) :: x(2), t(2)
        real(dp) :: n(curve%degree + 1), dn(curve%degree + 1)
        integer :: a

        call basis(curve, xi, n, dn)
        x = 0.0_dp
        t = 0.0_dp
        do a = 1, curve%degree + 1
            x = x + n(a)*curve%nodes(:, gs_sl_node(curve, p, a))
            t = t + dn(a)*curve%nodes(:, gs_sl_node(curve, p, a))
        end do
    end subroutine gs_sl_point

    subroutine gs_sl_plan_create(curve, flux_degree, plan, status)
        !! Pair kinds (0 separated, 1 same, 2 second follows first, 3 second
        !! precedes first) and the separated-pair subdivision m and Gauss
        !! order, chosen from the present geometry as in gs_bem_assemble.
        type(gs_sl_curve_t), intent(in) :: curve
        integer, intent(in) :: flux_degree
        type(gs_sl_plan_t), intent(out) :: plan
        integer, intent(out) :: status
        real(dp), allocatable :: rad_x(:), rad_w(:)
        real(dp) :: lens(curve%n_panels), x0(2), x1(2), t(2), hmax, dmin
        integer :: np, a, b

        status = 1
        if (flux_degree < 0 .or. flux_degree > 4) return
        np = curve%n_panels
        plan%n_panels = np
        plan%flux_degree = flux_degree
        allocate (plan%kind(np, np), plan%sub(np, np), plan%order(np, np))
        do a = 1, np
            call gs_sl_point(curve, a, 0.0_dp, x0, t)
            call gs_sl_point(curve, a, 1.0_dp, x1, t)
            lens(a) = norm2(x1 - x0)
        end do
        do b = 1, np
            do a = 1, np
                plan%sub(a, b) = 0
                plan%order(a, b) = 0
                if (b == a) then
                    plan%kind(a, b) = 1
                else if (b == mod(a, np) + 1) then
                    plan%kind(a, b) = 2
                else if (a == mod(b, np) + 1) then
                    plan%kind(a, b) = 3
                else
                    plan%kind(a, b) = 0
                    hmax = max(lens(a), lens(b))
                    dmin = pair_distance(curve, a, b)
                    plan%sub(a, b) = max(1, min(12, &
                        ceiling(2.5_dp*hmax/max(dmin, 1.0e-14_dp))))
                    plan%order(a, b) = 8
                    if (dmin > 8.0_dp*hmax) plan%order(a, b) = 6
                end if
            end do
        end do
        call graded_rule(rad_x, rad_w)
        call build_self_rule(rad_x, rad_w, plan%self_xi, plan%self_eta, &
            plan%self_w, plan%self_oa)
        call build_corner_rule(rad_x, rad_w, .true., plan%nxt_xi, plan%nxt_eta, &
            plan%nxt_w, plan%nxt_oa, plan%nxt_ob)
        call build_corner_rule(rad_x, rad_w, .false., plan%prv_xi, &
            plan%prv_eta, plan%prv_w, plan%prv_oa, plan%prv_ob)
        status = 0
    end subroutine gs_sl_plan_create

    subroutine gs_sl_assemble(curve, plan, v, tau, dvtau, gradient, hessian)
        !! V; with tau also dvtau = d(V tau)/dX, the gradient and the Hessian
        !! of tau' V tau / 2 in X.
        type(gs_sl_curve_t), intent(in) :: curve
        type(gs_sl_plan_t), intent(in) :: plan
        real(dp), intent(out) :: v(:, :)
        real(dp), intent(in), optional :: tau(:)
        real(dp), intent(out), optional :: dvtau(:, :), gradient(:), hessian(:, :)
        real(dp), allocatable :: fx(:), fw(:), lv(:, :), ldv(:, :), lg(:), lh(:, :)
        real(dp), allocatable :: tau_l(:)
        logical :: shape
        integer :: np, a, b, n, i, j, m, ng, nx, ny, ns

        np = curve%n_panels
        nx = 2*curve%n_nodes
        ny = size(v, 1)
        shape = present(tau)
        ns = 0
        if (shape) ns = nx
        allocate (tau_l(ny), source=0.0_dp)
        if (shape) tau_l = tau
        allocate (lv(ny, ny), source=0.0_dp)
        allocate (ldv(ny, ns), lg(ns), lh(ns, ns), source=0.0_dp)
        !$omp parallel do schedule(dynamic, 1) default(shared) &
        !$omp& private(a, b, n, i, j, m, ng, fx, fw) &
        !$omp& reduction(+:lv, ldv, lg, lh)
        do b = 1, np
            do a = 1, np
                select case (plan%kind(a, b))
                case (1)
                    do n = 1, size(plan%self_xi)
                        call add_pair(curve, plan, shape, tau_l, a, b, &
                            plan%self_xi(n), plan%self_eta(n), plan%self_w(n), &
                            1, plan%self_oa(n), 0.0_dp, lv, ldv, lg, lh)
                    end do
                case (2)
                    do n = 1, size(plan%nxt_xi)
                        call add_pair(curve, plan, shape, tau_l, a, b, &
                            plan%nxt_xi(n), plan%nxt_eta(n), plan%nxt_w(n), &
                            2, plan%nxt_oa(n), plan%nxt_ob(n), lv, ldv, lg, lh)
                    end do
                case (3)
                    do n = 1, size(plan%prv_xi)
                        call add_pair(curve, plan, shape, tau_l, a, b, &
                            plan%prv_xi(n), plan%prv_eta(n), plan%prv_w(n), &
                            3, plan%prv_oa(n), plan%prv_ob(n), lv, ldv, lg, lh)
                    end do
                case default
                    m = plan%sub(a, b)
                    ng = plan%order(a, b)
                    if (allocated(fx)) deallocate (fx, fw)
                    allocate (fx(ng), fw(ng))
                    call gauss_unit(ng, fx, fw)
                    do i = 1, m*ng
                        do j = 1, m*ng
                            call add_pair(curve, plan, shape, tau_l, a, b, &
                                (real((i - 1)/ng, dp) + fx(mod(i - 1, ng) + 1))/m, &
                                (real((j - 1)/ng, dp) + fx(mod(j - 1, ng) + 1))/m, &
                                fw(mod(i - 1, ng) + 1)*fw(mod(j - 1, ng) + 1) &
                                /real(m*m, dp), 0, 0.0_dp, 0.0_dp, lv, ldv, lg, lh)
                        end do
                    end do
                end select
            end do
        end do
        !$omp end parallel do
        v = lv
        if (shape) then
            dvtau = ldv
            gradient = lg
            hessian = lh
        end if
    end subroutine gs_sl_assemble

    subroutine add_pair(curve, plan, shape, tau, pa, pb, xi, eta, w, mode, oa, &
            ob, v, dvtau, gradient, hessian)
        !! One quadrature pair (xi on panel pa, eta on panel pb).
        type(gs_sl_curve_t), intent(in) :: curve
        type(gs_sl_plan_t), intent(in) :: plan
        logical, intent(in) :: shape
        real(dp), intent(in) :: tau(:)
        integer, intent(in) :: pa, pb, mode
        real(dp), intent(in) :: xi, eta, w, oa, ob
        real(dp), intent(inout) :: v(:, :), dvtau(:, :), gradient(:), hessian(:, :)
        real(dp) :: jac(n_in, 4*(curve%degree + 1)), u(n_in), ell
        real(dp) :: f, gr(n_in), he(n_in, n_in)
        real(dp) :: phi_a(plan%flux_degree + 1), phi_b(plan%flux_degree + 1)
        real(dp) :: alpha, beta, gl(4*(curve%degree + 1))
        real(dp) :: hl(4*(curve%degree + 1), 4*(curve%degree + 1))
        integer :: glob(4*(curve%degree + 1)), nl, ii, jj, i, j, k, l

        call pair_inputs(curve, pa, pb, xi, eta, mode, oa, ob, u, ell, jac, glob)
        call pair_jet(u, ell, f, gr, he)
        nl = plan%flux_degree + 1
        call flux_basis(plan%flux_degree, xi, phi_a)
        call flux_basis(plan%flux_degree, eta, phi_b)
        do j = 1, nl
            jj = (pb - 1)*nl + j
            do i = 1, nl
                ii = (pa - 1)*nl + i
                v(ii, jj) = v(ii, jj) + w*f*phi_a(i)*phi_b(j)
            end do
        end do
        if (.not. shape) return
        alpha = dot_product(tau((pa - 1)*nl + 1:pa*nl), phi_a)
        beta = dot_product(tau((pb - 1)*nl + 1:pb*nl), phi_b)
        gl = matmul(gr, jac)
        hl = matmul(transpose(jac), matmul(he, jac))
        do k = 1, size(glob)
            do i = 1, nl
                ii = (pa - 1)*nl + i
                dvtau(ii, glob(k)) = dvtau(ii, glob(k)) + w*phi_a(i)*beta*gl(k)
            end do
            gradient(glob(k)) = gradient(glob(k)) + 0.5_dp*w*alpha*beta*gl(k)
            do l = 1, size(glob)
                hessian(glob(l), glob(k)) = hessian(glob(l), glob(k)) &
                    + 0.5_dp*w*alpha*beta*hl(l, k)
            end do
        end do
    end subroutine add_pair

    subroutine gs_sl_load(curve, flux_degree, field, g, weights, dg, gradient, &
            hessian)
        !! g_i = int psi phi_i ds; with weights w (one per flux dof) also
        !! dg = dg/dX and the gradient and Hessian of w'g in X.  Ten Gauss
        !! points per panel, as the mass matrices of the operator module.
        type(gs_sl_curve_t), intent(in) :: curve
        integer, intent(in) :: flux_degree
        procedure(gs_sl_field) :: field
        real(dp), intent(out) :: g(:)
        real(dp), intent(in), optional :: weights(:)
        real(dp), intent(out), optional :: dg(:, :), gradient(:), hessian(:, :)
        real(dp) :: gx(10), gw(10), x(2), t(2), psi, dpsi(2), d2psi(3)
        real(dp) :: n(curve%degree + 1), dn(curve%degree + 1)
        real(dp) :: phi(flux_degree + 1), out(15), gr(4), he(4, 4)
        real(dp) :: jac(4, 2*(curve%degree + 1)), gl(2*(curve%degree + 1))
        real(dp) :: hl(2*(curve%degree + 1), 2*(curve%degree + 1)), wsum
        integer :: glob(2*(curve%degree + 1)), p, q, a, c, i, k, l, nl, ii
        logical :: shape

        shape = present(weights)
        nl = flux_degree + 1
        g = 0.0_dp
        if (shape) then
            dg = 0.0_dp
            gradient = 0.0_dp
            hessian = 0.0_dp
        end if
        call gauss_unit(10, gx, gw)
        do p = 1, curve%n_panels
            do q = 1, 10
                call basis(curve, gx(q), n, dn)
                call gs_sl_point(curve, p, gx(q), x, t)
                call field(x, psi, dpsi, d2psi)
                call generated_gs_boundary_load_jet(t(1), t(2), psi, dpsi(1), &
                    dpsi(2), d2psi(1), d2psi(2), d2psi(3), out(1), out(2), &
                    out(3), out(4), out(5), out(6), out(7), out(8), out(9), &
                    out(10), out(11), out(12), out(13), out(14), out(15))
                call flux_basis(flux_degree, gx(q), phi)
                do i = 1, nl
                    ii = (p - 1)*nl + i
                    g(ii) = g(ii) + gw(q)*out(1)*phi(i)
                end do
                if (.not. shape) cycle
                gr = out(2:5)
                call unpack_hessian(out(6:15), he)
                jac = 0.0_dp
                do a = 1, curve%degree + 1
                    do c = 1, 2
                        k = 2*(a - 1) + c
                        glob(k) = 2*(gs_sl_node(curve, p, a) - 1) + c
                        jac(c, k) = n(a)
                        jac(2 + c, k) = dn(a)
                    end do
                end do
                gl = matmul(gr, jac)
                hl = matmul(transpose(jac), matmul(he, jac))
                wsum = dot_product(weights((p - 1)*nl + 1:p*nl), phi)
                do k = 1, size(glob)
                    do i = 1, nl
                        ii = (p - 1)*nl + i
                        dg(ii, glob(k)) = dg(ii, glob(k)) + gw(q)*phi(i)*gl(k)
                    end do
                    gradient(glob(k)) = gradient(glob(k)) + gw(q)*wsum*gl(k)
                    do l = 1, size(glob)
                        hessian(glob(l), glob(k)) = hessian(glob(l), glob(k)) &
                            + gw(q)*wsum*hl(l, k)
                    end do
                end do
            end do
        end do
    end subroutine gs_sl_load

    ! ------------------------------------------------------------------
    ! Pair inputs and integrand jets
    ! ------------------------------------------------------------------

    pure subroutine pair_inputs(curve, pa, pb, xi, eta, mode, oa, ob, u, ell, &
            jac, glob)
        !! Integrand inputs u = (R, Z, R', Z', t_x, t_y, delta) at one point
        !! pair, their constant Jacobian in the local node coordinates
        !! (panel pa nodes, then panel pb nodes) and the global coordinates.
        !! delta = x - y, except on one panel, where delta = (x - y)/(xi -
        !! eta) and ell = ln|xi - eta|.  Near the shared corner of
        !! neighbouring panels delta is formed from offsets to the corner.
        type(gs_sl_curve_t), intent(in) :: curve
        integer, intent(in) :: pa, pb, mode
        real(dp), intent(in) :: xi, eta, oa, ob
        real(dp), intent(out) :: u(n_in), ell, jac(:, :)
        integer, intent(out) :: glob(:)
        real(dp) :: na(curve%degree + 1), dna(curve%degree + 1)
        real(dp) :: nb(curve%degree + 1), dnb(curve%degree + 1)
        real(dp) :: da(curve%degree + 1), db(curve%degree + 1), xa(2), xb(2)
        integer :: k1, a, c, ka, kb

        k1 = curve%degree + 1
        call basis(curve, xi, na, dna)
        call basis(curve, eta, nb, dnb)
        ell = 0.0_dp
        select case (mode)
        case (1)
            da = divided_difference(curve, xi, eta)
            db = 0.0_dp
            ell = log(abs(oa))
        case (2)
            da = offset_from(curve, 1.0_dp, -oa)
            db = offset_from(curve, 0.0_dp, ob)
        case (3)
            da = offset_from(curve, 0.0_dp, oa)
            db = offset_from(curve, 1.0_dp, -ob)
        case default
            da = na
            db = nb
        end select
        jac = 0.0_dp
        u = 0.0_dp
        do a = 1, k1
            xa = curve%nodes(:, gs_sl_node(curve, pa, a))
            xb = curve%nodes(:, gs_sl_node(curve, pb, a))
            do c = 1, 2
                ka = 2*(a - 1) + c
                kb = 2*k1 + ka
                glob(ka) = 2*(gs_sl_node(curve, pa, a) - 1) + c
                glob(kb) = 2*(gs_sl_node(curve, pb, a) - 1) + c
                jac(c, ka) = na(a)
                jac(2 + c, kb) = nb(a)
                jac(4 + c, ka) = dna(a)
                jac(6 + c, kb) = dnb(a)
                jac(8 + c, ka) = da(a)
                jac(8 + c, kb) = -db(a)
                u(c) = u(c) + na(a)*xa(c)
                u(2 + c) = u(2 + c) + nb(a)*xb(c)
                u(4 + c) = u(4 + c) + dna(a)*xa(c)
                u(6 + c) = u(6 + c) + dnb(a)*xb(c)
                u(8 + c) = u(8 + c) + da(a)*xa(c) - db(a)*xb(c)
            end do
        end do
    end subroutine pair_inputs

    pure subroutine pair_jet(u, ell, f, gr, he)
        !! Integrand G |t_x| |t_y| with gradient and Hessian in u; the
        !! log-split form near coincidence, the direct form elsewhere.
        real(dp), intent(in) :: u(n_in), ell
        real(dp), intent(out) :: f, gr(n_in), he(n_in, n_in)
        real(dp) :: out(66), d2, s2

        d2 = exp(2.0_dp*ell)*(u(9)**2 + u(10)**2)
        s2 = (u(1) + u(3))**2 + (u(2) - u(4))**2
        if (d2 < split_q_max*s2) then
            call generated_gs_single_layer_near_jet(u(1), u(2), u(3), u(4), &
                u(5), u(6), u(7), u(8), u(9), u(10), ell, out(1), out(2), &
                out(3), out(4), out(5), out(6), out(7), out(8), out(9), &
                out(10), out(11), out(12), out(13), out(14), out(15), out(16), &
                out(17), out(18), out(19), out(20), out(21), out(22), out(23), &
                out(24), out(25), out(26), out(27), out(28), out(29), out(30), &
                out(31), out(32), out(33), out(34), out(35), out(36), out(37), &
                out(38), out(39), out(40), out(41), out(42), out(43), out(44), &
                out(45), out(46), out(47), out(48), out(49), out(50), out(51), &
                out(52), out(53), out(54), out(55), out(56), out(57), out(58), &
                out(59), out(60), out(61), out(62), out(63), out(64), out(65), &
                out(66))
            f = out(1)
            gr = out(2:11)
            call unpack_hessian(out(12:66), he)
        else
            call generated_gs_single_layer_far_jet(u(1), u(2), u(3), u(4), &
                u(5), u(6), u(7), u(8), out(1), out(2), out(3), out(4), &
                out(5), out(6), out(7), out(8), out(9), out(10), out(11), &
                out(12), out(13), out(14), out(15), out(16), out(17), out(18), &
                out(19), out(20), out(21), out(22), out(23), out(24), out(25), &
                out(26), out(27), out(28), out(29), out(30), out(31), out(32), &
                out(33), out(34), out(35), out(36), out(37), out(38), out(39), &
                out(40), out(41), out(42), out(43), out(44), out(45))
            f = out(1)
            gr = 0.0_dp
            gr(1:8) = out(2:9)
            he = 0.0_dp
            call unpack_hessian(out(10:45), he(1:8, 1:8))
        end if
    end subroutine pair_jet

    pure subroutine unpack_hessian(packed, h)
        !! Row-wise upper triangle to a full symmetric matrix.
        real(dp), intent(in) :: packed(:)
        real(dp), intent(out) :: h(:, :)
        integer :: i, j, k

        k = 0
        do i = 1, size(h, 1)
            do j = i, size(h, 1)
                k = k + 1
                h(i, j) = packed(k)
                h(j, i) = packed(k)
            end do
        end do
    end subroutine unpack_hessian

    ! ------------------------------------------------------------------
    ! Panel basis
    ! ------------------------------------------------------------------

    pure subroutine basis(curve, xi, n, dn)
        type(gs_sl_curve_t), intent(in) :: curve
        real(dp), intent(in) :: xi
        real(dp), intent(out) :: n(:), dn(:)
        integer :: a, j

        n = 0.0_dp
        dn = 0.0_dp
        do a = 1, curve%degree + 1
            do j = curve%degree + 1, 1, -1
                n(a) = n(a)*xi + curve%lagrange(a, j)
                if (j > 1) dn(a) = dn(a)*xi + (j - 1)*curve%lagrange(a, j)
            end do
        end do
    end subroutine basis

    pure function divided_difference(curve, xi, eta) result(d)
        !! (N_a(xi) - N_a(eta))/(xi - eta) without cancellation:
        !! (xi**j - eta**j)/(xi - eta) = sum_l xi**l eta**(j-1-l).
        type(gs_sl_curve_t), intent(in) :: curve
        real(dp), intent(in) :: xi, eta
        real(dp) :: d(curve%degree + 1), h
        integer :: a, j, l

        d = 0.0_dp
        do a = 1, curve%degree + 1
            do j = 2, curve%degree + 1
                h = 0.0_dp
                do l = 0, j - 2
                    h = h + xi**l*eta**(j - 2 - l)
                end do
                d(a) = d(a) + curve%lagrange(a, j)*h
            end do
        end do
    end function divided_difference

    pure function offset_from(curve, base, h) result(d)
        !! N_a(base + h) - N_a(base) as a polynomial in h without constant
        !! term, so it keeps full relative accuracy as h -> 0.
        type(gs_sl_curve_t), intent(in) :: curve
        real(dp), intent(in) :: base, h
        real(dp) :: d(curve%degree + 1), coef(curve%degree + 1)
        integer :: a, j, l

        do a = 1, curve%degree + 1
            ! Taylor coefficients of N_a about base: coef(l+1) = N^(l)/l!.
            coef = 0.0_dp
            do j = 1, curve%degree + 1
                do l = 0, j - 1
                    coef(l + 1) = coef(l + 1) + curve%lagrange(a, j) &
                        *binomial(j - 1, l)*base**(j - 1 - l)
                end do
            end do
            d(a) = 0.0_dp
            do l = curve%degree, 1, -1
                d(a) = (d(a) + coef(l + 1))*h
            end do
        end do
    end function offset_from

    pure function binomial(n, k) result(c)
        integer, intent(in) :: n, k
        real(dp) :: c
        integer :: i

        c = 1.0_dp
        do i = 1, k
            c = c*real(n - k + i, dp)/real(i, dp)
        end do
    end function binomial

    function pair_distance(curve, pa, pb) result(dmin)
        type(gs_sl_curve_t), intent(in) :: curve
        integer, intent(in) :: pa, pb
        real(dp) :: dmin, xa(2), xb(2), t(2)
        integer :: i, j

        dmin = huge(1.0_dp)
        do i = 0, 4
            call gs_sl_point(curve, pa, 0.25_dp*i, xa, t)
            do j = 0, 4
                call gs_sl_point(curve, pb, 0.25_dp*j, xb, t)
                dmin = min(dmin, norm2(xa - xb))
            end do
        end do
    end function pair_distance

    pure subroutine invert_small(a, status)
        !! In-place Gauss-Jordan inverse with partial pivoting.
        real(dp), intent(inout) :: a(:, :)
        integer, intent(out) :: status
        real(dp) :: aug(size(a, 1), 2*size(a, 1)), row(2*size(a, 1)), piv
        integer :: n, i, k, p

        n = size(a, 1)
        aug = 0.0_dp
        aug(:, 1:n) = a
        do i = 1, n
            aug(i, n + i) = 1.0_dp
        end do
        status = 1
        do k = 1, n
            p = k - 1 + maxloc(abs(aug(k:n, k)), dim=1)
            if (aug(p, k) == 0.0_dp) return
            row = aug(k, :)
            aug(k, :) = aug(p, :)
            aug(p, :) = row
            piv = aug(k, k)
            aug(k, :) = aug(k, :)/piv
            do i = 1, n
                if (i /= k) aug(i, :) = aug(i, :) - aug(i, k)*aug(k, :)
            end do
        end do
        a = aug(:, n + 1:2*n)
        status = 0
    end subroutine invert_small

end module fortfem_gs_single_layer_shape_2d
