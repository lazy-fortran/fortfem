module fortfem_meridian_magnetic_coupling
    !! Coupling of the per-n FEEC complex of fortfem_meridian_feec_n to the
    !! toroidal-mode boundary operators of
    !! fortfem_toroidal_mode_boundary_operators_2d on the boundary of the
    !! meridional mesh (the magnetic FEM-BEM transmission conditions).
    !!
    !! Boundary curve.  The boundary edges, ordered counter-clockwise, are
    !! the panels of a gs_bem_curve_t: panel p is boundary edge p, its nodes
    !! 2p - 1 and 2p are the start vertex and the geometric mid node of the
    !! edge.  For geometry order 1 or 2 the panel polynomial is the edge
    !! geometry itself and the panel parameter xi in [0, 1] is the reference
    !! edge parameter, so traces live on one curve.  Geometry order 3 has no
    !! P2 panel and is rejected.
    !!
    !! Normal-flux pairing.  For the V2 field B = (q, s), q = R B_pol,
    !!
    !!     C2(i, j) = int_Gamma psi_i (q_j . n) ds,
    !!
    !! with the P2 trace basis psi_i of the curve and n the outward normal.
    !! The contravariant Piola map gives q . n ds = qhat . nhat dshat, so C2
    !! is assembled from reference RT values and is exact for curved edges.
    !! (q . n) = R B . n, hence int psi (q . n) ds = int psi B . n dS / (2 pi)
    !! is the per-n surface pairing (dS = R ds dphi, the 2 pi dropped as in
    !! the FEEC masses).
    !!
    !! Exterior closure.  Outside Gamma a current-free per-n field is
    !! B = grad_n u with L_n u = 0, u decaying.  With t = R d_n u the
    !! Steklov-Poincare map of the boundary-operator module gives
    !! S u = -t, so the transmission conditions n . curl_n A = d_n u and
    !! n x H continuous give, for the interior vector potential A in V1,
    !!
    !!     (curl_n A', curl_n A)_M2 + <C A', S^{-1} C A> = mu0 (A', J),
    !!
    !! with C = C2 curl_n: the exterior energy enters as the Hermitian
    !! positive Schur complement C^H S^{-1} C on the boundary trace.
    use fortfem_kinds, only: dp
    use fortfem_meridian_feec_n, only: meridian_feec_t, &
        meridian_point_basis_t, meridian_point_basis_init, meridian_point_basis
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gs_bem_curve_t, &
        gs_bem_curve_from_nodes, p2_basis
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private

    !> Counter-clockwise boundary chain: panel p is edge edge(p) of triangle
    !> tri(p) (local edge i joins local vertices i and i mod 3 + 1).
    type, public :: meridian_boundary_t
        integer :: n_panels = 0
        integer, allocatable :: tri(:), edge(:)
        type(gs_bem_curve_t) :: curve
    end type meridian_boundary_t

    public :: meridian_boundary_build
    public :: meridian_edge_point
    public :: meridian_flux_pairing

contains

    pure subroutine meridian_edge_point(edge, tau, xi, eta, nhat)
        !! Reference point (xi, eta) at parameter tau of local edge edge
        !! (vertex edge -> vertex edge mod 3 + 1, vertices (0,0), (1,0),
        !! (0,1)) and nhat = rot(-90) dxhat/dtau, the outward normal times
        !! dshat/dtau of a counter-clockwise triangle.
        integer, intent(in) :: edge
        real(dp), intent(in) :: tau
        real(dp), intent(out) :: xi, eta, nhat(2)

        select case (edge)
        case (1)
            xi = tau
            eta = 0.0_dp
            nhat = [0.0_dp, -1.0_dp]
        case (2)
            xi = 1.0_dp - tau
            eta = tau
            nhat = [1.0_dp, 1.0_dp]
        case default
            xi = 0.0_dp
            eta = 1.0_dp - tau
            nhat = [-1.0_dp, 0.0_dp]
        end select
    end subroutine meridian_edge_point

    subroutine meridian_boundary_build(fs, bd, status)
        !! Ordered boundary chain and its P2 panel curve.  The boundary must
        !! be one closed curve (the domain simply connected in the plane).
        type(meridian_feec_t), intent(in) :: fs
        type(meridian_boundary_t), intent(out) :: bd
        integer, intent(out) :: status

        type(meridian_point_basis_t) :: pb
        integer, allocatable :: next_t(:), next_e(:), start_v(:)
        real(dp), allocatable :: nodes(:, :), phi(:, :)
        real(dp) :: xi, eta, nh(2), x(2)
        integer :: nt, nv, t, i, v, nb, p, k, c, st

        status = 1
        if (fs%g > 2) return
        nt = fs%mm%topology%n_triangles
        nv = fs%mm%topology%n_vertices
        k = fs%k
        allocate (next_t(nv), next_e(nv))
        next_t = 0
        next_e = 0
        nb = 0
        do t = 1, nt
            do i = 1, 3
                if (.not. fs%sp(1)%bnd(fs%sp(1)%dof((i - 1)*k + 1, t))) cycle
                v = fs%mm%topology%triangles(i, t)
                if (next_t(v) /= 0) return
                next_t(v) = t
                next_e(v) = i
                nb = nb + 1
            end do
        end do
        if (nb < 4) return
        allocate (bd%tri(nb), bd%edge(nb), start_v(nb))
        v = 0
        do i = 1, nv
            if (next_t(i) /= 0) then
                v = i
                exit
            end if
        end do
        do p = 1, nb
            if (v == 0) return
            if (next_t(v) == 0) return
            start_v(p) = v
            bd%tri(p) = next_t(v)
            bd%edge(p) = next_e(v)
            c = mod(bd%edge(p), 3) + 1
            v = fs%mm%topology%triangles(c, bd%tri(p))
            if (p < nb .and. v == start_v(1)) return
        end do
        if (v /= start_v(1)) return
        bd%n_panels = nb
        call meridian_point_basis_init(fs, pb, st)
        if (st /= 0) return
        allocate (nodes(2, 2*nb), phi(3, fs%sp(0)%nloc))
        do p = 1, nb
            nodes(:, 2*p - 1) = fs%mm%topology%vertices(:, start_v(p))
            call meridian_edge_point(bd%edge(p), 0.5_dp, xi, eta, nh)
            call meridian_point_basis(fs, pb, 0, bd%tri(p), xi, eta, x, phi, st)
            if (st /= 0) return
            nodes(:, 2*p) = x
        end do
        call gs_bem_curve_from_nodes(nodes, bd%curve, status)
    end subroutine meridian_boundary_build

    subroutine meridian_flux_pairing(fs, bd, rows, cols, vals, status)
        !! Triplets of C2(i, j) = int psi_i (q_j . n) ds for the P2 trace
        !! dofs i = 1 .. 2 np of the boundary curve and the V2 dofs j.
        type(meridian_feec_t), intent(in) :: fs
        type(meridian_boundary_t), intent(in) :: bd
        integer, allocatable, intent(out) :: rows(:), cols(:)
        real(dp), allocatable, intent(out) :: vals(:)
        integer, intent(out) :: status

        type(meridian_point_basis_t) :: pb
        real(dp), allocatable :: tau(:), w(:), phi(:, :), ref(:, :), loc(:, :)
        real(dp) :: xi, eta, nh(2), x(2), psi(3), dpsi(3), flux
        integer :: nq, p, q, j, a, t, nln, np, pos, idx, st

        status = 1
        np = bd%n_panels
        nln = fs%nln
        nq = fs%k + 3
        allocate (tau(nq), w(nq), phi(3, fs%sp(2)%nloc), ref(2, nln))
        allocate (loc(3, nln))
        allocate (rows(3*nln*np), cols(3*nln*np), vals(3*nln*np))
        call gauss_legendre_ab(nq, 0.0_dp, 1.0_dp, tau, w)
        call meridian_point_basis_init(fs, pb, st)
        if (st /= 0) return
        pos = 0
        do p = 1, np
            t = bd%tri(p)
            loc = 0.0_dp
            do q = 1, nq
                call meridian_edge_point(bd%edge(p), tau(q), xi, eta, nh)
                call meridian_point_basis(fs, pb, 2, t, xi, eta, x, phi, st, &
                    ref)
                if (st /= 0) return
                call p2_basis(tau(q), psi, dpsi)
                do j = 1, nln
                    flux = ref(1, j)*nh(1) + ref(2, j)*nh(2)
                    loc(:, j) = loc(:, j) + w(q)*psi*flux
                end do
            end do
            do j = 1, nln
                do a = 1, 3
                    pos = pos + 1
                    idx = 2*p - 2 + a
                    if (idx > 2*np) idx = 1
                    rows(pos) = idx
                    cols(pos) = fs%sp(2)%dof(j, t)
                    vals(pos) = real(fs%sp(2)%sgn(j, t), dp)*loc(a, j)
                    if (abs(vals(pos)) < 1.0e-13_dp) vals(pos) = 0.0_dp
                end do
            end do
        end do
        status = 0
    end subroutine meridian_flux_pairing

end module fortfem_meridian_magnetic_coupling
