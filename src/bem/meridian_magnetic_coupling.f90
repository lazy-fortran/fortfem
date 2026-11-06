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
    !!
    !! Mode n = 0.  The exterior of a torus is not simply connected, so for
    !! n = 0 a single-valued scalar potential misses the poloidal field of
    !! a net toroidal current (C annihilates a constant a = R A_phi on
    !! Gamma).  The n = 0 poloidal field is closed instead by the
    !! Grad-Shafranov Steklov-Poincare map S_GS of
    !! fortfem_axisymmetric_gs_boundary_operators_2d acting on the Dirichlet
    !! trace of a: with the trace pairing T(i, j) = int psi_i a_j ds and the
    !! P2 trace mass M the trace is u = M^{-1} T A, and the exterior energy
    !! int |grad a|^2 / R dR dZ is (T A)^H M^{-1} S_GS M^{-1} (T A).  The
    !! n = 0 toroidal field has no decaying single-valued exterior
    !! extension; it vanishes outside, which is the natural condition
    !! R B_phi = 0 of the interior form.
    !!
    !! Coupled solve.  meridian_exterior_t holds the trace map T (T = C for
    !! n /= 0) and the dense exterior form E (S_n^{-1}, or M^{-1} S_GS M^{-1}
    !! for n = 0).  meridian_exterior_solve solves the Coulomb-gauge mixed
    !! problem
    !!
    !!     (curl_n A', curl_n A)_M2 + (T A')^H E (T A) + (A', grad_n p)_M1
    !!         = f(A'),        (A, grad_n q)_M1 = 0
    !!
    !! for all A' in V1 and q in V0.  The multiplier p vanishes when f
    !! annihilates gradients (a divergence-free source tangential to the
    !! boundary).  For n = 0 the constants leave the multiplier space
    !! (grad_0 1 = 0).
    use fortfem_kinds, only: dp
    use fortfem_meridian_feec_n, only: meridian_feec_t, &
        meridian_point_basis_t, meridian_point_basis_init, &
        meridian_point_basis, meridian_op_complex, meridian_zcsc_matmul
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gs_bem_curve_t, &
        gs_bem_curve_from_nodes, p2_basis, gs_bem_operators_t, &
        gs_bem_assemble, gs_bem_steklov_symmetric, gs_bem_panel_point
    use fortfem_toroidal_mode_boundary_operators_2d, only: tm_bem_assemble, &
        tm_bem_steklov_symmetric, tm_bem_neumann_to_dirichlet
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortfem_sparse_direct, only: sparse_direct_solve_once
    use fortsparse, only: csc_z_t, csc_from_triplet, fortsparse_status_t, &
        status_ok
    implicit none
    private

    !> Counter-clockwise boundary chain: panel p is edge edge(p) of triangle
    !> tri(p) (local edge i joins local vertices i and i mod 3 + 1).
    type, public :: meridian_boundary_t
        integer :: n_panels = 0
        integer, allocatable :: tri(:), edge(:)
        type(gs_bem_curve_t) :: curve
    end type meridian_boundary_t

    !> Exterior closure of one toroidal mode n on a boundary curve: the
    !> exterior field energy of the interior potential A is (T A)^H E (T A).
    type, public :: meridian_exterior_t
        integer :: n = 0
        type(gs_bem_curve_t) :: curve
        type(csc_z_t) :: trace !! T (2 n_panels x dim V1)
        real(dp), allocatable :: form(:, :) !! E (2 n_panels squared)
    end type meridian_exterior_t

    public :: meridian_boundary_build
    public :: meridian_edge_point
    public :: meridian_flux_pairing
    public :: meridian_exterior_form
    public :: meridian_exterior_build
    public :: meridian_exterior_solve
    public :: meridian_exterior_energy

    interface
        subroutine dposv(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: info
        end subroutine dposv
    end interface

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

    subroutine meridian_exterior_form(curve, n, form, status)
        !! Dense exterior form E on the P2 traces of the curve: S_n^{-1} for
        !! n /= 0 (acting on the normal-flux pairing), M^{-1} S_GS M^{-1} for
        !! n = 0 (acting on the Dirichlet-trace pairing of a = R A_phi).
        type(gs_bem_curve_t), intent(in) :: curve
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: form(:, :)
        integer, intent(out) :: status

        type(gs_bem_operators_t) :: ops
        real(dp), allocatable :: s(:, :), m(:, :), y(:, :)
        integer :: nx, info

        if (n /= 0) then
            call tm_bem_assemble(curve, n, 1, ops, status)
            if (status /= 0) return
            call tm_bem_steklov_symmetric(ops, s, status)
            if (status /= 0) return
            call tm_bem_neumann_to_dirichlet(s, form, status)
            return
        end if
        call gs_bem_assemble(curve, 1, ops, status)
        if (status /= 0) return
        call gs_bem_steklov_symmetric(ops, s, status)
        if (status /= 0) return
        status = 1
        nx = size(s, 1)
        m = ops%mxx
        y = s
        call dposv('L', nx, nx, m, nx, y, nx, info)
        if (info /= 0) return
        form = transpose(y)
        m = ops%mxx
        call dposv('L', nx, nx, m, nx, form, nx, info)
        if (info /= 0) return
        form = 0.5_dp*(form + transpose(form))
        status = 0
    end subroutine meridian_exterior_form

    subroutine meridian_exterior_build(fs, bd, n, ex, status)
        !! Exterior closure of mode n on the boundary chain of the curved
        !! Nedelec complex: T = C2 curl_n for n /= 0, the Dirichlet-trace
        !! pairing of a = R A_phi for n = 0.
        type(meridian_feec_t), intent(in) :: fs
        type(meridian_boundary_t), intent(in) :: bd
        integer, intent(in) :: n
        type(meridian_exterior_t), intent(out) :: ex
        integer, intent(out) :: status

        type(meridian_point_basis_t) :: pb
        type(csc_z_t) :: c2, curl
        type(fortsparse_status_t) :: fst
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:), tau(:), w(:), phi(:, :)
        real(dp) :: xi, eta, nh(2), x(2), xc(2), tg(2), speed, psi(3), dpsi(3)
        integer :: nb, np, nln, nl0, p, q, j, a, t, pos, st, nq, idx

        status = 1
        ex%n = n
        ex%curve = bd%curve
        np = bd%n_panels
        nb = 2*np
        if (n /= 0) then
            call meridian_flux_pairing(fs, bd, rows, cols, vals, st)
            if (st /= 0) return
            call csc_from_triplet(nb, fs%sp(2)%ndof, rows, cols, &
                cmplx(vals, 0.0_dp, dp), c2, fst)
            if (.not. status_ok(fst)) return
            call meridian_op_complex(fs%curl, n, curl, st)
            if (st /= 0) return
            call meridian_zcsc_matmul(c2, curl, ex%trace, st)
            if (st /= 0) return
        else
            nln = fs%nln
            nl0 = fs%nl0
            nq = fs%k0 + 3
            allocate (tau(nq), w(nq), phi(3, fs%sp(1)%nloc))
            allocate (rows(3*nl0*np), cols(3*nl0*np), vals(3*nl0*np))
            call gauss_legendre_ab(nq, 0.0_dp, 1.0_dp, tau, w)
            call meridian_point_basis_init(fs, pb, st)
            if (st /= 0) return
            vals = 0.0_dp
            pos = 0
            do p = 1, np
                t = bd%tri(p)
                do j = 1, nl0
                    do a = 1, 3
                        idx = 2*p - 2 + a
                        if (idx > nb) idx = 1
                        rows(pos + 3*(j - 1) + a) = idx
                        cols(pos + 3*(j - 1) + a) = fs%sp(1)%dof(nln + j, t)
                    end do
                end do
                do q = 1, nq
                    call meridian_edge_point(bd%edge(p), tau(q), xi, eta, nh)
                    call meridian_point_basis(fs, pb, 1, t, xi, eta, x, phi, &
                        st)
                    if (st /= 0) return
                    call gs_bem_panel_point(bd%curve, p, tau(q), xc, tg, speed)
                    call p2_basis(tau(q), psi, dpsi)
                    do j = 1, nl0
                        vals(pos + 3*(j - 1) + 1:pos + 3*j) = &
                            vals(pos + 3*(j - 1) + 1:pos + 3*j) &
                            + w(q)*speed*psi*phi(2, nln + j)*x(1) &
                            *real(fs%sp(1)%sgn(nln + j, t), dp)
                    end do
                end do
                pos = pos + 3*nl0
            end do
            call csc_from_triplet(nb, fs%sp(1)%ndof, rows, cols, &
                cmplx(vals, 0.0_dp, dp), ex%trace, fst)
            if (.not. status_ok(fst)) return
        end if
        call meridian_exterior_form(ex%curve, n, ex%form, status)
    end subroutine meridian_exterior_build

    subroutine meridian_exterior_solve(k, gauge, ex, f, u, p, status)
        !! Coulomb-gauge mixed FEM-BEM solve (see the module header).  k is
        !! the curl-curl form curl_n^H M2 curl_n, gauge = M1 grad_n and f the
        !! load f_i = f(A'_i).  u are the V1 coefficients, p the multiplier.
        type(csc_z_t), intent(in) :: k, gauge
        type(meridian_exterior_t), intent(in) :: ex
        complex(dp), intent(in) :: f(:)
        complex(dp), allocatable, intent(out) :: u(:), p(:)
        integer, intent(out) :: status

        type(csc_z_t) :: a
        type(fortsparse_status_t) :: fst
        integer, allocatable :: rows(:), cols(:), list(:), slot(:)
        complex(dp), allocatable :: vals(:), cd(:, :), d(:, :)
        complex(dp), allocatable :: rhs(:), x(:)
        integer :: n1, n0, nb, nj, j, e, i, pos, total
        logical :: pin

        status = 1
        n1 = k%nrow
        n0 = gauge%ncol
        nb = ex%trace%nrow
        if (k%ncol /= n1 .or. gauge%nrow /= n1) return
        if (ex%trace%ncol /= n1 .or. size(f) /= n1) return
        if (any(shape(ex%form) /= [nb, nb])) return
        ! n = 0: grad_0 1 = 0, so one multiplier is fixed to zero
        pin = ex%n == 0
        allocate (slot(n1))
        slot = 0
        nj = 0
        do j = 1, n1
            if (ex%trace%col_ptr(j + 1) > ex%trace%col_ptr(j)) then
                nj = nj + 1
                slot(j) = nj
            end if
        end do
        allocate (list(nj), cd(nb, nj))
        cd = (0.0_dp, 0.0_dp)
        do j = 1, n1
            if (slot(j) == 0) cycle
            list(slot(j)) = j
            do e = ex%trace%col_ptr(j), ex%trace%col_ptr(j + 1) - 1
                cd(ex%trace%row_idx(e), slot(j)) = &
                    cd(ex%trace%row_idx(e), slot(j)) + ex%trace%val(e)
            end do
        end do
        d = matmul(conjg(transpose(cd)), matmul(cmplx(ex%form, 0.0_dp, dp), &
            cd))
        total = k%nnz + nj*nj + 2*gauge%nnz + 1
        allocate (rows(total), cols(total), vals(total))
        pos = 0
        do j = 1, n1
            do e = k%col_ptr(j), k%col_ptr(j + 1) - 1
                pos = pos + 1
                rows(pos) = k%row_idx(e)
                cols(pos) = j
                vals(pos) = k%val(e)
            end do
        end do
        do j = 1, nj
            do i = 1, nj
                pos = pos + 1
                rows(pos) = list(i)
                cols(pos) = list(j)
                vals(pos) = d(i, j)
            end do
        end do
        do j = 1, n0
            if (pin .and. j == 1) cycle
            do e = gauge%col_ptr(j), gauge%col_ptr(j + 1) - 1
                pos = pos + 1
                rows(pos) = gauge%row_idx(e)
                cols(pos) = n1 + j
                vals(pos) = gauge%val(e)
                pos = pos + 1
                rows(pos) = n1 + j
                cols(pos) = gauge%row_idx(e)
                vals(pos) = conjg(gauge%val(e))
            end do
        end do
        if (pin) then
            pos = pos + 1
            rows(pos) = n1 + 1
            cols(pos) = n1 + 1
            vals(pos) = (1.0_dp, 0.0_dp)
        end if
        call csc_from_triplet(n1 + n0, n1 + n0, rows(:pos), cols(:pos), &
            vals(:pos), a, fst)
        if (.not. status_ok(fst)) return
        allocate (rhs(n1 + n0), x(n1 + n0))
        rhs = (0.0_dp, 0.0_dp)
        rhs(:n1) = f
        call sparse_direct_solve_once(a, rhs, x, fst)
        if (.not. status_ok(fst)) return
        u = x(:n1)
        p = x(n1 + 1:)
        status = 0
    end subroutine meridian_exterior_solve

    function meridian_exterior_energy(ex, u) result(energy)
        !! Exterior field energy (T u)^H E (T u), with the 2 pi of the phi
        !! integral and the 1/(2 mu0) dropped as in the FEEC masses.
        type(meridian_exterior_t), intent(in) :: ex
        complex(dp), intent(in) :: u(:)
        real(dp) :: energy

        complex(dp), allocatable :: tu(:)
        integer :: j, e

        allocate (tu(ex%trace%nrow))
        tu = (0.0_dp, 0.0_dp)
        do j = 1, ex%trace%ncol
            do e = ex%trace%col_ptr(j), ex%trace%col_ptr(j + 1) - 1
                tu(ex%trace%row_idx(e)) = tu(ex%trace%row_idx(e)) &
                    + ex%trace%val(e)*u(j)
            end do
        end do
        energy = real(dot_product(tu, matmul(cmplx(ex%form, 0.0_dp, dp), &
            tu)), dp)
    end function meridian_exterior_energy

end module fortfem_meridian_magnetic_coupling
