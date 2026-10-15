module fortfem_equilibrated_hcurl_meridian_n
    !! Patchwise equilibrated (Prager--Synge) energy bound for the per-mode
    !! meridian curl-curl problem of fortfem_meridian_feec_n,
    !!
    !!   (curl_n v, curl_n A)_M2 = (v, J)_M1  for all v = (v_pol, 0) with
    !!   v x n = 0,  A = (A_pol, 0) (gauge a = R A_phi = 0),  n /= 0,
    !!
    !! in physical components with the measure R dR dZ, J the already
    !! mu0-scaled current (only its poloidal part enters the gauge-fixed
    !! problem).  A_h is any V1 field with a = 0; the bound is for
    !! ||curl_n (A - A_h)||_M2 and holds for every A_h, Galerkin or not.
    !!
    !! Proof.  For v = (v_pol, 0), curl_n v = (-i n rot90 v_pol / R, -rot v_pol)
    !! (physical components), so for H = (H_pol, s)
    !!
    !!   (H, curl_n v) = int R conj(H_pol).(-i n rot90 v_pol / R)
    !!                   - R conj(s) rot v_pol .
    !!
    !! The H_pol part contains no derivative of v, so H_pol needs only L2; the
    !! s part is integrated by parts, which needs s in H^1 (continuous) and
    !! produces no boundary term because v x n = 0.  Hence for every
    !! continuous scalar s, the field
    !!
    !!   H(s) = ( (1/(i n)) (grad(R s) + R rot90 J), s ),   rot90 v = (-v_Z, v_R),
    !!
    !! satisfies the strong Ampere law curl_n^* H = J tested against the whole
    !! gauge-fixed space W0 = {(v_pol, 0), v x n = 0}: (curl_n v, H) = (v, J)
    !! for all v in W0.  With B = curl_n A the exact field this gives (curl_n e,
    !! H - B) = 0 for e = A - A_h in W0, so
    !!
    !!   ||curl_n e||^2 + ||H - B||^2 = ||H - B_h||^2,  B_h = curl_n A_h,
    !!
    !!   ||curl_n (A - A_h)||_M2 <= ||H(s) - curl_n A_h||_M2  for every
    !!   continuous s.
    !!
    !! No divergence condition on J and no Friedrichs-type constant enters
    !! (the gauge a = 0 removes ker curl_n; the zeroth-order term n^2/R of the
    !! form is what makes the dual field algebraic in the poloidal part).  The
    !! bound is thus guaranteed up to quadrature of the integrals and the
    !! exactness of the discrete problem (the exact A solves the same W0
    !! problem with the same J).
    !!
    !! Choice of s.  H(s) satisfies the constraint for every continuous s, so s
    !! only has to make ||H(s) - B_h|| small: s minimizes this quadratic
    !! functional over the continuous Lagrange space of the potential complex
    !! fsd (V0 of fsd, order m >= 1; efficiency ~ 1 needs m >= k + 1, the poloidal
    !! part of the residual contains grad(R s)).  Two constructions:
    !!
    !!  HCURL_GLOBAL (default): the global minimizer; Hermitian positive
    !!    definite sparse normal equations (the mass of s is part of the Gram
    !!    matrix) solved by sparse LU.  The minimizer is nonlocal at the scale
    !!    1/|n| (screened Poisson problem with screening length 1/|n|), which
    !!    is why a purely local construction cannot be efficient for n h << 1.
    !!    Measured efficiency index 1.000-1.01 (m = k + 1, fine meshes).
    !!
    !!  HCURL_PATCHWISE: s = sum_a s_a, s_a in the Lagrange space on the vertex
    !!    patch omega_a, zero at the nodes of the edge opposite to a (so s is
    !!    continuous; free on the domain boundary), minimizing
    !!      ||h_a - psi_a Bh_pol||^2 + ||s_a - psi_a sh||^2,
    !!      h_a = (1/(i n)) (grad(R s_a) + psi_a R rot90 J - R sh grad psi_a),
    !!    psi_a the P1 hat function and (Bh_pol, sh) = B_h.  sum psi_a = 1 and
    !!    sum grad psi_a = 0 make sum_a (h_a, s_a) = H(sum_a s_a) (the
    !!    commutator term R sh grad psi_a cancels in the sum).  O(N) work, one
    !!    small Cholesky solve per patch, OpenMP over patches; exactly
    !!    equilibrated B_h with s in P_{m-1} is reproduced to rounding.  The
    !!    bound is still guaranteed, but the commutator data carry the jump
    !!    error of sh (O(h^k) divided by the patch size), so the efficiency
    !!    index degrades like 1/(n h) when n h << 1 (measured: 6 -> 22 over
    !!    four refinements for n = 1, k = 2; 1.03 -> 1.4 for n = 10).  It is
    !!    useful for large n or as a cheap check, not as the default.
    !!
    !! In both cases the bound is evaluated on the assembled s, so solver
    !! errors only loosen it.  The Pythagoras identity gives
    !! bound^2 = ||curl_n e||^2 + ||H(s) - curl_n A||^2.  (A weaker variant
    !! with an approximate constraint would add (1/|n|) (int R^3 |J -
    !! curl_n^* H|^2)^(1/2), since ||e/sqrt(R)|| <= ||curl_n e||/|n| for a = 0;
    !! it is not needed because H(s) is exactly equilibrated.)
    !!
    !! Edge (sheet) loads and the free trace.  Let the load also contain
    !! L(v) = int_E G (v . t) dl over oriented mesh edges E (G per unit
    !! length, t the unit tangent of the edge orientation), e.g. a toroidal
    !! sheet current on an interior curve, R K_t = G.  For real v and a
    !! piecewise continuous potential, integration by parts cell by cell gives
    !!
    !!   (curl_n v, H(s)) = (v, J) - sum_K oint_dK (R s) v . t_K dl,
    !!
    !! t_K counter-clockwise for K, so the continuous construction carries
    !! over when sigma = R s jumps across the loaded interior edges by
    !! [sigma] = sigma_left - sigma_right = -G (left: the cell for which the
    !! edge runs counter-clockwise).  Without the trace condition v x n = 0
    !! (free_trace: a coupled exterior such as the Steklov-Poincare closure
    !! (curl_n v, curl_n A) + <C v, S^{-1} C A> = (v, J) + L(v), C = C2 curl_n)
    !! the boundary term must be cancelled by the exterior one: with the
    !! boundary loads G_b defined by int G_b (v . t) dl = -<C v, mu>, the
    !! prescribed trace sigma = -G_b on the boundary gives
    !!
    !!   |||A - A_h|||^2 <= ||H(sigma) - curl_n A_h||^2 + ||mu - S^{-1} C A_h||_S^2,
    !!
    !! |||v|||^2 = ||curl_n v||^2 + <C v, S^{-1} C v>; with mu = S^{-1} C A_h
    !! the exterior term vanishes and the bound refers to the problem with
    !! that (discrete) exterior operator.  The caller supplies G_b and owns
    !! the exterior term; any remainder of an approximated exterior operator
    !! is outside this bound.
    !!
    !! With edge loads or free_trace the potential is sigma = R s in V0 of fsd
    !! (sigma = sigma_0 + w, sigma_0 continuous and zero on the boundary with
    !! free_trace, w the fixed broken lift: the interpolant of the jump data
    !! on the side cells of the interior loaded edges plus the interpolant of
    !! the boundary data).  The jump and trace are then those of the
    !! interpolant I G of G at the V0 nodes of the loaded edges, so the bound
    !! is for the problem with the loads I G; est%load_defect reports the
    !! relative L2 (edge parameter) size of G - I G as an indicator (its
    !! H^(-1/2) effect is not bounded).  Only HCURL_GLOBAL supports loads.
    !!
    !! Requirements: fs and fsd are built on the same mesh (equal geometry
    !! order) with the same quadrature degree; the source callback is called
    !! from OpenMP threads and must be thread-safe.  The integrals use that
    !! rule; their quadrature error is not enclosed.
    use fortfem_kinds, only: dp
    use fortfem_generated_meridian_hcurl_equilibration, only: &
        generated_meridian_hcurl_residual
    use fortfem_generated_meridian_hcurl_scaled_potential, only: &
        generated_meridian_hcurl_scaled_potential
    use fortfem_generated_meridian_edge_coordinate, only: &
        generated_meridian_edge_coordinate
    use fortfem_generated_meridian_edge_lift, only: generated_meridian_edge_lift
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_field_i, &
        meridian_op_complex, meridian_qgrad, meridian_qpoint
    use fortfem_equilibrated_flux_topology, only: vertex_stars
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        initialize_triangle_lagrange_basis, triangle_lagrange_basis_t, &
        triangle_lagrange_nodes, evaluate_triangle_lagrange_basis
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_polynomial, only: lagrange_weights, lagrange_deriv_weights
    use fortsparse, only: csc_from_triplet, csc_matvec, csc_z_t, &
        fortsparse_status_t, sparse_solve_once, status_ok
    implicit none
    private

    public :: meridian_hcurl_estimate_t, meridian_hcurl_equilibrated_bound
    public :: meridian_edge_load_i
    public :: meridian_hcurl_edge_remainder

    integer, parameter, public :: HCURL_OK = 0
    integer, parameter, public :: HCURL_BAD_INPUT = 1
    integer, parameter, public :: HCURL_NOT_GAUGE_FIXED = 2
    integer, parameter, public :: HCURL_BAD_CELL = 3
    integer, parameter, public :: HCURL_SOLVE_FAILED = 4

    !> Choice of the potential: global least squares or patchwise sum.
    integer, parameter, public :: HCURL_GLOBAL = 1
    integer, parameter, public :: HCURL_PATCHWISE = 2

    type :: meridian_hcurl_estimate_t
        real(dp) :: bound = huge(1.0_dp) !! ||H(s) - curl_n A_h||_M2
        real(dp) :: field_norm = 0.0_dp !! ||curl_n A_h||_M2
        !> Per triangle: (int_T |H(s) - curl_n A_h|^2 R dR dZ)^(1/2).
        real(dp), allocatable :: cell_indicator(:)
        !> Coefficients of s = H_phi in V0 of fsd (with edge loads or a free
        !> trace: of the continuous part sigma_0 of sigma = R s).
        complex(dp), allocatable :: potential(:)
        !> Potential sigma = R s with a broken lift (edge loads, free trace).
        logical :: scaled_potential = .false.
        !> Boundary trace of sigma prescribed (free_trace).
        logical :: trace_fixed = .false.
        !> Relative L2 (edge parameter) size of G - I G over the loaded edges.
        real(dp) :: load_defect = 0.0_dp
        !> Largest relative disagreement of the load data at shared nodes.
        real(dp) :: lift_mismatch = 0.0_dp
    end type meridian_hcurl_estimate_t

    abstract interface
        subroutine meridian_edge_load_i(cell, edge, tau, g)
            !! Load density G at parameter tau in [0, 1] of local edge edge of
            !! triangle cell (from local vertex edge to vertex mod(edge, 3) +
            !! 1), L(v) = int G (v . t) dl with t along that orientation.
            import :: dp
            integer, intent(in) :: cell, edge
            real(dp), intent(in) :: tau
            complex(dp), intent(out) :: g
        end subroutine meridian_edge_load_i
    end interface

    interface
        subroutine zposv(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: info
            complex(dp), intent(inout) :: a(lda, *), b(ldb, *)
        end subroutine zposv
    end interface

contains

    subroutine meridian_hcurl_edge_remainder(fsd, n, cells, edges, load, &
            load_derivative, bound, status)
        !! Energy majorant of the correction from I_k G to the actual G,
        !! using the same gauge and metric as meridian_hcurl_equilibrated_bound.
        !! load_derivative returns dG/dtau, NOT the arc-length derivative.
        !! Each physical loaded edge occurs once, oriented CCW in its listed
        !! cell. Callbacks must be thread-safe and G must be H^1 on each edge.
        !!
        !! Proof: r=G-I_k G vanishes at both edge endpoints. In the listed
        !! cell put sigma=-(la+lb)^4*r(lb/(la+lb)), zero in the other cells.
        !! This H^1 lift is zero on the two other edges and jumps by -r.
        !! Sum the lifts before evaluating H=(grad(sigma)/(i*n),sigma/R).
        !! Integration by parts gives (curl_n v,H)=L_r(v). Thus
        !! |||A_G-A_IG||| <= ||H|| with constant ONE, also with any positive
        !! exterior energy form. No trace or mesh-dependent constant enters.
        !!
        !! The returned norm uses sampled volume quadrature and ordinary
        !! arithmetic: practical, NOT an enclosed/verified numerical radius.
        !! Cost O(N + nq*number_of_loaded_cells) at fixed order, OpenMP.
        type(meridian_feec_t), intent(in) :: fsd
        integer, intent(in) :: n, cells(:), edges(:)
        procedure(meridian_edge_load_i) :: load, load_derivative
        real(dp), intent(out) :: bound
        integer, intent(out) :: status
        integer, allocatable :: lookup(:, :), active(:)
        complex(dp), allocatable :: coeff(:, :)
        real(dp), allocatable :: lamv(:, :), val(:), grad(:, :), nodes(:), iv(:), di(:)
        real(dp) :: x(2), wq, rr, grr(2), lam(3), glam(2, 3), tau
        real(dp) :: sr, si, gsr(2), gsi(2), total
        complex(dp) :: g, dg, gi, dgi, r, dr, sig, gs(2), h(3), zero2(2)
        integer :: nt, ne, nl, m, ie, t, e, a, b, j, q, st, bad, ic, nactive

        status = HCURL_BAD_INPUT
        bound = huge(1.0_dp)
        if (n == 0 .or. fsd%k < 1) return
        ne = size(cells)
        if (size(edges) /= ne) return
        nt = fsd%mm%topology%n_triangles
        if (any(cells < 1) .or. any(cells > nt)) return
        if (any(edges < 1) .or. any(edges > 3)) return
        allocate(lookup(3, nt), active(ne))
        lookup = 0
        nactive = 0
        do ie = 1, ne
            if (lookup(edges(ie), cells(ie)) /= 0) return
            if (all(lookup(:, cells(ie)) == 0)) then
                nactive = nactive + 1
                active(nactive) = cells(ie)
            end if
            lookup(edges(ie), cells(ie)) = ie
        end do
        call vertex_weights(fsd, lamv)
        nl = fsd%nl0
        m = fsd%k + 1
        allocate(coeff(m, ne), nodes(m))
        nodes = [(real(j - 1, dp)/fsd%k, j=1, m)]
        coeff = (0.0_dp, 0.0_dp)
        bad = 0
        !$omp parallel do default(shared) private(ie,t,e,j,g,tau) &
        !$omp reduction(max:bad) schedule(static)
        do ie = 1, ne
            t = cells(ie)
            e = edges(ie)
            do j = 1, m
                tau = nodes(j)
                call load(t, e, tau, g)
                if (.not. finite_complex(g)) bad = 1
                coeff(j, ie) = g
            end do
        end do
        !$omp end parallel do
        if (bad /= 0) return
        total = 0.0_dp
        zero2 = (0.0_dp, 0.0_dp)
        !$omp parallel default(shared) &
        !$omp private(ic,t,q,e,ie,a,b,j,st,x,wq,val,grad,iv,di,rr,grr,lam,glam, &
        !$omp tau,g,dg,gi,dgi,r,dr,sr,si,gsr,gsi,sig,gs,h) &
        !$omp reduction(+:total) reduction(max:bad)
        allocate(val(nl), grad(2, nl), iv(m), di(m))
        !$omp do schedule(static)
        do ic = 1, nactive
            t = active(ic)
            do q = 1, fsd%nq
                call meridian_qgrad(fsd, t, q, x, wq, val, grad, st)
                if (st /= 0) then
                    bad = 1
                    cycle
                end if
                lam = matmul(lamv, val)
                glam = matmul(grad, transpose(lamv))
                sig = (0.0_dp, 0.0_dp)
                gs = (0.0_dp, 0.0_dp)
                do e = 1, 3
                    ie = lookup(e, t)
                    if (ie == 0) cycle
                    a = e
                    b = mod(e, 3) + 1
                    if (lam(a) + lam(b) <= 0.0_dp) then
                        bad = 1
                        cycle
                    end if
                    call generated_meridian_edge_coordinate(lam(a), lam(b), tau)
                    call lagrange_weights(m, tau, nodes, iv)
                    call lagrange_deriv_weights(m, tau, nodes, di)
                    gi = sum(coeff(:, ie)*iv)
                    dgi = sum(coeff(:, ie)*di)
                    call load(t, e, tau, g)
                    call load_derivative(t, e, tau, dg)
                    if (.not. finite_complex(g)) bad = 1
                    if (.not. finite_complex(dg)) bad = 1
                    r = g - gi
                    dr = dg - dgi
                    call generated_meridian_edge_lift(lam(a), lam(b), &
                        glam(1, a), glam(2, a), glam(1, b), glam(2, b), &
                        real(r, dp), real(dr, dp), sr, gsr(1), gsr(2))
                    call generated_meridian_edge_lift(lam(a), lam(b), &
                        glam(1, a), glam(2, a), glam(1, b), glam(2, b), &
                        aimag(r), aimag(dr), si, gsi(1), gsi(2))
                    sig = sig + cmplx(sr, si, dp)
                    gs = gs + cmplx(gsr, gsi, dp)
                end do
                rr = x(1)
                grr = [1.0_dp, 0.0_dp]
                if (fsd%cyl_r0 > 0.0_dp) then
                    rr = fsd%cyl_r0
                    grr = 0.0_dp
                end if
                call scale_potential(rr, grr, sig, gs)
                call residual(real(n, dp), rr, grr, 0.0_dp, [0.0_dp, 0.0_dp], &
                    sig, gs, zero2, (0.0_dp, 0.0_dp), zero2, h)
                total = total + wq*real(dot_product(h, h), dp)
            end do
        end do
        !$omp end do
        !$omp end parallel
        if (bad /= 0) return
        if (.not. ieee_is_finite(total)) return
        bound = sqrt(max(0.0_dp, total))
        status = HCURL_OK
    end subroutine meridian_hcurl_edge_remainder

    pure logical function finite_complex(z) result(ok)
        complex(dp), intent(in) :: z
        ok = ieee_is_finite(real(z, dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex

    subroutine meridian_hcurl_equilibrated_bound(fs, fsd, n, a_h, source, est, &
            status, method, edge_cell, edge_local, edge_load, free_trace, side)
        !! Bound on ||curl_n (A - A_h)||_M2 for the gauge-fixed problem with
        !! mu0-scaled current source(x) = J (components 1 = R, 3 = Z used).
        !! fs: primal complex (V1 holds a_h); fsd: complex whose Lagrange
        !! space V0 carries the potential s (same mesh and quadrature).
        !! method: HCURL_GLOBAL (default) or HCURL_PATCHWISE.
        !! Optional edge loads (all three together): the oriented edges
        !! (edge_cell(i), edge_local(i)) carry the density edge_load (called
        !! serially).  Interior loaded edges must separate the mesh: side(nt)
        !! marks the cells on one side (default: the cells connected to the
        !! first interior edge_cell without crossing a loaded edge).  Boundary
        !! loaded edges require free_trace (default .false.: v x n = 0, the
        !! tangential trace of A is fixed; .true.: the trace is free and sigma
        !! is prescribed on the whole boundary, zero on unloaded edges).
        type(meridian_feec_t), intent(in) :: fs, fsd
        integer, intent(in) :: n
        complex(dp), intent(in) :: a_h(:)
        procedure(meridian_field_i) :: source
        type(meridian_hcurl_estimate_t), intent(out) :: est
        integer, intent(out) :: status
        integer, intent(in), optional :: method
        integer, intent(in), optional :: edge_cell(:), edge_local(:)
        procedure(meridian_edge_load_i), optional :: edge_load
        logical, intent(in), optional :: free_trace
        logical, intent(in), optional :: side(:)

        complex(dp), allocatable :: c2(:), s(:), ke(:, :, :), be(:, :, :)
        complex(dp), allocatable :: wfix(:, :)
        real(dp), allocatable :: lamv(:, :), cell2(:)
        logical, allocatable :: fixed(:)
        integer, allocatable :: vstart(:), vcell(:), vslot(:)
        type(csc_z_t) :: op
        real(dp) :: amax, tot, fnorm
        integer :: nt, nl0, st, bad, mth
        logical :: loaded, ftrace, scaled

        status = HCURL_BAD_INPUT
        mth = HCURL_GLOBAL
        if (present(method)) mth = method
        if (mth /= HCURL_GLOBAL .and. mth /= HCURL_PATCHWISE) return
        if (n == 0) return
        loaded = present(edge_cell)
        if (present(edge_local) .neqv. loaded) return
        if (present(edge_load) .neqv. loaded) return
        ftrace = .false.
        if (present(free_trace)) ftrace = free_trace
        scaled = loaded .or. ftrace
        if (scaled .and. mth /= HCURL_GLOBAL) return
        call check_consistent(fs, fsd, status)
        if (status /= HCURL_OK) return
        status = HCURL_BAD_INPUT
        if (size(a_h) /= fs%sp(1)%ndof) return
        amax = maxval(abs(a_h))
        if (fs%n0 > 0) then
            if (maxval(abs(a_h(fs%nn + 1:))) > 1.0e-12_dp*max(amax, 1.0e-300_dp)) &
                then
                status = HCURL_NOT_GAUGE_FIXED
                return
            end if
        end if

        nt = fs%mm%topology%n_triangles
        nl0 = fsd%nl0
        call meridian_op_complex(fs%curl, n, op, st)
        if (st /= 0) return
        c2 = csc_matvec(op, a_h)
        call vertex_weights(fsd, lamv)
        allocate (wfix(nl0, nt), fixed(fsd%n0))
        wfix = (0.0_dp, 0.0_dp)
        fixed = .false.
        est%scaled_potential = scaled
        est%trace_fixed = ftrace
        if (ftrace) fixed = fsd%sp(0)%bnd
        if (loaded) then
            call build_lifts(edge_load, st)
            if (st /= 0) then
                status = HCURL_BAD_INPUT
                return
            end if
        end if

        allocate (ke(nl0, nl0, nt), be(nl0, 0:3, nt))
        call element_systems(bad)
        if (bad /= 0) then
            status = HCURL_BAD_CELL
            return
        end if
        allocate (s(fsd%n0))
        if (mth == HCURL_GLOBAL) then
            call solve_global(st)
        else
            call solve_patchwise(st)
        end if
        if (st /= 0) then
            status = HCURL_SOLVE_FAILED
            return
        end if

        allocate (cell2(nt), est%cell_indicator(nt))
        call evaluate_bound(tot, fnorm, bad)
        if (bad /= 0) then
            status = HCURL_BAD_CELL
            return
        end if
        est%cell_indicator = sqrt(cell2)
        est%bound = sqrt(tot)
        est%field_norm = sqrt(fnorm)
        call move_alloc(s, est%potential)
        status = HCURL_OK

    contains

        subroutine element_systems(bad_out)
            !! Element Gram matrices ke(:, :, t) of the residual map s ->
            !! (r_pol, r_phi) and right-hand sides be(:, 0, t) (psi = 1) and
            !! be(:, slot, t) (hat function of the triangle vertex slot).
            integer, intent(out) :: bad_out
            real(dp), allocatable :: phi(:, :), val(:), grad(:, :)
            complex(dp) :: g(3, nl0), d1(3), d2(3, 2), ds(3), jp(2), jv(3), bp(2), &
                bh(3), r(3), zero2(2), sf, gsf(2)
            real(dp) :: x(2), xq(2), wq, wq2, rr, grr(2), psi, gpsi(2)
            integer :: t, q, i, j, sl, stq

            bad_out = 0
            zero2 = (0.0_dp, 0.0_dp)
            !$omp parallel default(shared) &
            !$omp private(t, q, i, j, sl, stq, phi, val, grad, g, d1, d2, ds, jp, &
            !$omp jv, bp, bh, r, x, xq, wq, wq2, rr, grr, psi, gpsi, sf, gsf) &
            !$omp reduction(max:bad_out)
            allocate (phi(3, fs%sp(2)%nloc), val(nl0), grad(2, nl0))
            !$omp do schedule(static)
            do t = 1, nt
                ke(:, :, t) = (0.0_dp, 0.0_dp)
                be(:, :, t) = (0.0_dp, 0.0_dp)
                do q = 1, fsd%nq
                    call meridian_qgrad(fsd, t, q, x, wq, val, grad, stq)
                    if (stq /= 0) then
                        bad_out = 1
                        cycle
                    end if
                    call meridian_qpoint(fs, 2, t, q, xq, wq2, phi, stq)
                    if (stq /= 0) then
                        bad_out = 1
                        cycle
                    end if
                    call field_at(t, phi, bp, bh)
                    call source(x, jv)
                    jp = [jv(1), jv(3)]
                    call radius(x, rr, grr)
                    ! fixed lift of the potential (zero without loads)
                    sf = sum(wfix(:, t)*val)
                    gsf = matmul(grad, wfix(:, t))
                    if (scaled) then
                        call scale_potential(rr, grr, sf, gsf)
                        call scale_basis(rr, grr, val, grad)
                    end if
                    ! data part: d = psi D1 + d_R(psi) D2(:, 1) + d_Z(psi) D2(:, 2)
                    call residual(real(n, dp), rr, grr, 1.0_dp, [0.0_dp, 0.0_dp], &
                        sf, gsf, jp, bh(3), bp, d1)
                    call residual(real(n, dp), rr, grr, 0.0_dp, [1.0_dp, 0.0_dp], &
                        (0.0_dp, 0.0_dp), zero2, zero2, bh(3), zero2, d2(:, 1))
                    call residual(real(n, dp), rr, grr, 0.0_dp, [0.0_dp, 1.0_dp], &
                        (0.0_dp, 0.0_dp), zero2, zero2, bh(3), zero2, d2(:, 2))
                    do j = 1, nl0
                        call residual(real(n, dp), rr, grr, 0.0_dp, &
                            [0.0_dp, 0.0_dp], cmplx(val(j), 0.0_dp, dp), &
                            cmplx(grad(:, j), 0.0_dp, dp), zero2, (0.0_dp, 0.0_dp), &
                            zero2, r)
                        g(:, j) = r
                    end do
                    do j = 1, nl0
                        do i = 1, nl0
                            ke(i, j, t) = ke(i, j, t) &
                                + wq*dot_product(g(:, i), g(:, j))
                        end do
                        be(j, 0, t) = be(j, 0, t) - wq*dot_product(g(:, j), d1)
                    end do
                    do sl = 1, 3
                        psi = 0.0_dp
                        gpsi = 0.0_dp
                        do j = 1, nl0
                            psi = psi + lamv(sl, j)*val(j)
                            gpsi = gpsi + lamv(sl, j)*grad(:, j)
                        end do
                        ds = psi*d1 + gpsi(1)*d2(:, 1) + gpsi(2)*d2(:, 2)
                        do j = 1, nl0
                            be(j, sl, t) = be(j, sl, t) - wq*dot_product(g(:, j), ds)
                        end do
                    end do
                end do
            end do
            !$omp end do
            !$omp end parallel
        end subroutine element_systems

        subroutine solve_global(info)
            !! Minimizer of ||H(s) - B_h|| over the continuous potential space:
            !! sparse Hermitian positive definite normal equations.
            integer, intent(out) :: info
            type(csc_z_t) :: gm
            type(fortsparse_status_t) :: fst
            complex(dp), allocatable :: rhs(:), vals(:)
            integer, allocatable :: ri(:), ci(:)
            integer :: t, i, j, cnt

            allocate (ri(nl0*nl0*nt + fsd%n0), ci(nl0*nl0*nt + fsd%n0))
            allocate (vals(nl0*nl0*nt + fsd%n0))
            allocate (rhs(fsd%n0))
            rhs = (0.0_dp, 0.0_dp)
            cnt = 0
            do t = 1, nt
                do j = 1, nl0
                    if (fixed(fsd%sp(0)%dof(j, t))) cycle
                    rhs(fsd%sp(0)%dof(j, t)) = rhs(fsd%sp(0)%dof(j, t)) + be(j, 0, t)
                    do i = 1, nl0
                        if (fixed(fsd%sp(0)%dof(i, t))) cycle
                        cnt = cnt + 1
                        ri(cnt) = fsd%sp(0)%dof(i, t)
                        ci(cnt) = fsd%sp(0)%dof(j, t)
                        vals(cnt) = ke(i, j, t)
                    end do
                end do
            end do
            ! Fixed (boundary) dofs of sigma_0 vanish: identity rows.
            do i = 1, fsd%n0
                if (.not. fixed(i)) cycle
                cnt = cnt + 1
                ri(cnt) = i
                ci(cnt) = i
                vals(cnt) = (1.0_dp, 0.0_dp)
            end do
            info = 1
            call csc_from_triplet(fsd%n0, fsd%n0, ri(:cnt), ci(:cnt), vals(:cnt), &
                gm, fst)
            if (.not. status_ok(fst)) return
            call sparse_solve_once(gm, rhs, s, fst)
            if (.not. status_ok(fst)) return
            info = 0
        end subroutine solve_global

        subroutine solve_patchwise(info)
            !! Patch minimizers s_a (hat-function weighted data, zero at the
            !! nodes of the edge opposite to a), summed over the vertices a.
            integer, intent(out) :: info
            complex(dp), allocatable :: pcoef(:), gmat(:, :), rhs(:)
            integer, allocatable :: vstart(:), vcell(:), vslot(:), pstart(:), &
                psize(:), pdof(:), loc(:), lpd(:), lidx(:, :)
            integer :: nv, a, i, j, k, t, slot, np, ncell, base, ierr, ii, jj
            integer :: failed

            nv = fs%mm%topology%n_vertices
            call vertex_stars(nv, fs%mm%topology%triangles, vstart, vcell, &
                vslot)
            allocate (pstart(nv + 1), psize(nv))
            pstart(1) = 1
            do a = 1, nv
                pstart(a + 1) = pstart(a) + (vstart(a + 1) - vstart(a))*nl0
            end do
            allocate (pdof(pstart(nv + 1) - 1), pcoef(pstart(nv + 1) - 1))
            failed = 0
            !$omp parallel default(shared) &
            !$omp private(a, i, j, k, t, slot, np, ncell, base, ierr, ii, jj, loc, &
            !$omp lpd, lidx, gmat, rhs) reduction(max:failed)
            allocate (loc(fsd%n0))
            loc = 0
            !$omp do schedule(dynamic, 16)
            do a = 1, nv
                ncell = vstart(a + 1) - vstart(a)
                allocate (lpd(ncell*nl0), lidx(nl0, ncell))
                np = 0
                do k = 1, ncell
                    t = vcell(vstart(a) + k - 1)
                    slot = vslot(vstart(a) + k - 1)
                    do j = 1, nl0
                        lidx(j, k) = 0
                        if (lamv(slot, j) <= 0.0_dp) cycle
                        ii = fsd%sp(0)%dof(j, t)
                        if (loc(ii) == 0) then
                            np = np + 1
                            loc(ii) = np
                            lpd(np) = ii
                        end if
                        lidx(j, k) = loc(ii)
                    end do
                end do
                allocate (gmat(np, np), rhs(np))
                gmat = (0.0_dp, 0.0_dp)
                rhs = (0.0_dp, 0.0_dp)
                do k = 1, ncell
                    t = vcell(vstart(a) + k - 1)
                    slot = vslot(vstart(a) + k - 1)
                    do j = 1, nl0
                        if (lidx(j, k) == 0) cycle
                        rhs(lidx(j, k)) = rhs(lidx(j, k)) + be(j, slot, t)
                        do i = 1, nl0
                            if (lidx(i, k) == 0) cycle
                            gmat(lidx(i, k), lidx(j, k)) = &
                                gmat(lidx(i, k), lidx(j, k)) + ke(i, j, t)
                        end do
                    end do
                end do
                call zposv('U', np, 1, gmat, np, rhs, np, ierr)
                if (ierr /= 0) then
                    failed = 1
                    rhs = (0.0_dp, 0.0_dp)
                end if
                base = pstart(a) - 1
                psize(a) = np
                do j = 1, np
                    pdof(base + j) = lpd(j)
                    pcoef(base + j) = rhs(j)
                    loc(lpd(j)) = 0
                end do
                deallocate (lpd, lidx, gmat, rhs)
            end do
            !$omp end do
            !$omp end parallel
            info = failed
            if (failed /= 0) return
            s = (0.0_dp, 0.0_dp)
            do a = 1, nv
                do i = pstart(a), pstart(a) + psize(a) - 1
                    jj = pdof(i)
                    s(jj) = s(jj) + pcoef(i)
                end do
            end do
        end subroutine solve_patchwise

        subroutine evaluate_bound(total, bnorm, bad_out)
            !! ||H(s) - B_h||^2 and ||B_h||^2 on the assembled s.
            real(dp), intent(out) :: total, bnorm
            integer, intent(out) :: bad_out
            real(dp), allocatable :: phi(:, :), val(:), grad(:, :)
            complex(dp) :: r(3), jp(2), jv(3), bp(2), bh(3), sq, gsq(2)
            real(dp) :: x(2), xq(2), wq, wq2, rr, grr(2), acc
            integer :: t, q, j, stq

            total = 0.0_dp
            bnorm = 0.0_dp
            bad_out = 0
            !$omp parallel default(shared) &
            !$omp private(t, q, j, phi, val, grad, r, jp, jv, bp, bh, sq, gsq, &
            !$omp x, xq, wq, wq2, rr, grr, acc, stq) &
            !$omp reduction(+:total, bnorm) reduction(max:bad_out)
            allocate (phi(3, fs%sp(2)%nloc), val(nl0), grad(2, nl0))
            !$omp do schedule(static)
            do t = 1, nt
                acc = 0.0_dp
                do q = 1, fsd%nq
                    call meridian_qgrad(fsd, t, q, x, wq, val, grad, stq)
                    if (stq /= 0) then
                        bad_out = 1
                        cycle
                    end if
                    call meridian_qpoint(fs, 2, t, q, xq, wq2, phi, stq)
                    if (stq /= 0) then
                        bad_out = 1
                        cycle
                    end if
                    call field_at(t, phi, bp, bh)
                    call source(x, jv)
                    jp = [jv(1), jv(3)]
                    call radius(x, rr, grr)
                    sq = (0.0_dp, 0.0_dp)
                    gsq = (0.0_dp, 0.0_dp)
                    do j = 1, nl0
                        sq = sq + (s(fsd%sp(0)%dof(j, t)) + wfix(j, t))*val(j)
                        gsq = gsq + (s(fsd%sp(0)%dof(j, t)) + wfix(j, t))*grad(:, j)
                    end do
                    if (scaled) call scale_potential(rr, grr, sq, gsq)
                    call residual(real(n, dp), rr, grr, 1.0_dp, [0.0_dp, 0.0_dp], &
                        sq, gsq, jp, bh(3), bp, r)
                    acc = acc + wq*sum(abs(r)**2)
                    bnorm = bnorm + wq*(abs(bp(1))**2 + abs(bp(2))**2 &
                        + abs(bh(3))**2)
                end do
                cell2(t) = acc
                total = total + acc
            end do
            !$omp end do
            !$omp end parallel
        end subroutine evaluate_bound

        subroutine build_lifts(load, info)
            !! Broken lift wfix of the potential sigma = R s from the edge
            !! loads: -G (left) / +G (right) on the side cells of interior
            !! loaded edges, -G minus the side lift on boundary loaded edges;
            !! interpolated at the V0 nodes of the edges (see the header).
            procedure(meridian_edge_load_i) :: load
            integer, intent(out) :: info
            type(triangle_lagrange_basis_t) :: basis
            integer, allocatable :: nb(:), stack(:), hits(:)
            logical, allocatable :: region(:), cut(:, :)
            complex(dp), allocatable :: wint(:), wbnd(:), first(:)
            real(dp), allocatable :: tq(:), wtq(:), bval(:), bgrad(:, :)
            complex(dp) :: g, gi
            real(dp) :: xi, eta, defect2, norm2_, gmax
            integer :: ne, ie, t, e, j, nv, top, c, ee, cn, dof, start, pass, q
            integer :: ng, stb

            info = 1
            ne = size(edge_cell)
            if (size(edge_local) /= ne .or. ne < 1) return
            if (any(edge_cell < 1) .or. any(edge_cell > nt)) return
            if (any(edge_local < 1) .or. any(edge_local > 3)) return
            nv = fsd%mm%topology%n_vertices
            call vertex_stars(nv, fsd%mm%topology%triangles(1:3, :), vstart, &
                vcell, vslot)
            allocate (nb(ne), cut(3, nt))
            cut = .false.
            do ie = 1, ne
                nb(ie) = neighbour(edge_cell(ie), edge_local(ie), ee)
                cut(edge_local(ie), edge_cell(ie)) = .true.
                if (nb(ie) > 0) cut(ee, nb(ie)) = .true.
            end do
            if (.not. ftrace .and. any(nb == 0)) return
            allocate (region(nt))
            region = .false.
            if (any(nb > 0)) then
                if (present(side)) then
                    if (size(side) /= nt) return
                    region = side
                else
                    ! flood fill from the first interior loaded edge
                    start = edge_cell(findloc(nb > 0, .true., dim=1))
                    allocate (stack(nt))
                    top = 1
                    stack(1) = start
                    region(start) = .true.
                    do while (top > 0)
                        c = stack(top)
                        top = top - 1
                        do e = 1, 3
                            if (cut(e, c)) cycle
                            cn = neighbour(c, e, ee)
                            if (cn == 0) cycle
                            if (region(cn)) cycle
                            region(cn) = .true.
                            top = top + 1
                            stack(top) = cn
                        end do
                    end do
                end if
                do ie = 1, ne
                    if (nb(ie) == 0) cycle
                    if (region(edge_cell(ie)) .eqv. region(nb(ie))) return
                end do
            end if
            allocate (wint(fsd%n0), wbnd(fsd%n0), first(fsd%n0), hits(fsd%n0))
            wint = (0.0_dp, 0.0_dp)
            wbnd = (0.0_dp, 0.0_dp)
            gmax = 0.0_dp
            ! pass 1: interior jumps, pass 2: boundary traces
            do pass = 1, 2
                hits = 0
                do ie = 1, ne
                    if ((nb(ie) > 0) .neqv. (pass == 1)) cycle
                    t = edge_cell(ie)
                    e = edge_local(ie)
                    do j = 1, nl0
                        if (lamv(mod(e + 1, 3) + 1, j) > 0.0_dp) cycle
                        call load(t, e, edge_tau(e, j), g)
                        gmax = max(gmax, abs(g))
                        dof = fsd%sp(0)%dof(j, t)
                        if (pass == 1) then
                            if (.not. region(t)) g = -g
                        else if (region(t)) then
                            g = g + wint(dof)
                        end if
                        if (hits(dof) == 0) first(dof) = -g
                        est%lift_mismatch = max(est%lift_mismatch, &
                            abs(-g - first(dof)))
                        hits(dof) = hits(dof) + 1
                        if (pass == 1) then
                            wint(dof) = wint(dof) + (-g - wint(dof))/hits(dof)
                        else
                            wbnd(dof) = wbnd(dof) + (-g - wbnd(dof))/hits(dof)
                        end if
                    end do
                end do
            end do
            if (gmax > 0.0_dp) est%lift_mismatch = est%lift_mismatch/gmax
            do t = 1, nt
                do j = 1, nl0
                    dof = fsd%sp(0)%dof(j, t)
                    wfix(j, t) = wbnd(dof)
                    if (region(t)) wfix(j, t) = wfix(j, t) + wint(dof)
                end do
            end do
            ! load_defect: G - I G at Gauss points of every loaded edge
            call initialize_triangle_lagrange_basis(fsd%k, basis, stb)
            if (stb /= 0) return
            ng = fsd%k + 2
            allocate (tq(ng), wtq(ng), bval(nl0), bgrad(2, nl0))
            call gauss_legendre_ab(ng, 0.0_dp, 1.0_dp, tq, wtq)
            defect2 = 0.0_dp
            norm2_ = 0.0_dp
            do ie = 1, ne
                t = edge_cell(ie)
                e = edge_local(ie)
                do q = 1, ng
                    call edge_reference_point(e, tq(q), xi, eta)
                    call evaluate_triangle_lagrange_basis(basis, xi, eta, bval, &
                        bgrad, stb)
                    if (stb /= 0) return
                    gi = (0.0_dp, 0.0_dp)
                    do j = 1, nl0
                        if (lamv(mod(e + 1, 3) + 1, j) > 0.0_dp) cycle
                        call load(t, e, edge_tau(e, j), g)
                        gi = gi + g*bval(j)
                    end do
                    call load(t, e, tq(q), g)
                    defect2 = defect2 + wtq(q)*abs(g - gi)**2
                    norm2_ = norm2_ + wtq(q)*abs(g)**2
                end do
            end do
            if (norm2_ > 0.0_dp) est%load_defect = sqrt(defect2/norm2_)
            info = 0
        end subroutine build_lifts

        integer function neighbour(t, e, e_other) result(cn)
            !! Cell across local edge e of t (0 on the boundary) and the local
            !! index e_other of that edge in it.
            integer, intent(in) :: t, e
            integer, intent(out) :: e_other
            integer :: a, b, k2, c, sl

            a = fs%mm%topology%triangles(e, t)
            b = fs%mm%topology%triangles(mod(e, 3) + 1, t)
            cn = 0
            e_other = 0
            do k2 = vstart(b), vstart(b + 1) - 1
                c = vcell(k2)
                if (c == t) cycle
                sl = vslot(k2)
                ! edge b -> a in c: local edge sl ends at mod(sl, 3) + 1
                if (fs%mm%topology%triangles(mod(sl, 3) + 1, c) == a) then
                    cn = c
                    e_other = sl
                end if
            end do
        end function neighbour

        real(dp) function edge_tau(e, j) result(tau)
            !! Parameter of local node j on local edge e: the barycentric
            !! coordinate of the end vertex mod(e, 3) + 1.
            integer, intent(in) :: e, j
            tau = lamv(mod(e, 3) + 1, j)
        end function edge_tau

        subroutine field_at(t, phi, bp_out, bh_out)
            !! Poloidal (R, Z) components bp_out and toroidal component
            !! bh_out(3) of B_h = curl_n A_h (V2 coefficients c2) at a
            !! quadrature point with local basis phi (components R, phi, Z).
            integer, intent(in) :: t
            real(dp), intent(in) :: phi(:, :)
            complex(dp), intent(out) :: bp_out(2), bh_out(3)
            complex(dp) :: cl
            integer :: i2

            bh_out = (0.0_dp, 0.0_dp)
            do i2 = 1, fs%sp(2)%nloc
                cl = real(fs%sp(2)%sgn(i2, t), dp)*c2(fs%sp(2)%dof(i2, t))
                bh_out = bh_out + cl*phi(:, i2)
            end do
            bp_out = [bh_out(1), bh_out(3)]
            bh_out(3) = bh_out(2)
        end subroutine field_at

        subroutine radius(x, rr, grr)
            real(dp), intent(in) :: x(2)
            real(dp), intent(out) :: rr, grr(2)

            if (fs%cyl_r0 > 0.0_dp) then
                rr = fs%cyl_r0
                grr = 0.0_dp
            else
                rr = x(1)
                grr = [1.0_dp, 0.0_dp]
            end if
        end subroutine radius

    end subroutine meridian_hcurl_equilibrated_bound

    pure subroutine scale_potential(rr, grr, s, gs)
        !! (sigma, grad sigma) -> (s, grad s) for s = sigma/R (generated leaf,
        !! real and imaginary parts).
        real(dp), intent(in) :: rr, grr(2)
        complex(dp), intent(inout) :: s, gs(2)
        real(dp) :: re(3), im(3)

        call generated_meridian_hcurl_scaled_potential(rr, grr(1), grr(2), &
            real(s, dp), real(gs(1), dp), real(gs(2), dp), re(1), re(2), re(3))
        call generated_meridian_hcurl_scaled_potential(rr, grr(1), grr(2), &
            aimag(s), aimag(gs(1)), aimag(gs(2)), im(1), im(2), im(3))
        s = cmplx(re(1), im(1), dp)
        gs = cmplx(re(2:3), im(2:3), dp)
    end subroutine scale_potential

    pure subroutine scale_basis(rr, grr, val, grad)
        !! Basis values and gradients of sigma -> those of s = sigma/R.
        real(dp), intent(in) :: rr, grr(2)
        real(dp), intent(inout) :: val(:), grad(:, :)
        real(dp) :: o(3)
        integer :: j

        do j = 1, size(val)
            call generated_meridian_hcurl_scaled_potential(rr, grr(1), grr(2), &
                val(j), grad(1, j), grad(2, j), o(1), o(2), o(3))
            val(j) = o(1)
            grad(:, j) = o(2:3)
        end do
    end subroutine scale_basis

    pure subroutine edge_reference_point(edge, tau, xi, eta)
        !! Reference point at parameter tau of local edge edge (vertex edge ->
        !! vertex mod(edge, 3) + 1; vertices (0,0), (1,0), (0,1)).
        integer, intent(in) :: edge
        real(dp), intent(in) :: tau
        real(dp), intent(out) :: xi, eta

        select case (edge)
        case (1)
            xi = tau
            eta = 0.0_dp
        case (2)
            xi = 1.0_dp - tau
            eta = tau
        case default
            xi = 0.0_dp
            eta = 1.0_dp - tau
        end select
    end subroutine edge_reference_point

    pure subroutine residual(n, rr, grr, psi, gpsi, s, gs, jp, sh, bp, r)
        !! Generated pointwise residual (complex wrapper of the FortSym leaf).
        real(dp), intent(in) :: n, rr, grr(2), psi, gpsi(2)
        complex(dp), intent(in) :: s, gs(2), jp(2), sh, bp(2)
        complex(dp), intent(out) :: r(3)
        real(dp) :: o(6)

        call generated_meridian_hcurl_residual(n, rr, grr(1), grr(2), psi, &
            gpsi(1), gpsi(2), real(s, dp), aimag(s), real(gs(1), dp), &
            aimag(gs(1)), real(gs(2), dp), aimag(gs(2)), real(jp(1), dp), &
            aimag(jp(1)), real(jp(2), dp), aimag(jp(2)), real(sh, dp), &
            aimag(sh), real(bp(1), dp), aimag(bp(1)), real(bp(2), dp), &
            aimag(bp(2)), o(1), o(2), o(3), o(4), o(5), o(6))
        r(1) = cmplx(o(1), o(2), dp)
        r(2) = cmplx(o(3), o(4), dp)
        r(3) = cmplx(o(5), o(6), dp)
    end subroutine residual

    subroutine vertex_weights(fsd, lamv)
        !! P1 hat-function weights lamv(slot, j): the barycentric coordinate of
        !! triangle vertex slot at local Lagrange node j of the potential
        !! space (slot 1 = (0,0), 2 = (1,0), 3 = (0,1)).
        type(meridian_feec_t), intent(in) :: fsd
        real(dp), allocatable, intent(out) :: lamv(:, :)
        type(triangle_lagrange_basis_t) :: basis
        real(dp), allocatable :: nodes(:, :)
        integer :: status, j

        call initialize_triangle_lagrange_basis(fsd%k, basis, status)
        call triangle_lagrange_nodes(basis, nodes)
        allocate (lamv(3, size(nodes, 2)))
        do j = 1, size(nodes, 2)
            lamv(1, j) = 1.0_dp - nodes(1, j) - nodes(2, j)
            lamv(2, j) = nodes(1, j)
            lamv(3, j) = nodes(2, j)
        end do
        ! snap round-off so that nodes on the opposite edge are exactly zero
        where (abs(lamv) < 0.25_dp/real(fsd%k, dp)) lamv = 0.0_dp
    end subroutine vertex_weights

    subroutine check_consistent(fs, fsd, status)
        !! fs and fsd must share mesh, geometry order and quadrature rule.
        type(meridian_feec_t), intent(in) :: fs, fsd
        integer, intent(out) :: status

        status = HCURL_BAD_INPUT
        if (fs%mm%topology%n_triangles /= fsd%mm%topology%n_triangles) return
        if (fs%mm%topology%n_vertices /= fsd%mm%topology%n_vertices) return
        if (fs%g /= fsd%g) return
        if (fs%nq /= fsd%nq) return
        if (abs(fs%cyl_r0 - fsd%cyl_r0) > 0.0_dp) return
        if (fs%k < 1 .or. fsd%k < 1) return
        if (maxval(abs(fs%qw - fsd%qw)) > 1.0e-14_dp) return
        if (any(fs%mm%topology%triangles /= fsd%mm%topology%triangles)) return
        if (maxval(abs(fs%mm%nodes - fsd%mm%nodes)) > 0.0_dp) return
        status = HCURL_OK
    end subroutine check_consistent

end module fortfem_equilibrated_hcurl_meridian_n
