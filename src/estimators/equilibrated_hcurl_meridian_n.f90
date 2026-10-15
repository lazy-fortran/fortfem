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
    !! Consumers.  The Steklov (exterior) term is not part of this module; it
    !! adds a boundary functional to the right-hand side, which would
    !! constrain R s on the boundary (a prescribed boundary trace of s) and
    !! drop the free boundary treatment.  Both solvers could take prescribed
    !! values for the boundary nodes of s; that hook is not implemented.
    !!
    !! Requirements: fs and fsd are built on the same mesh (equal geometry
    !! order) with the same quadrature degree; the source callback is called
    !! from OpenMP threads and must be thread-safe.  The integrals use that
    !! rule; their quadrature error is not enclosed.
    use fortfem_kinds, only: dp
    use fortfem_generated_meridian_hcurl_equilibration, only: &
        generated_meridian_hcurl_residual
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_field_i, &
        meridian_op_complex, meridian_qgrad, meridian_qpoint
    use fortfem_equilibrated_flux_topology, only: vertex_stars
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        initialize_triangle_lagrange_basis, triangle_lagrange_basis_t, &
        triangle_lagrange_nodes
    use fortsparse, only: csc_from_triplet, csc_matvec, csc_z_t, &
        fortsparse_status_t, sparse_solve_once, status_ok
    implicit none
    private

    public :: meridian_hcurl_estimate_t, meridian_hcurl_equilibrated_bound

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
        !> Coefficients of s = H_phi in V0 of fsd.
        complex(dp), allocatable :: potential(:)
    end type meridian_hcurl_estimate_t

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

    subroutine meridian_hcurl_equilibrated_bound(fs, fsd, n, a_h, source, est, &
            status, method)
        !! Bound on ||curl_n (A - A_h)||_M2 for the gauge-fixed problem with
        !! mu0-scaled current source(x) = J (components 1 = R, 3 = Z used).
        !! fs: primal complex (V1 holds a_h); fsd: complex whose Lagrange
        !! space V0 carries the potential s (same mesh and quadrature).
        !! method: HCURL_GLOBAL (default) or HCURL_PATCHWISE.
        type(meridian_feec_t), intent(in) :: fs, fsd
        integer, intent(in) :: n
        complex(dp), intent(in) :: a_h(:)
        procedure(meridian_field_i) :: source
        type(meridian_hcurl_estimate_t), intent(out) :: est
        integer, intent(out) :: status
        integer, intent(in), optional :: method

        complex(dp), allocatable :: c2(:), s(:), ke(:, :, :), be(:, :, :)
        real(dp), allocatable :: lamv(:, :), cell2(:)
        type(csc_z_t) :: op
        real(dp) :: amax, tot, fnorm
        integer :: nt, nl0, st, bad, mth

        status = HCURL_BAD_INPUT
        mth = HCURL_GLOBAL
        if (present(method)) mth = method
        if (mth /= HCURL_GLOBAL .and. mth /= HCURL_PATCHWISE) return
        if (n == 0) return
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
                bh(3), r(3), zero2(2)
            real(dp) :: x(2), xq(2), wq, wq2, rr, grr(2), psi, gpsi(2)
            integer :: t, q, i, j, sl, stq

            bad_out = 0
            zero2 = (0.0_dp, 0.0_dp)
            !$omp parallel default(shared) &
            !$omp private(t, q, i, j, sl, stq, phi, val, grad, g, d1, d2, ds, jp, &
            !$omp jv, bp, bh, r, x, xq, wq, wq2, rr, grr, psi, gpsi) &
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
                    ! data part: d = psi D1 + d_R(psi) D2(:, 1) + d_Z(psi) D2(:, 2)
                    call residual(real(n, dp), rr, grr, 1.0_dp, [0.0_dp, 0.0_dp], &
                        (0.0_dp, 0.0_dp), zero2, jp, bh(3), bp, d1)
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

            allocate (ri(nl0*nl0*nt), ci(nl0*nl0*nt), vals(nl0*nl0*nt))
            allocate (rhs(fsd%n0))
            rhs = (0.0_dp, 0.0_dp)
            cnt = 0
            do t = 1, nt
                do j = 1, nl0
                    rhs(fsd%sp(0)%dof(j, t)) = rhs(fsd%sp(0)%dof(j, t)) + be(j, 0, t)
                    do i = 1, nl0
                        cnt = cnt + 1
                        ri(cnt) = fsd%sp(0)%dof(i, t)
                        ci(cnt) = fsd%sp(0)%dof(j, t)
                        vals(cnt) = ke(i, j, t)
                    end do
                end do
            end do
            info = 1
            call csc_from_triplet(fsd%n0, fsd%n0, ri, ci, vals, gm, fst)
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
                        sq = sq + s(fsd%sp(0)%dof(j, t))*val(j)
                        gsq = gsq + s(fsd%sp(0)%dof(j, t))*grad(:, j)
                    end do
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
