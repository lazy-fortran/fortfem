module fortfem_hcurl_equilibrated_2d
    !! Patchwise reaction-type a posteriori energy bound for the planar
    !! H(curl) problem curl(nu curl u) + kappa u = f on an affine CCW triangle
    !! mesh with u x n = 0, scalar nu > 0 and kappa > 0, for a discrete u_h in
    !! the order-one Nedelec space (cellwise constant curl).
    !!
    !! With s = nu curl u, rot h = (d_y h, -d_x h) and g = f - kappa u_h, the
    !! error e = u - u_h has zero tangential trace, so for every h in H^1
    !! (no boundary condition on h, since int curl(e) h = int e . rot h)
    !!
    !!   int nu curl e curl e + kappa |e|^2
    !!       = int (g - rot h) . e + int (h - nu curl u_h) curl e,
    !!
    !! and Cauchy--Schwarz gives the bound
    !!
    !!   |||u - u_h||| <= eta,
    !!   eta^2 = ||nu^(-1/2)(h - nu curl u_h)||^2 + ||kappa^(-1/2)(g - rot h)||^2,
    !!
    !! |||v|||^2 = int nu (curl v)^2 + kappa |v|^2. No Friedrichs constant is
    !! needed: kappa > 0 controls the L2 part. Any h gives a valid bound.
    !!
    !! The reconstruction h is continuous piecewise quadratic (Lagrange P2 =
    !! P_(k+1) for the order-one Nedelec space, so that rot h approximates
    !! rot(nu curl u) one order better than the energy error) with nodal values
    !! H = sum_a h_a, where h_a lives on the vertex patch of a and minimizes
    !!
    !!   J_a(h_a) = ||nu^(-1/2)(h_a - psi_a s_h)||^2
    !!            + ||kappa^(-1/2)(rot h_a - psi_a g - s_h rot psi_a)||^2,
    !!
    !! over the P2 functions of the patch whose nodes are not on an interior
    !! edge of the patch link (so h_a extends by zero). psi_a is the hat function
    !! of a and s_h = nu curl u_h. For exact data rot(psi_a s) = psi_a rot s +
    !! s rot psi_a, so h_a = psi_a s gives J_a = 0 whenever s is in P_1: the
    !! local problems are exact for polynomial solutions. In the targets s_h is
    !! replaced by one smooth P2 recovery s (see recover_flux): the piecewise
    !! constant s_h would make (s_h - s) rot psi_a an O(1) inconsistency.
    !! Summing over a, the
    !! targets add up to s_h and g because the hats sum to one. Each J_a is a
    !! symmetric positive definite system (the L2 term alone is coercive).
    !!
    !! Pointwise integrands are emitted by FortSym (src/generated, see
    !! tools/codegen/app/gen_hcurl_bound_products.f90); topology, gather and
    !! scatter, dense patch solves and the OpenMP schedule are handwritten.
    !! Integrals use the caller's reference rule; quadrature error is not
    !! enclosed.
    use fortfem_kinds, only: dp
    use fortfem_patch_topology_2d, only: topology_t, build_topology, dgesv
    use fortfem_generated_hcurl_p2_basis, only: generated_hcurl_p2_basis
    use fortfem_generated_hcurl_rot_leaf, only: generated_hcurl_rot_leaf
    use fortfem_generated_hcurl_patch_matrix, only: generated_hcurl_patch_matrix
    use fortfem_generated_hcurl_patch_rhs, only: generated_hcurl_patch_rhs
    use fortfem_generated_hcurl_defect, only: generated_hcurl_defect
    implicit none
    private

    public :: hcurl_equilibrated_estimate_t, hcurl_equilibrated_estimate
    public :: HCURL_EQUILIBRATION_OK, HCURL_EQUILIBRATION_BAD_INPUT
    public :: HCURL_EQUILIBRATION_BAD_CELL, HCURL_EQUILIBRATION_PATCH_FAILED

    integer, parameter :: HCURL_EQUILIBRATION_OK = 0
    integer, parameter :: HCURL_EQUILIBRATION_BAD_INPUT = 1
    integer, parameter :: HCURL_EQUILIBRATION_BAD_CELL = 2
    integer, parameter :: HCURL_EQUILIBRATION_PATCH_FAILED = 3

    type :: basis_t
        real(dp), allocatable :: hat(:, :), p2(:, :), gref(:, :, :)
    end type basis_t

    type :: hcurl_equilibrated_estimate_t
        real(dp) :: bound = huge(1.0_dp)     !< sqrt(flux_defect^2 + field_defect^2)
        real(dp) :: flux_defect = 0.0_dp     !< ||nu^(-1/2)(h - nu curl u_h)||
        real(dp) :: field_defect = 0.0_dp    !< ||kappa^(-1/2)(g - rot h)||
        !> Per cell squared contributions to the two defects.
        real(dp), allocatable :: cell_flux(:), cell_field(:)
        !> Nodal values of h: vertex values, then edge midpoint values.
        real(dp), allocatable :: reconstruction(:)
    end type hcurl_equilibrated_estimate_t

contains

    subroutine hcurl_equilibrated_estimate(vertices, triangles, xi, eta, weights, &
            nu, kappa, source, discrete_field, discrete_flux, estimate, status)
        !! Data are sampled at the affine images of the reference points (xi,
        !! eta) of each cell (weights sum to 1/2). nu, kappa: (nq, nt). source
        !! f and discrete_field u_h: (2, nq, nt). discrete_flux: nu curl u_h,
        !! (nq, nt). triangles must be counter-clockwise.
        real(dp), intent(in) :: vertices(:, :)          !< (2, nv)
        integer, intent(in) :: triangles(:, :)          !< (3, nt)
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        real(dp), intent(in) :: nu(:, :), kappa(:, :)
        real(dp), intent(in) :: source(:, :, :), discrete_field(:, :, :)
        real(dp), intent(in) :: discrete_flux(:, :)
        type(hcurl_equilibrated_estimate_t), intent(out) :: estimate
        integer, intent(out) :: status
        type(topology_t) :: topo
        type(basis_t) :: basis
        real(dp), allocatable :: det(:), jac(:, :), rotc(:, :, :)
        real(dp), allocatable :: nodal(:), patch_coef(:), smooth(:, :)
        integer, allocatable :: patch_dof(:), patch_count(:), patch_offset(:)
        integer :: nv, nt, nq, a, t, q, k, failed, total, nedges

        status = HCURL_EQUILIBRATION_BAD_INPUT
        nv = size(vertices, 2)
        nt = size(triangles, 2)
        nq = size(xi)
        if (nt < 1 .or. size(vertices, 1) /= 2 .or. size(triangles, 1) /= 3) return
        if (any(triangles < 1) .or. any(triangles > nv)) return
        if (size(eta) /= nq .or. size(weights) /= nq) return
        if (any(shape(nu) /= [nq, nt]) .or. any(shape(kappa) /= [nq, nt])) return
        if (any(shape(source) /= [2, nq, nt])) return
        if (any(shape(discrete_field) /= [2, nq, nt])) return
        if (any(shape(discrete_flux) /= [nq, nt])) return
        if (any(.not. (nu > 0.0_dp)) .or. any(.not. (kappa > 0.0_dp))) return

        status = HCURL_EQUILIBRATION_BAD_CELL
        allocate (det(nt), jac(4, nt), rotc(2, 3, nt))
        call tabulate_basis(xi, eta, basis)
        do t = 1, nt
            det(t) = cell_geometry_rot(vertices, triangles(:, t), jac(:, t), &
                rotc(:, :, t))
            if (.not. (det(t) > 0.0_dp)) return
        end do

        status = HCURL_EQUILIBRATION_OK
        call build_topology(nv, triangles, topo)
        nedges = maxval(topo%cell_edge)
        ! Patch a owns the slice patch_offset(a) + 1 .. + patch_count(a) of the
        ! flat arrays, with room for every P2 dof of its patch.
        allocate (patch_offset(nv), patch_count(nv), nodal(nv + nedges))
        patch_offset(1) = 0
        do a = 2, nv
            patch_offset(a) = patch_offset(a - 1) + 4*(topo%vertex_start(a) &
                - topo%vertex_start(a - 1)) + 2
        end do
        total = patch_offset(nv) + 4*(topo%vertex_start(nv + 1) &
            - topo%vertex_start(nv)) + 2
        allocate (patch_dof(total), patch_coef(total))
        call recover_flux(vertices, triangles, topo, weights, det, basis, &
            discrete_flux, smooth)
        patch_count = 0
        failed = 0
        !$omp parallel do schedule(dynamic, 16) reduction(max:failed)
        do a = 1, nv
            if (topo%vertex_start(a + 1) == topo%vertex_start(a)) cycle
            call solve_patch(a, nv, topo, triangles, weights, det, jac, basis, rotc, &
                nu, kappa, source, discrete_field, smooth, &
                patch_dof(patch_offset(a) + 1:), patch_coef(patch_offset(a) + 1:), &
                patch_count(a), failed)
        end do
        !$omp end parallel do
        if (failed /= 0) then
            status = HCURL_EQUILIBRATION_PATCH_FAILED
            return
        end if

        ! Nodal values of h: each dof sums its coefficient over the patches
        ! that contain it (each patch contributes once).
        nodal = 0.0_dp
        do a = 1, nv
            do k = 1, patch_count(a)
                nodal(patch_dof(patch_offset(a) + k)) = &
                    nodal(patch_dof(patch_offset(a) + k)) &
                    + patch_coef(patch_offset(a) + k)
            end do
        end do
        estimate%reconstruction = nodal
        call cell_defects_all(nv, triangles, topo, weights, det, jac, basis, nu, &
            kappa, source, discrete_field, discrete_flux, nodal, estimate)
        status = HCURL_EQUILIBRATION_OK
    end subroutine hcurl_equilibrated_estimate

    subroutine recover_flux(vertices, triangles, topo, weights, det, basis, &
            discrete_flux, smooth)
        !! Continuous P2 recovery s of the discrete flux: at each vertex a
        !! least-squares quadratic fitted to the cell means of its three-ring
        !! (area weighted, exact cell averages of the monomials; falls back to
        !! a line or a constant when the ring has too few cells), evaluated at
        !! the vertex and at the midpoints of the edges it owns (an edge is
        !! owned by its smaller endpoint). The recovery is one function shared by
        !! all patch targets, so the terms s rot psi_a still cancel in the sum
        !! over a.
        !! It is O(h^3) accurate for smooth s with O(h^2) gradient error, which
        !! keeps (s_h - s) rot psi_a and rot(s - h) of higher order than the
        !! energy error (the piecewise constant s_h would leave an O(1)
        !! inconsistency, a continuous P1 recovery an O(h) one).
        real(dp), intent(in) :: vertices(:, :), weights(:), det(:)
        integer, intent(in) :: triangles(:, :)
        type(topology_t), intent(in) :: topo
        type(basis_t), intent(in) :: basis
        real(dp), intent(in) :: discrete_flux(:, :)
        real(dp), allocatable, intent(out) :: smooth(:, :)
        real(dp), allocatable :: vval(:)
        real(dp) :: m(6, 6), r(6), area, mean, bavg(6), tmp(6), hs, mm(6, 6), cf(6)
        real(dp) :: xm(2), mono(6), vv
        logical, allocatable :: edge_done(:)
        integer :: nedges, e, v1, v2, ed
        real(dp) :: p(2), wsum
        integer, allocatable :: stamp_cells(:), cell_list(:)
        real(dp) :: lo, hi
        integer :: ncl, first, last, ring
        integer, parameter :: recovery_rings = 3
        integer :: nv, nt, nq, a, i, j, k, v, t, q, ipiv(6), info, n_in, nb

        nv = size(vertices, 2)
        nt = size(triangles, 2)
        nq = size(weights)
        wsum = sum(weights)
        nedges = maxval(topo%cell_edge)
        allocate (vval(nv + nedges), smooth(nq, nt), stamp_cells(nt), cell_list(nt), &
            edge_done(nedges))
        edge_done = .false.
        stamp_cells = 0
        do a = 1, nv
            vval(a) = 0.0_dp
            if (topo%vertex_start(a + 1) == topo%vertex_start(a)) cycle
            ! Length scale: mean cell diameter proxy of the patch.
            hs = 0.0_dp
            do i = topo%vertex_start(a), topo%vertex_start(a + 1) - 1
                hs = hs + sqrt(0.5_dp*det(topo%vertex_cell(i)))
            end do
            hs = hs/real(topo%vertex_start(a + 1) - topo%vertex_start(a), dp)
            m = 0.0_dp
            r = 0.0_dp
            n_in = 0
            lo = huge(1.0_dp)
            hi = -huge(1.0_dp)
            ! Cells of the three-ring around a, collected breadth first.
            ncl = 0
            do i = topo%vertex_start(a), topo%vertex_start(a + 1) - 1
                ncl = ncl + 1
                cell_list(ncl) = topo%vertex_cell(i)
                stamp_cells(cell_list(ncl)) = a
            end do
            first = 1
            do ring = 1, recovery_rings - 1
                last = ncl
                do i = first, last
                    do j = 1, 3
                        v = triangles(j, cell_list(i))
                        do k = topo%vertex_start(v), topo%vertex_start(v + 1) - 1
                            t = topo%vertex_cell(k)
                            if (stamp_cells(t) == a) cycle
                            stamp_cells(t) = a
                            ncl = ncl + 1
                            cell_list(ncl) = t
                        end do
                    end do
                end do
                first = last + 1
            end do
            do i = 1, ncl
                t = cell_list(i)
                area = 0.5_dp*det(t)
                mean = dot_product(weights, discrete_flux(:, t))/wsum
                lo = min(lo, mean)
                hi = max(hi, mean)
                bavg = 0.0_dp
                do q = 1, nq
                    p = (matmul(vertices(:, triangles(:, t)), basis%hat(:, q)) &
                        - vertices(:, a))/hs
                    bavg = bavg + weights(q)/wsum*[1.0_dp, p(1), p(2), p(1)**2, &
                        p(1)*p(2), p(2)**2]
                end do
                m = m + area*spread(bavg, 2, 6)*spread(bavg, 1, 6)
                r = r + area*mean*bavg
            end do
            n_in = ncl
            ! Constant fallback: area-weighted mean.
            cf = 0.0_dp
            cf(1) = r(1)/m(1, 1)
            nb = 1
            if (n_in >= 3) nb = 3
            if (n_in >= 8) nb = 6
            do while (nb > 1)
                mm(1:nb, 1:nb) = m(1:nb, 1:nb)
                tmp(1:nb) = r(1:nb)
                call dgesv(nb, 1, mm, 6, ipiv, tmp, 6, info)
                if (info == 0) then
                    cf(1:nb) = tmp(1:nb)
                    cf(nb + 1:6) = 0.0_dp
                    exit
                end if
                nb = merge(3, 1, nb == 6)
            end do
            ! Guard: a one-sided boundary fit must not leave the range of the
            ! data by more than half its width (smooth data never does).
            vv = hi - lo
            lo = lo - 0.5_dp*vv
            hi = hi + 0.5_dp*vv
            vval(a) = min(max(cf(1), lo), hi)
            do i = topo%vertex_start(a), topo%vertex_start(a + 1) - 1
                t = topo%vertex_cell(i)
                do e = 1, 3
                    v1 = triangles(e, t)
                    v2 = triangles(mod(e, 3) + 1, t)
                    ed = topo%cell_edge(e, t)
                    if (a /= min(v1, v2) .or. edge_done(ed)) cycle
                    edge_done(ed) = .true.
                    xm = (0.5_dp*(vertices(:, v1) + vertices(:, v2)) &
                        - vertices(:, a))/hs
                    mono = [1.0_dp, xm(1), xm(2), xm(1)**2, xm(1)*xm(2), xm(2)**2]
                    vv = dot_product(cf, mono)
                    vval(nv + ed) = min(max(vv, lo), hi)
                end do
            end do
        end do
        do t = 1, nt
            do i = 1, nq
                smooth(i, t) = 0.0_dp
                do j = 1, 6
                    smooth(i, t) = smooth(i, t) + basis%p2(j, i) &
                        *vval(cell_dof(triangles, topo, nv, t, j))
                end do
            end do
        end do
    end subroutine recover_flux

    subroutine tabulate_basis(xi, eta, basis)
        !! P1 hats and P2 basis values and reference gradients at the rule.
        real(dp), intent(in) :: xi(:), eta(:)
        type(basis_t), intent(out) :: basis
        integer :: q, nq

        nq = size(xi)
        allocate (basis%hat(3, nq), basis%p2(6, nq), basis%gref(2, 6, nq))
        do q = 1, nq
            basis%hat(:, q) = [1.0_dp - xi(q) - eta(q), xi(q), eta(q)]
            call generated_hcurl_p2_basis(xi(q), eta(q), basis%p2(1, q), &
                basis%p2(2, q), basis%p2(3, q), basis%p2(4, q), basis%p2(5, q), &
                basis%p2(6, q), basis%gref(1, 1, q), basis%gref(1, 2, q), &
                basis%gref(1, 3, q), basis%gref(1, 4, q), basis%gref(1, 5, q), &
                basis%gref(1, 6, q), basis%gref(2, 1, q), basis%gref(2, 2, q), &
                basis%gref(2, 3, q), basis%gref(2, 4, q), basis%gref(2, 5, q), &
                basis%gref(2, 6, q))
        end do
    end subroutine tabulate_basis

    real(dp) function cell_geometry_rot(vertices, cell, jac, rotc) result(det)
        !! Determinant and Jacobian of the affine map, and the physical rot of
        !! the three P1 hat functions (constant on the cell), via the generated
        !! leaf.
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: cell(3)
        real(dp), intent(out) :: jac(4), rotc(2, 3)
        real(dp) :: ref(2, 3)
        integer :: j

        jac(1) = vertices(1, cell(2)) - vertices(1, cell(1))
        jac(2) = vertices(2, cell(2)) - vertices(2, cell(1))
        jac(3) = vertices(1, cell(3)) - vertices(1, cell(1))
        jac(4) = vertices(2, cell(3)) - vertices(2, cell(1))
        det = jac(1)*jac(4) - jac(3)*jac(2)
        ref = reshape([-1.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        do j = 1, 3
            call generated_hcurl_rot_leaf(jac(1), jac(2), jac(3), jac(4), &
                ref(1, j), ref(2, j), rotc(1, j), rotc(2, j))
        end do
    end function cell_geometry_rot

    subroutine p2_rot(jac, gref, rot)
        !! Physical rot of the six P2 basis functions at one rule point.
        real(dp), intent(in) :: jac(4), gref(2, 6)
        real(dp), intent(out) :: rot(2, 6)
        integer :: j

        do j = 1, 6
            call generated_hcurl_rot_leaf(jac(1), jac(2), jac(3), jac(4), &
                gref(1, j), gref(2, j), rot(1, j), rot(2, j))
        end do
    end subroutine p2_rot

    pure integer function cell_dof(triangles, topo, nv, t, j) result(d)
        !! Global P2 dof of local function j: vertices 1-3, then midpoints of
        !! local edges 1-3 (numbered after the vertices).
        integer, intent(in) :: triangles(:, :), nv, t, j
        type(topology_t), intent(in) :: topo

        if (j <= 3) then
            d = triangles(j, t)
        else
            d = nv + topo%cell_edge(j - 3, t)
        end if
    end function cell_dof

    subroutine solve_patch(a, nv, topo, triangles, weights, det, jac, basis, rotc, &
            nu, kappa, source, discrete_field, discrete_flux, patch_dof, &
            patch_coef, count, failed)
        !! Local SPD problem for h_a on the patch of vertex a. Stores the
        !! nonzero coefficients (free dof, value) of h_a in the patch slice.
        integer, intent(in) :: a, nv, triangles(:, :)
        type(topology_t), intent(in) :: topo
        real(dp), intent(in) :: weights(:), det(:), jac(:, :), rotc(:, :, :)
        type(basis_t), intent(in) :: basis
        real(dp), intent(in) :: nu(:, :), kappa(:, :), source(:, :, :)
        real(dp), intent(in) :: discrete_field(:, :, :), discrete_flux(:, :)
        integer, intent(inout) :: patch_dof(:), count
        real(dp), intent(inout) :: patch_coef(:)
        integer, intent(inout) :: failed
        integer, allocatable :: dlist(:), freeindex(:), ipiv(:), fi(:)
        logical, allocatable :: constrained(:)
        real(dp), allocatable :: amat(:, :), rhs(:)
        real(dp) :: g(2), c, w, nuq, kq, psi, sflux, rpsi(2), rq(2, 6)
        integer :: m, ndl, nfree, i, t, s, l1, l2, e_opp, k, j, q, info, vs, pos

        vs = topo%vertex_start(a)
        m = topo%vertex_start(a + 1) - vs
        allocate (dlist(4*m + 2), constrained(4*m + 2))
        dlist = 0
        constrained = .false.
        ndl = 0
        do i = 1, m
            t = topo%vertex_cell(vs + i - 1)
            do j = 1, 6
                call add_dof(cell_dof(triangles, topo, nv, t, j), dlist, ndl)
            end do
        end do
        do i = 1, m
            t = topo%vertex_cell(vs + i - 1)
            s = topo%vertex_slot(vs + i - 1)
            l1 = mod(s, 3) + 1
            l2 = mod(l1, 3) + 1
            e_opp = l1
            ! The apex and the spokes are never on the link; they stay free.
            if (.not. topo%boundary_edge(topo%cell_edge(e_opp, t))) then
                constrained(position_of(triangles(l1, t), dlist, ndl)) = .true.
                constrained(position_of(triangles(l2, t), dlist, ndl)) = .true.
                constrained(position_of(nv + topo%cell_edge(e_opp, t), dlist, &
                    ndl)) = .true.
            end if
        end do
        allocate (freeindex(ndl))
        nfree = 0
        do k = 1, ndl
            freeindex(k) = 0
            if (.not. constrained(k)) then
                nfree = nfree + 1
                freeindex(k) = nfree
            end if
        end do
        allocate (amat(nfree, nfree), rhs(nfree), ipiv(nfree), fi(6))
        amat = 0.0_dp
        rhs = 0.0_dp
        do i = 1, m
            t = topo%vertex_cell(vs + i - 1)
            s = topo%vertex_slot(vs + i - 1)
            do j = 1, 6
                pos = position_of(cell_dof(triangles, topo, nv, t, j), dlist, ndl)
                fi(j) = freeindex(pos)
            end do
            do q = 1, size(weights)
                call p2_rot(jac(:, t), basis%gref(:, :, q), rq)
                w = weights(q)*det(t)
                nuq = nu(q, t)
                kq = kappa(q, t)
                psi = basis%hat(s, q)
                rpsi = rotc(:, s, t)
                g = source(:, q, t) - kq*discrete_field(:, q, t)
                sflux = discrete_flux(q, t)
                do j = 1, 6
                    if (fi(j) == 0) cycle
                    call generated_hcurl_patch_rhs(w, nuq, kq, psi, sflux, g(1), &
                        g(2), rpsi(1), rpsi(2), basis%p2(j, q), rq(1, j), &
                        rq(2, j), c)
                    rhs(fi(j)) = rhs(fi(j)) + c
                    do k = 1, 6
                        if (fi(k) == 0) cycle
                        call generated_hcurl_patch_matrix(w, nuq, kq, &
                            basis%p2(j, q), basis%p2(k, q), rq(1, j), rq(2, j), &
                            rq(1, k), rq(2, k), c)
                        amat(fi(j), fi(k)) = amat(fi(j), fi(k)) + c
                    end do
                end do
            end do
        end do
        call dgesv(nfree, 1, amat, nfree, ipiv, rhs, nfree, info)
        if (info /= 0) then
            failed = 1
            return
        end if
        count = 0
        do k = 1, ndl
            if (freeindex(k) == 0) cycle
            count = count + 1
            patch_dof(count) = dlist(k)
            patch_coef(count) = rhs(freeindex(k))
        end do
    end subroutine solve_patch

    subroutine cell_defects_all(nv, triangles, topo, weights, det, jac, basis, nu, &
            kappa, source, discrete_field, discrete_flux, nodal, estimate)
        integer, intent(in) :: nv, triangles(:, :)
        type(topology_t), intent(in) :: topo
        real(dp), intent(in) :: weights(:), det(:), jac(:, :)
        type(basis_t), intent(in) :: basis
        real(dp), intent(in) :: nu(:, :), kappa(:, :), source(:, :, :)
        real(dp), intent(in) :: discrete_field(:, :, :), discrete_flux(:, :)
        real(dp), intent(in) :: nodal(:)
        type(hcurl_equilibrated_estimate_t), intent(inout) :: estimate
        real(dp) :: hq, rh(2), g(2), flux_sq, field_sq, w, hn(6), rq(2, 6)
        integer :: t, q, nt, j

        nt = size(triangles, 2)
        allocate (estimate%cell_flux(nt), estimate%cell_field(nt))
        estimate%cell_flux = 0.0_dp
        estimate%cell_field = 0.0_dp
        !$omp parallel do schedule(static) private(hq, rh, g, flux_sq, field_sq, &
        !$omp& w, hn, rq, j, q)
        do t = 1, nt
            do j = 1, 6
                hn(j) = nodal(cell_dof(triangles, topo, nv, t, j))
            end do
            do q = 1, size(weights)
                call p2_rot(jac(:, t), basis%gref(:, :, q), rq)
                hq = dot_product(hn, basis%p2(:, q))
                rh = matmul(rq, hn)
                g = source(:, q, t) - kappa(q, t)*discrete_field(:, q, t)
                w = weights(q)*det(t)
                call generated_hcurl_defect(w, nu(q, t), kappa(q, t), hq, &
                    discrete_flux(q, t), g(1), g(2), rh(1), rh(2), flux_sq, field_sq)
                estimate%cell_flux(t) = estimate%cell_flux(t) + flux_sq
                estimate%cell_field(t) = estimate%cell_field(t) + field_sq
            end do
        end do
        !$omp end parallel do
        estimate%flux_defect = sqrt(sum(estimate%cell_flux))
        estimate%field_defect = sqrt(sum(estimate%cell_field))
        estimate%bound = sqrt(estimate%flux_defect**2 + estimate%field_defect**2)
    end subroutine cell_defects_all

    pure subroutine add_dof(v, dlist, ndl)
        integer, intent(in) :: v
        integer, intent(inout) :: dlist(:), ndl

        if (position_of(v, dlist, ndl) == 0) then
            ndl = ndl + 1
            dlist(ndl) = v
        end if
    end subroutine add_dof

    pure integer function position_of(v, dlist, ndl) result(p)
        integer, intent(in) :: v, dlist(:), ndl
        integer :: k

        p = 0
        do k = 1, ndl
            if (dlist(k) == v) then
                p = k
                return
            end if
        end do
    end function position_of

end module fortfem_hcurl_equilibrated_2d
