module fortfem_equilibrated_flux_2d
    !! Patchwise equilibrated-flux (Prager--Synge) energy-error bound for
    !! -div(A grad u) = f on an affine triangle mesh, u = 0 on the whole
    !! boundary, scalar A > 0, conforming discrete u_h with u_h = 0 there.
    !!
    !! For every vertex a with hat function psi_a, the local problem on its
    !! patch finds sigma_a in RT_k (zero normal trace on the patch boundary
    !! away from the domain boundary) minimizing ||A^(-1/2)(sigma_a -
    !! psi_a A grad u_h)|| subject to -div sigma_a = Pi_k(f psi_a - A grad u_h
    !! . grad psi_a) - c_a, with c_a the patch mean of that load on interior
    !! patches (zero for a Galerkin solution in exact arithmetic) and c_a = 0
    !! on boundary patches. sigma = sum_a sigma_a lies in H(div) and, for all
    !! v in H^1_0, a(u - u_h, v) = int (f - Pi_k f) v - int d v
    !! + int (sigma - A grad u_h) . grad v, d = div sigma + Pi_k f, hence
    !!
    !!   |||u - u_h||| <= ||A^(-1/2)(sigma - A grad u_h)||
    !!       + (sum_T (h_T/pi)^2 / A_T,min ||f - Pi_k f||_T^2)^(1/2)
    !!       + C_F ||d||,
    !!
    !! |||v|||^2 = int A |grad v|^2, h_T the diameter (Payne--Weinberger on
    !! the convex cell), A_T,min a caller-supplied lower bound of A on T and
    !! C_F a caller-supplied Friedrichs constant ||v|| <= C_F |||v||| on
    !! H^1_0. The bound holds for any sigma in H(div), so local solve errors
    !! only loosen it: d is measured from the assembled sigma. Integrals use
    !! the caller's reference rule; their quadrature error is not enclosed.
    use fortfem_kinds, only: dp
    use fortfem_triangle_rt_arbitrary_order, only: triangle_rt_basis_t, &
        initialize_triangle_raviart_thomas, evaluate_triangle_raviart_thomas, &
        triangle_rt_dof_count
    use fortfem_triangle_piola_maps, only: map_triangle_rt_contravariant
    implicit none
    private

    public :: equilibrated_flux_estimate_t, equilibrated_flux_estimate
    ! Shared with the non-conforming estimators (equilibrated_nonconforming_2d).
    public :: reference_t, build_reference, project_source

    integer, parameter, public :: EQUILIBRATION_OK = 0
    integer, parameter, public :: EQUILIBRATION_BAD_INPUT = 1
    integer, parameter, public :: EQUILIBRATION_BAD_CELL = 2
    integer, parameter, public :: EQUILIBRATION_PATCH_FAILED = 3

    type :: equilibrated_flux_estimate_t
        real(dp) :: bound = huge(1.0_dp)        !< sum of the three terms
        real(dp) :: flux_defect = 0.0_dp        !< ||A^(-1/2)(sigma - A grad u_h)||
        real(dp) :: oscillation = 0.0_dp        !< data-oscillation term
        real(dp) :: divergence_defect = 0.0_dp  !< C_F ||div sigma + Pi_k f||
        !> Per cell: flux-defect and oscillation contributions (not squared).
        real(dp), allocatable :: cell_flux(:), cell_oscillation(:)
        !> Element-local physical RT_k coefficients of sigma, (dofs, cells).
        real(dp), allocatable :: flux_coefficients(:, :)
    end type equilibrated_flux_estimate_t

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(*), info
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
        end subroutine dgesv
    end interface

    !> Reference data shared by all cells.
    type :: reference_t
        integer :: degree = 0, nrt = 0, nk = 0, nq = 0, nedge = 0
        real(dp), allocatable :: weights(:)
        real(dp), allocatable :: rt(:, :, :), rt_div(:, :)  !< (2,nrt,nq),(nrt,nq)
        real(dp), allocatable :: poly(:, :)                 !< (nk, nq), P_k
        real(dp), allocatable :: hat(:, :)                  !< (3, nq)
    end type reference_t

    !> Mesh topology: vertex -> cells, cell edges, boundary vertices.
    type :: topology_t
        integer, allocatable :: vertex_start(:), vertex_cell(:), vertex_slot(:)
        integer, allocatable :: cell_edge(:, :), cell_orientation(:, :)
        logical, allocatable :: boundary_vertex(:)
    end type topology_t

contains

    subroutine equilibrated_flux_estimate(vertices, triangles, degree, xi, eta, &
            weights, coefficient, source, discrete_flux, coefficient_lower, &
            friedrichs_constant, estimate, status, load_flux)
        !! Data arrays are sampled at the affine images of the reference
        !! points (xi, eta) of each cell (weights sum to 1/2); the rule must
        !! integrate polynomials of degree 2k + 2 exactly for the local
        !! problems to be well posed. triangles must be counter-clockwise.
        real(dp), intent(in) :: vertices(:, :)          !< (2, nv)
        integer, intent(in) :: triangles(:, :)          !< (3, nt)
        integer, intent(in) :: degree                   !< k >= 0
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        real(dp), intent(in) :: coefficient(:, :)       !< A, (nq, nt)
        real(dp), intent(in) :: source(:, :)            !< f, (nq, nt)
        real(dp), intent(in) :: discrete_flux(:, :, :)  !< A grad u_h, (2, nq, nt)
        real(dp), intent(in) :: coefficient_lower(:)    !< (nt)
        real(dp), intent(in) :: friedrichs_constant
        type(equilibrated_flux_estimate_t), intent(out) :: estimate
        integer, intent(out) :: status
        !> Optional flux (2, nq, nt) used instead of discrete_flux in the
        !> divergence constraints of the patch problems, i.e. in
        !> Pi_k(f psi_a - load_flux . grad psi_a) (DG residual lifting). The
        !> minimised target stays psi_a * discrete_flux.
        real(dp), intent(in), optional :: load_flux(:, :, :)
        type(reference_t) :: ref
        type(topology_t) :: topo
        real(dp), allocatable :: projection(:, :), contribution(:, :, :)
        real(dp), allocatable :: divergence_sq(:)
        integer :: nv, nt, a, t, failed

        status = EQUILIBRATION_BAD_INPUT
        nv = size(vertices, 2)
        nt = size(triangles, 2)
        if (degree < 0 .or. nt < 1 .or. size(vertices, 1) /= 2) return
        if (size(triangles, 1) /= 3) return
        if (any(triangles < 1) .or. any(triangles > nv)) return
        if (size(eta) /= size(xi) .or. size(weights) /= size(xi)) return
        if (size(coefficient, 1) /= size(xi) .or. size(coefficient, 2) /= nt) return
        if (size(source, 1) /= size(xi) .or. size(source, 2) /= nt) return
        if (size(discrete_flux, 1) /= 2) return
        if (size(discrete_flux, 2) /= size(xi)) return
        if (size(discrete_flux, 3) /= nt) return
        if (size(coefficient_lower) /= nt) return
        if (present(load_flux)) then
            if (size(load_flux, 1) /= 2) return
            if (size(load_flux, 2) /= size(xi)) return
            if (size(load_flux, 3) /= nt) return
        end if
        if (any(.not. (coefficient_lower > 0.0_dp))) return
        if (any(.not. (coefficient > 0.0_dp))) return
        if (.not. (friedrichs_constant >= 0.0_dp)) return

        call build_reference(degree, xi, eta, weights, ref, status)
        if (status /= 0) then
            status = EQUILIBRATION_BAD_INPUT
            return
        end if
        status = EQUILIBRATION_BAD_CELL
        do t = 1, nt
            if (.not. (cell_determinant(vertices, triangles(:, t)) > 0.0_dp)) return
        end do
        call build_topology(nv, triangles, topo)

        allocate (projection(ref%nk, nt), estimate%cell_oscillation(nt))
        allocate (contribution(ref%nrt, 3, nt))
        contribution = 0.0_dp
        failed = 0
        !$omp parallel do schedule(static)
        do t = 1, nt
            call project_source(ref, vertices, triangles(:, t), source(:, t), &
                coefficient_lower(t), projection(:, t), &
                estimate%cell_oscillation(t))
        end do
        !$omp end parallel do

        !$omp parallel do schedule(dynamic, 16) reduction(max:failed)
        do a = 1, nv
            if (topo%vertex_start(a + 1) == topo%vertex_start(a)) cycle
            if (present(load_flux)) then
                call solve_patch(ref, topo, a, vertices, triangles, coefficient, &
                    source, discrete_flux, load_flux, contribution, failed)
            else
                call solve_patch(ref, topo, a, vertices, triangles, coefficient, &
                    source, discrete_flux, discrete_flux, contribution, failed)
            end if
        end do
        !$omp end parallel do
        if (failed /= 0) then
            status = EQUILIBRATION_PATCH_FAILED
            return
        end if

        allocate (estimate%flux_coefficients(ref%nrt, nt), estimate%cell_flux(nt))
        allocate (divergence_sq(nt))
        !$omp parallel do schedule(static)
        do t = 1, nt
            estimate%flux_coefficients(:, t) = sum(contribution(:, :, t), dim=2)
            call cell_defects(ref, vertices, triangles(:, t), coefficient(:, t), &
                discrete_flux(:, :, t), projection(:, t), &
                estimate%flux_coefficients(:, t), estimate%cell_flux(t), &
                divergence_sq(t))
        end do
        !$omp end parallel do
        estimate%flux_defect = sqrt(sum(estimate%cell_flux**2))
        estimate%oscillation = sqrt(sum(estimate%cell_oscillation**2))
        estimate%divergence_defect = friedrichs_constant*sqrt(sum(divergence_sq))
        estimate%bound = estimate%flux_defect + estimate%oscillation &
            + estimate%divergence_defect
        status = EQUILIBRATION_OK
    end subroutine equilibrated_flux_estimate

    subroutine build_reference(degree, xi, eta, weights, ref, status)
        integer, intent(in) :: degree
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        type(reference_t), intent(out) :: ref
        integer, intent(out) :: status
        type(triangle_rt_basis_t) :: basis
        integer :: q, i, j, l

        call initialize_triangle_raviart_thomas(degree, basis, status)
        if (status /= 0) return
        ref%degree = degree
        ref%nrt = triangle_rt_dof_count(basis)
        ref%nk = (degree + 1)*(degree + 2)/2
        ref%nq = size(xi)
        ref%nedge = degree + 1
        ref%weights = weights
        allocate (ref%rt(2, ref%nrt, ref%nq), ref%rt_div(ref%nrt, ref%nq))
        allocate (ref%poly(ref%nk, ref%nq), ref%hat(3, ref%nq))
        do q = 1, ref%nq
            call evaluate_triangle_raviart_thomas(basis, xi(q), eta(q), &
                ref%rt(:, :, q), ref%rt_div(:, q), status)
            if (status /= 0) return
            l = 0
            do i = 0, degree
                do j = 0, degree - i
                    l = l + 1
                    ref%poly(l, q) = xi(q)**i*eta(q)**j
                end do
            end do
            ref%hat(:, q) = [1.0_dp - xi(q) - eta(q), xi(q), eta(q)]
        end do
        status = 0
    end subroutine build_reference

    subroutine build_topology(nv, triangles, topo)
        !! Edges are numbered by first appearance; local edge e of a cell runs
        !! from its local vertex e to vertex mod(e, 3) + 1 (the FortFEM RT
        !! reference-edge convention) and has orientation +1 when that start
        !! vertex has the smaller global index.
        integer, intent(in) :: nv, triangles(:, :)
        type(topology_t), intent(out) :: topo
        integer, allocatable :: fill(:)
        integer :: nt, t, s, e, a, b, k, other, oe, nedges

        nt = size(triangles, 2)
        allocate (topo%vertex_start(nv + 1), fill(nv))
        topo%vertex_start = 0
        do t = 1, nt
            do s = 1, 3
                a = triangles(s, t)
                topo%vertex_start(a + 1) = topo%vertex_start(a + 1) + 1
            end do
        end do
        topo%vertex_start(1) = 1
        do a = 1, nv
            topo%vertex_start(a + 1) = topo%vertex_start(a + 1) + topo%vertex_start(a)
        end do
        allocate (topo%vertex_cell(3*nt), topo%vertex_slot(3*nt))
        fill = topo%vertex_start(1:nv)
        do t = 1, nt
            do s = 1, 3
                a = triangles(s, t)
                topo%vertex_cell(fill(a)) = t
                topo%vertex_slot(fill(a)) = s
                fill(a) = fill(a) + 1
            end do
        end do

        allocate (topo%cell_edge(3, nt), topo%cell_orientation(3, nt))
        allocate (topo%boundary_vertex(nv))
        topo%cell_edge = 0
        topo%boundary_vertex = .false.
        nedges = 0
        do t = 1, nt
            do e = 1, 3
                a = triangles(e, t)
                b = triangles(mod(e, 3) + 1, t)
                topo%cell_orientation(e, t) = merge(1, -1, a < b)
                other = 0
                oe = 0
                do k = topo%vertex_start(a), topo%vertex_start(a + 1) - 1
                    if (topo%vertex_cell(k) == t) cycle
                    oe = shared_edge(triangles(:, topo%vertex_cell(k)), a, b)
                    if (oe /= 0) then
                        other = topo%vertex_cell(k)
                        exit
                    end if
                end do
                if (other == 0) then
                    topo%boundary_vertex(a) = .true.
                    topo%boundary_vertex(b) = .true.
                end if
                if (other /= 0 .and. other < t) then
                    topo%cell_edge(e, t) = topo%cell_edge(oe, other)
                else
                    nedges = nedges + 1
                    topo%cell_edge(e, t) = nedges
                end if
            end do
        end do
    end subroutine build_topology

    pure integer function shared_edge(cell, a, b) result(edge)
        !! Local edge of `cell` joining global vertices a and b, or 0.
        integer, intent(in) :: cell(3), a, b
        integer :: e, p, q

        edge = 0
        do e = 1, 3
            p = cell(e)
            q = cell(mod(e, 3) + 1)
            if ((p == a .and. q == b) .or. (p == b .and. q == a)) then
                edge = e
                return
            end if
        end do
    end function shared_edge

    pure real(dp) function cell_determinant(vertices, cell) result(det)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: cell(3)
        real(dp) :: jac(2, 2)

        jac(:, 1) = vertices(:, cell(2)) - vertices(:, cell(1))
        jac(:, 2) = vertices(:, cell(3)) - vertices(:, cell(1))
        det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
    end function cell_determinant

    pure subroutine cell_geometry(vertices, cell, jac, det, hat_grad, diameter)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: cell(3)
        real(dp), intent(out) :: jac(2, 2), det, hat_grad(2, 3), diameter
        real(dp) :: inv_t(2, 2)

        jac(:, 1) = vertices(:, cell(2)) - vertices(:, cell(1))
        jac(:, 2) = vertices(:, cell(3)) - vertices(:, cell(1))
        det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
        ! J^(-T) maps reference gradients to physical ones.
        inv_t(1, :) = [jac(2, 2), -jac(2, 1)]/det
        inv_t(2, :) = [-jac(1, 2), jac(1, 1)]/det
        hat_grad(:, 1) = matmul(inv_t, [-1.0_dp, -1.0_dp])
        hat_grad(:, 2) = matmul(inv_t, [1.0_dp, 0.0_dp])
        hat_grad(:, 3) = matmul(inv_t, [0.0_dp, 1.0_dp])
        diameter = max(norm2(jac(:, 1)), norm2(jac(:, 2)), &
            norm2(jac(:, 2) - jac(:, 1)))
    end subroutine cell_geometry

    subroutine physical_rt(ref, jac, q, values, divergences, status)
        type(reference_t), intent(in) :: ref
        real(dp), intent(in) :: jac(2, 2)
        integer, intent(in) :: q
        real(dp), intent(out) :: values(:, :), divergences(:)
        integer, intent(out) :: status

        call map_triangle_rt_contravariant(jac, ref%rt(:, :, q), ref%rt_div(:, q), &
            values, divergences, status)
    end subroutine physical_rt

    subroutine project_source(ref, vertices, cell, f, a_lower, coefficients, &
            oscillation)
        !! L2 projection of f onto P_k(T) and the oscillation term
        !! (h_T/pi) A_T,min^(-1/2) ||f - Pi_k f||_T.
        type(reference_t), intent(in) :: ref
        real(dp), intent(in) :: vertices(:, :), f(:), a_lower
        integer, intent(in) :: cell(3)
        real(dp), intent(out) :: coefficients(:), oscillation
        real(dp) :: jac(2, 2), det, hat_grad(2, 3), diameter, wdet, rest
        real(dp) :: mass(ref%nk, ref%nk), rhs(ref%nk, 1)
        integer :: q, l, ipiv(ref%nk), info

        call cell_geometry(vertices, cell, jac, det, hat_grad, diameter)
        mass = 0.0_dp
        rhs = 0.0_dp
        do q = 1, ref%nq
            wdet = ref%weights(q)*det
            do l = 1, ref%nk
                mass(:, l) = mass(:, l) + wdet*ref%poly(:, q)*ref%poly(l, q)
            end do
            rhs(:, 1) = rhs(:, 1) + wdet*f(q)*ref%poly(:, q)
        end do
        call dgesv(ref%nk, 1, mass, ref%nk, ipiv, rhs, ref%nk, info)
        coefficients = rhs(:, 1)
        rest = 0.0_dp
        do q = 1, ref%nq
            rest = rest + ref%weights(q)*det &
                *(f(q) - dot_product(coefficients, ref%poly(:, q)))**2
        end do
        oscillation = diameter/acos(-1.0_dp)*sqrt(rest/a_lower)
        if (info /= 0) oscillation = huge(1.0_dp)
    end subroutine project_source

    subroutine solve_patch(ref, topo, a, vertices, triangles, coefficient, source, &
            discrete_flux, load_flux, contribution, failed)
        type(reference_t), intent(in) :: ref
        type(topology_t), intent(in) :: topo
        integer, intent(in) :: a, triangles(:, :)
        real(dp), intent(in) :: vertices(:, :), coefficient(:, :), source(:, :)
        real(dp), intent(in) :: discrete_flux(:, :, :), load_flux(:, :, :)
        real(dp), intent(inout) :: contribution(:, :, :)
        integer, intent(inout) :: failed
        integer :: edges(2*(topo%vertex_start(a + 1) - topo%vertex_start(a)))
        integer, allocatable :: map(:, :), sgn(:, :), ipiv(:)
        real(dp), allocatable :: kkt(:, :), rhs(:, :), qint(:)
        real(dp) :: jac(2, 2), det, hat_grad(2, 3), diameter, wdet, ainv, psi
        real(dp) :: phi(2, ref%nrt), div(ref%nrt), load, mean, area
        integer :: m, nedges, nsigma, nlam, n, i, j, k, l, t, s, e, ee, moment
        integer :: q, cell_dofs, interior, slot_edge(2), info, stat, lam0
        logical :: interior_patch

        m = topo%vertex_start(a + 1) - topo%vertex_start(a)
        interior_patch = .not. topo%boundary_vertex(a)
        cell_dofs = ref%degree*(ref%degree + 1)
        ! Edges of the patch incident to a: local edges s and s - 1 of each cell.
        nedges = 0
        do i = 1, m
            t = topo%vertex_cell(topo%vertex_start(a) + i - 1)
            s = topo%vertex_slot(topo%vertex_start(a) + i - 1)
            do ee = 1, 2
                e = merge(s, mod(s + 1, 3) + 1, ee == 1)
                if (.not. any(edges(:nedges) == topo%cell_edge(e, t))) then
                    nedges = nedges + 1
                    edges(nedges) = topo%cell_edge(e, t)
                end if
            end do
        end do
        nsigma = nedges*ref%nedge + m*cell_dofs
        nlam = m*ref%nk
        n = nsigma + nlam + merge(1, 0, interior_patch)
        allocate (map(ref%nrt, m), sgn(ref%nrt, m), kkt(n, n), rhs(n, 1))
        allocate (ipiv(n), qint(nlam))
        map = 0
        sgn = 1
        kkt = 0.0_dp
        rhs = 0.0_dp
        qint = 0.0_dp
        do i = 1, m
            t = topo%vertex_cell(topo%vertex_start(a) + i - 1)
            s = topo%vertex_slot(topo%vertex_start(a) + i - 1)
            slot_edge = [s, mod(s + 1, 3) + 1]
            do ee = 1, 2
                e = slot_edge(ee)
                j = findloc(edges(:nedges), topo%cell_edge(e, t), dim=1)
                do moment = 1, ref%nedge
                    k = (e - 1)*ref%nedge + moment
                    map(k, i) = (j - 1)*ref%nedge + moment
                    if (topo%cell_orientation(e, t) == -1 .and. &
                        mod(moment, 2) == 1) sgn(k, i) = -1
                end do
            end do
            do interior = 1, cell_dofs
                map(3*ref%nedge + interior, i) = nedges*ref%nedge &
                    + (i - 1)*cell_dofs + interior
            end do

            call cell_geometry(vertices, triangles(:, t), jac, det, hat_grad, &
                diameter)
            lam0 = nsigma + (i - 1)*ref%nk
            do q = 1, ref%nq
                call physical_rt(ref, jac, q, phi, div, stat)
                if (stat /= 0) then
                    failed = 1
                    return
                end if
                wdet = ref%weights(q)*det
                ainv = 1.0_dp/coefficient(q, t)
                psi = ref%hat(s, q)
                load = source(q, t)*psi &
                    - dot_product(load_flux(:, q, t), hat_grad(:, s))
                do k = 1, ref%nrt
                    if (map(k, i) == 0) cycle
                    rhs(map(k, i), 1) = rhs(map(k, i), 1) + wdet*sgn(k, i)*psi &
                        *ainv*dot_product(discrete_flux(:, q, t), phi(:, k))
                    do j = 1, ref%nrt
                        if (map(j, i) == 0) cycle
                        kkt(map(j, i), map(k, i)) = kkt(map(j, i), map(k, i)) &
                            + wdet*sgn(j, i)*sgn(k, i)*ainv &
                            *dot_product(phi(:, j), phi(:, k))
                    end do
                    do l = 1, ref%nk
                        kkt(lam0 + l, map(k, i)) = kkt(lam0 + l, map(k, i)) &
                            - wdet*sgn(k, i)*div(k)*ref%poly(l, q)
                    end do
                end do
                do l = 1, ref%nk
                    rhs(lam0 + l, 1) = rhs(lam0 + l, 1) + wdet*load*ref%poly(l, q)
                    qint((i - 1)*ref%nk + l) = qint((i - 1)*ref%nk + l) &
                        + wdet*ref%poly(l, q)
                end do
            end do
        end do
        ! Symmetric saddle point: copy the constraint block above the diagonal.
        kkt(1:nsigma, nsigma + 1:nsigma + nlam) = &
            transpose(kkt(nsigma + 1:nsigma + nlam, 1:nsigma))
        if (interior_patch) then
            ! Remove the patch mean of the load (the Galerkin residual at psi_a)
            ! and fix the constant multiplier with one extra mean constraint.
            area = 0.0_dp
            mean = 0.0_dp
            do i = 1, m
                area = area + qint((i - 1)*ref%nk + 1)
                mean = mean + rhs(nsigma + (i - 1)*ref%nk + 1, 1)
            end do
            mean = mean/area
            rhs(nsigma + 1:nsigma + nlam, 1) = rhs(nsigma + 1:nsigma + nlam, 1) &
                - mean*qint
            kkt(nsigma + 1:nsigma + nlam, n) = qint
            kkt(n, nsigma + 1:nsigma + nlam) = qint
        end if
        call dgesv(n, 1, kkt, n, ipiv, rhs, n, info)
        if (info /= 0) then
            failed = 1
            return
        end if
        do i = 1, m
            t = topo%vertex_cell(topo%vertex_start(a) + i - 1)
            s = topo%vertex_slot(topo%vertex_start(a) + i - 1)
            do k = 1, ref%nrt
                if (map(k, i) == 0) then
                    contribution(k, s, t) = 0.0_dp
                else
                    contribution(k, s, t) = sgn(k, i)*rhs(map(k, i), 1)
                end if
            end do
        end do
    end subroutine solve_patch

    subroutine cell_defects(ref, vertices, cell, coefficient, discrete_flux, &
            projection, sigma, flux_defect, divergence_sq)
        type(reference_t), intent(in) :: ref
        real(dp), intent(in) :: vertices(:, :), coefficient(:)
        integer, intent(in) :: cell(3)
        real(dp), intent(in) :: discrete_flux(:, :), projection(:), sigma(:)
        real(dp), intent(out) :: flux_defect, divergence_sq
        real(dp) :: jac(2, 2), det, hat_grad(2, 3), diameter, wdet
        real(dp) :: phi(2, ref%nrt), div(ref%nrt), value(2)
        integer :: q, stat

        call cell_geometry(vertices, cell, jac, det, hat_grad, diameter)
        flux_defect = 0.0_dp
        divergence_sq = 0.0_dp
        do q = 1, ref%nq
            call physical_rt(ref, jac, q, phi, div, stat)
            wdet = ref%weights(q)*det
            value = matmul(phi, sigma) - discrete_flux(:, q)
            flux_defect = flux_defect + wdet*dot_product(value, value)/coefficient(q)
            divergence_sq = divergence_sq + wdet*(dot_product(div, sigma) &
                + dot_product(projection, ref%poly(:, q)))**2
        end do
        flux_defect = sqrt(flux_defect)
    end subroutine cell_defects

end module fortfem_equilibrated_flux_2d
