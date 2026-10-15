module fortfem_equilibrated_flux_2d
    !! Patchwise equilibrated-flux (Prager--Synge) energy-error bound for
    !! -div(A grad u) = f on a triangle mesh, u = 0 on the whole boundary,
    !! scalar A > 0, conforming discrete u_h with u_h = 0 there. Cells are
    !! affine or isoparametric (Lagrange P2/P3 geometry, see below).
    !!
    !! Every cell K = F_K(T^) carries RT_k through the contravariant Piola map
    !! sigma = J sigma^/det J, div sigma = div^ sigma^/det J, with the Jacobian J
    !! evaluated pointwise (constant for affine cells). Edge degrees of freedom
    !! are the moments of the flux per unit reference-edge parameter, which the
    !! Piola map leaves unchanged, so the normal trace is continuous across a
    !! (possibly curved) shared edge exactly as for straight edges.
    !!
    !! For every vertex a with hat function psi_a (the reference barycentric
    !! function pulled back by F_K, a conforming partition of unity), the local
    !! problem on its patch finds sigma_a in RT_k (zero normal trace on the
    !! patch boundary away from the domain boundary; free normal trace on
    !! domain-boundary edges at Dirichlet vertices) minimizing
    !! ||A^(-1/2)(sigma_a - psi_a A grad u_h)|| subject to the divergence
    !! constraint posed in the reference cell,
    !!
    !!   div^ sigma^_a = -Pi^_k(det J (f psi_a - A grad u_h . grad psi_a)) + c_a
    !!                                                         Pi^_k(det J),
    !!
    !! with Pi^_k the L2(T^) projection onto P_k, c_a the patch mean of the load
    !! on interior patches (zero for a Galerkin solution in exact arithmetic) and
    !! c_a = 0 on boundary patches. Then sigma = sum_a sigma_a lies in H(div),
    !! div sigma + Pi f = d with Pi f := Pi^_k(det J f)/det J, and for all
    !! v in H^1_0, a(u - u_h, v) = int (f - Pi f) v + int d v
    !! + int (sigma - A grad u_h) . grad v, hence
    !!
    !!   |||u - u_h||| <= ||A^(-1/2)(sigma - A grad u_h)|| + osc + C_F ||d||,
    !!
    !! |||v|||^2 = int A |grad v|^2 and C_F a caller-supplied Friedrichs constant
    !! ||v|| <= C_F |||v||| on H^1_0. The oscillation term is, per affine
    !! (hence convex) cell, (h_T/pi) A_T,min^(-1/2) ||f - Pi f||_T (Payne--
    !! Weinberger, A_T,min a caller-supplied lower bound of A on T), and for
    !! every isoparametric (possibly non-convex) cell the unconditional
    !! C_F ||f - Pi f||_K, from int_K (f - Pi f) v <= ||f - Pi f||_K ||v||_K;
    !! the two groups are summed in quadrature separately. The bound holds for
    !! any sigma in H(div), so local solve errors only loosen it: d is measured
    !! from the assembled sigma. Integrals use the caller's reference rule;
    !! their quadrature error, and positivity of det J between the quadrature
    !! points of isoparametric cells, are not enclosed.
    !!
    !! Pointwise leaves (geometry, hat gradients, P_k monomials and the
    !! integrands) come from the FortSym-generated module
    !! fortfem_equilibrated_flux_products; the cell diameter, patch topology,
    !! gather/scatter and the dense solves stay handwritten. Degrees up to 6
    !! are supported.
    use fortfem_kinds, only: dp
    use fortfem_triangle_rt_arbitrary_order, only: triangle_rt_basis_t, &
        initialize_triangle_raviart_thomas, evaluate_triangle_raviart_thomas, &
        triangle_rt_dof_count
    use fortfem_triangle_piola_maps, only: map_triangle_rt_contravariant
    use fortfem_triangle_lagrange_arbitrary_order, only: triangle_lagrange_basis_t, &
        initialize_triangle_lagrange_basis, evaluate_triangle_lagrange_basis, &
        triangle_lagrange_nodes
    use fortfem_generated_equilibrated_jacobian_p1, only: &
        generated_equilibrated_jacobian_p1
    use fortfem_generated_equilibrated_jacobian_p2, only: &
        generated_equilibrated_jacobian_p2
    use fortfem_generated_equilibrated_jacobian_p3, only: &
        generated_equilibrated_jacobian_p3
    use fortfem_generated_equilibrated_determinant, only: &
        generated_equilibrated_determinant
    use fortfem_generated_equilibrated_inverse_transpose, only: &
        generated_equilibrated_inverse_transpose
    use fortfem_generated_equilibrated_hat_gradients, only: &
        generated_equilibrated_hat_gradients
    use fortfem_generated_equilibrated_hat_values, only: &
        generated_equilibrated_hat_values
    use fortfem_generated_equilibrated_monomials_degree0, only: &
        generated_equilibrated_monomials_degree0
    use fortfem_generated_equilibrated_monomials_degree1, only: &
        generated_equilibrated_monomials_degree1
    use fortfem_generated_equilibrated_monomials_degree2, only: &
        generated_equilibrated_monomials_degree2
    use fortfem_generated_equilibrated_monomials_degree3, only: &
        generated_equilibrated_monomials_degree3
    use fortfem_generated_equilibrated_monomials_degree4, only: &
        generated_equilibrated_monomials_degree4
    use fortfem_generated_equilibrated_monomials_degree5, only: &
        generated_equilibrated_monomials_degree5
    use fortfem_generated_equilibrated_monomials_degree6, only: &
        generated_equilibrated_monomials_degree6
    use fortfem_generated_equilibrated_load_integrand, only: &
        generated_equilibrated_load_integrand
    use fortfem_generated_equilibrated_flux_defect_integrand, only: &
        generated_equilibrated_flux_defect_integrand
    use fortfem_generated_equilibrated_oscillation_integrand, only: &
        generated_equilibrated_oscillation_integrand
    use fortfem_generated_equilibrated_divergence_integrand, only: &
        generated_equilibrated_divergence_integrand
    implicit none
    private

    public :: equilibrated_flux_estimate_t, equilibrated_flux_estimate

    integer, parameter, public :: EQUILIBRATION_OK = 0
    integer, parameter, public :: EQUILIBRATION_BAD_INPUT = 1
    integer, parameter, public :: EQUILIBRATION_BAD_CELL = 2
    integer, parameter, public :: EQUILIBRATION_PATCH_FAILED = 3
    integer, parameter :: max_degree = 6

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
        integer :: degree = 0, nrt = 0, nk = 0, nq = 0, nedge = 0, nloc = 3
        real(dp), allocatable :: weights(:)
        real(dp), allocatable :: rt(:, :, :), rt_div(:, :)  !< (2,nrt,nq),(nrt,nq)
        real(dp), allocatable :: poly(:, :)                 !< (nk, nq), P_k
        real(dp), allocatable :: hat(:, :)                  !< (3, nq)
        !> Reference gradients of the geometry shape functions, (nloc, 2, nq).
        real(dp), allocatable :: geo_grad(:, :, :)
    end type reference_t

    !> Cell geometry: Lagrange nodes (vertices first) and an affinity flag.
    type :: problem_t
        real(dp), allocatable :: nodes(:, :, :)  !< (2, nloc, nt)
        logical, allocatable :: affine(:)        !< (nt)
    end type problem_t

    !> Mesh topology: vertex -> cells, cell edges, boundary vertices.
    type :: topology_t
        integer, allocatable :: vertex_start(:), vertex_cell(:), vertex_slot(:)
        integer, allocatable :: cell_edge(:, :), cell_orientation(:, :)
        logical, allocatable :: boundary_vertex(:)
    end type topology_t

contains

    subroutine equilibrated_flux_estimate(vertices, triangles, degree, xi, eta, &
            weights, coefficient, source, discrete_flux, coefficient_lower, &
            friedrichs_constant, estimate, status, geometry_nodes)
        !! Data arrays are sampled at the (mapped) reference points (xi, eta)
        !! of each cell (weights sum to 1/2); the rule must integrate
        !! polynomials of degree 2k + 2 exactly for the local problems to be
        !! well posed. triangles must be counter-clockwise.
        !! Optional geometry_nodes(2, nloc, nt), nloc = 3, 6 or 10, gives the
        !! Lagrange geometry nodes of each cell: the vertices
        !! vertices(:, triangles(:, t)) first, then the nloc = 6 edge midpoints
        !! or the nloc = 10 two nodes per edge (ordered from local vertex e to
        !! e + 1, edges 1-2, 2-3, 3-1) followed by the centroid node. Nodes
        !! shared by neighbouring cells must coincide. A cell whose nodes lie on
        !! the straight-edge Lagrange positions is treated as affine.
        real(dp), intent(in) :: vertices(:, :)          !< (2, nv)
        integer, intent(in) :: triangles(:, :)          !< (3, nt)
        integer, intent(in) :: degree                   !< 0 <= k <= 6
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        real(dp), intent(in) :: coefficient(:, :)       !< A, (nq, nt)
        real(dp), intent(in) :: source(:, :)            !< f, (nq, nt)
        real(dp), intent(in) :: discrete_flux(:, :, :)  !< A grad u_h, (2, nq, nt)
        real(dp), intent(in) :: coefficient_lower(:)    !< (nt)
        real(dp), intent(in) :: friedrichs_constant
        type(equilibrated_flux_estimate_t), intent(out) :: estimate
        integer, intent(out) :: status
        real(dp), intent(in), optional :: geometry_nodes(:, :, :)
        type(reference_t) :: ref
        type(topology_t) :: topo
        type(problem_t) :: prob
        real(dp), allocatable :: projection(:, :), contribution(:, :, :)
        real(dp), allocatable :: divergence_sq(:), source_sq(:)
        integer :: nv, nt, nloc, a, t, failed, nodes_status
        logical :: cells_ok

        status = EQUILIBRATION_BAD_INPUT
        nv = size(vertices, 2)
        nt = size(triangles, 2)
        if (degree < 0 .or. degree > max_degree) return
        if (nt < 1 .or. size(vertices, 1) /= 2) return
        if (size(triangles, 1) /= 3) return
        if (any(triangles < 1) .or. any(triangles > nv)) return
        if (size(eta) /= size(xi) .or. size(weights) /= size(xi)) return
        if (size(coefficient, 1) /= size(xi) .or. size(coefficient, 2) /= nt) return
        if (size(source, 1) /= size(xi) .or. size(source, 2) /= nt) return
        if (size(discrete_flux, 1) /= 2) return
        if (size(discrete_flux, 2) /= size(xi)) return
        if (size(discrete_flux, 3) /= nt) return
        if (size(coefficient_lower) /= nt) return
        if (any(.not. (coefficient_lower > 0.0_dp))) return
        if (any(.not. (coefficient > 0.0_dp))) return
        if (.not. (friedrichs_constant >= 0.0_dp)) return

        nloc = 3
        if (present(geometry_nodes)) nloc = size(geometry_nodes, 2)
        if (nloc /= 3 .and. nloc /= 6 .and. nloc /= 10) return
        call setup_geometry(vertices, triangles, geometry_nodes, nloc, prob, &
            nodes_status)
        if (nodes_status /= 0) return

        call build_reference(degree, nloc, xi, eta, weights, ref, status)
        if (status /= 0) then
            status = EQUILIBRATION_BAD_INPUT
            return
        end if
        status = EQUILIBRATION_BAD_CELL
        cells_ok = .true.
        do t = 1, nt
            if (.not. positive_determinant(ref, prob%nodes(:, :, t))) cells_ok = .false.
        end do
        if (.not. cells_ok) return
        call build_topology(nv, triangles, topo)

        allocate (projection(ref%nk, nt), estimate%cell_oscillation(nt))
        allocate (contribution(ref%nrt, 3, nt), source_sq(nt))
        contribution = 0.0_dp
        failed = 0
        !$omp parallel do schedule(static)
        do t = 1, nt
            call project_source(ref, prob%nodes(:, :, t), source(:, t), &
                projection(:, t), source_sq(t))
            if (prob%affine(t)) then
                estimate%cell_oscillation(t) = cell_diameter(prob%nodes(:, 1:3, t)) &
                    /acos(-1.0_dp)*sqrt(source_sq(t)/coefficient_lower(t))
            else
                estimate%cell_oscillation(t) = friedrichs_constant*sqrt(source_sq(t))
            end if
        end do
        !$omp end parallel do

        !$omp parallel do schedule(dynamic, 16) reduction(max:failed)
        do a = 1, nv
            if (topo%vertex_start(a + 1) == topo%vertex_start(a)) cycle
            call solve_patch(ref, topo, prob, a, triangles, coefficient, &
                source, discrete_flux, contribution, failed)
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
            call cell_defects(ref, prob%nodes(:, :, t), coefficient(:, t), &
                discrete_flux(:, :, t), projection(:, t), &
                estimate%flux_coefficients(:, t), estimate%cell_flux(t), &
                divergence_sq(t))
        end do
        !$omp end parallel do
        estimate%flux_defect = sqrt(sum(estimate%cell_flux**2))
        estimate%oscillation = sqrt(sum(estimate%cell_oscillation**2, &
            mask=prob%affine)) + sqrt(sum(estimate%cell_oscillation**2, &
            mask=.not. prob%affine))
        estimate%divergence_defect = friedrichs_constant*sqrt(sum(divergence_sq))
        estimate%bound = estimate%flux_defect + estimate%oscillation &
            + estimate%divergence_defect
        status = EQUILIBRATION_OK
    end subroutine equilibrated_flux_estimate

    subroutine setup_geometry(vertices, triangles, geometry_nodes, nloc, prob, &
            status)
        !! Lagrange nodes of every cell (affine nodes from the vertices when
        !! geometry_nodes is absent) and the affinity flags.
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), nloc
        real(dp), intent(in), optional :: geometry_nodes(:, :, :)
        type(problem_t), intent(out) :: prob
        integer, intent(out) :: status
        real(dp), allocatable :: straight(:, :)
        real(dp) :: reference(2, nloc), diameter
        integer :: nt, t, i

        status = 1
        nt = size(triangles, 2)
        allocate (prob%nodes(2, nloc, nt), prob%affine(nt))
        call reference_geometry_nodes(nloc, reference)
        if (present(geometry_nodes)) then
            if (size(geometry_nodes, 1) /= 2 .or. size(geometry_nodes, 3) /= nt) return
            prob%nodes = geometry_nodes
        else
            do t = 1, nt
                prob%nodes(:, :, t) = vertices(:, triangles(:, t))
            end do
        end if
        allocate (straight(2, nloc))
        do t = 1, nt
            if (any(abs(prob%nodes(:, 1:3, t) - vertices(:, triangles(:, t))) &
                > 0.0_dp)) return
            do i = 1, nloc
                straight(:, i) = prob%nodes(:, 1, t) &
                    + reference(1, i)*(prob%nodes(:, 2, t) - prob%nodes(:, 1, t)) &
                    + reference(2, i)*(prob%nodes(:, 3, t) - prob%nodes(:, 1, t))
            end do
            diameter = cell_diameter(prob%nodes(:, 1:3, t))
            prob%affine(t) = all(abs(prob%nodes(:, :, t) - straight) &
                <= 1.0e-13_dp*diameter)
        end do
        status = 0
    end subroutine setup_geometry

    pure subroutine reference_geometry_nodes(nloc, nodes)
        !! Reference positions of the Lagrange geometry nodes: vertices, then
        !! per edge e (local vertex e to e + 1) its interior nodes, then the
        !! centroid for P3.
        integer, intent(in) :: nloc
        real(dp), intent(out) :: nodes(2, nloc)
        real(dp), parameter :: third = 1.0_dp/3.0_dp

        nodes(:, 1:3) = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], &
            [2, 3])
        if (nloc == 6) then
            nodes(:, 4:6) = reshape([0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, &
                0.5_dp], [2, 3])
        else if (nloc == 10) then
            nodes(:, 4:9) = reshape([third, 0.0_dp, 2*third, 0.0_dp, &
                2*third, third, third, 2*third, 0.0_dp, 2*third, 0.0_dp, third], &
                [2, 6])
            nodes(:, 10) = [third, third]
        end if
    end subroutine reference_geometry_nodes

    logical function positive_determinant(ref, nodes) result(ok)
        !! det J > 0 at every quadrature point of the cell.
        type(reference_t), intent(in) :: ref
        real(dp), intent(in) :: nodes(:, :)
        real(dp) :: jac(2, 2), det
        integer :: q

        ok = .true.
        do q = 1, ref%nq
            call jacobian_at(ref, nodes, q, jac, det)
            if (.not. (det > 0.0_dp)) ok = .false.
        end do
    end function positive_determinant

    subroutine build_reference(degree, nloc, xi, eta, weights, ref, status)
        integer, intent(in) :: degree, nloc
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        type(reference_t), intent(out) :: ref
        integer, intent(out) :: status
        type(triangle_rt_basis_t) :: basis
        type(triangle_lagrange_basis_t) :: geometry_basis
        real(dp), allocatable :: lagrange_nodes(:, :), lagrange_value(:)
        real(dp), allocatable :: lagrange_grad(:, :)
        real(dp) :: targets(2, nloc)
        integer :: q, i, j, order(nloc), geometry_degree

        call initialize_triangle_raviart_thomas(degree, basis, status)
        if (status /= 0) return
        ref%degree = degree
        ref%nrt = triangle_rt_dof_count(basis)
        ref%nk = (degree + 1)*(degree + 2)/2
        ref%nq = size(xi)
        ref%nedge = degree + 1
        ref%nloc = nloc
        ref%weights = weights
        allocate (ref%rt(2, ref%nrt, ref%nq), ref%rt_div(ref%nrt, ref%nq))
        allocate (ref%poly(ref%nk, ref%nq), ref%hat(3, ref%nq))
        allocate (ref%geo_grad(nloc, 2, ref%nq))
        ! Geometry shape functions in the node order of the module header: the
        ! Lagrange basis uses its own ordering, matched here by node position.
        geometry_degree = merge(1, merge(2, 3, nloc == 6), nloc == 3)
        call initialize_triangle_lagrange_basis(geometry_degree, geometry_basis, &
            status)
        if (status /= 0) return
        call triangle_lagrange_nodes(geometry_basis, lagrange_nodes)
        call reference_geometry_nodes(nloc, targets)
        do i = 1, nloc
            order(i) = 0
            do j = 1, nloc
                if (maxval(abs(lagrange_nodes(:, j) - targets(:, i))) < 1.0e-12_dp) &
                    order(i) = j
            end do
            if (order(i) == 0) then
                status = 1
                return
            end if
        end do
        allocate (lagrange_value(nloc), lagrange_grad(2, nloc))
        do q = 1, ref%nq
            call evaluate_triangle_raviart_thomas(basis, xi(q), eta(q), &
                ref%rt(:, :, q), ref%rt_div(:, q), status)
            if (status /= 0) return
            call reference_monomials(degree, xi(q), eta(q), ref%poly(:, q))
            call generated_equilibrated_hat_values(xi(q), eta(q), ref%hat(:, q))
            call evaluate_triangle_lagrange_basis(geometry_basis, xi(q), eta(q), &
                lagrange_value, lagrange_grad, status)
            if (status /= 0) return
            do i = 1, nloc
                ref%geo_grad(i, :, q) = lagrange_grad(:, order(i))
            end do
        end do
        status = 0
    end subroutine build_reference

    pure subroutine reference_monomials(degree, xi, eta, values)
        !! xi**i eta**j, i = 0..degree, j = 0..degree-i, in that nesting order.
        integer, intent(in) :: degree
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(:)

        select case (degree)
        case (0)
            call generated_equilibrated_monomials_degree0(xi, eta, values)
        case (1)
            call generated_equilibrated_monomials_degree1(xi, eta, values)
        case (2)
            call generated_equilibrated_monomials_degree2(xi, eta, values)
        case (3)
            call generated_equilibrated_monomials_degree3(xi, eta, values)
        case (4)
            call generated_equilibrated_monomials_degree4(xi, eta, values)
        case (5)
            call generated_equilibrated_monomials_degree5(xi, eta, values)
        case default
            call generated_equilibrated_monomials_degree6(xi, eta, values)
        end select
    end subroutine reference_monomials

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

    pure subroutine jacobian_at(ref, nodes, q, jac, det)
        !! Jacobian of the cell map at quadrature point q and its determinant.
        type(reference_t), intent(in) :: ref
        real(dp), intent(in) :: nodes(:, :)
        integer, intent(in) :: q
        real(dp), intent(out) :: jac(2, 2), det

        select case (ref%nloc)
        case (3)
            call generated_equilibrated_jacobian_p1(nodes, ref%geo_grad(:, :, q), jac)
        case (6)
            call generated_equilibrated_jacobian_p2(nodes, ref%geo_grad(:, :, q), jac)
        case default
            call generated_equilibrated_jacobian_p3(nodes, ref%geo_grad(:, :, q), jac)
        end select
        call generated_equilibrated_determinant(jac, det)
    end subroutine jacobian_at

    pure subroutine geometry_at(ref, nodes, q, jac, det, hat_grad)
        !! Pointwise Jacobian, determinant and physical hat gradients.
        type(reference_t), intent(in) :: ref
        real(dp), intent(in) :: nodes(:, :)
        integer, intent(in) :: q
        real(dp), intent(out) :: jac(2, 2), det, hat_grad(2, 3)
        real(dp) :: inv_t(2, 2)

        call jacobian_at(ref, nodes, q, jac, det)
        ! J^(-T) maps reference gradients to physical ones.
        call generated_equilibrated_inverse_transpose(jac, det, inv_t)
        call generated_equilibrated_hat_gradients(inv_t, hat_grad)
    end subroutine geometry_at

    pure real(dp) function cell_diameter(vertex_nodes) result(diameter)
        !! Longest vertex distance of the (straight-edge) cell.
        real(dp), intent(in) :: vertex_nodes(:, :)

        diameter = max(norm2(vertex_nodes(:, 2) - vertex_nodes(:, 1)), &
            norm2(vertex_nodes(:, 3) - vertex_nodes(:, 1)), &
            norm2(vertex_nodes(:, 3) - vertex_nodes(:, 2)))
    end function cell_diameter

    subroutine physical_rt(ref, jac, q, values, divergences, status)
        type(reference_t), intent(in) :: ref
        real(dp), intent(in) :: jac(2, 2)
        integer, intent(in) :: q
        real(dp), intent(out) :: values(:, :), divergences(:)
        integer, intent(out) :: status

        call map_triangle_rt_contravariant(jac, ref%rt(:, :, q), ref%rt_div(:, q), &
            values, divergences, status)
    end subroutine physical_rt

    subroutine project_source(ref, nodes, f, coefficients, source_sq)
        !! Reference-cell L2 projection of det J f onto P_k(T^), so that
        !! Pi f = (coefficients . P_k)/det J, and ||f - Pi f||_K^2.
        type(reference_t), intent(in) :: ref
        real(dp), intent(in) :: nodes(:, :), f(:)
        real(dp), intent(out) :: coefficients(:), source_sq
        real(dp) :: jac(2, 2), det, mass(ref%nk, ref%nk), rhs(ref%nk, 1)
        real(dp) :: integrand
        integer :: q, l, ipiv(ref%nk), info

        mass = 0.0_dp
        rhs = 0.0_dp
        do q = 1, ref%nq
            call jacobian_at(ref, nodes, q, jac, det)
            do l = 1, ref%nk
                mass(:, l) = mass(:, l) + ref%weights(q)*ref%poly(:, q)*ref%poly(l, q)
            end do
            rhs(:, 1) = rhs(:, 1) + ref%weights(q)*det*f(q)*ref%poly(:, q)
        end do
        call dgesv(ref%nk, 1, mass, ref%nk, ipiv, rhs, ref%nk, info)
        coefficients = rhs(:, 1)
        source_sq = 0.0_dp
        do q = 1, ref%nq
            call jacobian_at(ref, nodes, q, jac, det)
            call generated_equilibrated_oscillation_integrand(f(q), &
                dot_product(coefficients, ref%poly(:, q)), det, integrand)
            source_sq = source_sq + ref%weights(q)*integrand
        end do
        if (info /= 0) source_sq = huge(1.0_dp)
    end subroutine project_source

    subroutine solve_patch(ref, topo, prob, a, triangles, coefficient, source, &
            discrete_flux, contribution, failed)
        type(reference_t), intent(in) :: ref
        type(topology_t), intent(in) :: topo
        type(problem_t), intent(in) :: prob
        integer, intent(in) :: a, triangles(:, :)
        real(dp), intent(in) :: coefficient(:, :), source(:, :)
        real(dp), intent(in) :: discrete_flux(:, :, :)
        real(dp), intent(inout) :: contribution(:, :, :)
        integer, intent(inout) :: failed
        integer :: edges(2*(topo%vertex_start(a + 1) - topo%vertex_start(a)))
        integer, allocatable :: map(:, :), sgn(:, :), ipiv(:)
        real(dp), allocatable :: kkt(:, :), rhs(:, :), qint(:)
        real(dp) :: jac(2, 2), det, hat_grad(2, 3), wdet, ainv, psi
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

            lam0 = nsigma + (i - 1)*ref%nk
            do q = 1, ref%nq
                call geometry_at(ref, prob%nodes(:, :, t), q, jac, det, hat_grad)
                call physical_rt(ref, jac, q, phi, div, stat)
                if (stat /= 0) then
                    failed = 1
                    return
                end if
                wdet = ref%weights(q)*det
                ainv = 1.0_dp/coefficient(q, t)
                psi = ref%hat(s, q)
                call generated_equilibrated_load_integrand(source(q, t), psi, &
                    discrete_flux(:, q, t), hat_grad(:, s), load)
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

    subroutine cell_defects(ref, nodes, coefficient, discrete_flux, projection, &
            sigma, flux_defect, divergence_sq)
        type(reference_t), intent(in) :: ref
        real(dp), intent(in) :: nodes(:, :), coefficient(:)
        real(dp), intent(in) :: discrete_flux(:, :), projection(:), sigma(:)
        real(dp), intent(out) :: flux_defect, divergence_sq
        real(dp) :: jac(2, 2), det, phi(2, ref%nrt), div(ref%nrt), value(2)
        real(dp) :: integrand
        integer :: q, stat

        flux_defect = 0.0_dp
        divergence_sq = 0.0_dp
        do q = 1, ref%nq
            call jacobian_at(ref, nodes, q, jac, det)
            call physical_rt(ref, jac, q, phi, div, stat)
            value = matmul(phi, sigma) - discrete_flux(:, q)
            call generated_equilibrated_flux_defect_integrand(value, coefficient(q), &
                det, integrand)
            flux_defect = flux_defect + ref%weights(q)*integrand
            call generated_equilibrated_divergence_integrand(dot_product(div, sigma), &
                dot_product(projection, ref%poly(:, q)), det, integrand)
            divergence_sq = divergence_sq + ref%weights(q)*integrand
        end do
        flux_defect = sqrt(flux_defect)
    end subroutine cell_defects

end module fortfem_equilibrated_flux_2d
