module fortfem_equilibrated_flux_2d
    !! Patchwise equilibrated-flux (Prager--Synge) energy-error bound for
    !! -div(A grad u) = f on a triangle mesh, u = 0 on the Dirichlet part of the
    !! boundary and A grad u . n = g on an optional Neumann part, scalar A > 0,
    !! conforming discrete u_h with u_h = 0 on the Dirichlet part. Cells are
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
    !! Dirichlet boundary edges; on Neumann edges the trace is prescribed as the
    !! moments of psi_a g, see below) minimizing
    !! ||A^(-1/2)(sigma_a - psi_a A grad u_h)|| subject to the divergence
    !! constraint posed in the reference cell,
    !!
    !!   div^ sigma^_a = -Pi^_k(det J (f psi_a - A grad u_h . grad psi_a)) + c_a
    !!                                                         Pi^_k(det J),
    !!
    !! with Pi^_k the L2(T^) projection onto P_k and c_a = (int det J load_a +
    !! int_GammaN psi_a g)/|omega_a| on patches touching no Dirichlet edge
    !! (interior vertices and vertices inside the Neumann part; zero for a
    !! Galerkin solution in exact arithmetic), c_a = 0 elsewhere. On a Neumann
    !! edge the normal flux of sigma_a is the L2(0,1) projection onto P_k of
    !! psi_a g |x'(s)| (edge parameter s), i.e. its Legendre moments are the
    !! prescribed edge degrees of freedom. Then sigma = sum_a sigma_a lies in
    !! H(div), div sigma + Pi f = d with Pi f := Pi^_k(det J f)/det J, and for
    !! all v in H^1 vanishing on the Dirichlet part,
    !! a(u - u_h, v) = int (f - Pi f) v + int d v + int (sigma - A grad u_h) . grad v
    !! + int_GammaN (g - sigma.n) v, hence
    !!
    !!   |||u - u_h||| <= ||A^(-1/2)(sigma - A grad u_h)|| + osc + C_F ||d||
    !!                    + neumann,
    !!
    !! |||v|||^2 = int A |grad v|^2 and C_F a caller-supplied Friedrichs constant
    !! ||v|| <= C_F |||v||| on that space. The oscillation term is, per affine
    !! (hence convex) cell, (h_T/pi) A_T,min^(-1/2) ||f - Pi f||_T (Payne--
    !! Weinberger, A_T,min a caller-supplied lower bound of A on T), and for
    !! every isoparametric (possibly non-convex) cell the unconditional
    !! C_F ||f - Pi f||_K, from int_K (f - Pi f) v <= ||f - Pi f||_K ||v||_K;
    !! the two groups are summed in quadrature separately. The Neumann term is
    !! sqrt(sum_T (C_T sum_e ||g - sigma.n||_e)^2) over the Neumann edges e of
    !! the (affine) cell T, with the guaranteed trace constant
    !! C_T = h_T sqrt((1 + 3/pi^2)/(r_T A_T,min)), r_T the inradius (derivation
    !! in neumann_defect). The bound holds for
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
    use fortfem_generated_equilibrated_edge_legendre, only: &
        generated_equilibrated_edge_legendre
    use fortfem_generated_equilibrated_edge_speed, only: &
        generated_equilibrated_edge_speed
    use fortfem_generated_equilibrated_reference_normal, only: &
        generated_equilibrated_reference_normal
    use fortfem_generated_equilibrated_flux_normal_component, only: &
        generated_equilibrated_flux_normal_component
    use fortfem_generated_equilibrated_neumann_target, only: &
        generated_equilibrated_neumann_target
    use fortfem_generated_equilibrated_neumann_defect_integrand, only: &
        generated_equilibrated_neumann_defect_integrand
    use fortfem_equilibrated_flux_topology, only: topology_t, build_topology, &
        dgesv, EQUILIBRATION_OK, EQUILIBRATION_BAD_INPUT, EQUILIBRATION_BAD_CELL, &
        EQUILIBRATION_PATCH_FAILED
    implicit none
    private

    public :: equilibrated_flux_estimate_t, equilibrated_flux_estimate
    ! Shared with the non-conforming estimators (equilibrated_nonconforming_2d).
    public :: reference_t, build_reference, project_source

    public :: EQUILIBRATION_OK, EQUILIBRATION_BAD_INPUT, EQUILIBRATION_BAD_CELL, &
        EQUILIBRATION_PATCH_FAILED
    integer, parameter :: max_degree = 6
    !> Reference gradients of the P1 geometry functions, grad(node, component).
    real(dp), parameter :: p1_gradient(3, 2) = reshape( &
        [-1.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp], [3, 2])
    !> Counter-clockwise reference edge tangents; edge e runs from local vertex
    !> e to e + 1 and is parametrized by s in [0, 1].
    real(dp), parameter :: edge_tangent(2, 3) = reshape( &
        [1.0_dp, 0.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, -1.0_dp], [2, 3])

    type :: equilibrated_flux_estimate_t
        real(dp) :: bound = huge(1.0_dp)        !< sum of the three terms
        real(dp) :: flux_defect = 0.0_dp        !< ||A^(-1/2)(sigma - A grad u_h)||
        real(dp) :: oscillation = 0.0_dp        !< data-oscillation term
        real(dp) :: divergence_defect = 0.0_dp  !< C_F ||div sigma + Pi_k f||
        !> Neumann term sqrt(sum_T (C_T sum_e ||g - sigma.n||_e)^2), zero without
        !> Neumann edges.
        real(dp) :: neumann_oscillation = 0.0_dp
        !> Per cell: flux-defect and oscillation contributions (not squared).
        real(dp), allocatable :: cell_flux(:), cell_oscillation(:)
        !> Element-local physical RT_k coefficients of sigma, (dofs, cells).
        real(dp), allocatable :: flux_coefficients(:, :)
        !> Per cell: C_T sum_e ||g - sigma.n||_e over its Neumann edges.
        real(dp), allocatable :: cell_neumann(:)
    end type equilibrated_flux_estimate_t

    !> Reference data shared by all cells.
    type :: reference_t
        integer :: degree = 0, nrt = 0, nk = 0, nq = 0, nedge = 0, nloc = 3
        integer :: nq1 = 0
        real(dp), allocatable :: weights(:)
        real(dp), allocatable :: rt(:, :, :), rt_div(:, :)  !< (2,nrt,nq),(nrt,nq)
        real(dp), allocatable :: poly(:, :)                 !< (nk, nq), P_k
        real(dp), allocatable :: hat(:, :)                  !< (3, nq)
        !> Reference gradients of the geometry shape functions, (nloc, 2, nq).
        real(dp), allocatable :: geo_grad(:, :, :)
        !> Edge rule on [0, 1], per reference edge: hat values (3, nq1), RT normal
        !> trace sigma^.n of every RT function (nrt, nq1), shifted Legendre
        !> values (nedge, nq1).
        real(dp), allocatable :: edge_weights(:), edge_hat(:, :, :)
        real(dp), allocatable :: edge_flux(:, :, :), edge_legendre(:, :)
    end type reference_t

    !> Cell geometry: Lagrange nodes (vertices first) and an affinity flag.
    type :: problem_t
        real(dp), allocatable :: nodes(:, :, :)  !< (2, nloc, nt)
        logical, allocatable :: affine(:)        !< (nt)
        logical, allocatable :: neumann(:, :)    !< (3, nt) flagged edges
        real(dp), allocatable :: flux(:, :, :)   !< g at the edge points (nq1, 3, nt)
    end type problem_t

contains

    subroutine equilibrated_flux_estimate(vertices, triangles, degree, xi, eta, &
            weights, coefficient, source, discrete_flux, coefficient_lower, &
            friedrichs_constant, estimate, status, geometry_nodes, neumann_edge, &
            neumann_flux, edge_xi, edge_weights, load_flux)
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
        !! Optional Neumann data (all four together): neumann_edge(3, nt) flags
        !! local edge e of cell t (from local vertex e to e + 1) as a boundary
        !! edge carrying the prescribed outward flux g = A grad u . n,
        !! neumann_flux(nq1, 3, nt) samples g at the points edge_xi(nq1) in
        !! [0, 1] of that edge (quadrature weights edge_weights, summing to 1;
        !! exact for polynomials of degree k + 1 + deg g). Flagged edges must be
        !! boundary edges of affine cells; the other boundary edges are Dirichlet,
        !! and friedrichs_constant must then bound ||v|| <= C_F |||v||| for the
        !! H^1 functions vanishing on the Dirichlet edges only.
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
        logical, intent(in), optional :: neumann_edge(:, :)
        real(dp), intent(in), optional :: neumann_flux(:, :, :)
        real(dp), intent(in), optional :: edge_xi(:), edge_weights(:)
        !> Optional flux (2, nq, nt) used instead of discrete_flux in the
        !> divergence constraints of the patch problems, i.e. in
        !> Pi_k(f psi_a - load_flux . grad psi_a) (DG residual lifting). The
        !> minimised target stays psi_a * discrete_flux.
        real(dp), intent(in), optional :: load_flux(:, :, :)
        type(reference_t) :: ref
        type(topology_t) :: topo
        type(problem_t) :: prob
        real(dp), allocatable :: projection(:, :), contribution(:, :, :)
        real(dp), allocatable :: divergence_sq(:), source_sq(:)
        real(dp), allocatable :: edge_points(:), edge_rule(:)
        integer :: nv, nt, nloc, a, t, failed, nodes_status, nq1
        logical :: cells_ok, has_neumann

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
        if (present(load_flux)) then
            if (size(load_flux, 1) /= 2) return
            if (size(load_flux, 2) /= size(xi)) return
            if (size(load_flux, 3) /= nt) return
        end if
        if (any(.not. (coefficient_lower > 0.0_dp))) return
        if (any(.not. (coefficient > 0.0_dp))) return
        if (.not. (friedrichs_constant >= 0.0_dp)) return

        nloc = 3
        if (present(geometry_nodes)) nloc = size(geometry_nodes, 2)
        if (nloc /= 3 .and. nloc /= 6 .and. nloc /= 10) return
        call setup_geometry(vertices, triangles, geometry_nodes, nloc, prob, &
            nodes_status)
        if (nodes_status /= 0) return
        has_neumann = present(neumann_edge)
        if (present(neumann_flux) .neqv. has_neumann) return
        if (present(edge_xi) .neqv. has_neumann) return
        if (present(edge_weights) .neqv. has_neumann) return
        allocate (prob%neumann(3, nt))
        prob%neumann = .false.
        nq1 = 0
        if (has_neumann) then
            nq1 = size(edge_xi)
            if (nq1 < 1 .or. size(edge_weights) /= nq1) return
            if (size(neumann_edge, 1) /= 3 .or. size(neumann_edge, 2) /= nt) return
            if (size(neumann_flux, 1) /= nq1) return
            if (size(neumann_flux, 2) /= 3 .or. size(neumann_flux, 3) /= nt) return
            prob%neumann = neumann_edge
            do t = 1, nt
                if (any(prob%neumann(:, t)) .and. .not. prob%affine(t)) return
            end do
            allocate (edge_points(nq1), edge_rule(nq1))
            edge_points = edge_xi
            edge_rule = edge_weights
            prob%flux = neumann_flux
        else
            allocate (edge_points(0), edge_rule(0), prob%flux(0, 3, nt))
        end if

        call build_reference(degree, nloc, xi, eta, weights, edge_points, &
            edge_rule, ref, status)
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
        call build_topology(nv, triangles, topo, nodes_status, prob%neumann)
        if (nodes_status /= 0) then
            status = EQUILIBRATION_BAD_INPUT
            return
        end if

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
            if (present(load_flux)) then
                call solve_patch(ref, topo, prob, a, coefficient, source, &
                    discrete_flux, load_flux, contribution, failed)
            else
                call solve_patch(ref, topo, prob, a, coefficient, source, &
                    discrete_flux, discrete_flux, contribution, failed)
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
        allocate (estimate%cell_neumann(nt))
        estimate%cell_neumann = 0.0_dp
        if (has_neumann) then
            !$omp parallel do schedule(static)
            do t = 1, nt
                if (.not. any(prob%neumann(:, t))) cycle
                call neumann_defect(ref, prob, t, estimate%flux_coefficients(:, t), &
                    coefficient_lower(t), estimate%cell_neumann(t))
            end do
            !$omp end parallel do
        end if
        estimate%neumann_oscillation = sqrt(sum(estimate%cell_neumann**2))
        estimate%bound = estimate%flux_defect + estimate%oscillation &
            + estimate%divergence_defect + estimate%neumann_oscillation
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

    subroutine build_reference(degree, nloc, xi, eta, weights, edge_xi, &
            edge_weights, ref, status)
        integer, intent(in) :: degree, nloc
        real(dp), intent(in) :: xi(:), eta(:), weights(:), edge_xi(:)
        real(dp), intent(in) :: edge_weights(:)
        type(reference_t), intent(out) :: ref
        integer, intent(out) :: status
        type(triangle_rt_basis_t) :: basis
        type(triangle_lagrange_basis_t) :: geometry_basis
        real(dp), allocatable :: lagrange_nodes(:, :), lagrange_value(:)
        real(dp), allocatable :: lagrange_grad(:, :)
        real(dp) :: targets(2, nloc), point(2), legendre(7), normal(2)
        real(dp), allocatable :: rt_edge(:, :), div_edge(:)
        integer :: q, i, j, e, order(nloc), geometry_degree

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
        ref%nq1 = size(edge_xi)
        ref%edge_weights = edge_weights
        allocate (ref%edge_hat(3, ref%nq1, 3), ref%edge_flux(ref%nrt, ref%nq1, 3))
        allocate (ref%edge_legendre(ref%nedge, ref%nq1))
        allocate (rt_edge(2, ref%nrt), div_edge(ref%nrt))
        do e = 1, 3
            call generated_equilibrated_reference_normal(edge_tangent(:, e), normal)
            do q = 1, ref%nq1
                select case (e)
                case (1)
                    point = [edge_xi(q), 0.0_dp]
                case (2)
                    point = [1.0_dp - edge_xi(q), edge_xi(q)]
                case default
                    point = [0.0_dp, 1.0_dp - edge_xi(q)]
                end select
                call generated_equilibrated_hat_values(point(1), point(2), &
                    ref%edge_hat(:, q, e))
                call evaluate_triangle_raviart_thomas(basis, point(1), point(2), &
                    rt_edge, div_edge, status)
                if (status /= 0) return
                do i = 1, ref%nrt
                    call generated_equilibrated_flux_normal_component(rt_edge(:, i), &
                        normal, ref%edge_flux(i, q, e))
                end do
                call generated_equilibrated_edge_legendre(edge_xi(q), legendre)
                ref%edge_legendre(:, q) = legendre(1:ref%nedge)
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

    subroutine solve_patch(ref, topo, prob, a, coefficient, source, &
            discrete_flux, load_flux, contribution, failed)
        type(reference_t), intent(in) :: ref
        type(topology_t), intent(in) :: topo
        type(problem_t), intent(in) :: prob
        integer, intent(in) :: a
        real(dp), intent(in) :: coefficient(:, :), source(:, :)
        real(dp), intent(in) :: discrete_flux(:, :, :), load_flux(:, :, :)
        real(dp), intent(inout) :: contribution(:, :, :)
        integer, intent(inout) :: failed
        integer :: edges(2*(topo%vertex_start(a + 1) - topo%vertex_start(a)))
        integer, allocatable :: map(:, :), sgn(:, :), ipiv(:)
        real(dp), allocatable :: kkt(:, :), rhs(:, :), qint(:)
        real(dp) :: jac(2, 2), det, hat_grad(2, 3), wdet, ainv, psi
        real(dp) :: phi(2, ref%nrt), div(ref%nrt), load, mean, area
        real(dp) :: moments(ref%nedge), neumann_total
        real(dp), allocatable :: fixed_value(:)
        logical, allocatable :: fixed(:)
        integer :: m, nedges, nsigma, nlam, n, i, j, k, l, t, s, e, ee, moment
        integer :: q, cell_dofs, interior, slot_edge(2), info, stat, lam0
        logical :: needs_mean

        m = topo%vertex_start(a + 1) - topo%vertex_start(a)
        ! Patches touching no Dirichlet edge (interior vertices and vertices
        ! inside a Neumann part) need the compatibility shift c_a.
        needs_mean = .not. topo%dirichlet_vertex(a)
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
        n = nsigma + nlam + merge(1, 0, needs_mean)
        allocate (map(ref%nrt, m), sgn(ref%nrt, m), kkt(n, n), rhs(n, 1))
        allocate (ipiv(n), qint(nlam), fixed(n), fixed_value(n))
        fixed = .false.
        fixed_value = 0.0_dp
        neumann_total = 0.0_dp
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
                if (prob%neumann(e, t)) then
                    ! Prescribed trace: the moments of psi_a g, in the patch
                    ! orientation (local = sgn * patch).
                    call neumann_moments(ref, prob, t, e, s, moments)
                    neumann_total = neumann_total + moments(1)
                    do moment = 1, ref%nedge
                        k = (e - 1)*ref%nedge + moment
                        fixed(map(k, i)) = .true.
                        fixed_value(map(k, i)) = sgn(k, i)*moments(moment)
                    end do
                end if
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
                    load_flux(:, q, t), hat_grad(:, s), load)
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
        if (needs_mean) then
            ! Remove the patch mean of the load plus the prescribed boundary flux
            ! (the Galerkin residual at psi_a) and fix the constant multiplier
            ! with one extra mean constraint.
            area = 0.0_dp
            mean = neumann_total
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
        do j = 1, nsigma
            if (.not. fixed(j)) cycle
            rhs(:, 1) = rhs(:, 1) - kkt(:, j)*fixed_value(j)
            kkt(j, :) = 0.0_dp
            kkt(:, j) = 0.0_dp
            kkt(j, j) = 1.0_dp
            rhs(j, 1) = fixed_value(j)
        end do
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

    subroutine neumann_moments(ref, prob, t, e, slot, moments)
        !! Moments int_0^1 psi(s) g(s) |x'(s)| P_j(s) ds, j = 0..k, of the flux
        !! per unit parameter on local edge e of cell t, psi the hat function of
        !! the cell's local vertex `slot`: the edge degrees of freedom of the
        !! L2 projection of psi g onto P_k.
        type(reference_t), intent(in) :: ref
        type(problem_t), intent(in) :: prob
        integer, intent(in) :: t, e, slot
        real(dp), intent(out) :: moments(:)
        real(dp) :: jac(2, 2), speed, target_value
        integer :: i

        call generated_equilibrated_jacobian_p1(prob%nodes(:, 1:3, t), p1_gradient, &
            jac)
        call generated_equilibrated_edge_speed(jac, edge_tangent(:, e), speed)
        moments = 0.0_dp
        do i = 1, ref%nq1
            call generated_equilibrated_neumann_target(ref%edge_hat(slot, i, e), &
                prob%flux(i, e, t), speed, target_value)
            moments = moments + ref%edge_weights(i)*target_value*ref%edge_legendre(:, i)
        end do
    end subroutine neumann_moments

    subroutine neumann_defect(ref, prob, t, sigma, a_lower, value)
        !! C_T sum_e ||g - sigma.n||_e over the flagged edges of the affine cell
        !! t. For rho = g - sigma.n orthogonal to constants on the edge and any
        !! constant c, int_e rho v = int_e rho (v - c) <= ||rho||_e ||v - c||_e.
        !! With c the cell mean of v, Payne--Weinberger gives ||v - c||_T <= (h_T/pi)
        !! ||grad v||_T, and the divergence theorem for w^2 (x - x0) with the
        !! incenter x0 (so (x - x0).n = r_T on every edge, |x - x0| <= h_T) gives
        !! r_T ||w||^2_dT <= 3 ||w||_T^2 + h_T^2 ||grad w||_T^2 for w = v - c, hence
        !!   ||v - c||_e^2 <= h_T^2 (1 + 3/pi^2) / (r_T A_T,min) |||v|||_T^2.
        !! (Cf. Carstensen--Funken 2000, Veeser--Verfuerth 2012 for such traces.)
        type(reference_t), intent(in) :: ref
        type(problem_t), intent(in) :: prob
        integer, intent(in) :: t
        real(dp), intent(in) :: sigma(:), a_lower
        real(dp), intent(out) :: value
        real(dp) :: jac(2, 2), det, speed, integrand, flux_value, defect, rho
        real(dp) :: diameter, inradius, perimeter, constant
        real(dp), parameter :: pi = acos(-1.0_dp)
        integer :: e, i

        call generated_equilibrated_jacobian_p1(prob%nodes(:, 1:3, t), p1_gradient, &
            jac)
        call generated_equilibrated_determinant(jac, det)
        diameter = cell_diameter(prob%nodes(:, 1:3, t))
        perimeter = norm2(prob%nodes(:, 2, t) - prob%nodes(:, 1, t)) &
            + norm2(prob%nodes(:, 3, t) - prob%nodes(:, 2, t)) &
            + norm2(prob%nodes(:, 1, t) - prob%nodes(:, 3, t))
        inradius = det/perimeter
        constant = diameter*sqrt((1.0_dp + 3.0_dp/pi**2)/(inradius*a_lower))
        rho = 0.0_dp
        do e = 1, 3
            if (.not. prob%neumann(e, t)) cycle
            call generated_equilibrated_edge_speed(jac, edge_tangent(:, e), speed)
            defect = 0.0_dp
            do i = 1, ref%nq1
                flux_value = dot_product(sigma, ref%edge_flux(:, i, e))
                call generated_equilibrated_neumann_defect_integrand( &
                    prob%flux(i, e, t), flux_value, speed, integrand)
                defect = defect + ref%edge_weights(i)*integrand
            end do
            rho = rho + sqrt(defect)
        end do
        value = constant*rho
    end subroutine neumann_defect

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
