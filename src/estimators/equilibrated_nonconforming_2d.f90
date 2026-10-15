module fortfem_equilibrated_nonconforming_2d
    !! Equilibrated-residual (Prager--Synge) bounds for the non-conforming
    !! diffusion solvers of FortFEM on affine triangles, homogeneous Dirichlet
    !! data, scalar diffusion A > 0:  -div(A grad u) = f.
    !!
    !! Common identity. For a flux sigma in H(div), any s in H^1_0, w = u - s,
    !! d = div sigma + Pi_k f (Pi_k the cellwise L2 projection onto P_k):
    !!   (A grad w, grad w) = (f - Pi_k f, w) + (d, w) + (sigma - A grad s, grad w)
    !! and |(f - Pi_k f, w) + (d, w)| <= E |||w|||, |||w|||^2 = int A |grad w|^2,
    !! E = osc + C_F ||d||, osc^2 = sum_T (h_T/pi)^2 / A_T,min ||f - Pi_k f||_T^2
    !! (Payne--Weinberger), C_F the Friedrichs constant of |||.|||. All
    !! quantities are measured on the discrete sigma, so a local solve error
    !! only loosens a bound; quadrature error of the integrals (the caller's
    !! rule) is not enclosed.
    !!
    !! 1. Mixed RT_k-DG_k (mixed_poisson_flux_bound), A = 1, sigma_h = -grad u
    !!    in RT_k with div sigma_h = Pi_k f (up to quadrature). Reconstruction
    !!    s: per cell the P_(k+1) function minimising int |sigma_h + grad s|^2
    !!    with mean s = mean p_h (Arbogast--Chen / Stenberg postprocessing),
    !!    made H^1_0-conforming by Oswald nodal averaging; local, O(N), no
    !!    global solve. With X = ||sigma_h + grad s|| (the flux defect):
    !!      ||sigma - sigma_h|| <= sqrt(X^2 + E^2)         (flux error)
    !!      |||u - s|||          <= X + sqrt(X^2 + E^2)    (potential error)
    !!    (error^2 = X^2 + 2 (f - div sigma_h, w) - |||w|||^2 <= X^2 + E^2).
    !!
    !! 2. SIPG (sipg_energy_bound): flux sigma from the vertex-patch RT
    !!    equilibration of the broken flux A grad_h u_h with the DG residual
    !!    lifting in the divergence constraints: with the cellwise lifting
    !!    l_T in P_k(T)^2, int_T A l_T . tau = sum_(F in dT) w_F int_F A_T [u_h]
    !!    (tau . n_T), w_F = 1/2 interior, 1 boundary, the SIPG equation tested
    !!    with the hat function psi_a gives int (f psi_a - A (grad_h u_h - l)
    !!    . grad psi_a) = 0, so the interior patch problems are compatible and
    !!    d ~ 0, while sigma stays the minimiser of ||A^(-1/2)(sigma - psi_a A
    !!    grad_h u_h)||. With s the Oswald average of u_h (H^1_0, same degree),
    !!    rho = ||A^(-1/2)(sigma - A grad_h u_h)||, nu = ||A^(1/2) grad_h(u_h-s)||,
    !!    e = u - u_h and x = ||A^(1/2) grad_h e||:
    !!      x^2 - (E + rho + nu) x - nu (E + rho) <= 0,
    !!    hence x <= (b + sqrt(b^2 + 4c))/2, b = E + rho + nu, c = nu (E + rho)
    !!    (derived from (z, z + n) = R + (g, z), z = A^(1/2) grad(u - s),
    !!    n = A^(1/2) grad_h(s - u_h), g = A^(-1/2)(sigma - A grad_h u_h)).
    !!    The SIPG norm adds the exactly computable penalty jump term
    !!    J^2 = sum_F int sigma_F [u_h]^2 (u has no jump), so
    !!    |||u - u_h|||_DG^2 = x^2 + J^2 <= x_bound^2 + J^2. The face flux
    !!    term of the usual DG norm is not controlled. Requires k >= 1 (the hat
    !!    functions must lie in the broken space: SIPG with k = 0 is not
    !!    consistent), geo_order 1 (affine cells; curved cells are rejected)
    !!    and u = g = 0 on the boundary.
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_dg_curved_2d, only: dg_coefficients_t, dg_curved_space_t, &
        dg_data_t, DG_DEFAULT_PENALTY, penalty_factor
    use fortfem_equilibrated_flux_2d, only: equilibrated_flux_estimate, &
        equilibrated_flux_estimate_t, reference_t, build_reference, &
        project_source, EQUILIBRATION_OK, EQUILIBRATION_BAD_INPUT, &
        EQUILIBRATION_BAD_CELL
    use fortfem_oswald_average_2d, only: lagrange_continuous_numbering, &
        oswald_average_homogeneous
    use fortfem_triangle_discontinuous_dof_map, only: &
        build_triangle_discontinuous_dof_map
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_triangle_global_dof_map, only: build_triangle_trimmed_dof_map
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        evaluate_triangle_lagrange_basis, initialize_triangle_lagrange_basis, &
        triangle_lagrange_basis_t
    use fortfem_triangle_piola_maps, only: map_triangle_rt_contravariant
    use fortfem_generated_eq_cell_geometry, only: generated_eq_cell_geometry
    use fortfem_generated_eq_physical_gradient, only: &
        generated_eq_physical_gradient
    use fortfem_generated_eq_stiffness_term, only: generated_eq_stiffness_term
    use fortfem_generated_eq_flux_load_term, only: generated_eq_flux_load_term
    use fortfem_generated_eq_mass_term, only: generated_eq_mass_term
    use fortfem_generated_eq_flux_defect_density, only: &
        generated_eq_flux_defect_density
    use fortfem_generated_eq_energy_density, only: generated_eq_energy_density
    use fortfem_generated_eq_divergence_defect_density, only: &
        generated_eq_divergence_defect_density
    use fortfem_generated_eq_dg_fluxes, only: generated_eq_dg_fluxes
    use fortfem_generated_eq_lifting_load, only: generated_eq_lifting_load
    use fortfem_generated_eq_jump_density, only: generated_eq_jump_density
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private

    public :: nonconforming_estimate_t
    public :: mixed_poisson_flux_bound, sipg_energy_bound
    public :: NONCONFORMING_OK, NONCONFORMING_BAD_INPUT, NONCONFORMING_BAD_CELL
    public :: NONCONFORMING_BAD_COEFFICIENT, NONCONFORMING_SOLVE_FAILED
    public :: MIXED_POTENTIAL_PATCH, MIXED_POTENTIAL_OSWALD

    integer, parameter :: NONCONFORMING_OK = EQUILIBRATION_OK
    integer, parameter :: NONCONFORMING_BAD_INPUT = EQUILIBRATION_BAD_INPUT
    integer, parameter :: NONCONFORMING_BAD_CELL = EQUILIBRATION_BAD_CELL
    integer, parameter :: NONCONFORMING_BAD_COEFFICIENT = 4
    integer, parameter :: NONCONFORMING_SOLVE_FAILED = 5
    integer, parameter :: MIXED_POTENTIAL_PATCH = 1
    integer, parameter :: MIXED_POTENTIAL_OSWALD = 2

    type :: nonconforming_estimate_t
        !> Guaranteed bound: mixed ||A^(-1/2)(sigma - sigma_h)||, SIPG
        !> ||A^(1/2) grad_h (u - u_h)||.
        real(dp) :: bound = huge(1.0_dp)
        !> Mixed: |||u - s||| <= potential_bound. SIPG: DG norm with the
        !> penalty jump term, sqrt(bound^2 + jump_seminorm^2).
        real(dp) :: potential_bound = huge(1.0_dp)
        real(dp) :: flux_defect = 0.0_dp        !< X (mixed) or rho (SIPG)
        real(dp) :: oscillation = 0.0_dp
        real(dp) :: divergence_defect = 0.0_dp  !< C_F ||div sigma + Pi_k f||
        real(dp) :: nonconformity = 0.0_dp      !< SIPG nu
        real(dp) :: jump_seminorm = 0.0_dp      !< SIPG J
        !> Nodal values (P_(k+1) mixed, P_k SIPG) of the conforming
        !> reconstruction s on every cell, (nodes, cells).
        real(dp), allocatable :: potential(:, :)
    end type nonconforming_estimate_t

    abstract interface
        pure function scalar_source_2d(x, y) result(value)
            import :: dp
            real(dp), intent(in) :: x, y
            real(dp) :: value
        end function scalar_source_2d
    end interface

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(*), info
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
        end subroutine dgesv
    end interface

contains

    subroutine mixed_poisson_flux_bound(mesh, degree, quadrature_degree, &
            flux_dofs, pressure_dofs, source, friedrichs_constant, estimate, &
            status, reconstruction)
        !! Flux-error bound of solve_mixed_poisson_rt (A = 1, flux sigma_h =
        !! -grad p_h) on the counter-clockwise mesh; flux_dofs/pressure_dofs
        !! are its outputs. The rule of degree quadrature_degree must
        !! integrate polynomials of degree 2k + 2 exactly. reconstruction
        !! selects the H^1_0 potential s (any s is admissible, the bound is
        !! always rigorous): MIXED_POTENTIAL_PATCH (default) the vertex-patch
        !! H^1 projection of grad(psi_a ptilde) into P_(k+2) (partition of
        !! unity), or MIXED_POTENTIAL_OSWALD the Oswald average of ptilde in
        !! P_(k+1), where ptilde is the elementwise P_(k+1) postprocessing of
        !! sigma_h and p_h. The patch variant is far tighter (efficiency
        !! 1.04-1.5 against 2-3).
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: flux_dofs(:), pressure_dofs(:)
        procedure(scalar_source_2d) :: source
        real(dp), intent(in) :: friedrichs_constant
        type(nonconforming_estimate_t), intent(out) :: estimate
        integer, intent(out) :: status
        integer, intent(in), optional :: reconstruction

        type(reference_t) :: ref
        type(triangle_lagrange_basis_t) :: pbasis, dbasis
        integer, allocatable :: flux_map(:, :), transforms(:, :), pressure_map(:, :)
        real(dp), allocatable :: xi(:), eta(:), wq(:)
        real(dp), allocatable :: pval(:, :), pgrad(:, :, :), dval(:, :)
        real(dp), allocatable :: sgrad(:, :, :), sval(:, :)
        type(triangle_lagrange_basis_t) :: sbasis
        integer :: ns
        real(dp), allocatable :: local_s(:, :), sigma(:, :, :), div(:, :)
        real(dp), allocatable :: osc(:), defect(:), divergence_sq(:)
        real(dp), allocatable :: projection(:, :)
        real(dp) :: dgrad(2, (degree + 1)*(degree + 2)/2)
        integer :: nt, nv, nq, nl, k, t, q, flux_count, pressure_count, st, failed
        real(dp) :: x_norm, e_total
        integer :: method

        status = NONCONFORMING_BAD_INPUT
        k = degree
        method = MIXED_POTENTIAL_PATCH
        if (present(reconstruction)) method = reconstruction
        if (method /= MIXED_POTENTIAL_PATCH .and. &
            method /= MIXED_POTENTIAL_OSWALD) return
        nt = mesh%n_triangles
        nv = mesh%n_vertices
        if (k < 0 .or. nt < 1 .or. quadrature_degree < 2*k + 2) return
        if (.not. (friedrichs_constant >= 0.0_dp)) return
        call build_triangle_trimmed_dof_map(mesh, k + 1, flux_map, transforms, &
            flux_count, st)
        if (st /= 0) return
        call build_triangle_discontinuous_dof_map(mesh, k, pressure_map, &
            pressure_count, st)
        if (st /= 0) return
        if (size(flux_dofs) /= flux_count) return
        if (size(pressure_dofs) /= pressure_count) return
        call triangle_duffy_quadrature(quadrature_degree, xi, eta, wq, st)
        if (st /= 0) return
        nq = size(wq)
        call build_reference(k, xi, eta, wq, ref, st)
        if (st /= 0) return
        call initialize_triangle_lagrange_basis(k + 1, pbasis, st)
        if (st /= 0) return
        call initialize_triangle_lagrange_basis(k, dbasis, st)
        if (st /= 0) return
        nl = (k + 2)*(k + 3)/2
        allocate (pval(nl, nq), pgrad(2, nl, nq), dval(ref%nk, nq))
        ! The potential s lives in P_(k+1) (Oswald) or P_(k+2) (patches: the
        ! product psi_a ptilde of a hat function and a P_(k+1) function).
        ns = nl
        if (method == MIXED_POTENTIAL_PATCH) ns = (k + 3)*(k + 4)/2
        call initialize_triangle_lagrange_basis(merge(k + 2, k + 1, &
            method == MIXED_POTENTIAL_PATCH), sbasis, st)
        if (st /= 0) return
        allocate (sval(ns, nq), sgrad(2, ns, nq))
        do q = 1, nq
            call evaluate_triangle_lagrange_basis(pbasis, xi(q), eta(q), &
                pval(:, q), pgrad(:, :, q), st)
            if (st /= 0) return
            call evaluate_triangle_lagrange_basis(sbasis, xi(q), eta(q), &
                sval(:, q), sgrad(:, :, q), st)
            if (st /= 0) return
            call evaluate_triangle_lagrange_basis(dbasis, xi(q), eta(q), &
                dval(:, q), dgrad, st)
            if (st /= 0) return
        end do
        status = NONCONFORMING_BAD_CELL
        allocate (local_s(nl, nt), sigma(2, nq, nt), div(nq, nt))
        allocate (osc(nt), defect(nt), divergence_sq(nt))
        allocate (projection(ref%nk, nt))
        failed = 0
        !$omp parallel do schedule(static) reduction(max:failed) &
        !$omp private(st)
        do t = 1, nt
            call postprocess_cell(t, st)
            failed = max(failed, st)
        end do
        !$omp end parallel do
        if (failed /= 0) return
        allocate (estimate%potential(ns, nt))
        if (method == MIXED_POTENTIAL_OSWALD) then
            call oswald_average_homogeneous(mesh%triangles(:, :nt), k + 1, &
                local_s, estimate%potential, st)
            if (st /= 0) then
                status = NONCONFORMING_SOLVE_FAILED
                return
            end if
        else
            call patch_potential(mesh%vertices, mesh%triangles(:, :nt), xi, eta, &
                wq, pval, sgrad, sigma, local_s, estimate%potential, st)
            if (st /= 0) then
                status = NONCONFORMING_SOLVE_FAILED
                return
            end if
        end if
        failed = 0
        !$omp parallel do schedule(static) reduction(max:failed) &
        !$omp private(st)
        do t = 1, nt
            call defect_cell(t, st)
            failed = max(failed, st)
        end do
        !$omp end parallel do
        if (failed /= 0) return
        x_norm = sqrt(sum(defect))
        estimate%flux_defect = x_norm
        estimate%oscillation = sqrt(sum(osc**2))
        estimate%divergence_defect = friedrichs_constant*sqrt(sum(divergence_sq))
        e_total = estimate%oscillation + estimate%divergence_defect
        estimate%bound = sqrt(x_norm**2 + e_total**2)
        estimate%potential_bound = x_norm + estimate%bound
        status = NONCONFORMING_OK

    contains

        subroutine cell_data(t, vertices, jac, det, jit, diameter)
            integer, intent(in) :: t
            real(dp), intent(out) :: vertices(2, 3), jac(2, 2), det, jit(2, 2)
            real(dp), intent(out) :: diameter
            real(dp) :: l12, l23, l31

            vertices = mesh%vertices(:, mesh%triangles(1:3, t))
            call generated_eq_cell_geometry(vertices(1, 1), vertices(2, 1), &
                vertices(1, 2), vertices(2, 2), vertices(1, 3), vertices(2, 3), &
                det, jit(1, 1), jit(1, 2), jit(2, 1), jit(2, 2), l12, l23, l31)
            jac(:, 1) = vertices(:, 2) - vertices(:, 1)
            jac(:, 2) = vertices(:, 3) - vertices(:, 1)
            diameter = sqrt(max(l12, l23, l31))
        end subroutine cell_data

        subroutine postprocess_cell(t, st)
            !! sigma_h, p_h samples and the local P_(k+1) minimiser of
            !! int |sigma_h + grad s|^2 subject to mean s = mean p_h.
            integer, intent(in) :: t
            integer, intent(out) :: st
            real(dp) :: vertices(2, 3), jac(2, 2), det, jit(2, 2), diameter
            real(dp) :: phi(2, ref%nrt), pdiv(ref%nrt), local_flux(ref%nrt)
            real(dp) :: grad(2, nl, nq), kkt(nl + 1, nl + 1), rhs(nl + 1, 1)
            real(dp) :: wdet, term, pmean, local_p(ref%nk)
            real(dp) :: f(nq), point(2)
            integer :: i, j, q, ipiv(nl + 1), info

            call cell_data(t, vertices, jac, det, jit, diameter)
            st = NONCONFORMING_BAD_CELL
            if (.not. (det > 0.0_dp)) return
            local_flux = real(transforms(:, t), dp)*flux_dofs(flux_map(:, t))
            local_p = pressure_dofs(pressure_map(:, t))
            kkt = 0.0_dp
            rhs = 0.0_dp
            pmean = 0.0_dp
            do q = 1, nq
                call map_triangle_rt_contravariant(jac, ref%rt(:, :, q), &
                    ref%rt_div(:, q), phi, pdiv, info)
                if (info /= 0) return
                sigma(:, q, t) = matmul(phi, local_flux)
                div(q, t) = dot_product(pdiv, local_flux)
                wdet = wq(q)*det
                do i = 1, nl
                    call generated_eq_physical_gradient(jit(1, 1), jit(1, 2), &
                        jit(2, 1), jit(2, 2), pgrad(1, i, q), pgrad(2, i, q), &
                        grad(1, i, q), grad(2, i, q))
                end do
                do i = 1, nl
                    call generated_eq_flux_load_term(wdet, sigma(1, q, t), &
                        sigma(2, q, t), grad(1, i, q), grad(2, i, q), term)
                    rhs(i, 1) = rhs(i, 1) - term
                    do j = 1, nl
                        call generated_eq_stiffness_term(wdet, 1.0_dp, &
                            grad(1, i, q), grad(2, i, q), grad(1, j, q), &
                            grad(2, j, q), term)
                        kkt(i, j) = kkt(i, j) + term
                    end do
                    call generated_eq_mass_term(wdet, 1.0_dp, pval(i, q), &
                        1.0_dp, term)
                    kkt(i, nl + 1) = kkt(i, nl + 1) + term
                    kkt(nl + 1, i) = kkt(nl + 1, i) + term
                end do
                pmean = pmean + wdet*dot_product(dval(:, q), local_p)
                point = vertices(:, 1) + xi(q)*jac(:, 1) + eta(q)*jac(:, 2)
                f(q) = source(point(1), point(2))
            end do
            rhs(nl + 1, 1) = pmean
            call dgesv(nl + 1, 1, kkt, nl + 1, ipiv, rhs, nl + 1, info)
            if (info /= 0) then
                st = NONCONFORMING_SOLVE_FAILED
                return
            end if
            local_s(:, t) = rhs(1:nl, 1)
            call project_source(ref, mesh%vertices, mesh%triangles(1:3, t), f, 1.0_dp, &
                projection(:, t), osc(t))
            st = NONCONFORMING_OK
        end subroutine postprocess_cell

        subroutine defect_cell(t, st)
            integer, intent(in) :: t
            integer, intent(out) :: st
            real(dp) :: vertices(2, 3), jac(2, 2), det, jit(2, 2), diameter
            real(dp) :: gs(2), gref(2), wdet, term, pf
            integer :: qq

            call cell_data(t, vertices, jac, det, jit, diameter)
            st = NONCONFORMING_OK
            defect(t) = 0.0_dp
            divergence_sq(t) = 0.0_dp
            do qq = 1, nq
                wdet = wq(qq)*det
                gref = matmul(sgrad(:, :, qq), estimate%potential(:, t))
                call generated_eq_physical_gradient(jit(1, 1), jit(1, 2), &
                    jit(2, 1), jit(2, 2), gref(1), gref(2), gs(1), gs(2))
                call generated_eq_flux_defect_density(wdet, 1.0_dp, &
                    sigma(1, qq, t), sigma(2, qq, t), gs(1), gs(2), term)
                defect(t) = defect(t) + term
                pf = dot_product(projection(:, t), ref%poly(:, qq))
                call generated_eq_divergence_defect_density(wdet, div(qq, t), &
                    pf, term)
                divergence_sq(t) = divergence_sq(t) + term
            end do
        end subroutine defect_cell

    end subroutine mixed_poisson_flux_bound

    subroutine patch_potential(vertices, triangles, xi, eta, wq, pval, pgrad, &
            flux, ptilde, potential, status)
        !! Vertex-patch H^1_0 potential reconstruction, continuous P_m. For
        !! every vertex a the function s_a in P_m vanishing on the boundary
        !! of the patch and on the domain boundary is the (dense SPD solve)
        !! H^1 projection of g_a = -psi_a flux + ptilde grad psi_a on the
        !! patch, (grad s_a, grad v) = (g_a, grad v) for all such v (Ern--
        !! Vohralik partition of unity): g_a is the broken gradient of psi_a
        !! ptilde for a broken potential ptilde with grad ptilde = -flux, and
        !! sum_a g_a = -flux. s = sum_a s_a lies in H^1_0; the degree m = k + 2
        !! makes psi_a ptilde (ptilde in P_(k+1)) exactly representable.
        !! pgrad(2, nl, nq): reference gradients of the P_m basis at the rule
        !! (xi, eta, wq); flux(2, nq, nt) physical samples; ptilde(:, nt)
        !! broken potential nodal values and pval(:, nq) the values at the rule
        !! of its (lower degree) basis.
        real(dp), intent(in) :: vertices(:, :), xi(:), eta(:), wq(:)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: pval(:, :), pgrad(:, :, :), flux(:, :, :)
        real(dp), intent(in) :: ptilde(:, :)
        real(dp), intent(out) :: potential(:, :)
        integer, intent(out) :: status
        integer, allocatable :: bary(:, :), node_id(:, :), vstart(:), vcell(:)
        integer, allocatable :: vslot(:)
        logical, allocatable :: boundary_node(:)
        real(dp), allocatable :: contribution(:, :, :)
        integer :: nt, nv, nl, nq, m, a, t, s, i, failed, st

        status = NONCONFORMING_BAD_INPUT
        nt = size(triangles, 2)
        nv = size(vertices, 2)
        nl = size(pgrad, 2)
        nq = size(wq)
        m = nint(0.5_dp*(sqrt(8.0_dp*nl + 1.0_dp) - 3.0_dp))
        call lagrange_continuous_numbering(triangles, m, bary, node_id, &
            boundary_node, st)
        if (st /= 0) return
        allocate (vstart(nv + 1), vcell(3*nt), vslot(3*nt))
        vstart = 0
        do t = 1, nt
            do s = 1, 3
                vstart(triangles(s, t) + 1) = vstart(triangles(s, t) + 1) + 1
            end do
        end do
        vstart(1) = 1
        do a = 1, nv
            vstart(a + 1) = vstart(a + 1) + vstart(a)
        end do
        block
            integer :: fill(nv)
            fill = vstart(1:nv)
            do t = 1, nt
                do s = 1, 3
                    a = triangles(s, t)
                    vcell(fill(a)) = t
                    vslot(fill(a)) = s
                    fill(a) = fill(a) + 1
                end do
            end do
        end block
        allocate (contribution(nl, 3, nt))
        contribution = 0.0_dp
        failed = 0
        !$omp parallel do schedule(dynamic, 16) reduction(max:failed)
        do a = 1, nv
            if (vstart(a + 1) == vstart(a)) cycle
            call solve_vertex(a, failed)
        end do
        !$omp end parallel do
        status = NONCONFORMING_SOLVE_FAILED
        if (failed /= 0) return
        do i = 1, nl
            potential(i, :) = sum(contribution(i, :, :), dim=1)
        end do
        status = NONCONFORMING_OK

    contains

        subroutine solve_vertex(a, failed)
            integer, intent(in) :: a
            integer, intent(inout) :: failed
            integer :: cells, free_count, i, j, l, ll, q, g, info, t, s
            integer :: free_id(nl*(vstart(a + 1) - vstart(a)))
            integer :: local_free(nl, vstart(a + 1) - vstart(a))
            integer, allocatable :: ipiv(:)
            real(dp), allocatable :: matrix(:, :), rhs(:, :)
            real(dp) :: v(2, 3), det, jit(2, 2), l12, l23, l31, wdet, term
            real(dp) :: grad(2, nl), hat(2), weight
            real(dp), parameter :: hat_reference(2, 3) = reshape( &
                [-1.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])

            cells = vstart(a + 1) - vstart(a)
            free_count = 0
            local_free = 0
            do i = 1, cells
                t = vcell(vstart(a) + i - 1)
                s = vslot(vstart(a) + i - 1)
                do l = 1, nl
                    if (bary(s, l) == 0) cycle
                    g = node_id(l, t)
                    if (boundary_node(g)) cycle
                    j = 0
                    do ll = 1, free_count
                        if (free_id(ll) == g) j = ll
                    end do
                    if (j == 0) then
                        free_count = free_count + 1
                        free_id(free_count) = g
                        j = free_count
                    end if
                    local_free(l, i) = j
                end do
            end do
            if (free_count == 0) return
            allocate (matrix(free_count, free_count), rhs(free_count, 1))
            allocate (ipiv(free_count))
            matrix = 0.0_dp
            rhs = 0.0_dp
            do i = 1, cells
                t = vcell(vstart(a) + i - 1)
                s = vslot(vstart(a) + i - 1)
                v = vertices(:, triangles(1:3, t))
                call generated_eq_cell_geometry(v(1, 1), v(2, 1), v(1, 2), &
                    v(2, 2), v(1, 3), v(2, 3), det, jit(1, 1), jit(1, 2), &
                    jit(2, 1), jit(2, 2), l12, l23, l31)
                call generated_eq_physical_gradient(jit(1, 1), jit(1, 2), &
                    jit(2, 1), jit(2, 2), hat_reference(1, s), &
                    hat_reference(2, s), hat(1), hat(2))
                do q = 1, nq
                    wdet = wq(q)*det
                    weight = merge(1.0_dp - xi(q) - eta(q), &
                        merge(xi(q), eta(q), s == 2), s == 1)
                    do l = 1, nl
                        call generated_eq_physical_gradient(jit(1, 1), jit(1, 2), &
                            jit(2, 1), jit(2, 2), pgrad(1, l, q), pgrad(2, l, q), &
                            grad(1, l), grad(2, l))
                    end do
                    do l = 1, nl
                        if (local_free(l, i) == 0) cycle
                        call generated_eq_flux_load_term(wdet*weight, &
                            flux(1, q, t), flux(2, q, t), grad(1, l), &
                            grad(2, l), term)
                        rhs(local_free(l, i), 1) = rhs(local_free(l, i), 1) - term
                        call generated_eq_stiffness_term(wdet, &
                            dot_product(pval(:, q), ptilde(:, t)), hat(1), &
                            hat(2), grad(1, l), grad(2, l), term)
                        rhs(local_free(l, i), 1) = rhs(local_free(l, i), 1) + term
                        do ll = 1, nl
                            if (local_free(ll, i) == 0) cycle
                            call generated_eq_stiffness_term(wdet, 1.0_dp, &
                                grad(1, l), grad(2, l), grad(1, ll), grad(2, ll), &
                                term)
                            matrix(local_free(l, i), local_free(ll, i)) = &
                                matrix(local_free(l, i), local_free(ll, i)) + term
                        end do
                    end do
                end do
            end do
            call dgesv(free_count, 1, matrix, free_count, ipiv, rhs, free_count, &
                info)
            if (info /= 0) then
                failed = 1
                return
            end if
            do i = 1, cells
                t = vcell(vstart(a) + i - 1)
                s = vslot(vstart(a) + i - 1)
                do l = 1, nl
                    if (local_free(l, i) == 0) cycle
                    contribution(l, s, t) = rhs(local_free(l, i), 1)
                end do
            end do
        end subroutine solve_vertex

    end subroutine patch_potential

    subroutine sipg_energy_bound(space, coef, data, u, quadrature_degree, &
            coefficient_lower, friedrichs_constant, estimate, status, penalty)
        !! Bound of the SIPG solution u (dg_curved_assemble/dg_curved_rhs with
        !! real isotropic K = A I, b = 0, c = 0, g = 0 on an affine mesh of
        !! counter-clockwise triangles, k >= 1); only the real part of u is
        !! used. penalty is the eta of the assembly (default DG_DEFAULT_PENALTY).
        !! coefficient_lower(t) <= A on cell t. The rule of degree
        !! quadrature_degree must integrate polynomials of degree 2k + 2 exactly.
        type(dg_curved_space_t), intent(in) :: space
        class(dg_coefficients_t), intent(in) :: coef
        class(dg_data_t), intent(in) :: data
        complex(dp), intent(in) :: u(:)
        integer, intent(in) :: quadrature_degree
        real(dp), intent(in) :: coefficient_lower(:), friedrichs_constant
        type(nonconforming_estimate_t), intent(out) :: estimate
        integer, intent(out) :: status
        real(dp), intent(in), optional :: penalty

        type(equilibrated_flux_estimate_t) :: flux
        type(triangle_lagrange_basis_t) :: basis
        real(dp), allocatable :: xi(:), eta(:), wq(:), vval(:, :), vgrad(:, :, :)
        real(dp), allocatable :: fnode(:), fw(:), fval(:, :, :)
        real(dp), allocatable :: acoef(:, :), fsrc(:, :), broken_flux(:, :, :)
        real(dp), allocatable :: load_flux(:, :, :), nodal(:, :), jump_sq(:)
        real(dp), allocatable :: energy_sq(:)
        real(dp), parameter :: vref(2, 3) = reshape( &
            [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        real(dp) :: eta_penalty, b, c, e_total, rho, nu, x_bound
        integer, allocatable :: triangles(:, :)
        integer :: nt, k, nl, nq, nfq, t, q, e, st, failed

        status = NONCONFORMING_BAD_INPUT
        k = space%k
        nt = space%nt
        nl = space%nl
        if (k < 1 .or. nt < 1 .or. quadrature_degree < 2*k + 2) return
        if (space%mm%geo_order /= 1) return
        if (.not. allocated(space%face_t)) return
        if (size(u) /= space%ndof) return
        if (size(coefficient_lower) /= nt) return
        if (any(.not. (coefficient_lower > 0.0_dp))) return
        if (.not. (friedrichs_constant >= 0.0_dp)) return
        eta_penalty = DG_DEFAULT_PENALTY
        if (present(penalty)) eta_penalty = penalty
        triangles = space%mm%topology%triangles(:, :nt)
        call triangle_duffy_quadrature(quadrature_degree, xi, eta, wq, st)
        if (st /= 0) return
        nq = size(wq)
        call initialize_triangle_lagrange_basis(k, basis, st)
        if (st /= 0) return
        nfq = quadrature_degree/2 + 1
        allocate (fnode(nfq), fw(nfq), fval(nl, nfq, 3))
        call gauss_legendre_ab(nfq, 0.0_dp, 1.0_dp, fnode, fw)
        allocate (vval(nl, nq), vgrad(2, nl, nq))
        block
            real(dp) :: dummy(2, nl)
            do q = 1, nq
                call evaluate_triangle_lagrange_basis(basis, xi(q), eta(q), &
                    vval(:, q), vgrad(:, :, q), st)
                if (st /= 0) return
            end do
            do e = 1, 3
                do q = 1, nfq
                    call evaluate_triangle_lagrange_basis(basis, &
                        vref(1, e) + fnode(q)*(vref(1, mod(e, 3) + 1) - vref(1, e)), &
                        vref(2, e) + fnode(q)*(vref(2, mod(e, 3) + 1) - vref(2, e)), &
                        fval(:, q, e), dummy, st)
                    if (st /= 0) return
                end do
            end do
        end block
        allocate (acoef(nq, nt), fsrc(nq, nt), broken_flux(2, nq, nt))
        allocate (load_flux(2, nq, nt), nodal(nl, nt), jump_sq(nt))
        allocate (energy_sq(nt))
        do t = 1, nt
            nodal(:, t) = real(u((t - 1)*nl + 1:t*nl), dp)
        end do
        status = NONCONFORMING_BAD_CELL
        failed = 0
        !$omp parallel do schedule(dynamic, 8) reduction(max:failed) private(st)
        do t = 1, nt
            call lift_cell(t, st)
            failed = max(failed, st)
        end do
        !$omp end parallel do
        status = failed
        if (failed /= 0) return

        call equilibrated_flux_estimate(space%mm%topology%vertices, triangles, k, &
            xi, eta, wq, acoef, fsrc, broken_flux, coefficient_lower, &
            friedrichs_constant, flux, st, load_flux)
        if (st /= EQUILIBRATION_OK) then
            status = st
            return
        end if
        allocate (estimate%potential(nl, nt))
        call oswald_average_homogeneous(triangles, k, nodal, &
            estimate%potential, st)
        if (st /= 0) then
            status = NONCONFORMING_SOLVE_FAILED
            return
        end if
        !$omp parallel do schedule(static)
        do t = 1, nt
            call energy_cell(t)
        end do
        !$omp end parallel do
        rho = flux%flux_defect
        nu = sqrt(sum(energy_sq))
        e_total = flux%oscillation + flux%divergence_defect
        b = e_total + rho + nu
        c = nu*(e_total + rho)
        x_bound = 0.5_dp*(b + sqrt(b*b + 4.0_dp*c))
        estimate%flux_defect = rho
        estimate%oscillation = flux%oscillation
        estimate%divergence_defect = flux%divergence_defect
        estimate%nonconformity = nu
        estimate%jump_seminorm = sqrt(sum(jump_sq))
        estimate%bound = x_bound
        estimate%potential_bound = sqrt(x_bound**2 + estimate%jump_seminorm**2)
        status = NONCONFORMING_OK

    contains

        subroutine lift_cell(t, st)
            !! Samples of A, f, the broken flux A grad u_h and the load flux
            !! A (grad u_h - l) of cell t, and the penalty jump term.
            integer, intent(in) :: t
            integer, intent(out) :: st
            real(dp) :: vertices(2, 3), jac(2, 2), det, jit(2, 2)
            real(dp) :: x(2, nq), grad(2, nq), lift(nq, 2)
            real(dp) :: mass(nl, nl), rhs(nl, 2), wdet, term, gref(2)
            real(dp) :: l12, l23, l31
            complex(dp) :: kk(2, 2, nq), cc(nq), fz(nq)
            real(dp) :: bb(2, nq)
            integer :: ipiv(nl), info, q, i, j, edge

            vertices = space%mm%topology%vertices(:, triangles(1:3, t))
            call generated_eq_cell_geometry(vertices(1, 1), vertices(2, 1), &
                vertices(1, 2), vertices(2, 2), vertices(1, 3), vertices(2, 3), &
                det, jit(1, 1), jit(1, 2), jit(2, 1), jit(2, 2), l12, l23, l31)
            jac(:, 1) = vertices(:, 2) - vertices(:, 1)
            jac(:, 2) = vertices(:, 3) - vertices(:, 1)
            st = NONCONFORMING_BAD_CELL
            if (.not. (det > 0.0_dp)) return
            do q = 1, nq
                x(:, q) = vertices(:, 1) + xi(q)*jac(:, 1) + eta(q)*jac(:, 2)
            end do
            call coef%evaluate(t, x, kk, bb, cc, info)
            if (info /= 0) return
            call data%source(t, x, fz, info)
            if (info /= 0) return
            st = NONCONFORMING_BAD_COEFFICIENT
            do q = 1, nq
                if (.not. isotropic_diffusion(kk(:, :, q), bb(:, q), cc(q))) return
                acoef(q, t) = real(kk(1, 1, q), dp)
                fsrc(q, t) = real(fz(q), dp)
            end do
            mass = 0.0_dp
            do q = 1, nq
                wdet = wq(q)*det
                gref = matmul(vgrad(:, :, q), nodal(:, t))
                call generated_eq_physical_gradient(jit(1, 1), jit(1, 2), &
                    jit(2, 1), jit(2, 2), gref(1), gref(2), grad(1, q), grad(2, q))
                do i = 1, nl
                    do j = 1, nl
                        call generated_eq_mass_term(wdet, acoef(q, t), &
                            vval(i, q), vval(j, q), term)
                        mass(i, j) = mass(i, j) + term
                    end do
                end do
            end do
            rhs = 0.0_dp
            jump_sq(t) = 0.0_dp
            do edge = 1, 3
                call face_terms(t, edge, vertices, rhs, jump_sq(t), st)
                if (st /= NONCONFORMING_OK) return
            end do
            call dgesv(nl, 2, mass, nl, ipiv, rhs, nl, info)
            st = NONCONFORMING_SOLVE_FAILED
            if (info /= 0) return
            lift = matmul(transpose(vval), rhs)
            do q = 1, nq
                call generated_eq_dg_fluxes(acoef(q, t), grad(1, q), grad(2, q), &
                    lift(q, 1), lift(q, 2), broken_flux(1, q, t), &
                    broken_flux(2, q, t), load_flux(1, q, t), load_flux(2, q, t))
            end do
            st = NONCONFORMING_OK
        end subroutine lift_cell

        subroutine face_terms(t, edge, vertices, rhs, jump_total, st)
            !! Lifting load of local edge `edge` of cell t and its share of
            !! the squared penalty jump seminorm.
            integer, intent(in) :: t, edge
            real(dp), intent(in) :: vertices(2, 3)
            real(dp), intent(inout) :: rhs(:, :), jump_total
            integer, intent(out) :: st
            real(dp) :: xf(2, nfq), tangent(2), length, normal(2), weight, w
            real(dp) :: a_here(nfq), a_there(nfq), jump, sigma, ratio, term(2)
            real(dp) :: bb(2, nfq), u_here, u_there
            complex(dp) :: kk(2, 2, nfq), cc(nfq)
            integer :: f, side, nb, e2, q, q2, l, info, v2
            logical :: interior

            st = NONCONFORMING_BAD_CELL
            f = space%elem_face(edge, t)
            side = space%elem_side(edge, t)
            interior = f <= space%nint
            v2 = mod(edge, 3) + 1
            tangent = vertices(:, v2) - vertices(:, edge)
            length = norm2(tangent)
            normal = [tangent(2), -tangent(1)]/length
            do q = 1, nfq
                xf(:, q) = vertices(:, edge) + fnode(q)*tangent
            end do
            call coef%evaluate(t, xf, kk, bb, cc, info)
            if (info /= 0) return
            a_here = real(kk(1, 1, :), dp)
            a_there = a_here
            ratio = space%flen(f)/space%area(t)
            weight = 1.0_dp
            nb = 0
            e2 = 0
            if (interior) then
                nb = space%face_t(3 - side, f)
                e2 = space%face_e(3 - side, f)
                call coef%evaluate(nb, xf, kk, bb, cc, info)
                if (info /= 0) return
                a_there = real(kk(1, 1, :), dp)
                ratio = max(ratio, space%flen(f)/space%area(nb))
                weight = 0.5_dp
            end if
            do q = 1, nfq
                u_here = dot_product(fval(:, q, edge), nodal(:, t))
                jump = u_here
                if (interior) then
                    q2 = q
                    if (space%face_rev(f)) q2 = nfq + 1 - q
                    u_there = dot_product(fval(:, q2, e2), nodal(:, nb))
                    jump = u_here - u_there
                end if
                w = fw(q)*length*weight
                do l = 1, nl
                    call generated_eq_lifting_load(w, a_here(q), jump, normal(1), &
                        normal(2), fval(l, q, edge), term(1), term(2))
                    rhs(l, :) = rhs(l, :) + term
                end do
                sigma = penalty_factor(space, eta_penalty, ratio) &
                    *max(a_here(q), a_there(q))
                call generated_eq_jump_density(w, sigma, jump, term(1))
                jump_total = jump_total + term(1)
            end do
            st = NONCONFORMING_OK
        end subroutine face_terms

        subroutine energy_cell(t)
            integer, intent(in) :: t
            real(dp) :: vertices(2, 3), det, jit(2, 2), g(2), gref(2), term
            real(dp) :: l1, l2, l3
            integer :: q

            vertices = space%mm%topology%vertices(:, triangles(1:3, t))
            call generated_eq_cell_geometry(vertices(1, 1), vertices(2, 1), &
                vertices(1, 2), vertices(2, 2), vertices(1, 3), vertices(2, 3), &
                det, jit(1, 1), jit(1, 2), jit(2, 1), jit(2, 2), l1, l2, l3)
            energy_sq(t) = 0.0_dp
            do q = 1, nq
                gref = matmul(vgrad(:, :, q), nodal(:, t) - estimate%potential(:, t))
                call generated_eq_physical_gradient(jit(1, 1), jit(1, 2), &
                    jit(2, 1), jit(2, 2), gref(1), gref(2), g(1), g(2))
                call generated_eq_energy_density(wq(q)*det, acoef(q, t), g(1), &
                    g(2), term)
                energy_sq(t) = energy_sq(t) + term
            end do
        end subroutine energy_cell

    end subroutine sipg_energy_bound

    pure logical function isotropic_diffusion(kk, b, c) result(ok)
        !! K = A I with real A > 0, no advection and no reaction.
        complex(dp), intent(in) :: kk(2, 2), c
        real(dp), intent(in) :: b(2)
        real(dp) :: scale

        scale = abs(kk(1, 1))
        ok = real(kk(1, 1), dp) > 0.0_dp
        if (.not. ok) return
        ok = abs(aimag(kk(1, 1))) <= 1.0e-12_dp*scale .and. &
            abs(kk(1, 2)) <= 1.0e-12_dp*scale .and. &
            abs(kk(2, 1)) <= 1.0e-12_dp*scale .and. &
            abs(kk(2, 2) - kk(1, 1)) <= 1.0e-12_dp*scale .and. &
            abs(c) == 0.0_dp .and. all(b == 0.0_dp)
    end function isotropic_diffusion

end module fortfem_equilibrated_nonconforming_2d
