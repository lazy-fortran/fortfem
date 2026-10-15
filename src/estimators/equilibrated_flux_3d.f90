module fortfem_equilibrated_flux_3d
    !! Patchwise equilibrated-flux (Prager--Synge) energy-error bound for
    !! -div(A grad u) = f on an affine tetrahedral mesh, u = 0 on the whole
    !! boundary, scalar A > 0, for any u_h in H^1_0 (here sampled as
    !! discrete_flux = A grad u_h). The 3D analogue of
    !! fortfem_equilibrated_flux_2d with the FortFEM tetrahedral RT_k space
    !! (dimension (k+1)(k+2)(k+4)/2, P_k face moments).
    !!
    !! For every vertex a with hat function psi_a, the local problem on its
    !! patch finds sigma_a in the RT_k global dofs of the patch: the face dofs
    !! of faces containing a and the interior dofs of the patch cells (zero
    !! normal trace on the patch boundary away from the domain boundary) that
    !! minimizes ||A^(-1/2)(sigma_a - psi_a A grad u_h)|| subject to
    !! -div sigma_a = Pi_k(f psi_a - A grad u_h . grad psi_a) - c_a, with c_a
    !! the patch mean of that load on interior patches. sigma = sum_a sigma_a
    !! is in H(div) and
    !!
    !!   |||u - u_h||| <= ||A^(-1/2)(sigma - A grad u_h)||
    !!       + (sum_T (h_T/pi)^2 / A_T,min ||f - Pi_k f||_T^2)^(1/2)
    !!       + C_F ||div sigma + Pi_k f||,
    !!
    !! with h_T the diameter (Payne--Weinberger, convex tetrahedron),
    !! A_T,min a caller-supplied lower bound of A on T and C_F a caller
    !! supplied Friedrichs constant ||v|| <= C_F |||v||| on H^1_0. The
    !! global dofs, their transforms (face orientation and permutation) come
    !! from FortFEM's tetrahedral RT dof map, so sigma_a is conforming by
    !! construction. Pointwise formulas are FortSym-generated
    !! (gen_equilibrated_flux_3d_products); topology, scatter, dense solves
    !! and OpenMP are handwritten. Quadrature error is not enclosed.
    use fortfem_kinds, only: dp
    use fortfem_tetra_rt_arbitrary_order, only: tetra_rt_t, initialize_tetra_rt, &
        evaluate_tetra_rt
    use fortfem_tetra_rt_global_dof_map, only: build_tetra_rt_dof_map, &
        build_tetra_rt_basis_transform
    use fortfem_tetra_piola_maps, only: map_tetra_rt_contravariant
    use fortfem_equilibrated_flux_topology, only: vertex_stars, &
        EQUILIBRATION_OK, EQUILIBRATION_BAD_INPUT, EQUILIBRATION_BAD_CELL, &
        EQUILIBRATION_PATCH_FAILED
    use fortfem_generated_tetra_simplex_geometry, only: &
        generated_tetra_simplex_geometry
    use fortfem_generated_tetra_monomials_degree0, only: &
        generated_tetra_monomials_degree0
    use fortfem_generated_tetra_monomials_degree1, only: &
        generated_tetra_monomials_degree1
    use fortfem_generated_tetra_monomials_degree2, only: &
        generated_tetra_monomials_degree2
    use fortfem_generated_tetra_monomials_degree3, only: &
        generated_tetra_monomials_degree3
    use fortfem_generated_tetra_monomials_degree4, only: &
        generated_tetra_monomials_degree4
    use fortfem_generated_tetra_edge_lengths_squared, only: &
        generated_tetra_edge_lengths_squared
    use fortfem_generated_tetra_patch_load, only: generated_tetra_patch_load
    use fortfem_generated_tetra_patch_gram, only: generated_tetra_patch_gram
    use fortfem_generated_tetra_patch_flux_rhs, only: &
        generated_tetra_patch_flux_rhs
    use fortfem_generated_tetra_patch_multiplier, only: &
        generated_tetra_patch_multiplier
    use fortfem_generated_tetra_projection_terms, only: &
        generated_tetra_projection_terms
    use fortfem_generated_tetra_flux_defect, only: generated_tetra_flux_defect
    use fortfem_generated_tetra_divergence_defect, only: &
        generated_tetra_divergence_defect
    use fortfem_generated_tetra_oscillation_defect, only: &
        generated_tetra_oscillation_defect
    implicit none
    private

    public :: equilibrated_flux_estimate_3d_t, equilibrated_flux_estimate_3d

    public :: EQUILIBRATION_OK, EQUILIBRATION_BAD_INPUT, EQUILIBRATION_BAD_CELL, &
        EQUILIBRATION_PATCH_FAILED

    type :: equilibrated_flux_estimate_3d_t
        real(dp) :: bound = huge(1.0_dp)        !< sum of the three terms
        real(dp) :: flux_defect = 0.0_dp        !< ||A^(-1/2)(sigma - A grad u_h)||
        real(dp) :: oscillation = 0.0_dp        !< data-oscillation term
        real(dp) :: divergence_defect = 0.0_dp  !< C_F ||div sigma + Pi_k f||
        !> Per cell: flux-defect and oscillation contributions (not squared).
        real(dp), allocatable :: cell_flux(:), cell_oscillation(:)
        !> Element-local physical RT_k coefficients of sigma, (dofs, cells).
        real(dp), allocatable :: flux_coefficients(:, :)
    end type equilibrated_flux_estimate_3d_t

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(*), info
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
        end subroutine dgesv
    end interface

    !> Local faces of the FortFEM tetrahedral RT dof map (vertices per face).
    integer, parameter :: LOCAL_FACE_VERTICES(3, 4) = &
        reshape([1, 4, 3, 1, 2, 4, 1, 3, 2, 2, 3, 4], [3, 4])

    !> Reference data shared by all cells.
    type :: reference_t
        integer :: degree = 0, nrt = 0, nk = 0, nq = 0, face_dofs = 0
        real(dp), allocatable :: weights(:), points(:, :)   !< (nq), (3, nq)
        real(dp), allocatable :: values(:, :, :), divergences(:, :)
        real(dp), allocatable :: poly(:, :)                 !< (nk, nq)
    end type reference_t

    !> Per-cell geometry and the global RT dof data.
    type :: mesh_t
        real(dp), allocatable :: jacobian(:, :, :), determinant(:)
        real(dp), allocatable :: hat_grad(:, :, :), diameter(:)
        real(dp), allocatable :: transform(:, :, :)        !< (nrt, nrt, nt)
        integer, allocatable :: global_dofs(:, :)          !< (nrt, nt)
        integer, allocatable :: face_orientations(:, :)    !< (4, nt)
        integer, allocatable :: face_permutations(:, :, :) !< (3, 4, nt)
        integer, allocatable :: vertex_start(:), vertex_cell(:), vertex_slot(:)
        logical, allocatable :: boundary_vertex(:)
    end type mesh_t

contains

    subroutine equilibrated_flux_estimate_3d(vertices, tetrahedra, degree, xi, eta, &
            zeta, weights, coefficient, source, discrete_flux, coefficient_lower, &
            friedrichs_constant, estimate, status)
        !! Data arrays are sampled at the affine images of the reference points
        !! (xi, eta, zeta) of each cell (weights sum to 1/6); the rule must
        !! integrate polynomials of degree 2k + 2 exactly. tetrahedra must be
        !! positively oriented (positive determinant). 0 <= degree <= 4.
        real(dp), intent(in) :: vertices(:, :)          !< (3, nv)
        integer, intent(in) :: tetrahedra(:, :)         !< (4, nt)
        integer, intent(in) :: degree                   !< k in 0..4
        real(dp), intent(in) :: xi(:), eta(:), zeta(:), weights(:)
        real(dp), intent(in) :: coefficient(:, :)       !< A, (nq, nt)
        real(dp), intent(in) :: source(:, :)            !< f, (nq, nt)
        real(dp), intent(in) :: discrete_flux(:, :, :)  !< A grad u_h, (3, nq, nt)
        real(dp), intent(in) :: coefficient_lower(:)    !< (nt)
        real(dp), intent(in) :: friedrichs_constant
        type(equilibrated_flux_estimate_3d_t), intent(out) :: estimate
        integer, intent(out) :: status
        type(reference_t) :: ref
        type(mesh_t) :: mesh
        real(dp), allocatable :: projection(:, :), contribution(:, :, :)
        real(dp), allocatable :: divergence_sq(:)
        integer :: nv, nt, a, t, failed

        status = EQUILIBRATION_BAD_INPUT
        nv = size(vertices, 2)
        nt = size(tetrahedra, 2)
        if (degree < 0 .or. degree > 4 .or. nt < 1) return
        if (size(vertices, 1) /= 3 .or. size(tetrahedra, 1) /= 4) return
        if (any(tetrahedra < 1) .or. any(tetrahedra > nv)) return
        if (size(eta) /= size(xi) .or. size(zeta) /= size(xi)) return
        if (size(weights) /= size(xi)) return
        if (size(coefficient, 1) /= size(xi) .or. size(coefficient, 2) /= nt) return
        if (size(source, 1) /= size(xi) .or. size(source, 2) /= nt) return
        if (size(discrete_flux, 1) /= 3 .or. size(discrete_flux, 2) /= size(xi)) return
        if (size(discrete_flux, 3) /= nt) return
        if (size(coefficient_lower) /= nt) return
        if (any(.not. (coefficient_lower > 0.0_dp))) return
        if (any(.not. (coefficient > 0.0_dp))) return
        if (.not. (friedrichs_constant >= 0.0_dp)) return

        call build_reference(degree, xi, eta, zeta, weights, ref, status)
        if (status /= 0) then
            status = EQUILIBRATION_BAD_INPUT
            return
        end if
        call build_mesh(vertices, tetrahedra, degree, ref, mesh, status)
        if (status /= 0) then
            status = EQUILIBRATION_BAD_CELL
            return
        end if

        allocate (projection(ref%nk, nt), estimate%cell_oscillation(nt))
        allocate (contribution(ref%nrt, 4, nt))
        contribution = 0.0_dp
        !$omp parallel do schedule(static)
        do t = 1, nt
            call project_source(ref, mesh, t, source(:, t), coefficient_lower(t), &
                projection(:, t), estimate%cell_oscillation(t))
        end do
        !$omp end parallel do

        failed = 0
        !$omp parallel do schedule(dynamic, 4) reduction(max:failed)
        do a = 1, nv
            if (mesh%vertex_start(a + 1) == mesh%vertex_start(a)) cycle
            call solve_patch(ref, mesh, a, coefficient, source, discrete_flux, &
                contribution, failed)
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
            call cell_defects(ref, mesh, t, coefficient(:, t), &
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
    end subroutine equilibrated_flux_estimate_3d

    subroutine build_reference(degree, xi, eta, zeta, weights, ref, status)
        integer, intent(in) :: degree
        real(dp), intent(in) :: xi(:), eta(:), zeta(:), weights(:)
        type(reference_t), intent(out) :: ref
        integer, intent(out) :: status
        type(tetra_rt_t) :: basis
        integer :: q

        call initialize_tetra_rt(degree, basis, status)
        if (status /= 0) return
        ref%degree = degree
        ref%nrt = basis%dof_count
        ref%nk = monomial_count(degree)
        ref%face_dofs = (degree + 1)*(degree + 2)/2
        ref%nq = size(xi)
        ref%weights = weights
        allocate (ref%points(3, ref%nq), ref%values(3, ref%nrt, ref%nq))
        allocate (ref%divergences(ref%nrt, ref%nq), ref%poly(ref%nk, ref%nq))
        do q = 1, ref%nq
            ref%points(:, q) = [xi(q), eta(q), zeta(q)]
            call evaluate_tetra_rt(basis, ref%points(:, q), ref%values(:, :, q), &
                ref%divergences(:, q), status)
            if (status /= 0) return
            call monomials(degree, ref%points(:, q), ref%poly(:, q))
        end do
        status = 0
    end subroutine build_reference

    subroutine build_mesh(vertices, tetrahedra, degree, ref, mesh, status)
        !! Geometry of every cell, the RT dof transforms and the vertex stars,
        !! and the boundary vertices (vertices of faces met by one cell only).
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree
        type(reference_t), intent(in) :: ref
        type(mesh_t), intent(out) :: mesh
        integer, intent(out) :: status
        integer, allocatable :: faces(:, :), face_count(:)
        real(dp) :: geometry(22), lengths(6), cell_vertices(3, 4)
        integer :: nv, nt, t, s, f, face, nfaces

        status = 1
        nv = size(vertices, 2)
        nt = size(tetrahedra, 2)
        allocate (mesh%jacobian(3, 3, nt), mesh%determinant(nt), mesh%diameter(nt))
        allocate (mesh%hat_grad(3, 4, nt), mesh%transform(ref%nrt, ref%nrt, nt))
        do t = 1, nt
            cell_vertices = vertices(:, tetrahedra(:, t))
            call generated_tetra_simplex_geometry(cell_vertices, geometry)
            mesh%determinant(t) = geometry(10)
            if (.not. (mesh%determinant(t) > 0.0_dp)) return
            mesh%jacobian(:, :, t) = reshape(geometry(1:9), [3, 3])
            do s = 1, 4
                mesh%hat_grad(:, s, t) = geometry(11 + 3*(s - 1):13 + 3*(s - 1))
            end do
            call generated_tetra_edge_lengths_squared(cell_vertices, lengths)
            mesh%diameter(t) = sqrt(maxval(lengths))
        end do

        call build_tetra_rt_dof_map(degree, tetrahedra, faces, mesh%global_dofs, &
            mesh%face_orientations, mesh%face_permutations, status)
        if (status /= 0) return
        do t = 1, nt
            call build_tetra_rt_basis_transform(degree, &
                mesh%face_orientations(:, t), mesh%face_permutations(:, :, t), &
                mesh%transform(:, :, t), status)
            if (status /= 0) return
        end do

        nfaces = size(faces, 2)
        allocate (face_count(nfaces), mesh%boundary_vertex(nv))
        face_count = 0
        do t = 1, nt
            do f = 1, 4
                face = (mesh%global_dofs((f - 1)*ref%face_dofs + 1, t) - 1) &
                    /ref%face_dofs + 1
                face_count(face) = face_count(face) + 1
            end do
        end do
        mesh%boundary_vertex = .false.
        do face = 1, nfaces
            if (face_count(face) == 1) mesh%boundary_vertex(faces(:, face)) = .true.
        end do
        call vertex_stars(nv, tetrahedra, mesh%vertex_start, mesh%vertex_cell, &
            mesh%vertex_slot)
        status = 0
    end subroutine build_mesh

    subroutine project_source(ref, mesh, t, f, a_lower, coefficients, oscillation)
        !! L2 projection of f onto P_k(T) and the oscillation term
        !! (h_T/pi) A_T,min^(-1/2) ||f - Pi_k f||_T.
        type(reference_t), intent(in) :: ref
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: t
        real(dp), intent(in) :: f(:), a_lower
        real(dp), intent(out) :: coefficients(:), oscillation
        real(dp) :: mass(ref%nk, ref%nk), rhs(ref%nk), pair(2), rest, wdet, pf
        real(dp) :: defect(1)
        integer :: q, i, j, info, ipiv(ref%nk)

        mass = 0.0_dp
        rhs = 0.0_dp
        do q = 1, ref%nq
            wdet = ref%weights(q)*mesh%determinant(t)
            do i = 1, ref%nk
                do j = 1, ref%nk
                    call generated_tetra_projection_terms(ref%poly(i, q), &
                        ref%poly(j, q), f(q), wdet, pair)
                    mass(j, i) = mass(j, i) + pair(1)
                end do
                call generated_tetra_projection_terms(ref%poly(i, q), &
                    ref%poly(i, q), f(q), wdet, pair)
                rhs(i) = rhs(i) + pair(2)
            end do
        end do
        call dgesv(ref%nk, 1, mass, ref%nk, ipiv, rhs, ref%nk, info)
        coefficients = rhs
        rest = 0.0_dp
        do q = 1, ref%nq
            wdet = ref%weights(q)*mesh%determinant(t)
            pf = dot_product(coefficients, ref%poly(:, q))
            call generated_tetra_oscillation_defect(f(q), pf, wdet, defect)
            rest = rest + defect(1)
        end do
        oscillation = mesh%diameter(t)/acos(-1.0_dp)*sqrt(rest/a_lower)
        if (info /= 0) oscillation = huge(1.0_dp)
    end subroutine project_source

    subroutine solve_patch(ref, mesh, a, coefficient, source, discrete_flux, &
            contribution, failed)
        !! Dense KKT solve for sigma_a. Unknowns are the distinct global RT dofs
        !! of the patch that are not on faces opposite to a; cell i contributes
        !! Ci (nrt x ns) mapping them to the oriented local dofs.
        type(reference_t), intent(in) :: ref
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: a
        real(dp), intent(in) :: coefficient(:, :), source(:, :)
        real(dp), intent(in) :: discrete_flux(:, :, :)
        real(dp), intent(inout) :: contribution(:, :, :)
        integer, intent(inout) :: failed
        integer :: m, ns, nl, n, i, j, q, k, l, t, s, f, d, p, info
        integer :: fdc, nrt, nk, lam0
        integer, allocatable :: cmap(:, :), gids(:), ipiv(:), cell(:), slot(:)
        real(dp), allocatable :: cmat(:, :, :), kkt(:, :), rhs(:), qint(:, :)
        real(dp), allocatable :: lrhs(:, :), gram(:, :), dpos(:, :), rloc(:)
        real(dp), allocatable :: bblock(:, :), gc(:, :), x(:)
        real(dp) :: phi(3, ref%nrt), dv(ref%nrt), load_values(8), wdet, ainv
        real(dp) :: flux(3), psi, load, pair(1), mult(2)
        real(dp) :: area, mean
        logical :: interior_patch

        m = mesh%vertex_start(a + 1) - mesh%vertex_start(a)
        interior_patch = .not. mesh%boundary_vertex(a)
        nrt = ref%nrt
        nk = ref%nk
        fdc = ref%face_dofs
        allocate (cmap(nrt, m), gids(m*nrt), cell(m), slot(m))
        cmap = 0
        ns = 0
        do i = 1, m
            t = mesh%vertex_cell(mesh%vertex_start(a) + i - 1)
            s = mesh%vertex_slot(mesh%vertex_start(a) + i - 1)
            cell(i) = t
            slot(i) = s
            do f = 1, 4
                if (.not. any(LOCAL_FACE_VERTICES(:, f) == s)) cycle
                do d = 1, fdc
                    j = (f - 1)*fdc + d
                    call find_or_add(mesh%global_dofs(j, t), gids, ns, p)
                    cmap(j, i) = p
                end do
            end do
            do j = 4*fdc + 1, nrt
                call find_or_add(mesh%global_dofs(j, t), gids, ns, p)
                cmap(j, i) = p
            end do
        end do

        nl = m*nk
        n = ns + nl + merge(1, 0, interior_patch)
        lam0 = ns
        allocate (cmat(nrt, ns, m), kkt(n, n), rhs(n), ipiv(n))
        allocate (qint(nk, m), lrhs(nk, m), gram(nrt, nrt), dpos(nrt, nk))
        allocate (rloc(nrt), bblock(nk, ns), gc(nrt, ns), x(ns))
        kkt = 0.0_dp
        rhs = 0.0_dp
        cmat = 0.0_dp
        do i = 1, m
            t = cell(i)
            do j = 1, nrt
                p = cmap(j, i)
                if (p > 0) cmat(:, p, i) = cmat(:, p, i) + mesh%transform(:, j, t)
            end do
        end do

        do i = 1, m
            t = cell(i)
            s = slot(i)
            gram = 0.0_dp
            dpos = 0.0_dp
            rloc = 0.0_dp
            lrhs(:, i) = 0.0_dp
            qint(:, i) = 0.0_dp
            do q = 1, ref%nq
                call physical_rt(ref, mesh, t, q, phi, dv, info)
                if (info /= 0) then
                    failed = 1
                    return
                end if
                wdet = ref%weights(q)*mesh%determinant(t)
                ainv = 1.0_dp/coefficient(q, t)
                flux = discrete_flux(:, q, t)
                call generated_tetra_patch_load(ref%points(:, q), source(q, t), &
                    flux, mesh%hat_grad(:, :, t), load_values)
                load = load_values(s)
                psi = load_values(4 + s)
                do k = 1, nrt
                    call generated_tetra_patch_flux_rhs(phi(:, k), flux, psi, &
                        ainv, wdet, pair)
                    rloc(k) = rloc(k) + pair(1)
                    do j = 1, nrt
                        call generated_tetra_patch_gram(phi(:, j), phi(:, k), &
                            ainv, wdet, pair)
                        gram(j, k) = gram(j, k) + pair(1)
                    end do
                    do l = 1, nk
                        call generated_tetra_patch_multiplier(dv(k), load, &
                            ref%poly(l, q), wdet, mult)
                        dpos(k, l) = dpos(k, l) + mult(1)
                    end do
                end do
                do l = 1, nk
                    call generated_tetra_patch_multiplier(1.0_dp, load, &
                        ref%poly(l, q), wdet, mult)
                    qint(l, i) = qint(l, i) + mult(1)
                    lrhs(l, i) = lrhs(l, i) + mult(2)
                end do
            end do
            ! Oriented local dofs -> patch unknowns: C_i^T (G C_i), C_i^T r_loc,
            ! and the constraint block -C_i^T D (rows of the multipliers).
            gc = matmul(gram, cmat(:, :, i))
            kkt(1:ns, 1:ns) = kkt(1:ns, 1:ns) + matmul(transpose(cmat(:, :, i)), gc)
            rhs(1:ns) = rhs(1:ns) + matmul(transpose(cmat(:, :, i)), rloc)
            bblock = -matmul(transpose(dpos), cmat(:, :, i))
            lam0 = ns + (i - 1)*nk
            kkt(lam0 + 1:lam0 + nk, 1:ns) = bblock
            rhs(lam0 + 1:lam0 + nk) = lrhs(:, i)
        end do
        kkt(1:ns, ns + 1:ns + nl) = transpose(kkt(ns + 1:ns + nl, 1:ns))

        if (interior_patch) then
            ! Remove the patch mean of the load and fix the constant multiplier
            ! with one extra mean constraint (as in the 2D module).
            area = sum(qint(1, :))
            mean = sum(lrhs(1, :))/area
            do i = 1, m
                lam0 = ns + (i - 1)*nk
                rhs(lam0 + 1:lam0 + nk) = rhs(lam0 + 1:lam0 + nk) - mean*qint(:, i)
                kkt(lam0 + 1:lam0 + nk, n) = qint(:, i)
                kkt(n, lam0 + 1:lam0 + nk) = qint(:, i)
            end do
        end if
        call dgesv(n, 1, kkt, n, ipiv, rhs, n, info)
        if (info /= 0) then
            failed = 1
            return
        end if
        x = rhs(1:ns)
        do i = 1, m
            t = cell(i)
            s = slot(i)
            contribution(:, s, t) = matmul(cmat(:, :, i), x)
        end do
        deallocate (cmap, gids, cell, slot, cmat, kkt, rhs, ipiv, qint, lrhs)
        deallocate (gram, dpos, rloc, bblock, gc, x)
    end subroutine solve_patch

    subroutine find_or_add(global, gids, count_, position)
        !! Position of a global dof in the patch list, appending it if new.
        integer, intent(in) :: global
        integer, intent(inout) :: gids(:), count_
        integer, intent(out) :: position
        integer :: k

        position = 0
        do k = 1, count_
            if (gids(k) == global) then
                position = k
                return
            end if
        end do
        count_ = count_ + 1
        gids(count_) = global
        position = count_
    end subroutine find_or_add

    subroutine physical_rt(ref, mesh, t, q, phi, dv, status)
        type(reference_t), intent(in) :: ref
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: t, q
        real(dp), intent(out) :: phi(3, ref%nrt), dv(ref%nrt)
        integer, intent(out) :: status

        call map_tetra_rt_contravariant(mesh%jacobian(:, :, t), ref%values(:, :, q), &
            ref%divergences(:, q), phi, dv, status)
    end subroutine physical_rt

    subroutine cell_defects(ref, mesh, t, coefficient, discrete_flux, projection, &
            sigma, flux_defect, divergence_sq)
        type(reference_t), intent(in) :: ref
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: t
        real(dp), intent(in) :: coefficient(:), discrete_flux(:, :), projection(:)
        real(dp), intent(in) :: sigma(:)
        real(dp), intent(out) :: flux_defect, divergence_sq
        real(dp) :: phi(3, ref%nrt), dv(ref%nrt), wdet, pf, div_sigma
        real(dp) :: defect(1), phys(3)
        integer :: q, stat

        flux_defect = 0.0_dp
        divergence_sq = 0.0_dp
        do q = 1, ref%nq
            call physical_rt(ref, mesh, t, q, phi, dv, stat)
            wdet = ref%weights(q)*mesh%determinant(t)
            phys = matmul(phi, sigma)
            div_sigma = dot_product(dv, sigma)
            pf = dot_product(projection, ref%poly(:, q))
            call generated_tetra_flux_defect(phys, discrete_flux(:, q), &
                coefficient(q), wdet, defect)
            flux_defect = flux_defect + defect(1)
            call generated_tetra_divergence_defect(div_sigma, pf, wdet, defect)
            divergence_sq = divergence_sq + defect(1)
        end do
        flux_defect = sqrt(flux_defect)
    end subroutine cell_defects

    pure integer function monomial_count(degree) result(count_)
        integer, intent(in) :: degree

        count_ = (degree + 1)*(degree + 2)*(degree + 3)/6
    end function monomial_count

    subroutine monomials(degree, point, poly)
        !! Generated monomial basis of P_degree at the reference point.
        integer, intent(in) :: degree
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: poly(:)

        select case (degree)
        case (0)
            call generated_tetra_monomials_degree0(point, poly)
        case (1)
            call generated_tetra_monomials_degree1(point, poly)
        case (2)
            call generated_tetra_monomials_degree2(point, poly)
        case (3)
            call generated_tetra_monomials_degree3(point, poly)
        case default
            call generated_tetra_monomials_degree4(point, poly)
        end select
    end subroutine monomials

end module fortfem_equilibrated_flux_3d
