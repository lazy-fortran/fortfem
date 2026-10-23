module fortfem_parametric_laplace_bem_3d
    !! Scalar Laplace Galerkin BEM on a caller-owned closed parametric surface.
    !! Normals point out of the bounded region. No axisymmetry is assumed.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_generated_bem_frame, only: bem_frame
    use fortfem_generated_bem_kernel, only: bem_kernel
    use fortfem_generated_bem_pair_product, only: bem_pair_product
    use fortfem_generated_bem_curl_density, only: bem_curl_density
    use fortfem_generated_bem_refinement_product, only: bem_refinement_product
    use fortfem_generated_p1_basis, only: generated_p1_jet
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_equilibrated_flux_topology, only: topology_t, build_topology
    use fortnum_linalg, only: dense_solve
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private

    type, abstract, public :: bem_surface_t
    contains
        procedure(surface_evaluate), deferred :: evaluate
    end type
    abstract interface
        pure subroutine surface_evaluate(self, panel, uv, x, a, b, status)
            import :: bem_surface_t, dp
            class(bem_surface_t), intent(in) :: self
            integer, intent(in) :: panel
            real(dp), intent(in) :: uv(2)
            real(dp), intent(out) :: x(3), a(3), b(3)
            integer, intent(out) :: status
        end subroutine
    end interface

    type, public :: bem_operators_t
        real(dp), allocatable :: v(:, :), k(:, :), adjoint(:, :), w(:, :)
        real(dp), allocatable :: mass(:, :), area(:)
        integer :: genus = -1
    end type
    type, public :: bem_exterior_dtn_t
        private
        real(dp), allocatable :: flux_map(:, :), steklov(:, :)
    contains
        procedure :: apply => apply_exterior_dtn
        procedure :: matrix => exterior_dtn_matrix
    end type
    type :: surface_point_t
        real(dp) :: x(3), a(3), b(3), normal(3), jacobian
        real(dp) :: basis(3), curl(3, 3)
    end type
    public :: assemble_parametric_laplace_bem, build_exterior_dtn
    public :: parametric_bem_residual, assignment(=)
    public :: parametric_bem_refinement_estimate
    interface assignment(=)
        module procedure assign_operators, assign_dtn
    end interface
    real(dp), parameter :: scale = 1.0_dp/(4.0_dp*acos(-1.0_dp))
    real(dp), parameter :: corners(2, 3) = reshape( &
        [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
contains
    subroutine assign_operators(lhs, rhs)
        type(bem_operators_t), intent(out) :: lhs
        type(bem_operators_t), intent(in) :: rhs
        if (allocated(rhs%v)) lhs%v = rhs%v
        if (allocated(rhs%k)) lhs%k = rhs%k
        if (allocated(rhs%adjoint)) lhs%adjoint = rhs%adjoint
        if (allocated(rhs%w)) lhs%w = rhs%w
        if (allocated(rhs%mass)) lhs%mass = rhs%mass
        if (allocated(rhs%area)) lhs%area = rhs%area
        lhs%genus = rhs%genus
    end subroutine

    subroutine assign_dtn(lhs, rhs)
        type(bem_exterior_dtn_t), intent(out) :: lhs
        type(bem_exterior_dtn_t), intent(in) :: rhs
        if (allocated(rhs%flux_map)) lhs%flux_map = rhs%flux_map
        if (allocated(rhs%steklov)) lhs%steklov = rhs%steklov
    end subroutine

    subroutine assemble_parametric_laplace_bem(surface, faces, degree, op, status)
        !! P1 continuous potentials / P0 discontinuous fluxes; geometry may
        !! have arbitrary order through evaluate. Edge orientation is checked.
        !! Coincident panels use a target-centred Duffy fan. Close panels use
        !! a fan about a closest point on the source; the radial square map
        !! clusters nodes near the apex. Increase degree to check quadrature.
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: faces(:, :), degree
        type(bem_operators_t), intent(out) :: op
        integer, intent(out) :: status
        type(surface_point_t), allocatable :: points(:, :)
        type(topology_t) :: topo
        real(dp), allocatable :: xi(:), eta(:), weights(:), nodes(:), lw(:)
        real(dp), allocatable :: centers(:, :), diameter(:)
        real(dp) :: vs, ks(3), ws(3, 3), target_weight
        integer, allocatable :: counts(:), signs(:)
        integer :: nv, nf, nq, nl, i, j, t, a, b, e, errors, local_status
        logical :: near

        status = 1
        if (size(faces, 1) /= 3 .or. size(faces, 2) < 4) return
        if (any(faces < 1) .or. degree < 2) return
        nv = maxval(faces)
        nf = size(faces, 2)
        do i = 1, nf
            if (faces(1, i) == faces(2, i) .or. &
                faces(2, i) == faces(3, i) .or. &
                faces(1, i) == faces(3, i)) return
        end do
        call build_topology(nv, faces, topo, local_status)
        if (local_status /= 0) return
        e = maxval(topo%cell_edge)
        allocate(counts(e), signs(e))
        counts = 0
        signs = 0
        do i = 1, nf
            do a = 1, 3
                j = topo%cell_edge(a, i)
                counts(j) = counts(j) + 1
                signs(j) = signs(j) + topo%cell_orientation(a, i)
            end do
        end do
        if (any(counts /= 2) .or. any(signs /= 0)) return
        ! Connectedness via the vertex adjacency, without a dense incidence.
        call connected_surface(faces, nv, local_status)
        if (local_status /= 0) return
        if (mod(2 - nv + e - nf, 2) /= 0) return
        op%genus = (2 - nv + e - nf)/2
        if (op%genus < 0) return
        call triangle_duffy_quadrature(degree, xi, eta, weights, local_status)
        if (local_status /= 0) return
        nq = size(weights)
        nl = max(3, (degree + 3)/2)
        allocate(nodes(nl), lw(nl), points(nq, nf))
        call gauss_legendre_ab(nl, 0.0_dp, 1.0_dp, nodes, lw)
        allocate(centers(3, nf), diameter(nf))
        allocate(op%v(nf, nf), op%k(nf, nv), op%mass(nf, nv), &
            op%w(nv, nv), op%area(nf))
        op%v = 0.0_dp
        op%k = 0.0_dp
        op%w = 0.0_dp
        op%mass = 0.0_dp
        errors = 0
        !$omp parallel do default(shared) private(i,t,local_status) &
        !$omp reduction(+:errors) schedule(static)
        do i = 1, nf
            do t = 1, nq
                call sample_surface(surface, i, [xi(t), eta(t)], &
                    points(t, i), local_status)
                if (local_status /= 0) errors = errors + 1
            end do
        end do
        !$omp end parallel do
        status = 2
        if (errors /= 0) return
        do i = 1, nf
            centers(:, i) = 0.0_dp
            do t = 1, nq
                centers(:, i) = centers(:, i) + points(t, i)%x/real(nq, dp)
            end do
            diameter(i) = 0.0_dp
            do t = 1, nq
                diameter(i) = max(diameter(i), norm2(points(t, i)%x-centers(:, i)))
                op%mass(i, faces(:, i)) = op%mass(i, faces(:, i)) + &
                    weights(t)*points(t, i)%jacobian*points(t, i)%basis
            end do
            op%area(i) = sum(op%mass(i, :))
        end do
        !$omp parallel do default(shared) schedule(dynamic,1) &
        !$omp private(i,j,t,a,b,near,vs,ks,ws,target_weight,local_status) &
        !$omp reduction(+:errors)
        do i = 1, nf
            do j = 1, nf
                near = norm2(centers(:, i)-centers(:, j)) < &
                    1.5_dp*(diameter(i)+diameter(j))
                do t = 1, nq
                    target_weight = weights(t)
                    if (near) then
                        call integrate_source_fan(surface, j, points(t, i), &
                            [xi(t), eta(t)], i == j, nodes, lw, vs, ks, ws, &
                            local_status)
                    else
                        call integrate_source_regular(points(t, i), &
                            points(:, j), weights, vs, ks, ws)
                        local_status = 0
                    end if
                    if (local_status /= 0) errors = errors + 1
                    op%v(i, j) = op%v(i, j) + target_weight*vs
                    op%k(i, faces(:, j)) = op%k(i, faces(:, j)) + &
                        target_weight*ks
                    ! Closed-surface identity: K u = integral d_nG (u_y-u_x)
                    ! - u_x/2. Subtract at each target before accumulation;
                    ! this removes the constant singular contribution exactly.
                    op%k(i, faces(:, i)) = op%k(i, faces(:, i)) - &
                        target_weight*sum(ks)*points(t, i)%basis
                    do b = 1, 3
                        do a = 1, 3
                            !$omp atomic update
                            op%w(faces(a, i), faces(b, j)) = &
                                op%w(faces(a, i), faces(b, j)) + &
                                target_weight*ws(a, b)
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        if (errors /= 0) return
        op%k = op%k - 0.5_dp*op%mass
        ! These Galerkin forms are symmetric; averaging removes directional
        ! quadrature error without altering the double-layer jump relation.
        op%v = 0.5_dp*(op%v + transpose(op%v))
        op%w = 0.5_dp*(op%w + transpose(op%w))
        op%adjoint = transpose(op%k)
        if (.not. all(ieee_is_finite(op%v))) return
        if (.not. all(ieee_is_finite(op%k))) return
        if (.not. all(ieee_is_finite(op%w))) return
        status = 0
    end subroutine

    subroutine connected_surface(faces, nv, status)
        integer, intent(in) :: faces(:, :), nv
        integer, intent(out) :: status
        integer :: labels(nv), i, a, lo, hi
        labels = [(i, i=1,nv)]
        do i = 1, size(faces, 2)
            do a = 2, 3
                lo = min(labels(faces(1, i)), labels(faces(a, i)))
                hi = max(labels(faces(1, i)), labels(faces(a, i)))
                where (labels == hi) labels = lo
            end do
        end do
        status = 1
        if (all(labels == 1)) status = 0
    end subroutine

    pure subroutine sample_surface(surface, panel, uv, point, status)
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: panel
        real(dp), intent(in) :: uv(2)
        type(surface_point_t), intent(out) :: point
        integer, intent(out) :: status
        real(dp) :: gradient(2), hessian(2, 2), norm_a, norm_b
        integer :: a
        call surface%evaluate(panel, uv, point%x, point%a, point%b, status)
        if (status /= 0) return
        status = 2
        if (.not. all(ieee_is_finite(point%x))) return
        if (.not. all(ieee_is_finite(point%a))) return
        if (.not. all(ieee_is_finite(point%b))) return
        norm_a = norm2(point%a)
        norm_b = norm2(point%b)
        if (min(norm_a, norm_b) <= tiny(1.0_dp)) return
        ! Reject a degenerate chart before evaluating the generated division.
        if (abs(dot_product(point%a, point%b)) >= &
            (1.0_dp-64.0_dp*epsilon(1.0_dp))*norm_a*norm_b) return
        call bem_frame(point%a(1), point%a(2), point%a(3), &
            point%b(1), point%b(2), point%b(3), point%jacobian, &
            point%normal(1), point%normal(2), point%normal(3))
        do a = 1, 3
            call generated_p1_jet(a, uv(1), uv(2), point%basis(a), gradient, hessian)
            call bem_curl_density(point%a(1), point%a(2), point%a(3), &
                point%b(1), point%b(2), point%b(3), gradient(1), gradient(2), &
                point%curl(1, a), point%curl(2, a), point%curl(3, a))
        end do
        status = 0
    end subroutine

    pure subroutine pair_integrand(target, source, v, k, w)
        type(surface_point_t), intent(in) :: target, source
        real(dp), intent(out) :: v, k(3), w(3, 3)
        real(dp) :: d(3), green, double, unused_v, unused_k
        integer :: a, b
        d = target%x - source%x
        call bem_kernel(d(1), d(2), d(3), source%normal(1), &
            source%normal(2), source%normal(3), green, double)
        do b = 1, 3
            do a = 1, 3
                call bem_pair_product(scale*green, scale*double, &
                    target%jacobian, source%jacobian, source%basis(b), &
                    target%curl(1,a), target%curl(2,a), target%curl(3,a), &
                    source%curl(1,b), source%curl(2,b), source%curl(3,b), &
                    unused_v, unused_k, w(a,b))
                if (a == 1) k(b) = unused_k
                if (a == 1 .and. b == 1) v = unused_v
            end do
        end do
    end subroutine

    pure subroutine integrate_source_regular(target, source, weights, v, k, w)
        type(surface_point_t), intent(in) :: target, source(:)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(out) :: v, k(3), w(3, 3)
        real(dp) :: vi, ki(3), wi(3, 3)
        integer :: i
        v = 0.0_dp
        k = 0.0_dp
        w = 0.0_dp
        do i = 1, size(weights)
            call pair_integrand(target, source(i), vi, ki, wi)
            v = v + weights(i)*vi
            k = k + weights(i)*ki
            w = w + weights(i)*wi
        end do
    end subroutine

    pure subroutine integrate_source_fan(surface, panel, target, uv, same, &
            nodes, weights, v, k, w, status)
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: panel
        type(surface_point_t), intent(in) :: target
        real(dp), intent(in) :: uv(2), nodes(:), weights(:)
        logical, intent(in) :: same
        real(dp), intent(out) :: v, k(3), w(3, 3)
        integer, intent(out) :: status
        type(surface_point_t) :: source
        real(dp) :: apex(2), e1(2), e2(2), det, local_uv(2), r, factor
        real(dp) :: vi, ki(3), wi(3, 3)
        integer :: edge, i, j
        v = 0.0_dp
        k = 0.0_dp
        w = 0.0_dp
        status = 0
        apex = uv
        if (.not. same) then
            call closest_reference_point(surface, panel, target%x, apex, status)
            if (status /= 0) return
        end if
        do edge = 1, 3
            e1 = corners(:, edge) - apex
            e2 = corners(:, mod(edge, 3)+1) - apex
            det = abs(e1(1)*e2(2)-e1(2)*e2(1))
            if (det < 64.0_dp*epsilon(1.0_dp)) cycle
            do i = 1, size(nodes)
                r = nodes(i)**2
                do j = 1, size(nodes)
                    local_uv = apex + r*((1.0_dp-nodes(j))*e1+nodes(j)*e2)
                    call sample_surface(surface, panel, local_uv, source, status)
                    if (status /= 0) return
                    call pair_integrand(target, source, vi, ki, wi)
                    factor = weights(i)*weights(j)*2.0_dp*nodes(i)**3*det
                    v = v + factor*vi
                    k = k + factor*ki
                    w = w + factor*wi
                end do
            end do
        end do
    end subroutine

    pure subroutine closest_reference_point(surface, panel, x, uv, status)
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: panel
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: uv(2)
        integer, intent(out) :: status
        real(dp) :: p(3), a(3), b(3), d(3), aa, ab, bb, det, step(2), t
        integer :: iteration
        uv = [1.0_dp/3.0_dp, 1.0_dp/3.0_dp]
        do iteration = 1, 8
            call surface%evaluate(panel, uv, p, a, b, status)
            if (status /= 0) return
            status = 2
            if (.not. all(ieee_is_finite(p))) return
            if (.not. all(ieee_is_finite(a))) return
            if (.not. all(ieee_is_finite(b))) return
            d = p-x
            aa = dot_product(a, a)
            ab = dot_product(a, b)
            bb = dot_product(b, b)
            det = aa*bb-ab*ab
            if (det <= tiny(1.0_dp)) then
                status = 2
                return
            end if
            step = [(bb*dot_product(a,d)-ab*dot_product(b,d))/det, &
                (aa*dot_product(b,d)-ab*dot_product(a,d))/det]
            uv = uv-step
            uv = max(uv, 0.0_dp)
            if (sum(uv) > 1.0_dp) then
                t = max(0.0_dp, min(1.0_dp, (uv(1)-uv(2)+1.0_dp)/2.0_dp))
                uv = [t, 1.0_dp-t]
            end if
        end do
        status = 0
    end subroutine

    subroutine build_exterior_dtn(op, exterior, status)
        !! One dense LU with all trace right-hand sides and one dgemm for the
        !! Schur product, both threaded by the linked BLAS/LAPACK. Subsequent
        !! applies use cached maps only. S = W + (M/2-K)^T V^{-1}(M/2-K) is positive
        !! exterior energy; physical outward flux q = -V^{-1}(M/2-K)u.
        type(bem_operators_t), intent(in) :: op
        type(bem_exterior_dtn_t), intent(out) :: exterior
        integer, intent(out) :: status
        real(dp), allocatable :: jump(:, :)
        real(dp) :: mean
        integer :: nf, nv, i, j
        interface
            subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                    beta, c, ldc)
                import :: dp
                character(len=1), intent(in) :: transa, transb
                integer, intent(in) :: m, n, k, lda, ldb, ldc
                real(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
                real(dp), intent(inout) :: c(ldc, *)
            end subroutine
        end interface
        status = 1
        if (.not. allocated(op%v)) return
        if (.not. allocated(op%k)) return
        if (.not. allocated(op%mass)) return
        if (.not. allocated(op%w)) return
        nf = size(op%v, 1)
        nv = size(op%k, 2)
        if (nf < 1 .or. nv < 1) return
        if (any(shape(op%v) /= [nf,nf])) return
        if (any(shape(op%k) /= [nf,nv])) return
        if (any(shape(op%mass) /= [nf,nv])) return
        if (any(shape(op%w) /= [nv,nv])) return
        if (.not. all(ieee_is_finite(op%v))) return
        if (.not. all(ieee_is_finite(op%k))) return
        if (.not. all(ieee_is_finite(op%mass))) return
        if (.not. all(ieee_is_finite(op%w))) return
        jump = 0.5_dp*op%mass-op%k
        allocate(exterior%flux_map(nf,nv))
        call dense_solve(op%v, -jump, exterior%flux_map, status)
        if (status /= 0) then
            deallocate(exterior%flux_map)
            return
        end if
        exterior%steklov = op%w
        call dgemm('T', 'N', nv, nv, nf, -1.0_dp, jump, nf, &
            exterior%flux_map, nf, 1.0_dp, exterior%steklov, nv)
        !$omp parallel do default(shared) private(i, j, mean) schedule(static)
        do j = 1, nv
            do i = 1, j - 1
                mean = 0.5_dp*(exterior%steklov(i, j) + exterior%steklov(j, i))
                exterior%steklov(i, j) = mean
                exterior%steklov(j, i) = mean
            end do
        end do
        !$omp end parallel do
    end subroutine

    subroutine exterior_dtn_matrix(self, matrix, status)
        class(bem_exterior_dtn_t), intent(in) :: self
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status
        status = 1
        if (.not. allocated(self%steklov)) return
        matrix = self%steklov
        status = 0
    end subroutine

    subroutine apply_exterior_dtn(self, trace, flux, load, status)
        class(bem_exterior_dtn_t), intent(in) :: self
        real(dp), intent(in) :: trace(:)
        real(dp), intent(out) :: flux(:), load(:)
        integer, intent(out) :: status
        flux = 0.0_dp
        load = 0.0_dp
        status = 1
        if (.not. allocated(self%flux_map)) return
        if (.not. allocated(self%steklov)) return
        if (size(trace) /= size(self%flux_map,2)) return
        if (size(flux) /= size(self%flux_map,1)) return
        if (size(load) /= size(trace)) return
        if (.not. all(ieee_is_finite(trace))) return
        flux = matmul(self%flux_map, trace)
        load = matmul(self%steklov, trace)
        status = 0
    end subroutine

    subroutine parametric_bem_refinement_estimate(fine, coarse_flux, fine_flux, &
            parent, panel_indicator, energy_estimate, status)
        !! Two-level P0 flux estimator on nested panels of the SAME surface.
        !! Energy: ||q_f - I q_c||_Vf. eta_T^2 = h_T ||q_f - I q_c||_L2(T)^2,
        !! with h_T=sqrt(area(T)). Reliability requires flux-error saturation;
        !! neither quadrature error nor the saturation factor is certified.
        type(bem_operators_t), intent(in) :: fine
        real(dp), intent(in) :: coarse_flux(:), fine_flux(:)
        integer, intent(in) :: parent(:)
        real(dp), intent(out) :: panel_indicator(:), energy_estimate
        integer, intent(out) :: status
        real(dp), allocatable :: difference(:), parent_area(:)
        real(dp) :: square, energy_squared, symmetry_scale, symmetry_defect
        real(dp) :: entry, reverse_entry
        integer :: nf, nc, panel, coarse, row
        logical :: invalid

        status = 1
        panel_indicator = 0.0_dp
        energy_estimate = 0.0_dp
        if (.not. allocated(fine%v)) return
        if (.not. allocated(fine%area)) return
        nf = size(fine_flux)
        nc = size(coarse_flux)
        if (min(nf, nc) < 1) return
        if (size(parent) /= nf .or. size(fine%area) /= nf) return
        if (size(panel_indicator) /= nc) return
        if (any(shape(fine%v) /= [nf, nf])) return
        if (any(parent < 1) .or. any(parent > nc)) return
        if (.not. all(ieee_is_finite(fine%area))) return
        if (any(fine%area <= 0.0_dp)) return
        if (.not. all(ieee_is_finite(fine_flux))) return
        if (.not. all(ieee_is_finite(coarse_flux))) return
        symmetry_scale = tiny(1.0_dp)
        symmetry_defect = 0.0_dp
        invalid = .false.
        !$omp parallel do default(shared) schedule(dynamic,16) &
        !$omp private(panel,row,entry,reverse_entry) &
        !$omp reduction(max:symmetry_scale,symmetry_defect) reduction(.or.:invalid)
        do panel = 1, nf
            do row = 1, panel
                entry = fine%v(row, panel)
                reverse_entry = fine%v(panel, row)
                if (.not. ieee_is_finite(entry) .or. &
                    .not. ieee_is_finite(reverse_entry)) then
                    invalid = .true.
                    cycle
                end if
                symmetry_scale = max(symmetry_scale, abs(entry), abs(reverse_entry))
                symmetry_defect = max(symmetry_defect, abs(entry-reverse_entry))
            end do
        end do
        !$omp end parallel do
        if (invalid) return
        if (symmetry_defect > 128*epsilon(1.0_dp)*symmetry_scale) return
        allocate(difference(nf), parent_area(nc))
        parent_area = 0.0_dp
        !$omp parallel do default(shared) private(panel,coarse,square)
        do panel = 1, nf
            coarse = parent(panel)
            call bem_refinement_product(fine_flux(panel), coarse_flux(coarse), &
                difference(panel), square)
            !$omp atomic update
            panel_indicator(coarse) = panel_indicator(coarse) + &
                fine%area(panel)*square
            !$omp atomic update
            parent_area(coarse) = parent_area(coarse) + fine%area(panel)
        end do
        !$omp end parallel do
        if (any(parent_area <= 0.0_dp)) then
            panel_indicator = 0.0_dp
            return
        end if
        energy_squared = 0.0_dp
        ! Column contractions preserve contiguous storage and need no workspace.
        !$omp parallel do default(shared) private(panel) &
        !$omp reduction(+:energy_squared) schedule(static)
        do panel = 1, nf
            energy_squared = energy_squared + difference(panel)* &
                dot_product(fine%v(:, panel), difference)
        end do
        !$omp end parallel do
        status = 2
        if (.not. ieee_is_finite(energy_squared)) return
        if (energy_squared < 0.0_dp) return
        panel_indicator = sqrt(sqrt(parent_area)*panel_indicator)
        if (.not. all(ieee_is_finite(panel_indicator))) return
        energy_estimate = sqrt(energy_squared)
        status = 0
    end subroutine

    subroutine parametric_bem_residual(op, trace, flux, panel_indicator, status)
        !! Mass-scaled equation defect. A practical quadrature/error indicator
        !! only when op is independently refined; no reliability certificate.
        type(bem_operators_t), intent(in) :: op
        real(dp), intent(in) :: trace(:), flux(:)
        real(dp), intent(out) :: panel_indicator(:)
        integer, intent(out) :: status
        panel_indicator = 0.0_dp
        status = 1
        if (.not. allocated(op%v)) return
        if (.not. allocated(op%k)) return
        if (.not. allocated(op%mass)) return
        if (.not. allocated(op%area)) return
        if (size(trace) /= size(op%k,2)) return
        if (size(flux) /= size(op%v,2)) return
        if (size(panel_indicator) /= size(op%v,1)) return
        if (any(op%area <= 0.0_dp)) return
        if (.not. all(ieee_is_finite(trace))) return
        if (.not. all(ieee_is_finite(flux))) return
        panel_indicator = abs(matmul(op%v,flux) - &
            matmul(op%k-0.5_dp*op%mass,trace))/sqrt(op%area)
        status = 0
    end subroutine
end module
