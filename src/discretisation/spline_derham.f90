module fortfem_spline_derham
    !! Discrete de Rham complex of a C0 multipatch spline section, optionally
    !! extended by a periodic Fourier direction t
    !! (doc/design/multipatch_spline_spaces.md, de Rham complex).
    !!
    !! On each block with B-splines N_i of degree p and the Curry-Schoenberg
    !! derivative splines D_i = p/(t_(i+p+1) - t_(i+1)) N_(i+1,p-1), so that
    !! N_i' = D_(i-1) - D_i, the forms are N_i N_j (0), D_i N_j dxi1 and
    !! N_i D_j dxi2 (1), D_i D_j dxi1 dxi2 (2). Their degrees of freedom are
    !! the nodes, edges and faces of the control mesh; interface edges are
    !! shared with an orientation sign (edges run from the lower to the
    !! higher global node), so grad and curl are incidence matrices.
    !! Along t the nodal basis L_k is the trigonometric cardinal basis at
    !! t_k = (k - 1)/n and the edge basis H_k the trigonometric histopolation
    !! basis on [t_k, t_(k+1)], n = 2M + 1, with L_k' = H_(k-1) - H_k.
    !!
    !! 3D spaces and global order (section index fastest):
    !!   V0 = N x T, V1 = [E x T ; N x T'], V2 = [F x T ; E x T'],
    !!   V3 = F x T';
    !!   grad = [G x I ; I x Dt], curl = [[C x I, 0] ; [I x Dt, -G x I]],
    !!   div = [I x Dt, -C x I].
    !! Proxy components in the parametric frame (xi1, xi2, t): a 1-form has
    !! (a1, a2, a3), a 2-form b1 dxi2 dt + b2 dt dxi1 + b3 dxi1 dxi2; the
    !! E x T' degrees of freedom are b2 on xi1-edges and -b1 on xi2-edges.
    !! The projections use the geometric degrees of freedom at the Greville
    !! points (values, edge integrals, cell integrals), so d Pi = Pi d.
    use fortfem_kinds, only: dp
    use fortfem_gll_line, only: gauss_rule
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_bspline_feec, only: evaluate_bspline_basis
    use fortfem_spline_layout, only: block_layout_t
    use fortfem_spline_space, only: spline_space_t, spline_greville
    use fortfem_spline_eval, only: spline_quadrature_t
    implicit none
    private

    type :: triplet_t
        integer :: count = 0
        integer, allocatable :: row(:), col(:)
        real(dp), allocatable :: val(:)
    contains
        procedure :: init => triplet_init
        procedure :: add => triplet_add
        procedure :: to_csr => triplet_to_csr
    end type triplet_t

    type :: block_ids_t
        integer, allocatable :: node(:, :), edge1(:, :), edge2(:, :), face(:, :)
        integer, allocatable :: sign1(:, :), sign2(:, :)
    end type block_ids_t

    type, public :: spline_derham_t
        integer :: degree = 0
        integer :: block_count = 0
        integer :: node_count = 0, edge_count = 0, face_count = 0
        integer :: plane_count = 1 ! n of the t direction, 1: section only
        logical :: toroidal = .false.
        integer :: dimension(0:3) = 0 ! dofs of V0..V3
        type(spline_space_t) :: space
        type(block_ids_t), allocatable :: ids(:)
        integer, allocatable :: edge_nodes(:, :) ! (2, section edge)
        type(csr_matrix_t) :: grad, curl, div ! section or 3D incidence
    end type spline_derham_t

    abstract interface
        subroutine derham_metric(block, xi, t, weight)
            !! Weight matrix of a mass matrix at (block, xi, t) for the
            !! parametric proxies, e.g. |J| G^-1 for 1-forms and G/|J| for
            !! 2-forms (G = J^T J), |J| for 0-forms and 1/|J| for 3-forms.
            import :: dp
            integer, intent(in) :: block
            real(dp), intent(in) :: xi(2), t
            real(dp), intent(out) :: weight(:, :)
        end subroutine derham_metric
        subroutine derham_jacobian(block, xi, t, jacobian)
            !! jacobian(:, j) = dx/dxi_j (and dx/dt as column 3 in 3D).
            import :: dp
            integer, intent(in) :: block
            real(dp), intent(in) :: xi(2), t
            real(dp), intent(out) :: jacobian(:, :)
        end subroutine derham_jacobian
        subroutine derham_field(block, xi, t, values)
            !! Parametric proxy components of the field at (block, xi, t).
            import :: dp
            integer, intent(in) :: block
            real(dp), intent(in) :: xi(2), t
            real(dp), intent(out) :: values(:)
        end subroutine derham_field
    end interface

    public :: derham_field
    public :: derham_metric
    public :: derham_jacobian
    public :: spline_derham_build
    public :: spline_derham_evaluate
    public :: spline_derham_project
    public :: spline_derham_t_basis
    public :: spline_derham_boundary_dofs
    public :: spline_derham_cell_basis
    public :: spline_derham_tree
    public :: spline_derham_harmonic
    public :: spline_derham_surface
    public :: spline_derham_mass
    public :: spline_derham_interpolation

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface

contains

    subroutine spline_derham_build(space, derham, status, toroidal_modes)
        !! Complex of a C0 section space whose extraction is a signed
        !! injection; with toroidal_modes = M the t direction has 2M + 1
        !! nodes. Status 1: unsupported space (C1, periodic block, THB or
        !! constrained extraction).
        type(spline_space_t), intent(in) :: space
        type(spline_derham_t), intent(out) :: derham
        integer, intent(out) :: status
        integer, intent(in), optional :: toroidal_modes

        integer, allocatable :: head(:), next(:)
        type(triplet_t) :: tg, tc
        integer :: k, i, j, n1, n2, nn, ne, nf, e, nt

        status = 1
        if (space%block_count /= space%layout_block_count) return
        if (any(space%layout%periodic)) return
        do j = 1, space%broken_count
            if (space%broken_ptr(j + 1) - space%broken_ptr(j) /= 1) return
            if (abs(abs(space%broken_weight(space%broken_ptr(j))) - 1.0_dp) > 0.0_dp) &
                return
        end do
        derham%space = space
        derham%degree = space%degree
        derham%block_count = space%block_count
        nn = space%dof_count
        allocate(derham%ids(space%block_count), derham%edge_nodes(2, 64), head(nn), &
            next(64))
        head = 0
        ne = 0
        nf = 0
        call tg%init()
        call tc%init()
        do k = 1, space%block_count
            associate (blk => space%blocks(k), ids => derham%ids(k))
                n1 = blk%n(1)
                n2 = blk%n(2)
                allocate(ids%node(n1, n2), ids%edge1(n1 - 1, n2), ids%sign1(n1 - 1, n2), &
                    ids%edge2(n1, n2 - 1), ids%sign2(n1, n2 - 1), &
                    ids%face(n1 - 1, n2 - 1))
                do j = 1, n2
                    do i = 1, n1
                        ids%node(i, j) = space%broken_dof(space%broken_ptr(blk%offset + &
                            i + n1*(j - 1)))
                    end do
                end do
                do j = 1, n2
                    do i = 1, n1 - 1
                        call find_edge(ids%node(i, j), ids%node(i + 1, j), &
                            ids%edge1(i, j), ids%sign1(i, j))
                    end do
                end do
                do j = 1, n2 - 1
                    do i = 1, n1
                        call find_edge(ids%node(i, j), ids%node(i, j + 1), &
                            ids%edge2(i, j), ids%sign2(i, j))
                    end do
                end do
                do j = 1, n2 - 1
                    do i = 1, n1 - 1
                        nf = nf + 1
                        ids%face(i, j) = nf
                        ! curl: e2(i+1,j) - e2(i,j) - e1(i,j+1) + e1(i,j)
                        call tc%add(nf, ids%edge2(i + 1, j), real(ids%sign2(i + 1, j), dp))
                        call tc%add(nf, ids%edge2(i, j), -real(ids%sign2(i, j), dp))
                        call tc%add(nf, ids%edge1(i, j + 1), -real(ids%sign1(i, j + 1), dp))
                        call tc%add(nf, ids%edge1(i, j), real(ids%sign1(i, j), dp))
                    end do
                end do
            end associate
        end do
        do e = 1, ne
            call tg%add(e, derham%edge_nodes(1, e), -1.0_dp)
            call tg%add(e, derham%edge_nodes(2, e), 1.0_dp)
        end do
        derham%node_count = nn
        derham%edge_count = ne
        derham%face_count = nf
        nt = 1
        if (present(toroidal_modes)) then
            if (toroidal_modes < 0) return
            nt = 2*toroidal_modes + 1
            derham%toroidal = .true.
        end if
        derham%plane_count = nt
        if (.not. derham%toroidal) then
            derham%dimension = [nn, ne, nf, 0]
            call tg%to_csr(ne, derham%grad)
            call tc%to_csr(nf, derham%curl)
        else
            derham%dimension = [nn*nt, (ne + nn)*nt, (nf + ne)*nt, nf*nt]
            call build_toroidal(derham, tg, tc)
        end if
        derham%edge_nodes = derham%edge_nodes(:, :ne)
        status = 0

    contains

        subroutine find_edge(na, nb, edge, sgn)
            integer, intent(in) :: na, nb
            integer, intent(out) :: edge, sgn
            integer :: lo, hi, q
            integer, allocatable :: grown(:, :)
            integer, allocatable :: grown_next(:)

            lo = min(na, nb)
            hi = max(na, nb)
            sgn = merge(1, -1, na < nb)
            q = head(lo)
            do while (q /= 0)
                if (derham%edge_nodes(2, q) == hi) then
                    edge = q
                    return
                end if
                q = next(q)
            end do
            if (ne == size(next)) then
                allocate(grown(2, 2*ne), grown_next(2*ne))
                grown(:, :ne) = derham%edge_nodes(:, :ne)
                grown_next(:ne) = next(:ne)
                call move_alloc(grown, derham%edge_nodes)
                call move_alloc(grown_next, next)
            end if
            ne = ne + 1
            derham%edge_nodes(:, ne) = [lo, hi]
            next(ne) = head(lo)
            head(lo) = ne
            edge = ne
        end subroutine find_edge

    end subroutine spline_derham_build

    subroutine build_toroidal(derham, tg, tc)
        !! 3D incidence from the section incidence and the cyclic difference.
        type(spline_derham_t), intent(inout) :: derham
        type(triplet_t), intent(in) :: tg, tc

        type(triplet_t) :: g3, c3, d3
        integer :: nn, ne, nf, nt, k, q, kn

        nn = derham%node_count
        ne = derham%edge_count
        nf = derham%face_count
        nt = derham%plane_count
        call g3%init()
        call c3%init()
        call d3%init()
        do k = 1, nt
            kn = modulo(k, nt) + 1
            do q = 1, tg%count
                ! G x I: E x T rows
                call g3%add(tg%row(q) + ne*(k - 1), tg%col(q) + nn*(k - 1), tg%val(q))
                ! -G x I on the a_t block of the E x T' rows of curl
                call c3%add(nf*nt + tg%row(q) + ne*(k - 1), &
                    ne*nt + tg%col(q) + nn*(k - 1), -tg%val(q))
            end do
            do q = 1, nn
                ! I x Dt: N x T' rows, t-edge k from node k to k + 1
                call g3%add(ne*nt + q + nn*(k - 1), q + nn*(kn - 1), 1.0_dp)
                call g3%add(ne*nt + q + nn*(k - 1), q + nn*(k - 1), -1.0_dp)
            end do
            do q = 1, tc%count
                call c3%add(tc%row(q) + nf*(k - 1), tc%col(q) + ne*(k - 1), tc%val(q))
                call d3%add(tc%row(q) + nf*(k - 1), nf*nt + tc%col(q) + ne*(k - 1), &
                    -tc%val(q))
            end do
            do q = 1, ne
                call c3%add(nf*nt + q + ne*(k - 1), q + ne*(kn - 1), 1.0_dp)
                call c3%add(nf*nt + q + ne*(k - 1), q + ne*(k - 1), -1.0_dp)
            end do
            do q = 1, nf
                call d3%add(q + nf*(k - 1), q + nf*(kn - 1), 1.0_dp)
                call d3%add(q + nf*(k - 1), q + nf*(k - 1), -1.0_dp)
            end do
        end do
        call g3%to_csr(derham%dimension(1), derham%grad)
        call c3%to_csr(derham%dimension(2), derham%curl)
        call d3%to_csr(derham%dimension(3), derham%div)
    end subroutine build_toroidal

    subroutine spline_derham_tree(derham, tree, order, parent_edge, roots)
        !! Spanning forest of the graph of grad (V0 dofs as vertices, V1 dofs
        !! as edges) by breadth-first search: tree(e) marks tree edges. The
        !! remaining edges form the cotree. order lists the vertices in
        !! search order, parent_edge(v) is the tree edge to v (0 at roots).
        !! With roots (e.g. the boundary nodes of Dirichlet n x A data) the
        !! search starts from all of them at once, as one super-root: no tree
        !! edge joins two roots, and the tree has (vertices - roots) edges on
        !! a connected graph.
        type(spline_derham_t), intent(in) :: derham
        logical, allocatable, intent(out) :: tree(:)
        integer, allocatable, intent(out), optional :: order(:), parent_edge(:)
        integer, intent(in), optional :: roots(:)

        type(csr_matrix_t) :: incidence
        integer, allocatable :: queue(:), ends(:, :), parent(:)
        logical, allocatable :: seen(:)
        integer :: nv, ne, e, k, head, tail, v, root, w, r

        nv = derham%dimension(0)
        ne = derham%dimension(1)
        allocate(ends(2, ne))
        do e = 1, ne
            ends(:, e) = 0
            do k = derham%grad%rowptr(e), derham%grad%rowptr(e + 1) - 1
                if (derham%grad%val(k) < 0.0_dp) then
                    ends(1, e) = derham%grad%colidx(k)
                else
                    ends(2, e) = derham%grad%colidx(k)
                end if
            end do
        end do
        call node_edges(ends, nv, incidence)
        allocate(tree(ne), seen(nv), queue(nv), parent(nv))
        tree = .false.
        seen = .false.
        parent = 0
        tail = 0
        if (present(roots)) then
            do r = 1, size(roots)
                if (seen(roots(r))) cycle
                seen(roots(r)) = .true.
                tail = tail + 1
                queue(tail) = roots(r)
            end do
            head = 1
            call search()
        end if
        do root = 1, nv
            if (seen(root)) cycle
            seen(root) = .true.
            tail = tail + 1
            head = tail
            queue(tail) = root
            call search()
        end do
        if (present(order)) order = queue
        if (present(parent_edge)) parent_edge = parent

    contains

        subroutine search()
            do while (head <= tail)
                v = queue(head)
                head = head + 1
                do r = incidence%rowptr(v), incidence%rowptr(v + 1) - 1
                    e = incidence%colidx(r)
                    w = sum(ends(:, e)) - v
                    if (seen(w)) cycle
                    seen(w) = .true.
                    tree(e) = .true.
                    parent(w) = e
                    tail = tail + 1
                    queue(tail) = w
                end do
            end do
        end subroutine search

    end subroutine spline_derham_tree

    subroutine node_edges(ends, nv, incidence)
        !! Edges at each vertex (CSR, values unused).
        integer, intent(in) :: ends(:, :), nv
        type(csr_matrix_t), intent(out) :: incidence

        integer, allocatable :: next(:)
        integer :: e, j, v

        allocate(incidence%rowptr(nv + 1), next(nv + 1))
        incidence%rowptr = 0
        do e = 1, size(ends, 2)
            do j = 1, 2
                incidence%rowptr(ends(j, e) + 1) = incidence%rowptr(ends(j, e) + 1) + 1
            end do
        end do
        incidence%rowptr(1) = 1
        do v = 2, nv + 1
            incidence%rowptr(v) = incidence%rowptr(v) + incidence%rowptr(v - 1)
        end do
        allocate(incidence%colidx(2*size(ends, 2)))
        next = incidence%rowptr
        do e = 1, size(ends, 2)
            do j = 1, 2
                v = ends(j, e)
                incidence%colidx(next(v)) = e
                next(v) = next(v) + 1
            end do
        end do
        incidence%n = nv
    end subroutine node_edges

    subroutine spline_derham_harmonic(derham, harmonic, count)
        !! Representatives h(:, 1:count) of the first cohomology: closed
        !! 1-forms (curl h = 0 exactly) that are not gradients, one per
        !! independent loop of the domain. They vanish on the tree edges.
        !! Section: every cotree edge left over by a dual spanning tree of
        !! the faces (rooted at the outside) generates one class, completed
        !! leaf-first along the dual tree. 3D (section x circle): the
        !! section classes on every t-node plus the toroidal class (one on
        !! every t-edge).
        type(spline_derham_t), intent(in) :: derham
        real(dp), allocatable, intent(out) :: harmonic(:, :)
        integer, intent(out) :: count

        real(dp), allocatable :: section(:, :)
        integer :: ns, nn, ne, nt, k, j

        call section_harmonic(derham, section, ns)
        if (.not. derham%toroidal) then
            harmonic = section
            count = ns
            return
        end if
        nn = derham%node_count
        ne = derham%edge_count
        nt = derham%plane_count
        count = ns + 1
        allocate(harmonic(derham%dimension(1), count))
        harmonic = 0.0_dp
        do j = 1, ns
            do k = 1, nt
                harmonic(ne*(k - 1) + 1:ne*k, j) = section(:, j)
            end do
        end do
        harmonic(ne*nt + 1:ne*nt + nn*nt, count) = 1.0_dp
        call vanish_on_tree(derham, harmonic)
    end subroutine spline_derham_harmonic

    subroutine vanish_on_tree(derham, harmonic)
        !! h <- h - grad phi with phi integrated along the tree, so that h is
        !! zero on the tree edges (same class, still closed).
        type(spline_derham_t), intent(in) :: derham
        real(dp), intent(inout) :: harmonic(:, :)

        logical, allocatable :: tree(:)
        integer, allocatable :: order(:), parent(:)
        real(dp), allocatable :: phi(:), gphi(:)
        integer :: j, q, v, e, k, r

        call spline_derham_tree(derham, tree, order, parent)
        allocate(phi(derham%dimension(0)), gphi(derham%dimension(1)))
        do j = 1, size(harmonic, 2)
            phi = 0.0_dp
            do q = 1, size(order)
                v = order(q)
                e = parent(v)
                if (e == 0) cycle
                do k = derham%grad%rowptr(e), derham%grad%rowptr(e + 1) - 1
                    if (derham%grad%colidx(k) == v) r = k
                end do
                ! (grad phi)(e) = h(e): phi(v) from the other end of e
                phi(v) = 0.0_dp
                do k = derham%grad%rowptr(e), derham%grad%rowptr(e + 1) - 1
                    if (k /= r) phi(v) = phi(v) - derham%grad%val(k)* &
                        phi(derham%grad%colidx(k))
                end do
                phi(v) = (phi(v) + harmonic(e, j))/derham%grad%val(r)
            end do
            gphi = 0.0_dp
            do e = 1, derham%dimension(1)
                do k = derham%grad%rowptr(e), derham%grad%rowptr(e + 1) - 1
                    gphi(e) = gphi(e) + derham%grad%val(k)*phi(derham%grad%colidx(k))
                end do
            end do
            harmonic(:, j) = harmonic(:, j) - gphi
        end do
    end subroutine vanish_on_tree

    subroutine spline_derham_surface(derham, tag, edges, faces, surface_curl, &
            harmonic)
        !! Boundary surface of the 3D complex on the sides with the given tag
        !! (boundary curve x circle): its 1-form dofs `edges` and 2-form dofs
        !! `faces` (indices into V1 and V2), the surface curl (rows `faces`,
        !! columns `edges`, the restriction of curl) and two closed,
        !! non-exact surface 1-forms in V1 indexing: the toroidal loop (one on
        !! the boundary t-edges) and the poloidal loop (one on a cut edge of
        !! the boundary curve on every t-node). The poloidal class belongs to
        !! the boundary curve of the first tagged edge.
        type(spline_derham_t), intent(in) :: derham
        integer, intent(in) :: tag
        integer, allocatable, intent(out) :: edges(:), faces(:)
        type(csr_matrix_t), intent(out) :: surface_curl
        real(dp), allocatable, intent(out) :: harmonic(:, :)

        type(triplet_t) :: tc
        integer, allocatable :: column(:), nodes(:)
        integer :: r, k, ne, nt, cut, j

        call spline_derham_boundary_dofs(derham, tag, 1, edges)
        call spline_derham_boundary_dofs(derham, tag, 2, faces)
        allocate(column(derham%dimension(1)))
        column = 0
        column(edges) = [(j, j = 1, size(edges))]
        call tc%init()
        do r = 1, size(faces)
            do k = derham%curl%rowptr(faces(r)), derham%curl%rowptr(faces(r) + 1) - 1
                if (column(derham%curl%colidx(k)) /= 0) call tc%add(r, &
                    column(derham%curl%colidx(k)), derham%curl%val(k))
            end do
        end do
        call tc%to_csr(size(faces), surface_curl)
        ne = derham%edge_count
        nt = derham%plane_count
        allocate(harmonic(derham%dimension(1), 2))
        harmonic = 0.0_dp
        call spline_derham_boundary_dofs(derham, tag, 0, nodes)
        do j = 1, size(nodes)
            ! t-edge of node n at t-node k has index ne nt + n + nn (k - 1)
            harmonic(ne*nt + nodes(j), 1) = 1.0_dp
        end do
        cut = edges(1)
        do k = 1, nt
            harmonic(cut + ne*(k - 1), 2) = 1.0_dp
        end do
    end subroutine spline_derham_surface

    subroutine spline_derham_mass(derham, quad, form, metric, matrix, status, &
            points)
        !! Mass matrix of V_form, int proxy_i . W proxy_j dxi1 dxi2 (dt) over
        !! the parametric cells (reference measure, without the section map
        !! Jacobian), with the weight W from `metric`, the section quadrature
        !! `quad` and `points` Gauss points per t-interval (default 2 M + 2).
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        integer, intent(in) :: form
        procedure(derham_metric) :: metric
        type(csr_matrix_t), intent(out) :: matrix
        integer, intent(out) :: status
        integer, intent(in), optional :: points

        type(triplet_t) :: tm
        real(dp), allocatable :: tx(:), tw(:), vals(:, :), w(:, :)
        real(dp), allocatable :: ev(:, :, :), nv(:, :, :), fv(:, :, :), m(:, :)
        integer, allocatable :: dofs(:), sgn(:), ed(:), es(:), nd(:), ns(:)
        integer, allocatable :: fd(:), fs(:)
        real(dp) :: lt(derham%plane_count, 2), ht(derham%plane_count, 2), wt, tt
        integer :: nt, nc, cell, q, j, kt, nq, ncomp, a, b, ne, nn, nf, point

        status = 1
        if (form < 0 .or. form > 3) return
        if (form == 3 .and. .not. derham%toroidal) return
        nt = derham%plane_count
        ne = derham%edge_count
        nn = derham%node_count
        nf = derham%face_count
        nq = 2*nt
        if (present(points)) nq = points
        if (derham%toroidal) then
            call gauss_rule(nq, tx, tw, status)
            if (status /= 0) return
            tx = 0.5_dp*(tx + 1.0_dp)/nt
            tw = 0.5_dp*tw/nt
        else
            tx = [0.0_dp]
            tw = [1.0_dp]
        end if
        ncomp = merge(3, merge(2, 1, form == 1), derham%toroidal .and. &
            (form == 1 .or. form == 2))
        call tm%init()
        allocate(w(ncomp, ncomp))
        do cell = 1, derham%space%cell_count
            call spline_derham_cell_basis(derham, quad, cell, 0, nv, nd, ns)
            call spline_derham_cell_basis(derham, quad, cell, 1, ev, ed, es)
            call spline_derham_cell_basis(derham, quad, cell, 2, fv, fd, fs)
            call local_dofs(dofs, sgn)
            nc = size(dofs)
            allocate(m(nc, nc), vals(ncomp, nc))
            m = 0.0_dp
            do q = 1, quad%points_per_cell
                point = quad%point_offset(cell) + q - 1
                do kt = 1, merge(nt, 1, derham%toroidal)
                    do j = 1, size(tx)
                        tt = (kt - 1)/real(nt, dp) + tx(j)
                        wt = reference_weight(cell, q)*tw(j)
                        if (derham%toroidal) call spline_derham_t_basis(nt, tt, lt, ht)
                        call local_values(q, vals)
                        call metric(quad%block(point), quad%xi(:, point), tt, w)
                        m = m + wt*matmul(transpose(vals), matmul(w, vals))
                    end do
                end do
            end do
            do b = 1, nc
                do a = 1, nc
                    if (m(a, b) /= 0.0_dp) call tm%add(dofs(a), dofs(b), &
                        sgn(a)*sgn(b)*m(a, b))
                end do
            end do
            deallocate(m, vals)
        end do
        call tm%to_csr(derham%dimension(form), matrix)
        status = 0

    contains

        real(dp) function reference_weight(cell, q) result(w0)
            !! Gauss weight of point q of the cell in xi (no det DG).
            integer, intent(in) :: cell, q
            integer :: k, e(2), q1, q2, nq1

            k = derham%space%cell_block(cell)
            e = derham%space%cell_element(:, cell)
            nq1 = quad%points_per_direction
            q1 = modulo(q - 1, nq1) + 1
            q2 = (q - 1)/nq1 + 1
            w0 = (derham%space%blocks(k)%breaks1(e(1) + 1) - &
                derham%space%blocks(k)%breaks1(e(1)))* &
                (derham%space%blocks(k)%breaks2(e(2) + 1) - &
                derham%space%blocks(k)%breaks2(e(2)))*quad%t_weight(q1)*quad%t_weight(q2)
        end function reference_weight

        subroutine local_dofs(dofs, sgn)
            !! Global dofs and signs of the cell functions of V_form.
            integer, allocatable, intent(out) :: dofs(:), sgn(:)
            integer :: k

            allocate(dofs(0), sgn(0))
            do k = 1, nt
                select case (form)
                case (0)
                    dofs = [dofs, nd + nn*(k - 1)]
                    sgn = [sgn, ns]
                case (1)
                    dofs = [dofs, ed + ne*(k - 1)]
                    sgn = [sgn, es]
                case (2)
                    if (derham%toroidal) then
                        dofs = [dofs, fd + nf*(k - 1)]
                    else
                        dofs = [dofs, fd]
                    end if
                    sgn = [sgn, fs]
                case default
                    dofs = [dofs, fd + nf*(k - 1)]
                    sgn = [sgn, fs]
                end select
            end do
            if (derham%toroidal .and. form == 1) then
                do k = 1, nt
                    dofs = [dofs, ne*nt + nd + nn*(k - 1)]
                    sgn = [sgn, ns]
                end do
            end if
            if (derham%toroidal .and. form == 2) then
                do k = 1, nt
                    dofs = [dofs, nf*nt + ed + ne*(k - 1)]
                    sgn = [sgn, es]
                end do
            end if
        end subroutine local_dofs

        subroutine local_values(q, vals)
            !! Proxy components of the cell functions at point q (and t).
            integer, intent(in) :: q
            real(dp), intent(out) :: vals(:, :)
            integer :: k, n, s

            vals = 0.0_dp
            n = 0
            do k = 1, nt
                select case (form)
                case (0)
                    do s = 1, size(nd)
                        n = n + 1
                        vals(1, n) = nv(1, s, q)*merge(lt(k, 1), 1.0_dp, &
                            derham%toroidal)
                    end do
                case (1)
                    do s = 1, size(ed)
                        n = n + 1
                        vals(1:2, n) = ev(:, s, q)*merge(lt(k, 1), 1.0_dp, &
                            derham%toroidal)
                    end do
                case (2)
                    do s = 1, size(fd)
                        n = n + 1
                        vals(ncomp, n) = fv(1, s, q)*merge(lt(k, 1), 1.0_dp, &
                            derham%toroidal)
                    end do
                case default
                    do s = 1, size(fd)
                        n = n + 1
                        vals(1, n) = fv(1, s, q)*ht(k, 1)
                    end do
                end select
            end do
            if (derham%toroidal .and. form == 1) then
                do k = 1, nt
                    do s = 1, size(nd)
                        n = n + 1
                        vals(3, n) = nv(1, s, q)*ht(k, 1)
                    end do
                end do
            end if
            if (derham%toroidal .and. form == 2) then
                ! E x T': b2 on xi1-edges, -b1 on xi2-edges
                do k = 1, nt
                    do s = 1, size(ed)
                        n = n + 1
                        vals(2, n) = ev(1, s, q)*ht(k, 1)
                        vals(1, n) = -ev(2, s, q)*ht(k, 1)
                    end do
                end do
            end if
        end subroutine local_values

    end subroutine spline_derham_mass

    subroutine spline_derham_interpolation(derham, jacobian, interpolation, status, &
            tolerance)
        !! Interpolation of nodal vector splines into V1 for auxiliary-space
        !! (Hiptmair-Xu) preconditioners: column node + n0 (d - 1) holds the
        !! commuting projection of the 1-form of phi_node e_d (proxy
        !! J^T e_d phi_node), d = 1..dim over Cartesian components. Per block
        !! the inverses of the collocation and histopolation matrices are
        !! truncated at `tolerance` (default 1e-6) relative to their largest
        !! entry, so the matrix stays sparse (1e-4 suffices for
        !! preconditioning).
        type(spline_derham_t), intent(in) :: derham
        procedure(derham_jacobian) :: jacobian
        type(csr_matrix_t), intent(out) :: interpolation
        integer, intent(out) :: status
        real(dp), intent(in), optional :: tolerance

        type(triplet_t) :: tp
        real(dp), allocatable :: x(:), w(:), g1(:), g2(:), mn1(:, :), md1(:, :)
        real(dp), allocatable :: mn2(:, :), md2(:, :), in1(:, :), id1(:, :), in2(:, :)
        real(dp), allocatable :: id2(:, :), n1g(:, :), n2g(:, :), q1(:, :, :, :)
        real(dp), allocatable :: q2(:, :, :, :), q3(:, :, :, :), s(:, :), c(:, :)
        real(dp), allocatable :: nv(:), dv(:)
        integer, allocatable :: owner(:)
        real(dp) :: tol, tk, h, jac(3, 3), lt(derham%plane_count, 2)
        real(dp) :: ht(derham%plane_count, 2), y
        integer :: dimx, nn, ne, nt, n0, k, kt, n1, n2, a, b, d, i, j, l, row, col
        integer :: nq, seg, m, qq

        status = 1
        tol = 1.0e-6_dp
        if (present(tolerance)) tol = tolerance
        dimx = merge(3, 2, derham%toroidal)
        nn = derham%node_count
        ne = derham%edge_count
        nt = derham%plane_count
        n0 = derham%dimension(0)
        h = 1.0_dp/nt
        nq = derham%degree + 3
        call gauss_rule(nq, x, w, status)
        if (status /= 0) return
        x = 0.5_dp*(x + 1.0_dp)
        w = 0.5_dp*w
        allocate(owner(derham%dimension(1)))
        owner = 0
        call tp%init()
        do k = 1, derham%block_count
            associate (blk => derham%space%blocks(k), ids => derham%ids(k))
                n1 = blk%n(1)
                n2 = blk%n(2)
                allocate(g1(n1), g2(n2))
                call spline_greville(blk%knots1, derham%degree, g1)
                call spline_greville(blk%knots2, derham%degree, g2)
                call direction_matrices(blk%knots1, derham%degree, g1, x, w, mn1, md1)
                call direction_matrices(blk%knots2, derham%degree, g2, x, w, mn2, md2)
                in1 = truncated_inverse(mn1, tol)
                id1 = truncated_inverse(md1, tol)
                in2 = truncated_inverse(mn2, tol)
                id2 = truncated_inverse(md2, tol)
                n1g = mn1
                n2g = mn2
                do kt = 1, nt
                    tk = (kt - 1)*h
                    ! q1(i', a, d, j'): int over xi1-edge i' at xi2 = g2(j')
                    allocate(q1(n1 - 1, n1, dimx, n2), q2(n2 - 1, n2, dimx, n1))
                    q1 = 0.0_dp
                    q2 = 0.0_dp
                    do j = 1, n2
                        do i = 1, n1 - 1
                            call edge_rule(blk%breaks1, g1(i), g1(i + 1), seg)
                            do m = 1, seg
                                do qq = 1, nq
                                    y = segment_point(blk%breaks1, g1(i), g1(i + 1), m, &
                                        x(qq))
                                    call direction_basis(blk%knots1, derham%degree, y, &
                                        nv, dv)
                                    call jacobian(k, [y, g2(j)], tk, jac(:dimx, :dimx))
                                    do d = 1, dimx
                                        q1(i, :, d, j) = q1(i, :, d, j) + &
                                            segment_weight(blk%breaks1, g1(i), &
                                            g1(i + 1), m, w(qq))*nv*jac(d, 1)
                                    end do
                                end do
                            end do
                        end do
                    end do
                    do i = 1, n1
                        do j = 1, n2 - 1
                            call edge_rule(blk%breaks2, g2(j), g2(j + 1), seg)
                            do m = 1, seg
                                do qq = 1, nq
                                    y = segment_point(blk%breaks2, g2(j), g2(j + 1), m, &
                                        x(qq))
                                    call direction_basis(blk%knots2, derham%degree, y, &
                                        nv, dv)
                                    call jacobian(k, [g1(i), y], tk, jac(:dimx, :dimx))
                                    do d = 1, dimx
                                        q2(j, :, d, i) = q2(j, :, d, i) + &
                                            segment_weight(blk%breaks2, g2(j), &
                                            g2(j + 1), m, w(qq))*nv*jac(d, 2)
                                    end do
                                end do
                            end do
                        end do
                    end do
                    if (derham%toroidal) then
                        allocate(q3(n1, n2, nt, dimx))
                        q3 = 0.0_dp
                        do l = 1, nt
                            do qq = 1, nq
                                call spline_derham_t_basis(nt, (l - 1)*h + h*x(qq), lt, ht)
                                do j = 1, n2
                                    do i = 1, n1
                                        call jacobian(k, [g1(i), g2(j)], &
                                            (l - 1)*h + h*x(qq), jac)
                                        q3(i, j, l, :) = q3(i, j, l, :) + h*w(qq)* &
                                            lt(kt, 1)*jac(:, 3)
                                    end do
                                end do
                            end do
                        end do
                    end if
                    do d = 1, dimx
                        do b = 1, n2
                            do a = 1, n1
                                col = ids%node(a, b) + nn*(kt - 1) + n0*(d - 1)
                                ! xi1-edges
                                s = q1(:, a, d, :)*spread(n2g(:, b), 1, n1 - 1)
                                c = matmul(id1, matmul(s, transpose(in2)))
                                call emit(c, ids%edge1, ids%sign1, ne*(kt - 1))
                                ! xi2-edges
                                s = transpose(q2(:, b, d, :))*spread(n1g(:, a), 2, &
                                    n2 - 1)
                                c = matmul(in1, matmul(s, transpose(id2)))
                                call emit(c, ids%edge2, ids%sign2, ne*(kt - 1))
                                if (derham%toroidal) then
                                    do l = 1, nt
                                        s = spread(n1g(:, a), 2, n2)*spread(n2g(:, b), &
                                            1, n1)*q3(:, :, l, d)
                                        c = matmul(in1, matmul(s, transpose(in2)))
                                        call emit_nodes(c, ids%node, ne*nt + nn*(l - 1))
                                    end do
                                end if
                            end do
                        end do
                    end do
                    deallocate(q1, q2)
                    if (allocated(q3)) deallocate(q3)
                end do
                deallocate(g1, g2)
            end associate
        end do
        call tp%to_csr(derham%dimension(1), interpolation)
        status = 0

    contains

        subroutine emit(c, id, sgn, offset)
            real(dp), intent(in) :: c(:, :)
            integer, intent(in) :: id(:, :), sgn(:, :), offset
            integer :: i0, j0, r

            do j0 = 1, size(c, 2)
                do i0 = 1, size(c, 1)
                    if (abs(c(i0, j0)) <= tol) cycle
                    r = offset + id(i0, j0)
                    if (owner(r) == 0) owner(r) = k
                    if (owner(r) /= k) cycle
                    call tp%add(r, col, sgn(i0, j0)*c(i0, j0))
                end do
            end do
        end subroutine emit

        subroutine emit_nodes(c, id, offset)
            real(dp), intent(in) :: c(:, :)
            integer, intent(in) :: id(:, :), offset
            integer :: i0, j0, r

            do j0 = 1, size(c, 2)
                do i0 = 1, size(c, 1)
                    if (abs(c(i0, j0)) <= tol) cycle
                    r = offset + id(i0, j0)
                    if (owner(r) == 0) owner(r) = k
                    if (owner(r) /= k) cycle
                    call tp%add(r, col, c(i0, j0))
                end do
            end do
        end subroutine emit_nodes

    end subroutine spline_derham_interpolation

    function truncated_inverse(a, tol) result(b)
        real(dp), intent(in) :: a(:, :), tol
        real(dp), allocatable :: b(:, :), m(:, :)
        integer, allocatable :: pivot(:)
        integer :: i, info

        m = a
        allocate(b(size(a, 1), size(a, 1)), pivot(size(a, 1)))
        b = 0.0_dp
        do i = 1, size(a, 1)
            b(i, i) = 1.0_dp
        end do
        call dgesv(size(a, 1), size(a, 1), m, size(a, 1), pivot, b, size(a, 1), info)
        where (abs(b) < tol*maxval(abs(b))) b = 0.0_dp
    end function truncated_inverse

    pure subroutine edge_rule(breaks, a, b, count)
        real(dp), intent(in) :: breaks(:), a, b
        integer, intent(out) :: count

        count = 1 + count_inside(breaks, a, b)
    end subroutine edge_rule

    pure integer function count_inside(breaks, a, b) result(n)
        real(dp), intent(in) :: breaks(:), a, b

        n = count(breaks > a .and. breaks < b)
    end function count_inside

    pure real(dp) function segment_point(breaks, a, b, m, x) result(y)
        !! Point x in [0, 1] of the m-th knot span piece of [a, b].
        real(dp), intent(in) :: breaks(:), a, b, x
        integer, intent(in) :: m
        real(dp) :: lo, hi

        call segment_ends(breaks, a, b, m, lo, hi)
        y = lo + (hi - lo)*x
    end function segment_point

    pure real(dp) function segment_weight(breaks, a, b, m, w) result(y)
        real(dp), intent(in) :: breaks(:), a, b, w
        integer, intent(in) :: m
        real(dp) :: lo, hi

        call segment_ends(breaks, a, b, m, lo, hi)
        y = (hi - lo)*w
    end function segment_weight

    pure subroutine segment_ends(breaks, a, b, m, lo, hi)
        real(dp), intent(in) :: breaks(:), a, b
        integer, intent(in) :: m
        real(dp), intent(out) :: lo, hi
        real(dp), allocatable :: cuts(:)

        cuts = [a, pack(breaks, breaks > a .and. breaks < b), b]
        lo = cuts(m)
        hi = cuts(m + 1)
    end subroutine segment_ends

    subroutine section_harmonic(derham, harmonic, count)
        type(spline_derham_t), intent(in) :: derham
        real(dp), allocatable, intent(out) :: harmonic(:, :)
        integer, intent(out) :: count

        type(spline_derham_t) :: plane
        type(csr_matrix_t) :: edge_faces
        logical, allocatable :: tree(:), dual(:), seen(:)
        integer, allocatable :: queue(:), parent_edge(:), free(:), face_ends(:, :)
        real(dp), allocatable :: value(:), sgn(:, :)
        integer :: ne, nf, e, k, f, g, head, tail, r, j, q, root
        real(dp) :: s

        plane = derham
        if (derham%toroidal) then
            call section_incidence(derham, plane)
        end if
        ne = plane%edge_count
        nf = plane%face_count
        call spline_derham_tree(plane, tree)
        ! faces at each cotree edge: dual vertices 1..nf, outside nf + 1
        allocate(face_ends(2, ne), sgn(2, ne))
        face_ends = nf + 1
        sgn = 0.0_dp
        do f = 1, nf
            do k = plane%curl%rowptr(f), plane%curl%rowptr(f + 1) - 1
                e = plane%curl%colidx(k)
                j = merge(1, 2, face_ends(1, e) == nf + 1 .and. sgn(1, e) == 0.0_dp)
                face_ends(j, e) = f
                sgn(j, e) = plane%curl%val(k)
            end do
        end do
        do e = 1, ne
            if (tree(e)) face_ends(:, e) = [nf + 2, nf + 2]
        end do
        call node_edges(min(face_ends, nf + 1), nf + 1, edge_faces)
        allocate(dual(ne), seen(nf + 1), queue(nf + 1), parent_edge(nf + 1))
        dual = .false.
        seen = .false.
        parent_edge = 0
        root = nf + 1
        seen(root) = .true.
        head = 1
        tail = 1
        queue(1) = root
        do while (head <= tail)
            f = queue(head)
            head = head + 1
            do r = edge_faces%rowptr(f), edge_faces%rowptr(f + 1) - 1
                e = edge_faces%colidx(r)
                if (tree(e)) cycle
                g = sum(face_ends(:, e)) - f
                if (seen(g)) cycle
                seen(g) = .true.
                dual(e) = .true.
                parent_edge(g) = e
                tail = tail + 1
                queue(tail) = g
            end do
        end do
        free = pack([(e, e = 1, ne)], .not. tree .and. .not. dual)
        count = size(free)
        allocate(harmonic(ne, count), value(ne))
        do j = 1, count
            value = 0.0_dp
            value(free(j)) = 1.0_dp
            do q = tail, 2, -1
                f = queue(q)
                e = parent_edge(f)
                s = 0.0_dp
                do k = plane%curl%rowptr(f), plane%curl%rowptr(f + 1) - 1
                    if (plane%curl%colidx(k) == e) then
                        r = k
                    else
                        s = s + plane%curl%val(k)*value(plane%curl%colidx(k))
                    end if
                end do
                value(e) = -s/plane%curl%val(r)
            end do
            harmonic(:, j) = value
        end do
    end subroutine section_harmonic

    subroutine section_incidence(derham, plane)
        !! Section complex of a 3D complex: the first E x N block of grad and
        !! the first F x E block of curl.
        type(spline_derham_t), intent(in) :: derham
        type(spline_derham_t), intent(inout) :: plane

        type(triplet_t) :: tg, tc
        integer :: r, k

        call tg%init()
        call tc%init()
        do r = 1, derham%edge_count
            do k = derham%grad%rowptr(r), derham%grad%rowptr(r + 1) - 1
                if (derham%grad%colidx(k) <= derham%node_count) &
                    call tg%add(r, derham%grad%colidx(k), derham%grad%val(k))
            end do
        end do
        do r = 1, derham%face_count
            do k = derham%curl%rowptr(r), derham%curl%rowptr(r + 1) - 1
                if (derham%curl%colidx(k) <= derham%edge_count) &
                    call tc%add(r, derham%curl%colidx(k), derham%curl%val(k))
            end do
        end do
        plane%toroidal = .false.
        plane%plane_count = 1
        plane%dimension = [derham%node_count, derham%edge_count, derham%face_count, 0]
        call tg%to_csr(derham%edge_count, plane%grad)
        call tc%to_csr(derham%face_count, plane%curl)
    end subroutine section_incidence

    subroutine spline_derham_cell_basis(derham, quad, cell, form, values, dofs, &
            signs)
        !! Section form basis on one cell at the quadrature points:
        !! values(comp, local, q) are parametric proxy components (form 0 and
        !! 2: one, form 1: two), dofs(local) the section degrees of freedom
        !! and signs(local) the orientation factors, so that a field with
        !! coefficients c has proxies sum_local signs c(dofs) values. In 3D
        !! multiply by spline_derham_t_basis along t.
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        integer, intent(in) :: cell, form
        real(dp), allocatable, intent(out) :: values(:, :, :)
        integer, allocatable, intent(out) :: dofs(:), signs(:)

        real(dp), allocatable :: nv1(:), dv1(:), nv2(:), dv2(:)
        integer :: k, e1, e2, f1, f2, p, nq, q, a, b, i, j, n

        p = derham%degree
        k = derham%space%cell_block(cell)
        e1 = derham%space%cell_element(1, cell)
        e2 = derham%space%cell_element(2, cell)
        f1 = derham%space%blocks(k)%first1(e1)
        f2 = derham%space%blocks(k)%first2(e2)
        nq = quad%points_per_cell
        select case (form)
        case (0)
            allocate(values(1, (p + 1)**2, nq), dofs((p + 1)**2), signs((p + 1)**2))
        case (1)
            allocate(values(2, 2*p*(p + 1), nq), dofs(2*p*(p + 1)), &
                signs(2*p*(p + 1)))
        case default
            allocate(values(1, p*p, nq), dofs(p*p), signs(p*p))
        end select
        values = 0.0_dp
        signs = 1
        associate (ids => derham%ids(k))
            do q = 1, nq
                call block_basis(derham, k, quad%xi(:, quad%point_offset(cell) + q - 1), &
                    nv1, dv1, nv2, dv2)
                n = 0
                select case (form)
                case (0)
                    do b = 0, p
                        do a = 0, p
                            n = n + 1
                            i = f1 + a
                            j = f2 + b
                            values(1, n, q) = nv1(i)*nv2(j)
                            dofs(n) = ids%node(i, j)
                        end do
                    end do
                case (1)
                    do b = 0, p
                        do a = 0, p - 1
                            n = n + 1
                            i = f1 + a
                            j = f2 + b
                            values(1, n, q) = dv1(i)*nv2(j)
                            dofs(n) = ids%edge1(i, j)
                            signs(n) = ids%sign1(i, j)
                        end do
                    end do
                    do b = 0, p - 1
                        do a = 0, p
                            n = n + 1
                            i = f1 + a
                            j = f2 + b
                            values(2, n, q) = nv1(i)*dv2(j)
                            dofs(n) = ids%edge2(i, j)
                            signs(n) = ids%sign2(i, j)
                        end do
                    end do
                case default
                    do b = 0, p - 1
                        do a = 0, p - 1
                            n = n + 1
                            i = f1 + a
                            j = f2 + b
                            values(1, n, q) = dv1(i)*dv2(j)
                            dofs(n) = ids%face(i, j)
                        end do
                    end do
                end select
            end do
        end associate
    end subroutine spline_derham_cell_basis

    subroutine spline_derham_boundary_dofs(derham, tag, form, dofs)
        !! Degrees of freedom of V_form with a trace on the sides with the
        !! given layout tag: nodes (form 0), tangential edges and, in 3D,
        !! t-edges at the boundary nodes (form 1), boundary edges times
        !! t-intervals (form 2, 3D only: the normal flux). Sorted, unique.
        type(spline_derham_t), intent(in) :: derham
        integer, intent(in) :: tag, form
        integer, allocatable, intent(out) :: dofs(:)

        logical, allocatable :: bnode(:), bedge(:)
        integer :: k, s, n1, n2, i, j, kt, nn, ne, nf, nt

        nn = derham%node_count
        ne = derham%edge_count
        nf = derham%face_count
        nt = derham%plane_count
        allocate(bnode(nn), bedge(ne))
        bnode = .false.
        bedge = .false.
        do k = 1, derham%block_count
            associate (ids => derham%ids(k))
                n1 = size(ids%node, 1)
                n2 = size(ids%node, 2)
                do s = 1, 4
                    if (derham%space%layout%side_tag(s, k) /= tag) cycle
                    select case (s)
                    case (1, 3)
                        j = merge(1, n2, s == 1)
                        bnode(ids%node(:, j)) = .true.
                        bedge(ids%edge1(:, j)) = .true.
                    case default
                        i = merge(n1, 1, s == 2)
                        bnode(ids%node(i, :)) = .true.
                        bedge(ids%edge2(i, :)) = .true.
                    end select
                end do
            end associate
        end do
        allocate(dofs(0))
        select case (form)
        case (0)
            do kt = 1, nt
                dofs = [dofs, pack([(i + nn*(kt - 1), i = 1, nn)], bnode)]
            end do
        case (1)
            do kt = 1, nt
                dofs = [dofs, pack([(i + ne*(kt - 1), i = 1, ne)], bedge)]
            end do
            if (derham%toroidal) then
                do kt = 1, nt
                    dofs = [dofs, pack([(ne*nt + i + nn*(kt - 1), i = 1, nn)], bnode)]
                end do
            end if
        case (2)
            if (derham%toroidal) then
                do kt = 1, nt
                    dofs = [dofs, pack([(nf*nt + i + ne*(kt - 1), i = 1, ne)], bedge)]
                end do
            end if
        end select
    end subroutine spline_derham_boundary_dofs

    subroutine triplet_init(self)
        class(triplet_t), intent(inout) :: self

        self%count = 0
        allocate(self%row(256), self%col(256), self%val(256))
    end subroutine triplet_init

    subroutine triplet_add(self, r, c, v)
        class(triplet_t), intent(inout) :: self
        integer, intent(in) :: r, c
        real(dp), intent(in) :: v
        integer, allocatable :: gr(:), gc(:)
        real(dp), allocatable :: gv(:)

        if (self%count == size(self%row)) then
            allocate(gr(2*self%count), gc(2*self%count), gv(2*self%count))
            gr(:self%count) = self%row
            gc(:self%count) = self%col
            gv(:self%count) = self%val
            call move_alloc(gr, self%row)
            call move_alloc(gc, self%col)
            call move_alloc(gv, self%val)
        end if
        self%count = self%count + 1
        self%row(self%count) = r
        self%col(self%count) = c
        self%val(self%count) = v
    end subroutine triplet_add

    subroutine triplet_to_csr(self, nrow, matrix)
        !! Rows sorted by column, duplicates summed, zeros dropped.
        class(triplet_t), intent(in) :: self
        integer, intent(in) :: nrow
        type(csr_matrix_t), intent(out) :: matrix

        integer, allocatable :: ptr(:), next(:), col(:)
        real(dp), allocatable :: val(:)
        integer :: q, r, i, j, n, x
        real(dp) :: y

        allocate(ptr(nrow + 1), next(nrow + 1), col(self%count), val(self%count))
        ptr = 0
        do q = 1, self%count
            ptr(self%row(q) + 1) = ptr(self%row(q) + 1) + 1
        end do
        ptr(1) = 1
        do r = 2, nrow + 1
            ptr(r) = ptr(r) + ptr(r - 1)
        end do
        next = ptr
        do q = 1, self%count
            col(next(self%row(q))) = self%col(q)
            val(next(self%row(q))) = self%val(q)
            next(self%row(q)) = next(self%row(q)) + 1
        end do
        matrix%n = nrow
        allocate(matrix%rowptr(nrow + 1), matrix%colidx(self%count), &
            matrix%val(self%count))
        n = 0
        matrix%rowptr(1) = 1
        do r = 1, nrow
            do i = ptr(r) + 1, ptr(r + 1) - 1
                x = col(i)
                y = val(i)
                j = i - 1
                do while (j >= ptr(r))
                    if (col(j) <= x) exit
                    col(j + 1) = col(j)
                    val(j + 1) = val(j)
                    j = j - 1
                end do
                col(j + 1) = x
                val(j + 1) = y
            end do
            do i = ptr(r), ptr(r + 1) - 1
                if (n >= matrix%rowptr(r)) then
                    if (matrix%colidx(n) == col(i)) then
                        matrix%val(n) = matrix%val(n) + val(i)
                        cycle
                    end if
                end if
                n = n + 1
                matrix%colidx(n) = col(i)
                matrix%val(n) = val(i)
            end do
            matrix%rowptr(r + 1) = n + 1
        end do
        matrix%colidx = matrix%colidx(:n)
        matrix%val = matrix%val(:n)
    end subroutine triplet_to_csr

    subroutine spline_derham_t_basis(nt, t, nodal, edge)
        !! Trigonometric cardinal (nodal) and histopolation (edge) bases of
        !! nt = 2M + 1 points at t, with their t-derivatives in column 2.
        integer, intent(in) :: nt
        real(dp), intent(in) :: t
        real(dp), intent(out) :: nodal(nt, 2), edge(nt, 2)

        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: h, x, s
        integer :: k, m

        h = 1.0_dp/nt
        do k = 1, nt
            x = t - (k - 1)*h
            nodal(k, :) = [1.0_dp, 0.0_dp]
            edge(k, :) = [1.0_dp, 0.0_dp]
            do m = 1, (nt - 1)/2
                s = sin(pi*m*h)/(pi*m*h)
                nodal(k, 1) = nodal(k, 1) + 2.0_dp*cos(2.0_dp*pi*m*x)
                nodal(k, 2) = nodal(k, 2) - 4.0_dp*pi*m*sin(2.0_dp*pi*m*x)
                edge(k, 1) = edge(k, 1) + 2.0_dp*cos(2.0_dp*pi*m*(x - 0.5_dp*h))/s
                edge(k, 2) = edge(k, 2) - 4.0_dp*pi*m*sin(2.0_dp*pi*m*(x - 0.5_dp*h))/s
            end do
            nodal(k, :) = nodal(k, :)/nt
        end do
    end subroutine spline_derham_t_basis

    subroutine block_basis(derham, k, xi, nv1, dv1, nv2, dv2)
        !! N and D values of block k at xi (all functions).
        type(spline_derham_t), intent(in) :: derham
        integer, intent(in) :: k
        real(dp), intent(in) :: xi(2)
        real(dp), allocatable, intent(out) :: nv1(:), dv1(:), nv2(:), dv2(:)

        call direction_basis(derham%space%blocks(k)%knots1, derham%degree, xi(1), &
            nv1, dv1)
        call direction_basis(derham%space%blocks(k)%knots2, derham%degree, xi(2), &
            nv2, dv2)
    end subroutine block_basis

    subroutine direction_basis(knots, p, x, nv, dv)
        real(dp), intent(in) :: knots(:), x
        integer, intent(in) :: p
        real(dp), allocatable, intent(out) :: nv(:), dv(:)

        real(dp), allocatable :: derivative(:), lower(:)
        integer :: n, i, status

        n = size(knots) - p - 1
        call evaluate_bspline_basis(knots, p, x, nv, derivative, status)
        call evaluate_bspline_basis(knots(2:size(knots) - 1), p - 1, x, lower, &
            derivative, status)
        allocate(dv(n - 1))
        do i = 1, n - 1
            dv(i) = real(p, dp)/(knots(i + p + 1) - knots(i + 1))*lower(i)
        end do
    end subroutine direction_basis

    subroutine spline_derham_evaluate(derham, form, coef, block, xi, t, values)
        !! Parametric proxy components of a discrete form at (block, xi, t):
        !! form 0: 1, form 1: 2 (section) or 3, form 2: 1 (section) or 3,
        !! form 3: 1. t is ignored for a section complex.
        type(spline_derham_t), intent(in) :: derham
        integer, intent(in) :: form, block
        real(dp), intent(in) :: coef(:), xi(2), t
        real(dp), intent(out) :: values(:)

        real(dp), allocatable :: nv1(:), dv1(:), nv2(:), dv2(:)
        real(dp) :: lt(derham%plane_count, 2), ht(derham%plane_count, 2)
        integer :: nn, ne, nf, nt, k, i, j, off

        call block_basis(derham, block, xi, nv1, dv1, nv2, dv2)
        nn = derham%node_count
        ne = derham%edge_count
        nf = derham%face_count
        nt = derham%plane_count
        lt = 0.0_dp
        ht = 0.0_dp
        if (derham%toroidal) then
            call spline_derham_t_basis(nt, t, lt, ht)
        else
            lt(1, 1) = 1.0_dp
        end if
        values = 0.0_dp
        associate (ids => derham%ids(block))
            do k = 1, nt
                select case (form)
                case (0)
                    do j = 1, size(nv2)
                        do i = 1, size(nv1)
                            values(1) = values(1) + coef(ids%node(i, j) + nn*(k - 1))* &
                                nv1(i)*nv2(j)*lt(k, 1)
                        end do
                    end do
                case (1)
                    do j = 1, size(nv2)
                        do i = 1, size(dv1)
                            values(1) = values(1) + ids%sign1(i, j)* &
                                coef(ids%edge1(i, j) + ne*(k - 1))*dv1(i)*nv2(j)*lt(k, 1)
                        end do
                    end do
                    do j = 1, size(dv2)
                        do i = 1, size(nv1)
                            values(2) = values(2) + ids%sign2(i, j)* &
                                coef(ids%edge2(i, j) + ne*(k - 1))*nv1(i)*dv2(j)*lt(k, 1)
                        end do
                    end do
                    if (derham%toroidal) then
                        off = ne*nt
                        do j = 1, size(nv2)
                            do i = 1, size(nv1)
                                values(3) = values(3) + coef(off + ids%node(i, j) + &
                                    nn*(k - 1))*nv1(i)*nv2(j)*ht(k, 1)
                            end do
                        end do
                    end if
                case (2)
                    if (derham%toroidal) then
                        do j = 1, size(dv2)
                            do i = 1, size(dv1)
                                values(3) = values(3) + coef(ids%face(i, j) + &
                                    nf*(k - 1))*dv1(i)*dv2(j)*lt(k, 1)
                            end do
                        end do
                        off = nf*nt
                        do j = 1, size(nv2)
                            do i = 1, size(dv1)
                                values(2) = values(2) + ids%sign1(i, j)*coef(off + &
                                    ids%edge1(i, j) + ne*(k - 1))*dv1(i)*nv2(j)*ht(k, 1)
                            end do
                        end do
                        do j = 1, size(dv2)
                            do i = 1, size(nv1)
                                values(1) = values(1) - ids%sign2(i, j)*coef(off + &
                                    ids%edge2(i, j) + ne*(k - 1))*nv1(i)*dv2(j)*ht(k, 1)
                            end do
                        end do
                    else
                        do j = 1, size(dv2)
                            do i = 1, size(dv1)
                                values(1) = values(1) + coef(ids%face(i, j))*dv1(i)*dv2(j)
                            end do
                        end do
                    end if
                case default
                    do j = 1, size(dv2)
                        do i = 1, size(dv1)
                            values(1) = values(1) + coef(ids%face(i, j) + nf*(k - 1))* &
                                dv1(i)*dv2(j)*ht(k, 1)
                        end do
                    end do
                end select
            end do
        end associate
    end subroutine spline_derham_evaluate

    subroutine spline_derham_project(derham, form, field, coef, status, points)
        !! Commuting projection of a field given by its parametric proxy
        !! components (derham_field) onto V_form, from point values, edge
        !! and cell integrals at the Greville points (Gauss rules with
        !! `points` nodes per knot span, default p + 3).
        type(spline_derham_t), intent(in) :: derham
        integer, intent(in) :: form
        procedure(derham_field) :: field
        real(dp), intent(out) :: coef(:)
        integer, intent(out) :: status
        integer, intent(in), optional :: points

        real(dp), allocatable :: x(:), w(:)
        integer :: nq, k, kt, nt

        status = 1
        if (form < 0 .or. form > 3) return
        if (.not. derham%toroidal .and. form == 3) return
        if (size(coef) /= derham%dimension(form)) return
        nq = derham%degree + 3
        if (present(points)) nq = points
        call gauss_rule(nq, x, w, status)
        if (status /= 0) return
        x = 0.5_dp*(x + 1.0_dp)
        w = 0.5_dp*w
        nt = derham%plane_count
        coef = 0.0_dp
        do k = 1, derham%block_count
            do kt = 1, nt
                call project_block(derham, k, kt, form, field, x, w, coef, status)
                if (status /= 0) return
            end do
        end do
        status = 0
    end subroutine spline_derham_project

    subroutine project_block(derham, k, kt, form, field, x, w, coef, status)
        !! Block k, t node or interval kt: geometric degrees of freedom of
        !! every component, then the tensor solves with the collocation and
        !! histopolation matrices of the block.
        type(spline_derham_t), intent(in) :: derham
        integer, intent(in) :: k, kt, form
        procedure(derham_field) :: field
        real(dp), intent(in) :: x(:), w(:)
        real(dp), intent(inout) :: coef(:)
        integer, intent(out) :: status

        real(dp), allocatable :: g1(:), g2(:), mn1(:, :), md1(:, :), mn2(:, :), md2(:, :)
        real(dp), allocatable :: sigma(:, :)
        real(dp) :: tk, h
        integer :: n1, n2, nn, ne, nf, nt, ncomp

        associate (blk => derham%space%blocks(k), ids => derham%ids(k))
            n1 = blk%n(1)
            n2 = blk%n(2)
            allocate(g1(n1), g2(n2))
            call spline_greville(blk%knots1, derham%degree, g1)
            call spline_greville(blk%knots2, derham%degree, g2)
            call direction_matrices(blk%knots1, derham%degree, g1, x, w, mn1, md1)
            call direction_matrices(blk%knots2, derham%degree, g2, x, w, mn2, md2)
            nn = derham%node_count
            ne = derham%edge_count
            nf = derham%face_count
            nt = derham%plane_count
            h = 1.0_dp/nt
            tk = (kt - 1)*h
            ncomp = 3
            status = 0
            select case (form)
            case (0)
                call dofs(1, .false., .false., .false., sigma)
                call solve2(mn1, mn2, sigma)
                call put(ids%node, sigma, (kt - 1)*nn)
            case (1)
                call dofs(1, .true., .false., .false., sigma)
                call solve2(md1, mn2, sigma)
                call put_signed(ids%edge1, ids%sign1, sigma, (kt - 1)*ne, 1.0_dp)
                call dofs(2, .false., .true., .false., sigma)
                call solve2(mn1, md2, sigma)
                call put_signed(ids%edge2, ids%sign2, sigma, (kt - 1)*ne, 1.0_dp)
                if (derham%toroidal) then
                    call dofs(3, .false., .false., .true., sigma)
                    call solve2(mn1, mn2, sigma)
                    call put(ids%node, sigma, ne*nt + (kt - 1)*nn)
                end if
            case (2)
                if (derham%toroidal) then
                    call dofs(3, .true., .true., .false., sigma)
                    call solve2(md1, md2, sigma)
                    call put(ids%face, sigma, (kt - 1)*nf)
                    call dofs(2, .true., .false., .true., sigma)
                    call solve2(md1, mn2, sigma)
                    call put_signed(ids%edge1, ids%sign1, sigma, nf*nt + (kt - 1)*ne, &
                        1.0_dp)
                    call dofs(1, .false., .true., .true., sigma)
                    call solve2(mn1, md2, sigma)
                    call put_signed(ids%edge2, ids%sign2, sigma, nf*nt + (kt - 1)*ne, &
                        -1.0_dp)
                else
                    call dofs(1, .true., .true., .false., sigma)
                    call solve2(md1, md2, sigma)
                    call put(ids%face, sigma, 0)
                end if
            case default
                call dofs(1, .true., .true., .true., sigma)
                call solve2(md1, md2, sigma)
                call put(ids%face, sigma, (kt - 1)*nf)
            end select
        end associate

    contains

        subroutine dofs(comp, int1, int2, int3, sigma)
            !! Values (or integrals over the Greville intervals in the flagged
            !! directions) of component comp.
            integer, intent(in) :: comp
            logical, intent(in) :: int1, int2, int3
            real(dp), allocatable, intent(out) :: sigma(:, :)

            real(dp), allocatable :: p1(:, :), p2(:, :), p3(:, :)
            real(dp) :: v(ncomp)
            integer :: i, j, a, b, c, m1, m2

            m1 = merge(n1 - 1, n1, int1)
            m2 = merge(n2 - 1, n2, int2)
            call nodes(g1, int1, p1)
            call nodes(g2, int2, p2)
            if (int3) then
                allocate(p3(2, size(x)))
                p3(1, :) = tk + h*x
                p3(2, :) = h*w
            else
                allocate(p3(2, 1))
                p3(:, 1) = [tk, 1.0_dp]
            end if
            allocate(sigma(m1, m2))
            sigma = 0.0_dp
            do j = 1, m2
                do i = 1, m1
                    do c = 1, size(p3, 2)
                        do b = 1, size(p2, 1)/2
                            do a = 1, size(p1, 1)/2
                                call field(k, [p1(2*a - 1, i), p2(2*b - 1, j)], p3(1, c), v)
                                sigma(i, j) = sigma(i, j) + p1(2*a, i)*p2(2*b, j)* &
                                    p3(2, c)*v(comp)
                            end do
                        end do
                    end do
                end do
            end do
        end subroutine dofs

        subroutine nodes(g, integrate, p)
            !! p(2a - 1, i), p(2a, i): points and weights of DOF i (one point
            !! of weight 1, or Gauss rules on the knot spans of [g_i, g_i+1]).
            real(dp), intent(in) :: g(:)
            logical, intent(in) :: integrate
            real(dp), allocatable, intent(out) :: p(:, :)

            real(dp), allocatable :: cuts(:), knots(:)
            integer :: i, s, q, nseg, maxseg

            if (.not. integrate) then
                allocate(p(2, size(g)))
                p(1, :) = g
                p(2, :) = 1.0_dp
                return
            end if
            knots = blk_breaks(g)
            maxseg = size(knots) + 1
            allocate(p(2*size(x)*maxseg, size(g) - 1))
            p = 0.0_dp
            do i = 1, size(g) - 1
                cuts = [g(i), pack(knots, knots > g(i) .and. knots < g(i + 1)), g(i + 1)]
                nseg = size(cuts) - 1
                do s = 1, nseg
                    do q = 1, size(x)
                        p(2*(q + size(x)*(s - 1)) - 1, i) = cuts(s) + &
                            (cuts(s + 1) - cuts(s))*x(q)
                        p(2*(q + size(x)*(s - 1)), i) = (cuts(s + 1) - cuts(s))*w(q)
                    end do
                end do
                ! unused slots keep weight 0 at point 0
                do q = 2*size(x)*nseg + 1, size(p, 1), 2
                    p(q, i) = g(i)
                end do
            end do
        end subroutine nodes

        function blk_breaks(g) result(b)
            real(dp), intent(in) :: g(:)
            real(dp), allocatable :: b(:)

            if (size(g) == n1) then
                b = derham%space%blocks(k)%breaks1
            else
                b = derham%space%blocks(k)%breaks2
            end if
        end function blk_breaks

        subroutine put(id, sigma, offset)
            integer, intent(in) :: id(:, :), offset
            real(dp), intent(in) :: sigma(:, :)

            coef(offset + reshape(id, [size(id)])) = reshape(sigma, [size(sigma)])
        end subroutine put

        subroutine put_signed(id, sgn, sigma, offset, factor)
            integer, intent(in) :: id(:, :), sgn(:, :), offset
            real(dp), intent(in) :: sigma(:, :), factor

            coef(offset + reshape(id, [size(id)])) = factor* &
                reshape(sigma*real(sgn, dp), [size(sigma)])
        end subroutine put_signed

    end subroutine project_block

    subroutine solve2(a, b, sigma)
        !! sigma <- a^(-1) sigma b^(-T).
        real(dp), intent(in) :: a(:, :), b(:, :)
        real(dp), intent(inout) :: sigma(:, :)

        real(dp), allocatable :: m(:, :), r(:, :)
        integer, allocatable :: pivot(:)
        integer :: info

        m = a
        allocate(pivot(max(size(a, 1), size(b, 1))))
        call dgesv(size(a, 1), size(sigma, 2), m, size(a, 1), pivot, sigma, &
            size(a, 1), info)
        r = transpose(sigma)
        m = b
        call dgesv(size(b, 1), size(r, 2), m, size(b, 1), pivot, r, size(b, 1), info)
        sigma = transpose(r)
    end subroutine solve2

    subroutine direction_matrices(knots, p, g, x, w, mn, md)
        !! mn(k, i) = N_i(g_k); md(k, i) = int_{g_k}^{g_k+1} D_i.
        real(dp), intent(in) :: knots(:), g(:), x(:), w(:)
        integer, intent(in) :: p
        real(dp), allocatable, intent(out) :: mn(:, :), md(:, :)

        real(dp), allocatable :: nv(:), dv(:), cuts(:), breaks(:)
        integer :: n, k, s, q

        n = size(g)
        allocate(mn(n, n), md(n - 1, n - 1))
        breaks = knots(p + 1:size(knots) - p)
        md = 0.0_dp
        do k = 1, n
            call direction_basis(knots, p, g(k), nv, dv)
            mn(k, :) = nv
        end do
        do k = 1, n - 1
            cuts = [g(k), pack(breaks, breaks > g(k) .and. breaks < g(k + 1)), g(k + 1)]
            do s = 1, size(cuts) - 1
                do q = 1, size(x)
                    call direction_basis(knots, p, cuts(s) + (cuts(s + 1) - cuts(s))*x(q), &
                        nv, dv)
                    md(k, :) = md(k, :) + (cuts(s + 1) - cuts(s))*w(q)*dv
                end do
            end do
        end do
    end subroutine direction_matrices

end module fortfem_spline_derham
