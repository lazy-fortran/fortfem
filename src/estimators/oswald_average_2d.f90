module fortfem_oswald_average_2d
    !! Oswald (nodal averaging) interpolation of a broken Lagrange P_m field on
    !! a conforming triangle mesh into the continuous Lagrange P_m space that
    !! vanishes on the domain boundary: the value at every shared Lagrange node
    !! is the arithmetic mean of the values of all cells containing it, and the
    !! nodes on boundary vertices and boundary edges are set to zero.
    !!
    !! Local nodes are ordered like fortfem_triangle_lagrange_arbitrary_order,
    !! so local vertex s is cell vertex s. A node is identified through its
    !! barycentric multi-index (integer, sum m): a vertex node by the global
    !! vertex, an edge node by the global edge and its distance index from the
    !! endpoint with the smaller global vertex number, an interior node by its
    !! cell. The averaging is a linear, local O(N) operator.
    use fortfem_kinds, only: dp
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        initialize_triangle_lagrange_basis, triangle_lagrange_basis_t, &
        triangle_lagrange_nodes
    implicit none
    private

    public :: oswald_average_homogeneous, lagrange_continuous_numbering

contains

    subroutine lagrange_continuous_numbering(triangles, degree, bary, node_id, &
            boundary_node, status)
        !! Global numbering of the nodes of continuous Lagrange P_m (m =
        !! degree >= 1): node_id(l, t) is the global node of local node l of
        !! cell t, bary(:, l) the integer barycentric multi-index (sum m) of
        !! local node l and boundary_node(g) is true for the nodes on the
        !! domain boundary (boundary vertices and boundary edges).
        integer, intent(in) :: triangles(:, :), degree
        integer, allocatable, intent(out) :: bary(:, :), node_id(:, :)
        logical, allocatable, intent(out) :: boundary_node(:)
        integer, intent(out) :: status
        type(triangle_lagrange_basis_t) :: basis
        real(dp), allocatable :: nodes(:, :)
        integer, allocatable :: cell_edge(:, :), edge_cells(:)
        integer :: nt, nv, ne, nl, ni, nnode, t, l, e, s, positive, lo
        integer :: ids(2), c(3), idx

        status = 1
        nt = size(triangles, 2)
        if (degree < 1 .or. nt < 1 .or. size(triangles, 1) /= 3) return
        call initialize_triangle_lagrange_basis(degree, basis, status)
        if (status /= 0) return
        status = 1
        nl = (degree + 1)*(degree + 2)/2
        call triangle_lagrange_nodes(basis, nodes)
        allocate (bary(3, nl))
        do l = 1, nl
            bary(:, l) = nint(real(degree, dp)*[1.0_dp - nodes(1, l) &
                - nodes(2, l), nodes(1, l), nodes(2, l)])
            if (sum(bary(:, l)) /= degree) return
        end do
        nv = maxval(triangles)
        call number_edges(nv, triangles, cell_edge, edge_cells, ne)
        ni = (degree - 1)*(degree - 2)/2
        nnode = nv + ne*(degree - 1) + nt*ni
        allocate (node_id(nl, nt), boundary_node(nnode))
        boundary_node = .false.
        do t = 1, nt
            idx = 0
            do l = 1, nl
                c = bary(:, l)
                positive = count_positive(c)
                if (positive == 1) then
                    s = maxloc(c, dim=1)
                    node_id(l, t) = triangles(s, t)
                else if (positive == 2) then
                    e = edge_of_node(c)
                    ids = [triangles(e, t), triangles(mod(e, 3) + 1, t)]
                    lo = merge(e, mod(e, 3) + 1, ids(1) < ids(2))
                    node_id(l, t) = nv + (cell_edge(e, t) - 1)*(degree - 1) &
                        + c(lo)
                    if (edge_cells(cell_edge(e, t)) == 1) &
                        boundary_node(node_id(l, t)) = .true.
                else
                    idx = idx + 1
                    node_id(l, t) = nv + ne*(degree - 1) + (t - 1)*ni + idx
                end if
            end do
            do e = 1, 3
                if (edge_cells(cell_edge(e, t)) /= 1) cycle
                boundary_node(triangles(e, t)) = .true.
                boundary_node(triangles(mod(e, 3) + 1, t)) = .true.
            end do
        end do
        status = 0
    end subroutine lagrange_continuous_numbering

    subroutine oswald_average_homogeneous(triangles, degree, local_values, &
            averaged, status)
        !! triangles(3, nt) counter-clockwise or not (only connectivity is
        !! used); local_values(nl, nt) nodal values of the broken field.
        integer, intent(in) :: triangles(:, :), degree
        real(dp), intent(in) :: local_values(:, :)
        real(dp), intent(out) :: averaged(:, :)
        integer, intent(out) :: status
        integer, allocatable :: bary(:, :), node_id(:, :), multiplicity(:)
        logical, allocatable :: boundary(:)
        real(dp), allocatable :: total(:)
        integer :: nt, nl, t, l, g

        averaged = 0.0_dp
        call lagrange_continuous_numbering(triangles, degree, bary, node_id, &
            boundary, status)
        if (status /= 0) return
        status = 1
        nt = size(triangles, 2)
        nl = size(bary, 2)
        if (size(local_values, 1) /= nl .or. size(local_values, 2) /= nt) return
        if (size(averaged, 1) /= nl .or. size(averaged, 2) /= nt) return
        allocate (total(size(boundary)), multiplicity(size(boundary)))
        total = 0.0_dp
        multiplicity = 0
        do t = 1, nt
            do l = 1, nl
                g = node_id(l, t)
                total(g) = total(g) + local_values(l, t)
                multiplicity(g) = multiplicity(g) + 1
            end do
        end do
        do t = 1, nt
            do l = 1, nl
                g = node_id(l, t)
                if (boundary(g)) cycle
                averaged(l, t) = total(g)/real(multiplicity(g), dp)
            end do
        end do
        status = 0
    end subroutine oswald_average_homogeneous

    pure integer function count_positive(c) result(n)
        integer, intent(in) :: c(3)

        n = count(c > 0)
    end function count_positive

    pure integer function edge_of_node(c) result(e)
        !! Local edge (e from vertex e to mod(e, 3) + 1) holding the node with
        !! one vanishing barycentric index.
        integer, intent(in) :: c(3)

        e = mod(minloc(c, dim=1), 3) + 1
    end function edge_of_node

    subroutine number_edges(nv, triangles, cell_edge, edge_cells, ne)
        !! Global edge numbers of every cell edge and the number of cells
        !! (1 = boundary edge, 2 = interior edge) of every edge.
        integer, intent(in) :: nv, triangles(:, :)
        integer, allocatable, intent(out) :: cell_edge(:, :), edge_cells(:)
        integer, intent(out) :: ne
        integer, allocatable :: start(:), fill(:), key_hi(:), key_edge(:)
        integer :: nt, t, e, a, b, lo, hi, k, slot

        nt = size(triangles, 2)
        allocate (start(nv + 1), fill(nv), cell_edge(3, nt))
        allocate (key_hi(3*nt), key_edge(3*nt))
        start = 0
        do t = 1, nt
            do e = 1, 3
                a = triangles(e, t)
                b = triangles(mod(e, 3) + 1, t)
                lo = min(a, b)
                start(lo + 1) = start(lo + 1) + 1
            end do
        end do
        start(1) = 1
        do a = 1, nv
            start(a + 1) = start(a + 1) + start(a)
        end do
        fill = 0
        ne = 0
        cell_edge = 0
        do t = 1, nt
            do e = 1, 3
                a = triangles(e, t)
                b = triangles(mod(e, 3) + 1, t)
                lo = min(a, b)
                hi = max(a, b)
                do k = start(lo), start(lo) + fill(lo) - 1
                    if (key_hi(k) == hi) then
                        cell_edge(e, t) = key_edge(k)
                        exit
                    end if
                end do
                if (cell_edge(e, t) == 0) then
                    ne = ne + 1
                    cell_edge(e, t) = ne
                    slot = start(lo) + fill(lo)
                    fill(lo) = fill(lo) + 1
                    key_hi(slot) = hi
                    key_edge(slot) = ne
                end if
            end do
        end do
        allocate (edge_cells(ne))
        edge_cells = 0
        do t = 1, nt
            do e = 1, 3
                edge_cells(cell_edge(e, t)) = edge_cells(cell_edge(e, t)) + 1
            end do
        end do
    end subroutine number_edges

end module fortfem_oswald_average_2d
