module fortfem_tensor_boundary_trace
    !! Boundary trace of a 3D tensor-product hexahedral mesh for BEM coupling.
    !! Boundary faces are the cell faces whose four corner vertices occur once.
    !! Each face becomes one quadrilateral BEM panel on the reference square
    !! [0,1]^2. The panel parametrisation is oriented so that a x b points out
    !! of the meshed region, and its geometry is the cell's own isoparametric
    !! map restricted to the face. Trace unknowns are the volume GLL face nodes.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortfem_tensor_mesh, only: tensor_mesh_t
    use fortfem_gll_line, only: gll_rule, lagrange_jet
    use fortfem_parametric_laplace_bem_3d, only: bem_surface_t
    implicit none
    private
    public :: tensor_boundary_surface_t, tensor_boundary_trace_t
    public :: build_tensor_boundary_trace

    type, extends(bem_surface_t) :: tensor_boundary_surface_t
        !! Panel geometry nodes on the face GLL grid of the geometry degree.
        real(dp), allocatable :: nodes(:), geometry(:, :, :)
    contains
        procedure :: evaluate => evaluate_face
    end type

    type :: tensor_boundary_trace_t
        integer :: degree = 0, npanel = 0, ntrace = 0, genus = -1
        !! panel_dof(i + (p+1)*(j-1), k): trace index of face node (i,j).
        integer, allocatable :: panel_cell(:), panel_dof(:, :), corners(:, :)
        !! trace_dof(t): volume dof of trace index t.
        integer, allocatable :: trace_dof(:)
        type(tensor_boundary_surface_t) :: surface
    end type
contains
    pure subroutine evaluate_face(self, panel, uv, x, a, b, status)
        class(tensor_boundary_surface_t), intent(in) :: self
        integer, intent(in) :: panel
        real(dp), intent(in) :: uv(2)
        real(dp), intent(out) :: x(3), a(3), b(3)
        integer, intent(out) :: status
        real(dp) :: bu(size(self%nodes)), du(size(self%nodes))
        real(dp) :: bv(size(self%nodes)), dv(size(self%nodes))
        integer :: i, j, n
        x = 0
        a = 0
        b = 0
        status = 1
        if (panel < 1 .or. panel > size(self%geometry, 3)) return
        n = size(self%nodes)
        call lagrange_jet(self%nodes, uv(1), bu, du)
        call lagrange_jet(self%nodes, uv(2), bv, dv)
        do j = 1, n
            do i = 1, n
                x = x + bu(i)*bv(j)*self%geometry(:, i + n*(j - 1), panel)
                a = a + du(i)*bv(j)*self%geometry(:, i + n*(j - 1), panel)
                b = b + bu(i)*dv(j)*self%geometry(:, i + n*(j - 1), panel)
            end do
        end do
        status = 0
    end subroutine

    pure function local_index(pos, n) result(index)
        integer, intent(in) :: pos(3), n
        integer :: index
        index = 1 + pos(1) + n*pos(2) + n*n*pos(3)
    end function

    pure subroutine face_frame(dir, side, du, dv)
        !! Tangent directions with e_du x e_dv along the outward reference normal.
        integer, intent(in) :: dir, side
        integer, intent(out) :: du, dv
        integer :: swap, sign_ref
        select case (dir)
        case (1)
            du = 2; dv = 3; sign_ref = 1
        case (2)
            du = 1; dv = 3; sign_ref = -1
        case default
            du = 1; dv = 2; sign_ref = 1
        end select
        if (sign_ref*(2*side - 1) < 0) then
            swap = du; du = dv; dv = swap
        end if
    end subroutine

    subroutine build_tensor_boundary_trace(mesh, trace, status)
        !! Status 1 rejects the mesh, 2 a non-manifold or inconsistently
        !! oriented boundary (a face used by more than two cells, or a
        !! boundary edge not shared by exactly two oppositely oriented panels).
        type(tensor_mesh_t), intent(in) :: mesh
        type(tensor_boundary_trace_t), intent(out) :: trace
        integer, intent(out) :: status
        integer, allocatable :: keys(:, :), table(:), count(:), owner(:, :)
        integer, allocatable :: map(:), edge_keys(:, :), edge_table(:), edge_sum(:)
        integer, allocatable :: edge_count(:)
        real(dp), allocatable :: z(:), w(:)
        integer :: e, dir, side, c, i, j, k, l, n, p, q, nkeys, hsize, slot, swap
        integer :: bits(3), pos(3), key(4), du, dv, panel, nv, nedge, vlo, vhi
        integer :: corner(4), corners(4), ekey(2)
        status = 1
        if (mesh%dim /= 3 .or. mesh%nelem < 1 .or. mesh%degree < 1) return
        if (.not. allocated(mesh%vertices) .or. .not. allocated(mesh%dof)) return
        p = mesh%degree
        q = mesh%geometry_degree
        hsize = 12*mesh%nelem + 1
        allocate (keys(4, 6*mesh%nelem), table(hsize), count(6*mesh%nelem), &
                  owner(2, 6*mesh%nelem))
        table = 0
        count = 0
        nkeys = 0
        do e = 1, mesh%nelem
            do dir = 1, 3
                do side = 0, 1
                    l = 0
                    do c = 1, 8
                        bits = [modulo(c - 1, 2), modulo((c - 1)/2, 2), (c - 1)/4]
                        if (bits(dir) /= side) cycle
                        l = l + 1
                        key(l) = mesh%vertices(c, e)
                    end do
                    call sort4(key)
                    slot = hash_slot(key, hsize)
                    do
                        if (table(slot) == 0) then
                            nkeys = nkeys + 1
                            table(slot) = nkeys
                            keys(:, nkeys) = key
                            exit
                        end if
                        if (all(keys(:, table(slot)) == key)) exit
                        slot = modulo(slot, hsize) + 1
                    end do
                    i = table(slot)
                    count(i) = count(i) + 1
                    if (count(i) == 1) owner(:, i) = [e, 2*dir - 1 + side]
                end do
            end do
        end do
        status = 2
        if (any(count(1:nkeys) > 2)) return
        trace%degree = p
        trace%npanel = count_one(count(1:nkeys))
        if (trace%npanel < 4) return
        n = p + 1
        allocate (trace%panel_cell(trace%npanel), trace%panel_dof(n*n, trace%npanel), &
                  trace%corners(4, trace%npanel), map(mesh%ndof), &
                  trace%surface%geometry(3, (q + 1)**2, trace%npanel))
        call gll_rule(q, z, w, status)
        if (status /= 0) return
        trace%surface%nodes = (z + 1)/2
        map = 0
        trace%ntrace = 0
        panel = 0
        do i = 1, nkeys
            if (count(i) /= 1) cycle
            panel = panel + 1
            e = owner(1, i)
            dir = (owner(2, i) + 1)/2
            side = owner(2, i) - 2*dir + 1
            call face_frame(dir, side, du, dv)
            trace%panel_cell(panel) = e
            do k = 0, p
                do j = 0, p
                    pos(dir) = side*p
                    pos(du) = j
                    pos(dv) = k
                    l = mesh%dof(local_index(pos, n), e)
                    if (map(l) == 0) then
                        trace%ntrace = trace%ntrace + 1
                        map(l) = trace%ntrace
                    end if
                    trace%panel_dof(1 + j + n*k, panel) = map(l)
                end do
            end do
            do k = 0, q
                do j = 0, q
                    pos(dir) = side*q
                    pos(du) = j
                    pos(dv) = k
                    trace%surface%geometry(:, 1 + j + (q + 1)*k, panel) = &
                        mesh%geometry(:, local_index(pos, q + 1), e)
                end do
            end do
            ! Corner vertices in counter-clockwise (u,v) order.
            do c = 1, 4
                pos(dir) = side
                pos(du) = merge(1, 0, c == 2 .or. c == 3)
                pos(dv) = merge(1, 0, c >= 3)
                trace%corners(c, panel) = mesh%vertices(local_index(pos, 2), e)
            end do
        end do
        allocate (trace%trace_dof(trace%ntrace))
        do l = 1, mesh%ndof
            if (map(l) > 0) trace%trace_dof(map(l)) = l
        end do
        ! Closed, consistently oriented boundary: each edge occurs once in each
        ! direction. Euler characteristic V - E + F = 2 - 2 genus.
        hsize = 8*trace%npanel + 1
        allocate (edge_keys(2, 4*trace%npanel), edge_table(hsize), &
                  edge_sum(4*trace%npanel), edge_count(4*trace%npanel))
        edge_table = 0
        edge_sum = 0
        edge_count = 0
        nedge = 0
        do panel = 1, trace%npanel
            corners = trace%corners(:, panel)
            do c = 1, 4
                vlo = corners(c)
                vhi = corners(modulo(c, 4) + 1)
                ekey = [min(vlo, vhi), max(vlo, vhi)]
                corner = [ekey, 0, 0]
                slot = hash_slot(corner, hsize)
                do
                    if (edge_table(slot) == 0) then
                        nedge = nedge + 1
                        edge_table(slot) = nedge
                        edge_keys(:, nedge) = ekey
                        exit
                    end if
                    if (all(edge_keys(:, edge_table(slot)) == ekey)) exit
                    slot = modulo(slot, hsize) + 1
                end do
                j = edge_table(slot)
                edge_count(j) = edge_count(j) + 1
                edge_sum(j) = edge_sum(j) + merge(1, -1, vlo < vhi)
            end do
        end do
        if (any(edge_count(1:nedge) /= 2)) return
        if (any(edge_sum(1:nedge) /= 0)) return
        nv = 0
        deallocate (map)
        allocate (map(maxval(trace%corners)))
        map = 0
        do panel = 1, trace%npanel
            do c = 1, 4
                swap = trace%corners(c, panel)
                if (map(swap) == 0) then
                    nv = nv + 1
                    map(swap) = 1
                end if
            end do
        end do
        if (mod(2 - (nv - nedge + trace%npanel), 2) /= 0) return
        trace%genus = (2 - (nv - nedge + trace%npanel))/2
        if (trace%genus < 0) return
        status = 0
    contains
        pure integer function count_one(values)
            integer, intent(in) :: values(:)
            integer :: m
            count_one = 0
            do m = 1, size(values)
                if (values(m) == 1) count_one = count_one + 1
            end do
        end function
    end subroutine

    pure subroutine sort4(key)
        integer, intent(inout) :: key(4)
        integer :: i, j, swap
        do i = 2, 4
            j = i
            do while (j > 1)
                if (key(j) >= key(j - 1)) exit
                swap = key(j); key(j) = key(j - 1); key(j - 1) = swap
                j = j - 1
            end do
        end do
    end subroutine

    pure integer function hash_slot(key, hsize)
        integer, intent(in) :: key(4), hsize
        integer(int64) :: hash
        integer :: i
        hash = 0_int64
        do i = 1, 4
            hash = modulo(hash*131_int64 + int(key(i), int64), int(hsize, int64))
        end do
        hash_slot = int(hash) + 1
    end function
end module
