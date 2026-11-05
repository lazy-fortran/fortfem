module fortfem_tensor_mesh
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortfem_gll_line, only: gll_rule, lagrange_table
    use fortfem_tensor_line_tables, only: tensor_line_tables_t
    implicit none
    private
    public :: tensor_mesh_t, build_tensor_mesh, structured_tensor_mesh
    public :: ogrid_tensor_mesh, geometry_map, ogrid_boundary
    type :: tensor_mesh_t
        integer :: dim = 0, degree = 0, geometry_degree = 0, nelem = 0, ndof = 0
        integer, allocatable :: vertices(:, :), dof(:, :)
        real(dp), allocatable :: geometry(:, :, :), coordinates(:, :)
        logical, allocatable :: boundary(:)
        ! Line-table discretisations (fortfem_line_mesh): per-direction 1D tables,
        ! the 1D cell of every cell and its patch. Unallocated `cell` selects the
        ! GLL p-nodal layout above, with geometry nodes of `geometry_degree`.
        type(tensor_line_tables_t) :: line(3)
        integer, allocatable :: cell(:, :), patch(:)
    end type
    abstract interface
        subroutine geometry_map(xi, x)
            import dp
            real(dp), intent(in) :: xi(:)
            real(dp), intent(out) :: x(:)
        end subroutine
        subroutine ogrid_boundary(angle, x)
            ! Physical outer-boundary point for a unit-disk polar angle.
            import dp
            real(dp), intent(in) :: angle
            real(dp), intent(out) :: x(2)
        end subroutine
    end interface
contains
    pure subroutine decode_index(index, base, position)
        integer, intent(in) :: index, base
        integer, intent(out) :: position(3)
        position(1) = modulo(index, base)
        position(2) = modulo(index/base, base)
        position(3) = index/base**2
    end subroutine

    subroutine build_tensor_mesh(dim, p, q, vertices, geometry, mesh, status, gather)
        integer, intent(in) :: dim, p, q, vertices(:, :)
        real(dp), intent(in) :: geometry(:, :, :)
        type(tensor_mesh_t), intent(out) :: mesh
        integer, intent(out) :: status
        integer, intent(in), optional :: gather(:, :)
        real(dp), allocatable :: x(:), w(:), z(:), v(:), b(:, :), db(:, :)
        real(dp) :: point(dim), factor
        integer :: e, a, i, j, k, r, s, t, ng, nl, idx, idxg, nd, nv, slot, hsize, l, c
        integer :: pos(3), bits(3), key(16), weight, vid(8), weights(8), swap
        integer, allocatable :: keys(:, :), table(:)
        integer(int64) :: hash
        status = 1
        if (dim < 2 .or. dim > 3 .or. p < 1 .or. q < 1) return
        nv = 2**dim
        if (size(vertices, 1) /= nv .or. size(vertices, 2) < 1) return
        if (any(vertices < 1)) return
        nl = (p + 1)**dim
        ng = (q + 1)**dim
        if (size(geometry, 1) /= dim .or. size(geometry, 2) /= ng) return
        if (size(geometry, 3) /= size(vertices, 2)) return
        mesh%dim = dim
        mesh%degree = p
        mesh%geometry_degree = q
        mesh%nelem = size(vertices, 2)
        mesh%vertices = vertices
        mesh%geometry = geometry
        allocate (mesh%dof(nl, mesh%nelem))
        if (present(gather)) then
            if (size(gather, 1) /= nl .or. size(gather, 2) /= mesh%nelem) return
            if (any(gather < 1)) return
            mesh%dof = gather
            nd = maxval(gather)
        else
            ! An orientation-independent entity key: sorted global corner IDs
            ! with integer tensor barycentric weights. Shared face rotations and
            ! edge reversals produce the same key. Interiors remain cell-local.
            hsize = 2*nl*mesh%nelem + 1
            allocate (keys(16, nl*mesh%nelem), table(hsize))
            table = 0
            nd = 0
            do e = 1, mesh%nelem
                do a = 1, nl
                    call decode_index(a - 1, p + 1, pos)
                    key = 0
                    l = 0
                    do c = 1, nv
                        call decode_index(c - 1, 2, bits)
                        weight = 1
                        do i = 1, dim
                            weight = weight*merge(pos(i), p - pos(i), bits(i) == 1)
                        end do
                        if (weight == 0) cycle
                        l = l + 1
                        vid(l) = vertices(c, e)
                        weights(l) = weight
                    end do
                    do i = 2, l
                        j = i
                        do while (j > 1)
                            if (vid(j) >= vid(j - 1)) exit
                            swap = vid(j); vid(j) = vid(j - 1); vid(j - 1) = swap
                            swap = weights(j); weights(j) = weights(j - 1)
                            weights(j - 1) = swap
                            j = j - 1
                        end do
                    end do
                    do i = 1, l
                        key(2*i - 1) = vid(i)
                        key(2*i) = weights(i)
                    end do
                    hash = 0_int64
                    do i = 1, 16
                        hash = modulo(hash*131_int64 + int(key(i), int64), &
                                      int(hsize, int64))
                    end do
                    slot = int(hash) + 1
                    do
                        if (table(slot) == 0) then
                            nd = nd + 1
                            table(slot) = nd
                            keys(:, nd) = key
                            exit
                        end if
                        if (all(keys(:, table(slot)) == key)) exit
                        slot = modulo(slot, hsize) + 1
                    end do
                    mesh%dof(a, e) = table(slot)
                end do
            end do
        end if
        mesh%ndof = nd
        allocate (mesh%coordinates(dim, nd), mesh%boundary(nd))
        mesh%boundary = .false.
        call gll_rule(p, x, w, status)
        if (status /= 0) return
        call gll_rule(q, z, v, status)
        if (status /= 0) return
        call lagrange_table(z, x, b, db)
        do e = 1, mesh%nelem
            do k = 1, merge(p + 1, 1, dim == 3)
                do j = 1, p + 1
                    do i = 1, p + 1
                        point = 0.0_dp
                        do t = 1, merge(q + 1, 1, dim == 3)
                            do s = 1, q + 1
                                do r = 1, q + 1
                                    idxg = r + (q + 1)*(s - 1) + (q + 1)**2*(t - 1)
                                    factor = b(i, r)*b(j, s)
                                    if (dim == 3) factor = factor*b(k, t)
                                    point = point + factor*geometry(:, idxg, e)
                                end do
                            end do
                        end do
                        idx = i + (p + 1)*(j - 1) + (p + 1)**2*(k - 1)
                        mesh%coordinates(:, mesh%dof(idx, e)) = point
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine

    subroutine structured_tensor_mesh(dim, cells, p, mesh, status, map, periodic, q, &
                                      dirichlet)
        integer, intent(in) :: dim, cells(3), p
        type(tensor_mesh_t), intent(out) :: mesh
        integer, intent(out) :: status
        procedure(geometry_map), optional :: map
        logical, intent(in), optional :: periodic(3)
        integer, intent(in), optional :: q
        !! Non-periodic directions without Dirichlet ends, e.g. an open seam.
        logical, intent(in), optional :: dirichlet(3)
        real(dp), allocatable :: x(:), w(:), geom(:, :, :)
        integer, allocatable :: vertices(:, :), gather(:, :)
        integer :: shape(3), nodes(3), vshape(3), degree, ng, e, i, a, b, k, c, &
                   ec(3), lc(3), vc(3), gc(3), ne, nl
        real(dp) :: xi(dim)
        logical :: per(3), fixed(3)
        status = 1
        if (dim < 2 .or. dim > 3 .or. p < 1) return
        if (any(cells(1:dim) < 1)) return
        fixed = .true.
        if (present(dirichlet)) fixed = dirichlet
        degree = p
        if (present(q)) degree = q
        if (degree < 1) return
        shape = cells
        if (dim == 2) shape(3) = 1
        per = .false.
        if (present(periodic)) per = periodic
        if (dim == 2) per(3) = .false.
        nodes = shape*p + 1
        vshape = shape + 1
        where (per)
            nodes = nodes - 1
            vshape = vshape - 1
        end where
        if (dim == 2) then
            nodes(3) = 1
            vshape(3) = 1
        end if
        ng = (degree + 1)**dim
        nl = (p + 1)**dim
        ne = product(shape)
        allocate (geom(dim, ng, ne), vertices(2**dim, ne), gather(nl, ne))
        call gll_rule(degree, x, w, status)
        if (status /= 0) return
        do e = 1, ne
            ec(1) = modulo(e - 1, shape(1))
            ec(2) = modulo((e - 1)/shape(1), shape(2))
            ec(3) = (e - 1)/(shape(1)*shape(2))
            do c = 1, 2**dim
                call decode_index(c - 1, 2, vc)
                vc = ec + vc
                vc = modulo(vc, vshape)
                vertices(c, e) = 1 + vc(1) + vshape(1)*(vc(2) + vshape(2)*vc(3))
            end do
            do a = 1, nl
                call decode_index(a - 1, p + 1, lc)
                gc = modulo(ec*p + lc, nodes)
                gather(a, e) = 1 + gc(1) + nodes(1)*(gc(2) + nodes(2)*gc(3))
            end do
            do a = 1, ng
                call decode_index(a - 1, degree + 1, lc)
                lc = lc + 1
                do i = 1, dim
                    xi(i) = (real(ec(i), dp) + (x(lc(i)) + 1)/2)/real(shape(i), dp)
                end do
                geom(:, a, e) = xi
                if (present(map)) call map(xi, geom(:, a, e))
            end do
        end do
        call build_tensor_mesh(dim, p, degree, vertices, geom, mesh, status, gather)
        if (status /= 0) return
        do k = 0, nodes(3) - 1
            do b = 0, nodes(2) - 1
                do a = 0, nodes(1) - 1
                    gc(1) = a
                    gc(2) = b
                    gc(3) = k
                    c = 1 + a + nodes(1)*(b + nodes(2)*k)
                    do i = 1, dim
                        if (per(i) .or. .not. fixed(i)) cycle
                    if (gc(i) == 0 .or. gc(i) == nodes(i) - 1) mesh%boundary(c) = .true.
                    end do
                end do
            end do
        end do
    end subroutine

    subroutine ogrid_tensor_mesh(dim, n, p, mesh, status, nphi, q, boundary, &
                                 center, core_half, map, periodic)
        ! Five-block O-grid through the axis. Without boundary: unit disk (2D)
        ! or the deformed reference torus (3D). With boundary (2D only): the
        ! core square is center +/- core_half and each outer block blends its
        ! core edge linearly into boundary(angle) along the radial direction.
        integer, intent(in) :: dim, n, p
        type(tensor_mesh_t), intent(out) :: mesh
        integer, intent(out) :: status
        integer, intent(in), optional :: nphi, q
        procedure(ogrid_boundary), optional :: boundary
        real(dp), intent(in), optional :: center(2), core_half(2)
        procedure(geometry_map), optional :: map
        logical, intent(in), optional :: periodic
        real(dp), allocatable :: x(:), w(:), geom(:, :, :), verts(:, :)
        integer, allocatable :: corners(:, :)
        real(dp) :: uv(3), point(3), angle, radial, tangent, rx, ry, tmp, phi, shape
        real(dp) :: inner(2), outer(2), origin(2), half(2), key(3), turn
        real(dp), parameter :: core = 0.32_dp, pi = acos(-1.0_dp)
        integer :: degree, np, ne, ng, e, block, ix, iy, iz, a, i, c, local(3), nv, v
        logical :: closed
        status = 1
        closed = .true.
        if (present(periodic)) closed = periodic
        if (dim < 2 .or. dim > 3 .or. n < 1 .or. p < 1) return
        degree = p
        if (present(q)) degree = q
        np = 1
        if (dim == 3) np = 8
        if (present(nphi)) then
            if (dim == 3) np = nphi
        end if
        if (degree < 1 .or. np < 1) return
        origin = 0
        half = core
        if (present(boundary)) then
            if (dim /= 2) return
            if (.not. present(center) .or. .not. present(core_half)) return
            if (any(core_half <= 0)) return
            origin = center
            half = core_half
        end if
        ng = (degree + 1)**dim
        ne = 5*n*n*np
        allocate (geom(dim, ng, ne), corners(2**dim, ne), verts(3, 2**dim*ne))
        call gll_rule(degree, x, w, status)
        if (status /= 0) return
        e = 0
        nv = 0
        do iz = 0, np - 1
            do block = 0, 4
                do iy = 0, n - 1
                    do ix = 0, n - 1
                        e = e + 1
                        do a = 1, ng + 2**dim
                            c = 0
                            if (a <= ng) then
                                call decode_index(a - 1, degree + 1, local)
                                do i = 1, 3
                                    uv(i) = (x(local(i) + 1) + 1)/2
                                end do
                            else
                                c = a - ng
                                call decode_index(c - 1, 2, local)
                                uv = real(local, dp)
                            end if
                            radial = (ix + uv(1))/n
                            tangent = 2*(iy + uv(2))/n - 1
                            phi = 2*pi*(iz + uv(3))/np
                            if (block == 0) then
                                rx = core*(2*radial - 1)
                                ry = core*tangent
                                if (present(boundary)) then
                                    rx = origin(1) + half(1)*rx/core
                                    ry = origin(2) + half(2)*ry/core
                                end if
                            else if (present(boundary)) then
                                angle = pi*tangent/4
                                inner = [1.0_dp, tangent]
                                do i = 2, block
                                    tmp = inner(1); inner(1) = -inner(2)
                                    inner(2) = tmp
                                end do
                                inner = origin + half*inner
                                call boundary(angle + (block - 1)*pi/2, outer)
                                rx = (1 - radial)*inner(1) + radial*outer(1)
                                ry = (1 - radial)*inner(2) + radial*outer(2)
                            else
                                angle = pi*tangent/4
                                rx = (1 - radial)*core + radial*cos(angle)
                                ry = (1 - radial)*core*tangent + radial*sin(angle)
                                do i = 2, block
                                    tmp = rx; rx = -ry; ry = tmp
                                end do
                            end if
                            point(1) = rx
                            point(2) = ry
                            point(3) = 0
                            turn = (iz + uv(3))/np
                            ! A mapped grid identifies corners by reference position.
                            key = [rx, ry, 0.0_dp]
                            if (dim == 3) key(3) = turn
                            if (closed .and. key(3) >= 1) key(3) = key(3) - 1
                            if (present(map)) then
                                if (dim == 3) then
                                    call map([rx, ry, turn], point)
                                else
                                    call map([rx, ry], point(1:2))
                                end if
                            else if (dim == 3) then
                                ! Smooth Cartesian deformation through the axis.
                                shape = 1 + 0.08_dp*cos(5*phi)
                                point(1) = (3 + shape*rx)*cos(phi)
                                point(2) = -(3 + shape*rx)*sin(phi)
                                point(3) = (1 - 0.06_dp*cos(5*phi))*ry + &
                                           0.04_dp*sin(5*phi)*rx
                            end if
                            if (a <= ng) then
                                geom(:, a, e) = point(1:dim)
                            else
                                if (.not. present(map)) key = point
                                v = 1
                                do while (v <= nv)
                                    if (maxval(abs(verts(:, v) - key)) &
                                        < 1.e-11_dp) exit
                                    v = v + 1
                                end do
                                if (v > nv) then
                                    nv = v
                                    verts(:, v) = key
                                end if
                                corners(c, e) = v
                            end if
                        end do
                    end do
                end do
            end do
        end do
        call build_tensor_mesh(dim, p, degree, corners, geom, mesh, status)
        if (status /= 0) return
        ! Only the outer ring is Dirichlet. Block seams and the axis are interior.
        e = 0
        do iz = 0, np - 1
            do block = 0, 4
                do iy = 0, n - 1
                    do ix = 0, n - 1
                        e = e + 1
                        if (block == 0 .or. ix /= n - 1) cycle
                        do a = 1, (p + 1)**dim
                            if (modulo(a - 1, p + 1) == p) &
                                mesh%boundary(mesh%dof(a, e)) = .true.
                        end do
                    end do
                end do
            end do
        end do
    end subroutine
end module
