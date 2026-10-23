module fortfem_tensor_mesh_tools
    ! Generic checks and transformations for unstructured curved quad/hex meshes
    ! given as corner IDs plus tensor GLL geometry nodes (fortfem_tensor_mesh).
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortfem_tensor_mesh, only: tensor_mesh_t, build_tensor_mesh, geometry_map
    use fortfem_gll_line, only: gll_rule, lagrange_table
    use fortfem_tensor_operator, only: determinant_at
    implicit none
    private
    public :: tensor_mesh_quality_t, orient_tensor_cells, mark_boundary_faces
    public :: tensor_mesh_check, extrude_tensor_mesh
    type :: tensor_mesh_quality_t
        real(dp) :: min_jacobian = 0, min_scaled_jacobian = 0, max_gap = 0
        integer :: boundary_faces = 0, inverted_cells = 0
    end type
contains
    pure subroutine lattice(index, base, position)
        integer, intent(in) :: index, base
        integer, intent(out) :: position(3)
        position(1) = modulo(index, base)
        position(2) = modulo(index/base, base)
        position(3) = index/base**2
    end subroutine

    subroutine jacobian_table(dim, q, points, jac)
        ! jac(comp, dir, point, node): d x_comp / d xi_dir contribution of a node.
        integer, intent(in) :: dim, q
        real(dp), intent(in) :: points(:)
        real(dp), allocatable, intent(out) :: jac(:, :, :)
        real(dp), allocatable :: z(:), w(:), b(:, :), d(:, :)
        integer :: status, np, n, a, s, pa(3), ps(3), dir, i
        real(dp) :: factor
        call gll_rule(q, z, w, status)
        call lagrange_table(z, points, b, d)
        np = size(points)**dim
        n = (q + 1)**dim
        allocate (jac(dim, np, n))
        do a = 1, n
            call lattice(a - 1, q + 1, pa)
            do s = 1, np
                call lattice(s - 1, size(points), ps)
                do dir = 1, dim
                    factor = 1
                    do i = 1, dim
                        if (i == dir) then
                            factor = factor*d(ps(i) + 1, pa(i) + 1)
                        else
                            factor = factor*b(ps(i) + 1, pa(i) + 1)
                        end if
                    end do
                    jac(dir, s, a) = factor
                end do
            end do
        end do
    end subroutine

    function cell_jacobian(dim, table, geometry, s) result(j)
        integer, intent(in) :: dim, s
        real(dp), intent(in) :: table(:, :, :), geometry(:, :)
        real(dp) :: j(3, 3)
        integer :: a, dir
        j = 0
        if (dim == 2) j(3, 3) = 1
        do a = 1, size(geometry, 2)
            do dir = 1, dim
                j(1:dim, dir) = j(1:dim, dir) + table(dir, s, a)*geometry(:, a)
            end do
        end do
    end function

    subroutine orient_tensor_cells(dim, q, corners, geometry, flipped)
        ! Swap the first two reference directions of every cell whose Jacobian
        ! is negative at the cell centre. The physical cell is unchanged.
        integer, intent(in) :: dim, q
        integer, intent(inout) :: corners(:, :)
        real(dp), intent(inout) :: geometry(:, :, :)
        integer, intent(out) :: flipped
        real(dp), allocatable :: table(:, :, :), old(:, :)
        integer, allocatable :: oldc(:)
        integer :: e, a, pa(3), c, pc(3)
        call jacobian_table(dim, q, [0.0_dp], table)
        flipped = 0
        do e = 1, size(geometry, 3)
            if (determinant_at(cell_jacobian(dim, table, geometry(:, :, e), 1)) &
                >= 0) cycle
            flipped = flipped + 1
            old = geometry(:, :, e)
            do a = 1, size(geometry, 2)
                call lattice(a - 1, q + 1, pa)
                geometry(:, a, e) = old(:, 1 + pa(2) + (q + 1)*pa(1) + &
                                        (q + 1)**2*pa(3))
            end do
            oldc = corners(:, e)
            do c = 1, size(corners, 1)
                call lattice(c - 1, 2, pc)
                corners(c, e) = oldc(1 + pc(2) + 2*pc(1) + 4*pc(3))
            end do
        end do
    end subroutine

    subroutine face_corners(dim, dir, side, list, count)
        integer, intent(in) :: dim, dir, side
        integer, intent(out) :: list(4), count
        integer :: c, pc(3)
        count = 0
        list = 0
        do c = 1, 2**dim
            call lattice(c - 1, 2, pc)
            if (pc(dir) /= side) cycle
            count = count + 1
            list(count) = c
        end do
    end subroutine

    subroutine mark_boundary_faces(mesh, faces)
        ! Geometric boundary: every cell face whose sorted corner IDs occur once.
        type(tensor_mesh_t), intent(inout) :: mesh
        integer, intent(out), optional :: faces
        integer, allocatable :: keys(:, :), table(:), count(:)
        integer :: e, dir, side, list(4), nc, key(4), i, j, swap, slot, hsize, nkey, &
                   a, pa(3), n, nboundary, id
        integer(int64) :: hash
        n = mesh%degree + 1
        hsize = 4*mesh%dim*mesh%nelem + 1
        allocate (keys(4, 2*mesh%dim*mesh%nelem), table(hsize), &
                  count(2*mesh%dim*mesh%nelem))
        table = 0
        count = 0
        nkey = 0
        nboundary = 0
        do i = 1, 2
            do e = 1, mesh%nelem
                do dir = 1, mesh%dim
                    do side = 0, 1
                        call face_corners(mesh%dim, dir, side, list, nc)
                        key = 0
                        key(1:nc) = mesh%vertices(list(1:nc), e)
                        do a = 2, nc
                            j = a
                            do while (j > 1)
                                if (key(j) >= key(j - 1)) exit
                                swap = key(j); key(j) = key(j - 1); key(j - 1) = swap
                                j = j - 1
                            end do
                        end do
                        hash = 0
                        do a = 1, 4
                            hash = modulo(hash*131_int64 + int(key(a), int64), &
                                          int(hsize, int64))
                        end do
                        slot = int(hash) + 1
                        do
                            id = table(slot)
                            if (id == 0) then
                                nkey = nkey + 1
                                keys(:, nkey) = key
                                table(slot) = nkey
                                id = nkey
                                exit
                            end if
                            if (all(keys(:, id) == key)) exit
                            slot = modulo(slot, hsize) + 1
                        end do
                        if (i == 1) then
                            count(id) = count(id) + 1
                        else if (count(id) == 1) then
                            nboundary = nboundary + 1
                            do a = 1, n**mesh%dim
                                call lattice(a - 1, n, pa)
                                if (pa(dir) == side*(n - 1)) &
                                    mesh%boundary(mesh%dof(a, e)) = .true.
                            end do
                        end if
                    end do
                end do
            end do
        end do
        if (present(faces)) faces = nboundary
    end subroutine

    subroutine tensor_mesh_check(mesh, quality)
        ! Jacobian sign/shape on a (q+2)-point GLL lattice per cell and conformity:
        ! every cell must reproduce the coordinates of its shared dof nodes.
        ! A periodic identification reports its period as max_gap.
        type(tensor_mesh_t), intent(in) :: mesh
        type(tensor_mesh_quality_t), intent(out) :: quality
        real(dp), allocatable :: z(:), w(:), table(:, :, :), x(:), wx(:), b(:, :), &
                                 d(:, :), zg(:)
        real(dp) :: j(3, 3), det, scale, point(3), factor, low
        integer :: e, s, status, q, dim, a, g, pa(3), pg(3), i
        dim = mesh%dim
        q = mesh%geometry_degree
        call gll_rule(q + 1, z, w, status)
        call jacobian_table(dim, q, z, table)
        quality%min_jacobian = huge(1.0_dp)
        quality%min_scaled_jacobian = huge(1.0_dp)
        quality%max_gap = 0
        call gll_rule(mesh%degree, x, wx, status)
        call gll_rule(q, zg, w, status)
        call lagrange_table(zg, x, b, d)
        do e = 1, mesh%nelem
            low = huge(1.0_dp)
            do s = 1, size(z)**dim
                j = cell_jacobian(dim, table, mesh%geometry(:, :, e), s)
                det = determinant_at(j)
                scale = 1
                do i = 1, dim
                    scale = scale*norm2(j(1:dim, i))
                end do
                low = min(low, det)
                quality%min_scaled_jacobian = min(quality%min_scaled_jacobian, &
                                                  det/max(scale, tiny(1.0_dp)))
            end do
            quality%min_jacobian = min(quality%min_jacobian, low)
            if (low <= 0) quality%inverted_cells = quality%inverted_cells + 1
            do a = 1, (mesh%degree + 1)**dim
                call lattice(a - 1, mesh%degree + 1, pa)
                point = 0
                do g = 1, (q + 1)**dim
                    call lattice(g - 1, q + 1, pg)
                    factor = 1
                    do i = 1, dim
                        factor = factor*b(pa(i) + 1, pg(i) + 1)
                    end do
                    point(1:dim) = point(1:dim) + factor*mesh%geometry(:, g, e)
                end do
                quality%max_gap = max(quality%max_gap, maxval(abs(point(1:dim) - &
                                                  mesh%coordinates(:, mesh%dof(a, e)))))
            end do
        end do
        if (mesh%nelem == 0) then
            quality%min_jacobian = 0
            quality%min_scaled_jacobian = 0
        end if
    end subroutine

    subroutine extrude_tensor_mesh(base, layers, p, mesh, status, map, periodic)
        ! Quad x phi: a 2D curved quad mesh times a 1D direction s in [0, 1]
        ! with `layers` cells, optionally periodic. The optional map receives
        ! (x, y, s) and returns the 3D point; default is the straight prism.
        type(tensor_mesh_t), intent(in) :: base
        integer, intent(in) :: layers, p
        type(tensor_mesh_t), intent(out) :: mesh
        integer, intent(out) :: status
        procedure(geometry_map), optional :: map
        logical, intent(in), optional :: periodic
        real(dp), allocatable :: z(:), w(:), geom(:, :, :)
        integer, allocatable :: corners(:, :)
        integer :: q, nv, ne, k, e, c, pc(3), a, pa(3), plane, flipped, layer
        real(dp) :: xi(3)
        logical :: closed
        status = 1
        if (base%dim /= 2 .or. layers < 1 .or. p < 1) return
        closed = .false.
        if (present(periodic)) closed = periodic
        if (closed .and. layers < 3) return
        q = base%geometry_degree
        nv = maxval(base%vertices)
        ne = base%nelem
        call gll_rule(q, z, w, status)
        if (status /= 0) return
        allocate (corners(8, ne*layers), geom(3, (q + 1)**3, ne*layers))
        do k = 0, layers - 1
            do e = 1, ne
                do c = 1, 8
                    call lattice(c - 1, 2, pc)
                    layer = k + pc(3)
                    if (closed) layer = modulo(layer, layers)
                    corners(c, e + ne*k) = base%vertices(1 + pc(1) + 2*pc(2), e) + &
                                           nv*layer
                end do
                do a = 1, (q + 1)**3
                    call lattice(a - 1, q + 1, pa)
                    plane = 1 + pa(1) + (q + 1)*pa(2)
                    xi(1:2) = base%geometry(:, plane, e)
                    xi(3) = (k + (z(pa(3) + 1) + 1)/2)/layers
                    geom(:, a, e + ne*k) = xi
                    if (present(map)) call map(xi, geom(:, a, e + ne*k))
                end do
            end do
        end do
        call orient_tensor_cells(3, q, corners, geom, flipped)
        call build_tensor_mesh(3, p, q, corners, geom, mesh, status)
        if (status /= 0) return
        call mark_boundary_faces(mesh)
    end subroutine
end module
