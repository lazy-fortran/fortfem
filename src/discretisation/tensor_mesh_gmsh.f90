module fortfem_tensor_mesh_gmsh
    ! Gmsh MSH 4.1 ASCII reader for curved quadrilateral (dim 2) and hexahedral
    ! (dim 3) meshes of geometry order 1..3 into the tensor gather-map mesh.
    ! Equispaced Gmsh nodes are converted exactly to the GLL geometry nodes of the
    ! same polynomial map, cells are oriented to a positive Jacobian, $Periodic
    ! node links identify corner vertices, and the Dirichlet set is the geometric
    ! boundary (faces used by one cell). Physical groups are not interpreted.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_mesh, only: tensor_mesh_t, build_tensor_mesh
    use fortfem_gll_line, only: gll_rule, lagrange_table
    use fortfem_tensor_mesh_tools, only: orient_tensor_cells, mark_boundary_faces
    implicit none
    private
    public :: read_gmsh_tensor_mesh, gmsh_lattice
    ! Lattice positions of Gmsh nodes, from gmsh.model.mesh.getElementProperties
    ! (Gmsh 4.15.2): reference coordinate (c + 1)/2*order per node.
    integer, parameter :: quad4(2, 4) = reshape([ &
        0, 0, 1, 0, 1, 1, 0, 1], [2, 4])
    integer, parameter :: quad9(2, 9) = reshape([ &
        0, 0, 2, 0, 2, 2, 0, 2, 1, 0, 2, 1, 1, 2, 0, 1, 1, 1], [2, 9])
    integer, parameter :: quad16(2, 16) = reshape([ &
        0, 0, 3, 0, 3, 3, 0, 3, 1, 0, 2, 0, 3, 1, 3, 2, 2, 3, 1, 3, 0, 2, 0, &
        1, 1, 1, 2, 1, 2, 2, 1, 2], [2, 16])
    integer, parameter :: hex8(3, 8) = reshape([ &
        0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1], [3, 8])
    integer, parameter :: hex27(3, 27) = reshape([ &
        0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 2, 0, 0, 0, 2, 2, 0, 2, 2, 2, 2, 0, 2, &
        2, 1, 0, 0, 0, 1, 0, 0, 0, 1, 2, 1, 0, 2, 0, 1, 1, 2, 0, 2, 2, 1, 0, &
        2, 1, 1, 0, 2, 0, 1, 2, 2, 1, 2, 1, 2, 2, 1, 1, 0, 1, 0, 1, 0, 1, 1, &
        2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1], [3, 27])
    integer, parameter :: hex64(3, 64) = reshape([ &
        0, 0, 0, 3, 0, 0, 3, 3, 0, 0, 3, 0, 0, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, &
        3, 1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 3, 1, 0, 3, &
        2, 0, 3, 0, 1, 3, 0, 2, 2, 3, 0, 1, 3, 0, 3, 3, 1, 3, 3, 2, 0, 3, 1, &
        0, 3, 2, 1, 0, 3, 2, 0, 3, 0, 1, 3, 0, 2, 3, 3, 1, 3, 3, 2, 3, 2, 3, &
        3, 1, 3, 3, 1, 1, 0, 1, 2, 0, 2, 2, 0, 2, 1, 0, 1, 0, 1, 2, 0, 1, 2, &
        0, 2, 1, 0, 2, 0, 1, 1, 0, 1, 2, 0, 2, 2, 0, 2, 1, 3, 1, 1, 3, 2, 1, &
        3, 2, 2, 3, 1, 2, 2, 3, 1, 1, 3, 1, 1, 3, 2, 2, 3, 2, 1, 1, 3, 2, 1, &
        3, 2, 2, 3, 1, 2, 3, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 1, 1, 2, 2, &
        1, 2, 2, 2, 2, 1, 2, 2], [3, 64])
contains
    subroutine gmsh_lattice(element_type, dim, order, position, status)
        integer, intent(in) :: element_type
        integer, intent(out) :: dim, order, status
        integer, allocatable, intent(out) :: position(:, :)
        status = 0
        select case (element_type)
        case (3)
            dim = 2; order = 1; position = quad4
        case (10)
            dim = 2; order = 2; position = quad9
        case (36)
            dim = 2; order = 3; position = quad16
        case (5)
            dim = 3; order = 1; position = hex8
        case (12)
            dim = 3; order = 2; position = hex27
        case (92)
            dim = 3; order = 3; position = hex64
        case default
            dim = 0; order = 0; status = 1
            allocate (position(0, 0))
        end select
    end subroutine

    subroutine read_gmsh_tensor_mesh(path, p, mesh, status, dim, message)
        character(len=*), intent(in) :: path
        integer, intent(in) :: p
        type(tensor_mesh_t), intent(out) :: mesh
        integer, intent(out) :: status
        integer, intent(in), optional :: dim
        character(len=:), allocatable, intent(out), optional :: message
        character(len=4096) :: line
        real(dp), allocatable :: coords(:, :), equi(:), gll(:), w(:), b(:, :), &
                                 d(:, :), geom(:, :, :)
        integer, allocatable :: tags(:), conn(:, :), lat(:, :), master(:), &
                                corners(:, :), vertex(:), position(:, :), types(:)
        real(dp) :: version, factor
        integer :: u, ios, nblock, nnode, minimum, maximum, blk, edim, etag, param, &
                   count, i, j, a, g, e, kind, nelem, target, order, edimension, ne, &
                   nl, nlink, naff, npair, tag, tagm, lorder, ldim, flipped, nv, &
                   pg(3), pl(3), filetype, datasize, nn, c, pc(3)
        character(len=128) :: why
        status = 1
        why = ''
        target = 0
        if (present(dim)) target = dim
        allocate (coords(3, 0), conn(0, 0), types(0), master(0))
        nelem = 0
        order = 0
        edimension = 0
        open (newunit=u, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            why = 'cannot open file'
            call finish()
            return
        end if
        do
            read (u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            select case (trim(adjustl(line)))
            case ('$MeshFormat')
                read (u, *, iostat=ios) version, filetype, datasize
                if (ios /= 0 .or. abs(version - 4.1_dp) > 1.e-6_dp &
                    .or. filetype /= 0) then
                    why = 'only MSH 4.1 ASCII is supported'
                    exit
                end if
            case ('$Nodes')
                read (u, *, iostat=ios) nblock, nnode, minimum, maximum
                if (ios /= 0 .or. minimum < 1) then
                    why = 'invalid $Nodes header'
                    exit
                end if
                deallocate (coords, master)
                allocate (coords(3, maximum), master(maximum))
                coords = 0
                do i = 1, maximum
                    master(i) = i
                end do
                do blk = 1, nblock
                    read (u, *, iostat=ios) edim, etag, param, count
                    if (ios /= 0) exit
                    if (allocated(tags)) deallocate (tags)
                    allocate (tags(count))
                    do i = 1, count
                        read (u, *, iostat=ios) tags(i)
                    end do
                    do i = 1, count
                        read (u, *, iostat=ios) coords(:, tags(i))
                    end do
                end do
                if (ios /= 0) then
                    why = 'invalid $Nodes block'
                    exit
                end if
            case ('$Elements')
                read (u, *, iostat=ios) nblock, ne, minimum, maximum
                if (ios /= 0) then
                    why = 'invalid $Elements header'
                    exit
                end if
                deallocate (conn, types)
                allocate (conn(64, ne), types(ne))
                conn = 0
                do blk = 1, nblock
                    read (u, *, iostat=ios) edim, etag, kind, count
                    if (ios /= 0) exit
                    call gmsh_lattice(kind, ldim, lorder, position, j)
                    nn = size(position, 2)
                    do i = 1, count
                        if (j == 0) then
                            read (u, *, iostat=ios) tag, conn(1:nn, nelem + 1)
                            nelem = nelem + 1
                            types(nelem) = kind
                        else
                            read (u, '(a)', iostat=ios) line
                        end if
                    end do
                end do
                if (ios /= 0) then
                    why = 'invalid $Elements block'
                    exit
                end if
            case ('$Periodic')
                read (u, *, iostat=ios) nlink
                do i = 1, nlink
                    read (u, *, iostat=ios) edim, etag, tag
                    read (u, *, iostat=ios) naff
                    read (u, *, iostat=ios) npair
                    do j = 1, npair
                        read (u, *, iostat=ios) tag, tagm
                        if (ios /= 0) exit
                        if (tag >= 1 .and. tag <= size(master) .and. tagm >= 1 &
                            .and. tagm <= size(master)) master(tag) = tagm
                    end do
                end do
                if (ios /= 0) then
                    why = 'invalid $Periodic block'
                    exit
                end if
            end select
        end do
        close (u)
        if (len_trim(why) > 0) then
            call finish()
            return
        end if
        ! Highest quad/hex dimension present, unless the caller selects one.
        do e = 1, nelem
            call gmsh_lattice(types(e), ldim, lorder, position, j)
            if (target == 0) edimension = max(edimension, ldim)
        end do
        if (target /= 0) edimension = target
        ne = 0
        do e = 1, nelem
            call gmsh_lattice(types(e), ldim, lorder, position, j)
            if (ldim /= edimension) cycle
            if (order == 0) order = lorder
            if (lorder /= order) then
                why = 'mixed geometry orders'
                call finish()
                return
            end if
            ne = ne + 1
        end do
        if (ne == 0 .or. edimension < 2) then
            why = 'no quadrilateral or hexahedral cells'
            call finish()
            return
        end if
        ! Resolve periodic chains so every slave names its final master.
        do i = 1, size(master)
            j = master(i)
            do a = 1, size(master)
                if (master(j) == j) exit
                j = master(j)
            end do
            master(i) = j
        end do
        nl = (order + 1)**edimension
        allocate (geom(edimension, nl, ne), corners(2**edimension, ne), &
                  vertex(size(master)), equi(order + 1))
        vertex = 0
        nv = 0
        do i = 0, order
            equi(i + 1) = -1 + 2*real(i, dp)/order
        end do
        call gll_rule(order, gll, w, status)
        if (status /= 0) return
        status = 1
        call lagrange_table(equi, gll, b, d)
        allocate (lat(edimension, nl))
        j = 0
        do e = 1, nelem
            call gmsh_lattice(types(e), ldim, lorder, position, a)
            if (ldim /= edimension) cycle
            j = j + 1
            lat = position
            geom(:, :, j) = 0
            do g = 1, nl
                pg = 0
                pg(1:edimension) = lat(:, g)
                ! Equispaced node value contributes to every GLL node.
                do a = 1, nl
                    pl = 0
                    call decode(a - 1, order + 1, pl)
                    factor = 1
                    do i = 1, edimension
                        factor = factor*b(pl(i) + 1, pg(i) + 1)
                    end do
                    geom(:, a, j) = geom(:, a, j) + &
                                    factor*coords(1:edimension, conn(g, e))
                end do
            end do
            do c = 1, 2**edimension
                call decode(c - 1, 2, pc)
                do g = 1, 2**edimension
                    if (all(lat(:, g) == pc(1:edimension)*order)) exit
                end do
                tag = master(conn(g, e))
                if (vertex(tag) == 0) then
                    nv = nv + 1
                    vertex(tag) = nv
                end if
                corners(c, j) = vertex(tag)
            end do
        end do
        call orient_tensor_cells(edimension, order, corners, geom, flipped)
        call build_tensor_mesh(edimension, p, order, corners, geom, mesh, status)
        if (status /= 0) then
            why = 'build_tensor_mesh rejected the cells'
            call finish()
            return
        end if
        call mark_boundary_faces(mesh)
        status = 0
        call finish()
    contains
        subroutine finish()
            if (present(message)) message = trim(why)
        end subroutine
        pure subroutine decode(index, base, position)
            integer, intent(in) :: index, base
            integer, intent(out) :: position(3)
            position(1) = modulo(index, base)
            position(2) = modulo(index/base, base)
            position(3) = index/base**2
        end subroutine
    end subroutine
end module
