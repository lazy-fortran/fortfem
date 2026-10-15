module fortfem_equilibrated_flux_topology
    !! Mesh topology shared by the patchwise a posteriori estimators on
    !! simplicial meshes (2D triangles, 3D tetrahedra, meridian and H(curl)
    !! estimators): vertex stars (vertex -> cells), planar edge numbering and
    !! boundary flags, status codes and the dense LAPACK solve of the local
    !! problems.
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: vertex_stars, topology_t, build_topology, shared_edge, dgesv

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(*), info
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
        end subroutine dgesv
    end interface

    !> Triangle mesh topology: vertex -> cells, cell edges, boundary flags.
    type :: topology_t
        integer, allocatable :: vertex_start(:), vertex_cell(:), vertex_slot(:)
        integer, allocatable :: cell_edge(:, :), cell_orientation(:, :)
        !> Vertex lies on a boundary edge that is not a Neumann edge.
        logical, allocatable :: dirichlet_vertex(:)
        !> Global edge e lies on the mesh boundary.
        logical, allocatable :: boundary_edge(:)
    end type topology_t

    !> Status codes shared by the 2D and 3D equilibrated-flux estimators.
    integer, parameter, public :: EQUILIBRATION_OK = 0
    integer, parameter, public :: EQUILIBRATION_BAD_INPUT = 1
    integer, parameter, public :: EQUILIBRATION_BAD_CELL = 2
    integer, parameter, public :: EQUILIBRATION_PATCH_FAILED = 3

contains

    subroutine vertex_stars(nv, cells, vertex_start, vertex_cell, vertex_slot)
        !! Compressed vertex -> (cell, local slot) lists: the entries of vertex
        !! a are vertex_start(a) .. vertex_start(a + 1) - 1, in cell order.
        integer, intent(in) :: nv, cells(:, :)
        integer, allocatable, intent(out) :: vertex_start(:), vertex_cell(:)
        integer, allocatable, intent(out) :: vertex_slot(:)
        integer, allocatable :: fill(:)
        integer :: nt, nlocal, t, s, a

        nt = size(cells, 2)
        nlocal = size(cells, 1)
        allocate (vertex_start(nv + 1), fill(nv))
        vertex_start = 0
        do t = 1, nt
            do s = 1, nlocal
                a = cells(s, t)
                vertex_start(a + 1) = vertex_start(a + 1) + 1
            end do
        end do
        vertex_start(1) = 1
        do a = 1, nv
            vertex_start(a + 1) = vertex_start(a + 1) + vertex_start(a)
        end do
        allocate (vertex_cell(nlocal*nt), vertex_slot(nlocal*nt))
        fill = vertex_start(1:nv)
        do t = 1, nt
            do s = 1, nlocal
                a = cells(s, t)
                vertex_cell(fill(a)) = t
                vertex_slot(fill(a)) = s
                fill(a) = fill(a) + 1
            end do
        end do
    end subroutine vertex_stars

    subroutine build_topology(nv, triangles, topo, status, neumann)
        !! Edges are numbered by first appearance; local edge e of a cell runs
        !! from its local vertex e to vertex mod(e, 3) + 1 (the FortFEM RT
        !! reference-edge convention) and has orientation +1 when that start
        !! vertex has the smaller global index.
        !! A vertex is a Dirichlet vertex when it lies on a boundary edge that is
        !! not flagged in the optional `neumann(3, nt)`; a flagged edge must be a
        !! boundary edge (status 1 otherwise).
        integer, intent(in) :: nv, triangles(:, :)
        type(topology_t), intent(out) :: topo
        integer, intent(out) :: status
        logical, intent(in), optional :: neumann(:, :)
        logical :: flagged
        integer :: nt, t, e, a, b, k, other, oe, nedges

        status = 1
        nt = size(triangles, 2)
        call vertex_stars(nv, triangles, topo%vertex_start, topo%vertex_cell, &
            topo%vertex_slot)

        allocate (topo%cell_edge(3, nt), topo%cell_orientation(3, nt))
        allocate (topo%dirichlet_vertex(nv), topo%boundary_edge(3*nt))
        topo%cell_edge = 0
        topo%dirichlet_vertex = .false.
        topo%boundary_edge = .false.
        nedges = 0
        do t = 1, nt
            do e = 1, 3
                a = triangles(e, t)
                b = triangles(mod(e, 3) + 1, t)
                topo%cell_orientation(e, t) = merge(1, -1, a < b)
                flagged = .false.
                if (present(neumann)) flagged = neumann(e, t)
                other = 0
                oe = 0
                do k = topo%vertex_start(a), topo%vertex_start(a + 1) - 1
                    if (topo%vertex_cell(k) == t) cycle
                    oe = shared_edge(triangles(:, topo%vertex_cell(k)), a, b)
                    if (oe /= 0) then
                        other = topo%vertex_cell(k)
                        exit
                    end if
                end do
                if (other == 0) then
                    if (.not. flagged) then
                        topo%dirichlet_vertex(a) = .true.
                        topo%dirichlet_vertex(b) = .true.
                    end if
                else
                    if (flagged) return
                end if
                if (other /= 0) then
                    if (other < t) then
                        topo%cell_edge(e, t) = topo%cell_edge(oe, other)
                        cycle
                    end if
                end if
                nedges = nedges + 1
                topo%cell_edge(e, t) = nedges
                topo%boundary_edge(nedges) = other == 0
            end do
        end do
        status = 0
    end subroutine build_topology

    pure integer function shared_edge(cell, a, b) result(edge)
        !! Local edge of `cell` joining global vertices a and b, or 0.
        integer, intent(in) :: cell(3), a, b
        integer :: e, p, q

        edge = 0
        do e = 1, 3
            p = cell(e)
            q = cell(mod(e, 3) + 1)
            if ((p == a .and. q == b) .or. (p == b .and. q == a)) then
                edge = e
                return
            end if
        end do
    end function shared_edge

end module fortfem_equilibrated_flux_topology
