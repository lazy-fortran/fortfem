module fortfem_patch_topology_2d
    !! Shared mesh topology of vertex patches for the planar patchwise
    !! a posteriori estimators (equilibrated fluxes, H(curl) reconstructions):
    !! vertex -> cells incidence, global edge numbering, boundary flags and
    !! the dense LAPACK solve used by the local problems.
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: topology_t, build_topology, shared_edge, dgesv

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(*), info
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
        end subroutine dgesv
    end interface

    !> Mesh topology: vertex -> cells, cell edges, boundary vertices and edges.
    type :: topology_t
        integer, allocatable :: vertex_start(:), vertex_cell(:), vertex_slot(:)
        integer, allocatable :: cell_edge(:, :), cell_orientation(:, :)
        logical, allocatable :: boundary_vertex(:), boundary_edge(:)
    end type topology_t

contains

    subroutine build_topology(nv, triangles, topo)
        !! Edges are numbered by first appearance; local edge e of a cell runs
        !! from its local vertex e to vertex mod(e, 3) + 1 (the FortFEM RT
        !! reference-edge convention) and has orientation +1 when that start
        !! vertex has the smaller global index.
        integer, intent(in) :: nv, triangles(:, :)
        type(topology_t), intent(out) :: topo
        integer, allocatable :: fill(:)
        integer :: nt, t, s, e, a, b, k, other, oe, nedges

        nt = size(triangles, 2)
        allocate (topo%vertex_start(nv + 1), fill(nv))
        topo%vertex_start = 0
        do t = 1, nt
            do s = 1, 3
                a = triangles(s, t)
                topo%vertex_start(a + 1) = topo%vertex_start(a + 1) + 1
            end do
        end do
        topo%vertex_start(1) = 1
        do a = 1, nv
            topo%vertex_start(a + 1) = topo%vertex_start(a + 1) + topo%vertex_start(a)
        end do
        allocate (topo%vertex_cell(3*nt), topo%vertex_slot(3*nt))
        fill = topo%vertex_start(1:nv)
        do t = 1, nt
            do s = 1, 3
                a = triangles(s, t)
                topo%vertex_cell(fill(a)) = t
                topo%vertex_slot(fill(a)) = s
                fill(a) = fill(a) + 1
            end do
        end do

        allocate (topo%cell_edge(3, nt), topo%cell_orientation(3, nt))
        allocate (topo%boundary_vertex(nv), topo%boundary_edge(3*nt))
        topo%cell_edge = 0
        topo%boundary_vertex = .false.
        topo%boundary_edge = .false.
        nedges = 0
        do t = 1, nt
            do e = 1, 3
                a = triangles(e, t)
                b = triangles(mod(e, 3) + 1, t)
                topo%cell_orientation(e, t) = merge(1, -1, a < b)
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
                    topo%boundary_vertex(a) = .true.
                    topo%boundary_vertex(b) = .true.
                end if
                if (other /= 0 .and. other < t) then
                    topo%cell_edge(e, t) = topo%cell_edge(oe, other)
                else
                    nedges = nedges + 1
                    topo%cell_edge(e, t) = nedges
                    topo%boundary_edge(nedges) = other == 0
                end if
            end do
        end do
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

end module fortfem_patch_topology_2d
