module fortfem_patch_topology
    !! Vertex patches of a triangle mesh in compressed-row form, shared by the
    !! patchwise equilibrated estimators.
    implicit none
    private

    public :: build_vertex_patches

contains

    subroutine build_vertex_patches(nv, triangles, vertex_start, vertex_cell, &
            vertex_slot)
        !! The triangles around vertex a are vertex_cell(vertex_start(a) :
        !! vertex_start(a + 1) - 1) in increasing triangle order; vertex_slot
        !! is the local index (1..3) of a in that triangle.
        integer, intent(in) :: nv, triangles(:, :)
        integer, allocatable, intent(out) :: vertex_start(:), vertex_cell(:), &
            vertex_slot(:)
        integer, allocatable :: fill(:)
        integer :: nt, t, s, a

        nt = size(triangles, 2)
        allocate (vertex_start(nv + 1), fill(nv))
        vertex_start = 0
        do t = 1, nt
            do s = 1, 3
                a = triangles(s, t)
                vertex_start(a + 1) = vertex_start(a + 1) + 1
            end do
        end do
        vertex_start(1) = 1
        do a = 1, nv
            vertex_start(a + 1) = vertex_start(a + 1) + vertex_start(a)
        end do
        allocate (vertex_cell(3*nt), vertex_slot(3*nt))
        fill = vertex_start(1:nv)
        do t = 1, nt
            do s = 1, 3
                a = triangles(s, t)
                vertex_cell(fill(a)) = t
                vertex_slot(fill(a)) = s
                fill(a) = fill(a) + 1
            end do
        end do
    end subroutine build_vertex_patches

end module fortfem_patch_topology
