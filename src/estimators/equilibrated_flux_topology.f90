module fortfem_equilibrated_flux_topology
    !! Mesh topology shared by the patchwise equilibrated-flux estimators on
    !! simplicial meshes (2D triangles and 3D tetrahedra).
    implicit none
    private

    public :: vertex_stars

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

end module fortfem_equilibrated_flux_topology
