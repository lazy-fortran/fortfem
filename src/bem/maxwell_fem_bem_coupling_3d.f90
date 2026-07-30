module fortfem_maxwell_fem_bem_coupling_3d
    !! Conforming pullback of the exterior RWG EFIE operator to the global
    !! lowest-order tetrahedral Nedelec edge space.
    use fortfem_kinds, only: dp
    use fortfem_maxwell_efie_rwg_3d, only: assemble_maxwell_efie_rwg_3d
    use fortfem_maxwell_rwg_surface, only: &
        build_maxwell_rwg_surface_space, &
        map_maxwell_rwg_to_tetra_nedelec_edges
    use fortfem_tetra_edge_dof_map, only: build_tetra_edge_dof_map
    implicit none
    private

    public :: assemble_maxwell_fem_bem_boundary_matrix_3d

contains

    subroutine assemble_maxwell_fem_bem_boundary_matrix_3d( &
            vertices, tetrahedra, boundary_triangles, wave_number, impedance, &
            quadrature_degree, tolerance, max_depth, matrix, status)
        real(dp), intent(in) :: vertices(:, :), wave_number, impedance
        integer, intent(in) :: tetrahedra(:, :), boundary_triangles(:, :)
        integer, intent(in) :: quadrature_degree, max_depth
        real(dp), intent(in) :: tolerance
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: efie(:, :)
        integer, allocatable :: edge_dofs(:, :), edge_orientations(:, :)
        integer, allocatable :: edges(:, :), rwg_dofs(:)
        integer, allocatable :: rwg_edge_triangles(:, :)
        integer, allocatable :: rwg_edges(:, :)
        real(dp), allocatable :: trace_scales(:)
        integer :: column, row

        status = 1
        call build_tetra_edge_dof_map( &
            tetrahedra, edges, edge_dofs, edge_orientations, status)
        if (status /= 0) return
        call build_maxwell_rwg_surface_space( &
            vertices, boundary_triangles, rwg_edges, rwg_edge_triangles, status)
        if (status /= 0 .or. size(rwg_edges, 2) == 0) return
        call map_maxwell_rwg_to_tetra_nedelec_edges( &
            vertices, tetrahedra, rwg_edges, rwg_dofs, trace_scales, status)
        if (status /= 0) return
        call assemble_maxwell_efie_rwg_3d( &
            vertices, boundary_triangles, wave_number, impedance, &
            quadrature_degree, tolerance, max_depth, efie, status)
        if (status /= 0) return
        allocate(matrix(size(edges, 2), size(edges, 2)))
        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        do column = 1, size(rwg_dofs)
            do row = 1, size(rwg_dofs)
                matrix(rwg_dofs(row), rwg_dofs(column)) = &
                    matrix(rwg_dofs(row), rwg_dofs(column)) + &
                    trace_scales(row)*trace_scales(column)*efie(row, column)
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_fem_bem_boundary_matrix_3d

end module fortfem_maxwell_fem_bem_coupling_3d
