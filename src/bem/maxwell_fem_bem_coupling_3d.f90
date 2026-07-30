module fortfem_maxwell_fem_bem_coupling_3d
    !! Conforming pullback of the exterior RWG EFIE operator to the global
    !! lowest-order tetrahedral Nedelec edge space.
    use fortfem_kinds, only: dp
    use fortfem_assembly_tetra_nedelec_3d, only: &
        assemble_tetra_nedelec_curl_mass_element
    use fortfem_maxwell_efie_rwg_3d, only: assemble_maxwell_efie_rwg_3d
    use fortfem_maxwell_rwg_surface, only: &
        assemble_maxwell_rwg_mass_matrix, &
        build_maxwell_rwg_surface_space, &
        map_maxwell_rwg_to_tetra_nedelec_edges
    use fortfem_tetra_edge_dof_map, only: build_tetra_edge_dof_map
    implicit none
    private

    public :: assemble_maxwell_fem_bem_boundary_matrix_3d
    public :: assemble_maxwell_fem_bem_system_3d
    public :: solve_maxwell_fem_bem_system_3d

    interface
        subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            complex(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*)
            complex(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine zgesv
    end interface

contains

    subroutine solve_maxwell_fem_bem_system_3d( &
            vertices, tetrahedra, boundary_triangles, curl_coefficient, &
            mass_coefficient, wave_number, impedance, quadrature_degree, &
            tolerance, max_depth, right_hand_side, field, current, status)
        real(dp), intent(in) :: vertices(:, :), curl_coefficient
        real(dp), intent(in) :: mass_coefficient, wave_number, impedance
        integer, intent(in) :: tetrahedra(:, :), boundary_triangles(:, :)
        integer, intent(in) :: quadrature_degree, max_depth
        real(dp), intent(in) :: tolerance
        complex(dp), intent(in) :: right_hand_side(:)
        complex(dp), allocatable, intent(out) :: field(:), current(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: matrix(:, :), solution(:, :)
        integer, allocatable :: edge_dofs(:, :), edge_orientations(:, :)
        integer, allocatable :: edges(:, :)
        integer, allocatable :: pivots(:)
        integer :: field_count, info

        status = 1
        call assemble_maxwell_fem_bem_system_3d( &
            vertices, tetrahedra, boundary_triangles, curl_coefficient, &
            mass_coefficient, wave_number, impedance, quadrature_degree, &
            tolerance, max_depth, matrix, status)
        if (status /= 0 .or. size(right_hand_side) /= size(matrix, 1)) return
        allocate(solution(size(matrix, 1), 1), pivots(size(matrix, 1)))
        solution(:, 1) = right_hand_side
        call zgesv( &
            size(matrix, 1), 1, matrix, size(matrix, 1), pivots, solution, &
            size(matrix, 1), info)
        if (info /= 0) then
            status = 2
            return
        end if
        call build_tetra_edge_dof_map( &
            tetrahedra, edges, edge_dofs, edge_orientations, status)
        if (status /= 0) return
        field_count = size(edges, 2)
        if (field_count < 1 .or. field_count >= size(matrix, 1)) return
        allocate(field(field_count), current(size(matrix, 1) - field_count))
        field = solution(:field_count, 1)
        current = solution(field_count + 1:, 1)
        status = 0
    end subroutine solve_maxwell_fem_bem_system_3d

    subroutine assemble_maxwell_fem_bem_system_3d( &
            vertices, tetrahedra, boundary_triangles, curl_coefficient, &
            mass_coefficient, wave_number, impedance, quadrature_degree, &
            tolerance, max_depth, matrix, status)
        real(dp), intent(in) :: vertices(:, :), curl_coefficient
        real(dp), intent(in) :: mass_coefficient, wave_number, impedance
        integer, intent(in) :: tetrahedra(:, :), boundary_triangles(:, :)
        integer, intent(in) :: quadrature_degree, max_depth
        real(dp), intent(in) :: tolerance
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: efie(:, :)
        integer, allocatable :: edge_dofs(:, :), edge_orientations(:, :)
        integer, allocatable :: edges(:, :), rwg_dofs(:)
        integer, allocatable :: rwg_edge_triangles(:, :), rwg_edges(:, :)
        real(dp), allocatable :: coupling(:, :), element_matrix(:, :)
        real(dp), allocatable :: rwg_mass(:, :), trace_scales(:)
        real(dp) :: element_vertices(3, 4)
        integer :: column, edge_count, local_column, local_row, node
        integer :: row, rwg_count, tetrahedron

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
        call assemble_maxwell_rwg_mass_matrix( &
            vertices, boundary_triangles, quadrature_degree, rwg_mass, status)
        if (status /= 0) return
        call assemble_maxwell_efie_rwg_3d( &
            vertices, boundary_triangles, wave_number, impedance, &
            quadrature_degree, tolerance, max_depth, efie, status)
        if (status /= 0) return
        edge_count = size(edges, 2)
        rwg_count = size(rwg_edges, 2)
        allocate(matrix(edge_count + rwg_count, edge_count + rwg_count))
        allocate(coupling(rwg_count, edge_count))
        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        coupling = 0.0_dp
        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                element_vertices(:, node) = &
                    vertices(:, tetrahedra(node, tetrahedron))
            end do
            call assemble_tetra_nedelec_curl_mass_element( &
                element_vertices, 1, quadrature_degree, element_matrix, status, &
                curl_coefficient=curl_coefficient, &
                mass_coefficient=mass_coefficient)
            if (status /= 0) return
            do local_column = 1, 6
                column = edge_dofs(local_column, tetrahedron)
                do local_row = 1, 6
                    row = edge_dofs(local_row, tetrahedron)
                    matrix(row, column) = matrix(row, column) + &
                        real(edge_orientations(local_row, tetrahedron)* &
                        edge_orientations(local_column, tetrahedron), dp)* &
                        element_matrix(local_row, local_column)
                end do
            end do
        end do
        do column = 1, rwg_count
            coupling(:, rwg_dofs(column)) = &
                coupling(:, rwg_dofs(column)) + &
                trace_scales(column)*rwg_mass(:, column)
        end do
        matrix(:edge_count, edge_count + 1:) = transpose(coupling)
        matrix(edge_count + 1:, :edge_count) = coupling
        matrix(edge_count + 1:, edge_count + 1:) = -efie
        status = 0
    end subroutine assemble_maxwell_fem_bem_system_3d

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
