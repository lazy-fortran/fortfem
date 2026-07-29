module fortfem_assembly_tetra_nedelec_3d
    use fortfem_kinds, only: dp
    use fortfem_tetra_edge_dof_map, only: build_tetra_edge_dof_map
    use fortfem_tetra_nedelec_first_order, only: &
        evaluate_tetra_nedelec_first_order
    use fortfem_tetra_piola_maps, only: map_tetra_nedelec_covariant
    use fortnum_linalg, only: det3
    use fortsparse, only: csc_from_triplet, csc_t, &
        FORTSPARSE_INVALID_MATRIX, fortsparse_status_t, status_set
    implicit none

    private

    public :: assemble_tetra_nedelec_curl_mass_csc

contains

    subroutine assemble_tetra_nedelec_curl_mass_csc( &
            mesh_vertices, tetrahedra, matrix, status, curl_coefficient, &
            mass_coefficient)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: curl_coefficient, mass_coefficient

        integer, allocatable :: columns(:), edges(:, :), global_dofs(:, :)
        integer, allocatable :: orientations(:, :), rows(:)
        real(dp), allocatable :: values(:)
        real(dp) :: curl_weight, element_matrix(6, 6), mass_weight
        real(dp) :: vertices(3, 4)
        integer :: column, entry, local_status, node, row, tetrahedron

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Tetrahedral Nedelec assembly failed")
        if (size(mesh_vertices, 1) /= 3) return
        if (size(tetrahedra, 1) /= 4) return
        if (size(tetrahedra, 2) < 1) return
        if (any(tetrahedra < 1)) return
        if (any(tetrahedra > size(mesh_vertices, 2))) return
        curl_weight = 1.0_dp
        mass_weight = 1.0_dp
        if (present(curl_coefficient)) curl_weight = curl_coefficient
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        call build_tetra_edge_dof_map( &
            tetrahedra, edges, global_dofs, orientations, local_status)
        if (local_status /= 0) return

        allocate(rows(36 * size(tetrahedra, 2)))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                vertices(:, node) = &
                    mesh_vertices(:, tetrahedra(node, tetrahedron))
            end do
            call assemble_tetra_nedelec_element( &
                vertices, curl_weight, mass_weight, element_matrix, &
                local_status)
            if (local_status /= 0) return
            do column = 1, 6
                do row = 1, 6
                    entry = entry + 1
                    rows(entry) = global_dofs(row, tetrahedron)
                    columns(entry) = global_dofs(column, tetrahedron)
                    values(entry) = real( &
                        orientations(row, tetrahedron) * &
                        orientations(column, tetrahedron), dp) * &
                        element_matrix(row, column)
                end do
            end do
        end do
        call csc_from_triplet( &
            size(edges, 2), size(edges, 2), rows, columns, values, &
            matrix, status)
    end subroutine assemble_tetra_nedelec_curl_mass_csc

    subroutine assemble_tetra_nedelec_element( &
            vertices, curl_coefficient, mass_coefficient, matrix, status)
        real(dp), intent(in) :: vertices(3, 4)
        real(dp), intent(in) :: curl_coefficient, mass_coefficient
        real(dp), intent(out) :: matrix(6, 6)
        integer, intent(out) :: status

        real(dp), parameter :: a = (5.0_dp + 3.0_dp * sqrt(5.0_dp)) / 20.0_dp
        real(dp), parameter :: b = (5.0_dp - sqrt(5.0_dp)) / 20.0_dp
        real(dp), parameter :: weight = 1.0_dp / 24.0_dp
        real(dp) :: determinant, jacobian(3, 3), points(3, 4)
        real(dp) :: physical_curls(3, 6), physical_values(3, 6)
        real(dp) :: reference_curls(3, 6), reference_values(3, 6)
        integer :: column, point, row

        status = 1
        matrix = 0.0_dp
        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
        determinant = det3(jacobian)
        if (determinant <= 64.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp, maxval(abs(jacobian))**3)) return
        points(:, 1) = [b, b, b]
        points(:, 2) = [a, b, b]
        points(:, 3) = [b, a, b]
        points(:, 4) = [b, b, a]
        do point = 1, 4
            call evaluate_tetra_nedelec_first_order( &
                points(:, point), reference_values, reference_curls, status)
            if (status /= 0) return
            call map_tetra_nedelec_covariant( &
                jacobian, reference_values, reference_curls, physical_values, &
                physical_curls, status)
            if (status /= 0) return
            do column = 1, 6
                do row = 1, 6
                    matrix(row, column) = matrix(row, column) + &
                        determinant * weight * ( &
                        curl_coefficient * dot_product( &
                        physical_curls(:, row), physical_curls(:, column)) + &
                        mass_coefficient * dot_product( &
                        physical_values(:, row), physical_values(:, column)))
                end do
            end do
        end do
        status = 0
    end subroutine assemble_tetra_nedelec_element

end module fortfem_assembly_tetra_nedelec_3d
