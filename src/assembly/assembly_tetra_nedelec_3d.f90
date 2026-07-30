module fortfem_assembly_tetra_nedelec_3d
    use fortfem_cartesian_helmholtz_pml, only: &
        cartesian_curl_curl_pml_coefficients
    use fortfem_kinds, only: dp
    use fortfem_tetra_edge_dof_map, only: build_tetra_edge_dof_map
    use fortfem_tetra_nedelec_first_order, only: &
        evaluate_tetra_nedelec_first_order
    use fortfem_tetra_nedelec_arbitrary_order, only: &
        evaluate_tetra_nedelec_first_kind, &
        initialize_tetra_nedelec_first_kind, tetra_nedelec_dof_count, &
        tetra_nedelec_first_kind_t
    use fortfem_tetra_nedelec_global_dof_map, only: &
        build_tetra_nedelec_basis_transform, build_tetra_nedelec_dof_map
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_tetra_piola_maps, only: map_tetra_nedelec_covariant
    use fortnum_linalg, only: det3
    use fortsparse, only: csc_from_triplet, csc_t, csc_z_t, &
        FORTSPARSE_INVALID_MATRIX, fortsparse_status_t, status_set
    implicit none

    private

    public :: assemble_tetra_nedelec_curl_mass_csc
    public :: assemble_tetra_nedelec_curl_mass_element
    public :: assemble_tetra_nedelec_pml_element
    public :: assemble_tetra_nedelec_pml_csc
    public :: assemble_tetra_nedelec_weighted_csc
    public :: assemble_tetra_nedelec_vector_load
    public :: assemble_tetra_nedelec_vector_load_order

    abstract interface
        pure subroutine tensor_coefficient_3d(x, y, z, value)
            import :: dp
            real(dp), intent(in) :: x, y, z
            real(dp), intent(out) :: value(3, 3)
        end subroutine tensor_coefficient_3d

        pure subroutine vector_source_3d(x, y, z, value)
            import :: dp
            real(dp), intent(in) :: x, y, z
            real(dp), intent(out) :: value(3)
        end subroutine vector_source_3d
    end interface

contains

    subroutine assemble_tetra_nedelec_vector_load_order( &
            mesh_vertices, tetrahedra, order, source, vector, status)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), order
        procedure(vector_source_3d) :: source
        real(dp), allocatable, intent(out) :: vector(:)
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: edge_orientations(:, :), edges(:, :)
        integer, allocatable :: face_permutations(:, :, :), faces(:, :)
        integer, allocatable :: global_dofs(:, :)
        real(dp), allocatable :: basis_transform(:, :), element_vector(:)
        real(dp), allocatable :: oriented_vector(:)
        real(dp) :: vertices(3, 4)
        integer :: dof, dof_count, local_status, node, tetrahedron

        if (order == 1) then
            call assemble_tetra_nedelec_vector_load( &
                mesh_vertices, tetrahedra, source, vector, status)
            return
        end if
        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Higher-order tetrahedral Nedelec vector load failed")
        if (.not. valid_tetra_mesh(mesh_vertices, tetrahedra)) return
        if (order < 2 .or. order > 4) return
        call build_tetra_nedelec_dof_map( &
            order, tetrahedra, edges, faces, global_dofs, &
            edge_orientations, face_permutations, local_status)
        if (local_status /= 0) return
        dof_count = size(global_dofs, 1)
        allocate(vector(maxval(global_dofs)))
        allocate(basis_transform(dof_count, dof_count))
        allocate(oriented_vector(dof_count))
        vector = 0.0_dp
        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                vertices(:, node) = &
                    mesh_vertices(:, tetrahedra(node, tetrahedron))
            end do
            call assemble_tetra_nedelec_load_element_order( &
                vertices, order, source, element_vector, local_status)
            if (local_status /= 0) return
            call build_tetra_nedelec_basis_transform( &
                order, edge_orientations(:, tetrahedron), &
                face_permutations(:, :, tetrahedron), basis_transform, &
                local_status)
            if (local_status /= 0) return
            oriented_vector = matmul(transpose(basis_transform), element_vector)
            do dof = 1, dof_count
                vector(global_dofs(dof, tetrahedron)) = &
                    vector(global_dofs(dof, tetrahedron)) + &
                    oriented_vector(dof)
            end do
        end do
        call status_set(status, 0, "")
    end subroutine assemble_tetra_nedelec_vector_load_order

    subroutine assemble_tetra_nedelec_load_element_order( &
            vertices, order, source, vector, status)
        real(dp), intent(in) :: vertices(3, 4)
        integer, intent(in) :: order
        procedure(vector_source_3d) :: source
        real(dp), allocatable, intent(out) :: vector(:)
        integer, intent(out) :: status

        type(tetra_nedelec_first_kind_t) :: basis
        real(dp), allocatable :: physical_curls(:, :), physical_values(:, :)
        real(dp), allocatable :: reference_curls(:, :)
        real(dp), allocatable :: reference_values(:, :)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: determinant, jacobian(3, 3), physical_point(3)
        real(dp) :: reference_point(3), source_value(3)
        integer :: dof, dof_count, point

        status = 1
        if (order < 2 .or. order > 4) return
        call initialize_tetra_nedelec_first_kind(order, basis, status)
        if (status /= 0) return
        call tetra_duffy_quadrature( &
            2 * order + 2, x, y, z, weights, status)
        if (status /= 0) return
        call tetra_geometry(vertices, jacobian, determinant, status)
        if (status /= 0) return
        dof_count = tetra_nedelec_dof_count(basis)
        allocate(vector(dof_count))
        allocate(reference_values(3, dof_count), reference_curls(3, dof_count))
        allocate(physical_values(3, dof_count), physical_curls(3, dof_count))
        vector = 0.0_dp
        do point = 1, size(weights)
            reference_point(1) = x(point)
            reference_point(2) = y(point)
            reference_point(3) = z(point)
            call evaluate_tetra_nedelec_first_kind( &
                basis, reference_point, reference_values, reference_curls, &
                status)
            if (status /= 0) return
            call map_tetra_nedelec_covariant( &
                jacobian, reference_values, reference_curls, physical_values, &
                physical_curls, status)
            if (status /= 0) return
            physical_point = vertices(:, 1)
            do dof = 1, 3
                physical_point = physical_point + &
                    jacobian(:, dof) * reference_point(dof)
            end do
            call source( &
                physical_point(1), physical_point(2), physical_point(3), &
                source_value)
            do dof = 1, dof_count
                vector(dof) = vector(dof) + determinant * weights(point) * &
                    dot_product(source_value, physical_values(:, dof))
            end do
        end do
        status = 0
    end subroutine assemble_tetra_nedelec_load_element_order

    subroutine assemble_tetra_nedelec_curl_mass_element( &
            vertices, order, quadrature_degree, matrix, status, &
            curl_coefficient, mass_coefficient)
        real(dp), intent(in) :: vertices(3, 4)
        integer, intent(in) :: order, quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: curl_coefficient, mass_coefficient

        type(tetra_nedelec_first_kind_t) :: basis
        real(dp), allocatable :: physical_curls(:, :), physical_values(:, :)
        real(dp), allocatable :: reference_curls(:, :)
        real(dp), allocatable :: reference_values(:, :)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: curl_weight, determinant, jacobian(3, 3), mass_weight
        real(dp) :: physical_weight, point(3)
        integer :: column, dof_count, point_index, row

        status = 1
        if (order < 1 .or. quadrature_degree < 0) return
        call initialize_tetra_nedelec_first_kind(order, basis, status)
        if (status /= 0) return
        call tetra_duffy_quadrature( &
            quadrature_degree, x, y, z, weights, status)
        if (status /= 0) return
        call tetra_geometry(vertices, jacobian, determinant, status)
        if (status /= 0) return

        curl_weight = 1.0_dp
        mass_weight = 1.0_dp
        if (present(curl_coefficient)) curl_weight = curl_coefficient
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        dof_count = tetra_nedelec_dof_count(basis)
        allocate(matrix(dof_count, dof_count))
        allocate( &
            reference_values(3, dof_count), reference_curls(3, dof_count), &
            physical_values(3, dof_count), physical_curls(3, dof_count))
        matrix = 0.0_dp
        do point_index = 1, size(weights)
            point = [x(point_index), y(point_index), z(point_index)]
            call evaluate_tetra_nedelec_first_kind( &
                basis, point, reference_values, reference_curls, status)
            if (status /= 0) return
            call map_tetra_nedelec_covariant( &
                jacobian, reference_values, reference_curls, physical_values, &
                physical_curls, status)
            if (status /= 0) return
            physical_weight = determinant * weights(point_index)
            do column = 1, dof_count
                do row = 1, dof_count
                    matrix(row, column) = matrix(row, column) + &
                        physical_weight * (curl_weight * dot_product( &
                        physical_curls(:, row), physical_curls(:, column)) + &
                        mass_weight * dot_product( &
                        physical_values(:, row), physical_values(:, column)))
                end do
            end do
        end do
        status = 0
    end subroutine assemble_tetra_nedelec_curl_mass_element

    subroutine assemble_tetra_nedelec_pml_element( &
            vertices, order, quadrature_degree, stretch, wave_number, matrix, &
            status)
        real(dp), intent(in) :: vertices(3, 4)
        integer, intent(in) :: order, quadrature_degree
        complex(dp), intent(in) :: stretch(3)
        real(dp), intent(in) :: wave_number
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        type(tetra_nedelec_first_kind_t) :: basis
        complex(dp) :: curl_coefficient(3), mass_coefficient(3)
        real(dp), allocatable :: physical_curls(:, :), physical_values(:, :)
        real(dp), allocatable :: reference_curls(:, :)
        real(dp), allocatable :: reference_values(:, :)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: determinant, jacobian(3, 3), physical_weight, point(3)
        integer :: column, dof_count, point_index, row

        status = 1
        if (allocated(matrix)) deallocate(matrix)
        if (order < 1 .or. order > 4) return
        if (quadrature_degree < 0 .or. wave_number <= 0.0_dp) return
        call cartesian_curl_curl_pml_coefficients( &
            stretch, curl_coefficient, mass_coefficient, status)
        if (status /= 0) return
        call initialize_tetra_nedelec_first_kind(order, basis, status)
        if (status /= 0) return
        call tetra_duffy_quadrature( &
            quadrature_degree, x, y, z, weights, status)
        if (status /= 0) return
        call tetra_geometry(vertices, jacobian, determinant, status)
        if (status /= 0) return

        dof_count = tetra_nedelec_dof_count(basis)
        allocate(matrix(dof_count, dof_count))
        allocate( &
            reference_values(3, dof_count), reference_curls(3, dof_count), &
            physical_values(3, dof_count), physical_curls(3, dof_count))
        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        do point_index = 1, size(weights)
            point = [x(point_index), y(point_index), z(point_index)]
            call evaluate_tetra_nedelec_first_kind( &
                basis, point, reference_values, reference_curls, status)
            if (status /= 0) return
            call map_tetra_nedelec_covariant( &
                jacobian, reference_values, reference_curls, physical_values, &
                physical_curls, status)
            if (status /= 0) return
            physical_weight = determinant*weights(point_index)
            do column = 1, dof_count
                do row = 1, dof_count
                    matrix(row, column) = matrix(row, column) + &
                        physical_weight*(sum(curl_coefficient* &
                        physical_curls(:, row)*physical_curls(:, column)) - &
                        wave_number**2*sum(mass_coefficient* &
                        physical_values(:, row)*physical_values(:, column)))
                end do
            end do
        end do
        status = 0
    end subroutine assemble_tetra_nedelec_pml_element

    subroutine assemble_tetra_nedelec_pml_csc( &
            mesh_vertices, tetrahedra, order, stretch, wave_number, matrix, &
            status)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), order
        complex(dp), intent(in) :: stretch(:, :)
        real(dp), intent(in) :: wave_number
        type(csc_z_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), edge_orientations(:, :)
        integer, allocatable :: edges(:, :), face_permutations(:, :, :)
        integer, allocatable :: faces(:, :), global_dofs(:, :), rows(:)
        complex(dp), allocatable :: element_matrix(:, :)
        complex(dp), allocatable :: oriented_matrix(:, :), values(:)
        real(dp), allocatable :: basis_transform(:, :)
        real(dp) :: vertices(3, 4)
        integer :: column, dof_count, entry, global_dof_count
        integer :: local_status, node, row, tetrahedron

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Complex tetrahedral Nedelec PML assembly failed")
        if (.not. valid_tetra_mesh(mesh_vertices, tetrahedra)) return
        if (order < 1 .or. order > 4) return
        if (size(stretch, 1) /= 3) return
        if (size(stretch, 2) /= size(tetrahedra, 2)) return
        if (wave_number <= 0.0_dp) return
        call build_tetra_nedelec_dof_map( &
            order, tetrahedra, edges, faces, global_dofs, &
            edge_orientations, face_permutations, local_status)
        if (local_status /= 0) return
        dof_count = size(global_dofs, 1)
        global_dof_count = maxval(global_dofs)
        allocate( &
            rows(dof_count*dof_count*size(tetrahedra, 2)), &
            columns(dof_count*dof_count*size(tetrahedra, 2)), &
            values(dof_count*dof_count*size(tetrahedra, 2)))
        allocate( &
            basis_transform(dof_count, dof_count), &
            oriented_matrix(dof_count, dof_count))

        entry = 0
        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                vertices(:, node) = &
                    mesh_vertices(:, tetrahedra(node, tetrahedron))
            end do
            call assemble_tetra_nedelec_pml_element( &
                vertices, order, 2*order + 2, stretch(:, tetrahedron), &
                wave_number, element_matrix, local_status)
            if (local_status /= 0) return
            call build_tetra_nedelec_basis_transform( &
                order, edge_orientations(:, tetrahedron), &
                face_permutations(:, :, tetrahedron), basis_transform, &
                local_status)
            if (local_status /= 0) return
            oriented_matrix = matmul( &
                transpose(basis_transform), &
                matmul(element_matrix, basis_transform))
            do column = 1, dof_count
                do row = 1, dof_count
                    entry = entry + 1
                    rows(entry) = global_dofs(row, tetrahedron)
                    columns(entry) = global_dofs(column, tetrahedron)
                    values(entry) = oriented_matrix(row, column)
                end do
            end do
        end do
        call csc_from_triplet( &
            global_dof_count, global_dof_count, rows, columns, values, &
            matrix, status)
    end subroutine assemble_tetra_nedelec_pml_csc

    subroutine assemble_tetra_nedelec_weighted_csc( &
            mesh_vertices, tetrahedra, coefficient, mass_coefficient, &
            matrix, status)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :)
        procedure(tensor_coefficient_3d) :: coefficient
        real(dp), intent(in) :: mass_coefficient
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), edges(:, :), global_dofs(:, :)
        integer, allocatable :: orientations(:, :), rows(:)
        real(dp), allocatable :: triplet_values(:)
        real(dp) :: element_matrix(6, 6), vertices(3, 4)
        integer :: column, entry, local_status, node, row, tetrahedron

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Weighted tetrahedral Nedelec assembly failed")
        if (.not. valid_tetra_mesh(mesh_vertices, tetrahedra)) return
        call build_tetra_edge_dof_map( &
            tetrahedra, edges, global_dofs, orientations, local_status)
        if (local_status /= 0) return
        allocate(rows(36 * size(tetrahedra, 2)))
        allocate(columns(size(rows)), triplet_values(size(rows)))
        entry = 0
        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                vertices(:, node) = &
                    mesh_vertices(:, tetrahedra(node, tetrahedron))
            end do
            call assemble_tetra_nedelec_weighted_element( &
                vertices, coefficient, mass_coefficient, element_matrix, &
                local_status)
            if (local_status /= 0) return
            do column = 1, 6
                do row = 1, 6
                    entry = entry + 1
                    rows(entry) = global_dofs(row, tetrahedron)
                    columns(entry) = global_dofs(column, tetrahedron)
                    triplet_values(entry) = real( &
                        orientations(row, tetrahedron) * &
                        orientations(column, tetrahedron), dp) * &
                        element_matrix(row, column)
                end do
            end do
        end do
        call csc_from_triplet( &
            size(edges, 2), size(edges, 2), rows, columns, triplet_values, &
            matrix, status)
    end subroutine assemble_tetra_nedelec_weighted_csc

    subroutine assemble_tetra_nedelec_vector_load( &
            mesh_vertices, tetrahedra, source, vector, status)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :)
        procedure(vector_source_3d) :: source
        real(dp), allocatable, intent(out) :: vector(:)
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: edges(:, :), global_dofs(:, :)
        integer, allocatable :: orientations(:, :)
        real(dp) :: element_vector(6), vertices(3, 4)
        integer :: dof, local_status, node, tetrahedron

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Tetrahedral Nedelec vector load assembly failed")
        if (.not. valid_tetra_mesh(mesh_vertices, tetrahedra)) return
        call build_tetra_edge_dof_map( &
            tetrahedra, edges, global_dofs, orientations, local_status)
        if (local_status /= 0) return
        allocate(vector(size(edges, 2)))
        vector = 0.0_dp
        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                vertices(:, node) = &
                    mesh_vertices(:, tetrahedra(node, tetrahedron))
            end do
            call assemble_tetra_nedelec_load_element( &
                vertices, source, element_vector, local_status)
            if (local_status /= 0) return
            do dof = 1, 6
                vector(global_dofs(dof, tetrahedron)) = &
                    vector(global_dofs(dof, tetrahedron)) + &
                    real(orientations(dof, tetrahedron), dp) * &
                    element_vector(dof)
            end do
        end do
        call status_set(status, 0, "")
    end subroutine assemble_tetra_nedelec_vector_load

    subroutine assemble_tetra_nedelec_curl_mass_csc( &
            mesh_vertices, tetrahedra, matrix, status, curl_coefficient, &
            mass_coefficient, order)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: curl_coefficient, mass_coefficient
        integer, intent(in), optional :: order

        integer, allocatable :: columns(:), edges(:, :), global_dofs(:, :)
        integer, allocatable :: orientations(:, :), rows(:)
        real(dp), allocatable :: values(:)
        real(dp) :: curl_weight, element_matrix(6, 6), mass_weight
        real(dp) :: vertices(3, 4)
        integer :: column, entry, local_status, node, polynomial_order
        integer :: row, tetrahedron

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Tetrahedral Nedelec assembly failed")
        if (size(mesh_vertices, 1) /= 3) return
        if (size(tetrahedra, 1) /= 4) return
        if (size(tetrahedra, 2) < 1) return
        if (any(tetrahedra < 1)) return
        if (any(tetrahedra > size(mesh_vertices, 2))) return
        polynomial_order = 1
        if (present(order)) polynomial_order = order
        if (polynomial_order < 1 .or. polynomial_order > 4) return
        curl_weight = 1.0_dp
        mass_weight = 1.0_dp
        if (present(curl_coefficient)) curl_weight = curl_coefficient
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        if (polynomial_order /= 1) then
            call assemble_arbitrary_order_curl_mass_csc( &
                mesh_vertices, tetrahedra, polynomial_order, curl_weight, &
                mass_weight, matrix, status)
            return
        end if
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

    subroutine assemble_arbitrary_order_curl_mass_csc( &
            mesh_vertices, tetrahedra, order, curl_coefficient, &
            mass_coefficient, matrix, status)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), order
        real(dp), intent(in) :: curl_coefficient, mass_coefficient
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), edge_orientations(:, :)
        integer, allocatable :: edges(:, :), face_permutations(:, :, :)
        integer, allocatable :: faces(:, :), global_dofs(:, :), rows(:)
        real(dp), allocatable :: basis_transform(:, :), element_matrix(:, :)
        real(dp), allocatable :: oriented_matrix(:, :), values(:)
        real(dp) :: vertices(3, 4)
        integer :: column, dof_count, entry, global_dof_count
        integer :: local_status, node, row, tetrahedron

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Higher-order tetrahedral Nedelec assembly failed")
        if (order < 2 .or. order > 4) return
        call build_tetra_nedelec_dof_map( &
            order, tetrahedra, edges, faces, global_dofs, &
            edge_orientations, face_permutations, local_status)
        if (local_status /= 0) return
        dof_count = size(global_dofs, 1)
        global_dof_count = maxval(global_dofs)
        allocate( &
            rows(dof_count * dof_count * size(tetrahedra, 2)), &
            columns(dof_count * dof_count * size(tetrahedra, 2)), &
            values(dof_count * dof_count * size(tetrahedra, 2)))
        allocate( &
            basis_transform(dof_count, dof_count), &
            oriented_matrix(dof_count, dof_count))

        entry = 0
        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                vertices(:, node) = &
                    mesh_vertices(:, tetrahedra(node, tetrahedron))
            end do
            call assemble_tetra_nedelec_curl_mass_element( &
                vertices, order, 2 * order, element_matrix, local_status, &
                curl_coefficient, mass_coefficient)
            if (local_status /= 0) return
            call build_tetra_nedelec_basis_transform( &
                order, edge_orientations(:, tetrahedron), &
                face_permutations(:, :, tetrahedron), basis_transform, &
                local_status)
            if (local_status /= 0) return
            oriented_matrix = matmul( &
                transpose(basis_transform), &
                matmul(element_matrix, basis_transform))
            do column = 1, dof_count
                do row = 1, dof_count
                    entry = entry + 1
                    rows(entry) = global_dofs(row, tetrahedron)
                    columns(entry) = global_dofs(column, tetrahedron)
                    values(entry) = oriented_matrix(row, column)
                end do
            end do
        end do
        call csc_from_triplet( &
            global_dof_count, global_dof_count, rows, columns, values, &
            matrix, status)
    end subroutine assemble_arbitrary_order_curl_mass_csc

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

    subroutine assemble_tetra_nedelec_weighted_element( &
            vertices, coefficient, mass_coefficient, matrix, status)
        real(dp), intent(in) :: vertices(3, 4)
        procedure(tensor_coefficient_3d) :: coefficient
        real(dp), intent(in) :: mass_coefficient
        real(dp), intent(out) :: matrix(6, 6)
        integer, intent(out) :: status

        real(dp) :: determinant, jacobian(3, 3), physical_point(3)
        real(dp) :: physical_curls(3, 6), physical_values(3, 6)
        real(dp) :: points(3, 4), reference_curls(3, 6)
        real(dp) :: reference_values(3, 6), tensor(3, 3), weights(4)
        integer :: column, point, row

        status = 1
        matrix = 0.0_dp
        call tetra_geometry(vertices, jacobian, determinant, status)
        if (status /= 0) return
        call tetra_degree_two_quadrature(points, weights)
        do point = 1, 4
            call evaluate_tetra_nedelec_first_order( &
                points(:, point), reference_values, reference_curls, status)
            if (status /= 0) return
            call map_tetra_nedelec_covariant( &
                jacobian, reference_values, reference_curls, physical_values, &
                physical_curls, status)
            if (status /= 0) return
            physical_point = vertices(:, 1) + &
                matmul(jacobian, points(:, point))
            call coefficient( &
                physical_point(1), physical_point(2), physical_point(3), &
                tensor)
            do column = 1, 6
                do row = 1, 6
                    matrix(row, column) = matrix(row, column) + &
                        determinant * weights(point) * (dot_product( &
                        physical_curls(:, row), &
                        matmul(tensor, physical_curls(:, column))) + &
                        mass_coefficient * dot_product( &
                        physical_values(:, row), physical_values(:, column)))
                end do
            end do
        end do
        status = 0
    end subroutine assemble_tetra_nedelec_weighted_element

    subroutine assemble_tetra_nedelec_load_element( &
            vertices, source, vector, status)
        real(dp), intent(in) :: vertices(3, 4)
        procedure(vector_source_3d) :: source
        real(dp), intent(out) :: vector(6)
        integer, intent(out) :: status

        real(dp) :: determinant, jacobian(3, 3), physical_point(3)
        real(dp) :: physical_curls(3, 6), physical_values(3, 6)
        real(dp) :: points(3, 4), reference_curls(3, 6)
        real(dp) :: reference_values(3, 6), source_value(3), weights(4)
        integer :: dof, point

        status = 1
        vector = 0.0_dp
        call tetra_geometry(vertices, jacobian, determinant, status)
        if (status /= 0) return
        call tetra_degree_two_quadrature(points, weights)
        do point = 1, 4
            call evaluate_tetra_nedelec_first_order( &
                points(:, point), reference_values, reference_curls, status)
            if (status /= 0) return
            call map_tetra_nedelec_covariant( &
                jacobian, reference_values, reference_curls, physical_values, &
                physical_curls, status)
            if (status /= 0) return
            physical_point = vertices(:, 1) + &
                matmul(jacobian, points(:, point))
            call source( &
                physical_point(1), physical_point(2), physical_point(3), &
                source_value)
            do dof = 1, 6
                vector(dof) = vector(dof) + determinant * weights(point) * &
                    dot_product(source_value, physical_values(:, dof))
            end do
        end do
        status = 0
    end subroutine assemble_tetra_nedelec_load_element

    pure subroutine tetra_geometry(vertices, jacobian, determinant, status)
        real(dp), intent(in) :: vertices(3, 4)
        real(dp), intent(out) :: jacobian(3, 3), determinant
        integer, intent(out) :: status

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
        determinant = det3(jacobian)
        status = 1
        if (determinant <= 64.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp, maxval(abs(jacobian))**3)) return
        status = 0
    end subroutine tetra_geometry

    pure subroutine tetra_degree_two_quadrature(points, weights)
        real(dp), intent(out) :: points(3, 4), weights(4)

        real(dp), parameter :: a = (5.0_dp + 3.0_dp * sqrt(5.0_dp)) / 20.0_dp
        real(dp), parameter :: b = (5.0_dp - sqrt(5.0_dp)) / 20.0_dp

        points(:, 1) = [b, b, b]
        points(:, 2) = [a, b, b]
        points(:, 3) = [b, a, b]
        points(:, 4) = [b, b, a]
        weights = 1.0_dp / 24.0_dp
    end subroutine tetra_degree_two_quadrature

    pure logical function valid_tetra_mesh(mesh_vertices, tetrahedra)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :)

        valid_tetra_mesh = .false.
        if (size(mesh_vertices, 1) /= 3) return
        if (size(tetrahedra, 1) /= 4) return
        if (size(tetrahedra, 2) < 1) return
        if (any(tetrahedra < 1)) return
        if (any(tetrahedra > size(mesh_vertices, 2))) return
        valid_tetra_mesh = .true.
    end function valid_tetra_mesh

end module fortfem_assembly_tetra_nedelec_3d
