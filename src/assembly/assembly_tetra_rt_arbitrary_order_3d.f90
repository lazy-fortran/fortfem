module fortfem_assembly_tetra_rt_arbitrary_order_3d
    use fortfem_kinds, only: dp
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_tetra_discontinuous_arbitrary_order, only: &
        evaluate_tetra_discontinuous, initialize_tetra_discontinuous, &
        tetra_discontinuous_dof_count, tetra_discontinuous_t
    use fortfem_tetra_piola_maps, only: map_tetra_rt_contravariant, &
        map_tetra_rt_contravariant_jvp, map_tetra_rt_contravariant_vjp
    use fortfem_tetra_rt_arbitrary_order, only: &
        evaluate_tetra_rt, initialize_tetra_rt, tetra_rt_dof_count, tetra_rt_t
    use fortfem_tetra_rt_global_dof_map, only: &
        build_tetra_rt_basis_transform, build_tetra_rt_dof_map
    use fortnum_linalg, only: det3, det3_jvp, det3_vjp
    use fortsparse, only: &
        csc_from_triplet, csc_t, FORTSPARSE_INVALID_MATRIX, &
        fortsparse_status_t, status_set
    implicit none

    private

    public :: assemble_tetra_rt_div_mass_csc
    public :: assemble_tetra_rt_div_mass_csc_jvp
    public :: assemble_tetra_rt_div_mass_csc_vjp
    public :: assemble_tetra_rt_div_mass_element
    public :: assemble_tetra_rt_div_mass_element_jvp
    public :: assemble_tetra_rt_div_mass_element_vjp
    public :: assemble_tetra_rt_divergence_csc
    public :: assemble_tetra_rt_vector_load

    abstract interface
        pure subroutine vector_source_3d(x, y, z, value)
            import :: dp
            real(dp), intent(in) :: x, y, z
            real(dp), intent(out) :: value(3)
        end subroutine vector_source_3d
    end interface

contains

    subroutine assemble_tetra_rt_vector_load( &
            mesh_vertices, tetrahedra, degree, quadrature_degree, source, &
            vector, status)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        procedure(vector_source_3d) :: source
        real(dp), allocatable, intent(out) :: vector(:)
        type(fortsparse_status_t), intent(out) :: status

        type(tetra_rt_t) :: basis
        integer, allocatable :: face_orientations(:, :)
        integer, allocatable :: face_permutations(:, :, :), faces(:, :)
        integer, allocatable :: global_dofs(:, :)
        real(dp), allocatable :: basis_transform(:, :), element_vector(:)
        real(dp), allocatable :: oriented_vector(:)
        real(dp), allocatable :: physical_divergences(:)
        real(dp), allocatable :: physical_values(:, :)
        real(dp), allocatable :: reference_divergences(:)
        real(dp), allocatable :: reference_values(:, :)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: determinant, jacobian(3, 3), physical_point(3)
        real(dp) :: point(3), source_value(3), vertices(3, 4)
        integer :: component, dof, dof_count, local_status, node, point_index
        integer :: tetrahedron

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Tetrahedral RT vector load failed")
        if (.not. valid_tetra_mesh(mesh_vertices, tetrahedra)) return
        if (degree < 0 .or. quadrature_degree < 0) return
        call initialize_tetra_rt(degree, basis, local_status)
        if (local_status /= 0) return
        call tetra_duffy_quadrature( &
            quadrature_degree, x, y, z, weights, local_status)
        if (local_status /= 0) return
        call build_tetra_rt_dof_map( &
            degree, tetrahedra, faces, global_dofs, face_orientations, &
            face_permutations, local_status)
        if (local_status /= 0) return

        dof_count = tetra_rt_dof_count(basis)
        allocate(vector(maxval(global_dofs)))
        allocate(element_vector(dof_count))
        allocate(oriented_vector(dof_count))
        allocate(basis_transform(dof_count, dof_count))
        allocate(reference_values(3, dof_count))
        allocate(reference_divergences(dof_count))
        allocate(physical_values(3, dof_count))
        allocate(physical_divergences(dof_count))
        vector = 0.0_dp
        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                vertices(:, node) = &
                    mesh_vertices(:, tetrahedra(node, tetrahedron))
            end do
            call tetra_geometry(vertices, jacobian, determinant, local_status)
            if (local_status /= 0) return
            element_vector = 0.0_dp
            do point_index = 1, size(weights)
                point(1) = x(point_index)
                point(2) = y(point_index)
                point(3) = z(point_index)
                call evaluate_tetra_rt( &
                    basis, point, reference_values, reference_divergences, &
                    local_status)
                if (local_status /= 0) return
                call map_tetra_rt_contravariant( &
                    jacobian, reference_values, reference_divergences, &
                    physical_values, physical_divergences, local_status)
                if (local_status /= 0) return
                physical_point = vertices(:, 1)
                do component = 1, 3
                    physical_point = physical_point + &
                        jacobian(:, component)*point(component)
                end do
                call source( &
                    physical_point(1), physical_point(2), physical_point(3), &
                    source_value)
                do dof = 1, dof_count
                    element_vector(dof) = element_vector(dof) + &
                        determinant*weights(point_index)*dot_product( &
                        source_value, physical_values(:, dof))
                end do
            end do
            call build_tetra_rt_basis_transform( &
                degree, face_orientations(:, tetrahedron), &
                face_permutations(:, :, tetrahedron), basis_transform, &
                local_status)
            if (local_status /= 0) return
            oriented_vector = &
                matmul(transpose(basis_transform), element_vector)
            do dof = 1, dof_count
                vector(global_dofs(dof, tetrahedron)) = &
                    vector(global_dofs(dof, tetrahedron)) + &
                    oriented_vector(dof)
            end do
        end do
        call status_set(status, 0, "")
    end subroutine assemble_tetra_rt_vector_load

    subroutine assemble_tetra_rt_divergence_csc( &
            mesh_vertices, tetrahedra, degree, quadrature_degree, matrix, &
            status)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :)
        integer, intent(in) :: degree, quadrature_degree
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        type(tetra_discontinuous_t) :: scalar_basis
        type(tetra_rt_t) :: rt_basis
        integer, allocatable :: columns(:), face_orientations(:, :)
        integer, allocatable :: face_permutations(:, :, :), faces(:, :)
        integer, allocatable :: global_dofs(:, :), rows(:)
        real(dp), allocatable :: basis_transform(:, :)
        real(dp), allocatable :: local_matrix(:, :), oriented_matrix(:, :)
        real(dp), allocatable :: physical_divergences(:)
        real(dp), allocatable :: physical_values(:, :)
        real(dp), allocatable :: reference_divergences(:)
        real(dp), allocatable :: reference_values(:, :)
        real(dp), allocatable :: scalar_values(:), values(:)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: determinant, jacobian(3, 3), point(3), vertices(3, 4)
        integer :: entry, local_status, node, point_index, rt_dof
        integer :: rt_dof_count, scalar_dof, scalar_dof_count
        integer :: scalar_global_count, tetrahedron

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Tetrahedral RT-DG divergence assembly failed")
        if (.not. valid_tetra_mesh(mesh_vertices, tetrahedra)) return
        if (degree < 0) return
        if (quadrature_degree < 0) return
        call initialize_tetra_rt(degree, rt_basis, local_status)
        if (local_status /= 0) return
        call initialize_tetra_discontinuous( &
            degree, scalar_basis, local_status)
        if (local_status /= 0) return
        call tetra_duffy_quadrature( &
            quadrature_degree, x, y, z, weights, local_status)
        if (local_status /= 0) return
        call build_tetra_rt_dof_map( &
            degree, tetrahedra, faces, global_dofs, face_orientations, &
            face_permutations, local_status)
        if (local_status /= 0) return

        rt_dof_count = tetra_rt_dof_count(rt_basis)
        scalar_dof_count = tetra_discontinuous_dof_count(scalar_basis)
        scalar_global_count = scalar_dof_count*size(tetrahedra, 2)
        allocate( &
            reference_values(3, rt_dof_count), &
            reference_divergences(rt_dof_count), &
            physical_values(3, rt_dof_count), &
            physical_divergences(rt_dof_count), &
            scalar_values(scalar_dof_count), &
            local_matrix(scalar_dof_count, rt_dof_count), &
            oriented_matrix(scalar_dof_count, rt_dof_count), &
            basis_transform(rt_dof_count, rt_dof_count))
        allocate(rows(scalar_global_count*rt_dof_count))
        allocate(columns(size(rows)), values(size(rows)))

        entry = 0
        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                vertices(:, node) = &
                    mesh_vertices(:, tetrahedra(node, tetrahedron))
            end do
            call tetra_geometry(vertices, jacobian, determinant, local_status)
            if (local_status /= 0) return
            local_matrix = 0.0_dp
            do point_index = 1, size(weights)
                point = [x(point_index), y(point_index), z(point_index)]
                call evaluate_tetra_rt( &
                    rt_basis, point, reference_values, &
                    reference_divergences, local_status)
                if (local_status /= 0) return
                call map_tetra_rt_contravariant( &
                    jacobian, reference_values, reference_divergences, &
                    physical_values, physical_divergences, local_status)
                if (local_status /= 0) return
                call evaluate_tetra_discontinuous( &
                    scalar_basis, point(1), point(2), point(3), &
                    scalar_values, local_status)
                if (local_status /= 0) return
                do rt_dof = 1, rt_dof_count
                    do scalar_dof = 1, scalar_dof_count
                        local_matrix(scalar_dof, rt_dof) = &
                            local_matrix(scalar_dof, rt_dof) + &
                            determinant*weights(point_index)* &
                            scalar_values(scalar_dof)* &
                            physical_divergences(rt_dof)
                    end do
                end do
            end do
            call build_tetra_rt_basis_transform( &
                degree, face_orientations(:, tetrahedron), &
                face_permutations(:, :, tetrahedron), basis_transform, &
                local_status)
            if (local_status /= 0) return
            oriented_matrix = matmul(local_matrix, basis_transform)
            do rt_dof = 1, rt_dof_count
                do scalar_dof = 1, scalar_dof_count
                    entry = entry + 1
                    rows(entry) = &
                        (tetrahedron - 1)*scalar_dof_count + scalar_dof
                    columns(entry) = global_dofs(rt_dof, tetrahedron)
                    values(entry) = oriented_matrix(scalar_dof, rt_dof)
                end do
            end do
        end do
        call csc_from_triplet( &
            scalar_global_count, maxval(global_dofs), rows, columns, values, &
            matrix, status)
    end subroutine assemble_tetra_rt_divergence_csc

    subroutine assemble_tetra_rt_div_mass_element( &
            vertices, degree, quadrature_degree, matrix, status, &
            divergence_coefficient, mass_coefficient)
        real(dp), intent(in) :: vertices(3, 4)
        integer, intent(in) :: degree, quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: divergence_coefficient
        real(dp), intent(in), optional :: mass_coefficient

        type(tetra_rt_t) :: basis
        real(dp), allocatable :: physical_divergences(:)
        real(dp), allocatable :: physical_values(:, :)
        real(dp), allocatable :: reference_divergences(:)
        real(dp), allocatable :: reference_values(:, :)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: determinant, divergence_weight, jacobian(3, 3)
        real(dp) :: mass_weight, physical_weight, point(3)
        integer :: column, dof_count, point_index, row

        status = 1
        if (degree < 0) return
        if (quadrature_degree < 0) return
        call tetra_geometry(vertices, jacobian, determinant, status)
        if (status /= 0) return
        call initialize_tetra_rt(degree, basis, status)
        if (status /= 0) return
        call tetra_duffy_quadrature( &
            quadrature_degree, x, y, z, weights, status)
        if (status /= 0) return

        divergence_weight = 1.0_dp
        mass_weight = 1.0_dp
        if (present(divergence_coefficient)) then
            divergence_weight = divergence_coefficient
        end if
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        dof_count = tetra_rt_dof_count(basis)
        allocate(matrix(dof_count, dof_count))
        allocate( &
            reference_values(3, dof_count), &
            reference_divergences(dof_count), &
            physical_values(3, dof_count), &
            physical_divergences(dof_count))
        matrix = 0.0_dp
        do point_index = 1, size(weights)
            point = [x(point_index), y(point_index), z(point_index)]
            call evaluate_tetra_rt( &
                basis, point, reference_values, reference_divergences, status)
            if (status /= 0) return
            call map_tetra_rt_contravariant( &
                jacobian, reference_values, reference_divergences, &
                physical_values, physical_divergences, status)
            if (status /= 0) return
            physical_weight = determinant*weights(point_index)
            do column = 1, dof_count
                do row = 1, dof_count
                    matrix(row, column) = matrix(row, column) + &
                        physical_weight*(divergence_weight* &
                        physical_divergences(row)* &
                        physical_divergences(column) + mass_weight* &
                        dot_product( &
                        physical_values(:, row), physical_values(:, column)))
                end do
            end do
        end do
        status = 0
    end subroutine assemble_tetra_rt_div_mass_element

    subroutine assemble_tetra_rt_div_mass_element_jvp( &
            vertices, degree, quadrature_degree, divergence_coefficient, &
            mass_coefficient, vertices_dot, divergence_coefficient_dot, &
            mass_coefficient_dot, matrix_dot, status)
        real(dp), intent(in) :: vertices(3, 4), vertices_dot(3, 4)
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: divergence_coefficient, mass_coefficient
        real(dp), intent(in) :: divergence_coefficient_dot
        real(dp), intent(in) :: mass_coefficient_dot
        real(dp), allocatable, intent(out) :: matrix_dot(:, :)
        integer, intent(out) :: status

        type(tetra_rt_t) :: basis
        real(dp), allocatable :: divergences(:), divergences_dot(:)
        real(dp), allocatable :: ref_divergences(:), ref_values(:, :)
        real(dp), allocatable :: values(:, :), values_dot(:, :)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp), allocatable :: zero_divergences(:), zero_values(:, :)
        real(dp) :: determinant, determinant_dot, divergence_energy
        real(dp) :: jacobian(3, 3), jacobian_dot(3, 3), mass_energy, point(3)
        integer :: column, dof_count, point_index, row

        status = 1
        if (allocated(matrix_dot)) deallocate(matrix_dot)
        if (degree < 0 .or. quadrature_degree < 0) return
        call initialize_tetra_rt(degree, basis, status)
        if (status /= 0) return
        call tetra_duffy_quadrature( &
            quadrature_degree, x, y, z, weights, status)
        if (status /= 0) return
        call tetra_geometry(vertices, jacobian, determinant, status)
        if (status /= 0) return
        call tetra_jacobian(vertices_dot, jacobian_dot)
        call det3_jvp(jacobian, jacobian_dot, determinant_dot)
        dof_count = tetra_rt_dof_count(basis)
        allocate(matrix_dot(dof_count, dof_count), source=0.0_dp)
        allocate(ref_values(3, dof_count), ref_divergences(dof_count))
        allocate(zero_values(3, dof_count), zero_divergences(dof_count), &
            source=0.0_dp)
        allocate(values(3, dof_count), divergences(dof_count))
        allocate(values_dot(3, dof_count), divergences_dot(dof_count))
        do point_index = 1, size(weights)
            point = [x(point_index), y(point_index), z(point_index)]
            call evaluate_tetra_rt( &
                basis, point, ref_values, ref_divergences, status)
            if (status /= 0) return
            call map_tetra_rt_contravariant( &
                jacobian, ref_values, ref_divergences, values, divergences, &
                status)
            if (status /= 0) return
            call map_tetra_rt_contravariant_jvp( &
                jacobian, ref_values, ref_divergences, jacobian_dot, &
                zero_values, zero_divergences, values_dot, divergences_dot, &
                status)
            if (status /= 0) return
            do column = 1, dof_count
                do row = 1, dof_count
                    divergence_energy = divergences(row)*divergences(column)
                    mass_energy = dot_product(values(:, row), values(:, column))
                    matrix_dot(row, column) = matrix_dot(row, column) + &
                        weights(point_index)*(determinant_dot*( &
                        divergence_coefficient*divergence_energy + &
                        mass_coefficient*mass_energy) + determinant*( &
                        divergence_coefficient_dot*divergence_energy + &
                        divergence_coefficient*( &
                        divergences_dot(row)*divergences(column) + &
                        divergences(row)*divergences_dot(column)) + &
                        mass_coefficient_dot*mass_energy + mass_coefficient*( &
                        dot_product(values_dot(:, row), values(:, column)) + &
                        dot_product(values(:, row), values_dot(:, column)))))
                end do
            end do
        end do
        status = 0
    end subroutine assemble_tetra_rt_div_mass_element_jvp

    subroutine assemble_tetra_rt_div_mass_element_vjp( &
            vertices, degree, quadrature_degree, divergence_coefficient, &
            mass_coefficient, matrix_bar, vertices_bar, &
            divergence_coefficient_bar, mass_coefficient_bar, status)
        real(dp), intent(in) :: vertices(3, 4)
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: divergence_coefficient, mass_coefficient
        real(dp), intent(in) :: matrix_bar(:, :)
        real(dp), intent(out) :: vertices_bar(3, 4)
        real(dp), intent(out) :: divergence_coefficient_bar
        real(dp), intent(out) :: mass_coefficient_bar
        integer, intent(out) :: status

        type(tetra_rt_t) :: basis
        real(dp), allocatable :: divergences(:), divergences_bar(:)
        real(dp), allocatable :: ref_divergences(:), ref_divergences_bar(:)
        real(dp), allocatable :: ref_values(:, :), ref_values_bar(:, :)
        real(dp), allocatable :: values(:, :), values_bar(:, :)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: determinant, determinant_bar, determinant_jacobian_bar(3, 3)
        real(dp) :: divergence_energy, jacobian(3, 3), jacobian_bar(3, 3)
        real(dp) :: local_jacobian_bar(3, 3), mass_energy, point(3), seed
        integer :: column, dof_count, point_index, row

        vertices_bar = 0.0_dp
        divergence_coefficient_bar = 0.0_dp
        mass_coefficient_bar = 0.0_dp
        status = 1
        if (degree < 0 .or. quadrature_degree < 0) return
        call initialize_tetra_rt(degree, basis, status)
        if (status /= 0) return
        dof_count = tetra_rt_dof_count(basis)
        if (size(matrix_bar, 1) /= dof_count .or. &
            size(matrix_bar, 2) /= dof_count) return
        call tetra_duffy_quadrature( &
            quadrature_degree, x, y, z, weights, status)
        if (status /= 0) return
        call tetra_geometry(vertices, jacobian, determinant, status)
        if (status /= 0) return
        allocate(ref_values(3, dof_count), ref_divergences(dof_count))
        allocate(ref_values_bar(3, dof_count), ref_divergences_bar(dof_count))
        allocate(values(3, dof_count), divergences(dof_count))
        allocate(values_bar(3, dof_count), divergences_bar(dof_count))
        jacobian_bar = 0.0_dp
        determinant_bar = 0.0_dp
        do point_index = 1, size(weights)
            point = [x(point_index), y(point_index), z(point_index)]
            call evaluate_tetra_rt( &
                basis, point, ref_values, ref_divergences, status)
            if (status /= 0) return
            call map_tetra_rt_contravariant( &
                jacobian, ref_values, ref_divergences, values, divergences, &
                status)
            if (status /= 0) return
            values_bar = 0.0_dp
            divergences_bar = 0.0_dp
            do column = 1, dof_count
                do row = 1, dof_count
                    seed = weights(point_index)*matrix_bar(row, column)
                    divergence_energy = divergences(row)*divergences(column)
                    mass_energy = dot_product(values(:, row), values(:, column))
                    determinant_bar = determinant_bar + seed*( &
                        divergence_coefficient*divergence_energy + &
                        mass_coefficient*mass_energy)
                    divergence_coefficient_bar = &
                        divergence_coefficient_bar + &
                        seed*determinant*divergence_energy
                    mass_coefficient_bar = mass_coefficient_bar + &
                        seed*determinant*mass_energy
                    divergences_bar(row) = divergences_bar(row) + &
                        seed*determinant*divergence_coefficient* &
                        divergences(column)
                    divergences_bar(column) = divergences_bar(column) + &
                        seed*determinant*divergence_coefficient*divergences(row)
                    values_bar(:, row) = values_bar(:, row) + &
                        seed*determinant*mass_coefficient*values(:, column)
                    values_bar(:, column) = values_bar(:, column) + &
                        seed*determinant*mass_coefficient*values(:, row)
                end do
            end do
            call map_tetra_rt_contravariant_vjp( &
                jacobian, ref_values, ref_divergences, values_bar, &
                divergences_bar, local_jacobian_bar, ref_values_bar, &
                ref_divergences_bar, status)
            if (status /= 0) return
            jacobian_bar = jacobian_bar + local_jacobian_bar
        end do
        call det3_vjp(jacobian, determinant_bar, determinant_jacobian_bar)
        jacobian_bar = jacobian_bar + determinant_jacobian_bar
        call tetra_jacobian_vjp(jacobian_bar, vertices_bar)
        status = 0
    end subroutine assemble_tetra_rt_div_mass_element_vjp

    subroutine assemble_tetra_rt_div_mass_csc( &
            mesh_vertices, tetrahedra, degree, quadrature_degree, matrix, &
            status, divergence_coefficient, mass_coefficient)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :)
        integer, intent(in) :: degree, quadrature_degree
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: divergence_coefficient
        real(dp), intent(in), optional :: mass_coefficient

        integer, allocatable :: columns(:), face_orientations(:, :)
        integer, allocatable :: face_permutations(:, :, :), faces(:, :)
        integer, allocatable :: global_dofs(:, :), rows(:)
        real(dp), allocatable :: basis_transform(:, :)
        real(dp), allocatable :: element_matrix(:, :), oriented_matrix(:, :)
        real(dp), allocatable :: values(:)
        real(dp) :: divergence_weight, mass_weight, vertices(3, 4)
        integer :: column, dof_count, entry, global_dof_count
        integer :: local_status, node, row, tetrahedron

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Tetrahedral RT sparse assembly failed")
        if (.not. valid_tetra_mesh(mesh_vertices, tetrahedra)) return
        if (degree < 0) return
        if (quadrature_degree < 0) return
        divergence_weight = 1.0_dp
        mass_weight = 1.0_dp
        if (present(divergence_coefficient)) then
            divergence_weight = divergence_coefficient
        end if
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        call build_tetra_rt_dof_map( &
            degree, tetrahedra, faces, global_dofs, face_orientations, &
            face_permutations, local_status)
        if (local_status /= 0) return

        dof_count = size(global_dofs, 1)
        global_dof_count = maxval(global_dofs)
        allocate( &
            rows(dof_count**2*size(tetrahedra, 2)), &
            columns(dof_count**2*size(tetrahedra, 2)), &
            values(dof_count**2*size(tetrahedra, 2)), &
            basis_transform(dof_count, dof_count), &
            oriented_matrix(dof_count, dof_count))
        entry = 0
        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                vertices(:, node) = &
                    mesh_vertices(:, tetrahedra(node, tetrahedron))
            end do
            call assemble_tetra_rt_div_mass_element( &
                vertices, degree, quadrature_degree, element_matrix, &
                local_status, divergence_weight, mass_weight)
            if (local_status /= 0) return
            call build_tetra_rt_basis_transform( &
                degree, face_orientations(:, tetrahedron), &
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
    end subroutine assemble_tetra_rt_div_mass_csc

    subroutine assemble_tetra_rt_div_mass_csc_jvp( &
            mesh_vertices, tetrahedra, degree, quadrature_degree, &
            divergence_coefficient, mass_coefficient, mesh_vertices_dot, &
            divergence_coefficient_dot, mass_coefficient_dot, matrix_dot, &
            status)
        real(dp), intent(in) :: mesh_vertices(:, :), mesh_vertices_dot(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        real(dp), intent(in) :: divergence_coefficient, mass_coefficient
        real(dp), intent(in) :: divergence_coefficient_dot
        real(dp), intent(in) :: mass_coefficient_dot
        type(csc_t), intent(out) :: matrix_dot
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), face_orientations(:, :)
        integer, allocatable :: face_permutations(:, :, :), faces(:, :)
        integer, allocatable :: global_dofs(:, :), rows(:)
        real(dp), allocatable :: basis_transform(:, :), element_dot(:, :)
        real(dp), allocatable :: oriented_dot(:, :), values(:)
        real(dp) :: vertices(3, 4), vertices_dot(3, 4)
        integer :: column, dof_count, entry, global_dof_count
        integer :: local_status, node, row, tetrahedron

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Tetrahedral RT div-mass JVP assembly failed")
        if (.not. valid_tetra_mesh(mesh_vertices, tetrahedra)) return
        if (any(shape(mesh_vertices_dot) /= shape(mesh_vertices))) return
        if (degree < 0 .or. quadrature_degree < 0) return
        call build_tetra_rt_dof_map( &
            degree, tetrahedra, faces, global_dofs, face_orientations, &
            face_permutations, local_status)
        if (local_status /= 0) return
        dof_count = size(global_dofs, 1)
        global_dof_count = maxval(global_dofs)
        allocate( &
            rows(dof_count*dof_count*size(tetrahedra, 2)), &
            columns(dof_count*dof_count*size(tetrahedra, 2)), &
            values(dof_count*dof_count*size(tetrahedra, 2)))
        allocate(basis_transform(dof_count, dof_count))
        allocate(oriented_dot(dof_count, dof_count))
        entry = 0
        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                vertices(:, node) = &
                    mesh_vertices(:, tetrahedra(node, tetrahedron))
                vertices_dot(:, node) = &
                    mesh_vertices_dot(:, tetrahedra(node, tetrahedron))
            end do
            call assemble_tetra_rt_div_mass_element_jvp( &
                vertices, degree, quadrature_degree, divergence_coefficient, &
                mass_coefficient, vertices_dot, divergence_coefficient_dot, &
                mass_coefficient_dot, element_dot, local_status)
            if (local_status /= 0) return
            call build_tetra_rt_basis_transform( &
                degree, face_orientations(:, tetrahedron), &
                face_permutations(:, :, tetrahedron), basis_transform, &
                local_status)
            if (local_status /= 0) return
            oriented_dot = matmul( &
                transpose(basis_transform), &
                matmul(element_dot, basis_transform))
            do column = 1, dof_count
                do row = 1, dof_count
                    entry = entry + 1
                    rows(entry) = global_dofs(row, tetrahedron)
                    columns(entry) = global_dofs(column, tetrahedron)
                    values(entry) = oriented_dot(row, column)
                end do
            end do
        end do
        call csc_from_triplet( &
            global_dof_count, global_dof_count, rows, columns, values, &
            matrix_dot, status)
    end subroutine assemble_tetra_rt_div_mass_csc_jvp

    subroutine assemble_tetra_rt_div_mass_csc_vjp( &
            mesh_vertices, tetrahedra, degree, quadrature_degree, &
            divergence_coefficient, mass_coefficient, matrix_values_bar, &
            mesh_vertices_bar, divergence_coefficient_bar, &
            mass_coefficient_bar, status)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        real(dp), intent(in) :: divergence_coefficient, mass_coefficient
        real(dp), intent(in) :: matrix_values_bar(:)
        real(dp), intent(out) :: mesh_vertices_bar(:, :)
        real(dp), intent(out) :: divergence_coefficient_bar
        real(dp), intent(out) :: mass_coefficient_bar
        type(fortsparse_status_t), intent(out) :: status

        type(csc_t) :: matrix
        integer, allocatable :: face_orientations(:, :), faces(:, :)
        integer, allocatable :: face_permutations(:, :, :)
        integer, allocatable :: global_dofs(:, :)
        real(dp), allocatable :: basis_transform(:, :), dense_bar(:, :)
        real(dp), allocatable :: element_bar(:, :), oriented_bar(:, :)
        real(dp) :: local_divergence_bar, local_mass_bar
        real(dp) :: local_vertices_bar(3, 4), vertices(3, 4)
        integer :: column, dof_count, entry, local_status, node, row
        integer :: tetrahedron

        mesh_vertices_bar = 0.0_dp
        divergence_coefficient_bar = 0.0_dp
        mass_coefficient_bar = 0.0_dp
        call assemble_tetra_rt_div_mass_csc( &
            mesh_vertices, tetrahedra, degree, quadrature_degree, matrix, &
            status, divergence_coefficient, mass_coefficient)
        if (status%code /= 0) return
        if (any(shape(mesh_vertices_bar) /= shape(mesh_vertices)) .or. &
            size(matrix_values_bar) /= matrix%nnz) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Tetrahedral RT div-mass VJP shapes differ")
            return
        end if
        allocate(dense_bar(matrix%nrow, matrix%ncol), source=0.0_dp)
        do column = 1, matrix%ncol
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                dense_bar(matrix%row_idx(entry), column) = &
                    matrix_values_bar(entry)
            end do
        end do
        call build_tetra_rt_dof_map( &
            degree, tetrahedra, faces, global_dofs, face_orientations, &
            face_permutations, local_status)
        if (local_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Tetrahedral RT div-mass VJP dof map failed")
            return
        end if
        dof_count = size(global_dofs, 1)
        allocate(basis_transform(dof_count, dof_count))
        allocate(oriented_bar(dof_count, dof_count))
        allocate(element_bar(dof_count, dof_count))
        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                vertices(:, node) = &
                    mesh_vertices(:, tetrahedra(node, tetrahedron))
            end do
            do column = 1, dof_count
                do row = 1, dof_count
                    oriented_bar(row, column) = dense_bar( &
                        global_dofs(row, tetrahedron), &
                        global_dofs(column, tetrahedron))
                end do
            end do
            call build_tetra_rt_basis_transform( &
                degree, face_orientations(:, tetrahedron), &
                face_permutations(:, :, tetrahedron), basis_transform, &
                local_status)
            if (local_status /= 0) then
                call status_set( &
                    status, FORTSPARSE_INVALID_MATRIX, &
                    "Tetrahedral RT div-mass VJP orientation failed")
                return
            end if
            element_bar = matmul( &
                basis_transform, &
                matmul(oriented_bar, transpose(basis_transform)))
            call assemble_tetra_rt_div_mass_element_vjp( &
                vertices, degree, quadrature_degree, divergence_coefficient, &
                mass_coefficient, element_bar, local_vertices_bar, &
                local_divergence_bar, local_mass_bar, local_status)
            if (local_status /= 0) then
                call status_set( &
                    status, FORTSPARSE_INVALID_MATRIX, &
                    "Tetrahedral RT div-mass element VJP failed")
                return
            end if
            do node = 1, 4
                mesh_vertices_bar(:, tetrahedra(node, tetrahedron)) = &
                    mesh_vertices_bar(:, tetrahedra(node, tetrahedron)) + &
                    local_vertices_bar(:, node)
            end do
            divergence_coefficient_bar = &
                divergence_coefficient_bar + local_divergence_bar
            mass_coefficient_bar = mass_coefficient_bar + local_mass_bar
        end do
    end subroutine assemble_tetra_rt_div_mass_csc_vjp

    pure subroutine tetra_jacobian(vertices, jacobian)
        real(dp), intent(in) :: vertices(3, 4)
        real(dp), intent(out) :: jacobian(3, 3)

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
    end subroutine tetra_jacobian

    pure subroutine tetra_jacobian_vjp(jacobian_bar, vertices_bar)
        real(dp), intent(in) :: jacobian_bar(3, 3)
        real(dp), intent(out) :: vertices_bar(3, 4)

        vertices_bar(:, 1) = -sum(jacobian_bar, dim=2)
        vertices_bar(:, 2) = jacobian_bar(:, 1)
        vertices_bar(:, 3) = jacobian_bar(:, 2)
        vertices_bar(:, 4) = jacobian_bar(:, 3)
    end subroutine tetra_jacobian_vjp

    pure subroutine tetra_geometry( &
            vertices, jacobian, determinant, status)
        real(dp), intent(in) :: vertices(3, 4)
        real(dp), intent(out) :: jacobian(3, 3), determinant
        integer, intent(out) :: status

        real(dp) :: tolerance

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
        determinant = det3(jacobian)
        tolerance = 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**3)
        status = merge(0, 1, determinant > tolerance)
    end subroutine tetra_geometry

    pure function valid_tetra_mesh(mesh_vertices, tetrahedra) result(valid)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :)
        logical :: valid

        valid = size(mesh_vertices, 1) == 3 .and. &
            size(tetrahedra, 1) == 4 .and. size(tetrahedra, 2) > 0
        if (.not. valid) return
        valid = all(tetrahedra >= 1) .and. &
            all(tetrahedra <= size(mesh_vertices, 2))
    end function valid_tetra_mesh

end module fortfem_assembly_tetra_rt_arbitrary_order_3d
