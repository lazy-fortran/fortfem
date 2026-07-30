module fortfem_assembly_rt_arbitrary_order_2d
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_triangle_discontinuous_dof_map, only: &
        build_triangle_discontinuous_dof_map
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_triangle_global_dof_map, only: &
        build_triangle_trimmed_dof_map
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        evaluate_triangle_lagrange_basis, initialize_triangle_lagrange_basis, &
        triangle_lagrange_basis_t
    use fortfem_triangle_piola_maps, only: map_triangle_rt_contravariant
    use fortfem_triangle_rt_arbitrary_order, only: &
        evaluate_triangle_raviart_thomas, initialize_triangle_raviart_thomas, &
        triangle_rt_basis_t, triangle_rt_dof_count
    use fortsparse, only: &
        csc_from_triplet, csc_t, FORTSPARSE_INVALID_MATRIX, &
        fortsparse_status_t, status_set
    implicit none

    private

    public :: assemble_triangle_rt_div_mass_element
    public :: assemble_triangle_rt_div_mass_csc
    public :: assemble_triangle_rt_divergence_csc
    public :: assemble_triangle_rt_cell_vector_load

contains

    subroutine assemble_triangle_rt_cell_vector_load( &
            mesh, degree, quadrature_degree, cell_values, vector, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: cell_values(:, :)
        real(dp), intent(out) :: vector(:)
        type(fortsparse_status_t), intent(out) :: status

        type(triangle_rt_basis_t) :: basis
        integer, allocatable :: global_dofs(:, :), transforms(:, :)
        real(dp), allocatable :: eta(:), local_vector(:)
        real(dp), allocatable :: physical_divergences(:)
        real(dp), allocatable :: physical_values(:, :)
        real(dp), allocatable :: reference_divergences(:)
        real(dp), allocatable :: reference_values(:, :), weights(:), xi(:)
        real(dp) :: determinant, jacobian(2, 2), physical_weight
        real(dp) :: vertices(2, 3)
        integer :: global_dof_count, local_dof, local_dof_count
        integer :: local_status, point, triangle

        vector = 0.0_dp
        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, "RT cell load assembly failed")
        if (degree < 0 .or. quadrature_degree < 0) return
        if (size(cell_values, 1) /= 2 .or. &
            size(cell_values, 2) /= mesh%n_triangles) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "RT cell load requires one vector per triangle")
            return
        end if
        call build_triangle_trimmed_dof_map( &
            mesh, degree + 1, global_dofs, transforms, global_dof_count, &
            local_status)
        if (local_status /= 0 .or. size(vector) /= global_dof_count) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "RT cell load requires a valid output space")
            return
        end if
        call initialize_triangle_raviart_thomas(degree, basis, local_status)
        if (local_status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, local_status)
        if (local_status /= 0) return

        local_dof_count = size(global_dofs, 1)
        allocate(local_vector(local_dof_count))
        allocate(reference_values(2, local_dof_count))
        allocate(reference_divergences(local_dof_count))
        allocate(physical_values(2, local_dof_count))
        allocate(physical_divergences(local_dof_count))
        do triangle = 1, mesh%n_triangles
            vertices(:, 1) = &
                mesh%vertices(:, mesh%triangles(1, triangle))
            vertices(:, 2) = &
                mesh%vertices(:, mesh%triangles(2, triangle))
            vertices(:, 3) = &
                mesh%vertices(:, mesh%triangles(3, triangle))
            jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
            jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
            determinant = jacobian(1, 1) * jacobian(2, 2) - &
                jacobian(1, 2) * jacobian(2, 1)
            if (determinant <= 64.0_dp * epsilon(1.0_dp) * &
                max(1.0_dp, maxval(abs(jacobian))**2)) return
            local_vector = 0.0_dp
            do point = 1, size(weights)
                call evaluate_triangle_raviart_thomas( &
                    basis, xi(point), eta(point), reference_values, &
                    reference_divergences, local_status)
                if (local_status /= 0) return
                call map_triangle_rt_contravariant( &
                    jacobian, reference_values, reference_divergences, &
                    physical_values, physical_divergences, local_status)
                if (local_status /= 0) return
                physical_weight = determinant * weights(point)
                do local_dof = 1, local_dof_count
                    local_vector(local_dof) = local_vector(local_dof) + &
                        physical_weight * dot_product( &
                        cell_values(:, triangle), &
                        physical_values(:, local_dof))
                end do
            end do
            do local_dof = 1, local_dof_count
                vector(global_dofs(local_dof, triangle)) = &
                    vector(global_dofs(local_dof, triangle)) + real( &
                    transforms(local_dof, triangle), dp) * &
                    local_vector(local_dof)
            end do
        end do
        call status_set(status, 0, "")
    end subroutine assemble_triangle_rt_cell_vector_load

    subroutine assemble_triangle_rt_div_mass_element( &
            vertices, degree, quadrature_degree, matrix, status, &
            divergence_coefficient, mass_coefficient)
        real(dp), intent(in) :: vertices(2, 3)
        integer, intent(in) :: degree, quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: divergence_coefficient
        real(dp), intent(in), optional :: mass_coefficient

        type(triangle_rt_basis_t) :: basis
        real(dp), allocatable :: eta(:), physical_divergences(:)
        real(dp), allocatable :: physical_values(:, :)
        real(dp), allocatable :: reference_divergences(:)
        real(dp), allocatable :: reference_values(:, :), weights(:), xi(:)
        real(dp) :: determinant, divergence_weight, jacobian(2, 2)
        real(dp) :: mass_weight, physical_weight
        integer :: basis_dof_count, column, point, row

        status = 1
        if (degree < 0 .or. quadrature_degree < 0) return
        divergence_weight = 1.0_dp
        mass_weight = 1.0_dp
        if (present(divergence_coefficient)) then
            divergence_weight = divergence_coefficient
        end if
        if (present(mass_coefficient)) mass_weight = mass_coefficient

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        determinant = jacobian(1, 1) * jacobian(2, 2) - &
            jacobian(1, 2) * jacobian(2, 1)
        if (determinant <= 64.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp, maxval(abs(jacobian))**2)) return

        call initialize_triangle_raviart_thomas(degree, basis, status)
        if (status /= 0) return
        basis_dof_count = triangle_rt_dof_count(basis)
        allocate(matrix(basis_dof_count, basis_dof_count))
        allocate(reference_values(2, basis_dof_count))
        allocate(reference_divergences(basis_dof_count))
        allocate(physical_values(2, basis_dof_count))
        allocate(physical_divergences(basis_dof_count))
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return

        matrix = 0.0_dp
        do point = 1, size(weights)
            call evaluate_triangle_raviart_thomas( &
                basis, xi(point), eta(point), reference_values, &
                reference_divergences, status)
            if (status /= 0) return
            call map_triangle_rt_contravariant( &
                jacobian, reference_values, reference_divergences, &
                physical_values, physical_divergences, status)
            if (status /= 0) return
            physical_weight = determinant * weights(point)
            do column = 1, basis_dof_count
                do row = 1, basis_dof_count
                    matrix(row, column) = matrix(row, column) + &
                        physical_weight * ( &
                        divergence_weight * physical_divergences(row) * &
                        physical_divergences(column) + mass_weight * &
                        dot_product( &
                        physical_values(:, row), physical_values(:, column)))
                end do
            end do
        end do
        status = 0
    end subroutine assemble_triangle_rt_div_mass_element

    subroutine assemble_triangle_rt_div_mass_csc( &
            mesh, degree, quadrature_degree, matrix, status, &
            divergence_coefficient, mass_coefficient)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: divergence_coefficient
        real(dp), intent(in), optional :: mass_coefficient

        integer, allocatable :: columns(:), global_dofs(:, :), rows(:)
        integer, allocatable :: transforms(:, :)
        real(dp), allocatable :: element_matrix(:, :), values(:)
        real(dp) :: divergence_weight, mass_weight, vertices(2, 3)
        integer :: column, entry, global_dof_count, local_dof_count
        integer :: local_status, row, triangle

        divergence_weight = 1.0_dp
        mass_weight = 1.0_dp
        if (present(divergence_coefficient)) then
            divergence_weight = divergence_coefficient
        end if
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        if (degree < 0 .or. quadrature_degree < 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "RT sparse assembly requires nonnegative degree")
            return
        end if
        call build_triangle_trimmed_dof_map( &
            mesh, degree + 1, global_dofs, transforms, global_dof_count, &
            local_status)
        if (local_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "RT sparse assembly requires a valid triangle mesh")
            return
        end if

        local_dof_count = size(global_dofs, 1)
        allocate(rows(mesh%n_triangles * local_dof_count**2))
        allocate(columns(mesh%n_triangles * local_dof_count**2))
        allocate(values(mesh%n_triangles * local_dof_count**2))
        entry = 0
        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            call assemble_triangle_rt_div_mass_element( &
                vertices, degree, quadrature_degree, element_matrix, &
                local_status, divergence_weight, mass_weight)
            if (local_status /= 0) then
                call status_set( &
                    status, FORTSPARSE_INVALID_MATRIX, &
                    "RT sparse assembly requires valid CCW triangles")
                return
            end if
            do column = 1, local_dof_count
                do row = 1, local_dof_count
                    entry = entry + 1
                    rows(entry) = global_dofs(row, triangle)
                    columns(entry) = global_dofs(column, triangle)
                    values(entry) = real( &
                        transforms(row, triangle) * &
                        transforms(column, triangle), dp) * &
                        element_matrix(row, column)
                end do
            end do
        end do
        call csc_from_triplet( &
            global_dof_count, global_dof_count, rows, columns, values, &
            matrix, status)
    end subroutine assemble_triangle_rt_div_mass_csc

    subroutine assemble_triangle_rt_divergence_csc( &
            mesh, degree, quadrature_degree, matrix, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        type(triangle_lagrange_basis_t) :: scalar_basis
        type(triangle_rt_basis_t) :: rt_basis
        integer, allocatable :: columns(:), rt_global_dofs(:, :), rows(:)
        integer, allocatable :: rt_transforms(:, :)
        integer, allocatable :: scalar_global_dofs(:, :)
        real(dp), allocatable :: eta(:), physical_divergences(:)
        real(dp), allocatable :: physical_values(:, :)
        real(dp), allocatable :: reference_divergences(:)
        real(dp), allocatable :: reference_values(:, :)
        real(dp), allocatable :: local_matrix(:, :)
        real(dp), allocatable :: scalar_gradients(:, :), scalar_values(:)
        real(dp), allocatable :: values(:), weights(:), xi(:)
        real(dp) :: determinant, jacobian(2, 2)
        real(dp) :: vertices(2, 3)
        integer :: entry, local_status, point, rt_dof, rt_dof_count
        integer :: rt_local_dof_count, scalar_dof, scalar_dof_count
        integer :: scalar_local_dof_count, triangle

        if (degree < 0 .or. quadrature_degree < 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "RT-DG divergence assembly requires nonnegative degree")
            return
        end if
        call build_triangle_trimmed_dof_map( &
            mesh, degree + 1, rt_global_dofs, rt_transforms, rt_dof_count, &
            local_status)
        if (local_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "RT-DG divergence assembly requires a valid triangle mesh")
            return
        end if
        call build_triangle_discontinuous_dof_map( &
            mesh, degree, scalar_global_dofs, scalar_dof_count, local_status)
        if (local_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "RT-DG divergence assembly could not map scalar degrees")
            return
        end if
        call initialize_triangle_raviart_thomas( &
            degree, rt_basis, local_status)
        if (local_status /= 0) return
        call initialize_triangle_lagrange_basis( &
            degree, scalar_basis, local_status)
        if (local_status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, local_status)
        if (local_status /= 0) return

        rt_local_dof_count = size(rt_global_dofs, 1)
        scalar_local_dof_count = size(scalar_global_dofs, 1)
        allocate(reference_values(2, rt_local_dof_count))
        allocate(reference_divergences(rt_local_dof_count))
        allocate(physical_values(2, rt_local_dof_count))
        allocate(physical_divergences(rt_local_dof_count))
        allocate(scalar_values(scalar_local_dof_count))
        allocate(scalar_gradients(2, scalar_local_dof_count))
        allocate(local_matrix(scalar_local_dof_count, rt_local_dof_count))
        allocate(rows( &
            mesh%n_triangles * rt_local_dof_count * scalar_local_dof_count))
        allocate(columns(size(rows)), values(size(rows)))

        entry = 0
        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
            jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
            determinant = jacobian(1, 1) * jacobian(2, 2) - &
                jacobian(1, 2) * jacobian(2, 1)
            if (determinant <= 0.0_dp) then
                call status_set( &
                    status, FORTSPARSE_INVALID_MATRIX, &
                    "RT-DG divergence assembly requires valid CCW triangles")
                return
            end if
            local_matrix = 0.0_dp
            do point = 1, size(weights)
                call evaluate_triangle_raviart_thomas( &
                    rt_basis, xi(point), eta(point), reference_values, &
                    reference_divergences, local_status)
                if (local_status /= 0) return
                call map_triangle_rt_contravariant( &
                    jacobian, reference_values, reference_divergences, &
                    physical_values, physical_divergences, local_status)
                if (local_status /= 0) return
                call evaluate_triangle_lagrange_basis( &
                    scalar_basis, xi(point), eta(point), scalar_values, &
                    scalar_gradients, local_status)
                if (local_status /= 0) return
                do rt_dof = 1, rt_local_dof_count
                    do scalar_dof = 1, scalar_local_dof_count
                        local_matrix(scalar_dof, rt_dof) = &
                            local_matrix(scalar_dof, rt_dof) + &
                            determinant * weights(point) * &
                            scalar_values(scalar_dof) * &
                            physical_divergences(rt_dof)
                    end do
                end do
            end do
            do rt_dof = 1, rt_local_dof_count
                do scalar_dof = 1, scalar_local_dof_count
                    entry = entry + 1
                    rows(entry) = scalar_global_dofs(scalar_dof, triangle)
                    columns(entry) = rt_global_dofs(rt_dof, triangle)
                    values(entry) = real( &
                        rt_transforms(rt_dof, triangle), dp) * &
                        local_matrix(scalar_dof, rt_dof)
                end do
            end do
        end do
        call csc_from_triplet( &
            scalar_dof_count, rt_dof_count, rows, columns, values, &
            matrix, status)
    end subroutine assemble_triangle_rt_divergence_csc

end module fortfem_assembly_rt_arbitrary_order_2d
