module fortfem_assembly_nedelec_arbitrary_order_2d
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
    use fortfem_triangle_nedelec_arbitrary_order, only: &
        evaluate_triangle_nedelec_first_kind, &
        initialize_triangle_nedelec_first_kind, triangle_nedelec_dof_count, &
        triangle_nedelec_first_kind_t
    use fortfem_triangle_piola_maps, only: &
        map_triangle_nedelec_covariant, &
        map_triangle_nedelec_covariant_jvp, &
        map_triangle_nedelec_covariant_vjp
    use fortnum_linalg, only: det2_jvp, det2_vjp
    use fortsparse, only: &
        csc_from_triplet, csc_t, FORTSPARSE_INVALID_MATRIX, &
        fortsparse_status_t, status_set
    implicit none

    private

    public :: assemble_triangle_nedelec_curl_mass_element
    public :: assemble_triangle_nedelec_curl_mass_element_jvp
    public :: assemble_triangle_nedelec_curl_mass_element_vjp
    public :: assemble_triangle_nedelec_curl_mass_csc
    public :: assemble_triangle_nedelec_curl_mass_csc_jvp
    public :: assemble_triangle_nedelec_curl_mass_csc_vjp
    public :: assemble_triangle_nedelec_cell_tensor_csc
    public :: assemble_triangle_nedelec_curl_csc
    public :: assemble_triangle_nedelec_cell_vector_load
    public :: assemble_triangle_nedelec_vector_load_samples
    public :: assemble_triangle_nedelec_vector_load_samples_jvp
    public :: assemble_triangle_nedelec_vector_load_samples_vjp

contains

    subroutine assemble_triangle_nedelec_cell_vector_load( &
            mesh, order, quadrature_degree, cell_values, vector, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: order, quadrature_degree
        real(dp), intent(in) :: cell_values(:, :)
        real(dp), intent(out) :: vector(:)
        type(fortsparse_status_t), intent(out) :: status

        type(triangle_nedelec_first_kind_t) :: basis
        integer, allocatable :: global_dofs(:, :), transforms(:, :)
        real(dp), allocatable :: eta(:), local_vector(:), physical_curls(:)
        real(dp), allocatable :: physical_values(:, :), reference_curls(:)
        real(dp), allocatable :: reference_values(:, :), weights(:), xi(:)
        real(dp) :: determinant, jacobian(2, 2), physical_weight
        real(dp) :: vertices(2, 3)
        integer :: global_dof_count, local_dof, local_dof_count
        integer :: local_status, point, triangle

        vector = 0.0_dp
        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Nedelec cell load assembly failed")
        if (order < 1 .or. quadrature_degree < 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Nedelec cell load requires positive order")
            return
        end if
        if (size(cell_values, 1) /= 2 .or. &
            size(cell_values, 2) /= mesh%n_triangles) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Nedelec cell load requires one vector per triangle")
            return
        end if
        call build_triangle_trimmed_dof_map( &
            mesh, order, global_dofs, transforms, global_dof_count, &
            local_status)
        if (local_status /= 0 .or. size(vector) /= global_dof_count) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Nedelec cell load requires a valid output space")
            return
        end if
        call initialize_triangle_nedelec_first_kind(order, basis, local_status)
        if (local_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Nedelec cell load could not initialize its basis")
            return
        end if
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, local_status)
        if (local_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Nedelec cell load could not build its quadrature")
            return
        end if

        local_dof_count = size(global_dofs, 1)
        allocate(local_vector(local_dof_count))
        allocate(reference_values(2, local_dof_count))
        allocate(reference_curls(local_dof_count))
        allocate(physical_values(2, local_dof_count))
        allocate(physical_curls(local_dof_count))
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
                max(1.0_dp, maxval(abs(jacobian))**2)) then
                call status_set( &
                    status, FORTSPARSE_INVALID_MATRIX, &
                    "Nedelec cell load requires valid CCW triangles")
                return
            end if
            local_vector = 0.0_dp
            do point = 1, size(weights)
                call evaluate_triangle_nedelec_first_kind( &
                    basis, xi(point), eta(point), reference_values, &
                    reference_curls, local_status)
                if (local_status /= 0) return
                call map_triangle_nedelec_covariant( &
                    jacobian, reference_values, reference_curls, &
                    physical_values, physical_curls, local_status)
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
    end subroutine assemble_triangle_nedelec_cell_vector_load

    subroutine assemble_triangle_nedelec_vector_load_samples( &
            mesh, order, quadrature_degree, source_values, vector, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: order, quadrature_degree
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), allocatable, intent(out) :: vector(:)
        type(fortsparse_status_t), intent(out) :: status

        type(triangle_nedelec_first_kind_t) :: basis
        integer, allocatable :: global_dofs(:, :), transforms(:, :)
        real(dp), allocatable :: curls(:), eta(:), local_vector(:)
        real(dp), allocatable :: ref_curls(:), ref_values(:, :), values(:, :)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp) :: determinant, jacobian(2, 2), vertices(2, 3)
        integer :: dof, dof_count, global_count, local_status, point, triangle

        call initialize_sampled_nedelec_load( &
            mesh, order, quadrature_degree, source_values, basis, &
            global_dofs, transforms, global_count, xi, eta, weights, status)
        if (status%code /= 0) return
        dof_count = size(global_dofs, 1)
        allocate(vector(global_count), source=0.0_dp)
        allocate(local_vector(dof_count))
        allocate(ref_values(2, dof_count), ref_curls(dof_count))
        allocate(values(2, dof_count), curls(dof_count))
        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            call triangle_jacobian(vertices, jacobian, determinant, local_status)
            if (local_status /= 0) return
            local_vector = 0.0_dp
            do point = 1, size(weights)
                call evaluate_triangle_nedelec_first_kind( &
                    basis, xi(point), eta(point), ref_values, ref_curls, &
                    local_status)
                if (local_status /= 0) return
                call map_triangle_nedelec_covariant( &
                    jacobian, ref_values, ref_curls, values, curls, &
                    local_status)
                if (local_status /= 0) return
                do dof = 1, dof_count
                    local_vector(dof) = local_vector(dof) + &
                        determinant*weights(point)*dot_product( &
                        source_values(:, point, triangle), values(:, dof))
                end do
            end do
            do dof = 1, dof_count
                vector(global_dofs(dof, triangle)) = &
                    vector(global_dofs(dof, triangle)) + &
                    real(transforms(dof, triangle), dp)*local_vector(dof)
            end do
        end do
        call status_set(status, 0, "")
    end subroutine assemble_triangle_nedelec_vector_load_samples

    subroutine assemble_triangle_nedelec_vector_load_samples_jvp( &
            mesh, order, quadrature_degree, source_values, source_gradients, &
            mesh_vertices_dot, source_parameter_dot, vector_dot, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: order, quadrature_degree
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :)
        real(dp), intent(in) :: mesh_vertices_dot(:, :)
        real(dp), intent(in) :: source_parameter_dot(:, :, :)
        real(dp), allocatable, intent(out) :: vector_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        type(triangle_nedelec_first_kind_t) :: basis
        integer, allocatable :: global_dofs(:, :), transforms(:, :)
        real(dp), allocatable :: curls(:), curls_dot(:), eta(:), local_dot(:)
        real(dp), allocatable :: ref_curls(:), ref_values(:, :), values(:, :)
        real(dp), allocatable :: values_dot(:, :), weights(:), xi(:)
        real(dp), allocatable :: zero_curls(:), zero_values(:, :)
        real(dp) :: determinant, determinant_dot, jacobian(2, 2)
        real(dp) :: jacobian_dot(2, 2), point_dot(2), reference_point(2)
        real(dp) :: source_dot(2), vertices(2, 3), vertices_dot(2, 3)
        integer :: dof, dof_count, global_count, local_status, point, triangle

        call initialize_sampled_nedelec_load( &
            mesh, order, quadrature_degree, source_values, basis, &
            global_dofs, transforms, global_count, xi, eta, weights, status)
        if (status%code /= 0) return
        if (any(shape(mesh_vertices_dot) /= shape(mesh%vertices))) return
        if (.not. valid_vector_sample_products( &
            source_values, source_gradients, source_parameter_dot, &
            size(weights), mesh%n_triangles)) return
        dof_count = size(global_dofs, 1)
        allocate(vector_dot(global_count), source=0.0_dp)
        allocate(local_dot(dof_count))
        allocate(ref_values(2, dof_count), ref_curls(dof_count))
        allocate(values(2, dof_count), curls(dof_count))
        allocate(values_dot(2, dof_count), curls_dot(dof_count))
        allocate(zero_values(2, dof_count), zero_curls(dof_count), &
            source=0.0_dp)
        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            vertices_dot = &
                mesh_vertices_dot(:, mesh%triangles(:, triangle))
            call triangle_jacobian(vertices, jacobian, determinant, local_status)
            if (local_status /= 0) return
            call triangle_jacobian_direction(vertices_dot, jacobian_dot)
            call det2_jvp(jacobian, jacobian_dot, determinant_dot)
            local_dot = 0.0_dp
            do point = 1, size(weights)
                reference_point = [xi(point), eta(point)]
                call evaluate_triangle_nedelec_first_kind( &
                    basis, xi(point), eta(point), ref_values, ref_curls, &
                    local_status)
                if (local_status /= 0) return
                call map_triangle_nedelec_covariant( &
                    jacobian, ref_values, ref_curls, values, curls, &
                    local_status)
                if (local_status /= 0) return
                call map_triangle_nedelec_covariant_jvp( &
                    jacobian, ref_values, ref_curls, jacobian_dot, zero_values, &
                    zero_curls, values_dot, curls_dot, local_status)
                if (local_status /= 0) return
                point_dot = vertices_dot(:, 1) + &
                    matmul(jacobian_dot, reference_point)
                source_dot = source_parameter_dot(:, point, triangle) + &
                    matmul(source_gradients(:, :, point, triangle), point_dot)
                do dof = 1, dof_count
                    local_dot(dof) = local_dot(dof) + weights(point)*( &
                        determinant_dot*dot_product( &
                        source_values(:, point, triangle), values(:, dof)) + &
                        determinant*(dot_product(source_dot, values(:, dof)) + &
                        dot_product(source_values(:, point, triangle), &
                        values_dot(:, dof))))
                end do
            end do
            do dof = 1, dof_count
                vector_dot(global_dofs(dof, triangle)) = &
                    vector_dot(global_dofs(dof, triangle)) + &
                    real(transforms(dof, triangle), dp)*local_dot(dof)
            end do
        end do
        call status_set(status, 0, "")
    end subroutine assemble_triangle_nedelec_vector_load_samples_jvp

    subroutine assemble_triangle_nedelec_vector_load_samples_vjp( &
            mesh, order, quadrature_degree, source_values, source_gradients, &
            vector_bar, mesh_vertices_bar, source_values_bar, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: order, quadrature_degree
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :), vector_bar(:)
        real(dp), intent(out) :: mesh_vertices_bar(:, :)
        real(dp), intent(out) :: source_values_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        type(triangle_nedelec_first_kind_t) :: basis
        integer, allocatable :: global_dofs(:, :), transforms(:, :)
        real(dp), allocatable :: curls(:), curls_bar(:), eta(:), local_bar(:)
        real(dp), allocatable :: ref_curls(:), ref_curls_bar(:)
        real(dp), allocatable :: ref_values(:, :), ref_values_bar(:, :)
        real(dp), allocatable :: values(:, :), values_bar(:, :)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp) :: determinant, determinant_bar, determinant_jacobian_bar(2, 2)
        real(dp) :: jacobian(2, 2), jacobian_bar(2, 2)
        real(dp) :: local_jacobian_bar(2, 2), local_vertices_bar(2, 3)
        real(dp) :: point_bar(2), reference_point(2), sample_bar(2), seed
        real(dp) :: vertices(2, 3)
        integer :: dof, dof_count, global_count, local_status, node
        integer :: point, triangle

        mesh_vertices_bar = 0.0_dp
        source_values_bar = 0.0_dp
        call initialize_sampled_nedelec_load( &
            mesh, order, quadrature_degree, source_values, basis, &
            global_dofs, transforms, global_count, xi, eta, weights, status)
        if (status%code /= 0) return
        if (any(shape(mesh_vertices_bar) /= shape(mesh%vertices))) return
        if (any(shape(source_values_bar) /= shape(source_values))) return
        if (.not. valid_vector_sample_gradients( &
            source_values, source_gradients, size(weights), &
            mesh%n_triangles)) return
        if (size(vector_bar) /= global_count) return
        dof_count = size(global_dofs, 1)
        allocate(local_bar(dof_count))
        allocate(ref_values(2, dof_count), ref_curls(dof_count))
        allocate(ref_values_bar(2, dof_count), ref_curls_bar(dof_count))
        allocate(values(2, dof_count), curls(dof_count))
        allocate(values_bar(2, dof_count), curls_bar(dof_count))
        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            call triangle_jacobian(vertices, jacobian, determinant, local_status)
            if (local_status /= 0) return
            do dof = 1, dof_count
                local_bar(dof) = real(transforms(dof, triangle), dp)* &
                    vector_bar(global_dofs(dof, triangle))
            end do
            local_vertices_bar = 0.0_dp
            jacobian_bar = 0.0_dp
            determinant_bar = 0.0_dp
            do point = 1, size(weights)
                reference_point = [xi(point), eta(point)]
                call evaluate_triangle_nedelec_first_kind( &
                    basis, xi(point), eta(point), ref_values, ref_curls, &
                    local_status)
                if (local_status /= 0) return
                call map_triangle_nedelec_covariant( &
                    jacobian, ref_values, ref_curls, values, curls, &
                    local_status)
                if (local_status /= 0) return
                values_bar = 0.0_dp
                curls_bar = 0.0_dp
                sample_bar = 0.0_dp
                do dof = 1, dof_count
                    seed = weights(point)*local_bar(dof)
                    determinant_bar = determinant_bar + seed*dot_product( &
                        source_values(:, point, triangle), values(:, dof))
                    sample_bar = sample_bar + seed*determinant*values(:, dof)
                    values_bar(:, dof) = values_bar(:, dof) + &
                        seed*determinant*source_values(:, point, triangle)
                end do
                source_values_bar(:, point, triangle) = sample_bar
                point_bar = matmul(transpose( &
                    source_gradients(:, :, point, triangle)), sample_bar)
                local_vertices_bar(:, 1) = &
                    local_vertices_bar(:, 1) + point_bar
                jacobian_bar = jacobian_bar + spread(point_bar, 2, 2)* &
                    spread(reference_point, 1, 2)
                call map_triangle_nedelec_covariant_vjp( &
                    jacobian, ref_values, ref_curls, values_bar, curls_bar, &
                    local_jacobian_bar, ref_values_bar, ref_curls_bar, &
                    local_status)
                if (local_status /= 0) return
                jacobian_bar = jacobian_bar + local_jacobian_bar
            end do
            call det2_vjp( &
                jacobian, determinant_bar, determinant_jacobian_bar)
            call triangle_jacobian_vjp_add( &
                jacobian_bar + determinant_jacobian_bar, local_vertices_bar)
            do node = 1, 3
                mesh_vertices_bar(:, mesh%triangles(node, triangle)) = &
                    mesh_vertices_bar(:, mesh%triangles(node, triangle)) + &
                    local_vertices_bar(:, node)
            end do
        end do
        call status_set(status, 0, "")
    end subroutine assemble_triangle_nedelec_vector_load_samples_vjp

    subroutine assemble_triangle_nedelec_curl_mass_element( &
            vertices, order, quadrature_degree, matrix, status, &
            curl_coefficient, mass_coefficient, mass_tensor)
        real(dp), intent(in) :: vertices(2, 3)
        integer, intent(in) :: order, quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: curl_coefficient, mass_coefficient
        real(dp), intent(in), optional :: mass_tensor(2, 2)

        type(triangle_nedelec_first_kind_t) :: basis
        real(dp), allocatable :: eta(:), physical_curls(:)
        real(dp), allocatable :: physical_values(:, :), reference_curls(:)
        real(dp), allocatable :: reference_values(:, :), weights(:), xi(:)
        real(dp) :: curl_weight, determinant, jacobian(2, 2)
        real(dp) :: mass_weight, physical_weight, tensor(2, 2)
        integer :: basis_dof_count, column, point, row

        status = 1
        if (order < 1 .or. quadrature_degree < 0) return
        curl_weight = 1.0_dp
        mass_weight = 1.0_dp
        if (present(curl_coefficient)) curl_weight = curl_coefficient
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        tensor = 0.0_dp
        tensor(1, 1) = mass_weight
        tensor(2, 2) = mass_weight
        if (present(mass_tensor)) tensor = mass_tensor

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        determinant = jacobian(1, 1) * jacobian(2, 2) - &
            jacobian(1, 2) * jacobian(2, 1)
        if (determinant <= 64.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp, maxval(abs(jacobian))**2)) return

        call initialize_triangle_nedelec_first_kind(order, basis, status)
        if (status /= 0) return
        basis_dof_count = triangle_nedelec_dof_count(basis)
        allocate(matrix(basis_dof_count, basis_dof_count))
        allocate(reference_values(2, basis_dof_count))
        allocate(reference_curls(basis_dof_count))
        allocate(physical_values(2, basis_dof_count))
        allocate(physical_curls(basis_dof_count))
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return

        matrix = 0.0_dp
        do point = 1, size(weights)
            call evaluate_triangle_nedelec_first_kind( &
                basis, xi(point), eta(point), reference_values, &
                reference_curls, status)
            if (status /= 0) return
            call map_triangle_nedelec_covariant( &
                jacobian, reference_values, reference_curls, physical_values, &
                physical_curls, status)
            if (status /= 0) return
            physical_weight = determinant * weights(point)
            do column = 1, basis_dof_count
                do row = 1, basis_dof_count
                    matrix(row, column) = matrix(row, column) + &
                        physical_weight * ( &
                        curl_weight * physical_curls(row) * &
                        physical_curls(column) + dot_product( &
                        physical_values(:, row), &
                        matmul(tensor, physical_values(:, column))))
                end do
            end do
        end do
        status = 0
    end subroutine assemble_triangle_nedelec_curl_mass_element

    subroutine assemble_triangle_nedelec_curl_mass_element_jvp( &
            vertices, order, quadrature_degree, curl_coefficient, mass_tensor, &
            vertices_dot, curl_coefficient_dot, mass_tensor_dot, matrix_dot, &
            status)
        real(dp), intent(in) :: vertices(2, 3), vertices_dot(2, 3)
        integer, intent(in) :: order, quadrature_degree
        real(dp), intent(in) :: curl_coefficient, mass_tensor(2, 2)
        real(dp), intent(in) :: curl_coefficient_dot, mass_tensor_dot(2, 2)
        real(dp), allocatable, intent(out) :: matrix_dot(:, :)
        integer, intent(out) :: status

        type(triangle_nedelec_first_kind_t) :: basis
        real(dp), allocatable :: curls(:), curls_dot(:), eta(:)
        real(dp), allocatable :: ref_curls(:), ref_values(:, :)
        real(dp), allocatable :: values(:, :), values_dot(:, :)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp), allocatable :: zero_curls(:), zero_values(:, :)
        real(dp) :: curl_energy, determinant, determinant_dot
        real(dp) :: jacobian(2, 2), jacobian_dot(2, 2)
        real(dp) :: mass_energy, mass_energy_dot
        integer :: column, dof_count, point, row

        status = 1
        if (allocated(matrix_dot)) deallocate(matrix_dot)
        if (order < 1 .or. quadrature_degree < 0) return
        call initialize_triangle_nedelec_first_kind(order, basis, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        call triangle_jacobian(vertices, jacobian, determinant, status)
        if (status /= 0) return
        call triangle_jacobian_direction(vertices_dot, jacobian_dot)
        call det2_jvp(jacobian, jacobian_dot, determinant_dot)
        dof_count = triangle_nedelec_dof_count(basis)
        allocate(matrix_dot(dof_count, dof_count), source=0.0_dp)
        allocate(ref_values(2, dof_count), ref_curls(dof_count))
        allocate(zero_values(2, dof_count), zero_curls(dof_count), &
            source=0.0_dp)
        allocate(values(2, dof_count), curls(dof_count))
        allocate(values_dot(2, dof_count), curls_dot(dof_count))
        do point = 1, size(weights)
            call evaluate_triangle_nedelec_first_kind( &
                basis, xi(point), eta(point), ref_values, ref_curls, status)
            if (status /= 0) return
            call map_triangle_nedelec_covariant( &
                jacobian, ref_values, ref_curls, values, curls, status)
            if (status /= 0) return
            call map_triangle_nedelec_covariant_jvp( &
                jacobian, ref_values, ref_curls, jacobian_dot, zero_values, &
                zero_curls, values_dot, curls_dot, status)
            if (status /= 0) return
            do column = 1, dof_count
                do row = 1, dof_count
                    curl_energy = curls(row)*curls(column)
                    mass_energy = dot_product( &
                        values(:, row), matmul(mass_tensor, values(:, column)))
                    mass_energy_dot = dot_product( &
                        values_dot(:, row), &
                        matmul(mass_tensor, values(:, column))) + dot_product( &
                        values(:, row), &
                        matmul(mass_tensor_dot, values(:, column)) + &
                        matmul(mass_tensor, values_dot(:, column)))
                    matrix_dot(row, column) = matrix_dot(row, column) + &
                        weights(point)*(determinant_dot*( &
                        curl_coefficient*curl_energy + mass_energy) + &
                        determinant*(curl_coefficient_dot*curl_energy + &
                        curl_coefficient*(curls_dot(row)*curls(column) + &
                        curls(row)*curls_dot(column)) + mass_energy_dot))
                end do
            end do
        end do
        status = 0
    end subroutine assemble_triangle_nedelec_curl_mass_element_jvp

    subroutine assemble_triangle_nedelec_curl_mass_element_vjp( &
            vertices, order, quadrature_degree, curl_coefficient, mass_tensor, &
            matrix_bar, vertices_bar, curl_coefficient_bar, mass_tensor_bar, &
            status)
        real(dp), intent(in) :: vertices(2, 3)
        integer, intent(in) :: order, quadrature_degree
        real(dp), intent(in) :: curl_coefficient, mass_tensor(2, 2)
        real(dp), intent(in) :: matrix_bar(:, :)
        real(dp), intent(out) :: vertices_bar(2, 3)
        real(dp), intent(out) :: curl_coefficient_bar, mass_tensor_bar(2, 2)
        integer, intent(out) :: status

        type(triangle_nedelec_first_kind_t) :: basis
        real(dp), allocatable :: curls(:), curls_bar(:), eta(:)
        real(dp), allocatable :: ref_curls(:), ref_curls_bar(:)
        real(dp), allocatable :: ref_values(:, :), ref_values_bar(:, :)
        real(dp), allocatable :: values(:, :), values_bar(:, :)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp) :: curl_energy, determinant, determinant_bar
        real(dp) :: determinant_jacobian_bar(2, 2), jacobian(2, 2)
        real(dp) :: jacobian_bar(2, 2), local_jacobian_bar(2, 2)
        real(dp) :: mass_energy, seed
        integer :: column, dof_count, point, row

        vertices_bar = 0.0_dp
        curl_coefficient_bar = 0.0_dp
        mass_tensor_bar = 0.0_dp
        status = 1
        if (order < 1 .or. quadrature_degree < 0) return
        call initialize_triangle_nedelec_first_kind(order, basis, status)
        if (status /= 0) return
        dof_count = triangle_nedelec_dof_count(basis)
        if (size(matrix_bar, 1) /= dof_count .or. &
            size(matrix_bar, 2) /= dof_count) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        call triangle_jacobian(vertices, jacobian, determinant, status)
        if (status /= 0) return
        allocate(ref_values(2, dof_count), ref_curls(dof_count))
        allocate(ref_values_bar(2, dof_count), ref_curls_bar(dof_count))
        allocate(values(2, dof_count), curls(dof_count))
        allocate(values_bar(2, dof_count), curls_bar(dof_count))
        jacobian_bar = 0.0_dp
        determinant_bar = 0.0_dp
        do point = 1, size(weights)
            call evaluate_triangle_nedelec_first_kind( &
                basis, xi(point), eta(point), ref_values, ref_curls, status)
            if (status /= 0) return
            call map_triangle_nedelec_covariant( &
                jacobian, ref_values, ref_curls, values, curls, status)
            if (status /= 0) return
            values_bar = 0.0_dp
            curls_bar = 0.0_dp
            do column = 1, dof_count
                do row = 1, dof_count
                    seed = weights(point)*matrix_bar(row, column)
                    curl_energy = curls(row)*curls(column)
                    mass_energy = dot_product( &
                        values(:, row), matmul(mass_tensor, values(:, column)))
                    determinant_bar = determinant_bar + &
                        seed*(curl_coefficient*curl_energy + mass_energy)
                    curl_coefficient_bar = curl_coefficient_bar + &
                        seed*determinant*curl_energy
                    mass_tensor_bar = mass_tensor_bar + &
                        seed*determinant*spread(values(:, row), 2, 2)* &
                        spread(values(:, column), 1, 2)
                    curls_bar(row) = curls_bar(row) + &
                        seed*determinant*curl_coefficient*curls(column)
                    curls_bar(column) = curls_bar(column) + &
                        seed*determinant*curl_coefficient*curls(row)
                    values_bar(:, row) = values_bar(:, row) + &
                        seed*determinant*matmul( &
                        mass_tensor, values(:, column))
                    values_bar(:, column) = values_bar(:, column) + &
                        seed*determinant*matmul( &
                        transpose(mass_tensor), values(:, row))
                end do
            end do
            call map_triangle_nedelec_covariant_vjp( &
                jacobian, ref_values, ref_curls, values_bar, curls_bar, &
                local_jacobian_bar, ref_values_bar, ref_curls_bar, status)
            if (status /= 0) return
            jacobian_bar = jacobian_bar + local_jacobian_bar
        end do
        call det2_vjp(jacobian, determinant_bar, determinant_jacobian_bar)
        call triangle_jacobian_vjp( &
            jacobian_bar + determinant_jacobian_bar, vertices_bar)
        status = 0
    end subroutine assemble_triangle_nedelec_curl_mass_element_vjp

    subroutine assemble_triangle_nedelec_curl_mass_csc( &
            mesh, order, quadrature_degree, matrix, status, &
            curl_coefficient, mass_coefficient, mass_tensor)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: order, quadrature_degree
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: curl_coefficient, mass_coefficient
        real(dp), intent(in), optional :: mass_tensor(2, 2)

        integer, allocatable :: columns(:), global_dofs(:, :), rows(:)
        integer, allocatable :: transforms(:, :)
        real(dp), allocatable :: element_matrix(:, :), values(:)
        real(dp) :: curl_weight, mass_weight, tensor(2, 2), vertices(2, 3)
        integer :: column, entry, global_dof_count, local_dof_count
        integer :: local_status, row, triangle

        curl_weight = 1.0_dp
        mass_weight = 1.0_dp
        if (present(curl_coefficient)) curl_weight = curl_coefficient
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        tensor = 0.0_dp
        tensor(1, 1) = mass_weight
        tensor(2, 2) = mass_weight
        if (present(mass_tensor)) tensor = mass_tensor
        if (order < 1 .or. quadrature_degree < 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Nedelec sparse assembly requires positive order")
            return
        end if
        call build_triangle_trimmed_dof_map( &
            mesh, order, global_dofs, transforms, global_dof_count, &
            local_status)
        if (local_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Nedelec sparse assembly requires a valid triangle mesh")
            return
        end if

        local_dof_count = size(global_dofs, 1)
        allocate(rows(mesh%n_triangles * local_dof_count**2))
        allocate(columns(mesh%n_triangles * local_dof_count**2))
        allocate(values(mesh%n_triangles * local_dof_count**2))
        entry = 0
        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            call assemble_triangle_nedelec_curl_mass_element( &
                vertices, order, quadrature_degree, element_matrix, &
                local_status, curl_weight, mass_tensor=tensor)
            if (local_status /= 0) then
                call status_set( &
                    status, FORTSPARSE_INVALID_MATRIX, &
                    "Nedelec sparse assembly requires valid CCW triangles")
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
    end subroutine assemble_triangle_nedelec_curl_mass_csc

    subroutine assemble_triangle_nedelec_curl_mass_csc_jvp( &
            mesh, order, quadrature_degree, curl_coefficient, mass_tensor, &
            mesh_vertices_dot, curl_coefficient_dot, mass_tensor_dot, &
            matrix_dot, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: order, quadrature_degree
        real(dp), intent(in) :: curl_coefficient, mass_tensor(2, 2)
        real(dp), intent(in) :: mesh_vertices_dot(:, :)
        real(dp), intent(in) :: curl_coefficient_dot, mass_tensor_dot(2, 2)
        type(csc_t), intent(out) :: matrix_dot
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), global_dofs(:, :), rows(:)
        integer, allocatable :: transforms(:, :)
        real(dp), allocatable :: element_dot(:, :), values(:)
        real(dp) :: vertices(2, 3), vertices_dot(2, 3)
        integer :: column, entry, global_dof_count, local_dof_count
        integer :: local_status, row, triangle

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Triangle Nedelec curl-mass JVP assembly failed")
        if (order < 1 .or. quadrature_degree < 0) return
        if (any(shape(mesh_vertices_dot) /= shape(mesh%vertices))) return
        call build_triangle_trimmed_dof_map( &
            mesh, order, global_dofs, transforms, global_dof_count, &
            local_status)
        if (local_status /= 0) return
        local_dof_count = size(global_dofs, 1)
        allocate(rows(mesh%n_triangles*local_dof_count**2))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            vertices_dot = &
                mesh_vertices_dot(:, mesh%triangles(:, triangle))
            call assemble_triangle_nedelec_curl_mass_element_jvp( &
                vertices, order, quadrature_degree, curl_coefficient, &
                mass_tensor, vertices_dot, curl_coefficient_dot, &
                mass_tensor_dot, element_dot, local_status)
            if (local_status /= 0) return
            do column = 1, local_dof_count
                do row = 1, local_dof_count
                    entry = entry + 1
                    rows(entry) = global_dofs(row, triangle)
                    columns(entry) = global_dofs(column, triangle)
                    values(entry) = real( &
                        transforms(row, triangle)* &
                        transforms(column, triangle), dp)* &
                        element_dot(row, column)
                end do
            end do
        end do
        call csc_from_triplet( &
            global_dof_count, global_dof_count, rows, columns, values, &
            matrix_dot, status)
    end subroutine assemble_triangle_nedelec_curl_mass_csc_jvp

    subroutine assemble_triangle_nedelec_curl_mass_csc_vjp( &
            mesh, order, quadrature_degree, curl_coefficient, mass_tensor, &
            matrix_values_bar, mesh_vertices_bar, curl_coefficient_bar, &
            mass_tensor_bar, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: order, quadrature_degree
        real(dp), intent(in) :: curl_coefficient, mass_tensor(2, 2)
        real(dp), intent(in) :: matrix_values_bar(:)
        real(dp), intent(out) :: mesh_vertices_bar(:, :)
        real(dp), intent(out) :: curl_coefficient_bar, mass_tensor_bar(2, 2)
        type(fortsparse_status_t), intent(out) :: status

        type(csc_t) :: matrix
        integer, allocatable :: global_dofs(:, :), transforms(:, :)
        real(dp), allocatable :: element_bar(:, :)
        real(dp) :: local_curl_bar, local_tensor_bar(2, 2)
        real(dp) :: local_vertices_bar(2, 3), vertices(2, 3)
        integer :: column, global_dof_count, local_dof_count, local_status
        integer :: node, row, triangle

        mesh_vertices_bar = 0.0_dp
        curl_coefficient_bar = 0.0_dp
        mass_tensor_bar = 0.0_dp
        call assemble_triangle_nedelec_curl_mass_csc( &
            mesh, order, quadrature_degree, matrix, status, curl_coefficient, &
            mass_tensor=mass_tensor)
        if (status%code /= 0) return
        if (any(shape(mesh_vertices_bar) /= shape(mesh%vertices)) .or. &
            size(matrix_values_bar) /= matrix%nnz) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Triangle Nedelec curl-mass VJP shapes differ")
            return
        end if
        call build_triangle_trimmed_dof_map( &
            mesh, order, global_dofs, transforms, global_dof_count, &
            local_status)
        if (local_status /= 0) return
        local_dof_count = size(global_dofs, 1)
        allocate(element_bar(local_dof_count, local_dof_count))
        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            do column = 1, local_dof_count
                do row = 1, local_dof_count
                    element_bar(row, column) = real( &
                        transforms(row, triangle)* &
                        transforms(column, triangle), dp)*csc_bar_at( &
                        matrix, matrix_values_bar, &
                        global_dofs(row, triangle), &
                        global_dofs(column, triangle))
                end do
            end do
            call assemble_triangle_nedelec_curl_mass_element_vjp( &
                vertices, order, quadrature_degree, curl_coefficient, &
                mass_tensor, element_bar, local_vertices_bar, local_curl_bar, &
                local_tensor_bar, local_status)
            if (local_status /= 0) return
            do node = 1, 3
                mesh_vertices_bar(:, mesh%triangles(node, triangle)) = &
                    mesh_vertices_bar(:, mesh%triangles(node, triangle)) + &
                    local_vertices_bar(:, node)
            end do
            curl_coefficient_bar = curl_coefficient_bar + local_curl_bar
            mass_tensor_bar = mass_tensor_bar + local_tensor_bar
        end do
        call status_set(status, 0, "")
    end subroutine assemble_triangle_nedelec_curl_mass_csc_vjp

    subroutine assemble_triangle_nedelec_cell_tensor_csc( &
            mesh, order, quadrature_degree, curl_values, mass_tensors, &
            matrix, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: order, quadrature_degree
        real(dp), intent(in) :: curl_values(:), mass_tensors(:, :, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), global_dofs(:, :), rows(:)
        integer, allocatable :: transforms(:, :)
        real(dp), allocatable :: element_matrix(:, :), values(:)
        real(dp) :: vertices(2, 3)
        integer :: column, entry, global_dof_count, local_dof_count
        integer :: local_status, row, triangle

        if (order < 1 .or. quadrature_degree < 0 .or. &
            size(curl_values) /= mesh%n_triangles .or. &
            size(mass_tensors, 1) /= 2 .or. &
            size(mass_tensors, 2) /= 2 .or. &
            size(mass_tensors, 3) /= mesh%n_triangles) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Nedelec tensor assembly requires one 2-by-2 tensor per cell")
            return
        end if
        call build_triangle_trimmed_dof_map( &
            mesh, order, global_dofs, transforms, global_dof_count, &
            local_status)
        if (local_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Nedelec tensor assembly requires a valid triangle mesh")
            return
        end if

        local_dof_count = size(global_dofs, 1)
        allocate(rows(mesh%n_triangles * local_dof_count**2))
        allocate(columns(mesh%n_triangles * local_dof_count**2))
        allocate(values(mesh%n_triangles * local_dof_count**2))
        entry = 0
        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            call assemble_triangle_nedelec_curl_mass_element( &
                vertices, order, quadrature_degree, element_matrix, &
                local_status, curl_values(triangle), &
                mass_tensor=mass_tensors(:, :, triangle))
            if (local_status /= 0) then
                call status_set( &
                    status, FORTSPARSE_INVALID_MATRIX, &
                    "Nedelec tensor element assembly failed")
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
    end subroutine assemble_triangle_nedelec_cell_tensor_csc

    subroutine assemble_triangle_nedelec_curl_csc( &
            mesh, order, quadrature_degree, matrix, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: order, quadrature_degree
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        type(triangle_lagrange_basis_t) :: scalar_basis
        type(triangle_nedelec_first_kind_t) :: nedelec_basis
        integer, allocatable :: columns(:), nedelec_global_dofs(:, :), rows(:)
        integer, allocatable :: nedelec_transforms(:, :)
        integer, allocatable :: scalar_global_dofs(:, :)
        real(dp), allocatable :: eta(:), local_matrix(:, :)
        real(dp), allocatable :: physical_curls(:), physical_values(:, :)
        real(dp), allocatable :: reference_curls(:), reference_values(:, :)
        real(dp), allocatable :: scalar_gradients(:, :), scalar_values(:)
        real(dp), allocatable :: values(:), weights(:), xi(:)
        real(dp) :: determinant, jacobian(2, 2), vertices(2, 3)
        integer :: entry, local_status, nedelec_dof, nedelec_dof_count
        integer :: nedelec_local_dof_count, point, scalar_dof
        integer :: scalar_dof_count, scalar_local_dof_count, triangle

        if (order < 1 .or. quadrature_degree < 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Nedelec-DG curl assembly requires positive order")
            return
        end if
        call build_triangle_trimmed_dof_map( &
            mesh, order, nedelec_global_dofs, nedelec_transforms, &
            nedelec_dof_count, local_status)
        if (local_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Nedelec-DG curl assembly requires a valid triangle mesh")
            return
        end if
        call build_triangle_discontinuous_dof_map( &
            mesh, order - 1, scalar_global_dofs, scalar_dof_count, &
            local_status)
        if (local_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Nedelec-DG curl assembly could not map scalar degrees")
            return
        end if
        call initialize_triangle_nedelec_first_kind( &
            order, nedelec_basis, local_status)
        if (local_status /= 0) return
        call initialize_triangle_lagrange_basis( &
            order - 1, scalar_basis, local_status)
        if (local_status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, local_status)
        if (local_status /= 0) return

        nedelec_local_dof_count = size(nedelec_global_dofs, 1)
        scalar_local_dof_count = size(scalar_global_dofs, 1)
        allocate(reference_values(2, nedelec_local_dof_count))
        allocate(reference_curls(nedelec_local_dof_count))
        allocate(physical_values(2, nedelec_local_dof_count))
        allocate(physical_curls(nedelec_local_dof_count))
        allocate(scalar_values(scalar_local_dof_count))
        allocate(scalar_gradients(2, scalar_local_dof_count))
        allocate(local_matrix( &
            scalar_local_dof_count, nedelec_local_dof_count))
        allocate(rows( &
            mesh%n_triangles * nedelec_local_dof_count * &
            scalar_local_dof_count))
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
                    "Nedelec-DG curl assembly requires valid CCW triangles")
                return
            end if
            local_matrix = 0.0_dp
            do point = 1, size(weights)
                call evaluate_triangle_nedelec_first_kind( &
                    nedelec_basis, xi(point), eta(point), reference_values, &
                    reference_curls, local_status)
                if (local_status /= 0) return
                call map_triangle_nedelec_covariant( &
                    jacobian, reference_values, reference_curls, &
                    physical_values, physical_curls, local_status)
                if (local_status /= 0) return
                call evaluate_triangle_lagrange_basis( &
                    scalar_basis, xi(point), eta(point), scalar_values, &
                    scalar_gradients, local_status)
                if (local_status /= 0) return
                do nedelec_dof = 1, nedelec_local_dof_count
                    do scalar_dof = 1, scalar_local_dof_count
                        local_matrix(scalar_dof, nedelec_dof) = &
                            local_matrix(scalar_dof, nedelec_dof) + &
                            determinant * weights(point) * &
                            scalar_values(scalar_dof) * &
                            physical_curls(nedelec_dof)
                    end do
                end do
            end do
            do nedelec_dof = 1, nedelec_local_dof_count
                do scalar_dof = 1, scalar_local_dof_count
                    entry = entry + 1
                    rows(entry) = scalar_global_dofs(scalar_dof, triangle)
                    columns(entry) = &
                        nedelec_global_dofs(nedelec_dof, triangle)
                    values(entry) = real( &
                        nedelec_transforms(nedelec_dof, triangle), dp) * &
                        local_matrix(scalar_dof, nedelec_dof)
                end do
            end do
        end do
        call csc_from_triplet( &
            scalar_dof_count, nedelec_dof_count, rows, columns, values, &
            matrix, status)
    end subroutine assemble_triangle_nedelec_curl_csc

    pure subroutine triangle_jacobian_direction(vertices, jacobian)
        real(dp), intent(in) :: vertices(2, 3)
        real(dp), intent(out) :: jacobian(2, 2)

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
    end subroutine triangle_jacobian_direction

    pure subroutine triangle_jacobian( &
            vertices, jacobian, determinant, status)
        real(dp), intent(in) :: vertices(2, 3)
        real(dp), intent(out) :: jacobian(2, 2), determinant
        integer, intent(out) :: status

        call triangle_jacobian_direction(vertices, jacobian)
        determinant = jacobian(1, 1)*jacobian(2, 2) - &
            jacobian(1, 2)*jacobian(2, 1)
        status = 1
        if (determinant <= 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**2)) return
        status = 0
    end subroutine triangle_jacobian

    pure subroutine triangle_jacobian_vjp(jacobian_bar, vertices_bar)
        real(dp), intent(in) :: jacobian_bar(2, 2)
        real(dp), intent(out) :: vertices_bar(2, 3)

        vertices_bar(:, 1) = -sum(jacobian_bar, dim=2)
        vertices_bar(:, 2) = jacobian_bar(:, 1)
        vertices_bar(:, 3) = jacobian_bar(:, 2)
    end subroutine triangle_jacobian_vjp

    pure subroutine triangle_jacobian_vjp_add(jacobian_bar, vertices_bar)
        real(dp), intent(in) :: jacobian_bar(2, 2)
        real(dp), intent(inout) :: vertices_bar(2, 3)

        vertices_bar(:, 1) = vertices_bar(:, 1) - &
            sum(jacobian_bar, dim=2)
        vertices_bar(:, 2) = vertices_bar(:, 2) + jacobian_bar(:, 1)
        vertices_bar(:, 3) = vertices_bar(:, 3) + jacobian_bar(:, 2)
    end subroutine triangle_jacobian_vjp_add

    subroutine initialize_sampled_nedelec_load( &
            mesh, order, quadrature_degree, source_values, basis, &
            global_dofs, transforms, global_count, xi, eta, weights, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: order, quadrature_degree
        real(dp), intent(in) :: source_values(:, :, :)
        type(triangle_nedelec_first_kind_t), intent(out) :: basis
        integer, allocatable, intent(out) :: global_dofs(:, :), transforms(:, :)
        integer, intent(out) :: global_count
        real(dp), allocatable, intent(out) :: xi(:), eta(:), weights(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: local_status

        global_count = 0
        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sampled triangle Nedelec vector load failed")
        if (order < 1 .or. quadrature_degree < 0) return
        call build_triangle_trimmed_dof_map( &
            mesh, order, global_dofs, transforms, global_count, local_status)
        if (local_status /= 0) return
        call initialize_triangle_nedelec_first_kind( &
            order, basis, local_status)
        if (local_status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, local_status)
        if (local_status /= 0) return
        if (.not. valid_vector_samples( &
            source_values, size(weights), mesh%n_triangles)) return
        call status_set(status, 0, "")
    end subroutine initialize_sampled_nedelec_load

    pure logical function valid_vector_samples( &
            source_values, point_count, triangle_count) result(valid)
        real(dp), intent(in) :: source_values(:, :, :)
        integer, intent(in) :: point_count, triangle_count

        valid = size(source_values, 1) == 2 .and. &
            size(source_values, 2) == point_count .and. &
            size(source_values, 3) == triangle_count
    end function valid_vector_samples

    pure logical function valid_vector_sample_gradients( &
            source_values, source_gradients, point_count, triangle_count) &
            result(valid)
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :)
        integer, intent(in) :: point_count, triangle_count

        valid = valid_vector_samples( &
            source_values, point_count, triangle_count)
        if (.not. valid) return
        valid = size(source_gradients, 1) == 2 .and. &
            size(source_gradients, 2) == 2 .and. &
            size(source_gradients, 3) == point_count .and. &
            size(source_gradients, 4) == triangle_count
    end function valid_vector_sample_gradients

    pure logical function valid_vector_sample_products( &
            source_values, source_gradients, source_parameter_dot, &
            point_count, triangle_count) result(valid)
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :)
        real(dp), intent(in) :: source_parameter_dot(:, :, :)
        integer, intent(in) :: point_count, triangle_count

        valid = valid_vector_sample_gradients( &
            source_values, source_gradients, point_count, triangle_count)
        if (.not. valid) return
        valid = all(shape(source_parameter_dot) == shape(source_values))
    end function valid_vector_sample_products

    pure real(dp) function csc_bar_at( &
            matrix, values_bar, row, column) result(value_bar)
        type(csc_t), intent(in) :: matrix
        real(dp), intent(in) :: values_bar(:)
        integer, intent(in) :: row, column

        integer :: entry

        value_bar = 0.0_dp
        do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
            if (matrix%row_idx(entry) == row) then
                value_bar = values_bar(entry)
                return
            end if
        end do
    end function csc_bar_at

end module fortfem_assembly_nedelec_arbitrary_order_2d
