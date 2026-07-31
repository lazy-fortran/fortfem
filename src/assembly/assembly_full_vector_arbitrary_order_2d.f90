module fortfem_assembly_full_vector_arbitrary_order_2d
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_triangle_bdm_arbitrary_order, only: &
        evaluate_triangle_bdm, initialize_triangle_bdm, triangle_bdm_basis_t
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_triangle_global_dof_map, only: &
        build_triangle_full_vector_dof_map
    use fortfem_triangle_nedelec_second_kind, only: &
        evaluate_triangle_nedelec_second_kind, &
        initialize_triangle_nedelec_second_kind, &
        triangle_nedelec_second_kind_t
    use fortfem_triangle_piola_maps, only: &
        map_triangle_nedelec_covariant, &
        map_triangle_nedelec_covariant_jvp, &
        map_triangle_nedelec_covariant_vjp, &
        map_triangle_rt_contravariant, map_triangle_rt_contravariant_jvp, &
        map_triangle_rt_contravariant_vjp
    use fortnum_linalg, only: det2_jvp, det2_vjp
    use fortsparse, only: &
        csc_from_triplet, csc_t, FORTSPARSE_INVALID_MATRIX, &
        fortsparse_status_t, status_set
    implicit none

    private

    public :: assemble_triangle_bdm_div_mass_csc
    public :: assemble_triangle_bdm_div_mass_csc_jvp
    public :: assemble_triangle_bdm_div_mass_csc_vjp
    public :: assemble_triangle_bdm_div_mass_element
    public :: assemble_triangle_bdm_div_mass_element_jvp
    public :: assemble_triangle_bdm_div_mass_element_vjp
    public :: assemble_triangle_bdm_cell_vector_load
    public :: assemble_triangle_bdm_vector_load_samples
    public :: assemble_triangle_bdm_vector_load_samples_jvp
    public :: assemble_triangle_bdm_vector_load_samples_vjp
    public :: assemble_triangle_nedelec_second_curl_mass_csc
    public :: assemble_triangle_nedelec_second_curl_mass_csc_jvp
    public :: assemble_triangle_nedelec_second_curl_mass_csc_vjp
    public :: assemble_triangle_nedelec_second_curl_mass_element
    public :: assemble_triangle_nedelec_second_curl_mass_element_jvp
    public :: assemble_triangle_nedelec_second_curl_mass_element_vjp
    public :: assemble_triangle_nedelec_second_cell_vector_load
    public :: assemble_triangle_nedelec_second_vector_load_samples
    public :: assemble_triangle_nedelec_second_vector_load_samples_jvp
    public :: assemble_triangle_nedelec_second_vector_load_samples_vjp

contains

    subroutine assemble_triangle_nedelec_second_cell_vector_load( &
            mesh, degree, quadrature_degree, cell_values, vector, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: cell_values(:, :)
        real(dp), intent(out) :: vector(:)
        type(fortsparse_status_t), intent(out) :: status

        type(triangle_nedelec_second_kind_t) :: basis
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
            "Second-kind Nedelec cell load assembly failed")
        if (degree < 1 .or. quadrature_degree < 0) return
        if (size(cell_values, 1) /= 2 .or. &
            size(cell_values, 2) /= mesh%n_triangles) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Second-kind cell load requires one vector per triangle")
            return
        end if
        call build_triangle_full_vector_dof_map( &
            mesh, degree, global_dofs, transforms, global_dof_count, &
            local_status)
        if (local_status /= 0 .or. size(vector) /= global_dof_count) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Second-kind cell load requires a valid output space")
            return
        end if
        call initialize_triangle_nedelec_second_kind( &
            degree, basis, local_status)
        if (local_status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, local_status)
        if (local_status /= 0) return

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
                max(1.0_dp, maxval(abs(jacobian))**2)) return
            local_vector = 0.0_dp
            do point = 1, size(weights)
                call evaluate_triangle_nedelec_second_kind( &
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
    end subroutine assemble_triangle_nedelec_second_cell_vector_load

    subroutine assemble_triangle_bdm_cell_vector_load( &
            mesh, degree, quadrature_degree, cell_values, vector, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: cell_values(:, :)
        real(dp), intent(out) :: vector(:)
        type(fortsparse_status_t), intent(out) :: status

        type(triangle_bdm_basis_t) :: basis
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
            status, FORTSPARSE_INVALID_MATRIX, &
            "BDM cell load assembly failed")
        if (degree < 1 .or. quadrature_degree < 0) return
        if (size(cell_values, 1) /= 2 .or. &
            size(cell_values, 2) /= mesh%n_triangles) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "BDM cell load requires one vector per triangle")
            return
        end if
        call build_triangle_full_vector_dof_map( &
            mesh, degree, global_dofs, transforms, global_dof_count, &
            local_status)
        if (local_status /= 0 .or. size(vector) /= global_dof_count) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "BDM cell load requires a valid output space")
            return
        end if
        call initialize_triangle_bdm(degree, basis, local_status)
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
                call evaluate_triangle_bdm( &
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
    end subroutine assemble_triangle_bdm_cell_vector_load

    subroutine assemble_triangle_bdm_vector_load_samples( &
            mesh, degree, quadrature_degree, source_values, vector, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), allocatable, intent(out) :: vector(:)
        type(fortsparse_status_t), intent(out) :: status

        call assemble_triangle_full_vector_load_samples( &
            mesh, degree, quadrature_degree, .true., source_values, vector, &
            status)
    end subroutine assemble_triangle_bdm_vector_load_samples

    subroutine assemble_triangle_nedelec_second_vector_load_samples( &
            mesh, degree, quadrature_degree, source_values, vector, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), allocatable, intent(out) :: vector(:)
        type(fortsparse_status_t), intent(out) :: status

        call assemble_triangle_full_vector_load_samples( &
            mesh, degree, quadrature_degree, .false., source_values, vector, &
            status)
    end subroutine assemble_triangle_nedelec_second_vector_load_samples

    subroutine assemble_triangle_bdm_vector_load_samples_jvp( &
            mesh, degree, quadrature_degree, source_values, source_gradients, &
            mesh_vertices_dot, source_parameter_dot, vector_dot, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :)
        real(dp), intent(in) :: mesh_vertices_dot(:, :)
        real(dp), intent(in) :: source_parameter_dot(:, :, :)
        real(dp), allocatable, intent(out) :: vector_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        call assemble_triangle_full_vector_load_samples_jvp( &
            mesh, degree, quadrature_degree, .true., source_values, &
            source_gradients, mesh_vertices_dot, source_parameter_dot, &
            vector_dot, status)
    end subroutine assemble_triangle_bdm_vector_load_samples_jvp

    subroutine assemble_triangle_nedelec_second_vector_load_samples_jvp( &
            mesh, degree, quadrature_degree, source_values, source_gradients, &
            mesh_vertices_dot, source_parameter_dot, vector_dot, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :)
        real(dp), intent(in) :: mesh_vertices_dot(:, :)
        real(dp), intent(in) :: source_parameter_dot(:, :, :)
        real(dp), allocatable, intent(out) :: vector_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        call assemble_triangle_full_vector_load_samples_jvp( &
            mesh, degree, quadrature_degree, .false., source_values, &
            source_gradients, mesh_vertices_dot, source_parameter_dot, &
            vector_dot, status)
    end subroutine assemble_triangle_nedelec_second_vector_load_samples_jvp

    subroutine assemble_triangle_bdm_vector_load_samples_vjp( &
            mesh, degree, quadrature_degree, source_values, source_gradients, &
            vector_bar, mesh_vertices_bar, source_values_bar, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :), vector_bar(:)
        real(dp), intent(out) :: mesh_vertices_bar(:, :)
        real(dp), intent(out) :: source_values_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        call assemble_triangle_full_vector_load_samples_vjp( &
            mesh, degree, quadrature_degree, .true., source_values, &
            source_gradients, vector_bar, mesh_vertices_bar, &
            source_values_bar, status)
    end subroutine assemble_triangle_bdm_vector_load_samples_vjp

    subroutine assemble_triangle_nedelec_second_vector_load_samples_vjp( &
            mesh, degree, quadrature_degree, source_values, source_gradients, &
            vector_bar, mesh_vertices_bar, source_values_bar, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :), vector_bar(:)
        real(dp), intent(out) :: mesh_vertices_bar(:, :)
        real(dp), intent(out) :: source_values_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        call assemble_triangle_full_vector_load_samples_vjp( &
            mesh, degree, quadrature_degree, .false., source_values, &
            source_gradients, vector_bar, mesh_vertices_bar, &
            source_values_bar, status)
    end subroutine assemble_triangle_nedelec_second_vector_load_samples_vjp

    subroutine assemble_triangle_full_vector_load_samples( &
            mesh, degree, quadrature_degree, normal_family, source_values, &
            vector, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), allocatable, intent(out) :: vector(:)
        type(fortsparse_status_t), intent(out) :: status

        type(triangle_bdm_basis_t) :: bdm_basis
        type(triangle_nedelec_second_kind_t) :: nedelec_basis
        integer, allocatable :: global_dofs(:, :), transforms(:, :)
        real(dp), allocatable :: derivatives(:), eta(:), local_vector(:)
        real(dp), allocatable :: ref_derivatives(:), ref_values(:, :)
        real(dp), allocatable :: values(:, :), weights(:), xi(:)
        real(dp) :: determinant, jacobian(2, 2), vertices(2, 3)
        integer :: dof, dof_count, global_count, local_status, point, triangle

        call initialize_sampled_full_vector_load( &
            mesh, degree, quadrature_degree, normal_family, source_values, &
            bdm_basis, nedelec_basis, global_dofs, transforms, global_count, &
            xi, eta, weights, status)
        if (status%code /= 0) return
        dof_count = size(global_dofs, 1)
        allocate(vector(global_count), source=0.0_dp)
        allocate(local_vector(dof_count))
        allocate(ref_values(2, dof_count), ref_derivatives(dof_count))
        allocate(values(2, dof_count), derivatives(dof_count))
        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            call triangle_jacobian(vertices, jacobian, determinant, local_status)
            if (local_status /= 0) return
            local_vector = 0.0_dp
            do point = 1, size(weights)
                call evaluate_full_vector_reference( &
                    normal_family, bdm_basis, nedelec_basis, xi(point), &
                    eta(point), ref_values, ref_derivatives, local_status)
                if (local_status /= 0) return
                call map_full_vector( &
                    normal_family, jacobian, ref_values, ref_derivatives, &
                    values, derivatives, local_status)
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
    end subroutine assemble_triangle_full_vector_load_samples

    subroutine assemble_triangle_full_vector_load_samples_jvp( &
            mesh, degree, quadrature_degree, normal_family, source_values, &
            source_gradients, mesh_vertices_dot, source_parameter_dot, &
            vector_dot, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :)
        real(dp), intent(in) :: mesh_vertices_dot(:, :)
        real(dp), intent(in) :: source_parameter_dot(:, :, :)
        real(dp), allocatable, intent(out) :: vector_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        type(triangle_bdm_basis_t) :: bdm_basis
        type(triangle_nedelec_second_kind_t) :: nedelec_basis
        integer, allocatable :: global_dofs(:, :), transforms(:, :)
        real(dp), allocatable :: derivatives(:), derivatives_dot(:), eta(:)
        real(dp), allocatable :: local_dot(:), ref_derivatives(:)
        real(dp), allocatable :: ref_values(:, :), values(:, :)
        real(dp), allocatable :: values_dot(:, :), weights(:), xi(:)
        real(dp), allocatable :: zero_derivatives(:), zero_values(:, :)
        real(dp) :: determinant, determinant_dot, jacobian(2, 2)
        real(dp) :: jacobian_dot(2, 2), point_dot(2), reference_point(2)
        real(dp) :: source_dot(2), vertices(2, 3), vertices_dot(2, 3)
        integer :: dof, dof_count, global_count, local_status, point, triangle

        call initialize_sampled_full_vector_load( &
            mesh, degree, quadrature_degree, normal_family, source_values, &
            bdm_basis, nedelec_basis, global_dofs, transforms, global_count, &
            xi, eta, weights, status)
        if (status%code /= 0) return
        if (any(shape(mesh_vertices_dot) /= shape(mesh%vertices))) return
        if (.not. valid_vector_sample_products( &
            source_values, source_gradients, source_parameter_dot, &
            size(weights), mesh%n_triangles)) return
        dof_count = size(global_dofs, 1)
        allocate(vector_dot(global_count), source=0.0_dp)
        allocate(local_dot(dof_count))
        allocate(ref_values(2, dof_count), ref_derivatives(dof_count))
        allocate(values(2, dof_count), derivatives(dof_count))
        allocate(values_dot(2, dof_count), derivatives_dot(dof_count))
        allocate(zero_values(2, dof_count), zero_derivatives(dof_count), &
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
                call evaluate_full_vector_reference( &
                    normal_family, bdm_basis, nedelec_basis, xi(point), &
                    eta(point), ref_values, ref_derivatives, local_status)
                if (local_status /= 0) return
                call map_full_vector( &
                    normal_family, jacobian, ref_values, ref_derivatives, &
                    values, derivatives, local_status)
                if (local_status /= 0) return
                call map_full_vector_jvp( &
                    normal_family, jacobian, ref_values, ref_derivatives, &
                    jacobian_dot, zero_values, zero_derivatives, values_dot, &
                    derivatives_dot, local_status)
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
    end subroutine assemble_triangle_full_vector_load_samples_jvp

    subroutine assemble_triangle_full_vector_load_samples_vjp( &
            mesh, degree, quadrature_degree, normal_family, source_values, &
            source_gradients, vector_bar, mesh_vertices_bar, source_values_bar, &
            status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :), vector_bar(:)
        real(dp), intent(out) :: mesh_vertices_bar(:, :)
        real(dp), intent(out) :: source_values_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        type(triangle_bdm_basis_t) :: bdm_basis
        type(triangle_nedelec_second_kind_t) :: nedelec_basis
        integer, allocatable :: global_dofs(:, :), transforms(:, :)
        real(dp), allocatable :: derivatives(:), derivatives_bar(:), eta(:)
        real(dp), allocatable :: local_bar(:), ref_derivatives(:)
        real(dp), allocatable :: ref_derivatives_bar(:), ref_values(:, :)
        real(dp), allocatable :: ref_values_bar(:, :), values(:, :)
        real(dp), allocatable :: values_bar(:, :), weights(:), xi(:)
        real(dp) :: determinant, determinant_bar, determinant_jacobian_bar(2, 2)
        real(dp) :: jacobian(2, 2), jacobian_bar(2, 2)
        real(dp) :: local_jacobian_bar(2, 2), local_vertices_bar(2, 3)
        real(dp) :: point_bar(2), reference_point(2), sample_bar(2), seed
        real(dp) :: vertices(2, 3)
        integer :: dof, dof_count, global_count, local_status, node
        integer :: point, triangle

        mesh_vertices_bar = 0.0_dp
        source_values_bar = 0.0_dp
        call initialize_sampled_full_vector_load( &
            mesh, degree, quadrature_degree, normal_family, source_values, &
            bdm_basis, nedelec_basis, global_dofs, transforms, global_count, &
            xi, eta, weights, status)
        if (status%code /= 0) return
        if (any(shape(mesh_vertices_bar) /= shape(mesh%vertices))) return
        if (any(shape(source_values_bar) /= shape(source_values))) return
        if (.not. valid_vector_sample_gradients( &
            source_values, source_gradients, size(weights), &
            mesh%n_triangles)) return
        if (size(vector_bar) /= global_count) return
        dof_count = size(global_dofs, 1)
        allocate(local_bar(dof_count))
        allocate(ref_values(2, dof_count), ref_derivatives(dof_count))
        allocate(ref_values_bar(2, dof_count))
        allocate(ref_derivatives_bar(dof_count))
        allocate(values(2, dof_count), derivatives(dof_count))
        allocate(values_bar(2, dof_count), derivatives_bar(dof_count))
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
                call evaluate_full_vector_reference( &
                    normal_family, bdm_basis, nedelec_basis, xi(point), &
                    eta(point), ref_values, ref_derivatives, local_status)
                if (local_status /= 0) return
                call map_full_vector( &
                    normal_family, jacobian, ref_values, ref_derivatives, &
                    values, derivatives, local_status)
                if (local_status /= 0) return
                values_bar = 0.0_dp
                derivatives_bar = 0.0_dp
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
                call map_full_vector_vjp( &
                    normal_family, jacobian, ref_values, ref_derivatives, &
                    values_bar, derivatives_bar, local_jacobian_bar, &
                    ref_values_bar, ref_derivatives_bar, local_status)
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
    end subroutine assemble_triangle_full_vector_load_samples_vjp

    subroutine assemble_triangle_nedelec_second_curl_mass_element( &
            vertices, degree, quadrature_degree, matrix, status, &
            curl_coefficient, mass_coefficient)
        real(dp), intent(in) :: vertices(2, 3)
        integer, intent(in) :: degree, quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: curl_coefficient, mass_coefficient

        real(dp) :: curl_weight, mass_weight

        curl_weight = 1.0_dp
        mass_weight = 1.0_dp
        if (present(curl_coefficient)) curl_weight = curl_coefficient
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        call assemble_triangle_full_vector_element( &
            vertices, degree, quadrature_degree, .false., curl_weight, &
            mass_weight, matrix, status)
    end subroutine assemble_triangle_nedelec_second_curl_mass_element

    subroutine assemble_triangle_bdm_div_mass_element( &
            vertices, degree, quadrature_degree, matrix, status, &
            divergence_coefficient, mass_coefficient)
        real(dp), intent(in) :: vertices(2, 3)
        integer, intent(in) :: degree, quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: divergence_coefficient
        real(dp), intent(in), optional :: mass_coefficient

        real(dp) :: divergence_weight, mass_weight

        divergence_weight = 1.0_dp
        mass_weight = 1.0_dp
        if (present(divergence_coefficient)) then
            divergence_weight = divergence_coefficient
        end if
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        call assemble_triangle_full_vector_element( &
            vertices, degree, quadrature_degree, .true., divergence_weight, &
            mass_weight, matrix, status)
    end subroutine assemble_triangle_bdm_div_mass_element

    subroutine assemble_triangle_nedelec_second_curl_mass_element_jvp( &
            vertices, degree, quadrature_degree, curl_coefficient, &
            mass_coefficient, vertices_dot, curl_coefficient_dot, &
            mass_coefficient_dot, matrix_dot, status)
        real(dp), intent(in) :: vertices(2, 3), vertices_dot(2, 3)
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: curl_coefficient, mass_coefficient
        real(dp), intent(in) :: curl_coefficient_dot, mass_coefficient_dot
        real(dp), allocatable, intent(out) :: matrix_dot(:, :)
        integer, intent(out) :: status

        call assemble_triangle_full_vector_element_jvp( &
            vertices, degree, quadrature_degree, .false., curl_coefficient, &
            mass_coefficient, vertices_dot, curl_coefficient_dot, &
            mass_coefficient_dot, matrix_dot, status)
    end subroutine assemble_triangle_nedelec_second_curl_mass_element_jvp

    subroutine assemble_triangle_bdm_div_mass_element_jvp( &
            vertices, degree, quadrature_degree, divergence_coefficient, &
            mass_coefficient, vertices_dot, divergence_coefficient_dot, &
            mass_coefficient_dot, matrix_dot, status)
        real(dp), intent(in) :: vertices(2, 3), vertices_dot(2, 3)
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: divergence_coefficient, mass_coefficient
        real(dp), intent(in) :: divergence_coefficient_dot
        real(dp), intent(in) :: mass_coefficient_dot
        real(dp), allocatable, intent(out) :: matrix_dot(:, :)
        integer, intent(out) :: status

        call assemble_triangle_full_vector_element_jvp( &
            vertices, degree, quadrature_degree, .true., &
            divergence_coefficient, mass_coefficient, vertices_dot, &
            divergence_coefficient_dot, mass_coefficient_dot, matrix_dot, &
            status)
    end subroutine assemble_triangle_bdm_div_mass_element_jvp

    subroutine assemble_triangle_nedelec_second_curl_mass_element_vjp( &
            vertices, degree, quadrature_degree, curl_coefficient, &
            mass_coefficient, matrix_bar, vertices_bar, curl_coefficient_bar, &
            mass_coefficient_bar, status)
        real(dp), intent(in) :: vertices(2, 3), matrix_bar(:, :)
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: curl_coefficient, mass_coefficient
        real(dp), intent(out) :: vertices_bar(2, 3)
        real(dp), intent(out) :: curl_coefficient_bar, mass_coefficient_bar
        integer, intent(out) :: status

        call assemble_triangle_full_vector_element_vjp( &
            vertices, degree, quadrature_degree, .false., curl_coefficient, &
            mass_coefficient, matrix_bar, vertices_bar, curl_coefficient_bar, &
            mass_coefficient_bar, status)
    end subroutine assemble_triangle_nedelec_second_curl_mass_element_vjp

    subroutine assemble_triangle_bdm_div_mass_element_vjp( &
            vertices, degree, quadrature_degree, divergence_coefficient, &
            mass_coefficient, matrix_bar, vertices_bar, &
            divergence_coefficient_bar, mass_coefficient_bar, status)
        real(dp), intent(in) :: vertices(2, 3), matrix_bar(:, :)
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: divergence_coefficient, mass_coefficient
        real(dp), intent(out) :: vertices_bar(2, 3)
        real(dp), intent(out) :: divergence_coefficient_bar
        real(dp), intent(out) :: mass_coefficient_bar
        integer, intent(out) :: status

        call assemble_triangle_full_vector_element_vjp( &
            vertices, degree, quadrature_degree, .true., &
            divergence_coefficient, mass_coefficient, matrix_bar, &
            vertices_bar, divergence_coefficient_bar, mass_coefficient_bar, &
            status)
    end subroutine assemble_triangle_bdm_div_mass_element_vjp

    subroutine assemble_triangle_nedelec_second_curl_mass_csc( &
            mesh, degree, quadrature_degree, matrix, status, &
            curl_coefficient, mass_coefficient)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: curl_coefficient, mass_coefficient

        real(dp) :: curl_weight, mass_weight

        curl_weight = 1.0_dp
        mass_weight = 1.0_dp
        if (present(curl_coefficient)) curl_weight = curl_coefficient
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        call assemble_triangle_full_vector_csc( &
            mesh, degree, quadrature_degree, .false., curl_weight, &
            mass_weight, matrix, status)
    end subroutine assemble_triangle_nedelec_second_curl_mass_csc

    subroutine assemble_triangle_bdm_div_mass_csc( &
            mesh, degree, quadrature_degree, matrix, status, &
            divergence_coefficient, mass_coefficient)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: divergence_coefficient
        real(dp), intent(in), optional :: mass_coefficient

        real(dp) :: divergence_weight, mass_weight

        divergence_weight = 1.0_dp
        mass_weight = 1.0_dp
        if (present(divergence_coefficient)) then
            divergence_weight = divergence_coefficient
        end if
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        call assemble_triangle_full_vector_csc( &
            mesh, degree, quadrature_degree, .true., divergence_weight, &
            mass_weight, matrix, status)
    end subroutine assemble_triangle_bdm_div_mass_csc

    subroutine assemble_triangle_nedelec_second_curl_mass_csc_jvp( &
            mesh, degree, quadrature_degree, curl_coefficient, &
            mass_coefficient, mesh_vertices_dot, curl_coefficient_dot, &
            mass_coefficient_dot, matrix_dot, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: curl_coefficient, mass_coefficient
        real(dp), intent(in) :: mesh_vertices_dot(:, :)
        real(dp), intent(in) :: curl_coefficient_dot, mass_coefficient_dot
        type(csc_t), intent(out) :: matrix_dot
        type(fortsparse_status_t), intent(out) :: status

        call assemble_triangle_full_vector_csc_jvp( &
            mesh, degree, quadrature_degree, .false., curl_coefficient, &
            mass_coefficient, mesh_vertices_dot, curl_coefficient_dot, &
            mass_coefficient_dot, matrix_dot, status)
    end subroutine assemble_triangle_nedelec_second_curl_mass_csc_jvp

    subroutine assemble_triangle_bdm_div_mass_csc_jvp( &
            mesh, degree, quadrature_degree, divergence_coefficient, &
            mass_coefficient, mesh_vertices_dot, divergence_coefficient_dot, &
            mass_coefficient_dot, matrix_dot, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: divergence_coefficient, mass_coefficient
        real(dp), intent(in) :: mesh_vertices_dot(:, :)
        real(dp), intent(in) :: divergence_coefficient_dot
        real(dp), intent(in) :: mass_coefficient_dot
        type(csc_t), intent(out) :: matrix_dot
        type(fortsparse_status_t), intent(out) :: status

        call assemble_triangle_full_vector_csc_jvp( &
            mesh, degree, quadrature_degree, .true., divergence_coefficient, &
            mass_coefficient, mesh_vertices_dot, divergence_coefficient_dot, &
            mass_coefficient_dot, matrix_dot, status)
    end subroutine assemble_triangle_bdm_div_mass_csc_jvp

    subroutine assemble_triangle_nedelec_second_curl_mass_csc_vjp( &
            mesh, degree, quadrature_degree, curl_coefficient, &
            mass_coefficient, matrix_values_bar, mesh_vertices_bar, &
            curl_coefficient_bar, mass_coefficient_bar, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: curl_coefficient, mass_coefficient
        real(dp), intent(in) :: matrix_values_bar(:)
        real(dp), intent(out) :: mesh_vertices_bar(:, :)
        real(dp), intent(out) :: curl_coefficient_bar, mass_coefficient_bar
        type(fortsparse_status_t), intent(out) :: status

        call assemble_triangle_full_vector_csc_vjp( &
            mesh, degree, quadrature_degree, .false., curl_coefficient, &
            mass_coefficient, matrix_values_bar, mesh_vertices_bar, &
            curl_coefficient_bar, mass_coefficient_bar, status)
    end subroutine assemble_triangle_nedelec_second_curl_mass_csc_vjp

    subroutine assemble_triangle_bdm_div_mass_csc_vjp( &
            mesh, degree, quadrature_degree, divergence_coefficient, &
            mass_coefficient, matrix_values_bar, mesh_vertices_bar, &
            divergence_coefficient_bar, mass_coefficient_bar, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: divergence_coefficient, mass_coefficient
        real(dp), intent(in) :: matrix_values_bar(:)
        real(dp), intent(out) :: mesh_vertices_bar(:, :)
        real(dp), intent(out) :: divergence_coefficient_bar
        real(dp), intent(out) :: mass_coefficient_bar
        type(fortsparse_status_t), intent(out) :: status

        call assemble_triangle_full_vector_csc_vjp( &
            mesh, degree, quadrature_degree, .true., divergence_coefficient, &
            mass_coefficient, matrix_values_bar, mesh_vertices_bar, &
            divergence_coefficient_bar, mass_coefficient_bar, status)
    end subroutine assemble_triangle_bdm_div_mass_csc_vjp

    subroutine assemble_triangle_full_vector_csc( &
            mesh, degree, quadrature_degree, normal_family, &
            derivative_coefficient, mass_coefficient, matrix, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: derivative_coefficient, mass_coefficient
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), global_dofs(:, :), rows(:)
        integer, allocatable :: transforms(:, :)
        real(dp), allocatable :: element_matrix(:, :), values(:)
        real(dp) :: vertices(2, 3)
        integer :: column, entry, global_dof_count, local_dof_count
        integer :: local_status, row, triangle

        if (degree < 1 .or. quadrature_degree < 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Full vector sparse assembly requires positive degree")
            return
        end if
        call build_triangle_full_vector_dof_map( &
            mesh, degree, global_dofs, transforms, global_dof_count, &
            local_status)
        if (local_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Full vector sparse assembly requires a valid triangle mesh")
            return
        end if

        local_dof_count = size(global_dofs, 1)
        allocate(rows(mesh%n_triangles * local_dof_count**2))
        allocate(columns(mesh%n_triangles * local_dof_count**2))
        allocate(values(mesh%n_triangles * local_dof_count**2))
        entry = 0
        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            call assemble_triangle_full_vector_element( &
                vertices, degree, quadrature_degree, normal_family, &
                derivative_coefficient, mass_coefficient, element_matrix, &
                local_status)
            if (local_status /= 0) then
                call status_set( &
                    status, FORTSPARSE_INVALID_MATRIX, &
                    "Full vector assembly requires valid CCW triangles")
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
    end subroutine assemble_triangle_full_vector_csc

    subroutine assemble_triangle_full_vector_csc_jvp( &
            mesh, degree, quadrature_degree, normal_family, &
            derivative_coefficient, mass_coefficient, mesh_vertices_dot, &
            derivative_coefficient_dot, mass_coefficient_dot, matrix_dot, &
            status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: derivative_coefficient, mass_coefficient
        real(dp), intent(in) :: mesh_vertices_dot(:, :)
        real(dp), intent(in) :: derivative_coefficient_dot
        real(dp), intent(in) :: mass_coefficient_dot
        type(csc_t), intent(out) :: matrix_dot
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), global_dofs(:, :), rows(:)
        integer, allocatable :: transforms(:, :)
        real(dp), allocatable :: element_dot(:, :), values(:)
        real(dp) :: vertices(2, 3), vertices_dot(2, 3)
        integer :: column, entry, global_count, local_count, local_status
        integer :: row, triangle

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Full-vector triangle JVP assembly failed")
        if (degree < 1 .or. quadrature_degree < 0) return
        if (any(shape(mesh_vertices_dot) /= shape(mesh%vertices))) return
        call build_triangle_full_vector_dof_map( &
            mesh, degree, global_dofs, transforms, global_count, local_status)
        if (local_status /= 0) return
        local_count = size(global_dofs, 1)
        allocate(rows(mesh%n_triangles*local_count**2))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            vertices_dot = &
                mesh_vertices_dot(:, mesh%triangles(:, triangle))
            call assemble_triangle_full_vector_element_jvp( &
                vertices, degree, quadrature_degree, normal_family, &
                derivative_coefficient, mass_coefficient, vertices_dot, &
                derivative_coefficient_dot, mass_coefficient_dot, element_dot, &
                local_status)
            if (local_status /= 0) return
            do column = 1, local_count
                do row = 1, local_count
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
            global_count, global_count, rows, columns, values, matrix_dot, &
            status)
    end subroutine assemble_triangle_full_vector_csc_jvp

    subroutine assemble_triangle_full_vector_csc_vjp( &
            mesh, degree, quadrature_degree, normal_family, &
            derivative_coefficient, mass_coefficient, matrix_values_bar, &
            mesh_vertices_bar, derivative_coefficient_bar, &
            mass_coefficient_bar, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: derivative_coefficient, mass_coefficient
        real(dp), intent(in) :: matrix_values_bar(:)
        real(dp), intent(out) :: mesh_vertices_bar(:, :)
        real(dp), intent(out) :: derivative_coefficient_bar
        real(dp), intent(out) :: mass_coefficient_bar
        type(fortsparse_status_t), intent(out) :: status

        type(csc_t) :: matrix
        integer, allocatable :: global_dofs(:, :), transforms(:, :)
        real(dp), allocatable :: element_bar(:, :)
        real(dp) :: local_derivative_bar, local_mass_bar
        real(dp) :: local_vertices_bar(2, 3), vertices(2, 3)
        integer :: column, global_count, local_count, local_status
        integer :: node, row, triangle

        mesh_vertices_bar = 0.0_dp
        derivative_coefficient_bar = 0.0_dp
        mass_coefficient_bar = 0.0_dp
        call assemble_triangle_full_vector_csc( &
            mesh, degree, quadrature_degree, normal_family, &
            derivative_coefficient, mass_coefficient, matrix, status)
        if (status%code /= 0) return
        if (any(shape(mesh_vertices_bar) /= shape(mesh%vertices)) .or. &
            size(matrix_values_bar) /= matrix%nnz) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Full-vector triangle VJP shapes differ")
            return
        end if
        call build_triangle_full_vector_dof_map( &
            mesh, degree, global_dofs, transforms, global_count, local_status)
        if (local_status /= 0) return
        local_count = size(global_dofs, 1)
        allocate(element_bar(local_count, local_count))
        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            do column = 1, local_count
                do row = 1, local_count
                    element_bar(row, column) = real( &
                        transforms(row, triangle)* &
                        transforms(column, triangle), dp)*csc_bar_at( &
                        matrix, matrix_values_bar, &
                        global_dofs(row, triangle), &
                        global_dofs(column, triangle))
                end do
            end do
            call assemble_triangle_full_vector_element_vjp( &
                vertices, degree, quadrature_degree, normal_family, &
                derivative_coefficient, mass_coefficient, element_bar, &
                local_vertices_bar, local_derivative_bar, local_mass_bar, &
                local_status)
            if (local_status /= 0) return
            do node = 1, 3
                mesh_vertices_bar(:, mesh%triangles(node, triangle)) = &
                    mesh_vertices_bar(:, mesh%triangles(node, triangle)) + &
                    local_vertices_bar(:, node)
            end do
            derivative_coefficient_bar = &
                derivative_coefficient_bar + local_derivative_bar
            mass_coefficient_bar = mass_coefficient_bar + local_mass_bar
        end do
        call status_set(status, 0, "")
    end subroutine assemble_triangle_full_vector_csc_vjp

    subroutine assemble_triangle_full_vector_element( &
            vertices, degree, quadrature_degree, normal_family, &
            derivative_coefficient, mass_coefficient, matrix, status)
        real(dp), intent(in) :: vertices(2, 3)
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: derivative_coefficient, mass_coefficient
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        type(triangle_bdm_basis_t) :: bdm_basis
        type(triangle_nedelec_second_kind_t) :: nedelec_basis
        real(dp), allocatable :: eta(:), physical_derivatives(:)
        real(dp), allocatable :: physical_values(:, :)
        real(dp), allocatable :: reference_derivatives(:)
        real(dp), allocatable :: reference_values(:, :), weights(:), xi(:)
        real(dp) :: determinant, jacobian(2, 2), physical_weight
        integer :: column, dof_count, point, row

        status = 1
        if (degree < 1 .or. quadrature_degree < 0) return
        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        determinant = jacobian(1, 1) * jacobian(2, 2) - &
            jacobian(1, 2) * jacobian(2, 1)
        if (determinant <= 64.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp, maxval(abs(jacobian))**2)) return

        dof_count = (degree + 1) * (degree + 2)
        if (normal_family) then
            call initialize_triangle_bdm(degree, bdm_basis, status)
        else
            call initialize_triangle_nedelec_second_kind( &
                degree, nedelec_basis, status)
        end if
        if (status /= 0) return
        allocate(matrix(dof_count, dof_count))
        allocate(reference_values(2, dof_count))
        allocate(reference_derivatives(dof_count))
        allocate(physical_values(2, dof_count))
        allocate(physical_derivatives(dof_count))
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return

        matrix = 0.0_dp
        do point = 1, size(weights)
            if (normal_family) then
                call evaluate_triangle_bdm( &
                    bdm_basis, xi(point), eta(point), reference_values, &
                    reference_derivatives, status)
                if (status /= 0) return
                call map_triangle_rt_contravariant( &
                    jacobian, reference_values, reference_derivatives, &
                    physical_values, physical_derivatives, status)
            else
                call evaluate_triangle_nedelec_second_kind( &
                    nedelec_basis, xi(point), eta(point), reference_values, &
                    reference_derivatives, status)
                if (status /= 0) return
                call map_triangle_nedelec_covariant( &
                    jacobian, reference_values, reference_derivatives, &
                    physical_values, physical_derivatives, status)
            end if
            if (status /= 0) return
            physical_weight = determinant * weights(point)
            do column = 1, dof_count
                do row = 1, dof_count
                    matrix(row, column) = matrix(row, column) + &
                        physical_weight * ( &
                        derivative_coefficient * physical_derivatives(row) * &
                        physical_derivatives(column) + &
                        mass_coefficient * dot_product( &
                        physical_values(:, row), physical_values(:, column)))
                end do
            end do
        end do
        status = 0
    end subroutine assemble_triangle_full_vector_element

    subroutine assemble_triangle_full_vector_element_jvp( &
            vertices, degree, quadrature_degree, normal_family, &
            derivative_coefficient, mass_coefficient, vertices_dot, &
            derivative_coefficient_dot, mass_coefficient_dot, matrix_dot, &
            status)
        real(dp), intent(in) :: vertices(2, 3), vertices_dot(2, 3)
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: derivative_coefficient, mass_coefficient
        real(dp), intent(in) :: derivative_coefficient_dot
        real(dp), intent(in) :: mass_coefficient_dot
        real(dp), allocatable, intent(out) :: matrix_dot(:, :)
        integer, intent(out) :: status

        type(triangle_bdm_basis_t) :: bdm_basis
        type(triangle_nedelec_second_kind_t) :: nedelec_basis
        real(dp), allocatable :: derivatives(:), derivatives_dot(:), eta(:)
        real(dp), allocatable :: ref_derivatives(:), ref_values(:, :)
        real(dp), allocatable :: values(:, :), values_dot(:, :)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp), allocatable :: zero_derivatives(:), zero_values(:, :)
        real(dp) :: determinant, determinant_dot, derivative_energy
        real(dp) :: jacobian(2, 2), jacobian_dot(2, 2), mass_energy
        integer :: column, dof_count, point, row

        status = 1
        if (degree < 1 .or. quadrature_degree < 0) return
        call triangle_jacobian(vertices, jacobian, determinant, status)
        if (status /= 0) return
        call triangle_jacobian_direction(vertices_dot, jacobian_dot)
        call det2_jvp(jacobian, jacobian_dot, determinant_dot)
        dof_count = (degree + 1)*(degree + 2)
        if (normal_family) then
            call initialize_triangle_bdm(degree, bdm_basis, status)
        else
            call initialize_triangle_nedelec_second_kind( &
                degree, nedelec_basis, status)
        end if
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate(matrix_dot(dof_count, dof_count), source=0.0_dp)
        allocate(ref_values(2, dof_count), ref_derivatives(dof_count))
        allocate(zero_values(2, dof_count), zero_derivatives(dof_count), &
            source=0.0_dp)
        allocate(values(2, dof_count), derivatives(dof_count))
        allocate(values_dot(2, dof_count), derivatives_dot(dof_count))
        do point = 1, size(weights)
            if (normal_family) then
                call evaluate_triangle_bdm( &
                    bdm_basis, xi(point), eta(point), ref_values, &
                    ref_derivatives, status)
                if (status /= 0) return
                call map_triangle_rt_contravariant( &
                    jacobian, ref_values, ref_derivatives, values, &
                    derivatives, status)
                if (status /= 0) return
                call map_triangle_rt_contravariant_jvp( &
                    jacobian, ref_values, ref_derivatives, jacobian_dot, &
                    zero_values, zero_derivatives, values_dot, derivatives_dot, &
                    status)
            else
                call evaluate_triangle_nedelec_second_kind( &
                    nedelec_basis, xi(point), eta(point), ref_values, &
                    ref_derivatives, status)
                if (status /= 0) return
                call map_triangle_nedelec_covariant( &
                    jacobian, ref_values, ref_derivatives, values, &
                    derivatives, status)
                if (status /= 0) return
                call map_triangle_nedelec_covariant_jvp( &
                    jacobian, ref_values, ref_derivatives, jacobian_dot, &
                    zero_values, zero_derivatives, values_dot, derivatives_dot, &
                    status)
            end if
            if (status /= 0) return
            do column = 1, dof_count
                do row = 1, dof_count
                    derivative_energy = &
                        derivatives(row)*derivatives(column)
                    mass_energy = dot_product(values(:, row), values(:, column))
                    matrix_dot(row, column) = matrix_dot(row, column) + &
                        weights(point)*(determinant_dot*( &
                        derivative_coefficient*derivative_energy + &
                        mass_coefficient*mass_energy) + determinant*( &
                        derivative_coefficient_dot*derivative_energy + &
                        derivative_coefficient*( &
                        derivatives_dot(row)*derivatives(column) + &
                        derivatives(row)*derivatives_dot(column)) + &
                        mass_coefficient_dot*mass_energy + mass_coefficient*( &
                        dot_product(values_dot(:, row), values(:, column)) + &
                        dot_product(values(:, row), values_dot(:, column)))))
                end do
            end do
        end do
        status = 0
    end subroutine assemble_triangle_full_vector_element_jvp

    subroutine assemble_triangle_full_vector_element_vjp( &
            vertices, degree, quadrature_degree, normal_family, &
            derivative_coefficient, mass_coefficient, matrix_bar, &
            vertices_bar, derivative_coefficient_bar, mass_coefficient_bar, &
            status)
        real(dp), intent(in) :: vertices(2, 3), matrix_bar(:, :)
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: derivative_coefficient, mass_coefficient
        real(dp), intent(out) :: vertices_bar(2, 3)
        real(dp), intent(out) :: derivative_coefficient_bar
        real(dp), intent(out) :: mass_coefficient_bar
        integer, intent(out) :: status

        type(triangle_bdm_basis_t) :: bdm_basis
        type(triangle_nedelec_second_kind_t) :: nedelec_basis
        real(dp), allocatable :: derivatives(:), derivatives_bar(:), eta(:)
        real(dp), allocatable :: ref_derivatives(:), ref_derivatives_bar(:)
        real(dp), allocatable :: ref_values(:, :), ref_values_bar(:, :)
        real(dp), allocatable :: values(:, :), values_bar(:, :)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp) :: determinant, determinant_bar, determinant_jacobian_bar(2, 2)
        real(dp) :: derivative_energy, jacobian(2, 2), jacobian_bar(2, 2)
        real(dp) :: local_jacobian_bar(2, 2), mass_energy, seed
        integer :: column, dof_count, point, row

        vertices_bar = 0.0_dp
        derivative_coefficient_bar = 0.0_dp
        mass_coefficient_bar = 0.0_dp
        status = 1
        if (degree < 1 .or. quadrature_degree < 0) return
        dof_count = (degree + 1)*(degree + 2)
        if (size(matrix_bar, 1) /= dof_count .or. &
            size(matrix_bar, 2) /= dof_count) return
        call triangle_jacobian(vertices, jacobian, determinant, status)
        if (status /= 0) return
        if (normal_family) then
            call initialize_triangle_bdm(degree, bdm_basis, status)
        else
            call initialize_triangle_nedelec_second_kind( &
                degree, nedelec_basis, status)
        end if
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate(ref_values(2, dof_count), ref_derivatives(dof_count))
        allocate(ref_values_bar(2, dof_count))
        allocate(ref_derivatives_bar(dof_count))
        allocate(values(2, dof_count), derivatives(dof_count))
        allocate(values_bar(2, dof_count), derivatives_bar(dof_count))
        jacobian_bar = 0.0_dp
        determinant_bar = 0.0_dp
        do point = 1, size(weights)
            if (normal_family) then
                call evaluate_triangle_bdm( &
                    bdm_basis, xi(point), eta(point), ref_values, &
                    ref_derivatives, status)
                if (status /= 0) return
                call map_triangle_rt_contravariant( &
                    jacobian, ref_values, ref_derivatives, values, &
                    derivatives, status)
            else
                call evaluate_triangle_nedelec_second_kind( &
                    nedelec_basis, xi(point), eta(point), ref_values, &
                    ref_derivatives, status)
                if (status /= 0) return
                call map_triangle_nedelec_covariant( &
                    jacobian, ref_values, ref_derivatives, values, &
                    derivatives, status)
            end if
            if (status /= 0) return
            values_bar = 0.0_dp
            derivatives_bar = 0.0_dp
            do column = 1, dof_count
                do row = 1, dof_count
                    seed = weights(point)*matrix_bar(row, column)
                    derivative_energy = &
                        derivatives(row)*derivatives(column)
                    mass_energy = dot_product(values(:, row), values(:, column))
                    determinant_bar = determinant_bar + seed*( &
                        derivative_coefficient*derivative_energy + &
                        mass_coefficient*mass_energy)
                    derivative_coefficient_bar = &
                        derivative_coefficient_bar + &
                        seed*determinant*derivative_energy
                    mass_coefficient_bar = mass_coefficient_bar + &
                        seed*determinant*mass_energy
                    derivatives_bar(row) = derivatives_bar(row) + &
                        seed*determinant*derivative_coefficient* &
                        derivatives(column)
                    derivatives_bar(column) = derivatives_bar(column) + &
                        seed*determinant*derivative_coefficient*derivatives(row)
                    values_bar(:, row) = values_bar(:, row) + &
                        seed*determinant*mass_coefficient*values(:, column)
                    values_bar(:, column) = values_bar(:, column) + &
                        seed*determinant*mass_coefficient*values(:, row)
                end do
            end do
            if (normal_family) then
                call map_triangle_rt_contravariant_vjp( &
                    jacobian, ref_values, ref_derivatives, values_bar, &
                    derivatives_bar, local_jacobian_bar, ref_values_bar, &
                    ref_derivatives_bar, status)
            else
                call map_triangle_nedelec_covariant_vjp( &
                    jacobian, ref_values, ref_derivatives, values_bar, &
                    derivatives_bar, local_jacobian_bar, ref_values_bar, &
                    ref_derivatives_bar, status)
            end if
            if (status /= 0) return
            jacobian_bar = jacobian_bar + local_jacobian_bar
        end do
        call det2_vjp(jacobian, determinant_bar, determinant_jacobian_bar)
        call triangle_jacobian_vjp( &
            jacobian_bar + determinant_jacobian_bar, vertices_bar)
        status = 0
    end subroutine assemble_triangle_full_vector_element_vjp

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

    subroutine initialize_sampled_full_vector_load( &
            mesh, degree, quadrature_degree, normal_family, source_values, &
            bdm_basis, nedelec_basis, global_dofs, transforms, global_count, &
            xi, eta, weights, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: source_values(:, :, :)
        type(triangle_bdm_basis_t), intent(out) :: bdm_basis
        type(triangle_nedelec_second_kind_t), intent(out) :: nedelec_basis
        integer, allocatable, intent(out) :: global_dofs(:, :), transforms(:, :)
        integer, intent(out) :: global_count
        real(dp), allocatable, intent(out) :: xi(:), eta(:), weights(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: local_status

        global_count = 0
        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sampled full-vector triangle load failed")
        if (degree < 1 .or. quadrature_degree < 0) return
        call build_triangle_full_vector_dof_map( &
            mesh, degree, global_dofs, transforms, global_count, local_status)
        if (local_status /= 0) return
        if (normal_family) then
            call initialize_triangle_bdm(degree, bdm_basis, local_status)
        else
            call initialize_triangle_nedelec_second_kind( &
                degree, nedelec_basis, local_status)
        end if
        if (local_status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, local_status)
        if (local_status /= 0) return
        if (.not. valid_vector_samples( &
            source_values, size(weights), mesh%n_triangles)) return
        call status_set(status, 0, "")
    end subroutine initialize_sampled_full_vector_load

    subroutine evaluate_full_vector_reference( &
            normal_family, bdm_basis, nedelec_basis, xi, eta, values, &
            derivatives, status)
        logical, intent(in) :: normal_family
        type(triangle_bdm_basis_t), intent(in) :: bdm_basis
        type(triangle_nedelec_second_kind_t), intent(in) :: nedelec_basis
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(:, :), derivatives(:)
        integer, intent(out) :: status

        if (normal_family) then
            call evaluate_triangle_bdm( &
                bdm_basis, xi, eta, values, derivatives, status)
        else
            call evaluate_triangle_nedelec_second_kind( &
                nedelec_basis, xi, eta, values, derivatives, status)
        end if
    end subroutine evaluate_full_vector_reference

    pure subroutine map_full_vector( &
            normal_family, jacobian, reference_values, reference_derivatives, &
            values, derivatives, status)
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: jacobian(2, 2), reference_values(:, :)
        real(dp), intent(in) :: reference_derivatives(:)
        real(dp), intent(out) :: values(:, :), derivatives(:)
        integer, intent(out) :: status

        if (normal_family) then
            call map_triangle_rt_contravariant( &
                jacobian, reference_values, reference_derivatives, values, &
                derivatives, status)
        else
            call map_triangle_nedelec_covariant( &
                jacobian, reference_values, reference_derivatives, values, &
                derivatives, status)
        end if
    end subroutine map_full_vector

    pure subroutine map_full_vector_jvp( &
            normal_family, jacobian, reference_values, reference_derivatives, &
            jacobian_dot, reference_values_dot, reference_derivatives_dot, &
            values_dot, derivatives_dot, status)
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: jacobian(2, 2), reference_values(:, :)
        real(dp), intent(in) :: reference_derivatives(:), jacobian_dot(2, 2)
        real(dp), intent(in) :: reference_values_dot(:, :)
        real(dp), intent(in) :: reference_derivatives_dot(:)
        real(dp), intent(out) :: values_dot(:, :), derivatives_dot(:)
        integer, intent(out) :: status

        if (normal_family) then
            call map_triangle_rt_contravariant_jvp( &
                jacobian, reference_values, reference_derivatives, &
                jacobian_dot, reference_values_dot, &
                reference_derivatives_dot, values_dot, derivatives_dot, status)
        else
            call map_triangle_nedelec_covariant_jvp( &
                jacobian, reference_values, reference_derivatives, &
                jacobian_dot, reference_values_dot, &
                reference_derivatives_dot, values_dot, derivatives_dot, status)
        end if
    end subroutine map_full_vector_jvp

    pure subroutine map_full_vector_vjp( &
            normal_family, jacobian, reference_values, reference_derivatives, &
            values_bar, derivatives_bar, jacobian_bar, reference_values_bar, &
            reference_derivatives_bar, status)
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: jacobian(2, 2), reference_values(:, :)
        real(dp), intent(in) :: reference_derivatives(:)
        real(dp), intent(in) :: values_bar(:, :), derivatives_bar(:)
        real(dp), intent(out) :: jacobian_bar(2, 2)
        real(dp), intent(out) :: reference_values_bar(:, :)
        real(dp), intent(out) :: reference_derivatives_bar(:)
        integer, intent(out) :: status

        if (normal_family) then
            call map_triangle_rt_contravariant_vjp( &
                jacobian, reference_values, reference_derivatives, values_bar, &
                derivatives_bar, jacobian_bar, reference_values_bar, &
                reference_derivatives_bar, status)
        else
            call map_triangle_nedelec_covariant_vjp( &
                jacobian, reference_values, reference_derivatives, values_bar, &
                derivatives_bar, jacobian_bar, reference_values_bar, &
                reference_derivatives_bar, status)
        end if
    end subroutine map_full_vector_vjp

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

end module fortfem_assembly_full_vector_arbitrary_order_2d
