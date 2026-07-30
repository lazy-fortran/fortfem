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
        map_triangle_nedelec_covariant, map_triangle_rt_contravariant
    use fortsparse, only: &
        csc_from_triplet, csc_t, FORTSPARSE_INVALID_MATRIX, &
        fortsparse_status_t, status_set
    implicit none

    private

    public :: assemble_triangle_bdm_div_mass_csc
    public :: assemble_triangle_bdm_div_mass_element
    public :: assemble_triangle_bdm_cell_vector_load
    public :: assemble_triangle_nedelec_second_curl_mass_csc
    public :: assemble_triangle_nedelec_second_curl_mass_element

contains

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

end module fortfem_assembly_full_vector_arbitrary_order_2d
