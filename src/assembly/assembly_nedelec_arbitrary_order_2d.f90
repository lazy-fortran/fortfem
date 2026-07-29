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
    use fortfem_triangle_piola_maps, only: map_triangle_nedelec_covariant
    use fortsparse, only: &
        csc_from_triplet, csc_t, FORTSPARSE_INVALID_MATRIX, &
        fortsparse_status_t, status_set
    implicit none

    private

    public :: assemble_triangle_nedelec_curl_mass_element
    public :: assemble_triangle_nedelec_curl_mass_csc
    public :: assemble_triangle_nedelec_curl_csc

contains

    subroutine assemble_triangle_nedelec_curl_mass_element( &
            vertices, order, quadrature_degree, matrix, status)
        real(dp), intent(in) :: vertices(2, 3)
        integer, intent(in) :: order, quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        type(triangle_nedelec_first_kind_t) :: basis
        real(dp), allocatable :: eta(:), physical_curls(:)
        real(dp), allocatable :: physical_values(:, :), reference_curls(:)
        real(dp), allocatable :: reference_values(:, :), weights(:), xi(:)
        real(dp) :: determinant, jacobian(2, 2), physical_weight
        integer :: basis_dof_count, column, point, row

        status = 1
        if (order < 1 .or. quadrature_degree < 0) return

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
                        physical_curls(row) * physical_curls(column) + &
                        dot_product( &
                        physical_values(:, row), physical_values(:, column)))
                end do
            end do
        end do
        status = 0
    end subroutine assemble_triangle_nedelec_curl_mass_element

    subroutine assemble_triangle_nedelec_curl_mass_csc( &
            mesh, order, quadrature_degree, matrix, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: order, quadrature_degree
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), global_dofs(:, :), rows(:)
        integer, allocatable :: transforms(:, :)
        real(dp), allocatable :: element_matrix(:, :), values(:)
        real(dp) :: vertices(2, 3)
        integer :: column, entry, global_dof_count, local_dof_count
        integer :: local_status, row, triangle

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
                local_status)
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

end module fortfem_assembly_nedelec_arbitrary_order_2d
