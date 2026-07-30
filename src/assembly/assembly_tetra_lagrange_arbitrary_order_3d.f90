module fortfem_assembly_tetra_lagrange_arbitrary_order_3d
    use fortfem_kinds, only: dp
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_tetra_lagrange_arbitrary_order, only: &
        evaluate_tetra_lagrange, initialize_tetra_lagrange, &
        tetra_lagrange_dof_count, tetra_lagrange_t
    use fortfem_tetra_lagrange_global_dof_map, only: &
        build_tetra_lagrange_dof_map
    use fortnum_linalg, only: det3, inv3
    use fortsparse, only: &
        csc_from_triplet, csc_t, FORTSPARSE_INVALID_MATRIX, &
        fortsparse_status_t, status_set
    implicit none

    private

    public :: assemble_tetra_lagrange_stiffness_csc
    public :: assemble_tetra_lagrange_stiffness_element
    public :: assemble_tetra_lagrange_scalar_load
    public :: scalar_source_3d

    abstract interface
        pure subroutine scalar_source_3d(x, y, z, value)
            import :: dp
            real(dp), intent(in) :: x, y, z
            real(dp), intent(out) :: value
        end subroutine scalar_source_3d
    end interface

contains

    subroutine assemble_tetra_lagrange_stiffness_element( &
            vertices, degree, quadrature_degree, matrix, status, &
            stiffness_coefficient, mass_coefficient)
        real(dp), intent(in) :: vertices(3, 4)
        integer, intent(in) :: degree, quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: stiffness_coefficient
        real(dp), intent(in), optional :: mass_coefficient

        type(tetra_lagrange_t) :: basis
        real(dp), allocatable :: gradients(:, :), values(:)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: determinant, inverse_jacobian(3, 3), jacobian(3, 3)
        real(dp) :: mass_weight, physical_gradient(3), physical_weight
        real(dp) :: stiffness_weight
        integer :: column, dof_count, inverse_status, point, row

        status = 1
        if (degree < 1 .or. degree > 4) return
        if (quadrature_degree < 0) return
        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
        determinant = det3(jacobian)
        if (determinant <= 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**3)) return
        call inv3(jacobian, inverse_jacobian, inverse_status)
        if (inverse_status /= 0) return
        call initialize_tetra_lagrange(degree, basis, status)
        if (status /= 0) return
        call tetra_duffy_quadrature( &
            quadrature_degree, x, y, z, weights, status)
        if (status /= 0) return

        stiffness_weight = 1.0_dp
        mass_weight = 0.0_dp
        if (present(stiffness_coefficient)) then
            stiffness_weight = stiffness_coefficient
        end if
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        dof_count = tetra_lagrange_dof_count(basis)
        allocate( &
            matrix(dof_count, dof_count), values(dof_count), &
            gradients(3, dof_count))
        matrix = 0.0_dp
        do point = 1, size(weights)
            call evaluate_tetra_lagrange( &
                basis, [x(point), y(point), z(point)], values, gradients, &
                status)
            if (status /= 0) return
            physical_weight = determinant*weights(point)
            do column = 1, dof_count
                do row = 1, dof_count
                    physical_gradient = matmul( &
                        transpose(inverse_jacobian), gradients(:, row))
                    matrix(row, column) = matrix(row, column) + &
                        physical_weight*(stiffness_weight*dot_product( &
                        physical_gradient, matmul( &
                        transpose(inverse_jacobian), gradients(:, column))) + &
                        mass_weight*values(row)*values(column))
                end do
            end do
        end do
        status = 0
    end subroutine assemble_tetra_lagrange_stiffness_element

    subroutine assemble_tetra_lagrange_stiffness_csc( &
            mesh_vertices, tetrahedra, degree, quadrature_degree, matrix, &
            status, stiffness_coefficient, mass_coefficient)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :)
        integer, intent(in) :: degree, quadrature_degree
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: stiffness_coefficient
        real(dp), intent(in), optional :: mass_coefficient

        integer, allocatable :: columns(:), global_dofs(:, :), rows(:)
        real(dp), allocatable :: element_matrix(:, :), triplet_values(:)
        real(dp) :: mass_weight, stiffness_weight, vertices(3, 4)
        integer :: column, dof_count, entry, global_count, local_status
        integer :: node, row, tetrahedron

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Tetrahedral H1 sparse assembly failed")
        if (size(mesh_vertices, 1) /= 3) return
        if (size(tetrahedra, 1) /= 4 .or. size(tetrahedra, 2) < 1) return
        if (any(tetrahedra < 1) .or. &
            any(tetrahedra > size(mesh_vertices, 2))) return
        stiffness_weight = 1.0_dp
        mass_weight = 0.0_dp
        if (present(stiffness_coefficient)) then
            stiffness_weight = stiffness_coefficient
        end if
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        call build_tetra_lagrange_dof_map( &
            degree, tetrahedra, global_dofs, global_count, local_status)
        if (local_status /= 0) return
        dof_count = size(global_dofs, 1)
        allocate(rows(dof_count**2*size(tetrahedra, 2)))
        allocate(columns(size(rows)), triplet_values(size(rows)))
        entry = 0
        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                vertices(:, node) = &
                    mesh_vertices(:, tetrahedra(node, tetrahedron))
            end do
            call assemble_tetra_lagrange_stiffness_element( &
                vertices, degree, quadrature_degree, element_matrix, &
                local_status, stiffness_weight, mass_weight)
            if (local_status /= 0) return
            do column = 1, dof_count
                do row = 1, dof_count
                    entry = entry + 1
                    rows(entry) = global_dofs(row, tetrahedron)
                    columns(entry) = global_dofs(column, tetrahedron)
                    triplet_values(entry) = element_matrix(row, column)
                end do
            end do
        end do
        call csc_from_triplet( &
            global_count, global_count, rows, columns, triplet_values, &
            matrix, status)
    end subroutine assemble_tetra_lagrange_stiffness_csc

    subroutine assemble_tetra_lagrange_scalar_load( &
            mesh_vertices, tetrahedra, degree, quadrature_degree, source, &
            right_hand_side, status)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :)
        integer, intent(in) :: degree, quadrature_degree
        procedure(scalar_source_3d) :: source
        real(dp), allocatable, intent(out) :: right_hand_side(:)
        type(fortsparse_status_t), intent(out) :: status

        type(tetra_lagrange_t) :: basis
        integer, allocatable :: global_dofs(:, :)
        real(dp), allocatable :: gradients(:, :), values(:)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: determinant, jacobian(3, 3), physical_point(3)
        real(dp) :: source_value, vertices(3, 4)
        integer :: dof, dof_count, global_count, local_status, node, point
        integer :: tetrahedron

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Tetrahedral H1 scalar load assembly failed")
        if (size(mesh_vertices, 1) /= 3) return
        if (size(tetrahedra, 1) /= 4 .or. size(tetrahedra, 2) < 1) return
        if (any(tetrahedra < 1)) return
        if (any(tetrahedra > size(mesh_vertices, 2))) return
        call initialize_tetra_lagrange(degree, basis, local_status)
        if (local_status /= 0) return
        call build_tetra_lagrange_dof_map( &
            degree, tetrahedra, global_dofs, global_count, local_status)
        if (local_status /= 0) return
        call tetra_duffy_quadrature( &
            quadrature_degree, x, y, z, weights, local_status)
        if (local_status /= 0) return
        dof_count = tetra_lagrange_dof_count(basis)
        allocate ( &
            right_hand_side(global_count), values(dof_count), &
            gradients(3, dof_count))
        right_hand_side = 0.0_dp

        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                vertices(:, node) = &
                    mesh_vertices(:, tetrahedra(node, tetrahedron))
            end do
            jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
            jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
            jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
            determinant = det3(jacobian)
            if (determinant <= 64.0_dp*epsilon(1.0_dp)* &
                max(1.0_dp, maxval(abs(jacobian))**3)) return
            do point = 1, size(weights)
                call evaluate_tetra_lagrange( &
                    basis, [x(point), y(point), z(point)], values, gradients, &
                    local_status)
                if (local_status /= 0) return
                physical_point = vertices(:, 1) + &
                    matmul(jacobian, [x(point), y(point), z(point)])
                call source( &
                    physical_point(1), physical_point(2), physical_point(3), &
                    source_value)
                do dof = 1, dof_count
                    right_hand_side(global_dofs(dof, tetrahedron)) = &
                        right_hand_side(global_dofs(dof, tetrahedron)) + &
                        determinant*weights(point)*source_value*values(dof)
                end do
            end do
        end do
        call status_set(status, 0, "")
    end subroutine assemble_tetra_lagrange_scalar_load

end module fortfem_assembly_tetra_lagrange_arbitrary_order_3d
