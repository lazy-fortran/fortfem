module fortfem_assembly_tetra_lagrange_arbitrary_order_3d
    use fortfem_kinds, only: dp
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_tetra_lagrange_arbitrary_order, only: &
        evaluate_tetra_lagrange, initialize_tetra_lagrange, &
        tetra_lagrange_dof_count, tetra_lagrange_t
    use fortfem_tetra_lagrange_global_dof_map, only: &
        build_tetra_lagrange_dof_map
    use fortnum_linalg, only: det3, det3_jvp, det3_vjp, inv3, inv3_jvp, &
        inv3_vjp
    use fortsparse, only: &
        csc_from_triplet, csc_t, FORTSPARSE_INVALID_MATRIX, &
        fortsparse_status_t, status_set
    implicit none

    private

    public :: assemble_tetra_lagrange_stiffness_csc
    public :: assemble_tetra_lagrange_stiffness_element
    public :: assemble_tetra_lagrange_stiffness_element_jvp
    public :: assemble_tetra_lagrange_stiffness_element_vjp
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

    subroutine assemble_tetra_lagrange_stiffness_element_jvp( &
            vertices, degree, quadrature_degree, stiffness_coefficient, &
            mass_coefficient, vertices_dot, stiffness_coefficient_dot, &
            mass_coefficient_dot, matrix_dot, status)
        real(dp), intent(in) :: vertices(3, 4), vertices_dot(3, 4)
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: stiffness_coefficient, mass_coefficient
        real(dp), intent(in) :: stiffness_coefficient_dot
        real(dp), intent(in) :: mass_coefficient_dot
        real(dp), allocatable, intent(out) :: matrix_dot(:, :)
        integer, intent(out) :: status

        type(tetra_lagrange_t) :: basis
        real(dp), allocatable :: gradients(:, :), physical_gradients(:, :)
        real(dp), allocatable :: physical_gradients_dot(:, :), values(:)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: determinant, determinant_dot, gradient_energy
        real(dp) :: inverse_jacobian(3, 3), inverse_jacobian_dot(3, 3)
        real(dp) :: jacobian(3, 3), jacobian_dot(3, 3), mass_energy
        integer :: column, dof_count, inverse_status, point, row

        status = 1
        if (allocated(matrix_dot)) deallocate(matrix_dot)
        if (degree < 1 .or. degree > 4 .or. quadrature_degree < 0) return
        call tetra_jacobian(vertices, jacobian)
        call tetra_jacobian(vertices_dot, jacobian_dot)
        determinant = det3(jacobian)
        if (.not. valid_jacobian(jacobian, determinant)) return
        call det3_jvp(jacobian, jacobian_dot, determinant_dot)
        call inv3_jvp( &
            jacobian, jacobian_dot, inverse_jacobian, inverse_jacobian_dot, &
            inverse_status)
        if (inverse_status /= 0) return
        call initialize_tetra_lagrange(degree, basis, status)
        if (status /= 0) return
        call tetra_duffy_quadrature( &
            quadrature_degree, x, y, z, weights, status)
        if (status /= 0) return
        dof_count = tetra_lagrange_dof_count(basis)
        allocate(matrix_dot(dof_count, dof_count), source=0.0_dp)
        allocate(values(dof_count), gradients(3, dof_count))
        allocate(physical_gradients(3, dof_count))
        allocate(physical_gradients_dot(3, dof_count))
        do point = 1, size(weights)
            call evaluate_tetra_lagrange( &
                basis, [x(point), y(point), z(point)], values, gradients, &
                status)
            if (status /= 0) return
            physical_gradients = &
                matmul(transpose(inverse_jacobian), gradients)
            physical_gradients_dot = &
                matmul(transpose(inverse_jacobian_dot), gradients)
            do column = 1, dof_count
                do row = 1, dof_count
                    gradient_energy = dot_product( &
                        physical_gradients(:, row), &
                        physical_gradients(:, column))
                    mass_energy = values(row)*values(column)
                    matrix_dot(row, column) = matrix_dot(row, column) + &
                        weights(point)*(determinant_dot*( &
                        stiffness_coefficient*gradient_energy + &
                        mass_coefficient*mass_energy) + determinant*( &
                        stiffness_coefficient_dot*gradient_energy + &
                        stiffness_coefficient*(dot_product( &
                        physical_gradients_dot(:, row), &
                        physical_gradients(:, column)) + dot_product( &
                        physical_gradients(:, row), &
                        physical_gradients_dot(:, column))) + &
                        mass_coefficient_dot*mass_energy))
                end do
            end do
        end do
        status = 0
    end subroutine assemble_tetra_lagrange_stiffness_element_jvp

    subroutine assemble_tetra_lagrange_stiffness_element_vjp( &
            vertices, degree, quadrature_degree, stiffness_coefficient, &
            mass_coefficient, matrix_bar, vertices_bar, &
            stiffness_coefficient_bar, mass_coefficient_bar, status)
        real(dp), intent(in) :: vertices(3, 4)
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: stiffness_coefficient, mass_coefficient
        real(dp), intent(in) :: matrix_bar(:, :)
        real(dp), intent(out) :: vertices_bar(3, 4)
        real(dp), intent(out) :: stiffness_coefficient_bar
        real(dp), intent(out) :: mass_coefficient_bar
        integer, intent(out) :: status

        type(tetra_lagrange_t) :: basis
        real(dp), allocatable :: gradients(:, :), physical_gradients(:, :)
        real(dp), allocatable :: physical_gradients_bar(:, :), values(:)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: determinant, determinant_bar
        real(dp) :: determinant_jacobian_bar(3, 3), gradient_energy
        real(dp) :: inverse_jacobian(3, 3), inverse_jacobian_bar(3, 3)
        real(dp) :: inverse_jacobian_jacobian_bar(3, 3)
        real(dp) :: jacobian(3, 3), mass_energy, seed
        integer :: column, dof_count, inverse_status, point, row

        vertices_bar = 0.0_dp
        stiffness_coefficient_bar = 0.0_dp
        mass_coefficient_bar = 0.0_dp
        status = 1
        if (degree < 1 .or. degree > 4 .or. quadrature_degree < 0) return
        call initialize_tetra_lagrange(degree, basis, status)
        if (status /= 0) return
        dof_count = tetra_lagrange_dof_count(basis)
        if (size(matrix_bar, 1) /= dof_count .or. &
            size(matrix_bar, 2) /= dof_count) return
        call tetra_jacobian(vertices, jacobian)
        determinant = det3(jacobian)
        if (.not. valid_jacobian(jacobian, determinant)) return
        call inv3(jacobian, inverse_jacobian, inverse_status)
        if (inverse_status /= 0) return
        call tetra_duffy_quadrature( &
            quadrature_degree, x, y, z, weights, status)
        if (status /= 0) return
        allocate(values(dof_count), gradients(3, dof_count))
        allocate(physical_gradients(3, dof_count))
        allocate(physical_gradients_bar(3, dof_count))
        determinant_bar = 0.0_dp
        inverse_jacobian_bar = 0.0_dp
        do point = 1, size(weights)
            call evaluate_tetra_lagrange( &
                basis, [x(point), y(point), z(point)], values, gradients, &
                status)
            if (status /= 0) return
            physical_gradients = &
                matmul(transpose(inverse_jacobian), gradients)
            physical_gradients_bar = 0.0_dp
            do column = 1, dof_count
                do row = 1, dof_count
                    seed = weights(point)*matrix_bar(row, column)
                    gradient_energy = dot_product( &
                        physical_gradients(:, row), &
                        physical_gradients(:, column))
                    mass_energy = values(row)*values(column)
                    determinant_bar = determinant_bar + seed*( &
                        stiffness_coefficient*gradient_energy + &
                        mass_coefficient*mass_energy)
                    stiffness_coefficient_bar = &
                        stiffness_coefficient_bar + &
                        seed*determinant*gradient_energy
                    mass_coefficient_bar = mass_coefficient_bar + &
                        seed*determinant*mass_energy
                    physical_gradients_bar(:, row) = &
                        physical_gradients_bar(:, row) + &
                        seed*determinant*stiffness_coefficient* &
                        physical_gradients(:, column)
                    physical_gradients_bar(:, column) = &
                        physical_gradients_bar(:, column) + &
                        seed*determinant*stiffness_coefficient* &
                        physical_gradients(:, row)
                end do
            end do
            inverse_jacobian_bar = inverse_jacobian_bar + &
                matmul(gradients, transpose(physical_gradients_bar))
        end do
        call inv3_vjp( &
            jacobian, inverse_jacobian_bar, inverse_jacobian, &
            inverse_jacobian_jacobian_bar, inverse_status)
        if (inverse_status /= 0) return
        call det3_vjp(jacobian, determinant_bar, determinant_jacobian_bar)
        call tetra_jacobian_vjp( &
            inverse_jacobian_jacobian_bar + determinant_jacobian_bar, &
            vertices_bar)
        status = 0
    end subroutine assemble_tetra_lagrange_stiffness_element_vjp

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

    pure logical function valid_jacobian( &
            jacobian, determinant) result(valid)
        real(dp), intent(in) :: jacobian(3, 3), determinant

        valid = determinant > 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**3)
    end function valid_jacobian

end module fortfem_assembly_tetra_lagrange_arbitrary_order_3d
