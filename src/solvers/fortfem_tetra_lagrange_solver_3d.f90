module fortfem_tetra_lagrange_solver_3d
    use fortfem_assembly_tetra_lagrange_arbitrary_order_3d, only: &
        assemble_tetra_lagrange_scalar_load, &
        assemble_tetra_lagrange_stiffness_csc, scalar_source_3d
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_solve_constrained, &
        sparse_direct_solve_csc
    use fortfem_tetra_lagrange_arbitrary_order, only: &
        evaluate_tetra_lagrange, initialize_tetra_lagrange, &
        tetra_lagrange_barycentric_indices, tetra_lagrange_dof_count, &
        tetra_lagrange_nodes, tetra_lagrange_t
    use fortfem_tetra_lagrange_global_dof_map, only: &
        build_tetra_lagrange_dof_map
    use fortnum_linalg, only: det3, inv3, inv3_jvp, inv3_vjp
    use fortsparse, only: csc_t, FORTSPARSE_INTERNAL_ERROR, &
        FORTSPARSE_INVALID_MATRIX, fortsparse_status_t, status_set
    implicit none
    private

    type :: tetra_lagrange_solution_evaluator_t
        type(tetra_lagrange_t) :: basis
        integer, allocatable :: global_dofs(:, :)
        real(dp), allocatable :: basis_gradients(:, :), basis_values(:)
        real(dp), allocatable :: inverse_jacobians(:, :, :)
        integer :: global_dof_count = 0
    end type tetra_lagrange_solution_evaluator_t

    interface assignment(=)
        module procedure assign_tetra_lagrange_solution_evaluator
    end interface

    public :: assignment(=)
    public :: evaluate_tetra_lagrange_solution
    public :: evaluate_tetra_lagrange_solution_jvp
    public :: evaluate_tetra_lagrange_solution_vjp
    public :: evaluate_tetra_lagrange_solution_prepared
    public :: initialize_tetra_lagrange_solution_evaluator
    public :: solve_tetra_lagrange_diffusion_reaction
    public :: solve_tetra_lagrange_poisson
    public :: tetra_lagrange_solution_evaluator_t

contains

    subroutine solve_tetra_lagrange_poisson( &
            vertices, tetrahedra, degree, source, boundary_value, solution, &
            status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree
        procedure(scalar_source_3d) :: source, boundary_value
        real(dp), allocatable, intent(out) :: solution(:)
        type(fortsparse_status_t), intent(out) :: status

        call solve_tetra_lagrange_diffusion_reaction( &
            vertices, tetrahedra, degree, source, 1.0_dp, 0.0_dp, solution, &
            status, boundary_value)
    end subroutine solve_tetra_lagrange_poisson

    subroutine solve_tetra_lagrange_diffusion_reaction( &
            vertices, tetrahedra, degree, source, stiffness_coefficient, &
            mass_coefficient, solution, status, boundary_value)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree
        procedure(scalar_source_3d) :: source
        real(dp), intent(in) :: stiffness_coefficient, mass_coefficient
        real(dp), allocatable, intent(out) :: solution(:)
        type(fortsparse_status_t), intent(out) :: status
        procedure(scalar_source_3d), optional :: boundary_value

        type(csc_t) :: matrix
        logical, allocatable :: constrained(:)
        real(dp), allocatable :: constrained_values(:), right_hand_side(:)
        integer :: solve_status

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Tetrahedral H1 diffusion-reaction solve failed")
        if (degree < 1 .or. degree > 4) return
        if (stiffness_coefficient < 0.0_dp) return
        if (mass_coefficient < 0.0_dp) return
        if (stiffness_coefficient == 0.0_dp) then
            if (mass_coefficient == 0.0_dp) return
        end if
        if (.not. present(boundary_value)) then
            if (mass_coefficient == 0.0_dp) return
        end if

        call assemble_tetra_lagrange_stiffness_csc( &
            vertices, tetrahedra, degree, 2*degree + 4, matrix, status, &
            stiffness_coefficient, mass_coefficient)
        if (status%code /= 0) return
        call assemble_tetra_lagrange_scalar_load( &
            vertices, tetrahedra, degree, 2*degree + 4, source, &
            right_hand_side, status)
        if (status%code /= 0) return
        if (size(right_hand_side) /= matrix%nrow) return
        allocate (solution(matrix%nrow))

        if (present(boundary_value)) then
            call build_dirichlet_data( &
                vertices, tetrahedra, degree, boundary_value, constrained, &
                constrained_values, solve_status)
            if (solve_status /= 0) then
                call status_set( &
                    status, FORTSPARSE_INTERNAL_ERROR, &
                    "Tetrahedral H1 boundary construction failed")
                return
            end if
            call sparse_direct_solve_constrained( &
                matrix, right_hand_side, constrained, constrained_values, &
                solution, solve_status)
        else
            call sparse_direct_solve_csc( &
                matrix%nrow, matrix%col_ptr, matrix%row_idx, matrix%val, &
                right_hand_side, solution, solve_status)
        end if
        if (solve_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INTERNAL_ERROR, &
                "Tetrahedral H1 FortSparse solve failed")
            return
        end if
        call status_set(status, 0, "")
    end subroutine solve_tetra_lagrange_diffusion_reaction

    subroutine evaluate_tetra_lagrange_solution( &
            vertices, tetrahedra, degree, solution, tetrahedron, &
            reference_point, value, gradient, status)
        real(dp), intent(in) :: vertices(:, :), solution(:)
        integer, intent(in) :: tetrahedra(:, :), degree, tetrahedron
        real(dp), intent(in) :: reference_point(3)
        real(dp), intent(out) :: value, gradient(3)
        integer, intent(out) :: status

        type(tetra_lagrange_solution_evaluator_t) :: evaluator

        value = 0.0_dp
        gradient = 0.0_dp
        call initialize_tetra_lagrange_solution_evaluator( &
            vertices, tetrahedra, degree, evaluator, status)
        if (status /= 0) return
        call evaluate_tetra_lagrange_solution_prepared( &
            evaluator, solution, tetrahedron, reference_point, value, &
            gradient, status)
    end subroutine evaluate_tetra_lagrange_solution

    subroutine evaluate_tetra_lagrange_solution_jvp( &
            vertices, tetrahedra, degree, solution, tetrahedron, &
            reference_point, vertices_dot, solution_dot, value_dot, &
            gradient_dot, status)
        real(dp), intent(in) :: vertices(:, :), solution(:)
        integer, intent(in) :: tetrahedra(:, :), degree, tetrahedron
        real(dp), intent(in) :: reference_point(3)
        real(dp), intent(in) :: vertices_dot(:, :), solution_dot(:)
        real(dp), intent(out) :: value_dot, gradient_dot(3)
        integer, intent(out) :: status

        type(tetra_lagrange_t) :: basis
        integer, allocatable :: global_dofs(:, :)
        real(dp), allocatable :: basis_gradients(:, :), basis_values(:)
        real(dp) :: inverse(3, 3), inverse_dot(3, 3)
        real(dp) :: jacobian(3, 3), jacobian_dot(3, 3)
        real(dp) :: reference_gradient(3), reference_gradient_dot(3)
        integer :: dof, global_dof_count

        value_dot = 0.0_dp
        gradient_dot = 0.0_dp
        status = 1
        if (.not. valid_observation_inputs( &
            vertices, tetrahedra, solution, tetrahedron)) return
        if (any(shape(vertices_dot) /= shape(vertices))) return
        if (size(solution_dot) /= size(solution)) return
        call initialize_tetra_lagrange(degree, basis, status)
        if (status /= 0) return
        call build_tetra_lagrange_dof_map( &
            degree, tetrahedra, global_dofs, global_dof_count, status)
        if (status /= 0) return
        if (size(solution) /= global_dof_count) return
        allocate(basis_values(tetra_lagrange_dof_count(basis)))
        allocate(basis_gradients(3, size(basis_values)))
        call evaluate_tetra_lagrange( &
            basis, reference_point, basis_values, basis_gradients, status)
        if (status /= 0) return
        call tetrahedron_jacobian( &
            vertices, tetrahedra(:, tetrahedron), jacobian)
        call tetrahedron_jacobian( &
            vertices_dot, tetrahedra(:, tetrahedron), jacobian_dot)
        call inv3_jvp( &
            jacobian, jacobian_dot, inverse, inverse_dot, status)
        if (status /= 0) return
        reference_gradient = 0.0_dp
        reference_gradient_dot = 0.0_dp
        do dof = 1, size(basis_values)
            value_dot = value_dot + &
                solution_dot(global_dofs(dof, tetrahedron))*basis_values(dof)
            reference_gradient = reference_gradient + &
                solution(global_dofs(dof, tetrahedron))* &
                basis_gradients(:, dof)
            reference_gradient_dot = reference_gradient_dot + &
                solution_dot(global_dofs(dof, tetrahedron))* &
                basis_gradients(:, dof)
        end do
        gradient_dot = matmul(transpose(inverse_dot), reference_gradient) + &
            matmul(transpose(inverse), reference_gradient_dot)
    end subroutine evaluate_tetra_lagrange_solution_jvp

    subroutine evaluate_tetra_lagrange_solution_vjp( &
            vertices, tetrahedra, degree, solution, tetrahedron, &
            reference_point, value_bar, gradient_bar, vertices_bar, &
            solution_bar, status)
        real(dp), intent(in) :: vertices(:, :), solution(:)
        integer, intent(in) :: tetrahedra(:, :), degree, tetrahedron
        real(dp), intent(in) :: reference_point(3), value_bar, gradient_bar(3)
        real(dp), intent(out) :: vertices_bar(:, :), solution_bar(:)
        integer, intent(out) :: status

        type(tetra_lagrange_t) :: basis
        integer, allocatable :: global_dofs(:, :)
        real(dp), allocatable :: basis_gradients(:, :), basis_values(:)
        real(dp) :: inverse(3, 3), inverse_bar(3, 3)
        real(dp) :: jacobian(3, 3), jacobian_bar(3, 3)
        real(dp) :: reference_gradient(3), reference_gradient_bar(3)
        integer :: dof, global_dof_count, vertex_ids(4)

        vertices_bar = 0.0_dp
        solution_bar = 0.0_dp
        status = 1
        if (.not. valid_observation_inputs( &
            vertices, tetrahedra, solution, tetrahedron)) return
        if (any(shape(vertices_bar) /= shape(vertices))) return
        if (size(solution_bar) /= size(solution)) return
        call initialize_tetra_lagrange(degree, basis, status)
        if (status /= 0) return
        call build_tetra_lagrange_dof_map( &
            degree, tetrahedra, global_dofs, global_dof_count, status)
        if (status /= 0) return
        if (size(solution) /= global_dof_count) return
        allocate(basis_values(tetra_lagrange_dof_count(basis)))
        allocate(basis_gradients(3, size(basis_values)))
        call evaluate_tetra_lagrange( &
            basis, reference_point, basis_values, basis_gradients, status)
        if (status /= 0) return
        vertex_ids = tetrahedra(:, tetrahedron)
        call tetrahedron_jacobian(vertices, vertex_ids, jacobian)
        reference_gradient = 0.0_dp
        do dof = 1, size(basis_values)
            reference_gradient = reference_gradient + &
                solution(global_dofs(dof, tetrahedron))* &
                basis_gradients(:, dof)
        end do
        call inv3(jacobian, inverse, status)
        if (status /= 0) return
        inverse_bar = spread(reference_gradient, 2, 3)* &
            spread(gradient_bar, 1, 3)
        call inv3_vjp( &
            jacobian, inverse_bar, inverse, jacobian_bar, status)
        if (status /= 0) return
        reference_gradient_bar = matmul(inverse, gradient_bar)
        do dof = 1, size(basis_values)
            solution_bar(global_dofs(dof, tetrahedron)) = &
                solution_bar(global_dofs(dof, tetrahedron)) + &
                value_bar*basis_values(dof) + &
                dot_product(reference_gradient_bar, basis_gradients(:, dof))
        end do
        vertices_bar(:, vertex_ids(1)) = &
            -jacobian_bar(:, 1) - jacobian_bar(:, 2) - jacobian_bar(:, 3)
        vertices_bar(:, vertex_ids(2)) = jacobian_bar(:, 1)
        vertices_bar(:, vertex_ids(3)) = jacobian_bar(:, 2)
        vertices_bar(:, vertex_ids(4)) = jacobian_bar(:, 3)
    end subroutine evaluate_tetra_lagrange_solution_vjp

    pure logical function valid_observation_inputs( &
            vertices, tetrahedra, solution, tetrahedron) result(valid)
        real(dp), intent(in) :: vertices(:, :), solution(:)
        integer, intent(in) :: tetrahedra(:, :), tetrahedron

        valid = size(vertices, 1) == 3
        if (.not. valid) return
        valid = size(tetrahedra, 1) == 4
        if (.not. valid) return
        valid = tetrahedron >= 1
        if (.not. valid) return
        valid = tetrahedron <= size(tetrahedra, 2)
        if (.not. valid) return
        valid = size(solution) >= 1
        if (.not. valid) return
        valid = all(tetrahedra >= 1)
        if (.not. valid) return
        valid = all(tetrahedra <= size(vertices, 2))
    end function valid_observation_inputs

    pure subroutine tetrahedron_jacobian(coordinates, vertex_ids, jacobian)
        real(dp), intent(in) :: coordinates(:, :)
        integer, intent(in) :: vertex_ids(4)
        real(dp), intent(out) :: jacobian(3, 3)

        jacobian(:, 1) = &
            coordinates(:, vertex_ids(2)) - coordinates(:, vertex_ids(1))
        jacobian(:, 2) = &
            coordinates(:, vertex_ids(3)) - coordinates(:, vertex_ids(1))
        jacobian(:, 3) = &
            coordinates(:, vertex_ids(4)) - coordinates(:, vertex_ids(1))
    end subroutine tetrahedron_jacobian

    subroutine initialize_tetra_lagrange_solution_evaluator( &
            vertices, tetrahedra, degree, evaluator, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree
        type(tetra_lagrange_solution_evaluator_t), intent(out) :: evaluator
        integer, intent(out) :: status

        real(dp) :: determinant, jacobian(3, 3)
        integer :: cell, inverse_status, local_dof_count

        evaluator%global_dof_count = 0
        status = 1
        if (size(vertices, 1) /= 3) return
        if (size(tetrahedra, 1) /= 4) return
        if (any(tetrahedra < 1)) return
        if (any(tetrahedra > size(vertices, 2))) return
        call initialize_tetra_lagrange(degree, evaluator%basis, status)
        if (status /= 0) return
        call build_tetra_lagrange_dof_map( &
            degree, tetrahedra, evaluator%global_dofs, &
            evaluator%global_dof_count, status)
        if (status /= 0) return
        local_dof_count = tetra_lagrange_dof_count(evaluator%basis)
        allocate ( &
            evaluator%basis_values(local_dof_count), &
            evaluator%basis_gradients(3, local_dof_count), &
            evaluator%inverse_jacobians(3, 3, size(tetrahedra, 2)))
        do cell = 1, size(tetrahedra, 2)
            jacobian(:, 1) = &
                vertices(:, tetrahedra(2, cell)) - &
                vertices(:, tetrahedra(1, cell))
            jacobian(:, 2) = &
                vertices(:, tetrahedra(3, cell)) - &
                vertices(:, tetrahedra(1, cell))
            jacobian(:, 3) = &
                vertices(:, tetrahedra(4, cell)) - &
                vertices(:, tetrahedra(1, cell))
            determinant = det3(jacobian)
            if (determinant <= 64.0_dp*epsilon(1.0_dp)* &
                max(1.0_dp, maxval(abs(jacobian))**3)) then
                status = 1
                return
            end if
            call inv3( &
                jacobian, evaluator%inverse_jacobians(:, :, cell), &
                inverse_status)
            if (inverse_status /= 0) then
                status = 1
                return
            end if
        end do
        status = 0
    end subroutine initialize_tetra_lagrange_solution_evaluator

    subroutine evaluate_tetra_lagrange_solution_prepared( &
            evaluator, solution, tetrahedron, reference_point, value, &
            gradient, status)
        type(tetra_lagrange_solution_evaluator_t), intent(inout) :: evaluator
        real(dp), intent(in) :: solution(:), reference_point(3)
        integer, intent(in) :: tetrahedron
        real(dp), intent(out) :: value, gradient(3)
        integer, intent(out) :: status

        real(dp) :: reference_gradient(3)
        integer :: dof, local_status

        value = 0.0_dp
        gradient = 0.0_dp
        status = 1
        if (tetrahedron < 1) return
        if (.not. allocated(evaluator%inverse_jacobians)) return
        if (tetrahedron > size(evaluator%inverse_jacobians, 3)) return
        if (size(solution) /= evaluator%global_dof_count) return
        call evaluate_tetra_lagrange( &
            evaluator%basis, reference_point, evaluator%basis_values, &
            evaluator%basis_gradients, local_status)
        if (local_status /= 0) return

        reference_gradient = 0.0_dp
        do dof = 1, size(evaluator%basis_values)
            value = value + &
                solution(evaluator%global_dofs(dof, tetrahedron))* &
                evaluator%basis_values(dof)
            reference_gradient = reference_gradient + &
                solution(evaluator%global_dofs(dof, tetrahedron))* &
                evaluator%basis_gradients(:, dof)
        end do
        gradient = matmul( &
            transpose(evaluator%inverse_jacobians(:, :, tetrahedron)), &
            reference_gradient)
        status = 0
    end subroutine evaluate_tetra_lagrange_solution_prepared

    subroutine assign_tetra_lagrange_solution_evaluator(left, right)
        type(tetra_lagrange_solution_evaluator_t), intent(out) :: left
        type(tetra_lagrange_solution_evaluator_t), intent(in) :: right

        left%basis = right%basis
        left%global_dof_count = right%global_dof_count
        if (allocated(right%global_dofs)) then
            allocate(left%global_dofs, source=right%global_dofs)
        end if
        if (allocated(right%basis_values)) then
            allocate(left%basis_values, source=right%basis_values)
        end if
        if (allocated(right%basis_gradients)) then
            allocate(left%basis_gradients, source=right%basis_gradients)
        end if
        if (allocated(right%inverse_jacobians)) then
            allocate( &
                left%inverse_jacobians, source=right%inverse_jacobians)
        end if
    end subroutine assign_tetra_lagrange_solution_evaluator

    subroutine build_dirichlet_data( &
            vertices, tetrahedra, degree, boundary_value, constrained, &
            constrained_values, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree
        procedure(scalar_source_3d) :: boundary_value
        logical, allocatable, intent(out) :: constrained(:)
        real(dp), allocatable, intent(out) :: constrained_values(:)
        integer, intent(out) :: status

        type(tetra_lagrange_t) :: basis
        integer, allocatable :: global_dofs(:, :), indices(:, :)
        real(dp), allocatable :: nodes(:, :)
        real(dp) :: jacobian(3, 3), physical_point(3)
        integer :: cell, dof, face(3), global_count, incidence
        integer :: local_status, node, opposite, other_cell

        status = 1
        call initialize_tetra_lagrange(degree, basis, local_status)
        if (local_status /= 0) return
        call tetra_lagrange_barycentric_indices(basis, indices)
        call tetra_lagrange_nodes(basis, nodes)
        call build_tetra_lagrange_dof_map( &
            degree, tetrahedra, global_dofs, global_count, local_status)
        if (local_status /= 0) return
        allocate (constrained(global_count), constrained_values(global_count))
        constrained = .false.
        constrained_values = 0.0_dp

        do cell = 1, size(tetrahedra, 2)
            jacobian(:, 1) = &
                vertices(:, tetrahedra(2, cell)) - vertices(:, tetrahedra(1, cell))
            jacobian(:, 2) = &
                vertices(:, tetrahedra(3, cell)) - vertices(:, tetrahedra(1, cell))
            jacobian(:, 3) = &
                vertices(:, tetrahedra(4, cell)) - vertices(:, tetrahedra(1, cell))
            do opposite = 1, 4
                face = pack( &
                    tetrahedra(:, cell), &
                    [(node /= opposite, node=1, 4)])
                incidence = 0
                do other_cell = 1, size(tetrahedra, 2)
                    if (all_vertices_in_cell(face, tetrahedra(:, other_cell))) then
                        incidence = incidence + 1
                    end if
                end do
                if (incidence /= 1) cycle
                do dof = 1, size(indices, 2)
                    if (indices(opposite, dof) /= 0) cycle
                    constrained(global_dofs(dof, cell)) = .true.
                    physical_point = vertices(:, tetrahedra(1, cell)) + &
                        matmul(jacobian, nodes(:, dof))
                    call boundary_value( &
                        physical_point(1), physical_point(2), physical_point(3), &
                        constrained_values(global_dofs(dof, cell)))
                end do
            end do
        end do
        status = 0
    end subroutine build_dirichlet_data

    pure logical function all_vertices_in_cell(face, cell) result(found)
        integer, intent(in) :: face(3), cell(4)

        integer :: vertex

        found = .true.
        do vertex = 1, 3
            found = found .and. any(cell == face(vertex))
        end do
    end function all_vertices_in_cell

end module fortfem_tetra_lagrange_solver_3d
