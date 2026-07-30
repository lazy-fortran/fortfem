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
    use fortnum_linalg, only: det3, inv3
    use fortsparse, only: csc_t, FORTSPARSE_INTERNAL_ERROR, &
        FORTSPARSE_INVALID_MATRIX, fortsparse_status_t, status_set
    implicit none
    private

    public :: evaluate_tetra_lagrange_solution
    public :: solve_tetra_lagrange_diffusion_reaction
    public :: solve_tetra_lagrange_poisson

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

        type(tetra_lagrange_t) :: basis
        integer, allocatable :: global_dofs(:, :)
        real(dp), allocatable :: basis_gradients(:, :), basis_values(:)
        real(dp) :: determinant, inverse_jacobian(3, 3), jacobian(3, 3)
        real(dp) :: reference_gradient(3)
        integer :: dof, global_count, inverse_status, local_status

        value = 0.0_dp
        gradient = 0.0_dp
        status = 1
        if (size(vertices, 1) /= 3) return
        if (size(tetrahedra, 1) /= 4) return
        if (tetrahedron < 1 .or. tetrahedron > size(tetrahedra, 2)) return
        call initialize_tetra_lagrange(degree, basis, local_status)
        if (local_status /= 0) return
        call build_tetra_lagrange_dof_map( &
            degree, tetrahedra, global_dofs, global_count, local_status)
        if (local_status /= 0) return
        if (size(solution) /= global_count) return
        allocate ( &
            basis_values(tetra_lagrange_dof_count(basis)), &
            basis_gradients(3, tetra_lagrange_dof_count(basis)))
        call evaluate_tetra_lagrange( &
            basis, reference_point, basis_values, basis_gradients, local_status)
        if (local_status /= 0) return
        jacobian(:, 1) = &
            vertices(:, tetrahedra(2, tetrahedron)) - &
            vertices(:, tetrahedra(1, tetrahedron))
        jacobian(:, 2) = &
            vertices(:, tetrahedra(3, tetrahedron)) - &
            vertices(:, tetrahedra(1, tetrahedron))
        jacobian(:, 3) = &
            vertices(:, tetrahedra(4, tetrahedron)) - &
            vertices(:, tetrahedra(1, tetrahedron))
        determinant = det3(jacobian)
        if (determinant <= 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**3)) return
        call inv3(jacobian, inverse_jacobian, inverse_status)
        if (inverse_status /= 0) return

        reference_gradient = 0.0_dp
        do dof = 1, size(basis_values)
            value = value + &
                solution(global_dofs(dof, tetrahedron))*basis_values(dof)
            reference_gradient = reference_gradient + &
                solution(global_dofs(dof, tetrahedron))* &
                basis_gradients(:, dof)
        end do
        gradient = matmul(transpose(inverse_jacobian), reference_gradient)
        status = 0
    end subroutine evaluate_tetra_lagrange_solution

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
