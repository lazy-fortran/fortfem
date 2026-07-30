module fortfem_api_solvers_vector
    use fortfem_kinds, only: dp
    use fortfem_api_types, only: vector_bc_t, vector_function_t
    use fortfem_api_forms, only: compile_vector_form_csc, &
        compile_vector_form_rhs, form_equation_t
    use fortfem_advanced_solvers, only: solver_options_t, solver_stats_t, &
        solver_options, cg_solve, pcg_solve, bicgstab_solve, gmres_solve
    use fortfem_basis_edge_2d, only: edge_basis_2d_t
    use fortfem_sparse_direct, only: sparse_direct_solve_csc
    use fortsparse, only: csc_t, fortsparse_status_t
    implicit none

    private

    public :: solve_vector
    public :: solve_curl_curl_problem
    public :: solve_generic_vector_problem

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface

contains

    subroutine solve_vector(equation, Eh, bc, solver_type, options, stats)
        type(form_equation_t), intent(in) :: equation
        type(vector_function_t), intent(inout) :: Eh
        type(vector_bc_t), intent(in) :: bc
        character(len=*), intent(in), optional :: solver_type
        type(solver_options_t), intent(in), optional :: options
        type(solver_stats_t), intent(out), optional :: stats

        character(len=32) :: solver
        type(solver_options_t) :: local_opts
        type(solver_stats_t) :: local_stats

        solver = "gmres"
        if (present(solver_type)) solver = solver_type

        local_opts = solver_options(method="gmres", tolerance=1.0e-6_dp, &
            max_iterations=100, restart=20)
        if (present(options)) local_opts = options

        local_stats%converged = .false.
        local_stats%iterations = 0
        local_stats%final_residual = 0.0_dp
        local_stats%restarts = 0
        local_stats%method_used = ""

        write(*,*) "Solving vector problem: ", &
            trim(equation%lhs%description), " == ", &
            trim(equation%rhs%description)
        write(*,*) "Using solver: ", trim(solver)

        if ((trim(Eh%space%element_family) == "Nedelec" .or. &
            trim(Eh%space%element_family) == "Nedelec1" .or. &
            trim(Eh%space%element_family) == "Edge") .and. &
            Eh%space%degree >= 1) then
            call solve_compiled_vector_problem( &
                equation, Eh, bc, solver, local_opts, local_stats, &
                "Nedelec", Eh%space%degree, .false.)
        else if ((trim(Eh%space%element_family) == "RT" .or. &
                trim(Eh%space%element_family) == "Raviart-Thomas") .and. &
                Eh%space%degree >= 0) then
            call solve_compiled_vector_problem( &
                equation, Eh, bc, solver, local_opts, local_stats, &
                "RT", Eh%space%degree + 1, .true.)
        else if (trim(Eh%space%element_family) == "BDM" .and. &
                Eh%space%degree >= 1) then
            call solve_compiled_vector_problem( &
                equation, Eh, bc, solver, local_opts, local_stats, &
                "BDM", Eh%space%degree + 1, .true.)
        else if ((trim(Eh%space%element_family) == "Nedelec2" .or. &
                trim(Eh%space%element_family) == "Nedelec-second") .and. &
                Eh%space%degree >= 1) then
            call solve_compiled_vector_problem( &
                equation, Eh, bc, solver, local_opts, local_stats, &
                "Nedelec2", Eh%space%degree + 1, .false.)
        else if (index(equation%lhs%description, "curl") > 0) then
            call solve_curl_curl_problem(Eh, bc, solver, local_opts, &
                local_stats)
        else
            call solve_generic_vector_problem(Eh, bc)
        end if

        if (present(stats)) then
            stats = local_stats
        end if
    end subroutine solve_vector

    subroutine solve_compiled_vector_problem( &
            equation, field, boundary_condition, solver_type, options, stats, &
            family, edge_moment_count, normal_family)
        type(form_equation_t), intent(in) :: equation
        type(vector_function_t), intent(inout) :: field
        type(vector_bc_t), intent(in) :: boundary_condition
        character(len=*), intent(in) :: solver_type
        type(solver_options_t), intent(in) :: options
        type(solver_stats_t), intent(out) :: stats
        character(len=*), intent(in) :: family
        integer, intent(in) :: edge_moment_count
        logical, intent(in) :: normal_family

        type(csc_t) :: sparse_matrix
        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: dense_matrix(:, :), prescribed_values(:)
        real(dp), allocatable :: right_hand_side(:), solution(:)
        real(dp) :: edge_component, edge_vector(2), prescribed_value
        logical, allocatable :: prescribed_mask(:)
        integer :: boundary_index, degree_of_freedom, edge, dof_count, moment
        integer :: assembly_degree

        assembly_degree = field%space%degree
        dof_count = field%space%ndof
        allocate(right_hand_side(dof_count), prescribed_values(dof_count))
        allocate(prescribed_mask(dof_count))
        allocate(solution(dof_count))
        call compile_vector_form_csc( &
            equation%lhs, field%space%mesh%data, family, assembly_degree, &
            2 * edge_moment_count + 2, &
            sparse_matrix, sparse_status)
        if (sparse_status%code /= 0) then
            error stop "solve: unsupported Nedelec bilinear form"
        end if
        call compile_vector_form_rhs( &
            equation%rhs, field%space%mesh%data, family, assembly_degree, &
            2 * edge_moment_count + 2, &
            right_hand_side, sparse_status)
        if (sparse_status%code /= 0) then
            error stop "solve: unsupported Nedelec linear form"
        end if
        if (sparse_matrix%nrow /= dof_count .or. &
            sparse_matrix%ncol /= dof_count) then
            error stop "solve: Nedelec function-space dimension mismatch"
        end if
        prescribed_values = 0.0_dp
        prescribed_mask = .false.
        do boundary_index = 1, &
                size(field%space%mesh%data%boundary_edges)
            edge = field%space%mesh%data%boundary_edges(boundary_index)
            edge_vector = field%space%mesh%data%vertices(:, &
                field%space%mesh%data%edges(2, edge)) - &
                field%space%mesh%data%vertices(:, &
                field%space%mesh%data%edges(1, edge))
            if (normal_family) then
                edge_component = edge_vector(1)
                edge_vector(1) = edge_vector(2)
                edge_vector(2) = -edge_component
            end if
            do moment = 1, edge_moment_count
                degree_of_freedom = &
                    field%space%mesh%data%edge_to_dof(edge) * &
                    edge_moment_count + moment
                if (allocated(boundary_condition%edge_values)) then
                    prescribed_value = boundary_condition%edge_values( &
                        degree_of_freedom)
                else if (moment == 1) then
                    prescribed_value = &
                        dot_product(boundary_condition%values, edge_vector)
                else
                    prescribed_value = 0.0_dp
                end if
                prescribed_values(degree_of_freedom) = prescribed_value
                prescribed_mask(degree_of_freedom) = .true.
            end do
        end do

        solution = 0.0_dp
        select case (trim(solver_type))
        case ("direct")
            call solve_sparse_dirichlet( &
                sparse_matrix, right_hand_side, prescribed_values, &
                prescribed_mask, solution, stats)
        case default
            if (edge_moment_count /= 1 .or. normal_family) then
                error stop "solve: arbitrary-order vector form needs direct"
            end if
            allocate(dense_matrix(dof_count, dof_count))
            call copy_csc_to_dense(sparse_matrix, dense_matrix)
            call apply_dense_dirichlet( &
                dense_matrix, right_hand_side, prescribed_values, &
                field%space%mesh%data%n_interior_dofs)
            call gmres_solve( &
                dense_matrix, right_hand_side, solution, options, stats)
        end select

        if (.not. allocated(field%values)) then
            allocate(field%values(dof_count, 2))
        end if
        field%values(:, 1) = solution
        field%values(:, 2) = 0.0_dp
    end subroutine solve_compiled_vector_problem

    subroutine solve_sparse_dirichlet( &
            matrix, rhs, prescribed, prescribed_mask, solution, stats)
        type(csc_t), intent(in) :: matrix
        real(dp), intent(in) :: rhs(:), prescribed(:)
        logical, intent(in) :: prescribed_mask(:)
        real(dp), intent(out) :: solution(:)
        type(solver_stats_t), intent(out) :: stats

        integer, allocatable :: column_pointers(:), free_dofs(:)
        integer, allocatable :: free_index(:), row_indices(:)
        real(dp), allocatable :: interior_rhs(:), interior_solution(:)
        real(dp), allocatable :: residual(:), values(:)
        integer :: column, entry, free_count, interior_entry, row, status

        solution = prescribed
        call initialize_direct_stats(stats)
        free_count = count(.not. prescribed_mask)
        if (free_count == 0) then
            stats%converged = .true.
            return
        end if

        allocate(free_dofs(free_count), free_index(matrix%nrow))
        free_count = 0
        do column = 1, matrix%ncol
            if (prescribed_mask(column)) cycle
            free_count = free_count + 1
            free_dofs(free_count) = column
        end do
        free_index = 0
        do column = 1, free_count
            free_index(free_dofs(column)) = column
        end do
        allocate(column_pointers(free_count + 1))
        interior_entry = 0
        do column = 1, free_count
            column_pointers(column) = interior_entry + 1
            do entry = matrix%col_ptr(free_dofs(column)), &
                    matrix%col_ptr(free_dofs(column) + 1) - 1
                if (.not. prescribed_mask(matrix%row_idx(entry))) then
                    interior_entry = interior_entry + 1
                end if
            end do
        end do
        column_pointers(free_count + 1) = interior_entry + 1
        allocate(row_indices(interior_entry), values(interior_entry))

        interior_entry = 0
        do column = 1, free_count
            do entry = matrix%col_ptr(free_dofs(column)), &
                    matrix%col_ptr(free_dofs(column) + 1) - 1
                row = matrix%row_idx(entry)
                if (.not. prescribed_mask(row)) then
                    interior_entry = interior_entry + 1
                    row_indices(interior_entry) = free_index(row)
                    values(interior_entry) = matrix%val(entry)
                end if
            end do
        end do

        allocate(interior_rhs(free_count))
        allocate(interior_solution(free_count), residual(free_count))
        interior_rhs = rhs(free_dofs)
        do column = 1, matrix%ncol
            if (.not. prescribed_mask(column)) cycle
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                row = matrix%row_idx(entry)
                if (.not. prescribed_mask(row)) then
                    interior_rhs(free_index(row)) = &
                        interior_rhs(free_index(row)) - &
                        matrix%val(entry) * prescribed(column)
                end if
            end do
        end do

        call sparse_direct_solve_csc( &
            free_count, column_pointers, row_indices, values, &
            interior_rhs, interior_solution, status)
        if (status /= 0) error stop "solve: sparse Nedelec solve failed"
        solution(free_dofs) = interior_solution

        residual = -interior_rhs
        do column = 1, free_count
            do entry = column_pointers(column), &
                    column_pointers(column + 1) - 1
                residual(row_indices(entry)) = &
                    residual(row_indices(entry)) + &
                    values(entry) * interior_solution(column)
            end do
        end do
        stats%final_residual = sqrt(sum(residual**2))
        stats%converged = stats%final_residual < 1.0e-10_dp
    end subroutine solve_sparse_dirichlet

    subroutine initialize_direct_stats(stats)
        type(solver_stats_t), intent(out) :: stats

        stats%converged = .false.
        stats%iterations = 1
        stats%final_residual = 0.0_dp
        stats%solve_time = 0.0_dp
        stats%memory_usage = 0
        stats%method_used = "fortsparse_direct"
        stats%restarts = 0
        stats%parallel_efficiency = 0.0_dp
        stats%condition_estimate = 0.0_dp
    end subroutine initialize_direct_stats

    subroutine apply_dense_dirichlet( &
            matrix, rhs, prescribed, interior_count)
        real(dp), intent(inout) :: matrix(:, :), rhs(:)
        real(dp), intent(in) :: prescribed(:)
        integer, intent(in) :: interior_count

        integer :: degree_of_freedom

        do degree_of_freedom = interior_count + 1, size(rhs)
            rhs = rhs - matrix(:, degree_of_freedom) * &
                prescribed(degree_of_freedom)
            matrix(:, degree_of_freedom) = 0.0_dp
            matrix(degree_of_freedom, :) = 0.0_dp
            matrix(degree_of_freedom, degree_of_freedom) = 1.0_dp
            rhs(degree_of_freedom) = prescribed(degree_of_freedom)
        end do
    end subroutine apply_dense_dirichlet

    pure subroutine copy_csc_to_dense(sparse_matrix, dense_matrix)
        type(csc_t), intent(in) :: sparse_matrix
        real(dp), intent(out) :: dense_matrix(:, :)

        integer :: column, entry

        dense_matrix = 0.0_dp
        do column = 1, sparse_matrix%ncol
            do entry = sparse_matrix%col_ptr(column), &
                    sparse_matrix%col_ptr(column + 1) - 1
                dense_matrix(sparse_matrix%row_idx(entry), column) = &
                    dense_matrix(sparse_matrix%row_idx(entry), column) + &
                    sparse_matrix%val(entry)
            end do
        end do
    end subroutine copy_csc_to_dense

    subroutine solve_curl_curl_problem(Eh, bc, solver_type, options, stats)
        type(vector_function_t), intent(inout) :: Eh
        type(vector_bc_t), intent(in) :: bc
        character(len=*), intent(in) :: solver_type
        type(solver_options_t), intent(in) :: options
        type(solver_stats_t), intent(out) :: stats

        real(dp), allocatable :: A(:, :), b(:), x(:)
        integer :: ndof

        ndof = Eh%space%ndof
        allocate(A(ndof, ndof), b(ndof), x(ndof))
        x = 0.0_dp

        call assemble_curl_curl_system(Eh, A, b)

        select case (trim(solver_type))
        case ("gmres")
            call gmres_solve(A, b, x, options, stats)
        case ("direct")
            call solve_direct_vector(A, b, x)
            stats%converged = .true.
            stats%iterations = 1
            stats%final_residual = sqrt(sum((matmul(A, x) - b)**2))
            stats%solve_time = 0.0_dp
            stats%memory_usage = 0
            stats%method_used = "lapack_lu"
            stats%restarts = 0
            stats%parallel_efficiency = 0.0_dp
            stats%condition_estimate = 0.0_dp
        case default
            call gmres_solve(A, b, x, options, stats)
        end select

        if (allocated(Eh%values)) then
            do ndof = 1, size(x)
                Eh%values(ndof, 1) = x(ndof)
                Eh%values(ndof, 2) = 0.0_dp
            end do
        end if

        deallocate(A, b, x)
    end subroutine solve_curl_curl_problem

    subroutine assemble_curl_curl_system(Eh, A, b)
        type(vector_function_t), intent(in) :: Eh
        real(dp), intent(inout) :: A(:, :), b(:)

        integer :: ndof, e, i
        type(edge_basis_2d_t) :: edge_basis

        ndof = Eh%space%ndof

        A = 0.0_dp
        b = 0.0_dp

        call edge_basis%init(Eh%space%mesh%data)

        do e = 1, Eh%space%mesh%data%n_triangles
            call add_curl_curl_triangle(Eh, ndof, e, A, b)
        end do

        do i = 1, ndof
            if (Eh%space%mesh%data%is_boundary_edge(i)) then
                A(i, :) = 0.0_dp
                A(i, i) = 1.0_dp
                b(i) = 0.0_dp
            end if
        end do
    end subroutine assemble_curl_curl_system

    subroutine add_curl_curl_triangle(Eh, ndof, triangle_id, A, b)
        type(vector_function_t), intent(in) :: Eh
        integer, intent(in) :: ndof, triangle_id
        real(dp), intent(inout) :: A(:, :), b(:)

        integer :: v1, v2, v3
        integer :: edge1, edge2, edge3
        real(dp) :: x1, y1, x2, y2, x3, y3
        real(dp) :: area

        v1 = Eh%space%mesh%data%triangles(1, triangle_id)
        v2 = Eh%space%mesh%data%triangles(2, triangle_id)
        v3 = Eh%space%mesh%data%triangles(3, triangle_id)

        edge1 = 3*(triangle_id - 1) + 1
        edge2 = 3*(triangle_id - 1) + 2
        edge3 = 3*(triangle_id - 1) + 3

        x1 = Eh%space%mesh%data%vertices(1, v1)
        y1 = Eh%space%mesh%data%vertices(2, v1)
        x2 = Eh%space%mesh%data%vertices(1, v2)
        y2 = Eh%space%mesh%data%vertices(2, v2)
        x3 = Eh%space%mesh%data%vertices(1, v3)
        y3 = Eh%space%mesh%data%vertices(2, v3)

        area = 0.5_dp*abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))

        call accumulate_curl_curl_for_edges(edge1, edge2, edge3, ndof, area, &
            A, b)
    end subroutine add_curl_curl_triangle

    subroutine accumulate_curl_curl_for_edges(edge1, edge2, edge3, ndof, &
            area, A, b)
        integer, intent(in) :: edge1, edge2, edge3, ndof
        real(dp), intent(in) :: area
        real(dp), intent(inout) :: A(:, :), b(:)

        integer :: i, j
        real(dp) :: curl_basis_i, curl_basis_j

        do i = 1, 3
            do j = 1, 3
                curl_basis_i = 1.0_dp/area
                curl_basis_j = 1.0_dp/area

                if (i == 1 .and. edge1 > 0 .and. edge1 <= ndof) then
                    if (j == 1 .and. edge1 > 0 .and. edge1 <= ndof) then
                        A(edge1, edge1) = A(edge1, edge1) + area* &
                            curl_basis_i*curl_basis_j
                    end if
                    if (j == 2 .and. edge2 > 0 .and. edge2 <= ndof) then
                        A(edge1, edge2) = A(edge1, edge2) + area* &
                            curl_basis_i*curl_basis_j
                    end if
                    if (j == 3 .and. edge3 > 0 .and. edge3 <= ndof) then
                        A(edge1, edge3) = A(edge1, edge3) + area* &
                            curl_basis_i*curl_basis_j
                    end if
                end if
            end do

            if (i == 1 .and. edge1 > 0 .and. edge1 <= ndof) then
                A(edge1, edge1) = A(edge1, edge1) + area/3.0_dp
            end if
        end do

        if (edge1 > 0 .and. edge1 <= ndof) then
            b(edge1) = b(edge1) + area/3.0_dp
        end if
        if (edge2 > 0 .and. edge2 <= ndof) then
            b(edge2) = b(edge2) + area/3.0_dp
        end if
        if (edge3 > 0 .and. edge3 <= ndof) then
            b(edge3) = b(edge3) + area/3.0_dp
        end if
    end subroutine accumulate_curl_curl_for_edges

    subroutine solve_direct_vector(A, b, x)
        real(dp), intent(in) :: A(:, :), b(:)
        real(dp), intent(out) :: x(:)

        real(dp), allocatable :: A_work(:, :), b_work(:)
        integer :: n, info, ipiv(size(A, 1))

        n = size(A, 1)
        allocate(A_work(n, n), b_work(n))

        A_work = A
        b_work = b

        call dgesv(n, 1, A_work, n, ipiv, b_work, n, info)

        if (info == 0) then
            x = b_work
        else
            write(*,*) "Warning: Direct vector solver failed with info =", info
            x = 0.0_dp
        end if

        deallocate(A_work, b_work)
    end subroutine solve_direct_vector

    subroutine solve_generic_vector_problem(Eh, bc)
        type(vector_function_t), intent(inout) :: Eh
        type(vector_bc_t), intent(in) :: bc

        if (allocated(Eh%values)) then
            Eh%values(:, 1) = bc%values(1)
            Eh%values(:, 2) = bc%values(2)
        end if
    end subroutine solve_generic_vector_problem

end module fortfem_api_solvers_vector
