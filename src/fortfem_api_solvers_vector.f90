module fortfem_api_solvers_vector
    use fortfem_kinds, only: dp
    use fortfem_api_types, only: vector_bc_t, vector_function_t
    use fortfem_api_forms, only: compile_vector_form_csc, &
        compile_vector_form_rhs, form_equation_t
    use fortfem_advanced_solvers, only: solver_options_t, solver_stats_t, &
        solver_options, cg_solve, pcg_solve, bicgstab_solve, gmres_solve
    use fortfem_basis_edge_2d, only: edge_basis_2d_t
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

        if (trim(Eh%space%element_family) == "Nedelec" .and. &
            Eh%space%degree == 1) then
            call solve_compiled_nedelec_problem( &
                equation, Eh, bc, solver, local_opts, local_stats)
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

    subroutine solve_compiled_nedelec_problem( &
            equation, field, boundary_condition, solver_type, options, stats)
        type(form_equation_t), intent(in) :: equation
        type(vector_function_t), intent(inout) :: field
        type(vector_bc_t), intent(in) :: boundary_condition
        character(len=*), intent(in) :: solver_type
        type(solver_options_t), intent(in) :: options
        type(solver_stats_t), intent(out) :: stats

        type(csc_t) :: sparse_matrix
        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: dense_matrix(:, :), right_hand_side(:)
        real(dp), allocatable :: solution(:)
        real(dp) :: edge_vector(2), prescribed_value
        integer :: boundary_index, degree_of_freedom, edge, dof_count

        dof_count = field%space%ndof
        allocate(right_hand_side(dof_count), solution(dof_count))
        call compile_vector_form_csc( &
            equation%lhs, field%space%mesh%data, "Nedelec", 1, 4, &
            sparse_matrix, sparse_status)
        if (sparse_status%code /= 0) then
            error stop "solve: unsupported Nedelec bilinear form"
        end if
        call compile_vector_form_rhs( &
            equation%rhs, field%space%mesh%data, "Nedelec", 1, 4, &
            right_hand_side, sparse_status)
        if (sparse_status%code /= 0) then
            error stop "solve: unsupported Nedelec linear form"
        end if
        if (sparse_matrix%nrow /= dof_count .or. &
            sparse_matrix%ncol /= dof_count) then
            error stop "solve: Nedelec function-space dimension mismatch"
        end if
        allocate(dense_matrix(dof_count, dof_count))
        call copy_csc_to_dense(sparse_matrix, dense_matrix)

        do boundary_index = 1, &
                size(field%space%mesh%data%boundary_edges)
            edge = field%space%mesh%data%boundary_edges(boundary_index)
            degree_of_freedom = &
                field%space%mesh%data%edge_to_dof(edge) + 1
            edge_vector = field%space%mesh%data%vertices(:, &
                field%space%mesh%data%edges(2, edge)) - &
                field%space%mesh%data%vertices(:, &
                field%space%mesh%data%edges(1, edge))
            prescribed_value = &
                dot_product(boundary_condition%values, edge_vector)
            right_hand_side = right_hand_side - &
                dense_matrix(:, degree_of_freedom) * prescribed_value
            dense_matrix(:, degree_of_freedom) = 0.0_dp
            dense_matrix(degree_of_freedom, :) = 0.0_dp
            dense_matrix(degree_of_freedom, degree_of_freedom) = 1.0_dp
            right_hand_side(degree_of_freedom) = prescribed_value
        end do

        solution = 0.0_dp
        select case (trim(solver_type))
        case ("direct")
            call solve_direct_vector( &
                dense_matrix, right_hand_side, solution)
            stats%converged = maxval(abs( &
                matmul(dense_matrix, solution) - right_hand_side)) < 1.0e-10_dp
            stats%iterations = 1
            stats%final_residual = sqrt(sum( &
                (matmul(dense_matrix, solution) - right_hand_side)**2))
            stats%solve_time = 0.0_dp
            stats%memory_usage = 0
            stats%method_used = "lapack_lu"
            stats%restarts = 0
            stats%parallel_efficiency = 0.0_dp
            stats%condition_estimate = 0.0_dp
        case default
            call gmres_solve( &
                dense_matrix, right_hand_side, solution, options, stats)
        end select

        if (.not. allocated(field%values)) then
            allocate(field%values(dof_count, 2))
        end if
        field%values(:, 1) = solution
        field%values(:, 2) = 0.0_dp
    end subroutine solve_compiled_nedelec_problem

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
