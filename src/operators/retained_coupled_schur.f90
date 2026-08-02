module fortfem_retained_coupled_schur
    !! Coupled Schur reduction over retained per-field factors.
    !!
    !! For a block system
    !!
    !!   [ E  C ] [ y ] = [ f ]
    !!   [ D  F ] [ x ]   [ g ]
    !!
    !! with block-diagonal (or otherwise caller-composed) D, this module
    !! returns E-C D^{-1} F and f-C D^{-1} g.  D is represented by the
    !! retained field split contract; all off-diagonal blocks stay caller
    !! owned.  The value, JVP, and VJP paths are therefore reusable by HDG,
    !! FEM/BEM, DtN, PML, wall, and multi-field MHD clients.
    use fortfem_kinds, only: dp
    use fortfem_retained_field_split, only: &
        apply_retained_complex_field_split, &
        apply_retained_complex_field_split_jvp, &
        apply_retained_complex_field_split_vjp, &
        apply_retained_field_split, &
        apply_retained_field_split_jvp, &
        apply_retained_field_split_vjp, &
        retained_complex_field_split_t, retained_field_split_t
    use fortsparse, only: csc_t, csc_z_t, FORTSPARSE_INVALID_MATRIX, &
        FORTSPARSE_OK
    implicit none
    private

    public :: assemble_retained_coupled_schur
    public :: assemble_retained_coupled_schur_jvp
    public :: assemble_retained_coupled_schur_vjp

    interface assemble_retained_coupled_schur
        module procedure assemble_retained_coupled_schur_real
        module procedure assemble_retained_coupled_schur_complex
    end interface assemble_retained_coupled_schur

    interface assemble_retained_coupled_schur_jvp
        module procedure assemble_retained_coupled_schur_jvp_real
        module procedure assemble_retained_coupled_schur_jvp_complex
    end interface assemble_retained_coupled_schur_jvp

    interface assemble_retained_coupled_schur_vjp
        module procedure assemble_retained_coupled_schur_vjp_real
        module procedure assemble_retained_coupled_schur_vjp_complex
    end interface assemble_retained_coupled_schur_vjp

    interface outer_product
        module procedure outer_product_real
        module procedure outer_product_complex
    end interface outer_product

contains

    subroutine assemble_retained_coupled_schur_real( &
            split, exterior_matrix, exterior_to_fields, fields_to_exterior, &
            exterior_rhs, fields_rhs, schur_matrix, schur_rhs, status)
        type(retained_field_split_t), intent(inout) :: split
        real(dp), intent(in) :: exterior_matrix(:, :), exterior_to_fields(:, :)
        real(dp), intent(in) :: fields_to_exterior(:, :), exterior_rhs(:)
        real(dp), intent(in) :: fields_rhs(:)
        real(dp), intent(out) :: schur_matrix(:, :), schur_rhs(:)
        integer, intent(out) :: status

        real(dp), allocatable :: fields_solution(:), fields_response(:, :)
        integer :: exterior_count, fields_count, column, solve_status

        schur_matrix = 0.0_dp
        schur_rhs = 0.0_dp
        status = FORTSPARSE_INVALID_MATRIX
        if (.not. valid_real_inputs( &
            split, exterior_matrix, exterior_to_fields, fields_to_exterior, &
            exterior_rhs, fields_rhs, schur_matrix, schur_rhs)) return
        exterior_count = size(exterior_matrix, 1)
        fields_count = size(fields_rhs)
        allocate(fields_solution(fields_count), fields_response( &
            fields_count, exterior_count))
        call apply_retained_field_split( &
            split, fields_rhs, fields_solution, solve_status)
        if (solve_status /= FORTSPARSE_OK) then
            status = solve_status
            return
        end if
        do column = 1, exterior_count
            call apply_retained_field_split( &
                split, fields_to_exterior(:, column), fields_response(:, column), &
                solve_status)
            if (solve_status /= FORTSPARSE_OK) then
                status = solve_status
                return
            end if
        end do
        schur_matrix = exterior_matrix - matmul( &
            exterior_to_fields, fields_response)
        schur_rhs = exterior_rhs - matmul(exterior_to_fields, fields_solution)
        status = FORTSPARSE_OK
    end subroutine assemble_retained_coupled_schur_real

    subroutine assemble_retained_coupled_schur_complex( &
            split, exterior_matrix, exterior_to_fields, fields_to_exterior, &
            exterior_rhs, fields_rhs, schur_matrix, schur_rhs, status)
        type(retained_complex_field_split_t), intent(inout) :: split
        complex(dp), intent(in) :: exterior_matrix(:, :), exterior_to_fields(:, :)
        complex(dp), intent(in) :: fields_to_exterior(:, :), exterior_rhs(:)
        complex(dp), intent(in) :: fields_rhs(:)
        complex(dp), intent(out) :: schur_matrix(:, :), schur_rhs(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: fields_solution(:), fields_response(:, :)
        integer :: exterior_count, fields_count, column, solve_status

        schur_matrix = cmplx(0.0_dp, 0.0_dp, dp)
        schur_rhs = cmplx(0.0_dp, 0.0_dp, dp)
        status = FORTSPARSE_INVALID_MATRIX
        if (.not. valid_complex_inputs( &
            split, exterior_matrix, exterior_to_fields, fields_to_exterior, &
            exterior_rhs, fields_rhs, schur_matrix, schur_rhs)) return
        exterior_count = size(exterior_matrix, 1)
        fields_count = size(fields_rhs)
        allocate(fields_solution(fields_count), fields_response( &
            fields_count, exterior_count))
        call apply_retained_complex_field_split( &
            split, fields_rhs, fields_solution, solve_status)
        if (solve_status /= FORTSPARSE_OK) then
            status = solve_status
            return
        end if
        do column = 1, exterior_count
            call apply_retained_complex_field_split( &
                split, fields_to_exterior(:, column), fields_response(:, column), &
                solve_status)
            if (solve_status /= FORTSPARSE_OK) then
                status = solve_status
                return
            end if
        end do
        schur_matrix = exterior_matrix - matmul( &
            exterior_to_fields, fields_response)
        schur_rhs = exterior_rhs - matmul(exterior_to_fields, fields_solution)
        status = FORTSPARSE_OK
    end subroutine assemble_retained_coupled_schur_complex

    subroutine assemble_retained_coupled_schur_jvp_real( &
            split, exterior_matrix, exterior_to_fields, fields_to_exterior, &
            exterior_rhs, fields_rhs, matrices_dot, exterior_matrix_dot, &
            exterior_to_fields_dot, fields_to_exterior_dot, exterior_rhs_dot, &
            fields_rhs_dot, schur_matrix_dot, schur_rhs_dot, status)
        type(retained_field_split_t), intent(inout) :: split
        real(dp), intent(in) :: exterior_matrix(:, :), exterior_to_fields(:, :)
        real(dp), intent(in) :: fields_to_exterior(:, :), exterior_rhs(:)
        real(dp), intent(in) :: fields_rhs(:)
        type(csc_t), intent(in) :: matrices_dot(:)
        real(dp), intent(in) :: exterior_matrix_dot(:, :)
        real(dp), intent(in) :: exterior_to_fields_dot(:, :)
        real(dp), intent(in) :: fields_to_exterior_dot(:, :)
        real(dp), intent(in) :: exterior_rhs_dot(:), fields_rhs_dot(:)
        real(dp), intent(out) :: schur_matrix_dot(:, :), schur_rhs_dot(:)
        integer, intent(out) :: status

        real(dp), allocatable :: fields_solution(:), fields_solution_dot(:)
        real(dp), allocatable :: fields_response(:, :), fields_response_dot(:, :)
        integer :: exterior_count, fields_count, column, solve_status

        schur_matrix_dot = 0.0_dp
        schur_rhs_dot = 0.0_dp
        status = FORTSPARSE_INVALID_MATRIX
        if (.not. valid_real_inputs( &
            split, exterior_matrix, exterior_to_fields, fields_to_exterior, &
            exterior_rhs, fields_rhs, schur_matrix_dot, schur_rhs_dot)) return
        if (.not. valid_real_directions( &
            split, matrices_dot, exterior_matrix_dot, exterior_to_fields_dot, &
            fields_to_exterior_dot, exterior_rhs_dot, fields_rhs_dot, &
            schur_matrix_dot, schur_rhs_dot)) return
        exterior_count = size(exterior_matrix, 1)
        fields_count = size(fields_rhs)
        allocate(fields_solution(fields_count), fields_solution_dot(fields_count), &
            fields_response(fields_count, exterior_count), &
            fields_response_dot(fields_count, exterior_count))
        call apply_retained_field_split( &
            split, fields_rhs, fields_solution, solve_status)
        if (solve_status /= FORTSPARSE_OK) then
            status = solve_status
            return
        end if
        call apply_retained_field_split_jvp( &
            split, matrices_dot, fields_solution, fields_rhs_dot, &
            fields_solution_dot, solve_status)
        if (solve_status /= FORTSPARSE_OK) then
            status = solve_status
            return
        end if
        do column = 1, exterior_count
            call apply_retained_field_split( &
                split, fields_to_exterior(:, column), fields_response(:, column), &
                solve_status)
            if (solve_status /= FORTSPARSE_OK) then
                status = solve_status
                return
            end if
            call apply_retained_field_split_jvp( &
                split, matrices_dot, fields_response(:, column), &
                fields_to_exterior_dot(:, column), fields_response_dot(:, column), &
                solve_status)
            if (solve_status /= FORTSPARSE_OK) then
                status = solve_status
                return
            end if
        end do
        schur_matrix_dot = exterior_matrix_dot - matmul( &
            exterior_to_fields_dot, fields_response) - matmul( &
            exterior_to_fields, fields_response_dot)
        schur_rhs_dot = exterior_rhs_dot - matmul( &
            exterior_to_fields_dot, fields_solution) - matmul( &
            exterior_to_fields, fields_solution_dot)
        status = FORTSPARSE_OK
    end subroutine assemble_retained_coupled_schur_jvp_real

    subroutine assemble_retained_coupled_schur_jvp_complex( &
            split, exterior_matrix, exterior_to_fields, fields_to_exterior, &
            exterior_rhs, fields_rhs, matrices_dot, exterior_matrix_dot, &
            exterior_to_fields_dot, fields_to_exterior_dot, exterior_rhs_dot, &
            fields_rhs_dot, schur_matrix_dot, schur_rhs_dot, status)
        type(retained_complex_field_split_t), intent(inout) :: split
        complex(dp), intent(in) :: exterior_matrix(:, :), exterior_to_fields(:, :)
        complex(dp), intent(in) :: fields_to_exterior(:, :), exterior_rhs(:)
        complex(dp), intent(in) :: fields_rhs(:)
        type(csc_z_t), intent(in) :: matrices_dot(:)
        complex(dp), intent(in) :: exterior_matrix_dot(:, :)
        complex(dp), intent(in) :: exterior_to_fields_dot(:, :)
        complex(dp), intent(in) :: fields_to_exterior_dot(:, :)
        complex(dp), intent(in) :: exterior_rhs_dot(:), fields_rhs_dot(:)
        complex(dp), intent(out) :: schur_matrix_dot(:, :), schur_rhs_dot(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: fields_solution(:), fields_solution_dot(:)
        complex(dp), allocatable :: fields_response(:, :), fields_response_dot(:, :)
        integer :: exterior_count, fields_count, column, solve_status

        schur_matrix_dot = cmplx(0.0_dp, 0.0_dp, dp)
        schur_rhs_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = FORTSPARSE_INVALID_MATRIX
        if (.not. valid_complex_inputs( &
            split, exterior_matrix, exterior_to_fields, fields_to_exterior, &
            exterior_rhs, fields_rhs, schur_matrix_dot, schur_rhs_dot)) return
        if (.not. valid_complex_directions( &
            split, matrices_dot, exterior_matrix_dot, exterior_to_fields_dot, &
            fields_to_exterior_dot, exterior_rhs_dot, fields_rhs_dot, &
            schur_matrix_dot, schur_rhs_dot)) return
        exterior_count = size(exterior_matrix, 1)
        fields_count = size(fields_rhs)
        allocate(fields_solution(fields_count), fields_solution_dot(fields_count), &
            fields_response(fields_count, exterior_count), &
            fields_response_dot(fields_count, exterior_count))
        call apply_retained_complex_field_split( &
            split, fields_rhs, fields_solution, solve_status)
        if (solve_status /= FORTSPARSE_OK) then
            status = solve_status
            return
        end if
        call apply_retained_complex_field_split_jvp( &
            split, matrices_dot, fields_solution, fields_rhs_dot, &
            fields_solution_dot, solve_status)
        if (solve_status /= FORTSPARSE_OK) then
            status = solve_status
            return
        end if
        do column = 1, exterior_count
            call apply_retained_complex_field_split( &
                split, fields_to_exterior(:, column), fields_response(:, column), &
                solve_status)
            if (solve_status /= FORTSPARSE_OK) then
                status = solve_status
                return
            end if
            call apply_retained_complex_field_split_jvp( &
                split, matrices_dot, fields_response(:, column), &
                fields_to_exterior_dot(:, column), fields_response_dot(:, column), &
                solve_status)
            if (solve_status /= FORTSPARSE_OK) then
                status = solve_status
                return
            end if
        end do
        schur_matrix_dot = exterior_matrix_dot - matmul( &
            exterior_to_fields_dot, fields_response) - matmul( &
            exterior_to_fields, fields_response_dot)
        schur_rhs_dot = exterior_rhs_dot - matmul( &
            exterior_to_fields_dot, fields_solution) - matmul( &
            exterior_to_fields, fields_solution_dot)
        status = FORTSPARSE_OK
    end subroutine assemble_retained_coupled_schur_jvp_complex

    subroutine assemble_retained_coupled_schur_vjp_real( &
            split, exterior_matrix, exterior_to_fields, fields_to_exterior, &
            exterior_rhs, fields_rhs, schur_matrix_bar, schur_rhs_bar, &
            exterior_matrix_bar, exterior_to_fields_bar, fields_to_exterior_bar, &
            exterior_rhs_bar, fields_rhs_bar, matrices_bar, status)
        type(retained_field_split_t), intent(inout) :: split
        real(dp), intent(in) :: exterior_matrix(:, :), exterior_to_fields(:, :)
        real(dp), intent(in) :: fields_to_exterior(:, :), exterior_rhs(:)
        real(dp), intent(in) :: fields_rhs(:), schur_matrix_bar(:, :)
        real(dp), intent(in) :: schur_rhs_bar(:)
        real(dp), intent(out) :: exterior_matrix_bar(:, :)
        real(dp), intent(out) :: exterior_to_fields_bar(:, :)
        real(dp), intent(out) :: fields_to_exterior_bar(:, :)
        real(dp), intent(out) :: exterior_rhs_bar(:), fields_rhs_bar(:)
        type(csc_t), intent(out) :: matrices_bar(:)
        integer, intent(out) :: status

        real(dp), allocatable :: fields_solution(:), fields_response(:, :)
        real(dp), allocatable :: fields_solution_bar(:), fields_response_bar(:, :)
        real(dp), allocatable :: rhs_bar(:)
        type(csc_t), allocatable :: temporary_matrices_bar(:)
        integer :: exterior_count, fields_count, column, solve_status

        exterior_matrix_bar = 0.0_dp
        exterior_to_fields_bar = 0.0_dp
        fields_to_exterior_bar = 0.0_dp
        exterior_rhs_bar = 0.0_dp
        fields_rhs_bar = 0.0_dp
        status = FORTSPARSE_INVALID_MATRIX
        if (.not. valid_real_inputs( &
            split, exterior_matrix, exterior_to_fields, fields_to_exterior, &
            exterior_rhs, fields_rhs, schur_matrix_bar, schur_rhs_bar)) return
        if (.not. valid_real_cotangents( &
            split, schur_matrix_bar, schur_rhs_bar, exterior_matrix_bar, &
            exterior_to_fields_bar, fields_to_exterior_bar, exterior_rhs_bar, &
            fields_rhs_bar, matrices_bar)) return
        exterior_count = size(exterior_matrix, 1)
        fields_count = size(fields_rhs)
        allocate(fields_solution(fields_count), fields_response( &
            fields_count, exterior_count), fields_solution_bar(fields_count), &
            fields_response_bar(fields_count, exterior_count), rhs_bar(fields_count), &
            temporary_matrices_bar(size(matrices_bar)))
        call apply_retained_field_split( &
            split, fields_rhs, fields_solution, solve_status)
        if (solve_status /= FORTSPARSE_OK) then
            status = solve_status
            return
        end if
        do column = 1, exterior_count
            call apply_retained_field_split( &
                split, fields_to_exterior(:, column), fields_response(:, column), &
                solve_status)
            if (solve_status /= FORTSPARSE_OK) then
                status = solve_status
                return
            end if
        end do
        exterior_matrix_bar = schur_matrix_bar
        exterior_rhs_bar = schur_rhs_bar
        exterior_to_fields_bar = -matmul( &
            schur_matrix_bar, transpose(fields_response)) - outer_product( &
            schur_rhs_bar, fields_solution)
        fields_solution_bar = -matmul( &
            transpose(exterior_to_fields), schur_rhs_bar)
        fields_response_bar = -matmul( &
            transpose(exterior_to_fields), schur_matrix_bar)
        call accumulate_real_split_vjp( &
            split, fields_solution, fields_solution_bar, fields_rhs_bar, &
            matrices_bar, temporary_matrices_bar, status)
        if (status /= FORTSPARSE_OK) return
        do column = 1, exterior_count
            call apply_retained_field_split_vjp( &
                split, fields_response(:, column), fields_response_bar(:, column), &
                rhs_bar, temporary_matrices_bar, solve_status)
            if (solve_status /= FORTSPARSE_OK) then
                status = solve_status
                return
            end if
            fields_to_exterior_bar(:, column) = rhs_bar
            call add_real_matrix_bars(matrices_bar, temporary_matrices_bar)
        end do
        status = FORTSPARSE_OK
    end subroutine assemble_retained_coupled_schur_vjp_real

    subroutine assemble_retained_coupled_schur_vjp_complex( &
            split, exterior_matrix, exterior_to_fields, fields_to_exterior, &
            exterior_rhs, fields_rhs, schur_matrix_bar, schur_rhs_bar, &
            exterior_matrix_bar, exterior_to_fields_bar, fields_to_exterior_bar, &
            exterior_rhs_bar, fields_rhs_bar, matrices_bar, status)
        type(retained_complex_field_split_t), intent(inout) :: split
        complex(dp), intent(in) :: exterior_matrix(:, :), exterior_to_fields(:, :)
        complex(dp), intent(in) :: fields_to_exterior(:, :), exterior_rhs(:)
        complex(dp), intent(in) :: fields_rhs(:), schur_matrix_bar(:, :)
        complex(dp), intent(in) :: schur_rhs_bar(:)
        complex(dp), intent(out) :: exterior_matrix_bar(:, :)
        complex(dp), intent(out) :: exterior_to_fields_bar(:, :)
        complex(dp), intent(out) :: fields_to_exterior_bar(:, :)
        complex(dp), intent(out) :: exterior_rhs_bar(:), fields_rhs_bar(:)
        type(csc_z_t), intent(out) :: matrices_bar(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: fields_solution(:), fields_response(:, :)
        complex(dp), allocatable :: fields_solution_bar(:), fields_response_bar(:, :)
        complex(dp), allocatable :: rhs_bar(:)
        type(csc_z_t), allocatable :: temporary_matrices_bar(:)
        integer :: exterior_count, fields_count, column, solve_status

        exterior_matrix_bar = cmplx(0.0_dp, 0.0_dp, dp)
        exterior_to_fields_bar = cmplx(0.0_dp, 0.0_dp, dp)
        fields_to_exterior_bar = cmplx(0.0_dp, 0.0_dp, dp)
        exterior_rhs_bar = cmplx(0.0_dp, 0.0_dp, dp)
        fields_rhs_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = FORTSPARSE_INVALID_MATRIX
        if (.not. valid_complex_inputs( &
            split, exterior_matrix, exterior_to_fields, fields_to_exterior, &
            exterior_rhs, fields_rhs, schur_matrix_bar, schur_rhs_bar)) return
        if (.not. valid_complex_cotangents( &
            split, schur_matrix_bar, schur_rhs_bar, exterior_matrix_bar, &
            exterior_to_fields_bar, fields_to_exterior_bar, exterior_rhs_bar, &
            fields_rhs_bar, matrices_bar)) return
        exterior_count = size(exterior_matrix, 1)
        fields_count = size(fields_rhs)
        allocate(fields_solution(fields_count), fields_response( &
            fields_count, exterior_count), fields_solution_bar(fields_count), &
            fields_response_bar(fields_count, exterior_count), rhs_bar(fields_count), &
            temporary_matrices_bar(size(matrices_bar)))
        call apply_retained_complex_field_split( &
            split, fields_rhs, fields_solution, solve_status)
        if (solve_status /= FORTSPARSE_OK) then
            status = solve_status
            return
        end if
        do column = 1, exterior_count
            call apply_retained_complex_field_split( &
                split, fields_to_exterior(:, column), fields_response(:, column), &
                solve_status)
            if (solve_status /= FORTSPARSE_OK) then
                status = solve_status
                return
            end if
        end do
        exterior_matrix_bar = schur_matrix_bar
        exterior_rhs_bar = schur_rhs_bar
        exterior_to_fields_bar = -matmul( &
            schur_matrix_bar, conjg(transpose(fields_response))) - outer_product( &
            schur_rhs_bar, conjg(fields_solution))
        fields_solution_bar = -matmul( &
            conjg(transpose(exterior_to_fields)), schur_rhs_bar)
        fields_response_bar = -matmul( &
            conjg(transpose(exterior_to_fields)), schur_matrix_bar)
        call accumulate_complex_split_vjp( &
            split, fields_solution, fields_solution_bar, fields_rhs_bar, &
            matrices_bar, temporary_matrices_bar, status)
        if (status /= FORTSPARSE_OK) return
        do column = 1, exterior_count
            call apply_retained_complex_field_split_vjp( &
                split, fields_response(:, column), fields_response_bar(:, column), &
                rhs_bar, temporary_matrices_bar, solve_status)
            if (solve_status /= FORTSPARSE_OK) then
                status = solve_status
                return
            end if
            fields_to_exterior_bar(:, column) = rhs_bar
            call add_complex_matrix_bars(matrices_bar, temporary_matrices_bar)
        end do
        status = FORTSPARSE_OK
    end subroutine assemble_retained_coupled_schur_vjp_complex

    subroutine accumulate_real_split_vjp( &
            split, solution, solution_bar, rhs_bar, total_bar, temporary_bar, status)
        type(retained_field_split_t), intent(inout) :: split
        real(dp), intent(in) :: solution(:), solution_bar(:)
        real(dp), intent(out) :: rhs_bar(:)
        type(csc_t), intent(inout) :: total_bar(:), temporary_bar(:)
        integer, intent(out) :: status

        call apply_retained_field_split_vjp( &
            split, solution, solution_bar, rhs_bar, temporary_bar, status)
        if (status /= FORTSPARSE_OK) return
        call initialize_real_matrix_bars(total_bar, split)
        call add_real_matrix_bars(total_bar, temporary_bar)
    end subroutine accumulate_real_split_vjp

    subroutine accumulate_complex_split_vjp( &
            split, solution, solution_bar, rhs_bar, total_bar, temporary_bar, status)
        type(retained_complex_field_split_t), intent(inout) :: split
        complex(dp), intent(in) :: solution(:), solution_bar(:)
        complex(dp), intent(out) :: rhs_bar(:)
        type(csc_z_t), intent(inout) :: total_bar(:), temporary_bar(:)
        integer, intent(out) :: status

        call apply_retained_complex_field_split_vjp( &
            split, solution, solution_bar, rhs_bar, temporary_bar, status)
        if (status /= FORTSPARSE_OK) return
        call initialize_complex_matrix_bars(total_bar, split)
        call add_complex_matrix_bars(total_bar, temporary_bar)
    end subroutine accumulate_complex_split_vjp

    subroutine initialize_real_matrix_bars(bars, split)
        type(csc_t), intent(out) :: bars(:)
        type(retained_field_split_t), intent(in) :: split
        integer :: field

        do field = 1, size(bars)
            bars(field) = split%matrices(field)
            bars(field)%val = 0.0_dp
        end do
    end subroutine initialize_real_matrix_bars

    subroutine initialize_complex_matrix_bars(bars, split)
        type(csc_z_t), intent(out) :: bars(:)
        type(retained_complex_field_split_t), intent(in) :: split
        integer :: field

        do field = 1, size(bars)
            bars(field) = split%matrices(field)
            bars(field)%val = cmplx(0.0_dp, 0.0_dp, dp)
        end do
    end subroutine initialize_complex_matrix_bars

    subroutine add_real_matrix_bars(total, increment)
        type(csc_t), intent(inout) :: total(:)
        type(csc_t), intent(in) :: increment(:)
        integer :: field

        do field = 1, size(total)
            total(field)%val = total(field)%val + increment(field)%val
        end do
    end subroutine add_real_matrix_bars

    subroutine add_complex_matrix_bars(total, increment)
        type(csc_z_t), intent(inout) :: total(:)
        type(csc_z_t), intent(in) :: increment(:)
        integer :: field

        do field = 1, size(total)
            total(field)%val = total(field)%val + increment(field)%val
        end do
    end subroutine add_complex_matrix_bars

    logical function valid_real_inputs( &
            split, exterior_matrix, exterior_to_fields, fields_to_exterior, &
            exterior_rhs, fields_rhs, schur_matrix, schur_rhs) result(valid)
        type(retained_field_split_t), intent(in) :: split
        real(dp), intent(in) :: exterior_matrix(:, :), exterior_to_fields(:, :)
        real(dp), intent(in) :: fields_to_exterior(:, :), exterior_rhs(:)
        real(dp), intent(in) :: fields_rhs(:), schur_matrix(:, :), schur_rhs(:)
        integer :: exterior_count, fields_count

        valid = valid_real_split_shape(split, fields_count)
        if (.not. valid) return
        exterior_count = size(exterior_matrix, 1)
        if (exterior_count < 1 .or. size(exterior_matrix, 2) /= exterior_count) then
            valid = .false.
            return
        end if
        if (size(exterior_to_fields, 1) /= exterior_count .or. &
            size(exterior_to_fields, 2) /= fields_count .or. &
            size(fields_to_exterior, 1) /= fields_count .or. &
            size(fields_to_exterior, 2) /= exterior_count .or. &
            size(exterior_rhs) /= exterior_count .or. size(fields_rhs) /= fields_count .or. &
            size(schur_matrix, 1) /= exterior_count .or. &
            size(schur_matrix, 2) /= exterior_count .or. &
            size(schur_rhs) /= exterior_count) valid = .false.
    end function valid_real_inputs

    logical function valid_complex_inputs( &
            split, exterior_matrix, exterior_to_fields, fields_to_exterior, &
            exterior_rhs, fields_rhs, schur_matrix, schur_rhs) result(valid)
        type(retained_complex_field_split_t), intent(in) :: split
        complex(dp), intent(in) :: exterior_matrix(:, :), exterior_to_fields(:, :)
        complex(dp), intent(in) :: fields_to_exterior(:, :), exterior_rhs(:)
        complex(dp), intent(in) :: fields_rhs(:), schur_matrix(:, :), schur_rhs(:)
        integer :: exterior_count, fields_count

        valid = valid_complex_split_shape(split, fields_count)
        if (.not. valid) return
        exterior_count = size(exterior_matrix, 1)
        if (exterior_count < 1 .or. size(exterior_matrix, 2) /= exterior_count) then
            valid = .false.
            return
        end if
        if (size(exterior_to_fields, 1) /= exterior_count .or. &
            size(exterior_to_fields, 2) /= fields_count .or. &
            size(fields_to_exterior, 1) /= fields_count .or. &
            size(fields_to_exterior, 2) /= exterior_count .or. &
            size(exterior_rhs) /= exterior_count .or. size(fields_rhs) /= fields_count .or. &
            size(schur_matrix, 1) /= exterior_count .or. &
            size(schur_matrix, 2) /= exterior_count .or. &
            size(schur_rhs) /= exterior_count) valid = .false.
    end function valid_complex_inputs

    logical function valid_real_directions( &
            split, matrices_dot, exterior_matrix_dot, exterior_to_fields_dot, &
            fields_to_exterior_dot, exterior_rhs_dot, fields_rhs_dot, &
            schur_matrix_dot, schur_rhs_dot) result(valid)
        type(retained_field_split_t), intent(in) :: split
        type(csc_t), intent(in) :: matrices_dot(:)
        real(dp), intent(in) :: exterior_matrix_dot(:, :), exterior_to_fields_dot(:, :)
        real(dp), intent(in) :: fields_to_exterior_dot(:, :), exterior_rhs_dot(:)
        real(dp), intent(in) :: fields_rhs_dot(:), schur_matrix_dot(:, :)
        real(dp), intent(in) :: schur_rhs_dot(:)
        integer :: exterior_count, fields_count

        valid = valid_real_split_direction(split, matrices_dot, fields_count)
        if (.not. valid) return
        exterior_count = size(schur_matrix_dot, 1)
        if (size(exterior_matrix_dot, 1) /= exterior_count .or. &
            size(exterior_matrix_dot, 2) /= exterior_count .or. &
            size(exterior_to_fields_dot, 1) /= exterior_count .or. &
            size(exterior_to_fields_dot, 2) /= fields_count .or. &
            size(fields_to_exterior_dot, 1) /= fields_count .or. &
            size(fields_to_exterior_dot, 2) /= exterior_count .or. &
            size(exterior_rhs_dot) /= exterior_count .or. &
            size(fields_rhs_dot) /= fields_count .or. &
            size(schur_rhs_dot) /= exterior_count) valid = .false.
    end function valid_real_directions

    logical function valid_complex_directions( &
            split, matrices_dot, exterior_matrix_dot, exterior_to_fields_dot, &
            fields_to_exterior_dot, exterior_rhs_dot, fields_rhs_dot, &
            schur_matrix_dot, schur_rhs_dot) result(valid)
        type(retained_complex_field_split_t), intent(in) :: split
        type(csc_z_t), intent(in) :: matrices_dot(:)
        complex(dp), intent(in) :: exterior_matrix_dot(:, :), exterior_to_fields_dot(:, :)
        complex(dp), intent(in) :: fields_to_exterior_dot(:, :), exterior_rhs_dot(:)
        complex(dp), intent(in) :: fields_rhs_dot(:), schur_matrix_dot(:, :)
        complex(dp), intent(in) :: schur_rhs_dot(:)
        integer :: exterior_count, fields_count

        valid = valid_complex_split_direction(split, matrices_dot, fields_count)
        if (.not. valid) return
        exterior_count = size(schur_matrix_dot, 1)
        if (size(exterior_matrix_dot, 1) /= exterior_count .or. &
            size(exterior_matrix_dot, 2) /= exterior_count .or. &
            size(exterior_to_fields_dot, 1) /= exterior_count .or. &
            size(exterior_to_fields_dot, 2) /= fields_count .or. &
            size(fields_to_exterior_dot, 1) /= fields_count .or. &
            size(fields_to_exterior_dot, 2) /= exterior_count .or. &
            size(exterior_rhs_dot) /= exterior_count .or. &
            size(fields_rhs_dot) /= fields_count .or. &
            size(schur_rhs_dot) /= exterior_count) valid = .false.
    end function valid_complex_directions

    logical function valid_real_cotangents( &
            split, schur_matrix_bar, schur_rhs_bar, exterior_matrix_bar, &
            exterior_to_fields_bar, fields_to_exterior_bar, exterior_rhs_bar, &
            fields_rhs_bar, matrices_bar) result(valid)
        type(retained_field_split_t), intent(in) :: split
        real(dp), intent(in) :: schur_matrix_bar(:, :), schur_rhs_bar(:)
        real(dp), intent(in) :: exterior_matrix_bar(:, :), exterior_to_fields_bar(:, :)
        real(dp), intent(in) :: fields_to_exterior_bar(:, :), exterior_rhs_bar(:)
        real(dp), intent(in) :: fields_rhs_bar(:)
        type(csc_t), intent(in) :: matrices_bar(:)
        integer :: exterior_count, fields_count

        valid = valid_real_split_shape(split, fields_count)
        if (.not. valid) return
        exterior_count = size(schur_matrix_bar, 1)
        valid = size(schur_matrix_bar, 2) == exterior_count .and. &
            size(schur_rhs_bar) == exterior_count .and. &
            size(exterior_matrix_bar, 1) == exterior_count .and. &
            size(exterior_matrix_bar, 2) == exterior_count .and. &
            size(exterior_to_fields_bar, 1) == exterior_count .and. &
            size(exterior_to_fields_bar, 2) == fields_count .and. &
            size(fields_to_exterior_bar, 1) == fields_count .and. &
            size(fields_to_exterior_bar, 2) == exterior_count .and. &
            size(exterior_rhs_bar) == exterior_count .and. &
            size(fields_rhs_bar) == fields_count .and. &
            size(matrices_bar) == size(split%matrices)
    end function valid_real_cotangents

    logical function valid_complex_cotangents( &
            split, schur_matrix_bar, schur_rhs_bar, exterior_matrix_bar, &
            exterior_to_fields_bar, fields_to_exterior_bar, exterior_rhs_bar, &
            fields_rhs_bar, matrices_bar) result(valid)
        type(retained_complex_field_split_t), intent(in) :: split
        complex(dp), intent(in) :: schur_matrix_bar(:, :), schur_rhs_bar(:)
        complex(dp), intent(in) :: exterior_matrix_bar(:, :), exterior_to_fields_bar(:, :)
        complex(dp), intent(in) :: fields_to_exterior_bar(:, :), exterior_rhs_bar(:)
        complex(dp), intent(in) :: fields_rhs_bar(:)
        type(csc_z_t), intent(in) :: matrices_bar(:)
        integer :: exterior_count, fields_count

        valid = valid_complex_split_shape(split, fields_count)
        if (.not. valid) return
        exterior_count = size(schur_matrix_bar, 1)
        valid = size(schur_matrix_bar, 2) == exterior_count .and. &
            size(schur_rhs_bar) == exterior_count .and. &
            size(exterior_matrix_bar, 1) == exterior_count .and. &
            size(exterior_matrix_bar, 2) == exterior_count .and. &
            size(exterior_to_fields_bar, 1) == exterior_count .and. &
            size(exterior_to_fields_bar, 2) == fields_count .and. &
            size(fields_to_exterior_bar, 1) == fields_count .and. &
            size(fields_to_exterior_bar, 2) == exterior_count .and. &
            size(exterior_rhs_bar) == exterior_count .and. &
            size(fields_rhs_bar) == fields_count .and. &
            size(matrices_bar) == size(split%matrices)
    end function valid_complex_cotangents

    logical function valid_real_split_shape(split, total) result(valid)
        type(retained_field_split_t), intent(in) :: split
        integer, intent(out) :: total
        integer :: field

        total = 0
        valid = allocated(split%matrices) .and. allocated(split%factors) .and. &
            allocated(split%transpose_factors)
        if (.not. valid) return
        if (split%factored_count /= size(split%matrices) .or. &
            size(split%factors) /= size(split%matrices)) then
            valid = .false.
            return
        end if
        do field = 1, size(split%matrices)
            if (split%matrices(field)%nrow < 1 .or. &
                split%matrices(field)%nrow /= split%matrices(field)%ncol) then
                valid = .false.
                return
            end if
            total = total + split%matrices(field)%nrow
        end do
        valid = total > 0
    end function valid_real_split_shape

    logical function valid_complex_split_shape(split, total) result(valid)
        type(retained_complex_field_split_t), intent(in) :: split
        integer, intent(out) :: total
        integer :: field

        total = 0
        valid = allocated(split%matrices) .and. allocated(split%factors) .and. &
            allocated(split%adjoint_factors)
        if (.not. valid) return
        if (split%factored_count /= size(split%matrices) .or. &
            size(split%factors) /= size(split%matrices)) then
            valid = .false.
            return
        end if
        do field = 1, size(split%matrices)
            if (split%matrices(field)%nrow < 1 .or. &
                split%matrices(field)%nrow /= split%matrices(field)%ncol) then
                valid = .false.
                return
            end if
            total = total + split%matrices(field)%nrow
        end do
        valid = total > 0
    end function valid_complex_split_shape

    logical function valid_real_split_direction( &
            split, matrices_dot, total) result(valid)
        type(retained_field_split_t), intent(in) :: split
        type(csc_t), intent(in) :: matrices_dot(:)
        integer, intent(out) :: total
        integer :: field

        total = 0
        valid = size(matrices_dot) == size(split%matrices)
        if (.not. valid) return
        do field = 1, size(matrices_dot)
            if (.not. allocated(matrices_dot(field)%col_ptr) .or. &
                .not. allocated(matrices_dot(field)%row_idx) .or. &
                .not. allocated(matrices_dot(field)%val)) then
                valid = .false.
                return
            end if
            if (matrices_dot(field)%nrow /= split%matrices(field)%nrow .or. &
                matrices_dot(field)%ncol /= split%matrices(field)%ncol .or. &
                matrices_dot(field)%nnz /= split%matrices(field)%nnz .or. &
                size(matrices_dot(field)%col_ptr) /= matrices_dot(field)%ncol + 1 .or. &
                size(matrices_dot(field)%row_idx) /= matrices_dot(field)%nnz .or. &
                size(matrices_dot(field)%val) /= matrices_dot(field)%nnz) then
                valid = .false.
                return
            end if
            if (.not. all(matrices_dot(field)%col_ptr == &
                split%matrices(field)%col_ptr) .or. &
                .not. all(matrices_dot(field)%row_idx == &
                split%matrices(field)%row_idx)) then
                valid = .false.
                return
            end if
            total = total + matrices_dot(field)%nrow
        end do
    end function valid_real_split_direction

    logical function valid_complex_split_direction( &
            split, matrices_dot, total) result(valid)
        type(retained_complex_field_split_t), intent(in) :: split
        type(csc_z_t), intent(in) :: matrices_dot(:)
        integer, intent(out) :: total
        integer :: field

        total = 0
        valid = size(matrices_dot) == size(split%matrices)
        if (.not. valid) return
        do field = 1, size(matrices_dot)
            if (.not. allocated(matrices_dot(field)%col_ptr) .or. &
                .not. allocated(matrices_dot(field)%row_idx) .or. &
                .not. allocated(matrices_dot(field)%val)) then
                valid = .false.
                return
            end if
            if (matrices_dot(field)%nrow /= split%matrices(field)%nrow .or. &
                matrices_dot(field)%ncol /= split%matrices(field)%ncol .or. &
                matrices_dot(field)%nnz /= split%matrices(field)%nnz .or. &
                size(matrices_dot(field)%col_ptr) /= matrices_dot(field)%ncol + 1 .or. &
                size(matrices_dot(field)%row_idx) /= matrices_dot(field)%nnz .or. &
                size(matrices_dot(field)%val) /= matrices_dot(field)%nnz) then
                valid = .false.
                return
            end if
            if (.not. all(matrices_dot(field)%col_ptr == &
                split%matrices(field)%col_ptr) .or. &
                .not. all(matrices_dot(field)%row_idx == &
                split%matrices(field)%row_idx)) then
                valid = .false.
                return
            end if
            total = total + matrices_dot(field)%nrow
        end do
    end function valid_complex_split_direction

    function outer_product_real(left, right) result(product)
        real(dp), intent(in) :: left(:), right(:)
        real(dp) :: product(size(left), size(right))
        product = spread(left, dim=2, ncopies=size(right))* &
            spread(right, dim=1, ncopies=size(left))
    end function outer_product_real

    function outer_product_complex(left, right) result(product)
        complex(dp), intent(in) :: left(:), right(:)
        complex(dp) :: product(size(left), size(right))
        product = spread(left, dim=2, ncopies=size(right))* &
            spread(right, dim=1, ncopies=size(left))
    end function outer_product_complex

end module fortfem_retained_coupled_schur
