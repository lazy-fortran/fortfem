module fortfem_hdg_static_condensation
    !! Differentiable local Schur complement for HDG and hybridized blocks.
    !!
    !! For a local block system
    !! [ A B ] [x] = [f], the interior unknown x is eliminated and the
    !! [ C D ] [lambda]   [g]
    !! trace system S lambda = r is returned with S=D-C A^-1 B and
    !! r=g-C A^-1 f.  The operation is neutral: FEEC trace spaces, physical
    !! fluxes, and global skeleton assembly remain caller-owned.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, FORTSPARSE_SINGULAR
    implicit none
    private

    public :: assemble_hdg_static_condensation
    public :: assemble_hdg_static_condensation_jvp
    public :: assemble_hdg_static_condensation_vjp

contains

    subroutine assemble_hdg_static_condensation( &
            interior_matrix, interior_to_trace, trace_to_interior, trace_matrix, &
            interior_rhs, trace_rhs, condensed_matrix, condensed_rhs, status)
        !! Form the local HDG Schur complement and condensed right-hand side.
        real(dp), intent(in) :: interior_matrix(:, :), interior_to_trace(:, :)
        real(dp), intent(in) :: trace_to_interior(:, :), trace_matrix(:, :)
        real(dp), intent(in) :: interior_rhs(:), trace_rhs(:)
        real(dp), intent(out) :: condensed_matrix(:, :), condensed_rhs(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: interior_count, trace_count, info
        real(dp), allocatable :: interior_to_trace_solve(:, :), interior_rhs_solve(:)

        condensed_matrix = 0.0_dp
        condensed_rhs = 0.0_dp
        call validate_hdg_inputs( &
            interior_matrix, interior_to_trace, trace_to_interior, trace_matrix, &
            interior_rhs, trace_rhs, condensed_matrix, condensed_rhs, status)
        if (status%code /= FORTSPARSE_OK) return
        interior_count = size(interior_matrix, 1)
        trace_count = size(trace_matrix, 1)
        allocate(interior_to_trace_solve(interior_count, trace_count), &
            interior_rhs_solve(interior_count))
        call dense_solve(interior_matrix, interior_to_trace, &
            interior_to_trace_solve, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "HDG static condensation interior block is singular")
            return
        end if
        call dense_solve(interior_matrix, interior_rhs, interior_rhs_solve, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "HDG static condensation interior block is singular")
            return
        end if
        condensed_matrix = trace_matrix - matmul( &
            trace_to_interior, interior_to_trace_solve)
        condensed_rhs = trace_rhs - matmul(trace_to_interior, interior_rhs_solve)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_hdg_static_condensation

    subroutine assemble_hdg_static_condensation_jvp( &
            interior_matrix, interior_to_trace, trace_to_interior, trace_matrix, &
            interior_rhs, trace_rhs, interior_matrix_dot, interior_to_trace_dot, &
            trace_to_interior_dot, trace_matrix_dot, interior_rhs_dot, trace_rhs_dot, &
            condensed_matrix_dot, condensed_rhs_dot, status)
        !! Apply the implicit product-rule JVP of the local Schur complement.
        real(dp), intent(in) :: interior_matrix(:, :), interior_to_trace(:, :)
        real(dp), intent(in) :: trace_to_interior(:, :), trace_matrix(:, :)
        real(dp), intent(in) :: interior_rhs(:), trace_rhs(:)
        real(dp), intent(in) :: interior_matrix_dot(:, :), interior_to_trace_dot(:, :)
        real(dp), intent(in) :: trace_to_interior_dot(:, :), trace_matrix_dot(:, :)
        real(dp), intent(in) :: interior_rhs_dot(:), trace_rhs_dot(:)
        real(dp), intent(out) :: condensed_matrix_dot(:, :), condensed_rhs_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: interior_count, trace_count, info
        real(dp), allocatable :: interior_to_trace_solve(:, :), interior_rhs_solve(:)
        real(dp), allocatable :: interior_to_trace_solve_dot(:, :), interior_rhs_solve_dot(:)

        condensed_matrix_dot = 0.0_dp
        condensed_rhs_dot = 0.0_dp
        call validate_hdg_inputs( &
            interior_matrix, interior_to_trace, trace_to_interior, trace_matrix, &
            interior_rhs, trace_rhs, condensed_matrix_dot, condensed_rhs_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        interior_count = size(interior_matrix, 1)
        trace_count = size(trace_matrix, 1)
        if (.not. valid_hdg_direction( &
            interior_matrix_dot, interior_to_trace_dot, trace_to_interior_dot, &
            trace_matrix_dot, interior_rhs_dot, trace_rhs_dot, interior_count, &
            trace_count)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "HDG static condensation JVP has incompatible increments")
            return
        end if
        allocate(interior_to_trace_solve(interior_count, trace_count), &
            interior_rhs_solve(interior_count), interior_to_trace_solve_dot( &
            interior_count, trace_count), interior_rhs_solve_dot(interior_count))
        call dense_solve(interior_matrix, interior_to_trace, &
            interior_to_trace_solve, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "HDG static condensation interior block is singular")
            return
        end if
        call dense_solve(interior_matrix, interior_rhs, interior_rhs_solve, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "HDG static condensation interior block is singular")
            return
        end if
        call dense_solve(interior_matrix, interior_to_trace_dot - matmul( &
            interior_matrix_dot, interior_to_trace_solve), interior_to_trace_solve_dot, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "HDG static condensation interior block is singular")
            return
        end if
        call dense_solve(interior_matrix, interior_rhs_dot - matmul( &
            interior_matrix_dot, interior_rhs_solve), interior_rhs_solve_dot, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "HDG static condensation interior block is singular")
            return
        end if
        condensed_matrix_dot = trace_matrix_dot - matmul(trace_to_interior_dot, &
            interior_to_trace_solve) - matmul(trace_to_interior, interior_to_trace_solve_dot)
        condensed_rhs_dot = trace_rhs_dot - matmul(trace_to_interior_dot, interior_rhs_solve) - &
            matmul(trace_to_interior, interior_rhs_solve_dot)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_hdg_static_condensation_jvp

    subroutine assemble_hdg_static_condensation_vjp( &
            interior_matrix, interior_to_trace, trace_to_interior, trace_matrix, &
            interior_rhs, trace_rhs, condensed_matrix_bar, condensed_rhs_bar, &
            interior_matrix_bar, interior_to_trace_bar, trace_to_interior_bar, &
            trace_matrix_bar, interior_rhs_bar, trace_rhs_bar, status)
        !! Apply the reverse product of the local Schur complement.
        real(dp), intent(in) :: interior_matrix(:, :), interior_to_trace(:, :)
        real(dp), intent(in) :: trace_to_interior(:, :), trace_matrix(:, :)
        real(dp), intent(in) :: interior_rhs(:), trace_rhs(:)
        real(dp), intent(in) :: condensed_matrix_bar(:, :), condensed_rhs_bar(:)
        real(dp), intent(out) :: interior_matrix_bar(:, :), interior_to_trace_bar(:, :)
        real(dp), intent(out) :: trace_to_interior_bar(:, :), trace_matrix_bar(:, :)
        real(dp), intent(out) :: interior_rhs_bar(:), trace_rhs_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: interior_count, trace_count, info
        real(dp), allocatable :: interior_to_trace_solve(:, :), interior_rhs_solve(:)
        real(dp), allocatable :: x_bar(:, :), y_bar(:)
        real(dp), allocatable :: x_transpose_bar(:, :)
        real(dp), allocatable :: y_transpose_bar(:)

        interior_matrix_bar = 0.0_dp
        interior_to_trace_bar = 0.0_dp
        trace_to_interior_bar = 0.0_dp
        trace_matrix_bar = 0.0_dp
        interior_rhs_bar = 0.0_dp
        trace_rhs_bar = 0.0_dp
        call validate_hdg_inputs( &
            interior_matrix, interior_to_trace, trace_to_interior, trace_matrix, &
            interior_rhs, trace_rhs, condensed_matrix_bar, condensed_rhs_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        interior_count = size(interior_matrix, 1)
        trace_count = size(trace_matrix, 1)
        if (size(interior_matrix_bar, 1) /= interior_count .or. &
            size(interior_matrix_bar, 2) /= interior_count .or. &
            size(interior_to_trace_bar, 1) /= interior_count .or. &
            size(interior_to_trace_bar, 2) /= trace_count .or. &
            size(trace_to_interior_bar, 1) /= trace_count .or. &
            size(trace_to_interior_bar, 2) /= interior_count .or. &
            size(trace_matrix_bar, 1) /= trace_count .or. &
            size(trace_matrix_bar, 2) /= trace_count .or. &
            size(interior_rhs_bar) /= interior_count .or. &
            size(trace_rhs_bar) /= trace_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "HDG static condensation VJP has incompatible cotangents")
            return
        end if
        allocate(interior_to_trace_solve(interior_count, trace_count), &
            interior_rhs_solve(interior_count), x_bar(interior_count, trace_count), &
            y_bar(interior_count), x_transpose_bar(interior_count, trace_count), &
            y_transpose_bar(interior_count))
        call dense_solve(interior_matrix, interior_to_trace, &
            interior_to_trace_solve, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "HDG static condensation interior block is singular")
            return
        end if
        call dense_solve(interior_matrix, interior_rhs, interior_rhs_solve, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "HDG static condensation interior block is singular")
            return
        end if
        trace_matrix_bar = condensed_matrix_bar
        trace_rhs_bar = condensed_rhs_bar
        trace_to_interior_bar = -matmul(condensed_matrix_bar, &
            transpose(interior_to_trace_solve)) - outer_product(condensed_rhs_bar, interior_rhs_solve)
        y_bar = -matmul(transpose(trace_to_interior), condensed_rhs_bar)
        x_bar = -matmul(transpose(trace_to_interior), condensed_matrix_bar)
        call dense_solve(transpose(interior_matrix), x_bar, x_transpose_bar, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "HDG static condensation transposed interior block is singular")
            return
        end if
        call dense_solve(transpose(interior_matrix), y_bar, y_transpose_bar, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "HDG static condensation transposed interior block is singular")
            return
        end if
        interior_to_trace_bar = x_transpose_bar
        interior_rhs_bar = y_transpose_bar
        interior_matrix_bar = -matmul(x_transpose_bar, &
            transpose(interior_to_trace_solve)) - &
            outer_product(y_transpose_bar, interior_rhs_solve)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_hdg_static_condensation_vjp

    function outer_product(left, right) result(product)
        real(dp), intent(in) :: left(:), right(:)
        real(dp) :: product(size(left), size(right))

        product = spread(left, dim=2, ncopies=size(right))* &
            spread(right, dim=1, ncopies=size(left))
    end function outer_product

    logical function valid_hdg_direction( &
            interior_matrix_dot, interior_to_trace_dot, trace_to_interior_dot, &
            trace_matrix_dot, interior_rhs_dot, trace_rhs_dot, interior_count, trace_count)
        real(dp), intent(in) :: interior_matrix_dot(:, :), interior_to_trace_dot(:, :)
        real(dp), intent(in) :: trace_to_interior_dot(:, :), trace_matrix_dot(:, :)
        real(dp), intent(in) :: interior_rhs_dot(:), trace_rhs_dot(:)
        integer, intent(in) :: interior_count, trace_count

        valid_hdg_direction = size(interior_matrix_dot, 1) == interior_count .and. &
            size(interior_matrix_dot, 2) == interior_count .and. &
            size(interior_to_trace_dot, 1) == interior_count .and. &
            size(interior_to_trace_dot, 2) == trace_count .and. &
            size(trace_to_interior_dot, 1) == trace_count .and. &
            size(trace_to_interior_dot, 2) == interior_count .and. &
            size(trace_matrix_dot, 1) == trace_count .and. &
            size(trace_matrix_dot, 2) == trace_count .and. &
            size(interior_rhs_dot) == interior_count .and. &
            size(trace_rhs_dot) == trace_count .and. &
            all(ieee_is_finite(interior_matrix_dot)) .and. &
            all(ieee_is_finite(interior_to_trace_dot)) .and. &
            all(ieee_is_finite(trace_to_interior_dot)) .and. &
            all(ieee_is_finite(trace_matrix_dot)) .and. &
            all(ieee_is_finite(interior_rhs_dot)) .and. &
            all(ieee_is_finite(trace_rhs_dot))
    end function valid_hdg_direction

    subroutine validate_hdg_inputs( &
            interior_matrix, interior_to_trace, trace_to_interior, trace_matrix, &
            interior_rhs, trace_rhs, condensed_matrix, condensed_rhs, status)
        real(dp), intent(in) :: interior_matrix(:, :), interior_to_trace(:, :)
        real(dp), intent(in) :: trace_to_interior(:, :), trace_matrix(:, :)
        real(dp), intent(in) :: interior_rhs(:), trace_rhs(:)
        real(dp), intent(in) :: condensed_matrix(:, :), condensed_rhs(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: interior_count, trace_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "HDG static condensation received incompatible blocks")
        interior_count = size(interior_matrix, 1)
        trace_count = size(trace_matrix, 1)
        if (interior_count < 1 .or. trace_count < 1) return
        if (size(interior_matrix, 2) /= interior_count .or. &
            size(interior_to_trace, 1) /= interior_count .or. &
            size(interior_to_trace, 2) /= trace_count .or. &
            size(trace_to_interior, 1) /= trace_count .or. &
            size(trace_to_interior, 2) /= interior_count .or. &
            size(trace_matrix, 2) /= trace_count .or. &
            size(interior_rhs) /= interior_count .or. size(trace_rhs) /= trace_count .or. &
            size(condensed_matrix, 1) /= trace_count .or. &
            size(condensed_matrix, 2) /= trace_count .or. size(condensed_rhs) /= trace_count) return
        if (any(.not. ieee_is_finite(interior_matrix)) .or. &
            any(.not. ieee_is_finite(interior_to_trace)) .or. &
            any(.not. ieee_is_finite(trace_to_interior)) .or. &
            any(.not. ieee_is_finite(trace_matrix)) .or. &
            any(.not. ieee_is_finite(interior_rhs)) .or. &
            any(.not. ieee_is_finite(trace_rhs)) .or. &
            any(.not. ieee_is_finite(condensed_matrix)) .or. &
            any(.not. ieee_is_finite(condensed_rhs))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_hdg_inputs

end module fortfem_hdg_static_condensation
