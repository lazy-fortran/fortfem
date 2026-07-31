module fortfem_fci_boundary_patch_mortar
    !! Conservative transfer between FCI and fitted boundary-patch traces.
    !!
    !! The cross-mass is assembled on an overlap quadrature.  Its row and
    !! column sums are the diagonal trace measures, so the two normalized
    !! transfer operators reproduce constants and are weighted adjoints.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_mortar_trace_coupling, only: assemble_mortar_trace_coupling
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_fci_boundary_patch_mortar
    public :: assemble_fci_boundary_patch_mortar_jvp
    public :: assemble_fci_boundary_patch_mortar_vjp

contains

    subroutine assemble_fci_boundary_patch_mortar( &
            background_trace, patch_trace, overlap_weights, ownership_count, &
            cross_mass, background_mass, patch_mass, background_from_patch, &
            patch_from_background, overlap_measure, status)
        !! Assemble cross-mass and conservative transfers on one overlap.
        !!
        !! `background_from_patch` maps patch coefficients to background
        !! coefficients.  `patch_from_background` is its weighted adjoint.
        !! Every overlap quadrature row must have exactly one owner.
        real(dp), intent(in) :: background_trace(:, :), patch_trace(:, :)
        real(dp), intent(in) :: overlap_weights(:)
        integer, intent(in) :: ownership_count(:)
        real(dp), intent(out) :: cross_mass(:, :), background_mass(:), patch_mass(:)
        real(dp), intent(out) :: background_from_patch(:, :)
        real(dp), intent(out) :: patch_from_background(:, :)
        real(dp), intent(out) :: overlap_measure
        type(fortsparse_status_t), intent(out) :: status

        integer :: background, quadrature_count, background_count, patch_count
        integer :: patch
        real(dp), allocatable :: background_partition(:), patch_partition(:)
        real(dp) :: scale

        cross_mass = 0.0_dp
        background_mass = 0.0_dp
        patch_mass = 0.0_dp
        background_from_patch = 0.0_dp
        patch_from_background = 0.0_dp
        overlap_measure = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI boundary-patch mortar received incompatible arrays")

        quadrature_count = size(background_trace, 1)
        background_count = size(background_trace, 2)
        patch_count = size(patch_trace, 2)
        if (quadrature_count < 1) return
        if (background_count < 1 .or. patch_count < 1) return
        if (size(patch_trace, 1) /= quadrature_count) return
        if (size(overlap_weights) /= quadrature_count) return
        if (size(ownership_count) /= quadrature_count) return
        if (size(cross_mass, 1) /= background_count) return
        if (size(cross_mass, 2) /= patch_count) return
        if (size(background_mass) /= background_count) return
        if (size(patch_mass) /= patch_count) return
        if (size(background_from_patch, 1) /= background_count) return
        if (size(background_from_patch, 2) /= patch_count) return
        if (size(patch_from_background, 1) /= patch_count) return
        if (size(patch_from_background, 2) /= background_count) return
        if (any(ownership_count /= 1)) return
        if (any(.not. ieee_is_finite(background_trace))) return
        if (any(.not. ieee_is_finite(patch_trace))) return
        if (any(.not. ieee_is_finite(overlap_weights))) return
        if (any(overlap_weights <= 0.0_dp)) return

        allocate(background_partition(quadrature_count))
        allocate(patch_partition(quadrature_count))
        background_partition = sum(background_trace, dim=2)
        patch_partition = sum(patch_trace, dim=2)
        scale = 100.0_dp*epsilon(1.0_dp)
        if (maxval(abs(background_partition - 1.0_dp)) > scale) return
        if (maxval(abs(patch_partition - 1.0_dp)) > scale) return

        call assemble_mortar_trace_coupling( &
            background_trace, patch_trace, overlap_weights, cross_mass, status)
        if (status%code /= FORTSPARSE_OK) return
        if (numerical_rank(cross_mass) < min(background_count, patch_count)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI boundary-patch mortar cross-mass is rank deficient")
            return
        end if

        background_mass = sum(cross_mass, dim=2)
        patch_mass = sum(cross_mass, dim=1)
        if (any(background_mass <= 0.0_dp) .or. any(patch_mass <= 0.0_dp)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI boundary-patch mortar has an uncovered trace degree of freedom")
            return
        end if
        overlap_measure = sum(overlap_weights)
        do background = 1, background_count
            do patch = 1, patch_count
                background_from_patch(background, patch) = &
                    cross_mass(background, patch)/background_mass(background)
                patch_from_background(patch, background) = &
                    cross_mass(background, patch)/patch_mass(patch)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_fci_boundary_patch_mortar

    subroutine assemble_fci_boundary_patch_mortar_jvp( &
            background_trace, patch_trace, overlap_weights, ownership_count, &
            background_trace_dot, patch_trace_dot, overlap_weights_dot, &
            cross_mass_dot, background_mass_dot, patch_mass_dot, &
            background_from_patch_dot, patch_from_background_dot, &
            overlap_measure_dot, status)
        !! Apply the fixed-topology JVP of the mortar transfer contract.
        real(dp), intent(in) :: background_trace(:, :), patch_trace(:, :)
        real(dp), intent(in) :: overlap_weights(:)
        integer, intent(in) :: ownership_count(:)
        real(dp), intent(in) :: background_trace_dot(:, :), patch_trace_dot(:, :)
        real(dp), intent(in) :: overlap_weights_dot(:)
        real(dp), intent(out) :: cross_mass_dot(:, :)
        real(dp), intent(out) :: background_mass_dot(:), patch_mass_dot(:)
        real(dp), intent(out) :: background_from_patch_dot(:, :)
        real(dp), intent(out) :: patch_from_background_dot(:, :)
        real(dp), intent(out) :: overlap_measure_dot
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: cross_mass(:, :), background_mass(:)
        real(dp), allocatable :: patch_mass(:), background_from_patch(:, :)
        real(dp), allocatable :: patch_from_background(:, :)
        integer :: background, quadrature, patch
        integer :: quadrature_count, background_count, patch_count

        cross_mass_dot = 0.0_dp
        background_mass_dot = 0.0_dp
        patch_mass_dot = 0.0_dp
        background_from_patch_dot = 0.0_dp
        patch_from_background_dot = 0.0_dp
        overlap_measure_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI boundary-patch mortar JVP received incompatible arrays")

        quadrature_count = size(background_trace, 1)
        background_count = size(background_trace, 2)
        patch_count = size(patch_trace, 2)
        if (quadrature_count < 1) return
        if (size(background_trace_dot, 1) /= quadrature_count) return
        if (size(background_trace_dot, 2) /= background_count) return
        if (size(patch_trace_dot, 1) /= quadrature_count) return
        if (size(patch_trace_dot, 2) /= patch_count) return
        if (size(overlap_weights_dot) /= quadrature_count) return
        if (size(cross_mass_dot, 1) /= background_count) return
        if (size(cross_mass_dot, 2) /= patch_count) return
        if (size(background_mass_dot) /= background_count) return
        if (size(patch_mass_dot) /= patch_count) return
        if (size(background_from_patch_dot, 1) /= background_count) return
        if (size(background_from_patch_dot, 2) /= patch_count) return
        if (size(patch_from_background_dot, 1) /= patch_count) return
        if (size(patch_from_background_dot, 2) /= background_count) return
        if (any(.not. ieee_is_finite(background_trace_dot))) return
        if (any(.not. ieee_is_finite(patch_trace_dot))) return
        if (any(.not. ieee_is_finite(overlap_weights_dot))) return

        allocate(cross_mass(background_count, patch_count))
        allocate(background_mass(background_count), patch_mass(patch_count))
        allocate(background_from_patch(background_count, patch_count))
        allocate(patch_from_background(patch_count, background_count))
        call assemble_fci_boundary_patch_mortar( &
            background_trace, patch_trace, overlap_weights, ownership_count, &
            cross_mass, background_mass, patch_mass, background_from_patch, &
            patch_from_background, overlap_measure_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        overlap_measure_dot = sum(overlap_weights_dot)
        do quadrature = 1, quadrature_count
            do background = 1, background_count
                do patch = 1, patch_count
                    cross_mass_dot(background, patch) = &
                        cross_mass_dot(background, patch) + &
                        overlap_weights_dot(quadrature)* &
                        background_trace(quadrature, background)* &
                        patch_trace(quadrature, patch) + &
                        overlap_weights(quadrature)* &
                        background_trace_dot(quadrature, background)* &
                        patch_trace(quadrature, patch) + &
                        overlap_weights(quadrature)* &
                        background_trace(quadrature, background)* &
                        patch_trace_dot(quadrature, patch)
                end do
            end do
        end do
        background_mass_dot = sum(cross_mass_dot, dim=2)
        patch_mass_dot = sum(cross_mass_dot, dim=1)
        do background = 1, background_count
            do patch = 1, patch_count
                background_from_patch_dot(background, patch) = &
                    cross_mass_dot(background, patch)/background_mass(background) - &
                    cross_mass(background, patch)*background_mass_dot(background) &
                    /background_mass(background)**2
                patch_from_background_dot(patch, background) = &
                    cross_mass_dot(background, patch)/patch_mass(patch) - &
                    cross_mass(background, patch)*patch_mass_dot(patch) &
                    /patch_mass(patch)**2
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_fci_boundary_patch_mortar_jvp

    subroutine assemble_fci_boundary_patch_mortar_vjp( &
            background_trace, patch_trace, overlap_weights, ownership_count, &
            cross_mass_bar, background_mass_bar, patch_mass_bar, &
            background_from_patch_bar, patch_from_background_bar, &
            overlap_measure_bar, background_trace_bar, patch_trace_bar, &
            overlap_weights_bar, status)
        !! Apply the real fixed-topology VJP of the mortar transfer contract.
        real(dp), intent(in) :: background_trace(:, :), patch_trace(:, :)
        real(dp), intent(in) :: overlap_weights(:)
        integer, intent(in) :: ownership_count(:)
        real(dp), intent(in) :: cross_mass_bar(:, :)
        real(dp), intent(in) :: background_mass_bar(:), patch_mass_bar(:)
        real(dp), intent(in) :: background_from_patch_bar(:, :)
        real(dp), intent(in) :: patch_from_background_bar(:, :)
        real(dp), intent(in) :: overlap_measure_bar
        real(dp), intent(out) :: background_trace_bar(:, :), patch_trace_bar(:, :)
        real(dp), intent(out) :: overlap_weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: cross_mass(:, :), background_mass(:)
        real(dp), allocatable :: patch_mass(:), background_from_patch(:, :)
        real(dp), allocatable :: patch_from_background(:, :), cross_bar(:, :)
        real(dp), allocatable :: background_mass_total_bar(:), patch_mass_total_bar(:)
        real(dp) :: overlap_measure, overlap_measure_cotangent
        integer :: background, quadrature, patch
        integer :: quadrature_count, background_count, patch_count

        background_trace_bar = 0.0_dp
        patch_trace_bar = 0.0_dp
        overlap_weights_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI boundary-patch mortar VJP received incompatible arrays")

        quadrature_count = size(background_trace, 1)
        background_count = size(background_trace, 2)
        patch_count = size(patch_trace, 2)
        if (quadrature_count < 1) return
        if (size(patch_trace, 1) /= quadrature_count) return
        if (size(overlap_weights) /= quadrature_count) return
        if (size(cross_mass_bar, 1) /= background_count) return
        if (size(cross_mass_bar, 2) /= patch_count) return
        if (size(background_mass_bar) /= background_count) return
        if (size(patch_mass_bar) /= patch_count) return
        if (size(background_from_patch_bar, 1) /= background_count) return
        if (size(background_from_patch_bar, 2) /= patch_count) return
        if (size(patch_from_background_bar, 1) /= patch_count) return
        if (size(patch_from_background_bar, 2) /= background_count) return
        if (size(background_trace_bar, 1) /= quadrature_count) return
        if (size(background_trace_bar, 2) /= background_count) return
        if (size(patch_trace_bar, 1) /= quadrature_count) return
        if (size(patch_trace_bar, 2) /= patch_count) return
        if (size(overlap_weights_bar) /= quadrature_count) return
        if (.not. ieee_is_finite(overlap_measure_bar)) return
        if (any(.not. ieee_is_finite(cross_mass_bar))) return
        if (any(.not. ieee_is_finite(background_mass_bar))) return
        if (any(.not. ieee_is_finite(patch_mass_bar))) return
        if (any(.not. ieee_is_finite(background_from_patch_bar))) return
        if (any(.not. ieee_is_finite(patch_from_background_bar))) return
        overlap_measure_cotangent = overlap_measure_bar

        allocate(cross_mass(background_count, patch_count))
        allocate(background_mass(background_count), patch_mass(patch_count))
        allocate(background_from_patch(background_count, patch_count))
        allocate(patch_from_background(patch_count, background_count))
        call assemble_fci_boundary_patch_mortar( &
            background_trace, patch_trace, overlap_weights, ownership_count, &
            cross_mass, background_mass, patch_mass, background_from_patch, &
            patch_from_background, overlap_measure, status)
        if (status%code /= FORTSPARSE_OK) return

        allocate(cross_bar(background_count, patch_count))
        allocate(background_mass_total_bar(background_count))
        allocate(patch_mass_total_bar(patch_count))
        cross_bar = cross_mass_bar
        background_mass_total_bar = background_mass_bar
        patch_mass_total_bar = patch_mass_bar
        do background = 1, background_count
            do patch = 1, patch_count
                cross_bar(background, patch) = cross_bar(background, patch) + &
                    background_from_patch_bar(background, patch) &
                    /background_mass(background) + &
                    patch_from_background_bar(patch, background) &
                    /patch_mass(patch)
                background_mass_total_bar(background) = &
                    background_mass_total_bar(background) - &
                    background_from_patch_bar(background, patch)* &
                    cross_mass(background, patch)/background_mass(background)**2
                patch_mass_total_bar(patch) = patch_mass_total_bar(patch) - &
                    patch_from_background_bar(patch, background)* &
                    cross_mass(background, patch)/patch_mass(patch)**2
            end do
        end do
        do background = 1, background_count
            do patch = 1, patch_count
                cross_bar(background, patch) = cross_bar(background, patch) + &
                    background_mass_total_bar(background) + &
                    patch_mass_total_bar(patch)
            end do
        end do
        overlap_weights_bar = overlap_measure_cotangent
        do quadrature = 1, quadrature_count
            do background = 1, background_count
                do patch = 1, patch_count
                    background_trace_bar(quadrature, background) = &
                        background_trace_bar(quadrature, background) + &
                        overlap_weights(quadrature)*patch_trace(quadrature, patch)* &
                        cross_bar(background, patch)
                    patch_trace_bar(quadrature, patch) = &
                        patch_trace_bar(quadrature, patch) + &
                        overlap_weights(quadrature)* &
                        background_trace(quadrature, background)* &
                        cross_bar(background, patch)
                    overlap_weights_bar(quadrature) = &
                        overlap_weights_bar(quadrature) + &
                        background_trace(quadrature, background)* &
                        patch_trace(quadrature, patch)*cross_bar(background, patch)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_fci_boundary_patch_mortar_vjp

    integer function numerical_rank(matrix) result(rank)
        !! Compute a scale-aware rank for the small dense cross-mass block.
        real(dp), intent(in) :: matrix(:, :)

        real(dp), allocatable :: work(:, :)
        real(dp) :: factor, matrix_scale, pivot, tolerance
        integer :: column, row, pivot_row, row_count, column_count

        row_count = size(matrix, 1)
        column_count = size(matrix, 2)
        rank = 0
        if (row_count < 1 .or. column_count < 1) return
        allocate(work(row_count, column_count))
        work = matrix
        matrix_scale = maxval(abs(matrix))
        if (matrix_scale <= 0.0_dp) return
        tolerance = 100.0_dp*epsilon(1.0_dp)*matrix_scale* &
            real(max(row_count, column_count), dp)
        do column = 1, column_count
            pivot_row = rank + 1
            if (pivot_row > row_count) exit
            pivot = 0.0_dp
            do row = rank + 1, row_count
                if (abs(work(row, column)) > pivot) then
                    pivot = abs(work(row, column))
                    pivot_row = row
                end if
            end do
            if (pivot <= tolerance) cycle
            if (pivot_row /= rank + 1) then
                call swap_rows(work, pivot_row, rank + 1)
            end if
            rank = rank + 1
            do row = rank + 1, row_count
                factor = work(row, column)/work(rank, column)
                work(row, column:column_count) = &
                    work(row, column:column_count) - &
                    factor*work(rank, column:column_count)
            end do
            if (rank == row_count) exit
        end do
    end function numerical_rank

    subroutine swap_rows(matrix, first, second)
        real(dp), intent(inout) :: matrix(:, :)
        integer, intent(in) :: first, second

        real(dp) :: temporary
        integer :: column

        do column = 1, size(matrix, 2)
            temporary = matrix(first, column)
            matrix(first, column) = matrix(second, column)
            matrix(second, column) = temporary
        end do
    end subroutine swap_rows

end module fortfem_fci_boundary_patch_mortar
