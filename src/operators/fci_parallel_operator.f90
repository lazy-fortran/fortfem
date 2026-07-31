module fortfem_fci_parallel_operator
    !! Sparse support-operator building blocks for field-coordinate-independent
    !! (FCI) parallel derivatives.
    !!
    !! The map arrays describe interpolation from the lower and upper
    !! poloidal planes of each toroidal segment to its staggered flux boxes.
    !! The gradient is assembled as
    !!
    !!   Q = (Q_plus - Q_minus) / ell,
    !!
    !! while the support divergence is the negative volume-weighted adjoint
    !!
    !!   P = -W_c^{-1} Q^T W_s.
    !!
    !! This is the conservative pairing used by PARALLAX-style FCI support
    !! operators.  Field-line tracing and interpolation coefficient
    !! construction remain separate geometry services; this module only owns
    !! the algebraic contract and its sparse representation.
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_is_valid, csc_t, &
        fortsparse_status_t, status_set, FORTSPARSE_INVALID_MATRIX, &
        FORTSPARSE_OK
    implicit none
    private

    public :: assemble_fci_parallel_gradient_csc
    public :: assemble_fci_parallel_support_divergence_csc

contains

    subroutine assemble_fci_parallel_gradient_csc( &
            forward_map, backward_map, line_lengths, gradient, status)
        !! Assemble Q from mapped upper/lower plane interpolation matrices.
        !!
        !! `forward_map(:, :, segment)` maps the upper plane of a segment and
        !! `backward_map(:, :, segment)` maps its lower plane.  The canonical
        !! unknown vector is ordered plane-by-plane, and the staggered result
        !! is ordered segment-by-segment.  Each map may be sparse in practice;
        !! zero entries are omitted from the resulting CSC matrix.
        real(dp), intent(in) :: forward_map(:, :, :)
        real(dp), intent(in) :: backward_map(:, :, :)
        real(dp), intent(in) :: line_lengths(:, :)
        type(csc_t), intent(out) :: gradient
        type(fortsparse_status_t), intent(out) :: status

        integer :: n_plane, n_segment, n_staggered
        integer :: segment, sample, plane_node, entry_count
        integer :: lower_column, upper_column, row_count, column_count
        integer, allocatable :: rows(:), columns(:)
        real(dp), allocatable :: values(:)
        real(dp) :: coefficient

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI gradient assembly received incompatible maps")
        n_staggered = size(forward_map, 1)
        n_plane = size(forward_map, 2)
        n_segment = size(forward_map, 3)
        if (n_staggered < 1 .or. n_plane < 1 .or. n_segment < 1) return
        if (any(shape(backward_map) /= shape(forward_map))) return
        if (size(line_lengths, 1) /= n_staggered) return
        if (size(line_lengths, 2) /= n_segment) return
        if (any(line_lengths <= 0.0_dp)) return

        row_count = n_staggered*n_segment
        column_count = n_plane*(n_segment + 1)
        allocate( &
            rows(2*n_staggered*n_plane*n_segment), &
            columns(2*n_staggered*n_plane*n_segment), &
            values(2*n_staggered*n_plane*n_segment))
        entry_count = 0

        do segment = 1, n_segment
            do sample = 1, n_staggered
                lower_column = (segment - 1)*n_plane
                upper_column = segment*n_plane
                do plane_node = 1, n_plane
                    coefficient = -backward_map(sample, plane_node, segment) / &
                        line_lengths(sample, segment)
                    if (coefficient /= 0.0_dp) then
                        entry_count = entry_count + 1
                        rows(entry_count) = &
                            sample + (segment - 1)*n_staggered
                        columns(entry_count) = lower_column + plane_node
                        values(entry_count) = coefficient
                    end if
                    coefficient = forward_map(sample, plane_node, segment) / &
                        line_lengths(sample, segment)
                    if (coefficient /= 0.0_dp) then
                        entry_count = entry_count + 1
                        rows(entry_count) = &
                            sample + (segment - 1)*n_staggered
                        columns(entry_count) = upper_column + plane_node
                        values(entry_count) = coefficient
                    end if
                end do
            end do
        end do

        call csc_from_triplet( &
            row_count, column_count, rows(:entry_count), columns(:entry_count), &
            values(:entry_count), gradient, status)
        if (status%code == FORTSPARSE_OK) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI gradient CSC construction failed")
    end subroutine assemble_fci_parallel_gradient_csc

    subroutine assemble_fci_parallel_support_divergence_csc( &
            gradient, canonical_volumes, staggered_volumes, divergence, status)
        !! Assemble P = -W_c^{-1} Q^T W_s from a FCI gradient Q.
        !!
        !! The returned operator satisfies, up to floating-point roundoff,
        !!
        !!   u^T W_c P f = -(Q u)^T W_s f,
        !!
        !! which is the support-operator conservation/adjointness contract.
        type(csc_t), intent(in) :: gradient
        real(dp), intent(in) :: canonical_volumes(:)
        real(dp), intent(in) :: staggered_volumes(:)
        type(csc_t), intent(out) :: divergence
        type(fortsparse_status_t), intent(out) :: status

        integer :: column, entry, row
        integer, allocatable :: rows(:), columns(:)
        real(dp), allocatable :: values(:)

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI support divergence received an invalid gradient")
        if (.not. csc_is_valid(gradient)) return
        if (size(canonical_volumes) /= gradient%ncol) return
        if (size(staggered_volumes) /= gradient%nrow) return
        if (any(canonical_volumes <= 0.0_dp)) return
        if (any(staggered_volumes <= 0.0_dp)) return

        allocate(rows(gradient%nnz), columns(gradient%nnz))
        allocate(values(gradient%nnz))
        entry = 0
        do column = 1, gradient%ncol
            do row = gradient%col_ptr(column), gradient%col_ptr(column + 1) - 1
                entry = entry + 1
                rows(entry) = column
                columns(entry) = gradient%row_idx(row)
                values(entry) = -staggered_volumes(gradient%row_idx(row)) / &
                    canonical_volumes(column)*gradient%val(row)
            end do
        end do

        call csc_from_triplet( &
            gradient%ncol, gradient%nrow, rows, columns, values, divergence, &
            status)
        if (status%code == FORTSPARSE_OK) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI support divergence CSC construction failed")
    end subroutine assemble_fci_parallel_support_divergence_csc

end module fortfem_fci_parallel_operator
