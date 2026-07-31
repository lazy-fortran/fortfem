module fortfem_fci_plane_multigrid
    !! Two-level plane solver building block for PARALLAX-style FCI splits.
    !!
    !! The fine operator is a positive plane action.  One Jacobi smooth is
    !! followed by restriction, a direct coarse solve, prolongation, and one
    !! post-smooth.  The coarse solve is deliberately an implementation detail
    !! of this small contract; callers can replace it with a retained factor or
    !! a production multigrid hierarchy without changing the FCI line action.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_solve_csc
    use fortsparse, only: csc_is_valid, csc_matvec, csc_t, &
        fortsparse_status_t, status_set, FORTSPARSE_INVALID_MATRIX, &
        FORTSPARSE_OK
    implicit none

    private

    public :: apply_fci_plane_two_level_vcycle
    public :: apply_fci_plane_two_level_vcycles

contains

    subroutine apply_fci_plane_two_level_vcycle( &
            fine_operator, coarse_operator, restriction, prolongation, &
            diagonal, residual, correction, status)
        !! Apply one two-level V(1,1) cycle to a positive plane operator.
        !!
        !! `restriction` and `prolongation` are rectangular CSC maps with
        !! dimensions `(n_coarse,n_fine)` and `(n_fine,n_coarse)`.  The supplied
        !! diagonal is the positive fine-level smoother diagonal.  This routine
        !! solves the coarse residual directly; repeated solves should retain a
        !! factor or use a stronger multilevel implementation upstream.
        type(csc_t), intent(in) :: fine_operator
        type(csc_t), intent(in) :: coarse_operator
        type(csc_t), intent(in) :: restriction
        type(csc_t), intent(in) :: prolongation
        real(dp), intent(in) :: diagonal(:)
        real(dp), intent(in) :: residual(:)
        real(dp), intent(out) :: correction(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: fine_size, coarse_size, solve_status
        real(dp), allocatable :: fine_residual(:), coarse_rhs(:)
        real(dp), allocatable :: coarse_correction(:)

        correction = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI plane V-cycle received incompatible operators")
        if (.not. csc_is_valid(fine_operator) .or. &
            .not. csc_is_valid(coarse_operator) .or. &
            .not. csc_is_valid(restriction) .or. &
            .not. csc_is_valid(prolongation)) return
        if (fine_operator%nrow /= fine_operator%ncol .or. &
            coarse_operator%nrow /= coarse_operator%ncol) return
        fine_size = fine_operator%nrow
        coarse_size = coarse_operator%nrow
        if (restriction%nrow /= coarse_size .or. &
            restriction%ncol /= fine_size .or. &
            prolongation%nrow /= fine_size .or. &
            prolongation%ncol /= coarse_size) return
        if (size(diagonal) /= fine_size .or. size(residual) /= fine_size .or. &
            size(correction) /= fine_size) return
        if (any(.not. ieee_is_finite(fine_operator%val)) .or. &
            any(.not. ieee_is_finite(coarse_operator%val)) .or. &
            any(.not. ieee_is_finite(restriction%val)) .or. &
            any(.not. ieee_is_finite(prolongation%val)) .or. &
            any(.not. ieee_is_finite(diagonal)) .or. &
            any(.not. ieee_is_finite(residual)) .or. &
            any(diagonal <= 0.0_dp)) return

        allocate(fine_residual(fine_size), coarse_rhs(coarse_size))
        allocate(coarse_correction(coarse_size))

        correction = residual/diagonal
        fine_residual = residual - csc_matvec(fine_operator, correction)
        coarse_rhs = csc_matvec(restriction, fine_residual)
        call sparse_direct_solve_csc( &
            coarse_size, coarse_operator%col_ptr, coarse_operator%row_idx, &
            coarse_operator%val, coarse_rhs, coarse_correction, solve_status)
        if (solve_status /= 0) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI plane V-cycle coarse solve failed")
            return
        end if
        correction = correction + csc_matvec(prolongation, coarse_correction)
        fine_residual = residual - csc_matvec(fine_operator, correction)
        correction = correction + fine_residual/diagonal
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_fci_plane_two_level_vcycle

    subroutine apply_fci_plane_two_level_vcycles( &
            fine_operators, coarse_operators, restrictions, prolongations, &
            diagonals, residual, correction, status)
        !! Apply independent two-level cycles to a homogeneous plane stack.
        !!
        !! The plane unknowns are stored contiguously, with the same fine
        !! plane size on every block.  This is the small PARALLAX-style
        !! field-split contract: each perpendicular plane solve is independent
        !! and the caller owns the coupling to the parallel line action.
        type(csc_t), intent(in) :: fine_operators(:)
        type(csc_t), intent(in) :: coarse_operators(:)
        type(csc_t), intent(in) :: restrictions(:)
        type(csc_t), intent(in) :: prolongations(:)
        real(dp), intent(in) :: diagonals(:, :)
        real(dp), intent(in) :: residual(:)
        real(dp), intent(out) :: correction(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: plane_count, fine_size, plane, offset
        real(dp), allocatable :: candidate(:), plane_correction(:)
        type(fortsparse_status_t) :: plane_status

        correction = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI batched plane V-cycle received incompatible operators")
        plane_count = size(fine_operators)
        if (plane_count < 1) return
        if (size(coarse_operators) /= plane_count .or. &
            size(restrictions) /= plane_count .or. &
            size(prolongations) /= plane_count) return
        if (.not. csc_is_valid(fine_operators(1))) return
        if (fine_operators(1)%nrow /= fine_operators(1)%ncol) return
        fine_size = fine_operators(1)%nrow
        if (fine_size < 1) return
        if (size(diagonals, 1) /= fine_size .or. &
            size(diagonals, 2) /= plane_count) return
        if (size(residual) /= fine_size*plane_count .or. &
            size(correction) /= fine_size*plane_count) return
        do plane = 1, plane_count
            if (.not. csc_is_valid(fine_operators(plane))) return
            if (fine_operators(plane)%nrow /= fine_size .or. &
                fine_operators(plane)%ncol /= fine_size) return
        end do

        allocate(candidate(size(correction)), plane_correction(fine_size))
        candidate = 0.0_dp
        do plane = 1, plane_count
            offset = (plane - 1)*fine_size
            call apply_fci_plane_two_level_vcycle( &
                fine_operators(plane), coarse_operators(plane), &
                restrictions(plane), prolongations(plane), diagonals(:, plane), &
                residual(offset + 1:offset + fine_size), plane_correction, &
                plane_status)
            if (plane_status%code /= FORTSPARSE_OK) then
                call status_set(status, plane_status%code, plane_status%msg)
                return
            end if
            candidate(offset + 1:offset + fine_size) = plane_correction
        end do
        correction = candidate
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_fci_plane_two_level_vcycles

end module fortfem_fci_plane_multigrid
