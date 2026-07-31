module fortfem_fci_field_split_preconditioner
    !! Additive PARALLAX-style parallel/perpendicular field split.
    !!
    !! The parallel block is represented by a cached positive FCI Jacobi
    !! diagonal.  The perpendicular block is an independently supplied ragged
    !! plane V-cycle.  This module only composes the two algebraic actions; it
    !! does not assume a particular physical coefficient or boundary model.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_fci_plane_multigrid, only: &
        apply_fci_plane_two_level_vcycles_ragged
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_t, fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none

    private

    public :: apply_fci_additive_field_split_preconditioner

contains

    subroutine apply_fci_additive_field_split_preconditioner( &
            parallel_diagonal, plane_diagonal, fine_operators, coarse_operators, &
            restrictions, prolongations, plane_offsets, parallel_weight, &
            plane_weight, residual, correction, status)
        !! Apply a weighted additive parallel/perpendicular preconditioner.
        !!
        !! For positive weights this evaluates
        !!
        !!   M^{-1}r = w_parallel D_parallel^{-1}r
        !!             + w_plane M_plane^{-1}r,
        !!
        !! where `M_plane^{-1}` is the ragged two-level plane V-cycle.  The
        !! diagonal arrays are deliberately caller-owned so repeated Krylov
        !! iterations can cache geometry-dependent quantities.
        real(dp), intent(in) :: parallel_diagonal(:), plane_diagonal(:)
        type(csc_t), intent(in) :: fine_operators(:), coarse_operators(:)
        type(csc_t), intent(in) :: restrictions(:), prolongations(:)
        integer, intent(in) :: plane_offsets(:)
        real(dp), intent(in) :: parallel_weight, plane_weight
        real(dp), intent(in) :: residual(:)
        real(dp), intent(out) :: correction(:)
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: parallel_correction(:), plane_correction(:)
        type(fortsparse_status_t) :: plane_status

        correction = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI additive field split received incompatible arrays")
        if (size(residual) < 1 .or. size(correction) /= size(residual)) return
        if (size(parallel_diagonal) /= size(residual) .or. &
            size(plane_diagonal) /= size(residual)) return
        if (.not. ieee_is_finite(parallel_weight) .or. &
            .not. ieee_is_finite(plane_weight)) return
        if (parallel_weight < 0.0_dp .or. plane_weight < 0.0_dp .or. &
            parallel_weight + plane_weight <= 0.0_dp) return
        if (any(.not. ieee_is_finite(parallel_diagonal)) .or. &
            any(.not. ieee_is_finite(plane_diagonal)) .or. &
            any(.not. ieee_is_finite(residual))) return
        if (any(parallel_diagonal <= 0.0_dp) .or. &
            any(plane_diagonal <= 0.0_dp)) return

        allocate(parallel_correction(size(residual)), &
            plane_correction(size(residual)))
        parallel_correction = residual/parallel_diagonal
        call apply_fci_plane_two_level_vcycles_ragged( &
            fine_operators, coarse_operators, restrictions, prolongations, &
            plane_offsets, plane_diagonal, residual, plane_correction, plane_status)
        if (plane_status%code /= FORTSPARSE_OK) then
            call status_set(status, plane_status%code, plane_status%msg)
            return
        end if
        correction = parallel_weight*parallel_correction + &
            plane_weight*plane_correction
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_fci_additive_field_split_preconditioner

end module fortfem_fci_field_split_preconditioner
