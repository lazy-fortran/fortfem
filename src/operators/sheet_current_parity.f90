module fortfem_sheet_current_parity
    !! Geometry-neutral parity contract for an explicit sheet and its
    !! resolved Gaussian volume representation.
    !!
    !! The normal quadrature is a one-dimensional slab cross-section.  A
    !! caller supplies a constant tangential K and a physical surface measure
    !! A.  The explicit representation is A*K, while the resolved
    !! representation is the normal integral of K*delta_epsilon.  Keeping the
    !! two ledgers together gives fitted, cut, DG, and IGA clients one small
    !! manufactured oracle without selecting a PDE or constitutive law.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_regularized_surface_current_layer, only: &
        evaluate_regularized_surface_current_integral
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    real(dp), parameter :: equality_tolerance = 2.0e-12_dp

    public :: evaluate_sheet_current_parity

contains

    subroutine evaluate_sheet_current_parity( &
            signed_distance, normal_weights, sheet_current, thickness, &
            surface_measure, explicit_current, regularized_integrated, &
            explicit_integrated, relative_error, status)
        !! Compare a resolved slab layer with the explicit surface ledger.
        real(dp), intent(in) :: signed_distance(:), normal_weights(:)
        real(dp), intent(in) :: sheet_current(:, :), thickness, surface_measure
        real(dp), intent(in) :: explicit_current(:)
        real(dp), intent(out) :: regularized_integrated(:), explicit_integrated(:)
        real(dp), intent(out) :: relative_error
        type(fortsparse_status_t), intent(out) :: status

        real(dp) :: normalization
        integer :: quadrature, component
        real(dp) :: scale

        regularized_integrated = 0.0_dp
        explicit_integrated = 0.0_dp
        relative_error = 0.0_dp
        call validate_inputs( &
            signed_distance, normal_weights, sheet_current, thickness, &
            surface_measure, explicit_current, regularized_integrated, &
            explicit_integrated, status)
        if (status%code /= FORTSPARSE_OK) return

        do quadrature = 1, size(signed_distance)
            if (maxval(abs(sheet_current(quadrature, :) - explicit_current)) > &
                equality_tolerance*max(1.0_dp, maxval(abs(explicit_current)))) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "sheet-current parity requires a constant tangential K")
                return
            end if
        end do

        call evaluate_regularized_surface_current_integral( &
            signed_distance, normal_weights, sheet_current, thickness, &
            normalization, regularized_integrated, status)
        if (status%code /= FORTSPARSE_OK) return
        regularized_integrated = surface_measure*regularized_integrated
        explicit_integrated = surface_measure*explicit_current
        scale = max(1.0_dp, sqrt(sum(explicit_integrated**2)))
        relative_error = sqrt(sum((regularized_integrated - &
            explicit_integrated)**2))/scale
        if (.not. ieee_is_finite(relative_error)) then
            regularized_integrated = 0.0_dp
            explicit_integrated = 0.0_dp
            relative_error = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "sheet-current parity produced a non-finite error")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_sheet_current_parity

    subroutine validate_inputs( &
            signed_distance, normal_weights, sheet_current, thickness, &
            surface_measure, explicit_current, regularized_integrated, &
            explicit_integrated, status)
        real(dp), intent(in) :: signed_distance(:), normal_weights(:)
        real(dp), intent(in) :: sheet_current(:, :), thickness, surface_measure
        real(dp), intent(in) :: explicit_current(:)
        real(dp), intent(in) :: regularized_integrated(:), explicit_integrated(:)
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "sheet-current parity received incompatible arrays")
        if (size(signed_distance) < 2 .or. size(normal_weights) /= &
            size(signed_distance) .or. size(sheet_current, 1) /= &
            size(signed_distance) .or. size(sheet_current, 2) /= 3 .or. &
            size(explicit_current) /= 3 .or. size(regularized_integrated) /= 3 .or. &
            size(explicit_integrated) /= 3) return
        if (.not. ieee_is_finite(thickness) .or. thickness <= 0.0_dp .or. &
            .not. ieee_is_finite(surface_measure) .or. surface_measure <= 0.0_dp) return
        if (.not. all(ieee_is_finite(signed_distance)) .or. &
            .not. all(ieee_is_finite(normal_weights)) .or. &
            .not. all(ieee_is_finite(sheet_current)) .or. &
            .not. all(ieee_is_finite(explicit_current)) .or. &
            any(normal_weights <= 0.0_dp)) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_inputs

end module fortfem_sheet_current_parity
