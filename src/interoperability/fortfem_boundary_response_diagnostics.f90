module fortfem_boundary_response_diagnostics
    !! Work-weighted reciprocity and passivity diagnostics for boundary maps.
    !!
    !! A boundary operator maps a trace to a work-conjugate response.  Given
    !! positive diagonal quadrature/physical work weights W and a square
    !! complex response A, this neutral diagnostic inspects W A.  It does not
    !! build a FEM, BEM, DtN, PML, wall, or free-boundary operator and does not
    !! select an equation or normalization convention.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: evaluate_weighted_boundary_response_diagnostics

contains

    subroutine evaluate_weighted_boundary_response_diagnostics( &
            response_matrix, work_weights, reciprocity_error, &
            passivity_lower_bound, status)
        !! Return normalized weighted-transpose and Hermitian certificates.
        !!
        !! The reciprocity defect is
        !!
        !!   || W A - (W A)^T ||_max / max(1, ||W A||_max),
        !!
        !! while the passivity value is the Gershgorin lower bound of
        !!
        !!   (W A + (W A)^H)/2.
        complex(dp), intent(in) :: response_matrix(:, :)
        real(dp), intent(in) :: work_weights(:)
        real(dp), intent(out) :: reciprocity_error, passivity_lower_bound
        integer, intent(out) :: status

        complex(dp), allocatable :: weighted_response(:, :), hermitian_part(:, :)
        real(dp) :: scale, radius
        integer :: row, count

        reciprocity_error = huge(1.0_dp)
        passivity_lower_bound = -huge(1.0_dp)
        status = 1
        count = size(response_matrix, 1)
        if (count < 1 .or. size(response_matrix, 2) /= count .or. &
            size(work_weights) /= count) return
        if (.not. all(ieee_is_finite(work_weights)) .or. &
            any(work_weights <= 0.0_dp) .or. &
            .not. finite_complex(response_matrix)) return

        allocate(weighted_response(count, count), hermitian_part(count, count))
        weighted_response = spread(work_weights, 2, count)*response_matrix
        scale = max(1.0_dp, maxval(abs(weighted_response)))
        reciprocity_error = maxval(abs( &
            weighted_response - transpose(weighted_response)))/scale
        hermitian_part = 0.5_dp*(weighted_response + &
            conjg(transpose(weighted_response)))
        passivity_lower_bound = huge(1.0_dp)
        do row = 1, count
            radius = sum(abs(hermitian_part(row, :))) - &
                abs(hermitian_part(row, row))
            passivity_lower_bound = min(passivity_lower_bound, &
                real(hermitian_part(row, row), dp) - radius)
        end do
        status = 0
    end subroutine evaluate_weighted_boundary_response_diagnostics

    logical function finite_complex(values) result(valid)
        complex(dp), intent(in) :: values(..)

        valid = .true.
        select rank (values)
        rank (0)
            valid = ieee_is_finite(real(values, dp)) .and. &
                ieee_is_finite(aimag(values))
        rank (1)
            valid = all(ieee_is_finite(real(values, dp))) .and. &
                all(ieee_is_finite(aimag(values)))
        rank (2)
            valid = all(ieee_is_finite(real(values, dp))) .and. &
                all(ieee_is_finite(aimag(values)))
        rank default
            valid = .false.
        end select
    end function finite_complex

end module fortfem_boundary_response_diagnostics
