module fortfem_source_trace_map
    !! Neutral dense source-to-trace map and its complex products.
    !!
    !! The map is the caller-supplied linear action
    !!
    !!     trace = source_matrix source_coefficients.
    !!
    !! No geometry, kernel, constitutive law, or source normalization is
    !! inferred here.  The reciprocity helper only compares two supplied
    !! source/target work pairings with a positive target-side weight.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_source_trace_map
    public :: evaluate_source_trace_map_jvp
    public :: evaluate_source_trace_map_vjp
    public :: evaluate_weighted_source_trace_reciprocity_defect

contains

    subroutine evaluate_source_trace_map( &
            source_matrix, source_coefficients, trace, status)
        !! Apply a caller-owned rectangular source-to-trace matrix.
        complex(dp), intent(in) :: source_matrix(:, :), source_coefficients(:)
        complex(dp), allocatable, intent(out) :: trace(:)
        type(fortsparse_status_t), intent(out) :: status

        if (allocated(trace)) deallocate(trace)
        allocate(trace(0))
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "source-trace map received incompatible arrays")
        if (.not. valid_map_inputs(source_matrix, source_coefficients)) return

        deallocate(trace)
        allocate(trace(size(source_matrix, 1)))
        trace = matmul(source_matrix, source_coefficients)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_source_trace_map

    subroutine evaluate_source_trace_map_jvp( &
            source_matrix, source_coefficients, source_matrix_dot, &
            source_coefficients_dot, trace_dot, status)
        !! Apply the exact product-rule JVP of the source-to-trace map.
        complex(dp), intent(in) :: source_matrix(:, :), source_coefficients(:)
        complex(dp), intent(in) :: source_matrix_dot(:, :)
        complex(dp), intent(in) :: source_coefficients_dot(:)
        complex(dp), allocatable, intent(out) :: trace_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        if (allocated(trace_dot)) deallocate(trace_dot)
        allocate(trace_dot(0))
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "source-trace map JVP received incompatible arrays")
        if (.not. valid_map_inputs(source_matrix, source_coefficients)) return
        if (any(shape(source_matrix_dot) /= shape(source_matrix)) .or. &
            size(source_coefficients_dot) /= size(source_coefficients) .or. &
            .not. finite_complex(source_matrix_dot) .or. &
            .not. finite_complex(source_coefficients_dot)) return

        deallocate(trace_dot)
        allocate(trace_dot(size(source_matrix, 1)))
        trace_dot = matmul(source_matrix_dot, source_coefficients) + &
            matmul(source_matrix, source_coefficients_dot)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_source_trace_map_jvp

    subroutine evaluate_source_trace_map_vjp( &
            source_matrix, source_coefficients, trace_bar, source_matrix_bar, &
            source_coefficients_bar, status)
        !! Apply the real-part complex VJP of the source-to-trace map.
        !!
        !! The convention is
        !! `Re(sum(conjg(trace_bar)*trace_dot))`, so the returned products are
        !! `source_matrix_bar = trace_bar conjg(source_coefficients)^T` and
        !! `source_coefficients_bar = source_matrix^H trace_bar`.
        complex(dp), intent(in) :: source_matrix(:, :), source_coefficients(:)
        complex(dp), intent(in) :: trace_bar(:)
        complex(dp), allocatable, intent(out) :: source_matrix_bar(:, :)
        complex(dp), allocatable, intent(out) :: source_coefficients_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        if (allocated(source_matrix_bar)) deallocate(source_matrix_bar)
        if (allocated(source_coefficients_bar)) deallocate(source_coefficients_bar)
        allocate(source_matrix_bar(0, 0), source_coefficients_bar(0))
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "source-trace map VJP received incompatible arrays")
        if (.not. valid_map_inputs(source_matrix, source_coefficients)) return
        if (size(trace_bar) /= size(source_matrix, 1) .or. &
            .not. finite_complex(trace_bar)) return

        deallocate(source_matrix_bar, source_coefficients_bar)
        allocate(source_matrix_bar(size(source_matrix, 1), size(source_matrix, 2)), &
            source_coefficients_bar(size(source_coefficients)))
        source_matrix_bar = spread(trace_bar, 2, size(source_matrix, 2))* &
            spread(conjg(source_coefficients), 1, size(source_matrix, 1))
        source_coefficients_bar = matmul( &
            conjg(transpose(source_matrix)), trace_bar)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_source_trace_map_vjp

    subroutine evaluate_weighted_source_trace_reciprocity_defect( &
            source_matrix, source_one, target_one, source_two, target_two, &
            work_weights, reciprocity_defect, status)
        !! Compare reciprocal work pairings for two source/target experiments.
        !!
        !! For `trace_k = source_matrix source_k`, this returns
        !!
        !!   | target_1^T W trace_2 - target_2^T W trace_1 |
        !!   ------------------------------------------------------------,
        !!   max(1, |target_1^T W trace_2|, |target_2^T W trace_1|).
        !!
        !! The transpose pairing is intentional: this is a reciprocity
        !! diagnostic, not a Hermitian-adjoint or passivity certificate.
        complex(dp), intent(in) :: source_matrix(:, :)
        complex(dp), intent(in) :: source_one(:), target_one(:)
        complex(dp), intent(in) :: source_two(:), target_two(:)
        real(dp), intent(in) :: work_weights(:)
        real(dp), intent(out) :: reciprocity_defect
        type(fortsparse_status_t), intent(out) :: status

        complex(dp), allocatable :: trace_one(:), trace_two(:)
        complex(dp) :: first_pairing, second_pairing
        real(dp) :: scale

        reciprocity_defect = huge(1.0_dp)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "weighted source-trace reciprocity received incompatible arrays")
        if (.not. valid_map_inputs(source_matrix, source_one)) return
        if (size(source_two) /= size(source_one) .or. &
            size(target_one) /= size(source_matrix, 1) .or. &
            size(target_two) /= size(target_one) .or. &
            size(work_weights) /= size(target_one)) return
        if (.not. finite_complex(source_two) .or. &
            .not. finite_complex(target_one) .or. &
            .not. finite_complex(target_two) .or. &
            .not. all(ieee_is_finite(work_weights)) .or. &
            any(work_weights <= 0.0_dp)) return

        allocate(trace_one(size(source_matrix, 1)), trace_two(size(source_matrix, 1)))
        trace_one = matmul(source_matrix, source_one)
        trace_two = matmul(source_matrix, source_two)
        first_pairing = sum(target_one*work_weights*trace_two)
        second_pairing = sum(target_two*work_weights*trace_one)
        scale = max(1.0_dp, abs(first_pairing), abs(second_pairing))
        reciprocity_defect = abs(first_pairing - second_pairing)/scale
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_weighted_source_trace_reciprocity_defect

    logical function valid_map_inputs(source_matrix, source_coefficients) result(valid)
        complex(dp), intent(in) :: source_matrix(:, :), source_coefficients(:)

        valid = size(source_matrix, 1) > 0 .and. size(source_matrix, 2) > 0 .and. &
            size(source_coefficients) == size(source_matrix, 2) .and. &
            finite_complex(source_matrix) .and. finite_complex(source_coefficients)
    end function valid_map_inputs

    logical function finite_complex(values) result(valid)
        complex(dp), intent(in) :: values(..)

        valid = .true.
        select rank (values)
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

end module fortfem_source_trace_map
