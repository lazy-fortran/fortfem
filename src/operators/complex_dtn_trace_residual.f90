module fortfem_complex_dtn_trace_residual
    !! Neutral complex FEM/DtN source-to-trace coupling residual.
    !!
    !! For a caller-owned compact source-to-trace map S, exterior source
    !! coefficients q, a FEM trace t, and positive work weights w, this module
    !! evaluates
    !!
    !!   r = w (t - S q).
    !!
    !! The ``DtN`` label identifies the coupling role only.  No geometry,
    !! differential equation, basis normalization, or exterior operator is
    !! selected here; callers may supply any compatible complex trace map.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: assemble_complex_dtn_trace_residual
    public :: assemble_complex_dtn_trace_residual_jvp
    public :: assemble_complex_dtn_trace_residual_vjp

    interface finite_complex
        module procedure finite_complex_vector
        module procedure finite_complex_matrix
    end interface finite_complex

contains

    subroutine assemble_complex_dtn_trace_residual( &
            source_to_trace, source_coefficients, fem_trace, work_weights, &
            residual, status)
        complex(dp), intent(in) :: source_to_trace(:, :)
        complex(dp), intent(in) :: source_coefficients(:), fem_trace(:)
        real(dp), intent(in) :: work_weights(:)
        complex(dp), intent(out) :: residual(:)
        type(fortsparse_status_t), intent(out) :: status

        residual = cmplx(0.0_dp, 0.0_dp, dp)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "complex FEM/DtN trace residual received incompatible arrays")
        if (.not. valid_primal( &
            source_to_trace, source_coefficients, fem_trace, work_weights, &
            residual)) return

        residual = work_weights*( &
            fem_trace - matmul(source_to_trace, source_coefficients))
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_complex_dtn_trace_residual

    subroutine assemble_complex_dtn_trace_residual_jvp( &
            source_to_trace, source_coefficients, fem_trace, work_weights, &
            source_to_trace_dot, source_coefficients_dot, fem_trace_dot, &
            work_weights_dot, residual_dot, status)
        complex(dp), intent(in) :: source_to_trace(:, :)
        complex(dp), intent(in) :: source_coefficients(:), fem_trace(:)
        real(dp), intent(in) :: work_weights(:)
        complex(dp), intent(in) :: source_to_trace_dot(:, :)
        complex(dp), intent(in) :: source_coefficients_dot(:), fem_trace_dot(:)
        real(dp), intent(in) :: work_weights_dot(:)
        complex(dp), intent(out) :: residual_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        complex(dp) :: raw_trace(size(fem_trace))

        residual_dot = cmplx(0.0_dp, 0.0_dp, dp)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "complex FEM/DtN trace residual JVP received incompatible arrays")
        if (.not. valid_primal( &
            source_to_trace, source_coefficients, fem_trace, work_weights, &
            residual_dot)) return
        if (.not. valid_direction( &
            source_to_trace, source_coefficients, fem_trace, work_weights, &
            source_to_trace_dot, source_coefficients_dot, fem_trace_dot, &
            work_weights_dot)) return

        raw_trace = fem_trace - matmul(source_to_trace, source_coefficients)
        residual_dot = work_weights_dot*raw_trace + work_weights*( &
            fem_trace_dot - matmul(source_to_trace_dot, source_coefficients) - &
            matmul(source_to_trace, source_coefficients_dot))
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_complex_dtn_trace_residual_jvp

    subroutine assemble_complex_dtn_trace_residual_vjp( &
            source_to_trace, source_coefficients, fem_trace, work_weights, &
            residual_bar, source_to_trace_bar, source_coefficients_bar, &
            fem_trace_bar, work_weights_bar, status)
        !! Return the VJP for the pairing Re(residual_bar**H residual_dot).
        complex(dp), intent(in) :: source_to_trace(:, :)
        complex(dp), intent(in) :: source_coefficients(:), fem_trace(:)
        real(dp), intent(in) :: work_weights(:)
        complex(dp), intent(in) :: residual_bar(:)
        complex(dp), intent(out) :: source_to_trace_bar(:, :)
        complex(dp), intent(out) :: source_coefficients_bar(:)
        complex(dp), intent(out) :: fem_trace_bar(:)
        real(dp), intent(out) :: work_weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        complex(dp) :: raw_trace(size(fem_trace))
        complex(dp) :: weighted_bar(size(fem_trace))
        integer :: column, row

        source_to_trace_bar = cmplx(0.0_dp, 0.0_dp, dp)
        source_coefficients_bar = cmplx(0.0_dp, 0.0_dp, dp)
        fem_trace_bar = cmplx(0.0_dp, 0.0_dp, dp)
        work_weights_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "complex FEM/DtN trace residual VJP received incompatible arrays")
        if (.not. valid_primal( &
            source_to_trace, source_coefficients, fem_trace, work_weights, &
            residual_bar)) return
        if (.not. valid_adjoint( &
            source_to_trace, source_coefficients, fem_trace, work_weights, &
            residual_bar, source_to_trace_bar, source_coefficients_bar, &
            fem_trace_bar, work_weights_bar)) return

        raw_trace = fem_trace - matmul(source_to_trace, source_coefficients)
        weighted_bar = work_weights*residual_bar
        fem_trace_bar = weighted_bar
        work_weights_bar = real(conjg(residual_bar)*raw_trace, dp)
        do column = 1, size(source_coefficients)
            do row = 1, size(fem_trace)
                source_to_trace_bar(row, column) = &
                    -weighted_bar(row)*conjg(source_coefficients(column))
            end do
        end do
        source_coefficients_bar = -matmul( &
            conjg(transpose(source_to_trace)), weighted_bar)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_complex_dtn_trace_residual_vjp

    logical function valid_primal( &
            source_to_trace, source_coefficients, fem_trace, work_weights, &
            residual) result(valid)
        complex(dp), intent(in) :: source_to_trace(:, :)
        complex(dp), intent(in) :: source_coefficients(:), fem_trace(:)
        real(dp), intent(in) :: work_weights(:)
        complex(dp), intent(in) :: residual(:)

        valid = .false.
        if (size(source_to_trace, 1) < 1 .or. &
            size(source_to_trace, 2) < 1) return
        if (size(source_coefficients) /= size(source_to_trace, 2)) return
        if (size(fem_trace) /= size(source_to_trace, 1)) return
        if (size(work_weights) /= size(fem_trace)) return
        if (size(residual) /= size(fem_trace)) return
        if (any(work_weights <= 0.0_dp)) return
        if (.not. all(ieee_is_finite(work_weights))) return
        if (.not. finite_complex(source_to_trace)) return
        if (.not. finite_complex(source_coefficients)) return
        if (.not. finite_complex(fem_trace)) return
        valid = .true.
    end function valid_primal

    logical function valid_direction( &
            source_to_trace, source_coefficients, fem_trace, work_weights, &
            source_to_trace_dot, source_coefficients_dot, fem_trace_dot, &
            work_weights_dot) result(valid)
        complex(dp), intent(in) :: source_to_trace(:, :)
        complex(dp), intent(in) :: source_coefficients(:), fem_trace(:)
        real(dp), intent(in) :: work_weights(:)
        complex(dp), intent(in) :: source_to_trace_dot(:, :)
        complex(dp), intent(in) :: source_coefficients_dot(:), fem_trace_dot(:)
        real(dp), intent(in) :: work_weights_dot(:)

        valid = .false.
        if (any(shape(source_to_trace_dot) /= shape(source_to_trace))) return
        if (size(source_coefficients_dot) /= size(source_coefficients)) return
        if (size(fem_trace_dot) /= size(fem_trace)) return
        if (size(work_weights_dot) /= size(work_weights)) return
        if (.not. finite_complex(source_to_trace_dot)) return
        if (.not. finite_complex(source_coefficients_dot)) return
        if (.not. finite_complex(fem_trace_dot)) return
        if (.not. all(ieee_is_finite(work_weights_dot))) return
        valid = .true.
    end function valid_direction

    logical function valid_adjoint( &
            source_to_trace, source_coefficients, fem_trace, work_weights, &
            residual_bar, source_to_trace_bar, source_coefficients_bar, &
            fem_trace_bar, work_weights_bar) result(valid)
        complex(dp), intent(in) :: source_to_trace(:, :)
        complex(dp), intent(in) :: source_coefficients(:), fem_trace(:)
        real(dp), intent(in) :: work_weights(:)
        complex(dp), intent(in) :: residual_bar(:)
        complex(dp), intent(in) :: source_to_trace_bar(:, :)
        complex(dp), intent(in) :: source_coefficients_bar(:), fem_trace_bar(:)
        real(dp), intent(in) :: work_weights_bar(:)

        valid = .false.
        if (.not. finite_complex(residual_bar)) return
        if (any(shape(source_to_trace_bar) /= shape(source_to_trace))) return
        if (size(source_coefficients_bar) /= size(source_coefficients)) return
        if (size(fem_trace_bar) /= size(fem_trace)) return
        if (size(work_weights_bar) /= size(work_weights)) return
        valid = .true.
    end function valid_adjoint

    logical function finite_complex_vector(values) result(valid)
        complex(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex_vector

    logical function finite_complex_matrix(values) result(valid)
        complex(dp), intent(in) :: values(:, :)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex_matrix

end module fortfem_complex_dtn_trace_residual
