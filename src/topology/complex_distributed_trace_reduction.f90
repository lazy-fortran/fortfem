module fortfem_complex_distributed_trace_reduction
    !! Complex value/products for the neutral distributed trace ledger.
    !!
    !! The fixed ownership topology is supplied by
    !! ``distributed_trace_layout_t``.  This module applies the same global
    !! sum and unique-owner aggregation as the real reduction to complex trace
    !! components, without introducing a communicator or duplicating layout
    !! metadata.  Its VJP uses the real-part complex adjoint convention.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_distributed_trace_ownership, only: &
        assemble_distributed_trace_reduction, &
        assemble_distributed_trace_reduction_vjp, &
        distributed_trace_layout_component_count, &
        distributed_trace_layout_global_count, &
        distributed_trace_layout_local_count, distributed_trace_layout_t, &
        validate_distributed_trace_layout
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: assemble_complex_distributed_trace_reduction
    public :: assemble_complex_distributed_trace_reduction_jvp
    public :: assemble_complex_distributed_trace_reduction_vjp

contains

    subroutine assemble_complex_distributed_trace_reduction( &
            layout, local_values, global_values, owner_values, status)
        !! Sum complex local rows and separately retain unique owner rows.
        type(distributed_trace_layout_t), intent(in) :: layout
        complex(dp), intent(in) :: local_values(:, :)
        complex(dp), intent(out) :: global_values(:, :), owner_values(:, :)
        type(fortsparse_status_t), intent(out) :: status

        real(dp) :: global_imag(size(global_values, 1), size(global_values, 2))
        real(dp) :: global_real(size(global_values, 1), size(global_values, 2))
        real(dp) :: owner_imag(size(owner_values, 1), size(owner_values, 2))
        real(dp) :: owner_real(size(owner_values, 1), size(owner_values, 2))

        global_values = cmplx(0.0_dp, 0.0_dp, dp)
        owner_values = cmplx(0.0_dp, 0.0_dp, dp)
        call validate_reduction_inputs( &
            layout, local_values, global_values, owner_values, status)
        if (status%code /= FORTSPARSE_OK) return

        call assemble_distributed_trace_reduction( &
            layout, real(local_values, dp), global_real, owner_real, status)
        if (status%code /= FORTSPARSE_OK) return
        call assemble_distributed_trace_reduction( &
            layout, aimag(local_values), global_imag, owner_imag, status)
        if (status%code /= FORTSPARSE_OK) return
        global_values = cmplx(global_real, global_imag, dp)
        owner_values = cmplx(owner_real, owner_imag, dp)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_complex_distributed_trace_reduction

    subroutine assemble_complex_distributed_trace_reduction_jvp( &
            layout, local_values_dot, global_values_dot, owner_values_dot, &
            status)
        !! Apply the fixed-topology complex JVP.
        type(distributed_trace_layout_t), intent(in) :: layout
        complex(dp), intent(in) :: local_values_dot(:, :)
        complex(dp), intent(out) :: global_values_dot(:, :)
        complex(dp), intent(out) :: owner_values_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call assemble_complex_distributed_trace_reduction( &
            layout, local_values_dot, global_values_dot, owner_values_dot, &
            status)
    end subroutine assemble_complex_distributed_trace_reduction_jvp

    subroutine assemble_complex_distributed_trace_reduction_vjp( &
            layout, global_values_bar, owner_values_bar, local_values_bar, &
            status)
        !! Apply the real-part complex transpose of both aggregation outputs.
        type(distributed_trace_layout_t), intent(in) :: layout
        complex(dp), intent(in) :: global_values_bar(:, :)
        complex(dp), intent(in) :: owner_values_bar(:, :)
        complex(dp), intent(out) :: local_values_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        real(dp) :: local_imag(size(local_values_bar, 1), size(local_values_bar, 2))
        real(dp) :: local_real(size(local_values_bar, 1), size(local_values_bar, 2))

        local_values_bar = cmplx(0.0_dp, 0.0_dp, dp)
        call validate_vjp_inputs( &
            layout, global_values_bar, owner_values_bar, local_values_bar, &
            status)
        if (status%code /= FORTSPARSE_OK) return

        call assemble_distributed_trace_reduction_vjp( &
            layout, real(global_values_bar, dp), &
            real(owner_values_bar, dp), local_real, status)
        if (status%code /= FORTSPARSE_OK) return
        call assemble_distributed_trace_reduction_vjp( &
            layout, aimag(global_values_bar), aimag(owner_values_bar), &
            local_imag, status)
        if (status%code /= FORTSPARSE_OK) return
        local_values_bar = cmplx(local_real, local_imag, dp)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_complex_distributed_trace_reduction_vjp

    subroutine validate_reduction_inputs( &
            layout, local_values, global_values, owner_values, status)
        type(distributed_trace_layout_t), intent(in) :: layout
        complex(dp), intent(in) :: local_values(:, :)
        complex(dp), intent(in) :: global_values(:, :), owner_values(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call validate_distributed_trace_layout(layout, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(local_values, 1) /= &
            distributed_trace_layout_local_count(layout) .or. &
            size(local_values, 2) /= &
            distributed_trace_layout_component_count(layout) .or. &
            size(global_values, 1) /= &
            distributed_trace_layout_global_count(layout) .or. &
            size(global_values, 2) /= &
            distributed_trace_layout_component_count(layout) .or. &
            size(owner_values, 1) /= &
            distributed_trace_layout_global_count(layout) .or. &
            size(owner_values, 2) /= &
            distributed_trace_layout_component_count(layout)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "complex distributed trace reduction received incompatible arrays")
            return
        end if
        if (.not. finite_complex_matrix(local_values)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "complex distributed trace reduction received non-finite values")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_reduction_inputs

    subroutine validate_vjp_inputs( &
            layout, global_values_bar, owner_values_bar, local_values_bar, &
            status)
        type(distributed_trace_layout_t), intent(in) :: layout
        complex(dp), intent(in) :: global_values_bar(:, :)
        complex(dp), intent(in) :: owner_values_bar(:, :)
        complex(dp), intent(in) :: local_values_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call validate_distributed_trace_layout(layout, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(global_values_bar, 1) /= &
            distributed_trace_layout_global_count(layout) .or. &
            size(global_values_bar, 2) /= &
            distributed_trace_layout_component_count(layout) .or. &
            size(owner_values_bar, 1) /= &
            distributed_trace_layout_global_count(layout) .or. &
            size(owner_values_bar, 2) /= &
            distributed_trace_layout_component_count(layout) .or. &
            size(local_values_bar, 1) /= &
            distributed_trace_layout_local_count(layout) .or. &
            size(local_values_bar, 2) /= &
            distributed_trace_layout_component_count(layout)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "complex distributed trace VJP received incompatible arrays")
            return
        end if
        if (.not. finite_complex_matrix(global_values_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "complex distributed trace VJP received non-finite adjoints")
            return
        end if
        if (.not. finite_complex_matrix(owner_values_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "complex distributed trace VJP received non-finite adjoints")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_vjp_inputs

    logical function finite_complex_matrix(values) result(valid)
        complex(dp), intent(in) :: values(:, :)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex_matrix

end module fortfem_complex_distributed_trace_reduction
