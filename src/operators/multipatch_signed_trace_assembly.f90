module fortfem_multipatch_signed_trace_assembly
    !! Validated multipatch wrapper for signed trace block assembly.
    !!
    !! The map is fixed topology: each patch contributes a local matrix and
    !! right-hand side through a signed local-to-global trace map.  Assembly is
    !! delegated to the dense HDG reference kernel after this layer checks the
    !! stronger multipatch invariant that one patch cannot list a global trace
    !! ID twice.  No communicator, mesh, or constitutive model is selected.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_hdg_global_skeleton, only: &
        assemble_hdg_global_skeleton, assemble_hdg_global_skeleton_jvp, &
        assemble_hdg_global_skeleton_vjp
    use fortfem_hdg_global_skeleton_csc, only: &
        assemble_hdg_global_skeleton_csc, assemble_hdg_global_skeleton_csc_jvp, &
        assemble_hdg_global_skeleton_csc_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        csc_is_valid, csc_t, FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_multipatch_signed_trace_assembly
    public :: assemble_multipatch_signed_trace_assembly_jvp
    public :: assemble_multipatch_signed_trace_assembly_vjp
    public :: assemble_multipatch_signed_trace_assembly_csc
    public :: assemble_multipatch_signed_trace_assembly_csc_jvp
    public :: assemble_multipatch_signed_trace_assembly_csc_vjp

contains

    subroutine assemble_multipatch_signed_trace_assembly( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix, global_rhs, status)
        real(dp), intent(in) :: local_matrix(:, :, :), local_rhs(:, :)
        integer, intent(in) :: local_to_global(:, :), local_sign(:, :)
        integer, intent(in) :: global_count
        real(dp), intent(out) :: global_matrix(:, :), global_rhs(:)
        type(fortsparse_status_t), intent(out) :: status

        global_matrix = 0.0_dp
        global_rhs = 0.0_dp
        if (.not. valid_signed_patch_inputs( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix, global_rhs, status)) return
        call assemble_hdg_global_skeleton( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix, global_rhs, status)
    end subroutine assemble_multipatch_signed_trace_assembly

    subroutine assemble_multipatch_signed_trace_assembly_jvp( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            local_matrix_dot, local_rhs_dot, global_matrix_dot, global_rhs_dot, status)
        real(dp), intent(in) :: local_matrix(:, :, :), local_rhs(:, :)
        integer, intent(in) :: local_to_global(:, :), local_sign(:, :)
        integer, intent(in) :: global_count
        real(dp), intent(in) :: local_matrix_dot(:, :, :), local_rhs_dot(:, :)
        real(dp), intent(out) :: global_matrix_dot(:, :), global_rhs_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        global_matrix_dot = 0.0_dp
        global_rhs_dot = 0.0_dp
        if (.not. valid_signed_patch_inputs( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix_dot, global_rhs_dot, status)) return
        call assemble_hdg_global_skeleton_jvp( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            local_matrix_dot, local_rhs_dot, global_matrix_dot, global_rhs_dot, status)
    end subroutine assemble_multipatch_signed_trace_assembly_jvp

    subroutine assemble_multipatch_signed_trace_assembly_vjp( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix_bar, global_rhs_bar, local_matrix_bar, local_rhs_bar, status)
        real(dp), intent(in) :: local_matrix(:, :, :), local_rhs(:, :)
        integer, intent(in) :: local_to_global(:, :), local_sign(:, :)
        integer, intent(in) :: global_count
        real(dp), intent(in) :: global_matrix_bar(:, :), global_rhs_bar(:)
        real(dp), intent(out) :: local_matrix_bar(:, :, :), local_rhs_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        local_matrix_bar = 0.0_dp
        local_rhs_bar = 0.0_dp
        if (.not. valid_signed_patch_inputs( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix_bar, global_rhs_bar, status)) return
        call assemble_hdg_global_skeleton_vjp( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix_bar, global_rhs_bar, local_matrix_bar, local_rhs_bar, status)
    end subroutine assemble_multipatch_signed_trace_assembly_vjp

    subroutine assemble_multipatch_signed_trace_assembly_csc( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix, global_rhs, status)
        real(dp), intent(in) :: local_matrix(:, :, :), local_rhs(:, :)
        integer, intent(in) :: local_to_global(:, :), local_sign(:, :)
        integer, intent(in) :: global_count
        type(csc_t), intent(out) :: global_matrix
        real(dp), intent(out) :: global_rhs(:)
        type(fortsparse_status_t), intent(out) :: status

        global_rhs = 0.0_dp
        if (.not. valid_signed_patch_csc_inputs( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_rhs, status)) return
        call assemble_hdg_global_skeleton_csc( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix, global_rhs, status)
    end subroutine assemble_multipatch_signed_trace_assembly_csc

    subroutine assemble_multipatch_signed_trace_assembly_csc_jvp( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            local_matrix_dot, local_rhs_dot, global_matrix_dot, global_rhs_dot, status)
        real(dp), intent(in) :: local_matrix(:, :, :), local_rhs(:, :)
        integer, intent(in) :: local_to_global(:, :), local_sign(:, :)
        integer, intent(in) :: global_count
        real(dp), intent(in) :: local_matrix_dot(:, :, :), local_rhs_dot(:, :)
        type(csc_t), intent(out) :: global_matrix_dot
        real(dp), intent(out) :: global_rhs_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        global_rhs_dot = 0.0_dp
        if (.not. valid_signed_patch_csc_inputs( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_rhs_dot, status)) return
        call assemble_hdg_global_skeleton_csc_jvp( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            local_matrix_dot, local_rhs_dot, global_matrix_dot, global_rhs_dot, status)
    end subroutine assemble_multipatch_signed_trace_assembly_csc_jvp

    subroutine assemble_multipatch_signed_trace_assembly_csc_vjp( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix_bar, global_rhs_bar, local_matrix_bar, local_rhs_bar, status)
        real(dp), intent(in) :: local_matrix(:, :, :), local_rhs(:, :)
        integer, intent(in) :: local_to_global(:, :), local_sign(:, :)
        integer, intent(in) :: global_count
        type(csc_t), intent(in) :: global_matrix_bar
        real(dp), intent(in) :: global_rhs_bar(:)
        real(dp), intent(out) :: local_matrix_bar(:, :, :), local_rhs_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        local_matrix_bar = 0.0_dp
        local_rhs_bar = 0.0_dp
        if (.not. valid_signed_patch_csc_inputs( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_rhs_bar, status)) return
        if (.not. csc_is_valid(global_matrix_bar) .or. &
            global_matrix_bar%nrow /= global_count .or. &
            global_matrix_bar%ncol /= global_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "multipatch signed CSC VJP received an incompatible matrix cotangent")
            return
        end if
        call assemble_hdg_global_skeleton_csc_vjp( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix_bar, global_rhs_bar, local_matrix_bar, local_rhs_bar, status)
    end subroutine assemble_multipatch_signed_trace_assembly_csc_vjp

    logical function valid_signed_patch_csc_inputs( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_rhs, status) result(valid)
        real(dp), intent(in) :: local_matrix(:, :, :), local_rhs(:, :)
        integer, intent(in) :: local_to_global(:, :), local_sign(:, :)
        integer, intent(in) :: global_count
        real(dp), intent(in) :: global_rhs(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: local_count, patch_count, patch, local, other

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "multipatch signed CSC assembly received incompatible arrays")
        local_count = size(local_matrix, 1)
        patch_count = size(local_matrix, 3)
        if (global_count < 1 .or. local_count < 1 .or. patch_count < 1 .or. &
            size(local_matrix, 2) /= local_count .or. &
            size(local_rhs, 1) /= local_count .or. size(local_rhs, 2) /= patch_count .or. &
            size(local_to_global, 1) /= local_count .or. &
            size(local_to_global, 2) /= patch_count .or. &
            size(local_sign, 1) /= local_count .or. size(local_sign, 2) /= patch_count .or. &
            size(global_rhs) /= global_count) return
        if (any(local_to_global < 1) .or. any(local_to_global > global_count) .or. &
            any(abs(local_sign) /= 1) .or. any(.not. ieee_is_finite(local_matrix)) .or. &
            any(.not. ieee_is_finite(local_rhs)) .or. any(.not. ieee_is_finite(global_rhs))) return
        do patch = 1, patch_count
            do local = 1, local_count
                do other = local + 1, local_count
                    if (local_to_global(local, patch) == local_to_global(other, patch)) return
                end do
            end do
        end do
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function valid_signed_patch_csc_inputs

    logical function valid_signed_patch_inputs( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix, global_rhs, status) result(valid)
        real(dp), intent(in) :: local_matrix(:, :, :), local_rhs(:, :)
        integer, intent(in) :: local_to_global(:, :), local_sign(:, :)
        integer, intent(in) :: global_count
        real(dp), intent(in) :: global_matrix(:, :), global_rhs(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: local_count, patch_count, patch, local, other

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "multipatch signed trace assembly received incompatible arrays")
        local_count = size(local_matrix, 1)
        patch_count = size(local_matrix, 3)
        if (global_count < 1 .or. local_count < 1 .or. patch_count < 1 .or. &
            size(local_matrix, 2) /= local_count .or. &
            size(local_rhs, 1) /= local_count .or. size(local_rhs, 2) /= patch_count .or. &
            size(local_to_global, 1) /= local_count .or. &
            size(local_to_global, 2) /= patch_count .or. &
            size(local_sign, 1) /= local_count .or. size(local_sign, 2) /= patch_count .or. &
            size(global_matrix, 1) /= global_count .or. &
            size(global_matrix, 2) /= global_count .or. size(global_rhs) /= global_count) return
        if (any(local_to_global < 1) .or. any(local_to_global > global_count) .or. &
            any(abs(local_sign) /= 1)) return
        if (any(.not. ieee_is_finite(local_matrix)) .or. &
            any(.not. ieee_is_finite(local_rhs)) .or. &
            any(.not. ieee_is_finite(global_matrix)) .or. &
            any(.not. ieee_is_finite(global_rhs))) return
        do patch = 1, patch_count
            do local = 1, local_count
                do other = local + 1, local_count
                    if (local_to_global(local, patch) == local_to_global(other, patch)) return
                end do
            end do
        end do
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function valid_signed_patch_inputs

end module fortfem_multipatch_signed_trace_assembly
