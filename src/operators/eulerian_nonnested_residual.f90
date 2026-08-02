module fortfem_eulerian_nonnested_residual
    !! Closure-neutral residual contract for fixed-domain, non-nested clients.
    !!
    !! The force and divergence vectors are caller-owned samples.  An optional
    !! stabilization vector may be produced by
    !! `assemble_pseudo_transient_residual`; this layer only adds it to the
    !! fixed-domain residual and therefore does not choose a mass matrix or a
    !! relaxation closure.  Optional signed margins are classified by the
    !! existing continuation-event diagnostic.  Event outputs are discrete
    !! metadata and are intentionally held fixed by JVP/VJP actions.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_continuation_event, only: classify_continuation_event
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_eulerian_nonnested_residual
    public :: assemble_eulerian_nonnested_residual_jvp
    public :: assemble_eulerian_nonnested_residual_vjp

contains

    subroutine assemble_eulerian_nonnested_residual( &
            force_residual, divergence_residual, residual, status, &
            stabilization_residual, previous_margin, current_margin, event_tolerance, &
            event_code, event_index, minimum_margin)
        real(dp), intent(in) :: force_residual(:), divergence_residual(:)
        real(dp), intent(out) :: residual(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: stabilization_residual(:)
        real(dp), intent(in), optional :: previous_margin(:), current_margin(:)
        real(dp), intent(in), optional :: event_tolerance
        integer, intent(out), optional :: event_code, event_index
        real(dp), intent(out), optional :: minimum_margin
        integer :: local_event_code, local_event_index
        real(dp) :: local_minimum_margin

        call initialize_event_outputs( &
            local_event_code, local_event_index, local_minimum_margin, &
            event_code, event_index, minimum_margin)
        residual = 0.0_dp
        if (.not. validate_value_inputs( &
            force_residual, divergence_residual, residual, stabilization_residual, &
            previous_margin, current_margin, event_tolerance, status)) return
        residual(1:size(force_residual)) = force_residual
        residual(size(force_residual)+1:) = divergence_residual
        if (present(stabilization_residual)) then
            residual = residual + stabilization_residual
        end if
        if (.not. classify_margins( &
            previous_margin, current_margin, event_tolerance, local_event_code, &
            local_event_index, local_minimum_margin, status)) then
            residual = 0.0_dp
            return
        end if
        call copy_event_outputs( &
            local_event_code, local_event_index, local_minimum_margin, &
            event_code, event_index, minimum_margin)
        if (.not. all(ieee_is_finite(residual))) then
            residual = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Eulerian non-nested residual is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_eulerian_nonnested_residual

    subroutine assemble_eulerian_nonnested_residual_jvp( &
            force_residual, divergence_residual, force_residual_dot, &
            divergence_residual_dot, residual_dot, status, stabilization_residual, &
            stabilization_residual_dot, previous_margin, current_margin, &
            event_tolerance, &
            event_code, event_index, minimum_margin)
        real(dp), intent(in) :: force_residual(:), divergence_residual(:)
        real(dp), intent(in) :: force_residual_dot(:), divergence_residual_dot(:)
        real(dp), intent(out) :: residual_dot(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: stabilization_residual(:)
        real(dp), intent(in), optional :: stabilization_residual_dot(:)
        real(dp), intent(in), optional :: previous_margin(:), current_margin(:)
        real(dp), intent(in), optional :: event_tolerance
        integer, intent(out), optional :: event_code, event_index
        real(dp), intent(out), optional :: minimum_margin
        integer :: local_event_code, local_event_index
        real(dp) :: local_minimum_margin

        call initialize_event_outputs( &
            local_event_code, local_event_index, local_minimum_margin, &
            event_code, event_index, minimum_margin)
        residual_dot = 0.0_dp
        if (.not. validate_jvp_inputs( &
            force_residual, divergence_residual, force_residual_dot, &
            divergence_residual_dot, residual_dot, stabilization_residual, &
            stabilization_residual_dot, previous_margin, current_margin, &
            event_tolerance, &
            status)) return
        residual_dot(1:size(force_residual)) = force_residual_dot
        residual_dot(size(force_residual)+1:) = divergence_residual_dot
        if (present(stabilization_residual_dot)) residual_dot = residual_dot + &
            stabilization_residual_dot
        if (.not. classify_margins( &
            previous_margin, current_margin, event_tolerance, local_event_code, &
            local_event_index, local_minimum_margin, status)) then
            residual_dot = 0.0_dp
            return
        end if
        call copy_event_outputs( &
            local_event_code, local_event_index, local_minimum_margin, &
            event_code, event_index, minimum_margin)
        if (.not. all(ieee_is_finite(residual_dot))) then
            residual_dot = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Eulerian non-nested residual JVP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_eulerian_nonnested_residual_jvp

    subroutine assemble_eulerian_nonnested_residual_vjp( &
            force_residual, divergence_residual, residual_bar, force_residual_bar, &
            divergence_residual_bar, status, stabilization_residual, &
            stabilization_residual_bar)
        real(dp), intent(in) :: force_residual(:), divergence_residual(:)
        real(dp), intent(in) :: residual_bar(:)
        real(dp), intent(out) :: force_residual_bar(:), divergence_residual_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: stabilization_residual(:)
        real(dp), intent(out), optional :: stabilization_residual_bar(:)
        integer :: force_count, total_count

        force_residual_bar = 0.0_dp
        divergence_residual_bar = 0.0_dp
        if (present(stabilization_residual_bar)) stabilization_residual_bar = 0.0_dp
        force_count = size(force_residual)
        total_count = force_count + size(divergence_residual)
        if (.not. validate_vjp_inputs( &
            force_residual, divergence_residual, residual_bar, force_residual_bar, &
            divergence_residual_bar, stabilization_residual, &
            stabilization_residual_bar, &
            status)) return
        force_residual_bar = residual_bar(:force_count)
        divergence_residual_bar = residual_bar(force_count+1:total_count)
        if (present(stabilization_residual_bar)) then
            stabilization_residual_bar = residual_bar
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_eulerian_nonnested_residual_vjp

    subroutine initialize_event_outputs( &
            event_code_local, event_index_local, minimum_margin_local, &
            event_code, event_index, minimum_margin)
        integer, intent(out) :: event_code_local, event_index_local
        real(dp), intent(out) :: minimum_margin_local
        integer, intent(out), optional :: event_code, event_index
        real(dp), intent(out), optional :: minimum_margin

        event_code_local = 0
        event_index_local = 0
        minimum_margin_local = 0.0_dp
        if (present(event_code)) event_code = event_code_local
        if (present(event_index)) event_index = event_index_local
        if (present(minimum_margin)) minimum_margin = minimum_margin_local
    end subroutine initialize_event_outputs

    subroutine copy_event_outputs( &
            event_code_local, event_index_local, minimum_margin_local, &
            event_code, event_index, minimum_margin)
        integer, intent(in) :: event_code_local, event_index_local
        real(dp), intent(in) :: minimum_margin_local
        integer, intent(out), optional :: event_code, event_index
        real(dp), intent(out), optional :: minimum_margin

        if (present(event_code)) event_code = event_code_local
        if (present(event_index)) event_index = event_index_local
        if (present(minimum_margin)) minimum_margin = minimum_margin_local
    end subroutine copy_event_outputs

    logical function classify_margins( &
            previous_margin, current_margin, event_tolerance, event_code, event_index, &
            minimum_margin, status) result(valid)
        real(dp), intent(in), optional :: previous_margin(:), current_margin(:)
        real(dp), intent(in), optional :: event_tolerance
        integer, intent(out) :: event_code, event_index
        real(dp), intent(out) :: minimum_margin
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: tolerance

        event_code = 0
        event_index = 0
        minimum_margin = 0.0_dp
        if (present(previous_margin) .neqv. present(current_margin)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Eulerian topology margins must be supplied as a pair")
            valid = .false.
            return
        end if
        if (.not. present(previous_margin)) then
            if (present(event_tolerance)) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "Eulerian event tolerance requires topology margins")
                valid = .false.
                return
            end if
            call status_set(status, FORTSPARSE_OK, "")
            valid = .true.
            return
        end if
        if (present(event_tolerance)) then
            tolerance = event_tolerance
        else
            tolerance = 0.0_dp
        end if
        call classify_continuation_event( &
            previous_margin, current_margin, tolerance, event_code, event_index, &
            minimum_margin, status)
        valid = status%code == FORTSPARSE_OK
    end function classify_margins

    logical function validate_value_inputs( &
            force_residual, divergence_residual, residual, stabilization_residual, &
            previous_margin, current_margin, event_tolerance, status) result(valid)
        real(dp), intent(in) :: force_residual(:), divergence_residual(:), residual(:)
        real(dp), intent(in), optional :: stabilization_residual(:)
        real(dp), intent(in), optional :: previous_margin(:), current_margin(:)
        real(dp), intent(in), optional :: event_tolerance
        type(fortsparse_status_t), intent(out) :: status
        integer :: total_count

        total_count = size(force_residual) + size(divergence_residual)
        valid = total_count > 0 .and. size(residual) == total_count
        if (valid) then
            if (present(stabilization_residual)) valid = &
                size(stabilization_residual) == total_count
        end if
        if (.not. valid) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Eulerian residual vectors have incompatible lengths")
            return
        end if
        valid = all(ieee_is_finite(force_residual)) .and. &
            all(ieee_is_finite(divergence_residual))
        if (valid) then
            if (present(stabilization_residual)) valid = &
                all(ieee_is_finite(stabilization_residual))
        end if
        if (.not. valid) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Eulerian residual vectors contain non-finite values")
            return
        end if
        if (present(event_tolerance)) then
            valid = ieee_is_finite(event_tolerance) .and. event_tolerance >= 0.0_dp
            if (.not. valid) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "Eulerian event tolerance is invalid")
                return
            end if
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_value_inputs

    logical function validate_jvp_inputs( &
            force_residual, divergence_residual, force_residual_dot, &
            divergence_residual_dot, residual_dot, stabilization_residual, &
            stabilization_residual_dot, previous_margin, current_margin, &
            event_tolerance, &
            status) result(valid)
        real(dp), intent(in) :: force_residual(:), divergence_residual(:)
        real(dp), intent(in) :: force_residual_dot(:), divergence_residual_dot(:)
        real(dp), intent(in) :: residual_dot(:)
        real(dp), intent(in), optional :: stabilization_residual(:)
        real(dp), intent(in), optional :: stabilization_residual_dot(:)
        real(dp), intent(in), optional :: previous_margin(:), current_margin(:)
        real(dp), intent(in), optional :: event_tolerance
        type(fortsparse_status_t), intent(out) :: status
        integer :: total_count

        total_count = size(force_residual) + size(divergence_residual)
        valid = total_count > 0 .and. size(residual_dot) == total_count .and. &
            size(force_residual_dot) == size(force_residual) .and. &
            size(divergence_residual_dot) == size(divergence_residual)
        if (valid .neqv. (present(stabilization_residual) .eqv. &
            present(stabilization_residual_dot))) valid = .false.
        if (valid) then
            if (present(stabilization_residual)) valid = &
                size(stabilization_residual) == total_count .and. &
                size(stabilization_residual_dot) == total_count
        end if
        if (.not. valid) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Eulerian residual JVP vectors have incompatible lengths")
            return
        end if
        valid = all(ieee_is_finite(force_residual)) .and. &
            all(ieee_is_finite(divergence_residual)) .and. &
            all(ieee_is_finite(force_residual_dot)) .and. &
            all(ieee_is_finite(divergence_residual_dot))
        if (valid) then
            if (present(stabilization_residual)) valid = &
                all(ieee_is_finite(stabilization_residual)) .and. &
                all(ieee_is_finite(stabilization_residual_dot))
        end if
        if (.not. valid) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Eulerian residual JVP vectors contain non-finite values")
            return
        end if
        if (present(event_tolerance)) then
            valid = ieee_is_finite(event_tolerance) .and. event_tolerance >= 0.0_dp
            if (.not. valid) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "Eulerian event tolerance is invalid")
                return
            end if
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_jvp_inputs

    logical function validate_vjp_inputs( &
            force_residual, divergence_residual, residual_bar, force_residual_bar, &
            divergence_residual_bar, stabilization_residual, &
            stabilization_residual_bar, &
            status) result(valid)
        real(dp), intent(in) :: force_residual(:), divergence_residual(:)
        real(dp), intent(in) :: residual_bar(:)
        real(dp), intent(in) :: force_residual_bar(:), divergence_residual_bar(:)
        real(dp), intent(in), optional :: stabilization_residual(:)
        real(dp), intent(in), optional :: stabilization_residual_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: total_count

        total_count = size(force_residual) + size(divergence_residual)
        valid = total_count > 0 .and. size(residual_bar) == total_count .and. &
            size(force_residual_bar) == size(force_residual) .and. &
            size(divergence_residual_bar) == size(divergence_residual)
        if (valid) then
            if (present(stabilization_residual)) valid = &
                size(stabilization_residual) == total_count
        end if
        if (valid) then
            if (present(stabilization_residual_bar)) valid = &
                size(stabilization_residual_bar) == total_count
        end if
        if (.not. valid) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Eulerian residual VJP vectors have incompatible lengths")
            return
        end if
        valid = all(ieee_is_finite(force_residual)) .and. &
            all(ieee_is_finite(divergence_residual)) .and. &
            all(ieee_is_finite(residual_bar))
        if (valid) then
            if (present(stabilization_residual)) valid = &
                all(ieee_is_finite(stabilization_residual))
        end if
        if (.not. valid) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Eulerian residual VJP vectors contain non-finite values")
            return
        end if
        if (present(stabilization_residual) .neqv. &
            present(stabilization_residual_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Eulerian stabilization VJP requires a matching cotangent")
            valid = .false.
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_vjp_inputs

end module fortfem_eulerian_nonnested_residual
