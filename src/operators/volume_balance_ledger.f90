module fortfem_volume_balance_ledger
    !! Physics-agnostic volume-source balance ledger.
    !!
    !! `source_rate(cell, component)` is a caller-owned density or rate.  The
    !! positive cell measure converts it to an integrated cell contribution;
    !! the global ledger is the component-wise sum.  Units, signs, and the
    !! interpretation of each component remain with the client.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_volume_balance_ledger
    public :: assemble_volume_balance_ledger_jvp
    public :: assemble_volume_balance_ledger_vjp

contains

    subroutine assemble_volume_balance_ledger( &
            cell_weights, source_rate, cell_ledger, global_ledger, status)
        real(dp), intent(in) :: cell_weights(:), source_rate(:, :)
        real(dp), intent(out) :: cell_ledger(:, :), global_ledger(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: cell, component

        cell_ledger = 0.0_dp
        global_ledger = 0.0_dp
        if (.not. valid_value_inputs(cell_weights, source_rate, cell_ledger, &
            global_ledger, status)) return

        do cell = 1, size(cell_weights)
            do component = 1, size(source_rate, 2)
                cell_ledger(cell, component) = &
                    cell_weights(cell)*source_rate(cell, component)
                global_ledger(component) = global_ledger(component) + &
                    cell_ledger(cell, component)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_volume_balance_ledger

    subroutine assemble_volume_balance_ledger_jvp( &
            cell_weights, source_rate, cell_weights_dot, source_rate_dot, &
            cell_ledger_dot, global_ledger_dot, status)
        real(dp), intent(in) :: cell_weights(:), source_rate(:, :)
        real(dp), intent(in) :: cell_weights_dot(:), source_rate_dot(:, :)
        real(dp), intent(out) :: cell_ledger_dot(:, :), global_ledger_dot(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: cell, component

        cell_ledger_dot = 0.0_dp
        global_ledger_dot = 0.0_dp
        if (.not. valid_value_inputs(cell_weights, source_rate, cell_ledger_dot, &
            global_ledger_dot, status)) return
        if (size(cell_weights_dot) /= size(cell_weights) .or. &
            size(source_rate_dot, 1) /= size(source_rate, 1) .or. &
            size(source_rate_dot, 2) /= size(source_rate, 2) .or. &
            any(.not. ieee_is_finite(cell_weights_dot)) .or. &
            any(.not. ieee_is_finite(source_rate_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "volume balance ledger JVP has incompatible tangents")
            return
        end if

        do cell = 1, size(cell_weights)
            do component = 1, size(source_rate, 2)
                cell_ledger_dot(cell, component) = &
                    cell_weights_dot(cell)*source_rate(cell, component) + &
                    cell_weights(cell)*source_rate_dot(cell, component)
                global_ledger_dot(component) = global_ledger_dot(component) + &
                    cell_ledger_dot(cell, component)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_volume_balance_ledger_jvp

    subroutine assemble_volume_balance_ledger_vjp( &
            cell_weights, source_rate, cell_ledger_bar, global_ledger_bar, &
            cell_weights_bar, source_rate_bar, status)
        real(dp), intent(in) :: cell_weights(:), source_rate(:, :)
        real(dp), intent(in) :: cell_ledger_bar(:, :), global_ledger_bar(:)
        real(dp), intent(out) :: cell_weights_bar(:), source_rate_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: cell, component
        real(dp) :: cotangent

        cell_weights_bar = 0.0_dp
        source_rate_bar = 0.0_dp
        if (.not. valid_value_inputs(cell_weights, source_rate, cell_ledger_bar, &
            global_ledger_bar, status)) return
        if (size(cell_weights_bar) /= size(cell_weights) .or. &
            size(source_rate_bar, 1) /= size(source_rate, 1) .or. &
            size(source_rate_bar, 2) /= size(source_rate, 2) .or. &
            any(.not. ieee_is_finite(cell_ledger_bar)) .or. &
            any(.not. ieee_is_finite(global_ledger_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "volume balance ledger VJP has incompatible cotangents")
            return
        end if

        do cell = 1, size(cell_weights)
            do component = 1, size(source_rate, 2)
                cotangent = cell_ledger_bar(cell, component) + &
                    global_ledger_bar(component)
                cell_weights_bar(cell) = cell_weights_bar(cell) + &
                    cotangent*source_rate(cell, component)
                source_rate_bar(cell, component) = &
                    source_rate_bar(cell, component) + cotangent*cell_weights(cell)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_volume_balance_ledger_vjp

    logical function valid_value_inputs( &
            cell_weights, source_rate, cell_ledger, global_ledger, status)
        real(dp), intent(in) :: cell_weights(:), source_rate(:, :)
        real(dp), intent(in) :: cell_ledger(:, :), global_ledger(:)
        type(fortsparse_status_t), intent(out) :: status

        valid_value_inputs = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "volume balance ledger has incompatible arrays")
        if (size(cell_weights) < 1 .or. size(source_rate, 1) /= size(cell_weights) .or. &
            size(source_rate, 2) < 1 .or. size(cell_ledger, 1) /= size(cell_weights) .or. &
            size(cell_ledger, 2) /= size(source_rate, 2) .or. &
            size(global_ledger) /= size(source_rate, 2)) return
        if (any(.not. ieee_is_finite(cell_weights)) .or. &
            any(.not. ieee_is_finite(source_rate)) .or. &
            any(.not. ieee_is_finite(cell_ledger)) .or. &
            any(.not. ieee_is_finite(global_ledger))) return
        if (any(cell_weights <= 0.0_dp)) return
        valid_value_inputs = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function valid_value_inputs

end module fortfem_volume_balance_ledger
