module fortfem_equation_objective_callbacks
    !! Callback dispatch for neutral equation/objective/constraint blocks.
    !!
    !! The metadata registry remains the source of deterministic ordering;
    !! this layer only invokes caller-owned value, tangent, and adjoint
    !! callbacks for each named block and packs their rows.  It owns no
    !! constitutive law, profile, coordinate convention, or optimizer.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_equation_objective_registry, only: &
        equation_objective_block_t, equation_objective_registry_t, &
        pack_equation_objective_values, &
        validate_equation_objective_registry
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    abstract interface
        subroutine equation_objective_value_callback(block, state, values, status)
            import dp, equation_objective_block_t
            type(equation_objective_block_t), intent(in) :: block
            real(dp), intent(in) :: state(:)
            real(dp), intent(out) :: values(:)
            integer, intent(out) :: status
        end subroutine equation_objective_value_callback

        subroutine equation_objective_jvp_callback( &
                block, state, state_dot, values_dot, status)
            import dp, equation_objective_block_t
            type(equation_objective_block_t), intent(in) :: block
            real(dp), intent(in) :: state(:), state_dot(:)
            real(dp), intent(out) :: values_dot(:)
            integer, intent(out) :: status
        end subroutine equation_objective_jvp_callback

        subroutine equation_objective_vjp_callback( &
                block, state, values_bar, state_bar, status)
            import dp, equation_objective_block_t
            type(equation_objective_block_t), intent(in) :: block
            real(dp), intent(in) :: state(:), values_bar(:)
            real(dp), intent(out) :: state_bar(:)
            integer, intent(out) :: status
        end subroutine equation_objective_vjp_callback
    end interface

    public :: evaluate_equation_objective_callbacks
    public :: evaluate_equation_objective_callbacks_jvp
    public :: evaluate_equation_objective_callbacks_vjp

contains

    subroutine evaluate_equation_objective_callbacks( &
            registry, state, packed, value_callback, status)
        type(equation_objective_registry_t), intent(in) :: registry
        real(dp), intent(in) :: state(:)
        real(dp), allocatable, intent(out) :: packed(:)
        procedure(equation_objective_value_callback) :: value_callback
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: block_values(:, :), values(:)
        integer :: block, max_rows, callback_status, row_count

        allocate(packed(0))
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "equation/objective callback evaluation received invalid inputs")
        if (.not. valid_registry_state(registry, state, status)) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "equation/objective value callback failed")
        max_rows = maxval(registry%blocks%row_count)
        allocate(block_values(max_rows, size(registry%blocks)))
        block_values = 0.0_dp
        do block = 1, size(registry%blocks)
            row_count = registry%blocks(block)%row_count
            if (row_count == 0) cycle
            allocate(values(row_count))
            values = 0.0_dp
            call value_callback(registry%blocks(block), state, values, callback_status)
            if (callback_status /= 0 .or. .not. finite_real(values)) then
                deallocate(values)
                return
            end if
            block_values(:row_count, block) = values
            deallocate(values)
        end do
        call pack_equation_objective_values(registry, block_values, packed, status)
    end subroutine evaluate_equation_objective_callbacks

    subroutine evaluate_equation_objective_callbacks_jvp( &
            registry, state, state_dot, packed_dot, jvp_callback, status)
        type(equation_objective_registry_t), intent(in) :: registry
        real(dp), intent(in) :: state(:), state_dot(:)
        real(dp), allocatable, intent(out) :: packed_dot(:)
        procedure(equation_objective_jvp_callback) :: jvp_callback
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: block_values_dot(:, :), values_dot(:)
        integer :: block, max_rows, callback_status, row_count

        allocate(packed_dot(0))
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "equation/objective callback JVP received invalid inputs")
        if (.not. valid_registry_state(registry, state, status)) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "equation/objective JVP callback failed")
        if (size(state_dot) /= size(state) .or. .not. finite_real(state_dot)) return
        max_rows = maxval(registry%blocks%row_count)
        allocate(block_values_dot(max_rows, size(registry%blocks)))
        block_values_dot = 0.0_dp
        do block = 1, size(registry%blocks)
            row_count = registry%blocks(block)%row_count
            if (row_count == 0) cycle
            allocate(values_dot(row_count))
            values_dot = 0.0_dp
            call jvp_callback( &
                registry%blocks(block), state, state_dot, values_dot, callback_status)
            if (callback_status /= 0 .or. .not. finite_real(values_dot)) then
                deallocate(values_dot)
                return
            end if
            block_values_dot(:row_count, block) = values_dot
            deallocate(values_dot)
        end do
        call pack_equation_objective_values( &
            registry, block_values_dot, packed_dot, status)
    end subroutine evaluate_equation_objective_callbacks_jvp

    subroutine evaluate_equation_objective_callbacks_vjp( &
            registry, state, packed_bar, state_bar, vjp_callback, status)
        type(equation_objective_registry_t), intent(in) :: registry
        real(dp), intent(in) :: state(:), packed_bar(:)
        real(dp), allocatable, intent(out) :: state_bar(:)
        procedure(equation_objective_vjp_callback) :: vjp_callback
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: block_bar(:), state_bar_local(:)
        integer :: block, callback_status, row_count

        allocate(state_bar(size(state)))
        state_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "equation/objective callback VJP received invalid inputs")
        if (.not. valid_registry_state(registry, state, status)) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "equation/objective VJP callback failed")
        if (size(packed_bar) /= registry%total_rows .or. &
            .not. finite_real(packed_bar)) return
        do block = 1, size(registry%blocks)
            row_count = registry%blocks(block)%row_count
            if (row_count == 0) cycle
            allocate(block_bar(row_count), state_bar_local(size(state)))
            block_bar = packed_bar( &
                registry%blocks(block)%row_offset:&
                registry%blocks(block)%row_offset + row_count - 1)
            state_bar_local = 0.0_dp
            call vjp_callback( &
                registry%blocks(block), state, block_bar, state_bar_local, &
                callback_status)
            if (callback_status /= 0 .or. .not. finite_real(state_bar_local)) then
                deallocate(block_bar, state_bar_local)
                return
            end if
            state_bar = state_bar + state_bar_local
            deallocate(block_bar, state_bar_local)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_equation_objective_callbacks_vjp

    logical function valid_registry_state(registry, state, status) result(valid)
        type(equation_objective_registry_t), intent(in) :: registry
        real(dp), intent(in) :: state(:)
        type(fortsparse_status_t), intent(inout) :: status

        valid = .false.
        if (.not. validate_equation_objective_registry(registry, status)) return
        if (size(state) < 1 .or. .not. finite_real(state)) return
        valid = .true.
    end function valid_registry_state

    pure logical function finite_real(values) result(finite)
        real(dp), intent(in) :: values(:)

        finite = all(ieee_is_finite(values))
    end function finite_real

end module fortfem_equation_objective_callbacks
