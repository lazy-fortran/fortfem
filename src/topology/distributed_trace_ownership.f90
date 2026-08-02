module fortfem_distributed_trace_ownership
    !! Serial ownership and reduction ledger for distributed trace rows.
    !!
    !! A ledger packs several communicator-free partition_layout_t instances.
    !! Each instance contributes local physical trace rows, while its layout
    !! records the global row ID and the owner/ghost decision.  The ledger is
    !! intentionally MPI-free: a later communicator backend can exchange the
    !! same packed rows without changing a residual or mortar kernel.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_partition_layout, only: &
        partition_layout_t, &
        partition_layout_maps, partition_layout_global_count, &
        validate_partition_layout
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    type, public :: distributed_trace_layout_t
        private
        integer :: partition_count = 0
        integer :: global_count = 0
        integer :: local_count = 0
        integer :: component_count = 0
        type(partition_layout_t), allocatable :: partitions(:)
        integer, allocatable :: row_offsets(:)
    end type distributed_trace_layout_t

    public :: initialize_distributed_trace_layout
    public :: validate_distributed_trace_layout
    public :: distributed_trace_layout_partition_count
    public :: distributed_trace_layout_global_count
    public :: distributed_trace_layout_local_count
    public :: distributed_trace_layout_component_count
    public :: assemble_distributed_trace_reduction
    public :: assemble_distributed_trace_reduction_jvp
    public :: assemble_distributed_trace_reduction_vjp

contains

    subroutine initialize_distributed_trace_layout( &
            layout, partitions, component_count, status)
        !! Build a fixed packed ledger from complete partition metadata.
        type(distributed_trace_layout_t), intent(inout) :: layout
        type(partition_layout_t), intent(in) :: partitions(:)
        integer, intent(in) :: component_count
        type(fortsparse_status_t), intent(out) :: status

        integer :: partition, local_count, global_count
        integer, allocatable :: local_to_global(:), owner_rank(:)
        logical, allocatable :: owned(:)
        integer :: rank

        call clear_distributed_trace_layout(layout)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "distributed trace layout received incompatible metadata")
        if (size(partitions) < 1 .or. component_count < 1) return
        global_count = partition_layout_global_count(partitions(1))
        if (global_count < 1) return
        do partition = 1, size(partitions)
            call validate_partition_layout(partitions(partition), status)
            if (status%code /= FORTSPARSE_OK) return
            if (partition_layout_global_count(partitions(partition)) /= &
                global_count) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "distributed partitions disagree on global row count")
                return
            end if
        end do

        layout%partition_count = size(partitions)
        layout%global_count = global_count
        layout%component_count = component_count
        allocate(layout%partitions(layout%partition_count))
        allocate(layout%row_offsets(layout%partition_count + 1))
        layout%partitions = partitions
        layout%row_offsets(1) = 1
        do partition = 1, layout%partition_count
            call partition_layout_maps( &
                layout%partitions(partition), local_to_global, owner_rank, &
                owned, rank, status)
            if (status%code /= FORTSPARSE_OK) return
            local_count = size(local_to_global)
            layout%row_offsets(partition + 1) = &
                layout%row_offsets(partition) + local_count
        end do
        layout%local_count = layout%row_offsets(layout%partition_count + 1) - 1
        call validate_distributed_trace_layout(layout, status)
    end subroutine initialize_distributed_trace_layout

    subroutine validate_distributed_trace_layout(layout, status)
        !! Validate IDs, owner ranks, masks, and the complete owner ledger.
        type(distributed_trace_layout_t), intent(in) :: layout
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: local_to_global(:), owner_rank(:)
        integer, allocatable :: expected_owner(:), owner_count(:)
        logical, allocatable :: owned(:), seen(:)
        integer :: partition, local, global_id, rank, local_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "distributed trace layout received incompatible metadata")
        if (layout%partition_count < 1 .or. layout%global_count < 1 .or. &
            layout%local_count < 1 .or. layout%component_count < 1) return
        if (.not. allocated(layout%partitions) .or. &
            .not. allocated(layout%row_offsets)) return
        if (size(layout%partitions) /= layout%partition_count .or. &
            size(layout%row_offsets) /= layout%partition_count + 1) return
        if (layout%row_offsets(1) /= 1 .or. &
            layout%row_offsets(layout%partition_count + 1) /= &
            layout%local_count + 1) return
        if (any(layout%row_offsets(2:) <= &
            layout%row_offsets(:layout%partition_count))) return

        allocate(expected_owner(layout%global_count), owner_count(layout%global_count))
        allocate(seen(layout%global_count))
        expected_owner = -1
        owner_count = 0
        seen = .false.
        do partition = 1, layout%partition_count
            call validate_partition_layout(layout%partitions(partition), status)
            if (status%code /= FORTSPARSE_OK) return
            if (partition_layout_global_count(layout%partitions(partition)) /= &
                layout%global_count) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "distributed partitions disagree on global row count")
                return
            end if
            call partition_layout_maps( &
                layout%partitions(partition), local_to_global, owner_rank, &
                owned, rank, status)
            if (status%code /= FORTSPARSE_OK) return
            local_count = size(local_to_global)
            if (layout%row_offsets(partition + 1) - &
                layout%row_offsets(partition) /= local_count) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "distributed trace row offsets are inconsistent")
                return
            end if
            do local = 1, local_count
                global_id = local_to_global(local)
                if (.not. seen(global_id)) then
                    seen(global_id) = .true.
                    expected_owner(global_id) = owner_rank(local)
                else if (expected_owner(global_id) /= owner_rank(local)) then
                    call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        "distributed trace rows disagree on an owner rank")
                    return
                end if
                if (owned(local)) owner_count(global_id) = &
                    owner_count(global_id) + 1
            end do
        end do
        do global_id = 1, layout%global_count
            if (seen(global_id) .and. owner_count(global_id) /= 1) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "distributed trace row does not have exactly one owner")
                return
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_distributed_trace_layout

    integer function distributed_trace_layout_partition_count(layout)
        type(distributed_trace_layout_t), intent(in) :: layout

        distributed_trace_layout_partition_count = layout%partition_count
    end function distributed_trace_layout_partition_count

    integer function distributed_trace_layout_global_count(layout)
        type(distributed_trace_layout_t), intent(in) :: layout

        distributed_trace_layout_global_count = layout%global_count
    end function distributed_trace_layout_global_count

    integer function distributed_trace_layout_local_count(layout)
        type(distributed_trace_layout_t), intent(in) :: layout

        distributed_trace_layout_local_count = layout%local_count
    end function distributed_trace_layout_local_count

    integer function distributed_trace_layout_component_count(layout)
        type(distributed_trace_layout_t), intent(in) :: layout

        distributed_trace_layout_component_count = layout%component_count
    end function distributed_trace_layout_component_count

    subroutine assemble_distributed_trace_reduction( &
            layout, local_values, global_values, owner_values, status)
        !! Sum all local rows and, separately, the unique owned row.
        type(distributed_trace_layout_t), intent(in) :: layout
        real(dp), intent(in) :: local_values(:, :)
        real(dp), intent(out) :: global_values(:, :), owner_values(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: partition, local, row, global_id, rank
        integer, allocatable :: local_to_global(:), owner_rank(:)
        logical, allocatable :: owned(:)

        call validate_reduction_inputs( &
            layout, local_values, global_values, owner_values, status)
        if (status%code /= FORTSPARSE_OK) return
        global_values = 0.0_dp
        owner_values = 0.0_dp
        do partition = 1, layout%partition_count
            call partition_layout_maps( &
                layout%partitions(partition), local_to_global, owner_rank, &
                owned, rank, status)
            if (status%code /= FORTSPARSE_OK) return
            do local = 1, size(local_to_global)
                row = layout%row_offsets(partition) + local - 1
                global_id = local_to_global(local)
                global_values(global_id, :) = global_values(global_id, :) + &
                    local_values(row, :)
                if (owned(local)) owner_values(global_id, :) = &
                    owner_values(global_id, :) + local_values(row, :)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_distributed_trace_reduction

    subroutine assemble_distributed_trace_reduction_jvp( &
            layout, local_values_dot, global_values_dot, owner_values_dot, status)
        !! Apply the fixed-topology JVP of the trace reduction.
        type(distributed_trace_layout_t), intent(in) :: layout
        real(dp), intent(in) :: local_values_dot(:, :)
        real(dp), intent(out) :: global_values_dot(:, :), owner_values_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call validate_reduction_inputs( &
            layout, local_values_dot, global_values_dot, owner_values_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        call assemble_distributed_trace_reduction( &
            layout, local_values_dot, global_values_dot, owner_values_dot, status)
    end subroutine assemble_distributed_trace_reduction_jvp

    subroutine assemble_distributed_trace_reduction_vjp( &
            layout, global_values_bar, owner_values_bar, local_values_bar, status)
        !! Apply the real transpose of the two-output trace reduction.
        type(distributed_trace_layout_t), intent(in) :: layout
        real(dp), intent(in) :: global_values_bar(:, :), owner_values_bar(:, :)
        real(dp), intent(out) :: local_values_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: partition, local, row, global_id, rank
        integer, allocatable :: local_to_global(:), owner_rank(:)
        logical, allocatable :: owned(:)

        local_values_bar = 0.0_dp
        call validate_vjp_inputs( &
            layout, global_values_bar, owner_values_bar, local_values_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        do partition = 1, layout%partition_count
            call partition_layout_maps( &
                layout%partitions(partition), local_to_global, owner_rank, &
                owned, rank, status)
            if (status%code /= FORTSPARSE_OK) return
            do local = 1, size(local_to_global)
                row = layout%row_offsets(partition) + local - 1
                global_id = local_to_global(local)
                local_values_bar(row, :) = global_values_bar(global_id, :)
                if (owned(local)) local_values_bar(row, :) = &
                    local_values_bar(row, :) + owner_values_bar(global_id, :)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_distributed_trace_reduction_vjp

    subroutine validate_reduction_inputs( &
            layout, local_values, global_values, owner_values, status)
        type(distributed_trace_layout_t), intent(in) :: layout
        real(dp), intent(in) :: local_values(:, :)
        real(dp), intent(in) :: global_values(:, :), owner_values(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call validate_distributed_trace_layout(layout, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(local_values, 1) /= layout%local_count .or. &
            size(local_values, 2) /= layout%component_count .or. &
            size(global_values, 1) /= layout%global_count .or. &
            size(global_values, 2) /= layout%component_count .or. &
            size(owner_values, 1) /= layout%global_count .or. &
            size(owner_values, 2) /= layout%component_count .or. &
            any(.not. ieee_is_finite(local_values))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "distributed trace reduction received incompatible arrays")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_reduction_inputs

    subroutine validate_vjp_inputs( &
            layout, global_values_bar, owner_values_bar, local_values_bar, status)
        type(distributed_trace_layout_t), intent(in) :: layout
        real(dp), intent(in) :: global_values_bar(:, :), owner_values_bar(:, :)
        real(dp), intent(in) :: local_values_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call validate_distributed_trace_layout(layout, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(global_values_bar, 1) /= layout%global_count .or. &
            size(global_values_bar, 2) /= layout%component_count .or. &
            size(owner_values_bar, 1) /= layout%global_count .or. &
            size(owner_values_bar, 2) /= layout%component_count .or. &
            size(local_values_bar, 1) /= layout%local_count .or. &
            size(local_values_bar, 2) /= layout%component_count .or. &
            any(.not. ieee_is_finite(global_values_bar)) .or. &
            any(.not. ieee_is_finite(owner_values_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "distributed trace VJP received incompatible arrays")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_vjp_inputs

    subroutine clear_distributed_trace_layout(layout)
        type(distributed_trace_layout_t), intent(inout) :: layout

        layout%partition_count = 0
        layout%global_count = 0
        layout%local_count = 0
        layout%component_count = 0
        if (allocated(layout%partitions)) deallocate(layout%partitions)
        if (allocated(layout%row_offsets)) deallocate(layout%row_offsets)
    end subroutine clear_distributed_trace_layout

end module fortfem_distributed_trace_ownership
