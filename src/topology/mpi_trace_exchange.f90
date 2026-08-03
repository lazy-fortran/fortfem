module fortfem_mpi_trace_exchange
    !! Communicator-free owner/ghost transport schedule for trace rows.
    !!
    !! The schedule records packed local rows, their global IDs, owner ranks,
    !! and deterministic owner/ghost buffer indices.  A future MPI backend can
    !! exchange the two packed buffers without changing any finite-element
    !! kernel.  Packing and unpacking are fixed-topology linear maps with
    !! exact JVP and VJP actions; no MPI dependency is introduced here.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_partition_layout, only: &
        partition_layout_maps, partition_layout_t, &
        partition_layout_global_count, validate_partition_layout
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    type, public :: mpi_trace_exchange_schedule_t
        private
        integer :: partition_count = 0
        integer :: global_count = 0
        integer :: component_count = 0
        integer :: local_count = 0
        integer :: owner_count = 0
        integer :: ghost_count = 0
        integer, allocatable :: row_offsets(:)
        integer, allocatable :: local_to_global(:), owner_rank(:)
        logical, allocatable :: owned(:)
        integer, allocatable :: owner_buffer_index(:), ghost_buffer_index(:)
        integer, allocatable :: global_owner_index(:)
    contains
        procedure, private :: assign_mpi_trace_exchange_schedule
        generic :: assignment(=) => assign_mpi_trace_exchange_schedule
    end type mpi_trace_exchange_schedule_t

    public :: initialize_mpi_trace_exchange_schedule
    public :: validate_mpi_trace_exchange_schedule
    public :: mpi_trace_exchange_schedule_maps
    public :: pack_mpi_trace_exchange
    public :: pack_mpi_trace_exchange_jvp
    public :: pack_mpi_trace_exchange_vjp
    public :: unpack_mpi_trace_exchange
    public :: unpack_mpi_trace_exchange_jvp
    public :: unpack_mpi_trace_exchange_vjp
    public :: pack_complex_mpi_trace_exchange
    public :: pack_complex_mpi_trace_exchange_jvp
    public :: pack_complex_mpi_trace_exchange_vjp
    public :: unpack_complex_mpi_trace_exchange
    public :: unpack_complex_mpi_trace_exchange_jvp
    public :: unpack_complex_mpi_trace_exchange_vjp

contains

    subroutine assign_mpi_trace_exchange_schedule(lhs, rhs)
        class(mpi_trace_exchange_schedule_t), intent(out) :: lhs
        type(mpi_trace_exchange_schedule_t), intent(in) :: rhs

        lhs%partition_count = rhs%partition_count
        lhs%global_count = rhs%global_count
        lhs%component_count = rhs%component_count
        lhs%local_count = rhs%local_count
        lhs%owner_count = rhs%owner_count
        lhs%ghost_count = rhs%ghost_count
        if (allocated(rhs%row_offsets)) allocate(lhs%row_offsets, source=rhs%row_offsets)
        if (allocated(rhs%local_to_global)) allocate(lhs%local_to_global, source=rhs%local_to_global)
        if (allocated(rhs%owner_rank)) allocate(lhs%owner_rank, source=rhs%owner_rank)
        if (allocated(rhs%owned)) allocate(lhs%owned, source=rhs%owned)
        if (allocated(rhs%owner_buffer_index)) allocate(lhs%owner_buffer_index, source=rhs%owner_buffer_index)
        if (allocated(rhs%ghost_buffer_index)) allocate(lhs%ghost_buffer_index, source=rhs%ghost_buffer_index)
        if (allocated(rhs%global_owner_index)) allocate(lhs%global_owner_index, source=rhs%global_owner_index)
    end subroutine assign_mpi_trace_exchange_schedule

    subroutine initialize_mpi_trace_exchange_schedule(schedule, partitions, component_count, status)
        type(mpi_trace_exchange_schedule_t), intent(inout) :: schedule
        type(partition_layout_t), intent(in) :: partitions(:)
        integer, intent(in) :: component_count
        type(fortsparse_status_t), intent(out) :: status
        integer :: partition, local, row, global_id, rank, global_count
        integer :: local_count, owner_count, ghost_count
        integer, allocatable :: local_to_global(:), owner_rank(:)
        logical, allocatable :: owned(:), seen(:), owner_seen(:), rank_seen(:)
        integer, allocatable :: expected_rank(:)

        call clear_schedule(schedule)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "MPI trace schedule received incompatible partitions")
        if (size(partitions) < 1 .or. component_count < 1) return
        global_count = partition_layout_global_count(partitions(1))
        if (global_count < 1) return
        do partition = 1, size(partitions)
            call validate_partition_layout(partitions(partition), status)
            if (status%code /= FORTSPARSE_OK .or. &
                partition_layout_global_count(partitions(partition)) /= global_count) return
        end do

        schedule%partition_count = size(partitions)
        schedule%global_count = global_count
        schedule%component_count = component_count
        allocate(schedule%row_offsets(schedule%partition_count + 1))
        schedule%row_offsets(1) = 1
        do partition = 1, schedule%partition_count
            call partition_layout_maps(partitions(partition), local_to_global, owner_rank, &
                owned, rank, status)
            if (status%code /= FORTSPARSE_OK) return
            schedule%row_offsets(partition + 1) = schedule%row_offsets(partition) + size(local_to_global)
        end do
        schedule%local_count = schedule%row_offsets(schedule%partition_count + 1) - 1
        allocate(schedule%local_to_global(schedule%local_count), schedule%owner_rank(schedule%local_count), &
            schedule%owned(schedule%local_count), schedule%owner_buffer_index(schedule%local_count), &
            schedule%ghost_buffer_index(schedule%local_count), schedule%global_owner_index(global_count), &
            seen(global_count), owner_seen(global_count), rank_seen(global_count), &
            expected_rank(global_count))
        schedule%owner_buffer_index = 0
        schedule%ghost_buffer_index = 0
        schedule%global_owner_index = 0
        seen = .false.
        owner_seen = .false.
        rank_seen = .false.
        expected_rank = -1
        owner_count = 0
        ghost_count = 0
        do partition = 1, schedule%partition_count
            call partition_layout_maps(partitions(partition), local_to_global, owner_rank, &
                owned, rank, status)
            if (status%code /= FORTSPARSE_OK) return
            do local = 1, size(local_to_global)
                row = schedule%row_offsets(partition) + local - 1
                global_id = local_to_global(local)
                schedule%local_to_global(row) = global_id
                schedule%owner_rank(row) = owner_rank(local)
                schedule%owned(row) = owned(local)
                if (rank_seen(global_id)) then
                    if (expected_rank(global_id) /= owner_rank(local)) then
                        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                            "MPI trace schedule found inconsistent owner ranks")
                        return
                    end if
                else
                    rank_seen(global_id) = .true.
                    expected_rank(global_id) = owner_rank(local)
                end if
                seen(global_id) = .true.
                if (owner_seen(global_id) .and. owned(local)) then
                    call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        "MPI trace schedule found duplicate owned rows")
                    return
                end if
                if (owner_rank(local) < 0) return
                if (owned(local)) then
                    owner_count = owner_count + 1
                    owner_seen(global_id) = .true.
                    schedule%global_owner_index(global_id) = owner_count
                    schedule%owner_buffer_index(row) = owner_count
                else
                    ghost_count = ghost_count + 1
                    schedule%ghost_buffer_index(row) = ghost_count
                end if
            end do
        end do
        if (any(.not. seen) .or. any(.not. owner_seen)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "MPI trace schedule requires one owner for every global row")
            return
        end if
        schedule%owner_count = owner_count
        schedule%ghost_count = ghost_count
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine initialize_mpi_trace_exchange_schedule

    subroutine validate_mpi_trace_exchange_schedule(schedule, status)
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        type(fortsparse_status_t), intent(out) :: status
        integer :: row

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "MPI trace schedule has invalid metadata")
        if (schedule%partition_count < 1 .or. schedule%global_count < 1 .or. &
            schedule%component_count < 1 .or. schedule%local_count < 1 .or. &
            schedule%owner_count < 1 .or. schedule%ghost_count < 0) return
        if (.not. allocated(schedule%row_offsets) .or. &
            .not. allocated(schedule%local_to_global) .or. .not. allocated(schedule%owner_rank) .or. &
            .not. allocated(schedule%owned) .or. .not. allocated(schedule%owner_buffer_index) .or. &
            .not. allocated(schedule%ghost_buffer_index) .or. .not. allocated(schedule%global_owner_index)) return
        if (size(schedule%row_offsets) /= schedule%partition_count + 1 .or. &
            schedule%row_offsets(1) /= 1 .or. schedule%row_offsets(schedule%partition_count + 1) /= schedule%local_count + 1 .or. &
            any(schedule%row_offsets(2:) <= schedule%row_offsets(:schedule%partition_count)) .or. &
            size(schedule%local_to_global) /= schedule%local_count .or. size(schedule%owner_rank) /= schedule%local_count .or. &
            size(schedule%owned) /= schedule%local_count .or. size(schedule%owner_buffer_index) /= schedule%local_count .or. &
            size(schedule%ghost_buffer_index) /= schedule%local_count .or. size(schedule%global_owner_index) /= schedule%global_count) return
        if (any(schedule%local_to_global < 1) .or. any(schedule%local_to_global > schedule%global_count) .or. &
            any(schedule%owner_rank < 0) .or. any(schedule%global_owner_index < 1)) return
        if (count(schedule%owned) /= schedule%owner_count .or. count(.not. schedule%owned) /= schedule%ghost_count) return
        do row = 1, schedule%local_count
            if (schedule%owned(row)) then
                if (schedule%owner_buffer_index(row) < 1 .or. schedule%owner_buffer_index(row) > schedule%owner_count .or. &
                    schedule%ghost_buffer_index(row) /= 0) return
            else
                if (schedule%ghost_buffer_index(row) < 1 .or. schedule%ghost_buffer_index(row) > schedule%ghost_count .or. &
                    schedule%owner_buffer_index(row) /= 0) return
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_mpi_trace_exchange_schedule

    subroutine mpi_trace_exchange_schedule_maps(schedule, row_offsets, local_to_global, owner_rank, owned, &
            owner_buffer_index, ghost_buffer_index, status)
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, allocatable, intent(out) :: row_offsets(:), local_to_global(:), owner_rank(:)
        logical, allocatable, intent(out) :: owned(:)
        integer, allocatable, intent(out) :: owner_buffer_index(:), ghost_buffer_index(:)
        type(fortsparse_status_t), intent(out) :: status

        if (allocated(row_offsets)) deallocate(row_offsets)
        if (allocated(local_to_global)) deallocate(local_to_global)
        if (allocated(owner_rank)) deallocate(owner_rank)
        if (allocated(owned)) deallocate(owned)
        if (allocated(owner_buffer_index)) deallocate(owner_buffer_index)
        if (allocated(ghost_buffer_index)) deallocate(ghost_buffer_index)
        call validate_mpi_trace_exchange_schedule(schedule, status)
        if (status%code /= FORTSPARSE_OK) return
        allocate(row_offsets, source=schedule%row_offsets)
        allocate(local_to_global, source=schedule%local_to_global)
        allocate(owner_rank, source=schedule%owner_rank)
        allocate(owned, source=schedule%owned)
        allocate(owner_buffer_index, source=schedule%owner_buffer_index)
        allocate(ghost_buffer_index, source=schedule%ghost_buffer_index)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine mpi_trace_exchange_schedule_maps

    subroutine pack_mpi_trace_exchange(schedule, offsets, local_values, owner_values, ghost_values, status)
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, intent(in) :: offsets(:)
        real(dp), intent(in) :: local_values(:, :)
        real(dp), intent(out) :: owner_values(:, :), ghost_values(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: row

        owner_values = 0.0_dp
        ghost_values = 0.0_dp
        if (.not. valid_pack_inputs(schedule, offsets, local_values, owner_values, ghost_values, status)) return
        do row = 1, schedule%local_count
            if (schedule%owned(row)) then
                owner_values(schedule%owner_buffer_index(row), :) = local_values(row, :)
            else
                ghost_values(schedule%ghost_buffer_index(row), :) = local_values(row, :)
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine pack_mpi_trace_exchange

    subroutine pack_mpi_trace_exchange_jvp(schedule, offsets, local_values_dot, owner_values_dot, ghost_values_dot, status)
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, intent(in) :: offsets(:)
        real(dp), intent(in) :: local_values_dot(:, :)
        real(dp), intent(out) :: owner_values_dot(:, :), ghost_values_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call pack_mpi_trace_exchange(schedule, offsets, local_values_dot, owner_values_dot, ghost_values_dot, status)
    end subroutine pack_mpi_trace_exchange_jvp

    subroutine pack_mpi_trace_exchange_vjp(schedule, offsets, owner_values_bar, ghost_values_bar, local_values_bar, status)
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, intent(in) :: offsets(:)
        real(dp), intent(in) :: owner_values_bar(:, :), ghost_values_bar(:, :)
        real(dp), intent(out) :: local_values_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: row

        local_values_bar = 0.0_dp
        if (.not. valid_pack_bar_inputs(schedule, offsets, owner_values_bar, ghost_values_bar, local_values_bar, status)) return
        do row = 1, schedule%local_count
            if (schedule%owned(row)) then
                local_values_bar(row, :) = owner_values_bar(schedule%owner_buffer_index(row), :)
            else
                local_values_bar(row, :) = ghost_values_bar(schedule%ghost_buffer_index(row), :)
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine pack_mpi_trace_exchange_vjp

    subroutine unpack_mpi_trace_exchange(schedule, offsets, owner_values, ghost_values, local_values, status)
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, intent(in) :: offsets(:)
        real(dp), intent(in) :: owner_values(:, :), ghost_values(:, :)
        real(dp), intent(out) :: local_values(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: row

        local_values = 0.0_dp
        if (.not. valid_unpack_inputs(schedule, offsets, owner_values, ghost_values, local_values, status)) return
        do row = 1, schedule%local_count
            if (schedule%owned(row)) then
                local_values(row, :) = owner_values(schedule%owner_buffer_index(row), :)
            else
                local_values(row, :) = ghost_values(schedule%ghost_buffer_index(row), :)
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine unpack_mpi_trace_exchange

    subroutine unpack_mpi_trace_exchange_jvp(schedule, offsets, owner_values_dot, ghost_values_dot, local_values_dot, status)
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, intent(in) :: offsets(:)
        real(dp), intent(in) :: owner_values_dot(:, :), ghost_values_dot(:, :)
        real(dp), intent(out) :: local_values_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call unpack_mpi_trace_exchange(schedule, offsets, owner_values_dot, ghost_values_dot, local_values_dot, status)
    end subroutine unpack_mpi_trace_exchange_jvp

    subroutine unpack_mpi_trace_exchange_vjp(schedule, offsets, local_values_bar, owner_values_bar, ghost_values_bar, status)
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, intent(in) :: offsets(:)
        real(dp), intent(in) :: local_values_bar(:, :)
        real(dp), intent(out) :: owner_values_bar(:, :), ghost_values_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: row

        owner_values_bar = 0.0_dp
        ghost_values_bar = 0.0_dp
        if (.not. valid_unpack_bar_inputs(schedule, offsets, local_values_bar, owner_values_bar, ghost_values_bar, status)) return
        do row = 1, schedule%local_count
            if (schedule%owned(row)) then
                owner_values_bar(schedule%owner_buffer_index(row), :) = local_values_bar(row, :)
            else
                ghost_values_bar(schedule%ghost_buffer_index(row), :) = local_values_bar(row, :)
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine unpack_mpi_trace_exchange_vjp

    subroutine pack_complex_mpi_trace_exchange( &
            schedule, offsets, local_values, owner_values, ghost_values, status)
        !! Pack complex trace rows without changing the real transport API.
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, intent(in) :: offsets(:)
        complex(dp), intent(in) :: local_values(:, :)
        complex(dp), intent(out) :: owner_values(:, :), ghost_values(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: row

        owner_values = cmplx(0.0_dp, 0.0_dp, dp)
        ghost_values = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. valid_complex_exchange_shapes( &
                schedule, offsets, local_values, owner_values, ghost_values, &
                status)) return
        if (.not. finite_complex_2d(local_values)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "complex MPI trace packing received non-finite values")
            return
        end if
        do row = 1, schedule%local_count
            if (schedule%owned(row)) then
                owner_values(schedule%owner_buffer_index(row), :) = local_values(row, :)
            else
                ghost_values(schedule%ghost_buffer_index(row), :) = local_values(row, :)
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine pack_complex_mpi_trace_exchange

    subroutine pack_complex_mpi_trace_exchange_jvp( &
            schedule, offsets, local_values_dot, owner_values_dot, &
            ghost_values_dot, status)
        !! Apply the exact fixed-schedule complex packing JVP.
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, intent(in) :: offsets(:)
        complex(dp), intent(in) :: local_values_dot(:, :)
        complex(dp), intent(out) :: owner_values_dot(:, :), ghost_values_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call pack_complex_mpi_trace_exchange( &
            schedule, offsets, local_values_dot, owner_values_dot, &
            ghost_values_dot, status)
    end subroutine pack_complex_mpi_trace_exchange_jvp

    subroutine pack_complex_mpi_trace_exchange_vjp( &
            schedule, offsets, owner_values_bar, ghost_values_bar, &
            local_values_bar, status)
        !! Apply the packing VJP under the real-part complex pairing.
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, intent(in) :: offsets(:)
        complex(dp), intent(in) :: owner_values_bar(:, :), ghost_values_bar(:, :)
        complex(dp), intent(out) :: local_values_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: row

        local_values_bar = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. valid_complex_exchange_shapes( &
                schedule, offsets, local_values_bar, owner_values_bar, &
                ghost_values_bar, status)) return
        if (.not. finite_complex_2d(owner_values_bar) .or. &
            .not. finite_complex_2d(ghost_values_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "complex MPI trace packing VJP received non-finite values")
            return
        end if
        do row = 1, schedule%local_count
            if (schedule%owned(row)) then
                local_values_bar(row, :) = &
                    owner_values_bar(schedule%owner_buffer_index(row), :)
            else
                local_values_bar(row, :) = &
                    ghost_values_bar(schedule%ghost_buffer_index(row), :)
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine pack_complex_mpi_trace_exchange_vjp

    subroutine unpack_complex_mpi_trace_exchange( &
            schedule, offsets, owner_values, ghost_values, local_values, status)
        !! Unpack complex owner and received ghost buffers into local rows.
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, intent(in) :: offsets(:)
        complex(dp), intent(in) :: owner_values(:, :), ghost_values(:, :)
        complex(dp), intent(out) :: local_values(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: row

        local_values = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. valid_complex_exchange_shapes( &
                schedule, offsets, local_values, owner_values, ghost_values, &
                status)) return
        if (.not. finite_complex_2d(owner_values) .or. &
            .not. finite_complex_2d(ghost_values)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "complex MPI trace unpacking received non-finite values")
            return
        end if
        do row = 1, schedule%local_count
            if (schedule%owned(row)) then
                local_values(row, :) = owner_values(schedule%owner_buffer_index(row), :)
            else
                local_values(row, :) = ghost_values(schedule%ghost_buffer_index(row), :)
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine unpack_complex_mpi_trace_exchange

    subroutine unpack_complex_mpi_trace_exchange_jvp( &
            schedule, offsets, owner_values_dot, ghost_values_dot, &
            local_values_dot, status)
        !! Apply the exact fixed-schedule complex unpacking JVP.
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, intent(in) :: offsets(:)
        complex(dp), intent(in) :: owner_values_dot(:, :), ghost_values_dot(:, :)
        complex(dp), intent(out) :: local_values_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call unpack_complex_mpi_trace_exchange( &
            schedule, offsets, owner_values_dot, ghost_values_dot, &
            local_values_dot, status)
    end subroutine unpack_complex_mpi_trace_exchange_jvp

    subroutine unpack_complex_mpi_trace_exchange_vjp( &
            schedule, offsets, local_values_bar, owner_values_bar, &
            ghost_values_bar, status)
        !! Apply the unpacking VJP under the real-part complex pairing.
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, intent(in) :: offsets(:)
        complex(dp), intent(in) :: local_values_bar(:, :)
        complex(dp), intent(out) :: owner_values_bar(:, :), ghost_values_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: row

        owner_values_bar = cmplx(0.0_dp, 0.0_dp, dp)
        ghost_values_bar = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. valid_complex_exchange_shapes( &
                schedule, offsets, local_values_bar, owner_values_bar, &
                ghost_values_bar, status)) return
        if (.not. finite_complex_2d(local_values_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "complex MPI trace unpacking VJP received non-finite values")
            return
        end if
        do row = 1, schedule%local_count
            if (schedule%owned(row)) then
                owner_values_bar(schedule%owner_buffer_index(row), :) = &
                    local_values_bar(row, :)
            else
                ghost_values_bar(schedule%ghost_buffer_index(row), :) = &
                    local_values_bar(row, :)
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine unpack_complex_mpi_trace_exchange_vjp

    logical function valid_pack_inputs(schedule, offsets, local_values, owner_values, ghost_values, status) result(valid)
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, intent(in) :: offsets(:)
        real(dp), intent(in) :: local_values(:, :), owner_values(:, :), ghost_values(:, :)
        type(fortsparse_status_t), intent(out) :: status

        valid = .false.
        call validate_mpi_trace_exchange_schedule(schedule, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. offsets_match(schedule, offsets) .or. size(local_values, 1) /= schedule%local_count .or. &
            size(local_values, 2) /= schedule%component_count .or. size(owner_values, 1) /= schedule%owner_count .or. &
            size(owner_values, 2) /= schedule%component_count .or. size(ghost_values, 1) /= schedule%ghost_count .or. &
            size(ghost_values, 2) /= schedule%component_count .or. any(.not. ieee_is_finite(local_values))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, "MPI trace packing received incompatible values")
            return
        end if
        valid = .true.
    end function valid_pack_inputs

    logical function valid_pack_bar_inputs(schedule, offsets, owner_values_bar, ghost_values_bar, local_values_bar, status) result(valid)
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, intent(in) :: offsets(:)
        real(dp), intent(in) :: owner_values_bar(:, :), ghost_values_bar(:, :), local_values_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        valid = .false.
        call validate_mpi_trace_exchange_schedule(schedule, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. offsets_match(schedule, offsets) .or. size(owner_values_bar, 1) /= schedule%owner_count .or. &
            size(owner_values_bar, 2) /= schedule%component_count .or. size(ghost_values_bar, 1) /= schedule%ghost_count .or. &
            size(ghost_values_bar, 2) /= schedule%component_count .or. size(local_values_bar, 1) /= schedule%local_count .or. &
            size(local_values_bar, 2) /= schedule%component_count .or. any(.not. ieee_is_finite(owner_values_bar)) .or. &
            any(.not. ieee_is_finite(ghost_values_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, "MPI trace packing VJP received incompatible values")
            return
        end if
        valid = .true.
    end function valid_pack_bar_inputs

    logical function valid_unpack_inputs(schedule, offsets, owner_values, ghost_values, local_values, status) result(valid)
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, intent(in) :: offsets(:)
        real(dp), intent(in) :: owner_values(:, :), ghost_values(:, :), local_values(:, :)
        type(fortsparse_status_t), intent(out) :: status

        valid = .false.
        call validate_mpi_trace_exchange_schedule(schedule, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. offsets_match(schedule, offsets) .or. size(owner_values, 1) /= schedule%owner_count .or. &
            size(owner_values, 2) /= schedule%component_count .or. size(ghost_values, 1) /= schedule%ghost_count .or. &
            size(ghost_values, 2) /= schedule%component_count .or. size(local_values, 1) /= schedule%local_count .or. &
            size(local_values, 2) /= schedule%component_count .or. any(.not. ieee_is_finite(owner_values)) .or. &
            any(.not. ieee_is_finite(ghost_values))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, "MPI trace unpacking received incompatible values")
            return
        end if
        valid = .true.
    end function valid_unpack_inputs

    logical function valid_unpack_bar_inputs(schedule, offsets, local_values_bar, owner_values_bar, ghost_values_bar, status) result(valid)
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, intent(in) :: offsets(:)
        real(dp), intent(in) :: local_values_bar(:, :), owner_values_bar(:, :), ghost_values_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        valid = .false.
        call validate_mpi_trace_exchange_schedule(schedule, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. offsets_match(schedule, offsets) .or. size(local_values_bar, 1) /= schedule%local_count .or. &
            size(local_values_bar, 2) /= schedule%component_count .or. size(owner_values_bar, 1) /= schedule%owner_count .or. &
            size(owner_values_bar, 2) /= schedule%component_count .or. size(ghost_values_bar, 1) /= schedule%ghost_count .or. &
            size(ghost_values_bar, 2) /= schedule%component_count .or. any(.not. ieee_is_finite(local_values_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, "MPI trace unpacking VJP received incompatible values")
            return
        end if
        valid = .true.
    end function valid_unpack_bar_inputs

    logical function valid_complex_exchange_shapes( &
            schedule, offsets, local_values, owner_values, ghost_values, &
            status) result(valid)
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, intent(in) :: offsets(:)
        complex(dp), intent(in) :: local_values(:, :), owner_values(:, :)
        complex(dp), intent(in) :: ghost_values(:, :)
        type(fortsparse_status_t), intent(out) :: status

        valid = .false.
        call validate_mpi_trace_exchange_schedule(schedule, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. offsets_match(schedule, offsets) .or. &
            size(local_values, 1) /= schedule%local_count .or. &
            size(local_values, 2) /= schedule%component_count .or. &
            size(owner_values, 1) /= schedule%owner_count .or. &
            size(owner_values, 2) /= schedule%component_count .or. &
            size(ghost_values, 1) /= schedule%ghost_count .or. &
            size(ghost_values, 2) /= schedule%component_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "complex MPI trace exchange received incompatible shapes")
            return
        end if
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function valid_complex_exchange_shapes

    logical function finite_complex_2d(values) result(finite)
        complex(dp), intent(in) :: values(:, :)

        finite = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex_2d

    logical function offsets_match(schedule, offsets) result(matches)
        type(mpi_trace_exchange_schedule_t), intent(in) :: schedule
        integer, intent(in) :: offsets(:)

        matches = size(offsets) == size(schedule%row_offsets)
        if (matches) matches = all(offsets == schedule%row_offsets)
    end function offsets_match

    subroutine clear_schedule(schedule)
        type(mpi_trace_exchange_schedule_t), intent(inout) :: schedule

        schedule%partition_count = 0
        schedule%global_count = 0
        schedule%component_count = 0
        schedule%local_count = 0
        schedule%owner_count = 0
        schedule%ghost_count = 0
        if (allocated(schedule%row_offsets)) deallocate(schedule%row_offsets)
        if (allocated(schedule%local_to_global)) deallocate(schedule%local_to_global)
        if (allocated(schedule%owner_rank)) deallocate(schedule%owner_rank)
        if (allocated(schedule%owned)) deallocate(schedule%owned)
        if (allocated(schedule%owner_buffer_index)) deallocate(schedule%owner_buffer_index)
        if (allocated(schedule%ghost_buffer_index)) deallocate(schedule%ghost_buffer_index)
        if (allocated(schedule%global_owner_index)) deallocate(schedule%global_owner_index)
    end subroutine clear_schedule

end module fortfem_mpi_trace_exchange
