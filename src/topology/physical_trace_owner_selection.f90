module fortfem_physical_trace_owner_selection
    !! Deterministic owner selection and gather for physical trace partitions.
    !!
    !! The contract is deliberately communicator-free.  It checks that every
    !! global trace ID has one owned copy, that duplicate physical coordinates
    !! agree, and records the unique partition/local row used by a later MPI
    !! transport layer.  Gather, JVP, and VJP are fixed-topology operations.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_physical_trace_ownership, only: &
        physical_trace_ownership_t, physical_trace_ownership_maps, &
        physical_trace_ownership_dimension
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    type, public :: physical_trace_owner_selection_t
        private
        integer :: partition_count = 0
        integer :: global_count = 0
        integer :: component_dimension = 0
        integer, allocatable :: owner_partition(:), owner_local(:), owner_rank(:)
    end type physical_trace_owner_selection_t

    public :: initialize_physical_trace_owner_selection
    public :: validate_physical_trace_owner_selection
    public :: physical_trace_owner_selection_maps
    public :: gather_physical_trace_values
    public :: gather_physical_trace_values_jvp
    public :: gather_physical_trace_values_vjp

contains

    subroutine initialize_physical_trace_owner_selection( &
            selection, ownerships, coordinate_tolerance, status)
        type(physical_trace_owner_selection_t), intent(inout) :: selection
        type(physical_trace_ownership_t), intent(in) :: ownerships(:)
        real(dp), intent(in) :: coordinate_tolerance
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: coordinates(:, :)
        integer, allocatable :: global_ids(:), owner_rank(:)
        logical, allocatable :: owned(:), seen(:)
        real(dp), allocatable :: canonical(:, :)
        integer :: partition, local, global_id, maximum_id, dimension, rank
        real(dp) :: distance

        call clear_selection(selection)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "physical trace owner selection received invalid metadata")
        if (size(ownerships) < 1 .or. .not. ieee_is_finite(coordinate_tolerance) .or. &
            coordinate_tolerance < 0.0_dp) return
        dimension = physical_trace_ownership_dimension(ownerships(1))
        if (dimension < 1) return
        maximum_id = 0
        do partition = 1, size(ownerships)
            call physical_trace_ownership_maps(ownerships(partition), coordinates, &
                global_ids, owner_rank, owned, rank, status)
            if (status%code /= FORTSPARSE_OK .or. size(coordinates, 1) /= dimension) return
            maximum_id = max(maximum_id, maxval(global_ids))
        end do
        if (maximum_id < 1) return
        allocate(selection%owner_partition(maximum_id), selection%owner_local(maximum_id), &
            selection%owner_rank(maximum_id), seen(maximum_id), canonical(dimension, maximum_id))
        selection%partition_count = size(ownerships)
        selection%global_count = maximum_id
        selection%component_dimension = dimension
        selection%owner_partition = 0
        selection%owner_local = 0
        selection%owner_rank = -1
        seen = .false.
        canonical = 0.0_dp
        do partition = 1, size(ownerships)
            call physical_trace_ownership_maps(ownerships(partition), coordinates, &
                global_ids, owner_rank, owned, rank, status)
            if (status%code /= FORTSPARSE_OK) return
            do local = 1, size(global_ids)
                global_id = global_ids(local)
                if (.not. seen(global_id)) then
                    seen(global_id) = .true.
                    canonical(:, global_id) = coordinates(:, local)
                    selection%owner_rank(global_id) = owner_rank(local)
                else
                    distance = sqrt(sum((canonical(:, global_id) - coordinates(:, local))**2))
                    if (distance > coordinate_tolerance .or. &
                        selection%owner_rank(global_id) /= owner_rank(local)) return
                end if
                if (owned(local)) then
                    if (selection%owner_partition(global_id) /= 0) return
                    selection%owner_partition(global_id) = partition
                    selection%owner_local(global_id) = local
                end if
            end do
        end do
        if (any(.not. seen) .or. any(selection%owner_partition == 0)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "physical trace IDs need exactly one owned copy")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine initialize_physical_trace_owner_selection

    subroutine validate_physical_trace_owner_selection(selection, status)
        type(physical_trace_owner_selection_t), intent(in) :: selection
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "physical trace owner selection has invalid metadata")
        if (selection%partition_count < 1 .or. selection%global_count < 1 .or. &
            selection%component_dimension < 1 .or. &
            .not. allocated(selection%owner_partition) .or. &
            .not. allocated(selection%owner_local) .or. .not. allocated(selection%owner_rank)) return
        if (size(selection%owner_partition) /= selection%global_count .or. &
            size(selection%owner_local) /= selection%global_count .or. &
            size(selection%owner_rank) /= selection%global_count .or. &
            any(selection%owner_partition < 1) .or. &
            any(selection%owner_partition > selection%partition_count) .or. &
            any(selection%owner_local < 1) .or. any(selection%owner_rank < 0)) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_physical_trace_owner_selection

    subroutine physical_trace_owner_selection_maps( &
            selection, owner_partition, owner_local, owner_rank, status)
        type(physical_trace_owner_selection_t), intent(in) :: selection
        integer, allocatable, intent(out) :: owner_partition(:), owner_local(:), owner_rank(:)
        type(fortsparse_status_t), intent(out) :: status

        if (allocated(owner_partition)) deallocate(owner_partition)
        if (allocated(owner_local)) deallocate(owner_local)
        if (allocated(owner_rank)) deallocate(owner_rank)
        call validate_physical_trace_owner_selection(selection, status)
        if (status%code /= FORTSPARSE_OK) return
        allocate(owner_partition, source=selection%owner_partition)
        allocate(owner_local, source=selection%owner_local)
        allocate(owner_rank, source=selection%owner_rank)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine physical_trace_owner_selection_maps

    subroutine gather_physical_trace_values( &
            selection, partition_offsets, local_values, global_values, status)
        type(physical_trace_owner_selection_t), intent(in) :: selection
        integer, intent(in) :: partition_offsets(:)
        real(dp), intent(in) :: local_values(:, :)
        real(dp), intent(out) :: global_values(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: global_id, row

        global_values = 0.0_dp
        if (.not. validate_gather_inputs(selection, partition_offsets, local_values, &
            global_values, status)) return
        do global_id = 1, selection%global_count
            row = partition_offsets(selection%owner_partition(global_id)) + &
                selection%owner_local(global_id) - 1
            global_values(global_id, :) = local_values(row, :)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine gather_physical_trace_values

    subroutine gather_physical_trace_values_jvp( &
            selection, partition_offsets, local_values_dot, global_values_dot, status)
        type(physical_trace_owner_selection_t), intent(in) :: selection
        integer, intent(in) :: partition_offsets(:)
        real(dp), intent(in) :: local_values_dot(:, :)
        real(dp), intent(out) :: global_values_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call gather_physical_trace_values(selection, partition_offsets, local_values_dot, &
            global_values_dot, status)
    end subroutine gather_physical_trace_values_jvp

    subroutine gather_physical_trace_values_vjp( &
            selection, partition_offsets, global_values_bar, local_values_bar, status)
        type(physical_trace_owner_selection_t), intent(in) :: selection
        integer, intent(in) :: partition_offsets(:)
        real(dp), intent(in) :: global_values_bar(:, :)
        real(dp), intent(out) :: local_values_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: global_id, row

        local_values_bar = 0.0_dp
        if (.not. validate_gather_vjp_inputs(selection, partition_offsets, global_values_bar, &
            local_values_bar, status)) return
        do global_id = 1, selection%global_count
            row = partition_offsets(selection%owner_partition(global_id)) + &
                selection%owner_local(global_id) - 1
            local_values_bar(row, :) = global_values_bar(global_id, :)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine gather_physical_trace_values_vjp

    logical function validate_gather_inputs(selection, offsets, local_values, global_values, status)
        type(physical_trace_owner_selection_t), intent(in) :: selection
        integer, intent(in) :: offsets(:)
        real(dp), intent(in) :: local_values(:, :), global_values(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: partition, global_id, row

        validate_gather_inputs = .false.
        call validate_physical_trace_owner_selection(selection, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(offsets) /= selection%partition_count + 1 .or. offsets(1) /= 1 .or. &
            any(offsets(2:) <= offsets(:selection%partition_count)) .or. &
            offsets(size(offsets)) - 1 /= size(local_values, 1) .or. &
            size(local_values, 2) /= selection%component_dimension .or. &
            size(global_values, 1) /= selection%global_count .or. &
            size(global_values, 2) /= selection%component_dimension .or. &
            any(.not. ieee_is_finite(local_values))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, "physical trace gather has incompatible values")
            return
        end if
        do global_id = 1, selection%global_count
            partition = selection%owner_partition(global_id)
            row = offsets(partition) + selection%owner_local(global_id) - 1
            if (row < offsets(partition) .or. row >= offsets(partition + 1)) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, "physical trace owner row is invalid")
                return
            end if
        end do
        validate_gather_inputs = .true.
    end function validate_gather_inputs

    logical function validate_gather_vjp_inputs(selection, offsets, global_values_bar, local_values_bar, status)
        type(physical_trace_owner_selection_t), intent(in) :: selection
        integer, intent(in) :: offsets(:)
        real(dp), intent(in) :: global_values_bar(:, :), local_values_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        validate_gather_vjp_inputs = .false.
        call validate_physical_trace_owner_selection(selection, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(offsets) /= selection%partition_count + 1 .or. offsets(1) /= 1 .or. &
            any(offsets(2:) <= offsets(:selection%partition_count)) .or. &
            size(global_values_bar, 1) /= selection%global_count .or. &
            size(global_values_bar, 2) /= selection%component_dimension .or. &
            offsets(size(offsets)) - 1 /= size(local_values_bar, 1) .or. &
            size(local_values_bar, 2) /= selection%component_dimension .or. &
            any(.not. ieee_is_finite(global_values_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, "physical trace gather VJP has incompatible values")
            return
        end if
        validate_gather_vjp_inputs = .true.
    end function validate_gather_vjp_inputs

    subroutine clear_selection(selection)
        type(physical_trace_owner_selection_t), intent(inout) :: selection

        selection%partition_count = 0
        selection%global_count = 0
        selection%component_dimension = 0
        if (allocated(selection%owner_partition)) deallocate(selection%owner_partition)
        if (allocated(selection%owner_local)) deallocate(selection%owner_local)
        if (allocated(selection%owner_rank)) deallocate(selection%owner_rank)
    end subroutine clear_selection

end module fortfem_physical_trace_owner_selection
