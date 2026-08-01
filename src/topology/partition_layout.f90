module fortfem_partition_layout
    !! Serial partition metadata and deterministic reduction contract.
    !!
    !! A layout records the local-to-global IDs and owner rank of each local
    !! entry.  It deliberately contains no communicator or MPI state: local
    !! kernels can therefore be tested and run unchanged before a distributed
    !! backend is introduced.  Reduction routines accumulate in local index
    !! order, giving the serial/no-op implementation a deterministic result.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    type, public :: partition_layout_t
        private
        integer :: global_count = 0
        integer :: local_count = 0
        integer :: owned_count = 0
        integer :: ghost_count = 0
        integer :: rank = -1
        integer, allocatable :: local_to_global(:)
        integer, allocatable :: owner_rank(:)
        logical, allocatable :: owned(:)
    end type partition_layout_t

    public :: initialize_partition_layout
    public :: validate_partition_layout
    public :: partition_layout_maps
    public :: partition_layout_global_count
    public :: partition_layout_owned_count
    public :: partition_layout_ghost_count
    public :: assemble_partitioned_sum
    public :: assemble_partitioned_sum_jvp
    public :: assemble_partitioned_sum_vjp

contains

    subroutine initialize_partition_layout( &
            layout, global_count, local_to_global, owner_rank, owned, rank, status)
        !! Initialize fixed local-to-global and owner/ghost metadata.
        type(partition_layout_t), intent(inout) :: layout
        integer, intent(in) :: global_count, local_to_global(:), owner_rank(:)
        logical, intent(in) :: owned(:)
        integer, intent(in) :: rank
        type(fortsparse_status_t), intent(out) :: status

        call clear_layout(layout)
        call validate_raw_layout( &
            global_count, local_to_global, owner_rank, owned, rank, status)
        if (status%code /= FORTSPARSE_OK) return
        layout%global_count = global_count
        layout%local_count = size(local_to_global)
        layout%owned_count = count(owned)
        layout%ghost_count = layout%local_count - layout%owned_count
        layout%rank = rank
        allocate(layout%local_to_global(layout%local_count))
        allocate(layout%owner_rank(layout%local_count))
        allocate(layout%owned(layout%local_count))
        layout%local_to_global = local_to_global
        layout%owner_rank = owner_rank
        layout%owned = owned
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine initialize_partition_layout

    subroutine validate_partition_layout(layout, status)
        !! Validate a previously initialized partition layout.
        type(partition_layout_t), intent(in) :: layout
        type(fortsparse_status_t), intent(out) :: status

        call validate_raw_layout( &
            layout%global_count, layout%local_to_global, layout%owner_rank, &
            layout%owned, layout%rank, status)
        if (status%code /= FORTSPARSE_OK) return
        if (layout%local_count /= size(layout%local_to_global) .or. &
            layout%owned_count /= count(layout%owned) .or. &
            layout%ghost_count /= layout%local_count - layout%owned_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "partition layout counts are inconsistent")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_partition_layout

    subroutine partition_layout_maps( &
            layout, local_to_global, owner_rank, owned, rank, status)
        !! Export caller-owned copies of the fixed partition metadata.
        type(partition_layout_t), intent(in) :: layout
        integer, allocatable, intent(out) :: local_to_global(:), owner_rank(:)
        logical, allocatable, intent(out) :: owned(:)
        integer, intent(out) :: rank
        type(fortsparse_status_t), intent(out) :: status

        if (allocated(local_to_global)) deallocate(local_to_global)
        if (allocated(owner_rank)) deallocate(owner_rank)
        if (allocated(owned)) deallocate(owned)
        rank = -1
        call validate_partition_layout(layout, status)
        if (status%code /= FORTSPARSE_OK) return
        allocate(local_to_global(layout%local_count))
        allocate(owner_rank(layout%local_count))
        allocate(owned(layout%local_count))
        local_to_global = layout%local_to_global
        owner_rank = layout%owner_rank
        owned = layout%owned
        rank = layout%rank
    end subroutine partition_layout_maps

    integer function partition_layout_global_count(layout)
        type(partition_layout_t), intent(in) :: layout

        partition_layout_global_count = layout%global_count
    end function partition_layout_global_count

    integer function partition_layout_owned_count(layout)
        type(partition_layout_t), intent(in) :: layout

        partition_layout_owned_count = layout%owned_count
    end function partition_layout_owned_count

    integer function partition_layout_ghost_count(layout)
        type(partition_layout_t), intent(in) :: layout

        partition_layout_ghost_count = layout%ghost_count
    end function partition_layout_ghost_count

    subroutine assemble_partitioned_sum(layout, local_values, global_values, status)
        !! Accumulate a local vector into the serial global reduction buffer.
        type(partition_layout_t), intent(in) :: layout
        real(dp), intent(in) :: local_values(:)
        real(dp), intent(out) :: global_values(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: local

        global_values = 0.0_dp
        call validate_reduction_inputs(layout, local_values, global_values, status)
        if (status%code /= FORTSPARSE_OK) return
        do local = 1, layout%local_count
            global_values(layout%local_to_global(local)) = &
                global_values(layout%local_to_global(local)) + local_values(local)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_partitioned_sum

    subroutine assemble_partitioned_sum_jvp(layout, local_values_dot, global_values_dot, status)
        !! Apply the linear JVP of the serial partition reduction.
        type(partition_layout_t), intent(in) :: layout
        real(dp), intent(in) :: local_values_dot(:)
        real(dp), intent(out) :: global_values_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        call assemble_partitioned_sum(layout, local_values_dot, global_values_dot, status)
    end subroutine assemble_partitioned_sum_jvp

    subroutine assemble_partitioned_sum_vjp(layout, global_values_bar, local_values_bar, status)
        !! Apply the transpose of the serial partition reduction.
        type(partition_layout_t), intent(in) :: layout
        real(dp), intent(in) :: global_values_bar(:)
        real(dp), intent(out) :: local_values_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: local

        local_values_bar = 0.0_dp
        call validate_reduction_inputs(layout, local_values_bar, global_values_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(local_values_bar) /= layout%local_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "partition reduction VJP has an incompatible local cotangent")
            return
        end if
        do local = 1, layout%local_count
            local_values_bar(local) = global_values_bar( &
                layout%local_to_global(local))
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_partitioned_sum_vjp

    subroutine validate_raw_layout( &
            global_count, local_to_global, owner_rank, owned, rank, status)
        integer, intent(in) :: global_count, local_to_global(:), owner_rank(:), rank
        logical, intent(in) :: owned(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: first, second

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "partition layout received incompatible metadata")
        if (global_count < 1 .or. rank < 0 .or. size(local_to_global) < 1 .or. &
            size(owner_rank) /= size(local_to_global) .or. &
            size(owned) /= size(local_to_global)) return
        if (any(local_to_global < 1) .or. &
            any(local_to_global > global_count) .or. any(owner_rank < 0)) return
        do first = 1, size(local_to_global)
            if (owned(first)) then
                if (owner_rank(first) /= rank) return
            else
                if (owner_rank(first) == rank) return
            end if
            do second = first + 1, size(local_to_global)
                if (local_to_global(first) == local_to_global(second)) return
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_raw_layout

    subroutine validate_reduction_inputs(layout, local_values, global_values, status)
        type(partition_layout_t), intent(in) :: layout
        real(dp), intent(in) :: local_values(:), global_values(:)
        type(fortsparse_status_t), intent(out) :: status

        call validate_partition_layout(layout, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(local_values) /= layout%local_count .or. &
            size(global_values) /= layout%global_count .or. &
            any(.not. ieee_is_finite(local_values)) .or. &
            any(.not. ieee_is_finite(global_values))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "partition reduction received incompatible vectors")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_reduction_inputs

    subroutine clear_layout(layout)
        type(partition_layout_t), intent(inout) :: layout

        layout%global_count = 0
        layout%local_count = 0
        layout%owned_count = 0
        layout%ghost_count = 0
        layout%rank = -1
        if (allocated(layout%local_to_global)) deallocate(layout%local_to_global)
        if (allocated(layout%owner_rank)) deallocate(layout%owner_rank)
        if (allocated(layout%owned)) deallocate(layout%owned)
    end subroutine clear_layout

end module fortfem_partition_layout
