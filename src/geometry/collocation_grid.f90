module fortfem_collocation_grid
    !! Neutral tensor-product collocation metadata.
    !!
    !! This contract describes a bounded two-coordinate grid without making
    !! assumptions about a PDE, a constitutive law, or a particular geometry.
    !! The first coordinate is normally radial and the second coordinate is
    !! normally periodic, but callers may use the same metadata for any tensor
    !! product sampling.  Flattening follows Fortran column-major order:
    !! ``flat = i + (j - 1)*n_first``.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    integer, parameter, public :: COLLOCATION_GRID_LINEAR = 1
    integer, parameter, public :: COLLOCATION_GRID_QUADRATURE = 2
    integer, parameter, public :: COLLOCATION_GRID_CONCENTRIC = 3

    type, public :: collocation_grid_t
        private
        integer :: grid_kind = 0
        integer :: first_count = 0
        integer :: second_count = 0
        integer :: point_count = 0
        integer :: chunk_size = 0
        integer :: chunk_count = 0
        real(dp), allocatable :: first_nodes(:)
        real(dp), allocatable :: second_nodes(:)
        real(dp), allocatable :: first_weights(:)
        real(dp), allocatable :: second_weights(:)
        real(dp), allocatable :: weights(:)
        integer, allocatable :: periodic_endpoint_map(:)
        integer, allocatable :: axis_map(:)
        integer, allocatable :: chunk_start(:)
        integer, allocatable :: chunk_end(:)
    contains
        procedure, private :: assign_collocation_grid
        generic :: assignment(=) => assign_collocation_grid
    end type collocation_grid_t

    public :: initialize_collocation_grid
    public :: validate_collocation_grid
    public :: collocation_grid_metadata
    public :: collocation_grid_flat_index
    public :: collocation_grid_unflatten_index
    public :: collocation_grid_chunk_bounds
    public :: collocation_grid_point_count

contains

    subroutine assign_collocation_grid(lhs, rhs)
        class(collocation_grid_t), intent(out) :: lhs
        type(collocation_grid_t), intent(in) :: rhs

        lhs%grid_kind = rhs%grid_kind
        lhs%first_count = rhs%first_count
        lhs%second_count = rhs%second_count
        lhs%point_count = rhs%point_count
        lhs%chunk_size = rhs%chunk_size
        lhs%chunk_count = rhs%chunk_count
        if (allocated(rhs%first_nodes)) allocate(lhs%first_nodes, source=rhs%first_nodes)
        if (allocated(rhs%second_nodes)) allocate(lhs%second_nodes, source=rhs%second_nodes)
        if (allocated(rhs%first_weights)) allocate(lhs%first_weights, source=rhs%first_weights)
        if (allocated(rhs%second_weights)) allocate(lhs%second_weights, source=rhs%second_weights)
        if (allocated(rhs%weights)) allocate(lhs%weights, source=rhs%weights)
        if (allocated(rhs%periodic_endpoint_map)) then
            allocate(lhs%periodic_endpoint_map, source=rhs%periodic_endpoint_map)
        end if
        if (allocated(rhs%axis_map)) allocate(lhs%axis_map, source=rhs%axis_map)
        if (allocated(rhs%chunk_start)) allocate(lhs%chunk_start, source=rhs%chunk_start)
        if (allocated(rhs%chunk_end)) allocate(lhs%chunk_end, source=rhs%chunk_end)
    end subroutine assign_collocation_grid

    subroutine initialize_collocation_grid( &
            grid, grid_kind, first_nodes, second_nodes, first_weights, &
            second_weights, periodic_endpoint_map, axis_map, chunk_size, status)
        type(collocation_grid_t), intent(out) :: grid
        integer, intent(in) :: grid_kind
        real(dp), intent(in) :: first_nodes(:), second_nodes(:)
        real(dp), intent(in) :: first_weights(:), second_weights(:)
        integer, intent(in), optional :: periodic_endpoint_map(:), axis_map(:)
        integer, intent(in), optional :: chunk_size
        type(fortsparse_status_t), intent(out) :: status

        integer :: first, second, point, selected_chunk_size

        call clear_collocation_grid(grid)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "collocation grid received incompatible metadata")
        selected_chunk_size = size(first_nodes)*size(second_nodes)
        if (present(chunk_size)) selected_chunk_size = chunk_size
        if (.not. valid_raw_grid( &
            grid_kind, first_nodes, second_nodes, first_weights, second_weights, &
            periodic_endpoint_map, axis_map, selected_chunk_size)) return

        grid%grid_kind = grid_kind
        grid%first_count = size(first_nodes)
        grid%second_count = size(second_nodes)
        grid%point_count = grid%first_count*grid%second_count
        grid%chunk_size = selected_chunk_size
        grid%chunk_count = (grid%point_count + selected_chunk_size - 1) / &
            selected_chunk_size
        allocate(grid%first_nodes(grid%first_count), &
            grid%second_nodes(grid%second_count), &
            grid%first_weights(grid%first_count), &
            grid%second_weights(grid%second_count), grid%weights(grid%point_count), &
            grid%periodic_endpoint_map(grid%second_count), &
            grid%axis_map(grid%first_count), &
            grid%chunk_start(grid%chunk_count), grid%chunk_end(grid%chunk_count))
        grid%first_nodes = first_nodes
        grid%second_nodes = second_nodes
        grid%first_weights = first_weights
        grid%second_weights = second_weights
        point = 0
        do second = 1, grid%second_count
            do first = 1, grid%first_count
                point = point + 1
                grid%weights(point) = first_weights(first)*second_weights(second)
            end do
        end do
        grid%periodic_endpoint_map = 0
        if (present(periodic_endpoint_map)) then
            grid%periodic_endpoint_map = periodic_endpoint_map
        else
            do second = 1, grid%second_count
                grid%periodic_endpoint_map(second) = second
            end do
        end if
        grid%axis_map = 0
        if (present(axis_map)) then
            grid%axis_map = axis_map
        else
            do first = 1, grid%first_count
                grid%axis_map(first) = first
            end do
        end if
        do point = 1, grid%chunk_count
            grid%chunk_start(point) = (point - 1)*selected_chunk_size + 1
            grid%chunk_end(point) = min(point*selected_chunk_size, grid%point_count)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine initialize_collocation_grid

    logical function validate_collocation_grid(grid, status) result(valid)
        type(collocation_grid_t), intent(in) :: grid
        type(fortsparse_status_t), intent(out) :: status

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "collocation grid has inconsistent metadata")
        if (.not. allocated(grid%first_nodes) .or. &
            .not. allocated(grid%second_nodes) .or. &
            .not. allocated(grid%first_weights) .or. &
            .not. allocated(grid%second_weights) .or. &
            .not. allocated(grid%periodic_endpoint_map) .or. &
            .not. allocated(grid%axis_map)) return
        if (.not. valid_raw_grid( &
            grid%grid_kind, grid%first_nodes, grid%second_nodes, &
            grid%first_weights, grid%second_weights, &
            grid%periodic_endpoint_map, grid%axis_map, grid%chunk_size)) return
        if (grid%first_count /= size(grid%first_nodes) .or. &
            grid%second_count /= size(grid%second_nodes) .or. &
            grid%point_count /= grid%first_count*grid%second_count .or. &
            .not. allocated(grid%weights) .or. size(grid%weights) /= grid%point_count .or. &
            any(.not. ieee_is_finite(grid%weights)) .or. any(grid%weights <= 0.0_dp)) return
        if (.not. allocated(grid%chunk_start) .or. &
            .not. allocated(grid%chunk_end) .or. &
            size(grid%chunk_start) /= grid%chunk_count .or. &
            size(grid%chunk_end) /= grid%chunk_count) return
        if (grid%chunk_count /= (grid%point_count + grid%chunk_size - 1) / &
            grid%chunk_size) return
        if (grid%chunk_start(1) /= 1 .or. &
            grid%chunk_end(grid%chunk_count) /= grid%point_count) return
        if (any(grid%chunk_start < 1) .or. any(grid%chunk_end > grid%point_count) .or. &
            any(grid%chunk_start > grid%chunk_end)) return
        if (grid%chunk_count > 1) then
            if (any(grid%chunk_start(2:) /= grid%chunk_end(:grid%chunk_count-1) + 1)) return
        end if
        if (any(grid%chunk_end - grid%chunk_start + 1 > grid%chunk_size)) return
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_collocation_grid

    subroutine collocation_grid_metadata( &
            grid, grid_kind, first_nodes, second_nodes, weights, &
            periodic_endpoint_map, axis_map, chunk_size, chunk_start, chunk_end, status)
        type(collocation_grid_t), intent(in) :: grid
        integer, intent(out) :: grid_kind, chunk_size
        real(dp), allocatable, intent(out) :: first_nodes(:), second_nodes(:), weights(:)
        integer, allocatable, intent(out) :: periodic_endpoint_map(:), axis_map(:)
        integer, allocatable, intent(out) :: chunk_start(:), chunk_end(:)
        type(fortsparse_status_t), intent(out) :: status

        grid_kind = 0
        chunk_size = 0
        if (allocated(first_nodes)) deallocate(first_nodes)
        if (allocated(second_nodes)) deallocate(second_nodes)
        if (allocated(weights)) deallocate(weights)
        if (allocated(periodic_endpoint_map)) deallocate(periodic_endpoint_map)
        if (allocated(axis_map)) deallocate(axis_map)
        if (allocated(chunk_start)) deallocate(chunk_start)
        if (allocated(chunk_end)) deallocate(chunk_end)
        if (.not. validate_collocation_grid(grid, status)) return
        grid_kind = grid%grid_kind
        chunk_size = grid%chunk_size
        allocate(first_nodes(size(grid%first_nodes)), second_nodes(size(grid%second_nodes)), &
            weights(size(grid%weights)), periodic_endpoint_map(size(grid%periodic_endpoint_map)), &
            axis_map(size(grid%axis_map)), chunk_start(size(grid%chunk_start)), &
            chunk_end(size(grid%chunk_end)))
        first_nodes = grid%first_nodes
        second_nodes = grid%second_nodes
        weights = grid%weights
        periodic_endpoint_map = grid%periodic_endpoint_map
        axis_map = grid%axis_map
        chunk_start = grid%chunk_start
        chunk_end = grid%chunk_end
    end subroutine collocation_grid_metadata

    integer function collocation_grid_flat_index(grid, first_index, second_index, status)
        type(collocation_grid_t), intent(in) :: grid
        integer, intent(in) :: first_index, second_index
        type(fortsparse_status_t), intent(out) :: status

        collocation_grid_flat_index = 0
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "collocation grid flatten index is out of bounds")
        if (.not. validate_collocation_grid(grid, status)) return
        if (first_index < 1 .or. first_index > grid%first_count .or. &
            second_index < 1 .or. second_index > grid%second_count) return
        collocation_grid_flat_index = first_index + &
            (second_index - 1)*grid%first_count
        call status_set(status, FORTSPARSE_OK, "")
    end function collocation_grid_flat_index

    subroutine collocation_grid_unflatten_index( &
            grid, flat_index, first_index, second_index, status)
        type(collocation_grid_t), intent(in) :: grid
        integer, intent(in) :: flat_index
        integer, intent(out) :: first_index, second_index
        type(fortsparse_status_t), intent(out) :: status

        first_index = 0
        second_index = 0
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "collocation grid unflatten index is out of bounds")
        if (.not. validate_collocation_grid(grid, status)) return
        if (flat_index < 1 .or. flat_index > grid%point_count) return
        second_index = (flat_index - 1)/grid%first_count + 1
        first_index = flat_index - (second_index - 1)*grid%first_count
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine collocation_grid_unflatten_index

    subroutine collocation_grid_chunk_bounds( &
            grid, chunk, first_point, last_point, status)
        type(collocation_grid_t), intent(in) :: grid
        integer, intent(in) :: chunk
        integer, intent(out) :: first_point, last_point
        type(fortsparse_status_t), intent(out) :: status

        first_point = 0
        last_point = 0
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "collocation grid chunk is out of bounds")
        if (.not. validate_collocation_grid(grid, status)) return
        if (chunk < 1 .or. chunk > grid%chunk_count) return
        first_point = grid%chunk_start(chunk)
        last_point = grid%chunk_end(chunk)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine collocation_grid_chunk_bounds

    integer function collocation_grid_point_count(grid)
        type(collocation_grid_t), intent(in) :: grid

        collocation_grid_point_count = grid%point_count
    end function collocation_grid_point_count

    logical function valid_raw_grid( &
            grid_kind, first_nodes, second_nodes, first_weights, second_weights, &
            periodic_endpoint_map, axis_map, chunk_size) result(valid)
        integer, intent(in) :: grid_kind, chunk_size
        real(dp), intent(in) :: first_nodes(:), second_nodes(:), first_weights(:), &
            second_weights(:)
        integer, intent(in), optional :: periodic_endpoint_map(:), axis_map(:)
        integer :: index

        valid = .false.
        if (grid_kind < COLLOCATION_GRID_LINEAR .or. &
            grid_kind > COLLOCATION_GRID_CONCENTRIC) return
        if (size(first_nodes) < 1 .or. size(second_nodes) < 1 .or. &
            size(first_weights) /= size(first_nodes) .or. &
            size(second_weights) /= size(second_nodes) .or. chunk_size < 1 .or. &
            chunk_size > size(first_nodes)*size(second_nodes)) return
        if (any(.not. ieee_is_finite(first_nodes)) .or. &
            any(.not. ieee_is_finite(second_nodes)) .or. &
            any(.not. ieee_is_finite(first_weights)) .or. &
            any(.not. ieee_is_finite(second_weights))) return
        if (any(first_weights <= 0.0_dp) .or. any(second_weights <= 0.0_dp)) return
        if (size(first_nodes) > 1 .and. any(first_nodes(2:) <= first_nodes(:size(first_nodes)-1))) return
        if (size(second_nodes) > 1 .and. any(second_nodes(2:) <= second_nodes(:size(second_nodes)-1))) return
        if (grid_kind == COLLOCATION_GRID_CONCENTRIC .and. any(first_nodes < 0.0_dp)) return
        if (present(periodic_endpoint_map)) then
            if (size(periodic_endpoint_map) /= size(second_nodes)) return
            if (.not. valid_idempotent_map(periodic_endpoint_map)) return
        end if
        if (present(axis_map)) then
            if (size(axis_map) /= size(first_nodes)) return
            if (.not. valid_idempotent_map(axis_map)) return
        end if
        do index = 1, size(first_nodes)
            if (first_weights(index) <= 0.0_dp) return
        end do
        valid = .true.
    end function valid_raw_grid

    logical function valid_idempotent_map(mapping) result(valid)
        integer, intent(in) :: mapping(:)
        integer :: index

        valid = size(mapping) > 0 .and. all(mapping >= 1) .and. &
            all(mapping <= size(mapping))
        if (.not. valid) return
        do index = 1, size(mapping)
            if (mapping(mapping(index)) /= mapping(index)) then
                valid = .false.
                return
            end if
        end do
    end function valid_idempotent_map

    subroutine clear_collocation_grid(grid)
        type(collocation_grid_t), intent(out) :: grid

        grid%grid_kind = 0
        grid%first_count = 0
        grid%second_count = 0
        grid%point_count = 0
        grid%chunk_size = 0
        grid%chunk_count = 0
        if (allocated(grid%first_nodes)) deallocate(grid%first_nodes)
        if (allocated(grid%second_nodes)) deallocate(grid%second_nodes)
        if (allocated(grid%first_weights)) deallocate(grid%first_weights)
        if (allocated(grid%second_weights)) deallocate(grid%second_weights)
        if (allocated(grid%weights)) deallocate(grid%weights)
        if (allocated(grid%periodic_endpoint_map)) deallocate(grid%periodic_endpoint_map)
        if (allocated(grid%axis_map)) deallocate(grid%axis_map)
        if (allocated(grid%chunk_start)) deallocate(grid%chunk_start)
        if (allocated(grid%chunk_end)) deallocate(grid%chunk_end)
    end subroutine clear_collocation_grid

end module fortfem_collocation_grid
