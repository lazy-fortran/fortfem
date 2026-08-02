module fortfem_physical_trace_ownership
    !! Physical-coordinate ownership ledger for distributed trace rows.
    !!
    !! The existing partition ledger identifies rows by integer IDs.  This
    !! companion contract carries the physical coordinates needed to verify
    !! that independently generated FEM/BEM/DtN/IGA partitions refer to the
    !! same trace point before an MPI exchange is introduced.  It deliberately
    !! owns no communicator, mesh, or geometry search and does not differentiate
    !! discrete owner selectors.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    type, public :: physical_trace_ownership_t
        private
        integer :: dimension = 0
        integer :: point_count = 0
        integer :: rank = -1
        real(dp) :: coordinate_tolerance = 0.0_dp
        real(dp), allocatable :: coordinates(:, :)
        integer, allocatable :: global_ids(:)
        integer, allocatable :: owner_rank(:)
        logical, allocatable :: owned(:)
    end type physical_trace_ownership_t

    public :: initialize_physical_trace_ownership
    public :: validate_physical_trace_ownership
    public :: physical_trace_ownership_maps
    public :: physical_trace_ownership_dimension
    public :: physical_trace_ownership_point_count
    public :: physical_trace_ownership_rank
    public :: compare_physical_trace_coordinates

contains

    subroutine initialize_physical_trace_ownership( &
            ownership, coordinates, global_ids, owner_rank, owned, rank, &
            coordinate_tolerance, status)
        type(physical_trace_ownership_t), intent(inout) :: ownership
        real(dp), intent(in) :: coordinates(:, :)
        integer, intent(in) :: global_ids(:), owner_rank(:)
        logical, intent(in) :: owned(:)
        integer, intent(in) :: rank
        real(dp), intent(in) :: coordinate_tolerance
        type(fortsparse_status_t), intent(out) :: status

        call clear_ownership(ownership)
        call validate_raw_ownership( &
            coordinates, global_ids, owner_rank, owned, rank, &
            coordinate_tolerance, status)
        if (status%code /= FORTSPARSE_OK) return
        ownership%dimension = size(coordinates, 1)
        ownership%point_count = size(coordinates, 2)
        ownership%rank = rank
        ownership%coordinate_tolerance = coordinate_tolerance
        allocate(ownership%coordinates, source=coordinates)
        allocate(ownership%global_ids, source=global_ids)
        allocate(ownership%owner_rank, source=owner_rank)
        allocate(ownership%owned, source=owned)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine initialize_physical_trace_ownership

    subroutine validate_physical_trace_ownership(ownership, status)
        type(physical_trace_ownership_t), intent(in) :: ownership
        type(fortsparse_status_t), intent(out) :: status

        call validate_raw_ownership( &
            ownership%coordinates, ownership%global_ids, ownership%owner_rank, &
            ownership%owned, ownership%rank, ownership%coordinate_tolerance, status)
        if (status%code /= FORTSPARSE_OK) return
        if (ownership%dimension /= size(ownership%coordinates, 1) .or. &
            ownership%point_count /= size(ownership%coordinates, 2)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "physical trace ownership counts are inconsistent")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_physical_trace_ownership

    subroutine physical_trace_ownership_maps( &
            ownership, coordinates, global_ids, owner_rank, owned, rank, status)
        type(physical_trace_ownership_t), intent(in) :: ownership
        real(dp), allocatable, intent(out) :: coordinates(:, :)
        integer, allocatable, intent(out) :: global_ids(:), owner_rank(:)
        logical, allocatable, intent(out) :: owned(:)
        integer, intent(out) :: rank
        type(fortsparse_status_t), intent(out) :: status

        if (allocated(coordinates)) deallocate(coordinates)
        if (allocated(global_ids)) deallocate(global_ids)
        if (allocated(owner_rank)) deallocate(owner_rank)
        if (allocated(owned)) deallocate(owned)
        rank = -1
        call validate_physical_trace_ownership(ownership, status)
        if (status%code /= FORTSPARSE_OK) return
        allocate(coordinates, source=ownership%coordinates)
        allocate(global_ids, source=ownership%global_ids)
        allocate(owner_rank, source=ownership%owner_rank)
        allocate(owned, source=ownership%owned)
        rank = ownership%rank
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine physical_trace_ownership_maps

    integer function physical_trace_ownership_dimension(ownership)
        type(physical_trace_ownership_t), intent(in) :: ownership

        physical_trace_ownership_dimension = ownership%dimension
    end function physical_trace_ownership_dimension

    integer function physical_trace_ownership_point_count(ownership)
        type(physical_trace_ownership_t), intent(in) :: ownership

        physical_trace_ownership_point_count = ownership%point_count
    end function physical_trace_ownership_point_count

    integer function physical_trace_ownership_rank(ownership)
        type(physical_trace_ownership_t), intent(in) :: ownership

        physical_trace_ownership_rank = ownership%rank
    end function physical_trace_ownership_rank

    subroutine compare_physical_trace_coordinates( &
            reference, candidate, coordinate_tolerance, mismatch_count, &
            maximum_distance, status)
        !! Compare two partitions by global ID, independent of local ordering.
        type(physical_trace_ownership_t), intent(in) :: reference, candidate
        real(dp), intent(in) :: coordinate_tolerance
        integer, intent(out) :: mismatch_count
        real(dp), intent(out) :: maximum_distance
        type(fortsparse_status_t), intent(out) :: status

        integer :: reference_point, candidate_point
        real(dp) :: distance

        mismatch_count = 0
        maximum_distance = 0.0_dp
        call validate_physical_trace_ownership(reference, status)
        if (status%code /= FORTSPARSE_OK) return
        call validate_physical_trace_ownership(candidate, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. ieee_is_finite(coordinate_tolerance) .or. &
            coordinate_tolerance < 0.0_dp .or. &
            physical_trace_ownership_dimension(reference) /= &
            physical_trace_ownership_dimension(candidate)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "physical trace coordinate comparison has incompatible inputs")
            return
        end if
        do reference_point = 1, reference%point_count
            candidate_point = find_global_id( &
                candidate%global_ids, reference%global_ids(reference_point))
            if (candidate_point == 0) then
                mismatch_count = mismatch_count + 1
                cycle
            end if
            distance = sqrt(sum((reference%coordinates(:, reference_point) - &
                candidate%coordinates(:, candidate_point))**2))
            maximum_distance = max(maximum_distance, distance)
            if (distance > coordinate_tolerance) mismatch_count = mismatch_count + 1
        end do
        do candidate_point = 1, candidate%point_count
            if (find_global_id(reference%global_ids, &
                candidate%global_ids(candidate_point)) == 0) mismatch_count = &
                mismatch_count + 1
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine compare_physical_trace_coordinates

    subroutine validate_raw_ownership( &
            coordinates, global_ids, owner_rank, owned, rank, tolerance, status)
        real(dp), intent(in) :: coordinates(:, :), tolerance
        integer, intent(in) :: global_ids(:), owner_rank(:), rank
        logical, intent(in) :: owned(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: first, second

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "physical trace ownership received incompatible metadata")
        if (size(coordinates, 1) < 1 .or. size(coordinates, 2) < 1 .or. &
            size(global_ids) /= size(coordinates, 2) .or. &
            size(owner_rank) /= size(global_ids) .or. size(owned) /= size(global_ids) .or. &
            rank < 0 .or. .not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) return
        if (.not. all(ieee_is_finite(coordinates)) .or. &
            any(global_ids < 1) .or. any(owner_rank < 0)) return
        do first = 1, size(global_ids)
            if (owned(first) .neqv. owner_rank(first) == rank) return
            do second = first + 1, size(global_ids)
                if (global_ids(first) == global_ids(second)) return
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_raw_ownership

    integer function find_global_id(global_ids, global_id) result(index)
        integer, intent(in) :: global_ids(:), global_id

        integer :: candidate

        index = 0
        do candidate = 1, size(global_ids)
            if (global_ids(candidate) == global_id) then
                index = candidate
                return
            end if
        end do
    end function find_global_id

    subroutine clear_ownership(ownership)
        type(physical_trace_ownership_t), intent(inout) :: ownership

        ownership%dimension = 0
        ownership%point_count = 0
        ownership%rank = -1
        ownership%coordinate_tolerance = 0.0_dp
        if (allocated(ownership%coordinates)) deallocate(ownership%coordinates)
        if (allocated(ownership%global_ids)) deallocate(ownership%global_ids)
        if (allocated(ownership%owner_rank)) deallocate(ownership%owner_rank)
        if (allocated(ownership%owned)) deallocate(ownership%owned)
    end subroutine clear_ownership

end module fortfem_physical_trace_ownership
