module fortfem_physical_trace_reconciliation
    !! Fixed-topology reconciliation of independently ordered physical traces.
    !!
    !! A trace row is identified by the global ID in a
    !! ``physical_trace_ownership_t`` record.  This contract builds the
    !! candidate-to-reference permutation and the sign induced by local trace
    !! orientations.  It then applies that linear map to values, tangents, or
    !! cotangents.  Coordinates are checked before the map is accepted; no
    !! communicator, mesh search, constitutive law, or physical convention is
    !! selected here.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_physical_trace_ownership, only: &
        physical_trace_ownership_t, physical_trace_ownership_maps, &
        physical_trace_ownership_point_count, compare_physical_trace_coordinates
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    type, public :: physical_trace_reconciliation_t
        private
        integer :: point_count = 0
        integer :: component_dimension = 0
        integer, allocatable :: candidate_index(:)
        integer, allocatable :: orientation(:)
        real(dp) :: maximum_distance = 0.0_dp
    contains
        procedure, private :: assign_physical_trace_reconciliation
        generic :: assignment(=) => assign_physical_trace_reconciliation
    end type physical_trace_reconciliation_t

    public :: initialize_physical_trace_reconciliation
    public :: validate_physical_trace_reconciliation
    public :: physical_trace_reconciliation_maps
    public :: reconcile_physical_trace_values
    public :: reconcile_physical_trace_values_jvp
    public :: reconcile_physical_trace_values_vjp

contains

    subroutine assign_physical_trace_reconciliation(lhs, rhs)
        class(physical_trace_reconciliation_t), intent(out) :: lhs
        type(physical_trace_reconciliation_t), intent(in) :: rhs

        lhs%point_count = rhs%point_count
        lhs%component_dimension = rhs%component_dimension
        lhs%maximum_distance = rhs%maximum_distance
        if (allocated(rhs%candidate_index)) &
            allocate(lhs%candidate_index, source=rhs%candidate_index)
        if (allocated(rhs%orientation)) &
            allocate(lhs%orientation, source=rhs%orientation)
    end subroutine assign_physical_trace_reconciliation

    subroutine initialize_physical_trace_reconciliation( &
            reconciliation, reference, candidate, reference_orientation, &
            candidate_orientation, coordinate_tolerance, status)
        type(physical_trace_reconciliation_t), intent(inout) :: reconciliation
        type(physical_trace_ownership_t), intent(in) :: reference, candidate
        integer, intent(in) :: reference_orientation(:), candidate_orientation(:)
        real(dp), intent(in) :: coordinate_tolerance
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: reference_coordinates(:, :), candidate_coordinates(:, :)
        integer, allocatable :: reference_ids(:), candidate_ids(:)
        integer, allocatable :: reference_owner(:), candidate_owner(:)
        logical, allocatable :: reference_owned(:), candidate_owned(:)
        integer :: reference_rank, candidate_rank, mismatch_count
        integer :: reference_point, candidate_point
        real(dp) :: maximum_distance

        call clear_reconciliation(reconciliation)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "physical trace reconciliation received invalid orientation metadata")
        if (.not. ieee_is_finite(coordinate_tolerance) .or. &
            coordinate_tolerance < 0.0_dp) return
        if (physical_trace_ownership_point_count(reference) < 1 .or. &
            physical_trace_ownership_point_count(candidate) < 1 .or. &
            size(reference_orientation) /= physical_trace_ownership_point_count(reference) .or. &
            size(candidate_orientation) /= physical_trace_ownership_point_count(candidate)) return
        if (any(abs(reference_orientation) /= 1) .or. &
            any(abs(candidate_orientation) /= 1)) return
        call physical_trace_ownership_maps( &
            reference, reference_coordinates, reference_ids, reference_owner, &
            reference_owned, reference_rank, status)
        if (status%code /= FORTSPARSE_OK) return
        call physical_trace_ownership_maps( &
            candidate, candidate_coordinates, candidate_ids, candidate_owner, &
            candidate_owned, candidate_rank, status)
        if (status%code /= FORTSPARSE_OK) return
        call compare_physical_trace_coordinates( &
            reference, candidate, coordinate_tolerance, mismatch_count, &
            maximum_distance, status)
        if (status%code /= FORTSPARSE_OK) return
        if (mismatch_count /= 0) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "physical trace reconciliation found coordinate or ID mismatch")
            return
        end if
        reconciliation%point_count = size(reference_ids)
        reconciliation%component_dimension = size(reference_coordinates, 1)
        reconciliation%maximum_distance = maximum_distance
        allocate(reconciliation%candidate_index(reconciliation%point_count))
        allocate(reconciliation%orientation(reconciliation%point_count))
        do reference_point = 1, reconciliation%point_count
            candidate_point = find_global_id( &
                candidate_ids, reference_ids(reference_point))
            if (candidate_point == 0) then
                call clear_reconciliation(reconciliation)
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "physical trace reconciliation found an unmapped global ID")
                return
            end if
            reconciliation%candidate_index(reference_point) = candidate_point
            reconciliation%orientation(reference_point) = &
                reference_orientation(reference_point)*candidate_orientation(candidate_point)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine initialize_physical_trace_reconciliation

    subroutine validate_physical_trace_reconciliation(reconciliation, status)
        type(physical_trace_reconciliation_t), intent(in) :: reconciliation
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "physical trace reconciliation has invalid metadata")
        if (reconciliation%point_count < 1 .or. &
            reconciliation%component_dimension < 1 .or. &
            .not. ieee_is_finite(reconciliation%maximum_distance) .or. &
            reconciliation%maximum_distance < 0.0_dp) return
        if (.not. allocated(reconciliation%candidate_index) .or. &
            .not. allocated(reconciliation%orientation)) return
        if (size(reconciliation%candidate_index) /= reconciliation%point_count .or. &
            size(reconciliation%orientation) /= reconciliation%point_count .or. &
            any(reconciliation%candidate_index < 1) .or. &
            any(abs(reconciliation%orientation) /= 1)) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_physical_trace_reconciliation

    subroutine physical_trace_reconciliation_maps( &
            reconciliation, candidate_index, orientation, maximum_distance, status)
        type(physical_trace_reconciliation_t), intent(in) :: reconciliation
        integer, allocatable, intent(out) :: candidate_index(:), orientation(:)
        real(dp), intent(out) :: maximum_distance
        type(fortsparse_status_t), intent(out) :: status

        if (allocated(candidate_index)) deallocate(candidate_index)
        if (allocated(orientation)) deallocate(orientation)
        maximum_distance = 0.0_dp
        call validate_physical_trace_reconciliation(reconciliation, status)
        if (status%code /= FORTSPARSE_OK) return
        allocate(candidate_index, source=reconciliation%candidate_index)
        allocate(orientation, source=reconciliation%orientation)
        maximum_distance = reconciliation%maximum_distance
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine physical_trace_reconciliation_maps

    subroutine reconcile_physical_trace_values( &
            reconciliation, candidate_values, reference_values, status)
        type(physical_trace_reconciliation_t), intent(in) :: reconciliation
        real(dp), intent(in) :: candidate_values(:, :)
        real(dp), intent(out) :: reference_values(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: reference_point

        reference_values = 0.0_dp
        if (.not. valid_value_shapes( &
                reconciliation, candidate_values, reference_values, status)) return
        do reference_point = 1, reconciliation%point_count
            reference_values(reference_point, :) = &
                reconciliation%orientation(reference_point)*candidate_values( &
                reconciliation%candidate_index(reference_point), :)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine reconcile_physical_trace_values

    subroutine reconcile_physical_trace_values_jvp( &
            reconciliation, candidate_values_dot, reference_values_dot, status)
        type(physical_trace_reconciliation_t), intent(in) :: reconciliation
        real(dp), intent(in) :: candidate_values_dot(:, :)
        real(dp), intent(out) :: reference_values_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call reconcile_physical_trace_values( &
            reconciliation, candidate_values_dot, reference_values_dot, status)
    end subroutine reconcile_physical_trace_values_jvp

    subroutine reconcile_physical_trace_values_vjp( &
            reconciliation, reference_values_bar, candidate_values_bar, status)
        type(physical_trace_reconciliation_t), intent(in) :: reconciliation
        real(dp), intent(in) :: reference_values_bar(:, :)
        real(dp), intent(out) :: candidate_values_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: reference_point

        candidate_values_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "physical trace reconciliation VJP received invalid values")
        if (.not. valid_vjp_shapes( &
                reconciliation, reference_values_bar, candidate_values_bar, status)) return
        do reference_point = 1, reconciliation%point_count
            candidate_values_bar(reconciliation%candidate_index(reference_point), :) = &
                reconciliation%orientation(reference_point)* &
                reference_values_bar(reference_point, :)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine reconcile_physical_trace_values_vjp

    logical function valid_value_shapes( &
            reconciliation, candidate_values, reference_values, status) result(valid)
        type(physical_trace_reconciliation_t), intent(in) :: reconciliation
        real(dp), intent(in) :: candidate_values(:, :), reference_values(:, :)
        type(fortsparse_status_t), intent(out) :: status

        valid = .false.
        call validate_physical_trace_reconciliation(reconciliation, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(candidate_values, 1) < 1 .or. size(candidate_values, 2) < 1 .or. &
            size(reference_values, 1) /= reconciliation%point_count .or. &
            size(candidate_values, 1) < maxval(reconciliation%candidate_index) .or. &
            size(reference_values, 2) /= size(candidate_values, 2) .or. &
            any(.not. ieee_is_finite(candidate_values))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "physical trace reconciliation received incompatible values")
            return
        end if
        valid = .true.
    end function valid_value_shapes

    logical function valid_vjp_shapes( &
            reconciliation, reference_values_bar, candidate_values_bar, status) result(valid)
        type(physical_trace_reconciliation_t), intent(in) :: reconciliation
        real(dp), intent(in) :: reference_values_bar(:, :)
        real(dp), intent(in) :: candidate_values_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        valid = .false.
        call validate_physical_trace_reconciliation(reconciliation, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(reference_values_bar, 1) /= reconciliation%point_count .or. &
            size(reference_values_bar, 2) < 1 .or. &
            size(candidate_values_bar, 1) < maxval(reconciliation%candidate_index) .or. &
            size(candidate_values_bar, 2) /= size(reference_values_bar, 2) .or. &
            any(.not. ieee_is_finite(reference_values_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "physical trace reconciliation VJP received incompatible values")
            return
        end if
        valid = .true.
    end function valid_vjp_shapes

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

    subroutine clear_reconciliation(reconciliation)
        type(physical_trace_reconciliation_t), intent(inout) :: reconciliation

        reconciliation%point_count = 0
        reconciliation%component_dimension = 0
        reconciliation%maximum_distance = 0.0_dp
        if (allocated(reconciliation%candidate_index)) deallocate( &
            reconciliation%candidate_index)
        if (allocated(reconciliation%orientation)) deallocate(reconciliation%orientation)
    end subroutine clear_reconciliation

end module fortfem_physical_trace_reconciliation
