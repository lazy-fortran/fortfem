module fortfem_cell_identification
    !! Signed representative maps for periodic or quotient cell IDs.
    !!
    !! `representative(i)` is the canonical cell ID for cell i and
    !! `orientation(i)` is its orientation relative to that representative.
    !! The representative map is intentionally metadata only: constructing a
    !! quotient mesh or transforming a chain complex remains a higher-level
    !! operation with geometry-dependent choices.
    implicit none
    private

    type, public :: cell_identification_t
        private
        integer, allocatable :: representative(:)
        integer, allocatable :: orientation(:)
    end type cell_identification_t

    public :: initialize_cell_identification
    public :: validate_cell_identification
    public :: cell_identification_classes

contains

    subroutine initialize_cell_identification( &
            identification, representative, orientation, status)
        type(cell_identification_t), intent(inout) :: identification
        integer, intent(in) :: representative(:), orientation(:)
        integer, intent(out) :: status
        integer :: cell_count

        call clear_cell_identification(identification)
        status = 1
        if (size(representative) /= size(orientation)) then
            status = 2
            return
        end if
        cell_count = size(representative)
        if (any(representative < 1) .or. &
            any(representative > cell_count)) then
            status = 3
            return
        end if
        if (any(abs(orientation) /= 1)) then
            status = 4
            return
        end if
        allocate(identification%representative(cell_count))
        allocate(identification%orientation(cell_count))
        identification%representative = representative
        identification%orientation = orientation
        status = 0
    end subroutine initialize_cell_identification

    subroutine validate_cell_identification(identification, status)
        type(cell_identification_t), intent(in) :: identification
        integer, intent(out) :: status
        integer :: cell, cell_count

        status = 1
        if (.not. allocated(identification%representative) .or. &
            .not. allocated(identification%orientation)) return
        cell_count = size(identification%representative)
        if (size(identification%orientation) /= cell_count) then
            status = 2
            return
        end if
        if (any(identification%representative < 1) .or. &
            any(identification%representative > cell_count)) then
            status = 3
            return
        end if
        if (any(abs(identification%orientation) /= 1)) then
            status = 4
            return
        end if
        do cell = 1, cell_count
            if (identification%representative( &
                identification%representative(cell)) /= &
                identification%representative(cell)) then
                status = 5
                return
            end if
            if (identification%representative(cell) == cell) then
                if (identification%orientation(cell) /= 1) then
                    status = 6
                    return
                end if
            end if
        end do
        status = 0
    end subroutine validate_cell_identification

    subroutine cell_identification_classes( &
            identification, classes, class_count, status)
        type(cell_identification_t), intent(in) :: identification
        integer, allocatable, intent(out) :: classes(:)
        integer, intent(out) :: class_count
        integer, intent(out) :: status
        integer, allocatable :: labels(:)
        integer :: cell, representative

        if (allocated(classes)) deallocate(classes)
        class_count = 0
        call validate_cell_identification(identification, status)
        if (status /= 0) return
        allocate(classes(size(identification%representative)))
        if (size(classes) == 0) then
            status = 0
            return
        end if
        allocate(labels(size(classes)))
        labels = 0
        do cell = 1, size(classes)
            representative = identification%representative(cell)
            if (labels(representative) == 0) then
                class_count = class_count + 1
                labels(representative) = class_count
            end if
            classes(cell) = labels(representative)
        end do
        status = 0
    end subroutine cell_identification_classes

    subroutine clear_cell_identification(identification)
        type(cell_identification_t), intent(inout) :: identification

        if (allocated(identification%representative)) then
            deallocate(identification%representative)
        end if
        if (allocated(identification%orientation)) then
            deallocate(identification%orientation)
        end if
    end subroutine clear_cell_identification

end module fortfem_cell_identification
