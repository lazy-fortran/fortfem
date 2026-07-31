program test_cell_identification
    use check, only: check_condition, check_summary
    use fortfem_api, only: cell_identification_classes, &
        cell_identification_t, initialize_cell_identification, &
        validate_cell_identification
    implicit none

    type(cell_identification_t) :: identity, periodic, invalid
    integer, allocatable :: representative(:), orientation(:), classes(:)
    integer :: class_count, status

    allocate(representative(3), orientation(3))
    representative = [1, 2, 3]
    orientation = 1
    call initialize_cell_identification( &
        identity, representative, orientation, status)
    call check_condition(status == 0, "identity cell map initializes")
    call validate_cell_identification(identity, status)
    call check_condition(status == 0, "identity cell map validates")
    call cell_identification_classes( &
        identity, classes, class_count, status)
    call check_condition(status == 0 .and. class_count == 3 .and. &
        all(classes == [1, 2, 3]), &
        "identity map has one class per cell")

    deallocate(representative, orientation)
    allocate(representative(4), orientation(4))
    representative = [1, 1, 3, 3]
    orientation = [1, -1, 1, -1]
    call initialize_cell_identification( &
        periodic, representative, orientation, status)
    call cell_identification_classes( &
        periodic, classes, class_count, status)
    call check_condition(status == 0 .and. class_count == 2 .and. &
        all(classes == [1, 1, 2, 2]), &
        "periodic signed map has the independent quotient classes")

    representative = [2, 1, 3, 3]
    orientation = 1
    call initialize_cell_identification( &
        invalid, representative, orientation, status)
    call check_condition(status == 0, &
        "cyclic representative map can be constructed for rejection")
    call validate_cell_identification(invalid, status)
    call check_condition(status /= 0, &
        "non-idempotent representative cycle is rejected")

    representative = 1
    orientation = [ -1, 1, 1 ]
    call initialize_cell_identification( &
        invalid, representative, orientation, status)
    call check_condition(status /= 0, &
        "canonical representative with reversed orientation is rejected")

    orientation = 0
    call initialize_cell_identification( &
        invalid, representative, orientation, status)
    call check_condition(status /= 0, "zero orientation is rejected")

    deallocate(representative, orientation)
    allocate(representative(1), orientation(2))
    call initialize_cell_identification( &
        invalid, representative, orientation, status)
    call check_condition(status /= 0, &
        "mismatched identification arrays are rejected")

    call check_summary("signed cell identification")
end program test_cell_identification
