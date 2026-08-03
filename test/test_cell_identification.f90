program test_cell_identification
    use check, only: check_condition, check_summary
    use fortfem_core, only: cell_complex_betti_numbers, cell_complex_t, &
        initialize_cell_complex, quotient_cell_complex, validate_cell_complex
    use fortfem_feec, only: cell_identification_classes, &
        identify_boundary_matrix, cell_identification_t, &
        initialize_cell_identification, validate_cell_identification
    implicit none

    type(cell_complex_t) :: interval, quotient, invalid_complex
    type(cell_identification_t) :: identity, periodic, lower, upper, invalid
    integer, allocatable :: representative(:), orientation(:), classes(:)
    integer, allocatable :: boundary(:, :), quotient_boundary(:, :)
    integer :: betti(4), class_count, status

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

    deallocate(representative, orientation)
    allocate(representative(2), orientation(2), boundary(2, 1))
    representative = [1, 1]
    orientation = 1
    boundary = reshape([-1, 1], [2, 1])
    call initialize_cell_identification( &
        lower, representative, orientation, status)
    call initialize_cell_identification( &
        upper, [1], [1], status)
    call identify_boundary_matrix( &
        lower, upper, boundary, quotient_boundary, status)
    call check_condition(status == 0 .and. &
        size(quotient_boundary, 1) == 1 .and. &
        size(quotient_boundary, 2) == 1 .and. &
        quotient_boundary(1, 1) == 0, &
        "endpoint identification maps an interval boundary to a circle")

    deallocate(representative, orientation, boundary)
    allocate(representative(1), orientation(1), boundary(1, 2))
    representative = 1
    orientation = 1
    boundary = reshape([1, -1], [1, 2])
    call initialize_cell_identification( &
        lower, representative, orientation, status)
    deallocate(representative, orientation)
    allocate(representative(2), orientation(2))
    representative = [1, 1]
    orientation = [1, -1]
    call initialize_cell_identification( &
        upper, representative, orientation, status)
    call identify_boundary_matrix( &
        lower, upper, boundary, quotient_boundary, status)
    call check_condition(status == 0 .and. &
        quotient_boundary(1, 1) == 1, &
        "signed upper identification preserves boundary orientation")

    boundary = reshape([1, 1], [1, 2])
    call identify_boundary_matrix( &
        lower, upper, boundary, quotient_boundary, status)
    call check_condition(status /= 0, &
        "inconsistent identified boundary columns are rejected")

    deallocate(representative, orientation, boundary)
    allocate(representative(4), orientation(4))
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

    allocate(boundary(2, 1))
    boundary = reshape([-1, 1], [2, 1])
    call initialize_cell_complex( &
        interval, 2, boundary, status=status)
    call check_condition(status == 0, &
        "interval for quotient-complex test initializes")
    call initialize_cell_identification( &
        lower, [1, 1], [1, 1], status)
    call initialize_cell_identification( &
        upper, [1], [1], status)
    call quotient_cell_complex( &
        interval, lower, upper, quotient=quotient, status=status)
    call check_condition(status == 0, &
        "signed endpoint identification builds quotient complex")
    call validate_cell_complex(quotient, status)
    call cell_complex_betti_numbers(quotient, betti, status)
    call check_condition(status == 0 .and. all(betti == [1, 1, 0, 0]), &
        "interval quotient has the circle homology oracle")
    deallocate(boundary)
    if (allocated(quotient_boundary)) deallocate(quotient_boundary)
    allocate(boundary(1, 2), quotient_boundary(2, 1))
    boundary = reshape([-1, 1], [1, 2])
    quotient_boundary = reshape([1, 0], [2, 1])
    call initialize_cell_complex( &
        invalid_complex, 1, boundary, quotient_boundary, status=status)
    call initialize_cell_identification( &
        lower, [1], [1], status)
    call initialize_cell_identification( &
        upper, [1, 2], [1, 1], status)
    call initialize_cell_identification( &
        identity, [1], [1], status)
    call quotient_cell_complex( &
        invalid_complex, lower, upper, quotient=quotient, status=status, &
        face_identification=identity)
    call check_condition(status /= 0, &
        "malformed source complex is rejected before quotienting")

    call check_summary("signed cell identification")
end program test_cell_identification
