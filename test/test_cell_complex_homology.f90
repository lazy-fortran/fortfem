program test_cell_complex_homology
    use check, only: check_condition, check_summary
    use fortfem_core, only: cell_complex_cohomology_cocycle_basis, &
        cell_complex_homology_cycle_basis, cell_complex_t, initialize_cell_complex
    use fortfem_kinds, only: dp
    implicit none

    type(cell_complex_t) :: disk, torus
    integer, allocatable :: boundary_1(:, :), boundary_2(:, :)
    real(dp), allocatable :: homology(:, :), cohomology(:, :)
    integer :: homology_count, cohomology_count, status

    allocate(boundary_1(3, 3), boundary_2(3, 1))
    boundary_1 = reshape([ &
        -1, 0, 1, &
         1, -1, 0, &
         0, 1, -1], [3, 3])
    boundary_2 = 1
    call initialize_cell_complex(disk, 3, boundary_1, boundary_2, status=status)
    call check_condition(status == 0, "oriented disk complex initializes")
    call cell_complex_homology_cycle_basis( &
        disk, homology, homology_count, status)
    call check_condition(status == 0 .and. homology_count == 0 .and. &
        size(homology, 1) == 3 .and. size(homology, 2) == 0, &
        "a filled loop has no one-dimensional homology")
    call cell_complex_cohomology_cocycle_basis( &
        disk, cohomology, cohomology_count, status)
    call check_condition(status == 0 .and. cohomology_count == 0 .and. &
        size(cohomology, 1) == 3 .and. size(cohomology, 2) == 0, &
        "a filled loop has no one-dimensional cohomology")

    deallocate(boundary_1, boundary_2)
    allocate(boundary_1(1, 2), boundary_2(2, 1))
    boundary_1 = 0
    boundary_2 = 0
    call initialize_cell_complex(torus, 1, boundary_1, boundary_2, status=status)
    call check_condition(status == 0, "periodic torus complex initializes")
    call cell_complex_homology_cycle_basis( &
        torus, homology, homology_count, status)
    call check_condition(status == 0 .and. homology_count == 2 .and. &
        size(homology, 1) == 2 .and. &
        maxval(abs(matmul(real(boundary_1, dp), homology))) < 1.0e-14_dp, &
        "the torus returns two independent closed one-cycle representatives")
    call cell_complex_cohomology_cocycle_basis( &
        torus, cohomology, cohomology_count, status)
    call check_condition(status == 0 .and. cohomology_count == 2 .and. &
        size(cohomology, 1) == 2 .and. &
        maxval(abs(matmul(real(transpose(boundary_2), dp), cohomology))) < &
        1.0e-14_dp, "the torus returns two independent cohomology representatives")

    call check_summary("cell-complex homology representatives")
end program test_cell_complex_homology
