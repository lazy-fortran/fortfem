program test_cell_complex
    use check, only: check_condition, check_summary
    use fortfem_api, only: cell_complex_betti_numbers, cell_complex_cocycle_basis, &
        cell_complex_cycle_basis, cell_complex_harmonic_one_forms, &
        cell_complex_euler_characteristic, cell_complex_t, &
        initialize_cell_complex, validate_cell_complex
    use fortfem_kinds, only: dp
    implicit none

    type(cell_complex_t) :: interval, circle, sphere, torus, invalid
    real(dp), allocatable :: cycles(:, :), harmonic(:, :), edge_metric(:, :)
    integer, allocatable :: boundary_1(:, :), boundary_2(:, :)
    integer :: betti(4), euler, status
    integer :: cocycle_count, cycle_count, harmonic_count

    allocate(boundary_1(2, 1))
    boundary_1 = reshape([-1, 1], shape(boundary_1))
    call initialize_cell_complex( &
        interval, 2, boundary_1, status=status)
    call check_condition(status == 0, "interval complex initializes")
    call validate_cell_complex(interval, status)
    call check_condition(status == 0, "interval complex validates")
    call cell_complex_euler_characteristic(interval, euler, status)
    call check_condition(status == 0 .and. euler == 1, &
        "interval Euler characteristic is one")
    call cell_complex_betti_numbers(interval, betti, status)
    call check_condition(status == 0 .and. all(betti == [1, 0, 0, 0]), &
        "interval Betti numbers are independently known")
    call cell_complex_cycle_basis(interval, cycles, cycle_count, status)
    call check_condition(status == 0 .and. cycle_count == 0 .and. &
        size(cycles, 1) == 1, &
        "interval has no nontrivial one-cycle")
    allocate(edge_metric(1, 1))
    edge_metric = 1.0_dp
    call cell_complex_harmonic_one_forms( &
        interval, edge_metric, harmonic, harmonic_count, status)
    call check_condition(status == 0 .and. harmonic_count == 0, &
        "interval has no metric harmonic one-form")
    deallocate(edge_metric)

    deallocate(boundary_1)
    allocate(boundary_1(1, 1))
    boundary_1 = 0
    call initialize_cell_complex(circle, 1, boundary_1, status=status)
    call check_condition(status == 0, "circle complex initializes")
    call cell_complex_euler_characteristic(circle, euler, status)
    call cell_complex_betti_numbers(circle, betti, status)
    call check_condition(status == 0 .and. euler == 0 .and. &
        all(betti == [1, 1, 0, 0]), &
        "circle Euler and Betti numbers match the loop oracle")
    call cell_complex_cycle_basis(circle, cycles, cycle_count, status)
    call check_condition(status == 0 .and. cycle_count == 1 .and. &
        size(cycles, 1) == 1 .and. abs(abs(cycles(1, 1)) - 1.0_dp) < 1.0e-14_dp .and. &
        maxval(abs(matmul(real(boundary_1, dp), cycles))) < 1.0e-14_dp, &
        "circle cycle basis spans the exact integer kernel")
    call cell_complex_cocycle_basis(circle, cycles, cocycle_count, status)
    call check_condition(status == 0 .and. cocycle_count == 1 .and. &
        size(cycles, 1) == 1, "circle has one independent one-cocycle")
    allocate(edge_metric(1, 1))
    edge_metric = 2.0_dp
    call cell_complex_harmonic_one_forms( &
        circle, edge_metric, harmonic, harmonic_count, status)
    call check_condition(status == 0 .and. harmonic_count == 1 .and. &
        size(harmonic, 1) == 1, "circle has one metric harmonic one-form")
    deallocate(edge_metric)

    deallocate(boundary_1)
    allocate(boundary_1(1, 0), boundary_2(0, 1))
    call initialize_cell_complex( &
        sphere, 1, boundary_1, boundary_2, status=status)
    call check_condition(status == 0, "sphere CW complex initializes")
    call cell_complex_euler_characteristic(sphere, euler, status)
    call cell_complex_betti_numbers(sphere, betti, status)
    call check_condition(status == 0 .and. euler == 2 .and. &
        all(betti == [1, 0, 1, 0]), &
        "zero-cell sphere boundary has the expected Euler oracle")

    deallocate(boundary_1, boundary_2)
    allocate(boundary_1(1, 2), boundary_2(2, 1))
    boundary_1 = 0
    boundary_2 = 0
    call initialize_cell_complex( &
        torus, 1, boundary_1, boundary_2, status=status)
    call check_condition(status == 0, "torus CW complex initializes")
    call cell_complex_euler_characteristic(torus, euler, status)
    call cell_complex_betti_numbers(torus, betti, status)
    call check_condition(status == 0 .and. euler == 0 .and. &
        all(betti == [1, 2, 1, 0]), &
        "torus CW complex matches the independent homology oracle")
    call cell_complex_cycle_basis(torus, cycles, cycle_count, status)
    call check_condition(status == 0 .and. cycle_count == 2 .and. &
        maxval(abs(matmul(real(boundary_1, dp), cycles))) < 1.0e-14_dp, &
        "torus cycle basis has two independent closed one-forms")
    call cell_complex_cocycle_basis(torus, cycles, cocycle_count, status)
    call check_condition(status == 0 .and. cocycle_count == 2 .and. &
        maxval(abs(matmul(real(transpose(boundary_2), dp), cycles))) < &
        1.0e-14_dp, "torus cocycle basis annihilates the face boundary")
    allocate(edge_metric(2, 2))
    edge_metric = reshape([2.0_dp, 0.5_dp, 0.5_dp, 3.0_dp], [2, 2])
    call cell_complex_harmonic_one_forms( &
        torus, edge_metric, harmonic, harmonic_count, status)
    call check_condition(status == 0 .and. harmonic_count == 2 .and. &
        maxval(abs(matmul(real(boundary_1, dp), matmul(edge_metric, harmonic)))) < &
        1.0e-14_dp, "torus metric harmonic basis satisfies the co-closed oracle")
    edge_metric(2, 2) = -1.0_dp
    call cell_complex_harmonic_one_forms( &
        torus, edge_metric, harmonic, harmonic_count, status)
    call check_condition(status /= 0, "harmonic basis rejects a non-positive metric")

    boundary_1 = reshape([-1, 1], [1, 2])
    boundary_2 = reshape([1, 0], [2, 1])
    call initialize_cell_complex( &
        invalid, 1, boundary_1, boundary_2, status=status)
    call check_condition(status == 0, &
        "invalid boundary algebra can be constructed for rejection")
    call validate_cell_complex(invalid, status)
    call check_condition(status /= 0, &
        "boundary-of-boundary violation is rejected")

    call check_summary("oriented cell complex")
end program test_cell_complex
