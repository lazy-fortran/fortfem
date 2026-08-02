program test_biro_tree_cotree_benchmark
    !! Independent oracle for the small Bíró/Preis/Richter gauge fixture.
    use check, only: check_condition, check_summary
    use fortfem_feec, only: apply_tree_cotree_prolongation, &
        apply_tree_cotree_restriction, build_tree_cotree_gauge, &
        reduce_tree_cotree_dense_system, tree_cotree_gauge_edges, &
        tree_cotree_gauge_t
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: incidence(4, 5) = reshape([ &
        -1, 1, 0, 0, 0, -1, 1, 0, 0, 0, -1, 1, 1, 0, 0, -1, &
        -1, 0, 1, 0], [4, 5])
    ! Independent topological curl-curl Gram matrix C^T C for two faces.
    real(dp), parameter :: full_matrix(5, 5) = reshape([ &
        1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, &
        1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
        -1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [5, 5])
    real(dp), parameter :: expected_full(5) = [-0.35_dp, 0.50_dp, -0.60_dp, &
        1.20_dp, -0.25_dp]
    real(dp), parameter :: expected_projected(5) = [0.0_dp, 0.0_dp, 0.0_dp, &
        1.20_dp, -0.25_dp]
    real(dp), parameter :: expected_gauge(5) = [0.0_dp, 0.0_dp, 0.0_dp, &
        0.75_dp, -0.40_dp]
    real(dp), parameter :: expected_rhs(5) = [0.4_dp, 0.4_dp, 0.35_dp, &
        0.35_dp, -0.05_dp]
    real(dp), parameter :: expected_gradient(5) = [-0.35_dp, 0.50_dp, &
        -0.60_dp, 0.45_dp, 0.15_dp]
    real(dp) :: full_rhs(5)
    real(dp), allocatable :: reduced_matrix(:, :), reduced_rhs(:)
    real(dp), allocatable :: restricted(:), prolonged(:)
    integer, allocatable :: tree_edges(:), cotree_edges(:)
    type(tree_cotree_gauge_t) :: gauge
    integer :: status

    full_rhs = expected_rhs
    call check_condition(maxval(abs(matmul(full_matrix, expected_full) - &
        expected_rhs)) < 1.0e-13_dp, "literal curl-curl source is consistent")
    call check_condition(maxval(abs(matmul(full_matrix, expected_gradient))) < &
        1.0e-13_dp, "discrete gradient lies in the curl-curl nullspace")
    call build_tree_cotree_gauge(incidence, gauge, status)
    call check_condition(status == 0, "Biro graph admits a tree-cotree gauge")
    if (status /= 0) then
        call check_summary("Biro tree-cotree benchmark")
        return
    end if
    call tree_cotree_gauge_edges(gauge, tree_edges, cotree_edges, status)
    call check_condition(status == 0, "Biro graph edge extraction succeeds")
    if (status /= 0) then
        call check_summary("Biro tree-cotree benchmark")
        return
    end if
    call check_condition(size(tree_edges) == 3 .and. size(cotree_edges) == 2, &
        "Biro graph has three tree and two cotree edges")
    call apply_tree_cotree_restriction(gauge, expected_full, restricted, status)
    call check_condition(status == 0, "tree gauge restriction succeeds")
    if (status /= 0) then
        call check_summary("Biro tree-cotree benchmark")
        return
    end if
    call apply_tree_cotree_prolongation(gauge, restricted, prolonged, status)
    call check_condition(status == 0, "tree gauge prolongation succeeds")
    if (status /= 0) then
        call check_summary("Biro tree-cotree benchmark")
        return
    end if
    call check_condition(maxval(abs(prolonged - expected_projected)) < 1.0e-14_dp, &
        "tree gauge projection removes tree coefficients")
    call reduce_tree_cotree_dense_system( &
        gauge, full_matrix, full_rhs, reduced_matrix, reduced_rhs, status)
    call check_condition(status == 0, "Biro reduced curl-curl system succeeds")
    if (status == 0) then
        call check_condition(size(reduced_matrix, 1) == 2 .and. &
            size(reduced_matrix, 2) == 2, "Biro reduced block has two cotree DOFs")
        if (size(reduced_matrix, 1) == 2 .and. size(reduced_matrix, 2) == 2) then
            restricted = expected_gauge(4:5)
            call check_condition(maxval(abs(reduced_matrix - reshape([ &
                1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2]))) < 1.0e-13_dp, &
                "Biro reduced curl-curl block matches independent oracle")
            call check_condition(maxval(abs(reduced_rhs - [0.35_dp, -0.05_dp])) < &
                1.0e-13_dp, "Biro reduced RHS matches independent face source")
            call check_condition(maxval(abs(reduced_rhs - &
                matmul(reduced_matrix, restricted))) < 1.0e-13_dp, &
                "Biro reduced RHS matches the reconstructed cotree solution")
        end if
    end if
    call check_summary("Biro tree-cotree benchmark")
end program test_biro_tree_cotree_benchmark
