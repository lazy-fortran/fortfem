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
    real(dp), parameter :: full_matrix(5, 5) = reshape([ &
        4.0_dp, 0.1_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        0.1_dp, 5.0_dp, 0.2_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.2_dp, 6.0_dp, 0.0_dp, 0.1_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, 7.0_dp, 0.3_dp, &
        0.0_dp, 0.0_dp, 0.1_dp, 0.3_dp, 8.0_dp], [5, 5])
    real(dp), parameter :: expected_full(5) = [0.0_dp, 0.0_dp, 0.0_dp, &
        0.75_dp, -0.4_dp]
    real(dp) :: full_rhs(5)
    real(dp), allocatable :: reduced_matrix(:, :), reduced_rhs(:)
    real(dp), allocatable :: restricted(:), prolonged(:)
    integer, allocatable :: tree_edges(:), cotree_edges(:)
    type(tree_cotree_gauge_t) :: gauge
    integer :: status

    full_rhs = matmul(full_matrix, expected_full)
    call build_tree_cotree_gauge(incidence, gauge, status)
    call check_condition(status == 0, "Biro graph admits a tree-cotree gauge")
    call tree_cotree_gauge_edges(gauge, tree_edges, cotree_edges, status)
    call check_condition(status == 0 .and. size(tree_edges) == 3 .and. &
        size(cotree_edges) == 2, "Biro graph has three tree and two cotree edges")
    call apply_tree_cotree_restriction(gauge, expected_full, restricted, status)
    call apply_tree_cotree_prolongation(gauge, restricted, prolonged, status)
    call check_condition(status == 0 .and. maxval(abs(prolonged - expected_full)) < &
        1.0e-14_dp, "tree gauge preserves the physical cotree coefficients")
    call reduce_tree_cotree_dense_system( &
        gauge, full_matrix, full_rhs, reduced_matrix, reduced_rhs, status)
    call check_condition(status == 0 .and. maxval(abs(reduced_rhs - &
        matmul(reduced_matrix, restricted))) < 1.0e-13_dp, &
        "Biro reduced curl-curl system matches independent manufactured oracle")
    call check_summary("Biro tree-cotree benchmark")
end program test_biro_tree_cotree_benchmark
