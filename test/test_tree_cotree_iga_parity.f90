program test_tree_cotree_iga_parity
    use check, only: check_condition, check_summary
    use fortfem_feec, only: diagnose_tree_cotree_iga_invariance, &
        tree_cotree_iga_parity_t
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: incidence(3, 3) = reshape([ &
        -1, 1, 0, 0, -1, 1, -1, 0, 1], [3, 3])
    integer, parameter :: control_edge_local(3) = [1, 2, 3]
    integer, parameter :: signed_map_a(5) = [1, -2, 3, 2, 4]
    integer, parameter :: signed_map_b(5) = [-1, 2, -3, -2, 4]
    real(dp), parameter :: local_matrix(5, 5) = reshape([ &
        4.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 5.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 6.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp], [5, 5])
    real(dp), parameter :: local_rhs(5) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    real(dp), parameter :: period_weights(4) = [0.7_dp, -1.2_dp, 2.0_dp, 0.4_dp]
    type(tree_cotree_iga_parity_t) :: diagnostics
    integer :: status

    call diagnose_tree_cotree_iga_invariance(incidence, control_edge_local, &
        signed_map_a, signed_map_b, local_matrix, local_rhs, period_weights, &
        diagnostics, status)
    call check_condition(status == 0, &
        "IGA signed-map tree-cotree parity accepts a fixed graph")
    call check_condition(diagnostics%global_dof_count == 4 .and. &
        diagnostics%free_dof_count == 2 .and. &
        diagnostics%matrix_error < 1.0e-14_dp .and. &
        diagnostics%rhs_error < 1.0e-14_dp, &
        "signed IGA assembly has the expected fixed-gauge direct system")
    call check_condition(diagnostics%solution_error < 1.0e-14_dp .and. &
        diagnostics%period_error < 1.0e-14_dp .and. &
        abs(diagnostics%period_value - 2.0_dp) < 1.0e-14_dp, &
        "tree-cotree solve and oriented period are invariant under IGA signs")

    call diagnose_tree_cotree_iga_invariance(incidence, control_edge_local, &
        signed_map_a, [1, 2, 3, 2, -4], local_matrix, local_rhs, &
        period_weights, diagnostics, status)
    call check_condition(status /= 0, &
        "IGA parity rejects a control-edge map with an ambiguous class")
    call check_summary("IGA tree-cotree parity")
end program test_tree_cotree_iga_parity
