program test_tree_cotree_gauge
    use check, only: check_condition, check_summary
    use fortfem_feec, only: apply_tree_cotree_prolongation, &
        apply_tree_cotree_restriction, build_tree_cotree_gauge, &
        build_tree_cotree_dof_map, &
        reduce_tree_cotree_dense_system, reduce_tree_cotree_dense_system_jvp, &
        reduce_tree_cotree_dense_system_vjp, tree_cotree_gauge_edges, &
        tree_cotree_gauge_components, &
        tree_cotree_gauge_t, validate_tree_cotree_gauge
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: triangle_incidence(3, 3) = reshape([ &
        -1, 1, 0, 0, -1, 1, -1, 0, 1], [3, 3])
    integer, parameter :: disconnected_incidence(4, 2) = reshape([ &
        -1, 1, 0, 0, 0, 0, -1, 1], [4, 2])
    integer, parameter :: invalid_incidence(3, 2) = reshape([ &
        -1, 1, 0, 1, 1, -1], [3, 2])
    real(dp), parameter :: full_vector(3) = [10.0_dp, 20.0_dp, 30.0_dp]
    real(dp), parameter :: full_matrix(3, 3) = reshape([ &
        2.0_dp, 0.1_dp, 0.2_dp, 0.1_dp, 3.0_dp, 0.3_dp, &
        0.2_dp, 0.3_dp, 4.0_dp], [3, 3])
    real(dp), parameter :: full_rhs(3) = [1.0_dp, 2.0_dp, 3.0_dp]
    type(tree_cotree_gauge_t) :: gauge, disconnected, invalid
    integer, allocatable :: tree_edges(:), cotree_edges(:)
    integer, allocatable :: components(:)
    real(dp), allocatable :: reduced(:), prolonged(:)
    real(dp), allocatable :: reduced_matrix(:, :), reduced_rhs(:)
    real(dp), allocatable :: reduced_matrix_dot(:, :), reduced_rhs_dot(:)
    real(dp), allocatable :: matrix_bar(:, :), rhs_bar(:)
    logical, allocatable :: constrained_dofs(:)
    integer, allocatable :: free_dofs(:)
    real(dp) :: lhs, rhs
    integer :: component_count, status

    call build_tree_cotree_gauge(triangle_incidence, gauge, status)
    call check_condition(status == 0, "triangle graph builds a tree-cotree gauge")
    call validate_tree_cotree_gauge(gauge, status)
    call tree_cotree_gauge_edges(gauge, tree_edges, cotree_edges, status)
    call check_condition(status == 0 .and. size(tree_edges) == 2 .and. &
        size(cotree_edges) == 1 .and. all(tree_edges == [1, 2]) .and. &
        all(cotree_edges == [3]), &
        "triangle gauge selects a spanning tree and one cotree edge")
    call tree_cotree_gauge_components( &
        gauge, components, component_count, status)
    call check_condition(status == 0 .and. component_count == 1 .and. &
        all(components == [1, 1, 1]) .and. &
        size(tree_edges) == size(components) - component_count .and. &
        size(cotree_edges) == size(triangle_incidence, 2) - &
            size(components) + component_count, &
        "triangle gauge exposes deterministic component and cycle ranks")

    call apply_tree_cotree_restriction(gauge, full_vector, reduced, status)
    call apply_tree_cotree_prolongation(gauge, reduced, prolonged, status)
    call check_condition(status == 0 .and. size(reduced) == 1 .and. &
        abs(reduced(1) - 30.0_dp) < 1.0e-14_dp .and. &
        maxval(abs(prolonged - [0.0_dp, 0.0_dp, 30.0_dp])) < 1.0e-14_dp, &
        "fixed tree gauge restricts and prolongs cotree coefficients")

    call reduce_tree_cotree_dense_system( &
        gauge, full_matrix, full_rhs, reduced_matrix, reduced_rhs, status)
    call check_condition(status == 0 .and. size(reduced_matrix, 1) == 1 .and. &
        abs(reduced_matrix(1, 1) - 4.0_dp) < 1.0e-14_dp .and. &
        abs(reduced_rhs(1) - 3.0_dp) < 1.0e-14_dp, &
        "tree-cotree reduction extracts the direct cotree system")

    call reduce_tree_cotree_dense_system_jvp( &
        gauge, full_matrix, full_rhs, reduced_matrix_dot, reduced_rhs_dot, status)
    call reduce_tree_cotree_dense_system_vjp( &
        gauge, reduced_matrix, reduced_rhs, matrix_bar, rhs_bar, status)
    lhs = sum(reduced_matrix*reduced_matrix_dot) + &
        dot_product(reduced_rhs, reduced_rhs_dot)
    rhs = sum(matrix_bar*full_matrix) + dot_product(rhs_bar, full_rhs)
    call check_condition(status == 0 .and. abs(lhs - rhs) < 1.0e-14_dp, &
        "fixed tree-cotree reduction has a real adjoint identity")

    call build_tree_cotree_gauge(disconnected_incidence, disconnected, status)
    call tree_cotree_gauge_edges(disconnected, tree_edges, cotree_edges, status)
    call tree_cotree_gauge_components( &
        disconnected, components, component_count, status)
    call check_condition(status == 0 .and. size(tree_edges) == 2 .and. &
        size(cotree_edges) == 0 .and. component_count == 2 .and. &
        all(components == [1, 1, 2, 2]) .and. &
        size(tree_edges) == size(components) - component_count, &
        "disconnected graph exposes a deterministic spanning forest")

    call build_tree_cotree_dof_map( &
        gauge, [2, 4, 5], 6, constrained_dofs, free_dofs, status)
    call check_condition(status == 0 .and. size(constrained_dofs) == 6 .and. &
        all(constrained_dofs .eqv. &
            [.false., .true., .false., .true., .false., .false.]) &
        .and. all(free_dofs == [1, 3, 5, 6]), &
        "tree-cotree map constrains control edges and retains IGA moments")

    call build_tree_cotree_dof_map( &
        gauge, [2, 2, 5], 6, constrained_dofs, free_dofs, status)
    call check_condition(status /= 0, &
        "tree-cotree DOF map rejects ambiguous control-edge numbering")

    call build_tree_cotree_gauge(invalid_incidence, invalid, status)
    call check_condition(status /= 0, "tree-cotree gauge rejects invalid incidence")

    call check_summary("tree-cotree gauge")
end program test_tree_cotree_gauge
