program biro_tree_cotree_benchmark
    !! Small, reproducible tree--cotree curl--curl fixture.
    !!
    !! This is a license-safe manufactured reproduction of the direct-gauge
    !! topology used by Biro, Preis, and Richter.  It intentionally exposes
    !! the physical edge potential and its loop-current proxy before any
    !! residual or convergence diagnostics.
    use fortfem_feec, only: apply_tree_cotree_prolongation, &
        apply_tree_cotree_restriction, build_tree_cotree_gauge, &
        reduce_tree_cotree_dense_system, tree_cotree_gauge_t
    use fortfem_kinds, only: dp
    use fortplot, only: add_scatter, figure, legend, plot, quiver, savefig, &
        title, xlabel, ylabel
    implicit none

    integer, parameter :: edge_count = 5
    integer, parameter :: incidence(4, edge_count) = reshape([ &
        -1, 1, 0, 0, 0, -1, 1, 0, 0, 0, -1, 1, 1, 0, 0, -1, &
        -1, 0, 1, 0], [4, edge_count])
    integer, parameter :: edge_nodes(2, edge_count) = reshape([ &
        1, 2, 2, 3, 3, 4, 4, 1, 1, 3], [2, edge_count])
    real(dp), parameter :: node_x(4) = [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
    real(dp), parameter :: node_y(4) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: mesh_color(3) = [0.65_dp, 0.65_dp, 0.65_dp]
    real(dp), parameter :: node_color(3) = [0.05_dp, 0.05_dp, 0.05_dp]
    real(dp), parameter :: arrow_color(3) = [0.12_dp, 0.35_dp, 0.75_dp]
    real(dp), parameter :: edge_solution(edge_count) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.75_dp, -0.4_dp]
    real(dp), parameter :: full_matrix(edge_count, edge_count) = reshape([ &
        4.0_dp, 0.1_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        0.1_dp, 5.0_dp, 0.2_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.2_dp, 6.0_dp, 0.0_dp, 0.1_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, 7.0_dp, 0.3_dp, &
        0.0_dp, 0.0_dp, 0.1_dp, 0.3_dp, 8.0_dp], &
        [edge_count, edge_count])
    character(*), parameter :: output_directory = &
        "output/example/biro_tree_cotree_benchmark"
    type(tree_cotree_gauge_t) :: gauge
    real(dp), allocatable :: restricted(:), prolonged(:)
    real(dp), allocatable :: reduced_matrix(:, :), reduced_rhs(:)
    real(dp) :: full_rhs(edge_count), residual
    real(dp) :: midpoint_x(edge_count), midpoint_y(edge_count)
    real(dp) :: arrow_u(edge_count), arrow_v(edge_count)
    integer :: edge, first_node, second_node, status, unit

    call execute_command_line("mkdir -p "//output_directory)
    full_rhs = matmul(full_matrix, edge_solution)
    call build_tree_cotree_gauge(incidence, gauge, status)
    if (status /= 0) error stop "Biro tree-cotree gauge construction failed"
    call apply_tree_cotree_restriction(gauge, edge_solution, restricted, status)
    if (status /= 0) error stop "Biro tree-cotree restriction failed"
    call apply_tree_cotree_prolongation(gauge, restricted, prolonged, status)
    if (status /= 0) error stop "Biro tree-cotree prolongation failed"
    call reduce_tree_cotree_dense_system( &
        gauge, full_matrix, full_rhs, reduced_matrix, reduced_rhs, status)
    if (status /= 0) error stop "Biro tree-cotree direct reduction failed"
    residual = maxval(abs(matmul(full_matrix, prolonged) - full_rhs))
    if (residual > 1.0e-12_dp) error stop "Biro manufactured residual failed"

    open (newunit=unit, file=output_directory//"/solution.csv", &
        status="replace", action="write")
    write (unit, "(a)") "edge,x0,y0,x1,y1,A_edge,loop_current_proxy"
    do edge = 1, edge_count
        first_node = edge_nodes(1, edge)
        second_node = edge_nodes(2, edge)
        midpoint_x(edge) = 0.5_dp*(node_x(first_node) + node_x(second_node))
        midpoint_y(edge) = 0.5_dp*(node_y(first_node) + node_y(second_node))
        arrow_u(edge) = 0.18_dp*edge_solution(edge)* &
            (node_x(second_node) - node_x(first_node))
        arrow_v(edge) = 0.18_dp*edge_solution(edge)* &
            (node_y(second_node) - node_y(first_node))
        write (unit, "(*(es24.16,:,','))") real(edge, dp), &
            node_x(first_node), node_y(first_node), node_x(second_node), &
            node_y(second_node), prolonged(edge), dot_product( &
            real(incidence(:, edge), dp), [0.0_dp, 0.0_dp, 1.0_dp, -1.0_dp])
    end do
    close (unit)
    open (newunit=unit, file=output_directory//"/diagnostics.csv", &
        status="replace", action="write")
    write (unit, "(a,es24.16)") "manufactured_residual,", residual
    write (unit, "(a,es24.16)") "direct_cotree_energy,", &
        0.5_dp*dot_product(restricted, matmul(reduced_matrix, restricted))
    close (unit)

    ! Physical solution first: edge potential and its oriented current proxy.
    call figure(figsize=[7.2_dp, 6.2_dp])
    do edge = 1, edge_count
        first_node = edge_nodes(1, edge)
        second_node = edge_nodes(2, edge)
        call plot([node_x(first_node), node_x(second_node)], &
            [node_y(first_node), node_y(second_node)], color=mesh_color, &
            linewidth=1.5_dp)
    end do
    call add_scatter(node_x, node_y, marker="o", color=node_color, &
        label="mesh vertices")
    call quiver(midpoint_x, midpoint_y, arrow_u, arrow_v, scale=1.0_dp, &
        scale_units="xy", angles="xy", color=arrow_color, width=0.004_dp, &
        headwidth=3.0_dp)
    call xlabel("x")
    call ylabel("y")
    call title("Bíró tree--cotree curl--curl solution")
    call legend()
    call savefig(output_directory//"/solution.png")
end program biro_tree_cotree_benchmark
