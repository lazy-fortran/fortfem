program biro_tree_cotree_3d_gallery
    !! Three-dimensional manufactured edge-potential tree--cotree gallery.
    !!
    !! The cube is a deliberately small topology/metric analogue of the
    !! direct curl--curl gauge in Biro, Preis, and Richter.  It is not a
    !! redistribution of their application geometry or source data.  Edge
    !! unknowns are line integrals of a smooth manufactured vector potential;
    !! the singular C^T C system is reduced on a spanning tree and solved on
    !! the cotree.  The first figure shows the reconstructed edge solution,
    !! rather than only a residual or convergence curve.
    use fortfem_feec, only: apply_tree_cotree_prolongation, &
        build_tree_cotree_gauge, reduce_tree_cotree_dense_system, &
        tree_cotree_gauge_edges, tree_cotree_gauge_t
    use fortfem_kinds, only: dp
    use fortplot, only: add_3d_plot, add_scatter, colorbar, figure, legend, &
        plot, savefig, title, view_init, xlabel, ylabel
    implicit none

    integer, parameter :: vertex_count = 8, edge_count = 12
    integer, parameter :: face_count = 6
    integer, parameter :: incidence(vertex_count, edge_count) = reshape([ &
        -1, 1, 0, 0, 0, 0, 0, 0, &
        0, -1, 1, 0, 0, 0, 0, 0, &
        0, 0, -1, 1, 0, 0, 0, 0, &
        1, 0, 0, -1, 0, 0, 0, 0, &
        0, 0, 0, 0, -1, 1, 0, 0, &
        0, 0, 0, 0, 0, -1, 1, 0, &
        0, 0, 0, 0, 0, 0, -1, 1, &
        0, 0, 0, 0, 1, 0, 0, -1, &
        -1, 0, 0, 0, 1, 0, 0, 0, &
        0, -1, 0, 0, 0, 1, 0, 0, &
        0, 0, -1, 0, 0, 0, 1, 0, &
        0, 0, 0, -1, 0, 0, 0, 1], [vertex_count, edge_count])
    integer, parameter :: edge_nodes(2, edge_count) = reshape([ &
        1, 2, 2, 3, 3, 4, 4, 1, 5, 6, 6, 7, 7, 8, 8, 5, &
        1, 5, 2, 6, 3, 7, 4, 8], [2, edge_count])
    integer, parameter :: face_edge(face_count, edge_count) = reshape([ &
        1, 0, 1, 0, 0, 0, &
        1, 0, 0, 1, 0, 0, &
        1, 0, 0, 0, 1, 0, &
        1, 0, 0, 0, 0, 1, &
        0, 1, -1, 0, 0, 0, &
        0, 1, 0, -1, 0, 0, &
        0, 1, 0, 0, -1, 0, &
        0, 1, 0, 0, 0, -1, &
        0, 0, -1, 0, 0, 1, &
        0, 0, 1, -1, 0, 0, &
        0, 0, 0, 1, -1, 0, &
        0, 0, 0, 0, 1, -1], &
        [face_count, edge_count])
    real(dp), parameter :: vertex_x(vertex_count) = &
        [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
    real(dp), parameter :: vertex_y(vertex_count) = &
        [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: vertex_z(vertex_count) = &
        [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: tree_color(3) = [0.10_dp, 0.35_dp, 0.76_dp]
    real(dp), parameter :: cotree_color(3) = [0.86_dp, 0.26_dp, 0.14_dp]
    real(dp), parameter :: node_color(3) = [0.08_dp, 0.08_dp, 0.08_dp]
    character(*), parameter :: output_directory = &
        "output/example/biro_tree_cotree_3d_gallery"

    type(tree_cotree_gauge_t) :: gauge
    integer, allocatable :: tree_edges(:), cotree_edges(:)
    real(dp) :: curl_incidence(face_count, edge_count)
    real(dp) :: full_matrix(edge_count, edge_count)
    real(dp) :: edge_target(edge_count), full_rhs(edge_count)
    real(dp) :: face_circulation(face_count), residual
    real(dp), allocatable :: reduced_matrix(:, :), reduced_rhs(:)
    real(dp), allocatable :: reduced_solution(:), edge_solution(:)
    real(dp) :: midpoint_x(edge_count), midpoint_y(edge_count)
    real(dp) :: midpoint_z(edge_count), arrow_x(edge_count)
    real(dp) :: arrow_y(edge_count), arrow_z(edge_count)
    real(dp) :: magnitude(edge_count), edge_length, midpoint(3), tangent(3)
    real(dp) :: vector_potential(3), scale
    integer :: edge, first_node, second_node, face, status, unit
    integer :: command_status

    command_status = 0
    call execute_command_line("mkdir -p "//output_directory, &
        exitstat=command_status)
    if (command_status /= 0) error stop "cannot create Biro 3-D output directory"

    curl_incidence = real(face_edge, dp)
    full_matrix = matmul(transpose(curl_incidence), curl_incidence)
    do edge = 1, edge_count
        first_node = edge_nodes(1, edge)
        second_node = edge_nodes(2, edge)
        midpoint = 0.5_dp*[ &
            vertex_x(first_node) + vertex_x(second_node), &
            vertex_y(first_node) + vertex_y(second_node), &
            vertex_z(first_node) + vertex_z(second_node)]
        tangent = [vertex_x(second_node) - vertex_x(first_node), &
            vertex_y(second_node) - vertex_y(first_node), &
            vertex_z(second_node) - vertex_z(first_node)]
        vector_potential = manufactured_vector_potential(midpoint)
        edge_target(edge) = dot_product(vector_potential, tangent)
    end do
    full_rhs = matmul(full_matrix, edge_target)

    call build_tree_cotree_gauge(incidence, gauge, status)
    if (status /= 0) error stop "Biro 3-D tree-cotree gauge construction failed"
    call tree_cotree_gauge_edges(gauge, tree_edges, cotree_edges, status)
    if (status /= 0 .or. size(tree_edges) /= 7 .or. size(cotree_edges) /= 5) &
        error stop "Biro 3-D cube topology has an unexpected cycle rank"
    call reduce_tree_cotree_dense_system( &
        gauge, full_matrix, full_rhs, reduced_matrix, reduced_rhs, status)
    if (status /= 0) error stop "Biro 3-D reduced system assembly failed"
    call solve_dense(reduced_matrix, reduced_rhs, reduced_solution, status)
    if (status /= 0) error stop "Biro 3-D cotree direct solve failed"
    call apply_tree_cotree_prolongation( &
        gauge, reduced_solution, edge_solution, status)
    if (status /= 0) error stop "Biro 3-D tree-cotree prolongation failed"
    face_circulation = matmul(curl_incidence, edge_solution)
    residual = maxval(abs(matmul(full_matrix, edge_solution) - full_rhs))
    if (residual > 5.0e-12_dp) error stop "Biro 3-D manufactured residual failed"

    do edge = 1, edge_count
        first_node = edge_nodes(1, edge)
        second_node = edge_nodes(2, edge)
        midpoint_x(edge) = 0.5_dp*(vertex_x(first_node) + vertex_x(second_node))
        midpoint_y(edge) = 0.5_dp*(vertex_y(first_node) + vertex_y(second_node))
        midpoint_z(edge) = 0.5_dp*(vertex_z(first_node) + vertex_z(second_node))
        tangent = [vertex_x(second_node) - vertex_x(first_node), &
            vertex_y(second_node) - vertex_y(first_node), &
            vertex_z(second_node) - vertex_z(first_node)]
        edge_length = sqrt(dot_product(tangent, tangent))
        ! Keep the solved edge-potential direction legible in the first
        ! physical view; the line length is a display scale, not a change to
        ! the edge integral stored in ``solution.csv``.
        scale = 0.68_dp*edge_solution(edge)/max(edge_length, 1.0e-12_dp)
        arrow_x(edge) = scale*tangent(1)
        arrow_y(edge) = scale*tangent(2)
        arrow_z(edge) = scale*tangent(3)
        magnitude(edge) = abs(edge_solution(edge))
    end do

    open (newunit=unit, file=output_directory//"/solution.csv", &
        status="replace", action="write")
    write (unit, "(a)") "edge,x0,y0,z0,x1,y1,z1,A_edge,tree_cotree"
    do edge = 1, edge_count
        first_node = edge_nodes(1, edge)
        second_node = edge_nodes(2, edge)
        write (unit, "(*(es24.16,:,','))") real(edge, dp), &
            vertex_x(first_node), vertex_y(first_node), vertex_z(first_node), &
            vertex_x(second_node), vertex_y(second_node), vertex_z(second_node), &
            edge_solution(edge), real(edge_class(edge, tree_edges), dp)
    end do
    close (unit)
    open (newunit=unit, file=output_directory//"/diagnostics.csv", &
        status="replace", action="write")
    write (unit, "(a,es24.16)") "manufactured_residual,", residual
    do face = 1, face_count
        write (unit, "(a,i0,a,es24.16)") "face_circulation_", face, ",", &
            face_circulation(face)
    end do
    write (unit, "(a,es24.16)") "cotree_energy,", &
        0.5_dp*dot_product(reduced_solution, &
            matmul(reduced_matrix, reduced_solution))
    close (unit)

    ! Physical solution first: 3-D edge-element potential on the cube.
    call figure(figsize=[8.5_dp, 7.2_dp])
    do edge = 1, edge_count
        first_node = edge_nodes(1, edge)
        second_node = edge_nodes(2, edge)
        if (edge_class(edge, tree_edges) == 1) then
            if (edge == tree_edges(1)) then
                call add_3d_plot( &
                    [vertex_x(first_node), vertex_x(second_node)], &
                    [vertex_y(first_node), vertex_y(second_node)], &
                    [vertex_z(first_node), vertex_z(second_node)], &
                    color=tree_color, linewidth=2.2_dp, label="tree edges")
            else
                call add_3d_plot( &
                    [vertex_x(first_node), vertex_x(second_node)], &
                    [vertex_y(first_node), vertex_y(second_node)], &
                    [vertex_z(first_node), vertex_z(second_node)], &
                    color=tree_color, linewidth=2.2_dp)
            end if
        else
            if (edge == cotree_edges(1)) then
                call add_3d_plot( &
                    [vertex_x(first_node), vertex_x(second_node)], &
                    [vertex_y(first_node), vertex_y(second_node)], &
                    [vertex_z(first_node), vertex_z(second_node)], &
                    color=cotree_color, linewidth=2.8_dp, label="cotree edges")
            else
                call add_3d_plot( &
                    [vertex_x(first_node), vertex_x(second_node)], &
                    [vertex_y(first_node), vertex_y(second_node)], &
                    [vertex_z(first_node), vertex_z(second_node)], &
                    color=cotree_color, linewidth=2.8_dp)
            end if
        end if
        call add_3d_plot( &
            [midpoint_x(edge), midpoint_x(edge) + arrow_x(edge)], &
            [midpoint_y(edge), midpoint_y(edge) + arrow_y(edge)], &
            [midpoint_z(edge), midpoint_z(edge) + arrow_z(edge)], &
            color="black", linewidth=1.5_dp)
    end do
    call add_scatter(vertex_x, vertex_y, vertex_z, marker="o", color=node_color, &
        label="cube vertices")
    call add_scatter(midpoint_x, midpoint_y, midpoint_z, c=magnitude, &
        cmap="viridis", marker=".", markersize=7.0_dp, label="edge |A| DOFs")
    call colorbar(label="|edge potential|")
    call xlabel("x")
    call ylabel("y")
    call title("Bíró tree--cotree curl--curl solution on a 3-D cube")
    call legend()
    call view_init(elev=24.0_dp, azim=-48.0_dp)
    call savefig(output_directory//"/solution.png")
    call savefig(output_directory//"/solution_3d.png")

    call figure(figsize=[8.0_dp, 5.2_dp])
    call plot([(real(face, dp), face=1, face_count)], face_circulation, &
        marker="o", linewidth=2.0_dp, color=cotree_color)
    call xlabel("oriented cube face")
    call ylabel("discrete curl circulation")
    call title("Bíró tree--cotree face-circulation diagnostic")
    call savefig(output_directory//"/diagnostics_1d.png")

contains

    pure function manufactured_vector_potential(point) result(value)
        real(dp), intent(in) :: point(3)
        real(dp) :: value(3)

        value = [-0.30_dp*point(2) + 0.10_dp*point(3), &
            0.22_dp*point(1) + 0.05_dp*point(3), &
            0.16_dp*point(2) - 0.08_dp*point(1)]
    end function manufactured_vector_potential

    integer function edge_class(edge, tree) result(class)
        integer, intent(in) :: edge
        integer, intent(in) :: tree(:)

        if (any(tree == edge)) then
            class = 1
        else
            class = 2
        end if
    end function edge_class

    subroutine solve_dense(matrix, rhs, solution, status)
        real(dp), intent(in) :: matrix(:, :), rhs(:)
        real(dp), allocatable, intent(out) :: solution(:)
        integer, intent(out) :: status
        real(dp), allocatable :: work(:, :), work_rhs(:), row_buffer(:)
        real(dp) :: pivot, factor
        integer :: n, row, column, pivot_row

        status = 1
        allocate(solution(0))
        n = size(rhs)
        if (size(matrix, 1) /= n .or. size(matrix, 2) /= n .or. n < 1) return
        allocate(work(n, n), work_rhs(n), row_buffer(n))
        work = matrix
        work_rhs = rhs
        do column = 1, n - 1
            pivot_row = column
            do row = column + 1, n
                if (abs(work(row, column)) > abs(work(pivot_row, column))) &
                    pivot_row = row
            end do
            pivot = work(pivot_row, column)
            if (abs(pivot) < 1.0e-13_dp) return
            if (pivot_row /= column) then
                row_buffer = work(column, :)
                work(column, :) = work(pivot_row, :)
                work(pivot_row, :) = row_buffer
                factor = work_rhs(column)
                work_rhs(column) = work_rhs(pivot_row)
                work_rhs(pivot_row) = factor
            end if
            do row = column + 1, n
                factor = work(row, column)/work(column, column)
                work(row, column:n) = work(row, column:n) - &
                    factor*work(column, column:n)
                work_rhs(row) = work_rhs(row) - factor*work_rhs(column)
            end do
        end do
        if (abs(work(n, n)) < 1.0e-13_dp) return
        deallocate(solution)
        allocate(solution(n))
        solution(n) = work_rhs(n)/work(n, n)
        do row = n - 1, 1, -1
            solution(row) = (work_rhs(row) - &
                dot_product(work(row, row + 1:n), solution(row + 1:n)))/ &
                work(row, row)
        end do
        status = 0
    end subroutine solve_dense

end program biro_tree_cotree_3d_gallery
