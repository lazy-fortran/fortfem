program bem_sphere_3d
    use fortfem_api, only: &
        apply_laplace_single_layer_p0_hierarchical_3d, &
        assemble_laplace_single_layer_p0_3d, &
        solve_laplace_dirichlet_p0_3d
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, title, xlabel, ylabel
    implicit none

    character(*), parameter :: output_directory = &
        "output/example/bem_sphere_3d"
    real(dp), allocatable :: dense_action(:), density(:), fast_action(:)
    real(dp), allocatable :: matrix(:, :), vertices(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp) :: capacities(0:2), dense_seconds, exact(0:2)
    real(dp) :: fast_error, fast_seconds, panel_axis(0:2), seconds(0:2)
    integer :: interaction_count, level, status, unit

    call execute_command_line("mkdir -p "//output_directory)
    do level = 0, 2
        call build_sphere(level, vertices, triangles)
        panel_axis(level) = real(size(triangles, 2), dp)
        call cpu_time(dense_seconds)
        call solve_laplace_dirichlet_p0_3d( &
            vertices, triangles, 1.0_dp, 8, density, capacities(level), status)
        call cpu_time(seconds(level))
        seconds(level) = seconds(level) - dense_seconds
        if (status /= 0) error stop "sphere Galerkin solve failed"
    end do
    exact = 4.0_dp*acos(-1.0_dp)

    call assemble_laplace_single_layer_p0_3d( &
        vertices, triangles, 8, matrix, status)
    if (status /= 0) error stop "dense single-layer assembly failed"
    dense_action = matmul(matrix, density)
    call cpu_time(fast_seconds)
    call apply_laplace_single_layer_p0_hierarchical_3d( &
        vertices, triangles, density, 0.6_dp, 6, fast_action, status, &
        interaction_count)
    call cpu_time(dense_seconds)
    fast_seconds = dense_seconds - fast_seconds
    if (status /= 0) error stop "hierarchical single-layer apply failed"
    fast_error = norm2(fast_action - dense_action)/norm2(dense_action)

    call figure(figsize=[8.0_dp, 5.0_dp])
    call plot(panel_axis, exact, label="analytical 4 pi", linestyle="-")
    call plot( &
        panel_axis, capacities, label="P0 Galerkin BEM", linestyle="--", &
        marker="o")
    call xlabel("surface triangles")
    call ylabel("capacitance")
    call title("Three-dimensional sphere capacitance")
    call legend()
    call savefig(output_directory//"/sphere_capacitance.png")

    call figure(figsize=[8.0_dp, 5.0_dp])
    call plot( &
        [(real(level, dp), level=1, size(dense_action))], dense_action, &
        label="dense Galerkin", linestyle="-")
    call plot( &
        [(real(level, dp), level=1, size(fast_action))], fast_action, &
        label="hierarchical", linestyle="--")
    call xlabel("panel index")
    call ylabel("single-layer action")
    call title("Dense versus hierarchical 3D BEM")
    call legend()
    call savefig(output_directory//"/sphere_hierarchical_action.png")

    open( &
        newunit=unit, file=output_directory//"/benchmark.txt", &
        status="replace", action="write")
    do level = 0, 2
        write (unit, '(A,I0,A,ES14.6,A,ES14.6)') &
            "panels=", nint(panel_axis(level)), &
            " capacity=", capacities(level), " solve_seconds=", seconds(level)
    end do
    write (unit, '(A,ES14.6)') "hierarchical_seconds=", fast_seconds
    write (unit, '(A,ES14.6)') "hierarchical_relative_error=", fast_error
    write (unit, '(A,I0)') "hierarchical_interactions=", interaction_count
    write (unit, '(A,I0)') "dense_interactions=", size(triangles, 2)**2
    close(unit)

contains

    subroutine build_sphere(level, vertices, triangles)
        integer, intent(in) :: level
        real(dp), allocatable, intent(out) :: vertices(:, :)
        integer, allocatable, intent(out) :: triangles(:, :)

        integer :: refinement

        allocate(vertices(3, 6), triangles(3, 8))
        vertices(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
        vertices(:, 2) = [-1.0_dp, 0.0_dp, 0.0_dp]
        vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
        vertices(:, 4) = [0.0_dp, -1.0_dp, 0.0_dp]
        vertices(:, 5) = [0.0_dp, 0.0_dp, 1.0_dp]
        vertices(:, 6) = [0.0_dp, 0.0_dp, -1.0_dp]
        triangles(:, 1) = [1, 3, 5]
        triangles(:, 2) = [3, 2, 5]
        triangles(:, 3) = [2, 4, 5]
        triangles(:, 4) = [4, 1, 5]
        triangles(:, 5) = [3, 1, 6]
        triangles(:, 6) = [2, 3, 6]
        triangles(:, 7) = [4, 2, 6]
        triangles(:, 8) = [1, 4, 6]
        do refinement = 1, level
            call refine_sphere(vertices, triangles)
        end do
    end subroutine build_sphere

    subroutine refine_sphere(vertices, triangles)
        real(dp), allocatable, intent(inout) :: vertices(:, :)
        integer, allocatable, intent(inout) :: triangles(:, :)

        real(dp), allocatable :: expanded_vertices(:, :)
        integer, allocatable :: edge_midpoints(:), edge_vertices(:, :)
        integer, allocatable :: refined(:, :)
        integer :: edge_count, midpoint(3), old_vertex_count
        integer :: triangle, vertex_count

        old_vertex_count = size(vertices, 2)
        allocate(expanded_vertices(3, old_vertex_count + 3*size(triangles, 2)))
        expanded_vertices(:, :old_vertex_count) = vertices
        allocate(edge_vertices(2, 3*size(triangles, 2)))
        allocate(edge_midpoints(3*size(triangles, 2)))
        allocate(refined(3, 4*size(triangles, 2)))
        vertex_count = old_vertex_count
        edge_count = 0
        do triangle = 1, size(triangles, 2)
            midpoint(1) = midpoint_vertex( &
                triangles(1, triangle), triangles(2, triangle), vertices, &
                expanded_vertices, edge_vertices, edge_midpoints, &
                edge_count, vertex_count)
            midpoint(2) = midpoint_vertex( &
                triangles(2, triangle), triangles(3, triangle), vertices, &
                expanded_vertices, edge_vertices, edge_midpoints, &
                edge_count, vertex_count)
            midpoint(3) = midpoint_vertex( &
                triangles(3, triangle), triangles(1, triangle), vertices, &
                expanded_vertices, edge_vertices, edge_midpoints, &
                edge_count, vertex_count)
            refined(:, 4*triangle - 3) = &
                [triangles(1, triangle), midpoint(1), midpoint(3)]
            refined(:, 4*triangle - 2) = &
                [midpoint(1), triangles(2, triangle), midpoint(2)]
            refined(:, 4*triangle - 1) = &
                [midpoint(3), midpoint(2), triangles(3, triangle)]
            refined(:, 4*triangle) = midpoint
        end do
        vertices = expanded_vertices(:, :vertex_count)
        call move_alloc(refined, triangles)
    end subroutine refine_sphere

    function midpoint_vertex( &
            first, second, vertices, expanded_vertices, edge_vertices, &
            edge_midpoints, edge_count, vertex_count) result(midpoint)
        integer, intent(in) :: first, second
        real(dp), intent(in) :: vertices(:, :)
        real(dp), intent(inout) :: expanded_vertices(:, :)
        integer, intent(inout) :: edge_vertices(:, :), edge_midpoints(:)
        integer, intent(inout) :: edge_count, vertex_count
        integer :: midpoint

        integer :: edge, ordered(2)
        real(dp) :: midpoint_coordinate(3)

        ordered = [min(first, second), max(first, second)]
        do edge = 1, edge_count
            if (all(edge_vertices(:, edge) == ordered)) then
                midpoint = edge_midpoints(edge)
                return
            end if
        end do
        edge_count = edge_count + 1
        edge_vertices(:, edge_count) = ordered
        vertex_count = vertex_count + 1
        midpoint_coordinate = vertices(:, first) + vertices(:, second)
        expanded_vertices(:, vertex_count) = &
            midpoint_coordinate/norm2(midpoint_coordinate)
        edge_midpoints(edge_count) = vertex_count
        midpoint = vertex_count
    end function midpoint_vertex

end program bem_sphere_3d
