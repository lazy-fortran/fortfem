program bem_sphere_3d
    use fortfem_api, only: &
        apply_laplace_single_layer_p0_hierarchical_3d, &
        assemble_laplace_single_layer_p0_3d, &
        generate_sphere_surface_mesh, &
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
        call generate_sphere_surface_mesh(1.0_dp, level, vertices, triangles)
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

end program bem_sphere_3d
