program bem_sphere_3d
    use fortfem_api, only: &
        apply_helmholtz_cfie_p0_hierarchical_3d, &
        apply_helmholtz_single_layer_p0_hierarchical_3d, &
        apply_laplace_single_layer_p0_hierarchical_3d, &
        assemble_helmholtz_single_layer_p0_3d, &
        assemble_laplace_single_layer_p0_3d, &
        evaluate_helmholtz_cfie_p0_3d, &
        evaluate_helmholtz_representation_triangles_3d, &
        generate_sphere_surface_mesh, &
        solve_helmholtz_cfie_p0_3d, &
        solve_helmholtz_cfie_p0_hierarchical_3d, &
        solve_helmholtz_dirichlet_p0_3d, &
        solve_helmholtz_dirichlet_p0_hierarchical_3d, &
        solve_laplace_dirichlet_p0_3d
    use fortfem_kinds, only: dp
    use fortplot, only: add_scatter, figure, legend, plot, savefig, title, &
        xlabel, xscale, ylabel, yscale
    implicit none

    character(*), parameter :: output_directory = &
        "output/example/bem_sphere_3d"
    real(dp), allocatable :: dense_action(:), density(:), fast_action(:)
    real(dp), allocatable :: matrix(:, :), vertices(:, :)
    real(dp), allocatable :: panel_centers(:, :)
    complex(dp), allocatable :: complex_density(:), helmholtz_dense(:)
    complex(dp), allocatable :: helmholtz_dense_density(:)
    complex(dp), allocatable :: helmholtz_dirichlet(:)
    complex(dp), allocatable :: helmholtz_fast(:), helmholtz_matrix(:, :)
    complex(dp), allocatable :: helmholtz_fast_density(:)
    complex(dp), allocatable :: cfie_dense_density(:), cfie_fast_density(:)
    complex(dp), allocatable :: scaling_action(:), scaling_density(:)
    integer, allocatable :: triangles(:, :)
    integer, allocatable :: scaling_triangles(:, :)
    real(dp), allocatable :: scaling_vertices(:, :)
    complex(dp) :: analytical_field(16), dense_field(16), fast_field(16)
    complex(dp) :: cfie_analytical_field(16), cfie_dense_field(16)
    complex(dp) :: cfie_fast_field(16)
    real(dp) :: capacities(0:2), dense_seconds, exact(0:2)
    real(dp) :: field_radius(16)
    real(dp) :: fast_error, fast_seconds, panel_axis(0:2), seconds(0:2)
    real(dp) :: helmholtz_dense_field_error, helmholtz_density_error
    real(dp) :: helmholtz_dense_seconds, helmholtz_dense_solve_seconds
    real(dp) :: helmholtz_error, helmholtz_fast_field_error
    real(dp) :: helmholtz_fast_seconds, helmholtz_fast_solve_seconds
    real(dp) :: helmholtz_residual, start_time
    real(dp) :: cfie_dense_field_error, cfie_dense_solve_seconds
    real(dp) :: cfie_fast_field_error, cfie_fast_solve_seconds, cfie_residual
    integer :: cfie_interactions, cfie_iterations
    integer :: scaling_interactions(3)
    real(dp) :: scaling_dense_interactions(3), scaling_panels(3)
    real(dp) :: scaling_seconds(3)
    integer :: field_point, helmholtz_interactions, laplace_interactions
    integer :: level, solve_interactions, solve_iterations, status, unit

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
        laplace_interactions)
    call cpu_time(dense_seconds)
    fast_seconds = dense_seconds - fast_seconds
    if (status /= 0) error stop "hierarchical single-layer apply failed"
    fast_error = norm2(fast_action - dense_action)/norm2(dense_action)

    complex_density = cmplx(density, 0.25_dp*density, dp)
    call assemble_helmholtz_single_layer_p0_3d( &
        vertices, triangles, 0.7_dp, 8, helmholtz_matrix, status)
    if (status /= 0) error stop "dense Helmholtz assembly failed"
    call cpu_time(start_time)
    helmholtz_dense = matmul(helmholtz_matrix, complex_density)
    call cpu_time(helmholtz_dense_seconds)
    helmholtz_dense_seconds = helmholtz_dense_seconds - start_time
    call cpu_time(start_time)
    call apply_helmholtz_single_layer_p0_hierarchical_3d( &
        vertices, triangles, complex_density, 0.7_dp, 0.45_dp, 6, &
        helmholtz_fast, status, helmholtz_interactions)
    call cpu_time(helmholtz_fast_seconds)
    helmholtz_fast_seconds = helmholtz_fast_seconds - start_time
    if (status /= 0) error stop "hierarchical Helmholtz apply failed"
    helmholtz_error = norm2(abs(helmholtz_fast - helmholtz_dense))/ &
        norm2(abs(helmholtz_dense))

    call cpu_time(start_time)
    call solve_helmholtz_dirichlet_p0_3d( &
        vertices, triangles, cmplx(1.0_dp, 0.0_dp, dp), 0.7_dp, 8, &
        helmholtz_dense_density, status)
    call cpu_time(helmholtz_dense_solve_seconds)
    helmholtz_dense_solve_seconds = &
        helmholtz_dense_solve_seconds - start_time
    if (status /= 0) error stop "dense Helmholtz solve failed"
    call cpu_time(start_time)
    call solve_helmholtz_dirichlet_p0_hierarchical_3d( &
        vertices, triangles, cmplx(1.0_dp, 0.0_dp, dp), 0.7_dp, &
        0.3_dp, 6, 1.0e-10_dp, 80, 20, helmholtz_fast_density, status, &
        solve_iterations, helmholtz_residual, solve_interactions)
    call cpu_time(helmholtz_fast_solve_seconds)
    helmholtz_fast_solve_seconds = &
        helmholtz_fast_solve_seconds - start_time
    if (status /= 0) error stop "hierarchical Helmholtz solve failed"
    helmholtz_density_error = &
        norm2(abs(helmholtz_fast_density - helmholtz_dense_density))/ &
        norm2(abs(helmholtz_dense_density))
    allocate(helmholtz_dirichlet(size(vertices, 2)))
    helmholtz_dirichlet = cmplx(0.0_dp, 0.0_dp, dp)
    do field_point = 1, size(field_radius)
        field_radius(field_point) = 1.1_dp + &
            1.9_dp*real(field_point - 1, dp)/real(size(field_radius) - 1, dp)
        analytical_field(field_point) = &
            exp(cmplx(0.0_dp, 0.7_dp*(field_radius(field_point) - 1.0_dp), &
            dp))/field_radius(field_point)
        call evaluate_helmholtz_representation_triangles_3d( &
            vertices, triangles, helmholtz_dirichlet, &
            -helmholtz_dense_density, [0.0_dp, 0.0_dp, field_radius(field_point)], &
            0.7_dp, 8, dense_field(field_point), status)
        if (status /= 0) error stop "dense Helmholtz field evaluation failed"
        call evaluate_helmholtz_representation_triangles_3d( &
            vertices, triangles, helmholtz_dirichlet, &
            -helmholtz_fast_density, [0.0_dp, 0.0_dp, field_radius(field_point)], &
            0.7_dp, 8, fast_field(field_point), status)
        if (status /= 0) then
            error stop "hierarchical Helmholtz field evaluation failed"
        end if
    end do
    helmholtz_dense_field_error = maxval(abs(dense_field - analytical_field))
    helmholtz_fast_field_error = maxval(abs(fast_field - analytical_field))

    call cpu_time(start_time)
    call solve_helmholtz_cfie_p0_3d( &
        vertices, triangles, cmplx(1.0_dp, 0.0_dp, dp), acos(-1.0_dp), &
        acos(-1.0_dp), 8, cfie_dense_density, status)
    call cpu_time(cfie_dense_solve_seconds)
    cfie_dense_solve_seconds = cfie_dense_solve_seconds - start_time
    if (status /= 0) error stop "dense Helmholtz CFIE solve failed"
    call cpu_time(start_time)
    call solve_helmholtz_cfie_p0_hierarchical_3d( &
        vertices, triangles, cmplx(1.0_dp, 0.0_dp, dp), acos(-1.0_dp), &
        acos(-1.0_dp), 0.45_dp, 6, 1.0e-10_dp, 80, 20, &
        cfie_fast_density, status, cfie_iterations, cfie_residual, &
        cfie_interactions)
    call cpu_time(cfie_fast_solve_seconds)
    cfie_fast_solve_seconds = cfie_fast_solve_seconds - start_time
    if (status /= 0) error stop "hierarchical Helmholtz CFIE solve failed"
    do field_point = 1, size(field_radius)
        cfie_analytical_field(field_point) = exp(cmplx( &
            0.0_dp, acos(-1.0_dp)*(field_radius(field_point) - 1.0_dp), dp))/ &
            field_radius(field_point)
        call evaluate_helmholtz_cfie_p0_3d( &
            vertices, triangles, cfie_dense_density, &
            [0.0_dp, 0.0_dp, field_radius(field_point)], acos(-1.0_dp), &
            acos(-1.0_dp), 8, cfie_dense_field(field_point), status)
        if (status /= 0) error stop "dense Helmholtz CFIE evaluation failed"
        call evaluate_helmholtz_cfie_p0_3d( &
            vertices, triangles, cfie_fast_density, &
            [0.0_dp, 0.0_dp, field_radius(field_point)], acos(-1.0_dp), &
            acos(-1.0_dp), 8, cfie_fast_field(field_point), status)
        if (status /= 0) then
            error stop "hierarchical Helmholtz CFIE evaluation failed"
        end if
    end do
    cfie_dense_field_error = &
        maxval(abs(cfie_dense_field - cfie_analytical_field))
    cfie_fast_field_error = &
        maxval(abs(cfie_fast_field - cfie_analytical_field))

    allocate(panel_centers(3, size(triangles, 2)))
    do field_point = 1, size(triangles, 2)
        panel_centers(:, field_point) = sum( &
            vertices(:, triangles(:, field_point)), dim=2)/3.0_dp
    end do
    call figure(figsize=[7.5_dp, 6.0_dp])
    call add_scatter( &
        panel_centers(1, :), panel_centers(2, :), panel_centers(3, :), &
        c=density, &
        cmap="viridis", marker=".", markersize=5.0_dp, &
        label="Laplace P0 surface density")
    call title("Three-dimensional Laplace BEM density on the sphere")
    call savefig(output_directory//"/sphere_density_3d.png")

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

    call figure(figsize=[8.0_dp, 5.0_dp])
    call plot( &
        [(real(level, dp), level=1, size(helmholtz_dense))], &
        real(helmholtz_dense, dp), label="dense real", linestyle="-")
    call plot( &
        [(real(level, dp), level=1, size(helmholtz_fast))], &
        real(helmholtz_fast, dp), label="hierarchical real", linestyle="--")
    call plot( &
        [(real(level, dp), level=1, size(helmholtz_dense))], &
        aimag(helmholtz_dense), label="dense imaginary", linestyle="-")
    call plot( &
        [(real(level, dp), level=1, size(helmholtz_fast))], &
        aimag(helmholtz_fast), label="hierarchical imaginary", linestyle="--")
    call xlabel("panel index")
    call ylabel("Helmholtz single-layer action")
    call title("Dense versus hierarchical outgoing Helmholtz BEM")
    call legend()
    call savefig(output_directory//"/sphere_helmholtz_hierarchical_action.png")

    call figure(figsize=[8.0_dp, 5.0_dp])
    call plot( &
        field_radius, real(analytical_field, dp), label="analytical real", &
        linestyle="-")
    call plot( &
        field_radius, real(dense_field, dp), label="dense BEM real", &
        linestyle="--", marker="o")
    call plot( &
        field_radius, real(fast_field, dp), label="hierarchical BEM real", &
        linestyle=":")
    call plot( &
        field_radius, aimag(analytical_field), label="analytical imaginary", &
        linestyle="-")
    call plot( &
        field_radius, aimag(dense_field), label="dense BEM imaginary", &
        linestyle="--", marker="o")
    call plot( &
        field_radius, aimag(fast_field), &
        label="hierarchical BEM imaginary", linestyle=":")
    call xlabel("radial coordinate")
    call ylabel("outgoing Helmholtz field")
    call title("Analytical versus dense and hierarchical sphere BEM")
    call legend()
    call savefig(output_directory//"/sphere_helmholtz_solve.png")

    call figure(figsize=[8.0_dp, 5.0_dp])
    call plot( &
        field_radius, real(cfie_analytical_field, dp), &
        label="analytical real", linestyle="-")
    call plot( &
        field_radius, real(cfie_dense_field, dp), label="dense CFIE real", &
        linestyle="--", marker="o")
    call plot( &
        field_radius, real(cfie_fast_field, dp), &
        label="hierarchical CFIE real", linestyle=":")
    call plot( &
        field_radius, aimag(cfie_analytical_field), &
        label="analytical imaginary", linestyle="-")
    call plot( &
        field_radius, aimag(cfie_dense_field), &
        label="dense CFIE imaginary", linestyle="--", marker="o")
    call plot( &
        field_radius, aimag(cfie_fast_field), &
        label="hierarchical CFIE imaginary", linestyle=":")
    call xlabel("radial coordinate")
    call ylabel("outgoing Helmholtz field")
    call title("Resonance-safe sphere CFIE at k = pi")
    call legend()
    call savefig(output_directory//"/sphere_helmholtz_cfie_resonance.png")

    do level = 2, 4
        call generate_sphere_surface_mesh( &
            1.0_dp, level, scaling_vertices, scaling_triangles)
        scaling_panels(level - 1) = real(size(scaling_triangles, 2), dp)
        scaling_dense_interactions(level - 1) = scaling_panels(level - 1)**2
        allocate(scaling_density(size(scaling_triangles, 2)))
        scaling_density = cmplx(1.0_dp, 0.25_dp, dp)
        call cpu_time(start_time)
        call apply_helmholtz_cfie_p0_hierarchical_3d( &
            scaling_vertices, scaling_triangles, scaling_density, &
            acos(-1.0_dp), acos(-1.0_dp), 0.45_dp, 6, scaling_action, &
            status, scaling_interactions(level - 1))
        call cpu_time(scaling_seconds(level - 1))
        scaling_seconds(level - 1) = &
            scaling_seconds(level - 1) - start_time
        if (status /= 0) error stop "hierarchical CFIE scaling apply failed"
        deallocate(scaling_density)
    end do
    call figure(figsize=[8.0_dp, 5.0_dp])
    call plot( &
        scaling_panels, scaling_dense_interactions, &
        label="dense pair interactions", linestyle="--", marker="o")
    call plot( &
        scaling_panels, real(scaling_interactions, dp), &
        label="hierarchical interactions", linestyle="-", marker="o")
    call xscale("log")
    call yscale("log")
    call xlabel("surface triangles")
    call ylabel("operator interactions")
    call title("Helmholtz CFIE interaction scaling")
    call legend()
    call savefig(output_directory//"/sphere_helmholtz_cfie_scaling.png")

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
    write (unit, '(A,I0)') &
        "hierarchical_interactions=", laplace_interactions
    write (unit, '(A,I0)') "dense_interactions=", size(triangles, 2)**2
    write (unit, '(A,ES14.6)') &
        "helmholtz_dense_apply_seconds=", helmholtz_dense_seconds
    write (unit, '(A,ES14.6)') &
        "helmholtz_hierarchical_seconds=", helmholtz_fast_seconds
    write (unit, '(A,ES14.6)') &
        "helmholtz_hierarchical_relative_error=", helmholtz_error
    write (unit, '(A,I0)') &
        "helmholtz_hierarchical_interactions=", helmholtz_interactions
    write (unit, '(A,ES14.6)') &
        "helmholtz_dense_solve_seconds=", helmholtz_dense_solve_seconds
    write (unit, '(A,ES14.6)') &
        "helmholtz_hierarchical_solve_seconds=", helmholtz_fast_solve_seconds
    write (unit, '(A,I0)') &
        "helmholtz_hierarchical_solve_iterations=", solve_iterations
    write (unit, '(A,ES14.6)') &
        "helmholtz_hierarchical_solve_residual=", helmholtz_residual
    write (unit, '(A,I0)') &
        "helmholtz_hierarchical_solve_interactions=", solve_interactions
    write (unit, '(A,ES14.6)') &
        "helmholtz_hierarchical_density_error=", helmholtz_density_error
    write (unit, '(A,ES14.6)') &
        "helmholtz_dense_field_max_error=", helmholtz_dense_field_error
    write (unit, '(A,ES14.6)') &
        "helmholtz_hierarchical_field_max_error=", helmholtz_fast_field_error
    write (unit, '(A,ES14.6)') &
        "helmholtz_cfie_dense_solve_seconds=", cfie_dense_solve_seconds
    write (unit, '(A,ES14.6)') &
        "helmholtz_cfie_hierarchical_solve_seconds=", cfie_fast_solve_seconds
    write (unit, '(A,I0)') &
        "helmholtz_cfie_hierarchical_iterations=", cfie_iterations
    write (unit, '(A,ES14.6)') &
        "helmholtz_cfie_hierarchical_residual=", cfie_residual
    write (unit, '(A,I0)') &
        "helmholtz_cfie_hierarchical_interactions=", cfie_interactions
    write (unit, '(A,ES14.6)') &
        "helmholtz_cfie_dense_field_max_error=", cfie_dense_field_error
    write (unit, '(A,ES14.6)') &
        "helmholtz_cfie_hierarchical_field_max_error=", cfie_fast_field_error
    do level = 1, 3
        write (unit, '(A,I0,A,I0,A,ES14.6)') &
            "helmholtz_cfie_scaling_panels=", nint(scaling_panels(level)), &
            " interactions=", scaling_interactions(level), &
            " apply_seconds=", scaling_seconds(level)
    end do
    close(unit)

contains

end program bem_sphere_3d
