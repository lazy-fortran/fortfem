program adaptive_bem_prolate
    use fortfem_api, only: &
        estimate_laplace_p0_two_level_residual_3d, &
        generate_sphere_surface_mesh, mark_bem_dorfler, &
        refine_surface_mesh_marked, solve_laplace_dirichlet_p0_3d
    use fortfem_kinds, only: dp
    use fortplot, only: add_scatter, figure, legend, plot, savefig, title, &
        xlabel, ylabel, yscale
    implicit none

    integer, parameter :: step_count = 5
    character(*), parameter :: output_directory = &
        "output/example/adaptive_bem_prolate"
    real(dp), parameter :: major_axis = 2.0_dp, minor_axis = 1.0_dp
    integer, allocatable :: parent(:), refined_triangles(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: density(:), indicators(:)
    real(dp), allocatable :: refined_vertices(:, :), vertices(:, :)
    real(dp), allocatable :: panel_centers(:, :)
    logical, allocatable :: marked(:)
    real(dp) :: capacity, exact_capacity, error(step_count)
    real(dp) :: estimator(step_count), panels(step_count)
    integer :: status, step

    call execute_command_line("mkdir -p "//output_directory)
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    vertices(1:2, :) = minor_axis*vertices(1:2, :)
    vertices(3, :) = major_axis*vertices(3, :)
    exact_capacity = 4.0_dp*acos(-1.0_dp)* &
        sqrt(major_axis**2 - minor_axis**2)/ &
        acosh(major_axis/minor_axis)

    do step = 1, step_count
        call solve_laplace_dirichlet_p0_3d( &
            vertices, triangles, 1.0_dp, 8, density, capacity, status)
        if (status /= 0) error stop "adaptive BEM solve failed"
        call estimate_laplace_p0_two_level_residual_3d( &
            vertices, triangles, density, 1.0_dp, 6, indicators, status)
        if (status /= 0) error stop "adaptive BEM estimate failed"
        panels(step) = real(size(triangles, 2), dp)
        error(step) = abs(capacity - exact_capacity)/exact_capacity
        estimator(step) = norm2(indicators)
        print "(i4,2es14.5)", size(triangles, 2), error(step), estimator(step)
        if (step == step_count) exit
        call mark_bem_dorfler(indicators, 0.5_dp, marked, status)
        if (status /= 0) error stop "adaptive BEM marking failed"
        call refine_surface_mesh_marked( &
            vertices, triangles, marked, refined_vertices, &
            refined_triangles, parent, status)
        if (status /= 0) error stop "adaptive BEM refinement failed"
        call project_to_prolate(refined_vertices)
        call move_alloc(refined_vertices, vertices)
        call move_alloc(refined_triangles, triangles)
    end do

    allocate(panel_centers(3, size(triangles, 2)))
    do step = 1, size(triangles, 2)
        panel_centers(:, step) = sum( &
            vertices(:, triangles(:, step)), dim=2)/3.0_dp
    end do
    call figure(figsize=[7.5_dp, 6.0_dp])
    call add_scatter( &
        panel_centers(1, :), panel_centers(2, :), panel_centers(3, :), &
        c=density, &
        cmap="viridis", marker=".", markersize=5.0_dp, &
        label="P0 surface density")
    call title("Adaptive Laplace BEM density on a prolate spheroid")
    call savefig(output_directory//"/prolate_density_3d.png")

    call figure()
    call plot(panels, error, label="capacity error")
    call plot(panels, estimator, label="residual estimator")
    call yscale("log")
    call xlabel("surface panels")
    call ylabel("relative error / estimator")
    call title("Adaptive P0 BEM on a prolate spheroid")
    call legend()
    call savefig(output_directory//"/adaptive_convergence.png")

contains

    subroutine project_to_prolate(points)
        real(dp), intent(inout) :: points(:, :)
        real(dp) :: scale
        integer :: point

        do point = 1, size(points, 2)
            scale = sqrt( &
                sum((points(1:2, point)/minor_axis)**2) + &
                (points(3, point)/major_axis)**2)
            points(:, point) = points(:, point)/scale
        end do
    end subroutine project_to_prolate

end program adaptive_bem_prolate
