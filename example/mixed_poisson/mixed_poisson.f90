program mixed_poisson
    use fortfem_feec, only: &
        build_triangle_discontinuous_dof_map, &
        build_triangle_trimmed_dof_map, div, dx, &
        evaluate_triangle_lagrange_basis, evaluate_triangle_rt_interpolant, &
        form_expr_t, function_space, init_measures, inner, &
        initialize_triangle_lagrange_basis, &
        initialize_triangle_raviart_thomas, operator(*), &
        solve_symbolic_mixed_poisson_rt, test_function, &
        test_function_t, triangle_duffy_quadrature, triangle_lagrange_basis_t, &
        triangle_rt_basis_t, trial_function, trial_function_t, &
        vector_function_space, vector_test_function, &
        vector_test_function_t, vector_trial_function, vector_trial_function_t
    use fortfem_core, only: mesh_t, function_space_t, vector_function_space_t, &
        rectangle_mesh
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortplot, only: colorbar, figure, legend, pcolormesh, plot, savefig, &
        set_xscale, set_yscale, title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    character(*), parameter :: output_directory = &
        "output/example/mixed_poisson"
    integer, parameter :: mesh_sizes(3) = [5, 9, 17]
    real(dp) :: errors(2, 3), flux_errors(3), h(3), pressure_errors(3)
    type(mesh_2d_t) :: plotted_mesh
    real(dp), allocatable :: plotted_pressure_dofs(:)
    integer :: level

    call init_measures()
    call execute_command_line("mkdir -p "//output_directory)
    call initialize_gallery_sequence()
    do level = 1, size(mesh_sizes)
        h(level) = 1.0_dp/real(mesh_sizes(level) - 1, dp)
        if (level == size(mesh_sizes)) then
            call solve_and_measure( &
                mesh_sizes(level), errors(:, level), plotted_pressure_dofs, &
                plotted_mesh)
        else
            call solve_and_measure(mesh_sizes(level), errors(:, level))
        end if
    end do
    if (minval(log(errors(:, :2)/errors(:, 2:))/ &
        spread(log(h(:2)/h(2:)), 1, 2)) <= 1.45_dp) then
        error stop "symbolic mixed Poisson convergence regression"
    end if
    flux_errors = errors(1, :)
    pressure_errors = errors(2, :)

    call render_solution(plotted_mesh, plotted_pressure_dofs)

    call figure(figsize=[8.5_dp, 5.5_dp])
    call plot(h, flux_errors, label="RT1 flux L2 error", marker="o")
    call plot(h, pressure_errors, label="DG1 pressure L2 error", marker="s")
    call set_xscale("log")
    call set_yscale("log")
    call xlabel("mesh spacing h")
    call ylabel("L2 error")
    call title("Symbolic RT-DG mixed Poisson convergence")
    call legend()
    call savefig(output_directory//"/mixed_poisson_convergence_1d.png")
    call record_gallery_stage("diagnostics")

contains

    subroutine render_solution(mesh, pressure_dofs)
        integer, parameter :: nx = 32, ny = 32
        type(mesh_2d_t), intent(inout) :: mesh
        real(dp), intent(in) :: pressure_dofs(:)
        real(dp) :: x_edges(nx + 1), y_edges(ny + 1), values(nx, ny)
        real(dp) :: exact_values(nx, ny), point(2)
        type(triangle_lagrange_basis_t) :: pressure_basis
        integer, allocatable :: pressure_map(:, :)
        integer :: i, j, status, pressure_count, unit

        call initialize_triangle_lagrange_basis(1, pressure_basis, status)
        if (status /= 0) error stop "DG1 basis failed for gallery plot"
        call build_triangle_discontinuous_dof_map( &
            mesh, 1, pressure_map, pressure_count, status)
        if (status /= 0) error stop "DG1 map failed for gallery plot"
        if (size(pressure_dofs) /= pressure_count) &
            error stop "DG1 gallery coefficient size mismatch"
        do i = 1, nx + 1
            x_edges(i) = real(i - 1, dp)/real(nx, dp)
        end do
        do j = 1, ny + 1
            y_edges(j) = real(j - 1, dp)/real(ny, dp)
        end do
        do j = 1, ny
            do i = 1, nx
                point = [ &
                    0.5_dp*(x_edges(i) + x_edges(i + 1)), &
                    0.5_dp*(y_edges(j) + y_edges(j + 1))]
                call evaluate_pressure_at_point( &
                    mesh, pressure_basis, pressure_map, pressure_dofs, point, &
                    values(i, j), status)
                if (status /= 0) error stop "DG1 gallery evaluation failed"
                exact_values(i, j) = exact_pressure(point)
            end do
        end do
        call figure(figsize=[8.5_dp, 5.5_dp])
        call pcolormesh(x_edges, y_edges, values, cmap="viridis")
        call colorbar(label="solved DG1 pressure p_h")
        call xlabel("x")
        call ylabel("y")
        call title("Solved RT-DG mixed Poisson pressure field")
        call savefig(output_directory//"/mixed_poisson_solution_2d.png")
        call record_gallery_stage("physical_solution")

        open (newunit=unit, file=output_directory//"/mixed_poisson_solution.csv", &
            status="replace", action="write")
        write (unit, "(a)") "x,y,numerical,exact"
        do j = 1, ny
            do i = 1, nx
                write (unit, "(4(es24.16,:,','))") &
                    0.5_dp*(x_edges(i) + x_edges(i + 1)), &
                    0.5_dp*(y_edges(j) + y_edges(j + 1)), values(i, j), &
                    exact_values(i, j)
            end do
        end do
        close (unit)
    end subroutine render_solution

    subroutine evaluate_pressure_at_point( &
            mesh, basis, pressure_map, pressure_dofs, point, value, status)
        type(mesh_2d_t), intent(inout) :: mesh
        type(triangle_lagrange_basis_t), intent(in) :: basis
        integer, intent(in) :: pressure_map(:, :)
        real(dp), intent(in) :: pressure_dofs(:), point(2)
        real(dp), intent(out) :: value
        integer, intent(out) :: status

        real(dp) :: vertices(2, 3), reference_values(3), gradients(2, 3)
        real(dp) :: determinant, xi, eta, dx, dy
        integer :: triangle
        logical :: found

        found = .false.
        value = 0.0_dp
        status = 0
        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            determinant = (vertices(1, 2) - vertices(1, 1))* &
                (vertices(2, 3) - vertices(2, 1)) - &
                (vertices(1, 3) - vertices(1, 1))* &
                (vertices(2, 2) - vertices(2, 1))
            dx = point(1) - vertices(1, 1)
            dy = point(2) - vertices(2, 1)
            xi = (dx*(vertices(2, 3) - vertices(2, 1)) - &
                (vertices(1, 3) - vertices(1, 1))*dy)/determinant
            eta = ((vertices(1, 2) - vertices(1, 1))*dy - &
                dx*(vertices(2, 2) - vertices(2, 1)))/determinant
            if (xi >= -1.0e-12_dp .and. eta >= -1.0e-12_dp .and. &
                xi + eta <= 1.0_dp + 1.0e-12_dp) then
                call evaluate_triangle_lagrange_basis( &
                    basis, xi, eta, reference_values, gradients, status)
                if (status /= 0) return
                value = dot_product(reference_values, &
                    pressure_dofs(pressure_map(:, triangle)))
                found = .true.
                exit
            end if
        end do
        if (.not. found) status = 1
    end subroutine evaluate_pressure_at_point

    subroutine solve_and_measure( &
            vertex_count, errors, pressure_out, mesh_out)
        integer, intent(in) :: vertex_count
        real(dp), intent(out) :: errors(2)
        real(dp), allocatable, intent(out), optional :: pressure_out(:)
        type(mesh_2d_t), intent(out), optional :: mesh_out

        type(form_expr_t) :: balance_form, flux_form, pressure_flux_form
        type(fortsparse_status_t) :: status
        type(function_space_t), target :: pressure_space
        type(mesh_t), target :: mesh
        type(test_function_t) :: pressure_test
        type(trial_function_t) :: pressure_trial
        type(vector_function_space_t), target :: flux_space
        type(vector_test_function_t) :: flux_test
        type(vector_trial_function_t) :: flux_trial
        real(dp), allocatable :: flux_dofs(:), pressure_dofs(:)

        mesh = rectangle_mesh( &
            vertex_count, vertex_count, [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp])
        flux_space = vector_function_space(mesh, "RT", 1)
        pressure_space = function_space(mesh, "DG", 1)
        flux_trial = vector_trial_function(flux_space)
        flux_test = vector_test_function(flux_space)
        pressure_trial = trial_function(pressure_space)
        pressure_test = test_function(pressure_space)
        flux_form = inner(flux_trial, flux_test)*dx
        pressure_flux_form = &
            (-1.0_dp)*inner(pressure_trial, div(flux_test))*dx
        balance_form = inner(div(flux_trial), pressure_test)*dx
        call solve_symbolic_mixed_poisson_rt( &
            mesh%data, 1, 8, flux_form, pressure_flux_form, balance_form, &
            exact_source, flux_dofs, pressure_dofs, status)
        if (status%code /= 0) error stop "symbolic mixed solve failed"
        call compute_errors( &
            mesh%data, flux_dofs, pressure_dofs, errors)
        if (present(pressure_out)) then
            allocate(pressure_out(size(pressure_dofs)))
            pressure_out = pressure_dofs
        end if
        if (present(mesh_out)) mesh_out = mesh%data
    end subroutine solve_and_measure

    subroutine compute_errors(mesh, flux_dofs, pressure_dofs, errors)
        type(mesh_2d_t), intent(inout) :: mesh
        real(dp), intent(in) :: flux_dofs(:), pressure_dofs(:)
        real(dp), intent(out) :: errors(2)

        type(triangle_lagrange_basis_t) :: pressure_basis
        type(triangle_rt_basis_t) :: flux_basis
        integer, allocatable :: flux_map(:, :), pressure_map(:, :)
        integer, allocatable :: transforms(:, :)
        real(dp), allocatable :: eta(:), local_flux(:), local_pressure(:)
        real(dp), allocatable :: pressure_gradients(:, :)
        real(dp), allocatable :: pressure_values(:), weights(:), xi(:)
        real(dp) :: determinant, divergence, exact_flux(2)
        real(dp) :: numerical_flux(2), numerical_pressure, point(2)
        real(dp) :: vertices(2, 3)
        integer :: flux_count, local_dof, point_index, pressure_count
        integer :: status, triangle

        call build_triangle_trimmed_dof_map( &
            mesh, 2, flux_map, transforms, flux_count, status)
        if (status /= 0) error stop "RT1 map failed"
        call build_triangle_discontinuous_dof_map( &
            mesh, 1, pressure_map, pressure_count, status)
        if (status /= 0) error stop "DG1 map failed"
        call initialize_triangle_raviart_thomas(1, flux_basis, status)
        if (status /= 0) error stop "RT1 basis failed"
        call initialize_triangle_lagrange_basis(1, pressure_basis, status)
        if (status /= 0) error stop "DG1 basis failed"
        call triangle_duffy_quadrature(10, xi, eta, weights, status)
        if (status /= 0) error stop "mixed quadrature failed"
        allocate(local_flux(size(flux_map, 1)))
        allocate(local_pressure(size(pressure_map, 1)))
        allocate(pressure_values(size(pressure_map, 1)))
        allocate(pressure_gradients(2, size(pressure_map, 1)))
        errors = 0.0_dp
        do triangle = 1, mesh%n_triangles
            do local_dof = 1, size(flux_map, 1)
                local_flux(local_dof) = &
                    real(transforms(local_dof, triangle), dp)* &
                    flux_dofs(flux_map(local_dof, triangle))
            end do
            local_pressure = pressure_dofs(pressure_map(:, triangle))
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            determinant = &
                (vertices(1, 2) - vertices(1, 1))* &
                (vertices(2, 3) - vertices(2, 1)) - &
                (vertices(2, 2) - vertices(2, 1))* &
                (vertices(1, 3) - vertices(1, 1))
            do point_index = 1, size(weights)
                point = vertices(:, 1) + &
                    xi(point_index)*(vertices(:, 2) - vertices(:, 1)) + &
                    eta(point_index)*(vertices(:, 3) - vertices(:, 1))
                call evaluate_triangle_rt_interpolant( &
                    vertices, flux_basis, local_flux, xi(point_index), &
                    eta(point_index), numerical_flux, divergence, status)
                if (status /= 0) error stop "RT1 evaluation failed"
                call evaluate_triangle_lagrange_basis( &
                    pressure_basis, xi(point_index), eta(point_index), &
                    pressure_values, pressure_gradients, status)
                if (status /= 0) error stop "DG1 evaluation failed"
                numerical_pressure = dot_product( &
                    pressure_values, local_pressure)
                exact_flux = exact_darcy_flux(point)
                errors(1) = errors(1) + determinant*weights(point_index)* &
                    sum((numerical_flux - exact_flux)**2)
                errors(2) = errors(2) + determinant*weights(point_index)* &
                    (numerical_pressure - exact_pressure(point))**2
            end do
        end do
        errors = sqrt(errors)
    end subroutine compute_errors

    pure real(dp) function exact_pressure(point) result(value)
        real(dp), intent(in) :: point(2)
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = sin(pi*point(1))*sin(pi*point(2))
    end function exact_pressure

    pure function exact_darcy_flux(point) result(value)
        real(dp), intent(in) :: point(2)
        real(dp) :: value(2)
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = [ &
            -pi*cos(pi*point(1))*sin(pi*point(2)), &
            -pi*sin(pi*point(1))*cos(pi*point(2))]
    end function exact_darcy_flux

    pure real(dp) function exact_source(x, y) result(value)
        real(dp), intent(in) :: x, y
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = 2.0_dp*pi**2*sin(pi*x)*sin(pi*y)
    end function exact_source

    subroutine initialize_gallery_sequence()
        integer :: unit

        open (newunit=unit, file=output_directory//"/gallery_sequence.txt", &
            status="replace", action="write")
        close (unit)
    end subroutine initialize_gallery_sequence

    subroutine record_gallery_stage(stage)
        character(*), intent(in) :: stage
        integer :: unit

        open (newunit=unit, file=output_directory//"/gallery_sequence.txt", &
            status="old", position="append", action="write")
        write (unit, "(a)") stage
        close (unit)
    end subroutine record_gallery_stage

end program mixed_poisson
