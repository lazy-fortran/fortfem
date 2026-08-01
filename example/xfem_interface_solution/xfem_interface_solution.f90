program xfem_interface_solution
    !! Manufactured scalar and vector shifted-XFEM interface fields.
    use fortfem_api, only: &
        evaluate_shifted_enriched_space, &
        evaluate_shifted_vector_enriched_space
    use fortfem_kinds, only: dp
    use fortplot, only: colorbar, contourf, figure, pcolormesh, plot, quiver, savefig, &
        title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: nx = 49, ny = 49
    integer, parameter :: point_count = nx*ny, basis_count = 3
    integer, parameter :: component_count = 2
    real(dp), parameter :: interface_offset = 1.03123_dp
    real(dp), parameter :: scalar_coefficients(basis_count) = [ &
        0.25_dp, -0.40_dp, 0.30_dp]
    real(dp), parameter :: vector_coefficients(component_count, basis_count) = &
        reshape([0.30_dp, 0.10_dp, -0.20_dp, 0.40_dp, 0.10_dp, -0.25_dp], &
        [component_count, basis_count])
    real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
    character(*), parameter :: output_directory = &
        "output/example/xfem_interface_solution"
    real(dp) :: base_values(basis_count, point_count)
    real(dp) :: vector_base(component_count, basis_count, point_count)
    real(dp) :: level_values(point_count), anchor_values(basis_count)
    real(dp) :: enriched_values(basis_count, point_count)
    real(dp) :: vector_enriched(component_count, basis_count, point_count)
    real(dp) :: scalar_map(nx, ny), jump_map(nx, ny)
    real(dp) :: vector_magnitude(nx, ny)
    real(dp) :: x_grid(nx), y_grid(ny), x_edges(nx + 1), y_edges(ny + 1)
    real(dp) :: interface_x(2), interface_y(2)
    real(dp) :: quiver_x(9*9), quiver_y(9*9), quiver_u(9*9), quiver_v(9*9)
    real(dp) :: x, y, scalar_value, exact_scalar, vector_value(component_count)
    real(dp) :: exact_vector(component_count), maximum_error, vector_error
    real(dp) :: field_norm, start_time, elapsed
    integer :: index, ix, iy, point, basis, component, quiver_count, unit
    type(fortsparse_status_t) :: status

    call execute_command_line("mkdir -p "//output_directory)
    anchor_values = [-0.4_dp, -0.6_dp, -0.8_dp]
    do iy = 1, ny
        y = real(iy - 1, dp)/real(ny - 1, dp)
        y_grid(iy) = y
        do ix = 1, nx
            x = real(ix - 1, dp)/real(nx - 1, dp)
            x_grid(ix) = x
            point = point_index(ix, iy)
            level_values(point) = x + y - interface_offset
            base_values(:, point) = [1.0_dp, x, y]
            vector_base(:, :, point) = reshape([ &
                1.0_dp, y, x, 1.0_dp, x, y], &
                [component_count, basis_count])
        end do
    end do
    x_edges(2:nx) = 0.5_dp*(x_grid(:nx - 1) + x_grid(2:nx))
    y_edges(2:ny) = 0.5_dp*(y_grid(:ny - 1) + y_grid(2:ny))
    x_edges(1) = x_grid(1) - 0.5_dp*(x_grid(2) - x_grid(1))
    x_edges(nx + 1) = x_grid(nx) + 0.5_dp*(x_grid(nx) - x_grid(nx - 1))
    y_edges(1) = y_grid(1) - 0.5_dp*(y_grid(2) - y_grid(1))
    y_edges(ny + 1) = y_grid(ny) + 0.5_dp*(y_grid(ny) - y_grid(ny - 1))

    call cpu_time(start_time)
    call evaluate_shifted_enriched_space( &
        base_values, level_values, anchor_values, enriched_values, status)
    if (status%code /= 0) error stop "scalar shifted-XFEM construction failed"
    call evaluate_shifted_vector_enriched_space( &
        vector_base, level_values, anchor_values, vector_enriched, status)
    if (status%code /= 0) error stop "vector shifted-XFEM construction failed"
    call cpu_time(elapsed)
    elapsed = elapsed - start_time

    maximum_error = 0.0_dp
    vector_error = 0.0_dp
    do iy = 1, ny
        do ix = 1, nx
            point = point_index(ix, iy)
            scalar_value = dot_product(scalar_coefficients, base_values(:, point)) + &
                dot_product(scalar_coefficients, enriched_values(:, point))
            exact_scalar = dot_product(scalar_coefficients, base_values(:, point)) + &
                dot_product(scalar_coefficients, base_values(:, point))* &
                heaviside(level_values(point))
            scalar_map(ix, iy) = scalar_value
            jump_map(ix, iy) = dot_product( &
                scalar_coefficients, enriched_values(:, point))
            maximum_error = max(maximum_error, abs(scalar_value - exact_scalar))
            do component = 1, component_count
                vector_value(component) = dot_product( &
                    vector_coefficients(component, :), &
                    vector_base(component, :, point)) + &
                    dot_product(vector_coefficients(component, :), &
                    vector_enriched(component, :, point))
                exact_vector(component) = dot_product( &
                    vector_coefficients(component, :), &
                    vector_base(component, :, point))* &
                    (1.0_dp + heaviside(level_values(point)))
            end do
            vector_magnitude(ix, iy) = sqrt(dot_product(vector_value, vector_value))
            vector_error = max(vector_error, maxval(abs(vector_value - exact_vector)))
        end do
    end do
    if (maximum_error > 2.0e-14_dp) error stop "scalar shifted-XFEM oracle failed"
    if (vector_error > 2.0e-14_dp) error stop "vector shifted-XFEM oracle failed"

    quiver_count = 0
    do iy = 3, ny - 2, 5
        do ix = 3, nx - 2, 5
            quiver_count = quiver_count + 1
            point = point_index(ix, iy)
            quiver_x(quiver_count) = x_grid(ix)
            quiver_y(quiver_count) = y_grid(iy)
            quiver_u(quiver_count) = dot_product( &
                vector_coefficients(1, :), vector_base(1, :, point)) + &
                dot_product(vector_coefficients(1, :), vector_enriched(1, :, point))
            quiver_v(quiver_count) = dot_product( &
                vector_coefficients(2, :), vector_base(2, :, point)) + &
                dot_product(vector_coefficients(2, :), vector_enriched(2, :, point))
            field_norm = sqrt(quiver_u(quiver_count)**2 + &
                quiver_v(quiver_count)**2)
            if (field_norm > epsilon(1.0_dp)) then
                quiver_u(quiver_count) = &
                    0.075_dp*quiver_u(quiver_count)/field_norm
                quiver_v(quiver_count) = &
                    0.075_dp*quiver_v(quiver_count)/field_norm
            end if
        end do
    end do

    interface_x = [interface_offset - 1.0_dp, 1.0_dp]
    interface_y = [1.0_dp, interface_offset - 1.0_dp]
    call render_plots(quiver_count)

    open (newunit=unit, file=output_directory//"/benchmark.txt", &
        status="replace", action="write")
    write (unit, "(a)") "quantity,value"
    write (unit, "(a,i0)") "grid_points,", point_count
    write (unit, "(a,es24.16)") "maximum_scalar_oracle_error,", maximum_error
    write (unit, "(a,es24.16)") "maximum_vector_oracle_error,", vector_error
    write (unit, "(a,es24.16)") "space_construction_seconds,", elapsed
    close (unit)

contains

    integer function point_index(i, j)
        integer, intent(in) :: i, j

        point_index = i + (j - 1)*nx
    end function point_index

    pure real(dp) function heaviside(value)
        real(dp), intent(in) :: value

        if (value > 0.0_dp) then
            heaviside = 1.0_dp
        else
            heaviside = 0.0_dp
        end if
    end function heaviside

    subroutine render_plots(count)
        integer, intent(in) :: count

        call figure(figsize=[8.0_dp, 6.5_dp])
        call contourf(x_grid, y_grid, scalar_map, cmap="viridis", &
            show_colorbar=.true.)
        call colorbar(label="scalar enriched solution u")
        call plot(interface_x, interface_y, color=orange, linewidth=3.0_dp)
        call xlabel("x")
        call ylabel("y")
        call title("Shifted XFEM manufactured interface solution")
        call savefig(output_directory//"/xfem_interface_solution_2d.png")

        call figure(figsize=[8.0_dp, 6.5_dp])
        call pcolormesh(x_edges, y_edges, transpose(vector_magnitude), cmap="plasma")
        call plot(interface_x, interface_y, color=orange, linewidth=3.0_dp)
        call quiver(quiver_x(:count), quiver_y(:count), quiver_u(:count), &
            quiver_v(:count), scale=1.0_dp, scale_units="xy", &
            color="black", angles="xy", width=0.003_dp, headwidth=3.0_dp)
        call colorbar(label="vector field magnitude")
        call xlabel("x")
        call ylabel("y")
        call title("Vector shifted-XFEM field and interface")
        call savefig(output_directory//"/xfem_vector_field_2d.png")

        call figure(figsize=[8.0_dp, 6.5_dp])
        call contourf(x_grid, y_grid, jump_map, cmap="coolwarm", &
            show_colorbar=.true.)
        call colorbar(label="enriched jump contribution")
        call plot(interface_x, interface_y, color=orange, linewidth=3.0_dp)
        call xlabel("x")
        call ylabel("y")
        call title("Shifted-XFEM jump contribution")
        call savefig(output_directory//"/xfem_interface_jump_2d.png")
    end subroutine render_plots

end program xfem_interface_solution
