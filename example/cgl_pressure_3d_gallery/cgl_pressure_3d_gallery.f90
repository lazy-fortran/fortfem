program cgl_pressure_3d_gallery
    !! Physical-first 3-D CGL pressure/stress manufactured gallery.
    !!
    !! This is deliberately a neutral constitutive fixture.  It exercises the
    !! generated CGL tensor, tensor-normal traction, and pressure-work
    !! contractions on an oblique, spatially varying unit direction.  No
    !! equilibrium closure, plasma reader, or application-specific geometry is
    !! assumed.
    use fortfem_api, only: &
        evaluate_cgl_pressure_tensor, evaluate_cgl_pressure_tensor_jvp, &
        evaluate_cgl_pressure_tensor_vjp, evaluate_cgl_pressure_traction, &
        evaluate_cgl_pressure_traction_jvp, evaluate_cgl_pressure_traction_vjp, &
        evaluate_cgl_pressure_work, evaluate_cgl_pressure_work_jvp, &
        evaluate_cgl_pressure_work_vjp
    use fortfem_kinds, only: dp
    use fortplot, only: add_3d_plot, add_scatter, colorbar, figure, legend, &
        plot, savefig, title, view_init, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: points_per_axis = 4
    integer, parameter :: node_count = points_per_axis**3
    integer, parameter :: timing_repetitions = 250
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(*), parameter :: output_directory = &
        "output/example/cgl_pressure_3d_gallery"
    real(dp), parameter :: blue(3) = [0.05_dp, 0.20_dp, 0.70_dp]
    real(dp), parameter :: orange(3) = [0.90_dp, 0.35_dp, 0.05_dp]
    real(dp), parameter :: edge_color(3) = [0.35_dp, 0.35_dp, 0.35_dp]

    real(dp) :: coordinates(3, node_count)
    real(dp) :: p_parallel_values(node_count), p_perpendicular_values(node_count)
    real(dp) :: direction_values(3, node_count), normal_values(3, node_count)
    real(dp) :: pressure_values(3, 3, node_count)
    real(dp) :: traction_values(3, node_count)
    real(dp) :: gradient_values(3, 3, node_count)
    real(dp) :: work_values(node_count), trace_values(node_count)
    real(dp) :: tensor(3, 3), tensor_dot(3, 3), tensor_expected(3, 3)
    real(dp) :: traction(3), traction_dot(3), traction_expected(3)
    real(dp) :: work, work_dot, work_expected
    real(dp) :: point(3), direction_dot(3), normal_dot(3), gradient_dot(3, 3)
    real(dp) :: p_parallel_dot, p_perpendicular_dot
    real(dp) :: tensor_bar(3, 3), traction_bar(3), work_bar
    real(dp) :: p_parallel_bar, p_perpendicular_bar
    real(dp) :: direction_bar(3), normal_bar(3), gradient_bar(3, 3)
    real(dp) :: direction_dot_tangent(3), identity(3, 3), delta
    real(dp) :: tensor_error, traction_error, work_error
    real(dp) :: tensor_jvp_error, traction_jvp_error, work_jvp_error
    real(dp) :: tensor_vjp_error, traction_vjp_error, work_vjp_error
    real(dp) :: maximum_tensor_error, maximum_traction_error
    real(dp) :: maximum_work_error, maximum_tensor_jvp_error
    real(dp) :: maximum_traction_jvp_error, maximum_work_jvp_error
    real(dp) :: maximum_tensor_vjp_error, maximum_traction_vjp_error
    real(dp) :: maximum_work_vjp_error, minimum_trace, kernel_seconds
    real(dp) :: start_time, arrow_scale, traction_scale, magnitude
    real(dp) :: node_x(node_count), node_y(node_count), node_z(node_count)
    real(dp) :: sample_index(node_count)
    integer :: node, ix, iy, iz, repetition, unit, command_status
    type(fortsparse_status_t) :: status

    command_status = 0
    call execute_command_line("mkdir -p "//output_directory, &
        exitstat=command_status)
    if (command_status /= 0) error stop "cannot create CGL 3-D output directory"
    call initialize_gallery_sequence()
    identity = 0.0_dp
    identity(1, 1) = 1.0_dp
    identity(2, 2) = 1.0_dp
    identity(3, 3) = 1.0_dp
    tensor_bar = reshape([ &
        0.2_dp, -0.4_dp, 0.5_dp, 0.1_dp, 0.7_dp, -0.3_dp, &
        -0.6_dp, 0.8_dp, 0.9_dp], [3, 3])
    traction_bar = [0.6_dp, -0.2_dp, 0.8_dp]
    work_bar = 0.73_dp

    maximum_tensor_error = 0.0_dp
    maximum_traction_error = 0.0_dp
    maximum_work_error = 0.0_dp
    maximum_tensor_jvp_error = 0.0_dp
    maximum_traction_jvp_error = 0.0_dp
    maximum_work_jvp_error = 0.0_dp
    maximum_tensor_vjp_error = 0.0_dp
    maximum_traction_vjp_error = 0.0_dp
    maximum_work_vjp_error = 0.0_dp
    minimum_trace = huge(1.0_dp)
    call cpu_time(start_time)
    node = 0
    do iz = 0, points_per_axis - 1
        do iy = 0, points_per_axis - 1
            do ix = 0, points_per_axis - 1
                node = node + 1
                point = [ &
                    -1.0_dp + 2.0_dp*real(ix, dp)/ &
                    real(points_per_axis - 1, dp), &
                    -1.0_dp + 2.0_dp*real(iy, dp)/ &
                    real(points_per_axis - 1, dp), &
                    -1.0_dp + 2.0_dp*real(iz, dp)/ &
                    real(points_per_axis - 1, dp)]
                coordinates(:, node) = point
                sample_index(node) = real(node, dp)
                call sample_inputs(point, p_parallel_values(node), &
                    p_perpendicular_values(node), direction_values(:, node), &
                    normal_values(:, node), gradient_values(:, :, node), &
                    direction_dot, normal_dot, gradient_dot, &
                    p_parallel_dot, p_perpendicular_dot)
                call evaluate_cgl_pressure_tensor( &
                    p_parallel_values(node), p_perpendicular_values(node), &
                    direction_values(:, node), pressure_values(:, :, node), status)
                if (status%code /= 0) error stop "CGL tensor evaluation failed"
                delta = p_parallel_values(node) - p_perpendicular_values(node)
                tensor_expected = p_perpendicular_values(node)*identity + &
                    delta*outer_product(direction_values(:, node), &
                    direction_values(:, node))
                tensor_error = maxval(abs(pressure_values(:, :, node) - &
                    tensor_expected))
                maximum_tensor_error = max(maximum_tensor_error, tensor_error)

                call evaluate_cgl_pressure_traction( &
                    p_parallel_values(node), p_perpendicular_values(node), &
                    direction_values(:, node), normal_values(:, node), traction, &
                    status)
                if (status%code /= 0) error stop "CGL traction evaluation failed"
                traction_expected = matmul(tensor_expected, normal_values(:, node))
                traction_error = maxval(abs(traction - traction_expected))
                maximum_traction_error = max(maximum_traction_error, &
                    traction_error)
                traction_values(:, node) = traction

                call evaluate_cgl_pressure_work( &
                    p_parallel_values(node), p_perpendicular_values(node), &
                    direction_values(:, node), gradient_values(:, :, node), work, &
                    status)
                if (status%code /= 0) error stop "CGL work evaluation failed"
                work_expected = sum(tensor_expected*gradient_values(:, :, node))
                work_error = abs(work - work_expected)
                maximum_work_error = max(maximum_work_error, work_error)
                work_values(node) = work
                trace_values(node) = sum(pressure_values(:, :, node))
                minimum_trace = min(minimum_trace, trace_values(node))

                direction_dot_tangent = direction_dot - &
                    direction_values(:, node)*dot_product( &
                    direction_values(:, node), direction_dot)
                call evaluate_cgl_pressure_tensor_jvp( &
                    p_parallel_values(node), p_perpendicular_values(node), &
                    direction_values(:, node), p_parallel_dot, &
                    p_perpendicular_dot, direction_dot_tangent, tensor_dot, status)
                tensor_expected = p_perpendicular_dot*identity + &
                    (p_parallel_dot - p_perpendicular_dot)* &
                    outer_product(direction_values(:, node), direction_values(:, node)) + &
                    delta*(outer_product(direction_dot_tangent, &
                    direction_values(:, node)) + outer_product( &
                    direction_values(:, node), direction_dot_tangent))
                tensor_jvp_error = maxval(abs(tensor_dot - tensor_expected))
                maximum_tensor_jvp_error = max(maximum_tensor_jvp_error, &
                    tensor_jvp_error)

                call evaluate_cgl_pressure_traction_jvp( &
                    p_parallel_values(node), p_perpendicular_values(node), &
                    direction_values(:, node), normal_values(:, node), &
                    p_parallel_dot, p_perpendicular_dot, direction_dot_tangent, &
                    normal_dot, traction_dot, status)
                traction_expected = matmul(tensor_dot, normal_values(:, node)) + &
                    matmul(pressure_values(:, :, node), normal_dot)
                traction_jvp_error = maxval(abs(traction_dot - traction_expected))
                maximum_traction_jvp_error = max(maximum_traction_jvp_error, &
                    traction_jvp_error)

                call evaluate_cgl_pressure_work_jvp( &
                    p_parallel_values(node), p_perpendicular_values(node), &
                    direction_values(:, node), gradient_values(:, :, node), &
                    p_parallel_dot, p_perpendicular_dot, direction_dot_tangent, &
                    gradient_dot, work_dot, status)
                work_expected = sum(tensor_dot*gradient_values(:, :, node)) + &
                    sum(pressure_values(:, :, node)*gradient_dot)
                work_jvp_error = abs(work_dot - work_expected)
                maximum_work_jvp_error = max(maximum_work_jvp_error, work_jvp_error)

                call evaluate_cgl_pressure_tensor_vjp( &
                    p_parallel_values(node), p_perpendicular_values(node), &
                    direction_values(:, node), tensor_bar, p_parallel_bar, &
                    p_perpendicular_bar, direction_bar, status)
                tensor_vjp_error = abs(sum(tensor_bar*tensor_dot) - &
                    (p_parallel_bar*p_parallel_dot + p_perpendicular_bar* &
                    p_perpendicular_dot + dot_product(direction_bar, &
                    direction_dot_tangent)))
                maximum_tensor_vjp_error = max(maximum_tensor_vjp_error, &
                    tensor_vjp_error)

                call evaluate_cgl_pressure_traction_vjp( &
                    p_parallel_values(node), p_perpendicular_values(node), &
                    direction_values(:, node), normal_values(:, node), &
                    traction_bar, p_parallel_bar, p_perpendicular_bar, &
                    direction_bar, normal_bar, status)
                traction_vjp_error = abs(dot_product(traction_bar, traction_dot) - &
                    (p_parallel_bar*p_parallel_dot + p_perpendicular_bar* &
                    p_perpendicular_dot + dot_product(direction_bar, &
                    direction_dot_tangent) + dot_product(normal_bar, normal_dot)))
                maximum_traction_vjp_error = max(maximum_traction_vjp_error, &
                    traction_vjp_error)

                call evaluate_cgl_pressure_work_vjp( &
                    p_parallel_values(node), p_perpendicular_values(node), &
                    direction_values(:, node), gradient_values(:, :, node), &
                    work_bar, p_parallel_bar, p_perpendicular_bar, direction_bar, &
                    gradient_bar, status)
                work_vjp_error = abs(work_bar*work_dot - (p_parallel_bar* &
                    p_parallel_dot + p_perpendicular_bar*p_perpendicular_dot + &
                    dot_product(direction_bar, direction_dot_tangent) + &
                    sum(gradient_bar*gradient_dot)))
                maximum_work_vjp_error = max(maximum_work_vjp_error, &
                    work_vjp_error)
            end do
        end do
    end do
    call cpu_time(kernel_seconds)
    kernel_seconds = kernel_seconds - start_time
    call check_oracle_bounds()
    call measure_kernel_time()
    call write_solution_csv()
    call write_diagnostics()
    call render_solution()
    call record_gallery_stage("physical_solution")
    call render_diagnostics()
    call record_gallery_stage("diagnostics")

contains

    pure function outer_product(left, right) result(product)
        real(dp), intent(in) :: left(3), right(3)
        real(dp) :: product(3, 3)

        product = spread(left, dim=2, ncopies=3)* &
            spread(right, dim=1, ncopies=3)
    end function outer_product

    pure real(dp) function vector_norm(vector)
        real(dp), intent(in) :: vector(3)

        vector_norm = sqrt(dot_product(vector, vector))
    end function vector_norm

    pure subroutine sample_inputs(point_value, p_parallel, p_perpendicular, &
            direction, normal, gradient, direction_dot, normal_dot, &
            gradient_dot, p_parallel_dot, p_perpendicular_dot)
        real(dp), intent(in) :: point_value(3)
        real(dp), intent(out) :: p_parallel, p_perpendicular, direction(3)
        real(dp), intent(out) :: normal(3), gradient(3, 3)
        real(dp), intent(out) :: direction_dot(3), normal_dot(3)
        real(dp), intent(out) :: gradient_dot(3, 3)
        real(dp), intent(out) :: p_parallel_dot, p_perpendicular_dot
        real(dp) :: raw_direction(3), raw_normal(3)

        p_parallel = 3.6_dp + 0.4_dp*point_value(1) - &
            0.25_dp*point_value(2) + 0.18_dp*point_value(3)
        p_perpendicular = 1.2_dp + 0.15_dp*cos(pi*point_value(1))* &
            cos(0.5_dp*pi*point_value(2))*cos(0.75_dp*pi*point_value(3))
        raw_direction = [ &
            1.0_dp + 0.25_dp*sin(pi*point_value(1)), &
            2.0_dp + 0.20_dp*cos(pi*point_value(2)), &
            2.0_dp + 0.15_dp*sin(pi*point_value(3))]
        direction = raw_direction/vector_norm(raw_direction)
        raw_normal = [ &
            1.0_dp + 0.35_dp*point_value(1), &
            -0.8_dp + 0.25_dp*point_value(2), &
            0.6_dp + 0.20_dp*point_value(3)]
        normal = raw_normal/vector_norm(raw_normal)

        gradient = 0.0_dp
        gradient(1, 1) = pi*cos(pi*point_value(1))
        gradient(1, 2) = 0.20_dp
        gradient(2, 2) = -0.40_dp*pi*sin(pi*point_value(2))
        gradient(2, 3) = 0.10_dp
        gradient(3, 1) = 0.15_dp
        gradient(3, 3) = 0.30_dp*pi*cos(pi*point_value(3))

        direction_dot = [0.04_dp, -0.03_dp, 0.02_dp]
        normal_dot = [-0.02_dp, 0.05_dp, 0.03_dp]
        gradient_dot = reshape([ &
            -0.30_dp + 0.01_dp*point_value(1), &
            0.15_dp, 0.20_dp, -0.10_dp, &
            0.25_dp + 0.01_dp*point_value(2), 0.40_dp, &
            0.05_dp, -0.20_dp, 0.10_dp + 0.01_dp*point_value(3)], [3, 3])
        p_parallel_dot = 0.11_dp + 0.01_dp*point_value(1)
        p_perpendicular_dot = -0.07_dp + 0.02_dp*point_value(2)
    end subroutine sample_inputs

    subroutine check_oracle_bounds()
        if (maximum_tensor_error > 2.0e-13_dp) &
            error stop "CGL 3-D tensor contraction oracle failed"
        if (maximum_traction_error > 2.0e-13_dp) &
            error stop "CGL 3-D traction oracle failed"
        if (maximum_work_error > 2.0e-13_dp) &
            error stop "CGL 3-D pressure-work oracle failed"
        if (maximum_tensor_jvp_error > 3.0e-13_dp) &
            error stop "CGL 3-D tensor JVP oracle failed"
        if (maximum_traction_jvp_error > 3.0e-13_dp) &
            error stop "CGL 3-D traction JVP oracle failed"
        if (maximum_work_jvp_error > 3.0e-13_dp) &
            error stop "CGL 3-D pressure-work JVP oracle failed"
        if (maximum_tensor_vjp_error > 3.0e-12_dp) &
            error stop "CGL 3-D tensor VJP oracle failed"
        if (maximum_traction_vjp_error > 3.0e-12_dp) &
            error stop "CGL 3-D traction VJP oracle failed"
        if (maximum_work_vjp_error > 3.0e-12_dp) &
            error stop "CGL 3-D pressure-work VJP oracle failed"
        if (minimum_trace <= 0.0_dp) error stop "CGL pressure trace is nonpositive"
    end subroutine check_oracle_bounds

    subroutine measure_kernel_time()
        integer :: repetition, local_node
        real(dp) :: local_tensor(3, 3), local_traction(3), local_work
        real(dp) :: local_start, local_seconds

        call cpu_time(local_start)
        do repetition = 1, timing_repetitions
            do local_node = 1, node_count
                call evaluate_cgl_pressure_tensor( &
                    p_parallel_values(local_node), &
                    p_perpendicular_values(local_node), &
                    direction_values(:, local_node), local_tensor, status)
                call evaluate_cgl_pressure_traction( &
                    p_parallel_values(local_node), &
                    p_perpendicular_values(local_node), &
                    direction_values(:, local_node), normal_values(:, local_node), &
                    local_traction, status)
                call evaluate_cgl_pressure_work( &
                    p_parallel_values(local_node), &
                    p_perpendicular_values(local_node), &
                    direction_values(:, local_node), gradient_values(:, :, local_node), &
                    local_work, status)
            end do
        end do
        call cpu_time(local_seconds)
        kernel_seconds = kernel_seconds + &
            (local_seconds - local_start)/real(timing_repetitions, dp)
    end subroutine measure_kernel_time

    subroutine write_solution_csv()
        open (newunit=unit, file=output_directory//"/solution.csv", &
            status="replace", action="write")
        write (unit, "(a)") &
            "node,x,y,z,p_parallel,p_perpendicular,bx,by,bz,P11,P12,P13,"// &
            "P22,P23,P33,nx,ny,nz,t1,t2,t3,g11,g12,g13,g21,g22,g23,"// &
            "g31,g32,g33,pressure_work"
        do node = 1, node_count
            write (unit, "(i0,30(',',es24.16))") node, coordinates(:, node), &
                p_parallel_values(node), p_perpendicular_values(node), &
                direction_values(:, node), pressure_values(1, 1, node), &
                pressure_values(1, 2, node), pressure_values(1, 3, node), &
                pressure_values(2, 2, node), pressure_values(2, 3, node), &
                pressure_values(3, 3, node), normal_values(:, node), &
                traction_values(:, node), gradient_values(:, :, node), &
                work_values(node)
        end do
        close (unit)
    end subroutine write_solution_csv

    subroutine write_diagnostics()
        open (newunit=unit, file=output_directory//"/diagnostics.csv", &
            status="replace", action="write")
        write (unit, "(a)") "quantity,value"
        write (unit, "(a,es24.16)") "maximum_tensor_error,", &
            maximum_tensor_error
        write (unit, "(a,es24.16)") "maximum_traction_error,", &
            maximum_traction_error
        write (unit, "(a,es24.16)") "maximum_work_error,", maximum_work_error
        write (unit, "(a,es24.16)") "maximum_tensor_jvp_error,", &
            maximum_tensor_jvp_error
        write (unit, "(a,es24.16)") "maximum_traction_jvp_error,", &
            maximum_traction_jvp_error
        write (unit, "(a,es24.16)") "maximum_work_jvp_error,", &
            maximum_work_jvp_error
        write (unit, "(a,es24.16)") "maximum_tensor_vjp_error,", &
            maximum_tensor_vjp_error
        write (unit, "(a,es24.16)") "maximum_traction_vjp_error,", &
            maximum_traction_vjp_error
        write (unit, "(a,es24.16)") "maximum_work_vjp_error,", &
            maximum_work_vjp_error
        write (unit, "(a,es24.16)") "minimum_pressure_trace,", minimum_trace
        write (unit, "(a,es24.16)") "kernel_seconds,", kernel_seconds
        write (unit, "(a,i0)") "points,", node_count
        close (unit)

        open (newunit=unit, file=output_directory//"/benchmark.txt", &
            status="replace", action="write")
        write (unit, "(a,i0)") "points: ", node_count
        write (unit, "(a,i0)") "timing repetitions: ", timing_repetitions
        write (unit, "(a,es24.16)") "kernel seconds: ", kernel_seconds
        write (unit, "(a,es24.16)") "maximum tensor error: ", &
            maximum_tensor_error
        write (unit, "(a,es24.16)") "maximum traction error: ", &
            maximum_traction_error
        close (unit)
    end subroutine write_diagnostics

    subroutine render_solution()
        integer :: line_node
        logical :: boundary_node
        real(dp) :: b_start(3), b_end(3), t_start(3), t_end(3)

        do line_node = 1, node_count
            node_x(line_node) = coordinates(1, line_node)
            node_y(line_node) = coordinates(2, line_node)
            node_z(line_node) = coordinates(3, line_node)
        end do
        arrow_scale = 0.18_dp
        magnitude = maxval(sqrt(sum(traction_values**2, dim=1)))
        traction_scale = 0.25_dp/max(1.0_dp, magnitude)
        call figure(figsize=[8.8_dp, 7.2_dp])
        call add_scatter(node_x, node_y, node_z, c=trace_values, cmap="viridis", &
            marker="o", markersize=46.0_dp, label="pressure trace")
        call colorbar(label="tr(P)")
        do line_node = 1, node_count
            b_start = coordinates(:, line_node)
            b_end = b_start + arrow_scale*direction_values(:, line_node)
            if (line_node == 1) then
                call add_3d_plot([b_start(1), b_end(1)], &
                    [b_start(2), b_end(2)], [b_start(3), b_end(3)], &
                    color=blue, linewidth=1.1_dp, label="oblique B direction")
            else
                call add_3d_plot([b_start(1), b_end(1)], &
                    [b_start(2), b_end(2)], [b_start(3), b_end(3)], &
                    color=blue, linewidth=1.1_dp)
            end if
            boundary_node = maxval(abs(coordinates(:, line_node))) > &
                1.0_dp - 1.0e-12_dp
            if (boundary_node) then
                t_start = b_start
                t_end = t_start + traction_scale*traction_values(:, line_node)
                if (line_node == 1) then
                    call add_3d_plot([t_start(1), t_end(1)], &
                        [t_start(2), t_end(2)], [t_start(3), t_end(3)], &
                        color=orange, linewidth=1.4_dp, label="surface traction")
                else
                    call add_3d_plot([t_start(1), t_end(1)], &
                        [t_start(2), t_end(2)], [t_start(3), t_end(3)], &
                        color=orange, linewidth=1.4_dp)
                end if
            end if
        end do
        call draw_box_edges()
        call xlabel("x")
        call ylabel("y")
        call title("3-D CGL pressure tensor: trace, B direction, and traction")
        call legend()
        call view_init(elev=24.0_dp, azim=-52.0_dp)
        call savefig(output_directory//"/solution.png")
        call savefig(output_directory//"/solution_3d.png")
    end subroutine render_solution

    subroutine draw_box_edges()
        integer :: edge
        real(dp), parameter :: endpoints(3, 2, 12) = reshape([ &
            -1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp, &
            -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, &
            -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp, &
            1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, &
            1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, &
            -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, &
            -1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, &
            -1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, &
            -1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, &
            1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
            -1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
            1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
            [3, 2, 12])

        do edge = 1, 12
            call add_3d_plot(endpoints(1, :, edge), endpoints(2, :, edge), &
                endpoints(3, :, edge), color=edge_color, linewidth=0.8_dp)
        end do
    end subroutine draw_box_edges

    subroutine render_diagnostics()
        call figure(figsize=[8.5_dp, 5.5_dp])
        call plot(sample_index, p_parallel_values, label="p_parallel")
        call plot(sample_index, p_perpendicular_values, label="p_perpendicular")
        call plot(sample_index, trace_values/3.0_dp, label="tr(P)/3")
        call xlabel("sample index")
        call ylabel("pressure")
        call title("CGL pressure components (3-D samples)")
        call legend()
        call savefig(output_directory//"/pressure_components.png")
    end subroutine render_diagnostics

    subroutine initialize_gallery_sequence()
        open (newunit=unit, file=output_directory//"/gallery_sequence.txt", &
            status="replace", action="write")
        close (unit)
    end subroutine initialize_gallery_sequence

    subroutine record_gallery_stage(stage)
        character(*), intent(in) :: stage

        open (newunit=unit, file=output_directory//"/gallery_sequence.txt", &
            status="old", position="append", action="write")
        write (unit, "(a)") trim(stage)
        close (unit)
    end subroutine record_gallery_stage

end program cgl_pressure_3d_gallery
