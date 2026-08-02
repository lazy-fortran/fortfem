program anisotropic_tensor_diffusion
    !! Manufactured P1 diffusion with a strongly anisotropic tensor.
    use fortfem_feec, only: assemble_tensor_diffusion_matrix
    use fortfem_kinds, only: dp
    use fortplot, only: colorbar, contourf, figure, legend, plot, savefig, &
        title, xlabel, ylabel
    use fortnum_linalg, only: dense_solve
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: divisions = 12
    integer, parameter :: nodes_per_side = divisions + 1
    integer, parameter :: node_count = nodes_per_side*nodes_per_side
    integer, parameter :: triangle_count = 2*divisions*divisions
    real(dp), parameter :: anisotropy = 1000.0_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: tensor(2, 2) = reshape([ &
        anisotropy, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
    real(dp), parameter :: blue(3) = [0.0_dp, 114.0_dp, 178.0_dp]/255.0_dp
    real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
    real(dp), parameter :: gray(3) = [0.35_dp, 0.35_dp, 0.35_dp]
    character(*), parameter :: output_directory = &
        "output/example/anisotropic_tensor_diffusion"
    real(dp) :: coordinates(2, node_count), matrix(node_count, node_count)
    real(dp) :: right_hand_side(node_count), solution(node_count)
    real(dp) :: exact_solution(node_count), local_matrix(3, 3)
    real(dp) :: triangle_coordinates(2, 3), basis_gradients(1, 3, 2)
    real(dp) :: tensor_at_quadrature(1, 2, 2), quadrature_weight(1)
    real(dp) :: grid(nodes_per_side), solution_map(nodes_per_side, nodes_per_side)
    real(dp) :: exact_map(nodes_per_side, nodes_per_side)
    real(dp) :: error_map(nodes_per_side, nodes_per_side)
    real(dp) :: centerline(nodes_per_side), exact_centerline(nodes_per_side)
    real(dp) :: x, y, centroid_x, centroid_y, determinant, area
    real(dp) :: maximum_error, energy, start_time, elapsed
    integer :: ix, iy, node, triangle, local, info, unit
    integer :: nodes(3), lower_left, lower_right, upper_left, upper_right
    type(fortsparse_status_t) :: status

    call execute_command_line("mkdir -p "//output_directory)
    matrix = 0.0_dp
    right_hand_side = 0.0_dp
    do iy = 0, divisions
        y = real(iy, dp)/real(divisions, dp)
        grid(iy + 1) = y
        do ix = 0, divisions
            x = real(ix, dp)/real(divisions, dp)
            coordinates(:, node_index(ix, iy)) = [x, y]
            exact_solution(node_index(ix, iy)) = manufactured_solution(x, y)
        end do
    end do

    call cpu_time(start_time)
    triangle = 0
    do iy = 0, divisions - 1
        do ix = 0, divisions - 1
            lower_left = node_index(ix, iy)
            lower_right = node_index(ix + 1, iy)
            upper_left = node_index(ix, iy + 1)
            upper_right = node_index(ix + 1, iy + 1)
            triangle = triangle + 1
            call assemble_triangle( &
                [lower_left, lower_right, upper_right], triangle)
            triangle = triangle + 1
            call assemble_triangle( &
                [lower_left, upper_right, upper_left], triangle)
        end do
    end do

    do iy = 0, divisions
        do ix = 0, divisions
            node = node_index(ix, iy)
            if (ix == 0 .or. ix == divisions .or. &
                iy == 0 .or. iy == divisions) then
                matrix(node, :) = 0.0_dp
                matrix(:, node) = 0.0_dp
                matrix(node, node) = 1.0_dp
                right_hand_side(node) = 0.0_dp
            end if
        end do
    end do
    call dense_solve(matrix, right_hand_side, solution, info)
    call cpu_time(elapsed)
    elapsed = elapsed - start_time
    if (info /= 0) error stop "anisotropic tensor solve failed"

    maximum_error = 0.0_dp
    do iy = 0, divisions
        do ix = 0, divisions
            node = node_index(ix, iy)
            solution_map(ix + 1, iy + 1) = solution(node)
            exact_map(ix + 1, iy + 1) = exact_solution(node)
            error_map(ix + 1, iy + 1) = abs( &
                solution(node) - exact_solution(node))
            maximum_error = max(maximum_error, abs( &
                solution(node) - exact_solution(node)))
        end do
    end do
    do ix = 0, divisions
        centerline(ix + 1) = solution(node_index(ix, divisions/2))
        exact_centerline(ix + 1) = exact_solution(node_index(ix, divisions/2))
    end do
    energy = dot_product(solution, matmul(matrix, solution))
    if (maximum_error > 8.0e-2_dp) &
        error stop "anisotropic manufactured error oracle failed"
    if (energy <= 0.0_dp) error stop "anisotropic energy positivity failed"
    call render_plots()

    open (newunit=unit, file=output_directory//"/benchmark.txt", &
        status="replace", action="write")
    write (unit, "(a)") "quantity,value"
    write (unit, "(a,i0)") "divisions,", divisions
    write (unit, "(a,es24.16)") "anisotropy_ratio,", anisotropy
    write (unit, "(a,es24.16)") "maximum_nodal_error,", maximum_error
    write (unit, "(a,es24.16)") "discrete_energy,", energy
    write (unit, "(a,es24.16)") "assembly_and_solve_seconds,", elapsed
    close (unit)

contains

    integer function node_index(i, j)
        integer, intent(in) :: i, j

        node_index = i + 1 + j*nodes_per_side
    end function node_index

    real(dp) function manufactured_solution(x_value, y_value)
        real(dp), intent(in) :: x_value, y_value

        manufactured_solution = sin(pi*x_value)*sin(pi*y_value)
    end function manufactured_solution

    real(dp) function manufactured_source(x_value, y_value)
        real(dp), intent(in) :: x_value, y_value

        manufactured_source = pi*pi*(anisotropy + 1.0_dp)* &
            manufactured_solution(x_value, y_value)
    end function manufactured_source

    subroutine assemble_triangle(triangle_nodes, triangle_number)
        integer, intent(in) :: triangle_nodes(3), triangle_number

        integer :: first_basis, second_basis

        triangle_coordinates = coordinates(:, triangle_nodes)
        determinant = (triangle_coordinates(1, 2) - triangle_coordinates(1, 1))* &
            (triangle_coordinates(2, 3) - triangle_coordinates(2, 1)) - &
            (triangle_coordinates(1, 3) - triangle_coordinates(1, 1))* &
            (triangle_coordinates(2, 2) - triangle_coordinates(2, 1))
        area = 0.5_dp*abs(determinant)
        if (determinant <= 0.0_dp .or. area <= 0.0_dp) &
            error stop "anisotropic triangle orientation failed"
        basis_gradients = 0.0_dp
        basis_gradients(1, :, 1) = [ &
            triangle_coordinates(2, 2) - triangle_coordinates(2, 3), &
            triangle_coordinates(2, 3) - triangle_coordinates(2, 1), &
            triangle_coordinates(2, 1) - triangle_coordinates(2, 2)]/determinant
        basis_gradients(1, :, 2) = [ &
            triangle_coordinates(1, 3) - triangle_coordinates(1, 2), &
            triangle_coordinates(1, 1) - triangle_coordinates(1, 3), &
            triangle_coordinates(1, 2) - triangle_coordinates(1, 1)]/determinant
        tensor_at_quadrature(1, :, :) = tensor
        quadrature_weight(1) = area
        call assemble_tensor_diffusion_matrix( &
            basis_gradients, tensor_at_quadrature, quadrature_weight, &
            local_matrix, status)
        if (status%code /= 0) error stop "anisotropic local assembly failed"
        do second_basis = 1, 3
            do first_basis = 1, 3
                matrix(triangle_nodes(first_basis), triangle_nodes(second_basis)) = &
                    matrix(triangle_nodes(first_basis), triangle_nodes(second_basis)) + &
                    local_matrix(first_basis, second_basis)
            end do
        end do
        centroid_x = sum(triangle_coordinates(1, :))/3.0_dp
        centroid_y = sum(triangle_coordinates(2, :))/3.0_dp
        do local = 1, 3
            right_hand_side(triangle_nodes(local)) = &
                right_hand_side(triangle_nodes(local)) + area/3.0_dp* &
                manufactured_source(centroid_x, centroid_y)
        end do
        if (triangle_number < 1 .or. triangle_number > triangle_count) &
            error stop "anisotropic triangle count failed"
    end subroutine assemble_triangle

    subroutine render_plots()
        real(dp) :: line_x(nodes_per_side), line_y(nodes_per_side)
        integer :: line

        call figure(figsize=[8.0_dp, 6.0_dp])
        call contourf(grid, grid, solution_map, cmap="viridis", &
            show_colorbar=.true.)
        call colorbar(label="computed u")
        call xlabel("x")
        call ylabel("y")
        call title("Anisotropic tensor diffusion solution")
        call savefig(output_directory//"/anisotropic_solution_2d.png")

        call figure(figsize=[8.5_dp, 5.5_dp])
        call plot(grid, centerline, label="P1 tensor FEM", color=blue)
        call plot(grid, exact_centerline, label="manufactured solution", &
            color=orange, linestyle="--")
        call xlabel("x at y = 0.5")
        call ylabel("u")
        call title("Anisotropic diffusion centerline")
        call legend()
        call savefig(output_directory//"/anisotropic_centerline_1d.png")

        call figure(figsize=[8.0_dp, 6.0_dp])
        call contourf(grid, grid, error_map, cmap="magma", &
            show_colorbar=.true.)
        call colorbar(label="absolute nodal error")
        call xlabel("x")
        call ylabel("y")
        call title("Anisotropic diffusion absolute error")
        call savefig(output_directory//"/anisotropic_error_2d.png")

        call figure(figsize=[7.0_dp, 6.0_dp])
        do line = 1, nodes_per_side
            line_x = grid
            line_y = grid(line)
            call plot(line_x, line_y, color=gray, linewidth=0.5_dp)
            line_y = grid
            line_x = grid(line)
            call plot(line_x, line_y, color=gray, linewidth=0.5_dp)
        end do
        call xlabel("x")
        call ylabel("y")
        call title("Anisotropic diffusion structured P1 mesh")
        call savefig(output_directory//"/anisotropic_mesh_2d.png")
    end subroutine render_plots

end program anisotropic_tensor_diffusion
