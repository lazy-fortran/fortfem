---
title: field_aligned_anisotropic_3d_gallery Example
---

# field_aligned_anisotropic_3d_gallery Example

# Field-aligned anisotropic 3-D diffusion

This neutral manufactured-solution gallery exercises the public three-dimensional
tensor-diffusion contraction with a strongly anisotropic field-aligned tensor.
The unit direction is a fixed oblique vector, so the large coefficient acts
along a direction that is not aligned with any coordinate axis:

\[
K = k_\perp I + (k_\parallel-k_\perp)\,b b^T,
\qquad b=(1,2,2)/3,
\qquad -\nabla\cdot(K\nabla u)=f.
\]

The exact solution is `u = sin(pi*x) sin(pi*y) sin(pi*z)` on the unit cube
with homogeneous Dirichlet data.  The tetrahedral P1 assembly evaluates `K`
through `evaluate_field_aligned_constitutive_tensor` and calls
`assemble_tensor_diffusion_matrix_3d`; no plasma model, closure, geometry
reader, or external data is involved.

The first figure is the computed 3-D solution with flux-direction arrows and
tetrahedral edges.  A subsequent mid-plane slice shows the scalar solution and
the in-plane anisotropic flux.  The executable writes an independent analytic
source/error/energy oracle to `diagnostics.csv`, plus a CSV nodal solution for
downstream plotting.

The constitutive construction follows the standard field-aligned projector
used in neutral anisotropic transport formulations.  The gallery is an
implementation fixture for future field-aligned FEM/IGA and not a plasma
equilibrium calculation.

## Usage

```bash
fpm run --example field_aligned_anisotropic_3d_gallery
```

## Source Code

```fortran
program field_aligned_anisotropic_3d_gallery
    !! Physical-first 3-D field-aligned anisotropic diffusion gallery.
    use fortfem_feec, only: &
        assemble_tensor_diffusion_matrix_3d, &
        evaluate_field_aligned_constitutive_tensor
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve, det3, inv3
    use fortplot, only: add_3d_plot, add_scatter, colorbar, contourf, figure, &
        legend, quiver, savefig, title, view_init, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: cells_per_axis = 4
    integer, parameter :: nodes_per_axis = cells_per_axis + 1
    integer, parameter :: node_count = nodes_per_axis**3
    integer, parameter :: tetrahedra_per_cube = 6
    integer, parameter :: cube_count = cells_per_axis**3
    integer, parameter :: tetrahedron_count = tetrahedra_per_cube*cube_count
    integer, parameter :: slice_count = nodes_per_axis**2
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: k_parallel = 5000.0_dp
    real(dp), parameter :: k_perpendicular = 1.0_dp
    real(dp), parameter :: field_direction(3) = [1.0_dp, 2.0_dp, 2.0_dp]/3.0_dp
    character(*), parameter :: output_directory = &
        "output/example/field_aligned_anisotropic_3d_gallery"
    real(dp), parameter :: arrow_color(3) = [0.08_dp, 0.08_dp, 0.08_dp]
    real(dp), parameter :: edge_color(3) = [0.42_dp, 0.42_dp, 0.42_dp]

    real(dp) :: coordinates(3, node_count)
    real(dp) :: system_matrix(node_count, node_count)
    real(dp) :: right_hand_side(node_count), solution(node_count)
    real(dp) :: exact_solution(node_count), nodal_error(node_count)
    real(dp) :: constitutive_tensor(3, 3)
    real(dp) :: basis_gradients(1, 4, 3), local_matrix(4, 4)
    real(dp) :: tensor_at_quadrature(1, 3, 3), quadrature_weight(1)
    real(dp) :: jacobian(3, 3), inverse_jacobian(3, 3)
    real(dp) :: tetra_vertices(3, 4), reference_gradients(4, 3)
    real(dp) :: centroid(3), determinant, volume
    real(dp) :: solution_slice(nodes_per_axis, nodes_per_axis)
    real(dp) :: slice_x(nodes_per_axis), slice_y(nodes_per_axis)
    real(dp) :: flux_x(slice_count), flux_y(slice_count)
    real(dp) :: flux_u(slice_count), flux_v(slice_count)
    real(dp) :: node_x(node_count), node_y(node_count), node_z(node_count)
    real(dp) :: flux_at_node(3), point(3), scale
    real(dp) :: maximum_error, energy, residual, elapsed, start_time
    real(dp) :: source_value
    integer :: tetrahedra(4, tetrahedron_count)
    integer :: nodes(8), local_nodes(4), tetrahedron
    integer :: j, cube, local, node, info, unit, command_status
    integer :: first, second, third, fourth, fifth, sixth, seventh, eighth
    type(fortsparse_status_t) :: status

    command_status = 0
    call execute_command_line("mkdir -p "//output_directory, &
        exitstat=command_status)
    if (command_status /= 0) error stop "cannot create anisotropic gallery output"
    call initialize_gallery_sequence()
    call initialize_constitutive_tensor()
    call initialize_mesh()
    system_matrix = 0.0_dp
    right_hand_side = 0.0_dp

    call cpu_time(start_time)
    do tetrahedron = 1, tetrahedron_count
        local_nodes = tetrahedra(:, tetrahedron)
        tetra_vertices = coordinates(:, local_nodes)
        jacobian(:, 1) = tetra_vertices(:, 2) - tetra_vertices(:, 1)
        jacobian(:, 2) = tetra_vertices(:, 3) - tetra_vertices(:, 1)
        jacobian(:, 3) = tetra_vertices(:, 4) - tetra_vertices(:, 1)
        determinant = det3(jacobian)
        if (determinant <= 0.0_dp) error stop "tetrahedron orientation failed"
        volume = determinant/6.0_dp
        call inv3(jacobian, inverse_jacobian, info)
        if (info /= 0) error stop "tetrahedron inverse failed"
        reference_gradients(1, :) = [-1.0_dp, -1.0_dp, -1.0_dp]
        reference_gradients(2, :) = [1.0_dp, 0.0_dp, 0.0_dp]
        reference_gradients(3, :) = [0.0_dp, 1.0_dp, 0.0_dp]
        reference_gradients(4, :) = [0.0_dp, 0.0_dp, 1.0_dp]
        do local = 1, 4
            basis_gradients(1, local, :) = matmul( &
                transpose(inverse_jacobian), reference_gradients(local, :))
        end do
        tensor_at_quadrature(1, :, :) = constitutive_tensor
        quadrature_weight(1) = volume
        call assemble_tensor_diffusion_matrix_3d( &
            basis_gradients, tensor_at_quadrature, quadrature_weight, &
            local_matrix, status)
        if (status%code /= 0) error stop "3-D anisotropic local assembly failed"
        centroid = sum(tetra_vertices, dim=2)/4.0_dp
        source_value = manufactured_source(centroid)
        do local = 1, 4
            node = local_nodes(local)
            right_hand_side(node) = right_hand_side(node) + &
                volume*source_value/4.0_dp
            do j = 1, 4
                system_matrix(node, local_nodes(j)) = &
                    system_matrix(node, local_nodes(j)) + local_matrix(local, j)
            end do
        end do
    end do
    call impose_dirichlet_data()
    call dense_solve(system_matrix, right_hand_side, solution, info)
    elapsed = 0.0_dp
    call cpu_time(elapsed)
    elapsed = elapsed - start_time
    if (info /= 0) error stop "field-aligned anisotropic solve failed"

    maximum_error = 0.0_dp
    exact_solution = 0.0_dp
    do node = 1, node_count
        exact_solution(node) = manufactured_solution(coordinates(:, node))
        nodal_error(node) = abs(solution(node) - exact_solution(node))
        maximum_error = max(maximum_error, nodal_error(node))
    end do
    residual = maxval(abs(matmul(system_matrix, solution) - right_hand_side))
    energy = dot_product(solution, matmul(system_matrix, solution))
    if (maximum_error > 0.15_dp) error stop "anisotropic solution oracle failed"
    if (residual > 5.0e-10_dp) error stop "anisotropic residual oracle failed"
    if (energy <= 0.0_dp) error stop "anisotropic energy oracle failed"

    call write_solution_csv()
    call write_diagnostics(maximum_error, residual, energy, elapsed)
    call render_solution()
    call record_gallery_stage("physical_solution")
    call render_slice()
    call render_diagnostics(maximum_error)
    call record_gallery_stage("diagnostics")

contains

    subroutine initialize_gallery_sequence()
        integer :: sequence_unit

        open (newunit=sequence_unit, file=output_directory// &
            "/gallery_sequence.txt", status="replace", action="write")
        close (sequence_unit)
    end subroutine initialize_gallery_sequence

    subroutine record_gallery_stage(stage)
        character(len=*), intent(in) :: stage
        integer :: sequence_unit

        open (newunit=sequence_unit, file=output_directory// &
            "/gallery_sequence.txt", status="old", position="append", &
            action="write")
        write (sequence_unit, "(a)") trim(stage)
        close (sequence_unit)
    end subroutine record_gallery_stage

    subroutine initialize_constitutive_tensor()
        call evaluate_field_aligned_constitutive_tensor( &
            k_parallel, k_perpendicular, field_direction, constitutive_tensor, &
            status)
        if (status%code /= 0) error stop "field-aligned tensor construction failed"
        if (maxval(abs(constitutive_tensor - transpose(constitutive_tensor))) > &
            1.0e-13_dp) error stop "anisotropic tensor lost symmetry"
    end subroutine initialize_constitutive_tensor

    subroutine initialize_mesh()
        integer :: ix, iy, iz, tetra_offset

        do iz = 0, cells_per_axis
            do iy = 0, cells_per_axis
                do ix = 0, cells_per_axis
                    node = node_index(ix, iy, iz)
                    coordinates(:, node) = [ &
                        real(ix, dp)/real(cells_per_axis, dp), &
                        real(iy, dp)/real(cells_per_axis, dp), &
                        real(iz, dp)/real(cells_per_axis, dp)]
                end do
            end do
        end do
        tetrahedron = 0
        do iz = 0, cells_per_axis - 1
            do iy = 0, cells_per_axis - 1
                do ix = 0, cells_per_axis - 1
                    cube = cube_node(ix, iy, iz)
                    first = cube
                    second = cube_node(ix + 1, iy, iz)
                    third = cube_node(ix + 1, iy + 1, iz)
                    fourth = cube_node(ix, iy + 1, iz)
                    fifth = cube_node(ix, iy, iz + 1)
                    sixth = cube_node(ix + 1, iy, iz + 1)
                    seventh = cube_node(ix + 1, iy + 1, iz + 1)
                    eighth = cube_node(ix, iy + 1, iz + 1)
                    nodes = [first, second, third, fourth, fifth, sixth, &
                        seventh, eighth]
                    tetra_offset = tetrahedron
                    tetrahedra(:, tetra_offset + 1) = nodes([1, 2, 3, 7])
                    tetrahedra(:, tetra_offset + 2) = nodes([1, 3, 4, 7])
                    tetrahedra(:, tetra_offset + 3) = nodes([1, 4, 8, 7])
                    tetrahedra(:, tetra_offset + 4) = nodes([1, 8, 5, 7])
                    tetrahedra(:, tetra_offset + 5) = nodes([1, 5, 6, 7])
                    tetrahedra(:, tetra_offset + 6) = nodes([1, 6, 2, 7])
                    tetrahedron = tetrahedron + tetrahedra_per_cube
                end do
            end do
        end do
        if (tetrahedron /= tetrahedron_count) &
            error stop "anisotropic mesh count failed"
    end subroutine initialize_mesh

    integer function node_index(ix, iy, iz)
        integer, intent(in) :: ix, iy, iz

        node_index = 1 + ix + nodes_per_axis*(iy + nodes_per_axis*iz)
    end function node_index

    integer function cube_node(ix, iy, iz)
        integer, intent(in) :: ix, iy, iz

        cube_node = node_index(ix, iy, iz)
    end function cube_node

    real(dp) function manufactured_solution(point_value)
        real(dp), intent(in) :: point_value(3)

        manufactured_solution = sin(pi*point_value(1))* &
            sin(pi*point_value(2))*sin(pi*point_value(3))
    end function manufactured_solution

    real(dp) function manufactured_source(point_value)
        real(dp), intent(in) :: point_value(3)
        real(dp) :: hessian(3, 3), value

        value = manufactured_solution(point_value)
        hessian = 0.0_dp
        hessian(1, 1) = -pi*pi*value
        hessian(2, 2) = -pi*pi*value
        hessian(3, 3) = -pi*pi*value
        hessian(1, 2) = pi*pi*cos(pi*point_value(1))* &
            cos(pi*point_value(2))*sin(pi*point_value(3))
        hessian(2, 1) = hessian(1, 2)
        hessian(1, 3) = pi*pi*cos(pi*point_value(1))* &
            sin(pi*point_value(2))*cos(pi*point_value(3))
        hessian(3, 1) = hessian(1, 3)
        hessian(2, 3) = pi*pi*sin(pi*point_value(1))* &
            cos(pi*point_value(2))*cos(pi*point_value(3))
        hessian(3, 2) = hessian(2, 3)
        manufactured_source = -sum(constitutive_tensor*hessian)
    end function manufactured_source

    subroutine impose_dirichlet_data()
        do node = 1, node_count
            if (minval([coordinates(1, node), coordinates(2, node), &
                coordinates(3, node), 1.0_dp - coordinates(1, node), &
                1.0_dp - coordinates(2, node), &
                1.0_dp - coordinates(3, node)]) < 1.0e-13_dp) then
                system_matrix(node, :) = 0.0_dp
                system_matrix(:, node) = 0.0_dp
                system_matrix(node, node) = 1.0_dp
                right_hand_side(node) = 0.0_dp
            end if
        end do
    end subroutine impose_dirichlet_data

    subroutine write_solution_csv()
        open (newunit=unit, file=output_directory//"/solution.csv", &
            status="replace", action="write")
        write (unit, "(a)") "node,x,y,z,computed,exact,absolute_error"
        do node = 1, node_count
            write (unit, "(i0,6(',',es24.16))") node, coordinates(:, node), &
                solution(node), exact_solution(node), nodal_error(node)
        end do
        close (unit)
    end subroutine write_solution_csv

    subroutine write_diagnostics(maximum_error_value, residual_value, &
            energy_value, elapsed_value)
        real(dp), intent(in) :: maximum_error_value, residual_value
        real(dp), intent(in) :: energy_value, elapsed_value

        open (newunit=unit, file=output_directory//"/diagnostics.csv", &
            status="replace", action="write")
        write (unit, "(a)") "quantity,value"
        write (unit, "(a,es24.16)") "k_parallel,", k_parallel
        write (unit, "(a,es24.16)") "k_perpendicular,", k_perpendicular
        write (unit, "(a,es24.16)") "anisotropy_ratio,", &
            k_parallel/k_perpendicular
        write (unit, "(a,es24.16)") "maximum_nodal_error,", &
            maximum_error_value
        write (unit, "(a,es24.16)") "linear_residual_inf,", residual_value
        write (unit, "(a,es24.16)") "discrete_energy,", energy_value
        write (unit, "(a,es24.16)") "assembly_and_solve_seconds,", &
            elapsed_value
        close (unit)
    end subroutine write_diagnostics

    subroutine render_solution()
        real(dp) :: arrow_x(node_count), arrow_y(node_count)
        real(dp) :: arrow_z(node_count), magnitude(node_count)
        real(dp) :: local_point(3)

        do node = 1, node_count
            node_x(node) = coordinates(1, node)
            node_y(node) = coordinates(2, node)
            node_z(node) = coordinates(3, node)
            local_point = coordinates(:, node)
            call analytical_flux(local_point, flux_at_node)
            magnitude(node) = sqrt(sum(flux_at_node**2))
            scale = 0.10_dp/max(1.0_dp, magnitude(node))
            arrow_x(node) = scale*flux_at_node(1)
            arrow_y(node) = scale*flux_at_node(2)
            arrow_z(node) = scale*flux_at_node(3)
        end do
        call figure(figsize=[8.5_dp, 7.0_dp])
        call add_scatter(node_x, node_y, node_z, c=solution, cmap="viridis", &
            marker="o", markersize=24.0_dp, label="computed P1 solution")
        call colorbar(label="u_h")
        do node = 1, node_count
            if (magnitude(node) > 1.0e-12_dp) call add_3d_plot( &
                [node_x(node), node_x(node) + arrow_x(node)], &
                [node_y(node), node_y(node) + arrow_y(node)], &
                [node_z(node), node_z(node) + arrow_z(node)], &
                color=arrow_color, linewidth=0.8_dp)
        end do
        call draw_box_edges()
        call xlabel("x")
        call ylabel("y")
        call title("Field-aligned anisotropic 3-D diffusion solution")
        call legend()
        call view_init(elev=24.0_dp, azim=-48.0_dp)
        call savefig(output_directory//"/solution.png")
        call savefig(output_directory//"/solution_3d.png")
    end subroutine render_solution

    subroutine draw_box_edges()
        integer :: line
        real(dp), parameter :: endpoints(3, 2, 12) = reshape([ &
            0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
            1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, &
            1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
            0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, &
            0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, &
            0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
            0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, &
            1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
            0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
            1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
            [3, 2, 12])

        do line = 1, 12
            call add_3d_plot(endpoints(1, :, line), endpoints(2, :, line), &
                endpoints(3, :, line), color=edge_color, linewidth=0.7_dp)
        end do
    end subroutine draw_box_edges

    subroutine render_slice()
        integer :: ix, iy, index, middle_plane
        real(dp) :: point_value(3), flux(3)

        middle_plane = cells_per_axis/2
        do iy = 0, cells_per_axis
            slice_y(iy + 1) = real(iy, dp)/real(cells_per_axis, dp)
            do ix = 0, cells_per_axis
                slice_x(ix + 1) = real(ix, dp)/real(cells_per_axis, dp)
                node = node_index(ix, iy, middle_plane)
                solution_slice(ix + 1, iy + 1) = solution(node)
                index = ix + 1 + nodes_per_axis*iy
                point_value = coordinates(:, node)
                call analytical_flux(point_value, flux)
                flux_x(index) = point_value(1)
                flux_y(index) = point_value(2)
                flux_u(index) = 0.16_dp*flux(1)/max(1.0_dp, sqrt(sum(flux**2)))
                flux_v(index) = 0.16_dp*flux(2)/max(1.0_dp, sqrt(sum(flux**2)))
            end do
        end do
        call figure(figsize=[8.0_dp, 6.5_dp])
        call contourf(slice_x, slice_y, solution_slice, cmap="viridis", &
            show_colorbar=.true.)
        call colorbar(label="u_h at z = 0.5")
        call quiver(flux_x, flux_y, flux_u, flux_v, scale=1.0_dp, &
            scale_units="xy", angles="xy", color="white", width=0.003_dp, &
            headwidth=3.0_dp)
        call xlabel("x")
        call ylabel("y")
        call title("Computed solution and field-aligned flux slice")
        call legend()
        call savefig(output_directory//"/solution_slice.png")
    end subroutine render_slice

    subroutine analytical_flux(point_value, flux)
        real(dp), intent(in) :: point_value(3)
        real(dp), intent(out) :: flux(3)
        real(dp) :: gradient_value(3)

        gradient_value = [ &
            pi*cos(pi*point_value(1))*sin(pi*point_value(2))* &
            sin(pi*point_value(3)), &
            pi*sin(pi*point_value(1))*cos(pi*point_value(2))* &
            sin(pi*point_value(3)), &
            pi*sin(pi*point_value(1))*sin(pi*point_value(2))* &
            cos(pi*point_value(3))]
        flux = -matmul(constitutive_tensor, gradient_value)
    end subroutine analytical_flux

    subroutine render_diagnostics(maximum_error_value)
        real(dp), intent(in) :: maximum_error_value
        real(dp) :: values(2)

        values = [maximum_error_value, residual]
        call figure(figsize=[7.0_dp, 5.0_dp])
        call contourf(slice_x, slice_y, solution_slice, cmap="viridis", &
            show_colorbar=.true.)
        call title("Anisotropic solution diagnostic slice")
        call xlabel("x")
        call ylabel("y")
        call savefig(output_directory//"/diagnostics_slice.png")
        if (any(values < 0.0_dp)) error stop "diagnostic values invalid"
    end subroutine render_diagnostics

end program field_aligned_anisotropic_3d_gallery
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/field_aligned_anisotropic_3d_gallery/primary.png)

### solution.png

![solution.png](../../media/examples/field_aligned_anisotropic_3d_gallery/solution.png)

---

[← Back to all examples](../index.html)
