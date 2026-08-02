---
title: curl_curl Example
---

# curl_curl Example

# Curl-Curl Electromagnetic Example

This example solves a manufactured three-dimensional curl-curl problem with
FortFEM's verified, arbitrary-order tetrahedral Nédélec elements and the
FortSparse direct solver.

## Problem

On the unit cube, solve

```text
curl(curl(E)) + E = J,
n × E = 0 on the boundary,
```

with

```text
E = (0, 0, sin(pi x) sin(pi y)),
J = (0, 0, (2 pi² + 1) sin(pi x) sin(pi y)).
```

Both the field and its analytical curl are sampled after covariant Piola
mapping. Their errors must decrease for Nédélec orders one through five; the
program stops with an error if that independent convergence oracle fails.

## Numerics

- first-kind tetrahedral Nédélec spaces of orders 1–5;
- globally oriented edge, face, and volume degrees of freedom;
- sparse curl-curl-plus-mass assembly and solution through FortSparse;
- homogeneous PEC tangential boundary elimination.

## Output

- `curl_curl_field_slice_2d.png`: the order-five solved Nédélec field
  magnitude with the reconstructed in-plane `curl(E_h)` arrows on a physical
  mid-plane; this is the primary visual result.
- `curl_curl_p_convergence.png`: field and curl errors by polynomial order;
- `convergence.csv`: the numerical values used to generate the plot.

Generated output is written under `output/example/curl_curl/` and is not
checked into the repository.

## Usage

```bash
fpm run --example curl_curl
```

## Source Code

```fortran
program curl_curl_example
    !! Verified arbitrary-order tetrahedral Nedelec curl-curl solve.
    use fortfem_feec, only: solve_tetra_nedelec_curl_mass, &
        evaluate_tetra_nedelec_interpolant_at_point
    use fortfem_core, only: invert_tetra_affine_map
    use fortfem_kinds, only: dp
    use fortfem_tetra_nedelec_arbitrary_order, only: &
        evaluate_tetra_nedelec_first_kind, &
        initialize_tetra_nedelec_first_kind, tetra_nedelec_first_kind_t
    use fortfem_tetra_nedelec_global_dof_map, only: &
        build_tetra_nedelec_basis_transform, build_tetra_nedelec_dof_map
    use fortfem_tetra_piola_maps, only: map_tetra_nedelec_covariant
    use fortnum_linalg, only: det3
    use fortplot, only: add_scatter, colorbar, figure, legend, plot, quiver, &
        savefig, title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: maximum_order = 5
    character(*), parameter :: output_directory = "output/example/curl_curl"
    type(fortsparse_status_t) :: status
    real(dp), allocatable :: solution(:)
    real(dp) :: errors(2, maximum_order), orders(maximum_order)
    real(dp) :: vertices(3, 8)
    integer :: order, tetrahedra(4, 6), unit

    call execute_command_line("mkdir -p "//output_directory)
    call build_cube_mesh(vertices, tetrahedra)
    do order = 1, maximum_order
        call solve_tetra_nedelec_curl_mass( &
            vertices, tetrahedra, order, manufactured_source, 1.0_dp, &
            1.0_dp, solution, status, zero_tangential_boundary=.true.)
        if (status%code /= 0) error stop "Nedelec curl-curl solve failed"
        call measure_error( &
            vertices, tetrahedra, order, solution, errors(:, order))
        orders(order) = real(order, dp)
        write (*, "(a,i0,a,2(es12.4,1x))") &
            "Nedelec order ", order, " field/curl errors ", errors(:, order)
        if (order < maximum_order) deallocate(solution)
    end do
    if (.not. all(errors(:, 2:) < errors(:, :maximum_order - 1))) &
        error stop "curl-curl p-convergence regression"

    call render_field_slice(solution)
    deallocate(solution)

    call figure(figsize=[8.0_dp, 5.0_dp])
    call plot(orders, errors(1, :), label="field error", marker="o")
    call plot(orders, errors(2, :), label="curl error", marker="s")
    call xlabel("Nedelec polynomial order")
    call ylabel("maximum sampled error")
    call title("Manufactured three-dimensional curl-curl convergence")
    call legend()
    call savefig(output_directory//"/curl_curl_p_convergence.png")

    open (newunit=unit, file=output_directory//"/convergence.csv", &
        status="replace", action="write")
    write (unit, "(a)") "order,field_error,curl_error"
    do order = 1, maximum_order
        write (unit, "(*(es24.16,:,','))") orders(order), errors(:, order)
    end do
    close (unit)

contains

    subroutine render_field_slice(coefficients)
        integer, parameter :: nx = 17, ny = 17
        real(dp), intent(in) :: coefficients(:)
        real(dp) :: x_grid(nx), y_grid(ny), u_grid(nx, ny), v_grid(nx, ny)
        real(dp) :: x_flat(nx*ny), y_flat(nx*ny), u_flat(nx*ny), v_flat(nx*ny)
        real(dp) :: magnitude(nx*ny), point(3), value(3), curl_value(3)
        real(dp), allocatable :: basis_transform(:, :), local_dofs(:)
        integer, allocatable :: edge_orientations(:, :), edges(:, :)
        integer, allocatable :: face_permutations(:, :, :), faces(:, :)
        integer, allocatable :: global_dofs(:, :)
        type(tetra_nedelec_first_kind_t) :: basis
        integer :: dof_count, local_status
        integer :: i, j, index

        call initialize_tetra_nedelec_first_kind( &
            maximum_order, basis, local_status)
        if (local_status /= 0) error stop "Nedelec plot basis failed"
        call build_tetra_nedelec_dof_map( &
            maximum_order, tetrahedra, edges, faces, global_dofs, &
            edge_orientations, face_permutations, local_status)
        if (local_status /= 0) error stop "Nedelec plot DOF map failed"
        dof_count = size(global_dofs, 1)
        allocate( &
            basis_transform(dof_count, dof_count), local_dofs(dof_count))

        do i = 1, nx
            x_grid(i) = real(i - 1, dp)/real(nx - 1, dp)
        end do
        do j = 1, ny
            y_grid(j) = real(j - 1, dp)/real(ny - 1, dp)
        end do
        do j = 1, ny
            do i = 1, nx
                point = [x_grid(i), y_grid(j), 0.5_dp]
                call evaluate_global_field( &
                    point, coefficients, basis, global_dofs, &
                    edge_orientations, face_permutations, basis_transform, &
                    local_dofs, value, curl_value, local_status)
                if (local_status /= 0) error stop "Nedelec plot field failed"
                magnitude(index_point(i, j, nx)) = norm2(value)
                u_grid(i, j) = curl_value(1)
                v_grid(i, j) = curl_value(2)
            end do
        end do

        index = 0
        do j = 1, ny
            do i = 1, nx
                index = index + 1
                x_flat(index) = x_grid(i)
                y_flat(index) = y_grid(j)
                u_flat(index) = u_grid(i, j)
                v_flat(index) = v_grid(i, j)
            end do
        end do
        call figure(figsize=[7.5_dp, 6.0_dp])
        call add_scatter( &
            x_flat, y_flat, c=magnitude, cmap="viridis", marker=".", &
            markersize=8.0_dp, label="solved Nedelec field magnitude")
        call quiver( &
            x_flat, y_flat, u_flat, v_flat, scale=1.5_dp, &
            color="black", width=0.0035_dp)
        call colorbar(label="|E_h|")
        call xlabel("x")
        call ylabel("y")
        call title("Solved curl-curl field and curl(E_h) at z = 0.5")
        call savefig(output_directory//"/curl_curl_field_slice_2d.png")
    end subroutine render_field_slice

    subroutine evaluate_global_field( &
            point, coefficients, basis, global_dofs, edge_orientations, &
            face_permutations, basis_transform, local_dofs, value, curl_value, &
            local_status)
        real(dp), intent(in) :: point(3), coefficients(:)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        integer, intent(in) :: global_dofs(:, :), edge_orientations(:, :)
        integer, intent(in) :: face_permutations(:, :, :)
        real(dp), intent(inout) :: basis_transform(:, :), local_dofs(:)
        real(dp), intent(out) :: value(3), curl_value(3)
        integer, intent(out) :: local_status

        real(dp) :: reference(3)
        integer :: cell

        value = 0.0_dp
        curl_value = 0.0_dp
        local_status = 1
        do cell = 1, size(tetrahedra, 2)
            call invert_tetra_affine_map( &
                vertices(:, tetrahedra(:, cell)), point, reference, &
                local_status)
            if (local_status /= 0) cycle
            if (minval(reference) < -1.0e-10_dp .or. &
                sum(reference) > 1.0_dp + 1.0e-10_dp) cycle
            call build_tetra_nedelec_basis_transform( &
                maximum_order, edge_orientations(:, cell), &
                face_permutations(:, :, cell), basis_transform, local_status)
            if (local_status /= 0) return
            local_dofs = matmul( &
                basis_transform, coefficients(global_dofs(:, cell)))
            call evaluate_tetra_nedelec_interpolant_at_point( &
                vertices(:, tetrahedra(:, cell)), basis, local_dofs, point, &
                value, curl_value, local_status)
            return
        end do
    end subroutine evaluate_global_field

    pure integer function index_point(i, j, width) result(index)
        integer, intent(in) :: i, j, width

        index = (j - 1)*width + i
    end function index_point

    pure subroutine manufactured_source(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value(3)
        real(dp), parameter :: pi = acos(-1.0_dp)

        associate(unused => z)
        end associate
        value = [ &
            0.0_dp, 0.0_dp, &
            (2.0_dp*pi*pi + 1.0_dp)*sin(pi*x)*sin(pi*y)]
    end subroutine manufactured_source

    pure function exact_field(point) result(value)
        real(dp), intent(in) :: point(3)
        real(dp) :: value(3)
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = [ &
            0.0_dp, 0.0_dp, &
            sin(pi*point(1))*sin(pi*point(2))]
    end function exact_field

    pure function exact_curl(point) result(value)
        real(dp), intent(in) :: point(3)
        real(dp) :: value(3)
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = [ &
            pi*sin(pi*point(1))*cos(pi*point(2)), &
            -pi*cos(pi*point(1))*sin(pi*point(2)), 0.0_dp]
    end function exact_curl

    subroutine build_cube_mesh(mesh_vertices, cells)
        real(dp), intent(out) :: mesh_vertices(3, 8)
        integer, intent(out) :: cells(4, 6)

        real(dp) :: jacobian(3, 3)
        integer :: cell, temporary

        mesh_vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
        mesh_vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
        mesh_vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
        mesh_vertices(:, 4) = [1.0_dp, 1.0_dp, 0.0_dp]
        mesh_vertices(:, 5) = [0.0_dp, 0.0_dp, 1.0_dp]
        mesh_vertices(:, 6) = [1.0_dp, 0.0_dp, 1.0_dp]
        mesh_vertices(:, 7) = [0.0_dp, 1.0_dp, 1.0_dp]
        mesh_vertices(:, 8) = [1.0_dp, 1.0_dp, 1.0_dp]
        cells(:, 1) = [1, 2, 4, 8]
        cells(:, 2) = [1, 2, 6, 8]
        cells(:, 3) = [1, 3, 4, 8]
        cells(:, 4) = [1, 3, 7, 8]
        cells(:, 5) = [1, 5, 6, 8]
        cells(:, 6) = [1, 5, 7, 8]
        do cell = 1, size(cells, 2)
            jacobian(:, 1) = &
                mesh_vertices(:, cells(2, cell)) - &
                mesh_vertices(:, cells(1, cell))
            jacobian(:, 2) = &
                mesh_vertices(:, cells(3, cell)) - &
                mesh_vertices(:, cells(1, cell))
            jacobian(:, 3) = &
                mesh_vertices(:, cells(4, cell)) - &
                mesh_vertices(:, cells(1, cell))
            if (det3(jacobian) < 0.0_dp) then
                temporary = cells(3, cell)
                cells(3, cell) = cells(4, cell)
                cells(4, cell) = temporary
            end if
        end do
    end subroutine build_cube_mesh

    subroutine measure_error( &
            mesh_vertices, cells, basis_order, coefficients, error)
        real(dp), intent(in) :: mesh_vertices(:, :), coefficients(:)
        integer, intent(in) :: cells(:, :), basis_order
        real(dp), intent(out) :: error(2)

        type(tetra_nedelec_first_kind_t) :: basis
        integer, allocatable :: edge_orientations(:, :), edges(:, :)
        integer, allocatable :: face_permutations(:, :, :), faces(:, :)
        integer, allocatable :: global_dofs(:, :)
        real(dp), allocatable :: basis_transform(:, :), local_dofs(:)
        real(dp), allocatable :: physical_curls(:, :), physical_values(:, :)
        real(dp), allocatable :: reference_curls(:, :)
        real(dp), allocatable :: reference_values(:, :)
        real(dp) :: curl_value(3), jacobian(3, 3), point(3)
        real(dp) :: reference_point(3), value(3)
        integer :: dof_count, local_status, tetrahedron

        call build_tetra_nedelec_dof_map( &
            basis_order, cells, edges, faces, global_dofs, edge_orientations, &
            face_permutations, local_status)
        if (local_status /= 0) error stop "Nedelec DOF map failed"
        call initialize_tetra_nedelec_first_kind( &
            basis_order, basis, local_status)
        if (local_status /= 0) error stop "Nedelec basis failed"
        dof_count = size(global_dofs, 1)
        allocate( &
            basis_transform(dof_count, dof_count), local_dofs(dof_count), &
            reference_values(3, dof_count), reference_curls(3, dof_count), &
            physical_values(3, dof_count), physical_curls(3, dof_count))
        reference_point = [0.21_dp, 0.24_dp, 0.19_dp]
        error = 0.0_dp
        do tetrahedron = 1, size(cells, 2)
            call build_tetra_nedelec_basis_transform( &
                basis_order, edge_orientations(:, tetrahedron), &
                face_permutations(:, :, tetrahedron), basis_transform, &
                local_status)
            if (local_status /= 0) error stop "Nedelec transform failed"
            local_dofs = matmul( &
                basis_transform, coefficients(global_dofs(:, tetrahedron)))
            jacobian(:, 1) = &
                mesh_vertices(:, cells(2, tetrahedron)) - &
                mesh_vertices(:, cells(1, tetrahedron))
            jacobian(:, 2) = &
                mesh_vertices(:, cells(3, tetrahedron)) - &
                mesh_vertices(:, cells(1, tetrahedron))
            jacobian(:, 3) = &
                mesh_vertices(:, cells(4, tetrahedron)) - &
                mesh_vertices(:, cells(1, tetrahedron))
            point = mesh_vertices(:, cells(1, tetrahedron)) + &
                matmul(jacobian, reference_point)
            call evaluate_tetra_nedelec_first_kind( &
                basis, reference_point, reference_values, reference_curls, &
                local_status)
            if (local_status /= 0) error stop "Nedelec evaluation failed"
            call map_tetra_nedelec_covariant( &
                jacobian, reference_values, reference_curls, physical_values, &
                physical_curls, local_status)
            if (local_status /= 0) error stop "Nedelec Piola map failed"
            value = matmul(physical_values, local_dofs)
            curl_value = matmul(physical_curls, local_dofs)
            error(1) = max(error(1), maxval(abs(value - exact_field(point))))
            error(2) = max( &
                error(2), maxval(abs(curl_value - exact_curl(point))))
        end do
    end subroutine measure_error

end program curl_curl_example
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/curl_curl/primary.png)

### curl_curl_field_slice_2d.png

![curl_curl_field_slice_2d.png](../../media/examples/curl_curl/curl_curl_field_slice_2d.png)

### curl_curl_p_convergence.png

![curl_curl_p_convergence.png](../../media/examples/curl_curl/curl_curl_p_convergence.png)

---

[← Back to all examples](../index.html)
