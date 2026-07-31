---
title: laplace_symmetric_transmission Example
---

# laplace_symmetric_transmission Example

Manufactured transmission pair:

## Usage

```bash
fpm run --example laplace_symmetric_transmission
```

## Source Code

```fortran
program laplace_symmetric_transmission
    use fortfem_api, only: solve_laplace_symmetric_coupling_p1_p0
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, title, xlabel, ylabel
    implicit none

    integer, parameter :: points_per_side = 5
    integer, parameter :: panel_count = 4 * (points_per_side - 1)
    integer, parameter :: triangle_count = &
        2 * (points_per_side - 1)**2
    integer, parameter :: vertex_count = points_per_side**2
    real(dp), parameter :: side_length = 0.25_dp

    real(dp) :: dirichlet_jump(vertex_count), exterior_flux(panel_count)
    real(dp) :: interior_solution(vertex_count), neumann_jump(panel_count)
    real(dp) :: panel_end(2, panel_count), panel_start(2, panel_count)
    real(dp) :: vertices(2, vertex_count), volume_load(vertex_count)
    integer :: panel_nodes(2, panel_count), triangles(3, triangle_count)
    integer :: command_status, panel, status
    real(dp) :: error_flux, error_interior, length, normal(2)

    call execute_command_line( &
        "mkdir -p output/example/laplace_symmetric_transmission", &
        exitstat=command_status)
    if (command_status /= 0) error stop "cannot create example output directory"

    call build_square_mesh(vertices, triangles)
    call build_square_boundary(vertices, panel_nodes, panel_start, panel_end)

    ! Manufactured transmission pair:
    ! u_inside = 3/4 + x + 2 y, u_outside = 0.
    dirichlet_jump = 0.75_dp + vertices(1, :) + 2.0_dp * vertices(2, :)
    volume_load = 0.0_dp
    do panel = 1, panel_count
        length = norm2(panel_end(:, panel) - panel_start(:, panel))
        normal = [panel_end(2, panel) - panel_start(2, panel), &
            panel_start(1, panel) - panel_end(1, panel)] / length
        neumann_jump(panel) = normal(1) + 2.0_dp * normal(2)
    end do

    call solve_laplace_symmetric_coupling_p1_p0( &
        vertices, triangles, panel_start, panel_end, panel_nodes, 20, &
        volume_load, dirichlet_jump, neumann_jump, interior_solution, &
        exterior_flux, status)
    if (status /= 0) error stop "Laplace FEM-BEM transmission solve failed"

    error_interior = maxval(abs(interior_solution - dirichlet_jump))
    error_flux = maxval(abs(exterior_flux))
    print '(a,es12.4)', "maximum interior error: ", error_interior
    print '(a,es12.4)', "maximum exterior flux error: ", error_flux
    if (error_interior > 2.0e-11_dp .or. error_flux > 2.0e-11_dp) then
        error stop "Manufactured transmission solution was not recovered"
    end if

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(real([(panel, panel=1, vertex_count)], dp), &
        interior_solution, label="interior solution", marker="o")
    call plot(real([(panel, panel=1, vertex_count)], dp), &
        dirichlet_jump, label="manufactured solution", linestyle="--")
    call xlabel("vertex index")
    call ylabel("solution value")
    call title("Laplace symmetric FEM-BEM transmission solution")
    call legend()
    call savefig( &
        "output/example/laplace_symmetric_transmission/solution_1d.png")

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(real([(panel, panel=1, panel_count)], dp), exterior_flux, &
        label="exterior flux jump", marker="s")
    call xlabel("boundary panel index")
    call ylabel("flux jump")
    call title("Manufactured transmission flux residual")
    call savefig( &
        "output/example/laplace_symmetric_transmission/flux_1d.png")

contains

    subroutine build_square_mesh(vertices, triangles)
        real(dp), intent(out) :: vertices(2, vertex_count)
        integer, intent(out) :: triangles(3, triangle_count)

        real(dp) :: spacing
        integer :: column, lower_left, row, triangle, vertex

        spacing = side_length / real(points_per_side - 1, dp)
        vertex = 0
        do row = 0, points_per_side - 1
            do column = 0, points_per_side - 1
                vertex = vertex + 1
                vertices(:, vertex) = spacing * real([column, row], dp)
            end do
        end do

        triangle = 0
        do row = 1, points_per_side - 1
            do column = 1, points_per_side - 1
                lower_left = column + (row - 1) * points_per_side
                triangle = triangle + 1
                triangles(:, triangle) = [ &
                    lower_left, lower_left + 1, &
                    lower_left + points_per_side + 1]
                triangle = triangle + 1
                triangles(:, triangle) = [ &
                    lower_left, lower_left + points_per_side + 1, &
                    lower_left + points_per_side]
            end do
        end do
    end subroutine build_square_mesh

    subroutine build_square_boundary( &
            vertices, panel_nodes, panel_start, panel_end)
        real(dp), intent(in) :: vertices(2, vertex_count)
        integer, intent(out) :: panel_nodes(2, panel_count)
        real(dp), intent(out) :: panel_start(2, panel_count)
        real(dp), intent(out) :: panel_end(2, panel_count)

        integer :: boundary_nodes(panel_count)
        integer :: index, panel

        index = 0
        do panel = 1, points_per_side
            index = index + 1
            boundary_nodes(index) = panel
        end do
        do panel = 2, points_per_side
            index = index + 1
            boundary_nodes(index) = &
                panel * points_per_side
        end do
        do panel = points_per_side - 1, 1, -1
            index = index + 1
            boundary_nodes(index) = &
                (points_per_side - 1) * points_per_side + panel
        end do
        do panel = points_per_side - 1, 2, -1
            index = index + 1
            boundary_nodes(index) = &
                (panel - 1) * points_per_side + 1
        end do

        do panel = 1, panel_count
            panel_nodes(1, panel) = boundary_nodes(panel)
            panel_nodes(2, panel) = boundary_nodes(mod(panel, panel_count) + 1)
            panel_start(:, panel) = vertices(:, panel_nodes(1, panel))
            panel_end(:, panel) = vertices(:, panel_nodes(2, panel))
        end do
    end subroutine build_square_boundary

end program laplace_symmetric_transmission
```

## Generated Plots

### flux_1d.png

![flux_1d.png](../../media/examples/laplace_symmetric_transmission/flux_1d.png)

### solution_1d.png

![solution_1d.png](../../media/examples/laplace_symmetric_transmission/solution_1d.png)

---

[← Back to all examples](../index.html)
