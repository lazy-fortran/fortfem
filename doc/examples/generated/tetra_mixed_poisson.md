---
title: tetra_mixed_poisson Example
---

# tetra_mixed_poisson Example

# Tetrahedral symbolic mixed Poisson

This example extends the symbolic Darcy form from triangles to a unit cube
split into six tetrahedra:

\[
 (q,\tau) - (p,\nabla\!\cdot\tau)=0,\qquad
 (\nabla\!\cdot q,v)=(1,v).
\]

It runs every implemented tetrahedral pair from RT0/DG0 through RT5/DG5.
The RT mass and rectangular divergence expressions compile directly into
FortSparse CSC blocks. For each tetrahedron, the constant DG test moment of
the numerical divergence must equal the independently computed geometric
volume. A second plot records the coupled system size as order increases.

CI generates `tetra_mixed_conservation_1d.png` and
`tetra_mixed_dofs_1d.png`; generated media are not committed.

## Usage

```bash
fpm run --example tetra_mixed_poisson
```

## Source Code

```fortran
program tetra_mixed_poisson
    use fortfem_api, only: &
        compile_tetra_mixed_form_csc, div, dx, &
        evaluate_tetra_discontinuous, initialize_tetra_discontinuous, &
        generate_structured_tetra_box_mesh, init_measures, inner, operator(*), &
        solve_symbolic_tetra_mixed_poisson_rt, test_function_t, &
        tetra_discontinuous_t, trial_function_t, vector_test_function_t, &
        vector_trial_function_t
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det3
    use fortplot, only: add_scatter, figure, legend, plot, savefig, set_yscale, &
        title, xlabel, ylabel
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    character(*), parameter :: output_directory = &
        "output/example/tetra_mixed_poisson"
    type(csc_t) :: divergence
    type(fortsparse_status_t) :: status
    type(test_function_t) :: pressure_test
    type(trial_function_t) :: pressure_trial
    type(vector_test_function_t) :: flux_test
    type(vector_trial_function_t) :: flux_trial
    integer, allocatable :: tetrahedra(:, :)
    real(dp), allocatable :: balance(:), flux(:), pressure(:), vertices(:, :)
    real(dp), allocatable :: pressure_plot(:)
    real(dp) :: balance_error(6), bounds(3, 2), cell_balance
    real(dp) :: dofs(6), jacobian(3, 3), orders(6), volume
    integer :: degree, dg_count, local_status, tetrahedron

    call init_measures()
    call execute_command_line("mkdir -p "//output_directory)
    bounds(:, 1) = 0.0_dp
    bounds(:, 2) = 1.0_dp
    call generate_structured_tetra_box_mesh( &
        bounds, [1, 1, 1], vertices, tetrahedra, local_status)
    if (local_status /= 0) error stop "tetrahedral box mesh failed"

    do degree = 0, 5
        orders(degree + 1) = real(degree, dp)
        call solve_symbolic_tetra_mixed_poisson_rt( &
            vertices, tetrahedra, degree, 2*degree + 4, &
            inner(flux_trial, flux_test)*dx, &
            (-1.0_dp)*inner(pressure_trial, div(flux_test))*dx, &
            inner(div(flux_trial), pressure_test)*dx, unit_source, &
            flux, pressure, status)
        if (status%code /= 0) error stop "symbolic tetra solve failed"
        call compile_tetra_mixed_form_csc( &
            inner(div(flux_trial), pressure_test)*dx, vertices, tetrahedra, &
            degree, 2*degree + 4, divergence, status)
        if (status%code /= 0) error stop "symbolic tetra block failed"
        balance = csc_matvec(divergence, flux)
        dg_count = (degree + 1)*(degree + 2)*(degree + 3)/6
        balance_error(degree + 1) = 0.0_dp
        do tetrahedron = 1, size(tetrahedra, 2)
            jacobian(:, 1) = &
                vertices(:, tetrahedra(2, tetrahedron)) - &
                vertices(:, tetrahedra(1, tetrahedron))
            jacobian(:, 2) = &
                vertices(:, tetrahedra(3, tetrahedron)) - &
                vertices(:, tetrahedra(1, tetrahedron))
            jacobian(:, 3) = &
                vertices(:, tetrahedra(4, tetrahedron)) - &
                vertices(:, tetrahedra(1, tetrahedron))
            volume = abs(det3(jacobian))/6.0_dp
            cell_balance = balance((tetrahedron - 1)*dg_count + 1)
            balance_error(degree + 1) = max( &
                balance_error(degree + 1), abs(cell_balance - volume))
        end do
        dofs(degree + 1) = real(size(flux) + size(pressure), dp)
        if (balance_error(degree + 1) >= 2.0e-9_dp) &
            error stop "tetrahedral cell balance regression"
        if (degree == 5) then
            allocate(pressure_plot(size(pressure)))
            pressure_plot = pressure
        end if
        deallocate(balance, flux, pressure)
    end do

    call render_solution(pressure_plot)

    call figure(figsize=[8.5_dp, 5.5_dp])
    call plot(orders, max(balance_error, epsilon(1.0_dp)), &
        label="maximum cell balance error", marker="o")
    call set_yscale("log")
    call xlabel("RT-DG polynomial degree")
    call ylabel("absolute balance error")
    call title("Symbolic tetrahedral mixed Poisson conservation")
    call legend()
    call savefig(output_directory//"/tetra_mixed_conservation_1d.png")

    call figure(figsize=[8.5_dp, 5.5_dp])
    call plot(orders, dofs, label="coupled unknowns", marker="s")
    call xlabel("RT-DG polynomial degree")
    call ylabel("degrees of freedom")
    call title("Tetrahedral mixed-system order growth")
    call legend()
    call savefig(output_directory//"/tetra_mixed_dofs_1d.png")

contains

    subroutine render_solution(pressure_coefficients)
        real(dp), intent(in) :: pressure_coefficients(:)

        integer, parameter :: sample_side = 16
        type(tetra_discontinuous_t) :: basis
        real(dp), allocatable :: x_plot(:), y_plot(:), z_plot(:), values(:)
        real(dp), allocatable :: basis_values(:)
        real(dp) :: point(3), reference_point(3), value
        integer :: cell, count, i, j, k, local_status

        call initialize_tetra_discontinuous(5, basis, local_status)
        if (local_status /= 0) error stop "tetrahedral plot basis failed"
        if (size(pressure_coefficients) /= &
            size(tetrahedra, 2)*size(basis%exponents, 2)) then
            error stop "tetrahedral pressure plot shape mismatch"
        end if
        allocate( &
            x_plot(size(tetrahedra, 2)*(sample_side + 1)**3), &
            y_plot(size(tetrahedra, 2)*(sample_side + 1)**3), &
            z_plot(size(tetrahedra, 2)*(sample_side + 1)**3), &
            values(size(tetrahedra, 2)*(sample_side + 1)**3), &
            basis_values(size(basis%exponents, 2)))
        count = 0
        do cell = 1, size(tetrahedra, 2)
            do k = 0, sample_side
                do j = 0, sample_side
                    do i = 0, sample_side
                        if (i + j + k > sample_side) cycle
                        reference_point = [ &
                            real(i, dp)/real(sample_side, dp), &
                            real(j, dp)/real(sample_side, dp), &
                            real(k, dp)/real(sample_side, dp)]
                        call evaluate_tetra_discontinuous( &
                            basis, reference_point(1), reference_point(2), &
                            reference_point(3), basis_values, local_status)
                        if (local_status /= 0) cycle
                        value = dot_product(basis_values, pressure_coefficients( &
                            (cell - 1)*size(basis%exponents, 2) + 1: &
                            cell*size(basis%exponents, 2)))
                        point = vertices(:, tetrahedra(1, cell)) + &
                            reference_point(1)*(vertices(:, tetrahedra(2, cell)) &
                            - vertices(:, tetrahedra(1, cell))) + &
                            reference_point(2)*(vertices(:, tetrahedra(3, cell)) &
                            - vertices(:, tetrahedra(1, cell))) + &
                            reference_point(3)*(vertices(:, tetrahedra(4, cell)) &
                            - vertices(:, tetrahedra(1, cell)))
                        count = count + 1
                        x_plot(count) = point(1)
                        y_plot(count) = point(2)
                        z_plot(count) = point(3)
                        values(count) = value
                    end do
                end do
            end do
        end do
        call figure(figsize=[7.5_dp, 6.0_dp])
        call add_scatter( &
            x_plot(:count), y_plot(:count), z_plot(:count), &
            c=values(:count), cmap="viridis", marker=".", &
            markersize=4.0_dp, label="computed DG pressure samples")
        call title("RT-DG mixed Poisson computed pressure")
        call savefig(output_directory//"/tetra_mixed_solution_3d.png")
    end subroutine render_solution

    pure real(dp) function unit_source(x, y, z) result(value)
        real(dp), intent(in) :: x, y, z

        associate(unused => x + y + z)
        end associate
        value = 1.0_dp
    end function unit_source

end program tetra_mixed_poisson
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/tetra_mixed_poisson/primary.png)

### tetra_mixed_conservation_1d.png

![tetra_mixed_conservation_1d.png](../../media/examples/tetra_mixed_poisson/tetra_mixed_conservation_1d.png)

### tetra_mixed_dofs_1d.png

![tetra_mixed_dofs_1d.png](../../media/examples/tetra_mixed_poisson/tetra_mixed_dofs_1d.png)

### tetra_mixed_solution_3d.png

![tetra_mixed_solution_3d.png](../../media/examples/tetra_mixed_poisson/tetra_mixed_solution_3d.png)

---

[← Back to all examples](../index.html)
