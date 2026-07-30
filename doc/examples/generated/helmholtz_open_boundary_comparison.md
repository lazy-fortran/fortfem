---
title: helmholtz_open_boundary_comparison Example
---

# helmholtz_open_boundary_comparison Example

# Helmholtz open-boundary comparison

This example compares three truncations of the same one-dimensional outgoing
Helmholtz wave:

- the exact modal Dirichlet-to-Neumann condition \(u'=iku\);
- a quadratic complex-stretched perfectly matched layer;
- a substantially larger ordinary domain terminated by homogeneous
  Dirichlet data.

The analytical oracle is \(u(x)=e^{ikx}\). The larger-domain result
deliberately demonstrates that distance alone is not a nonreflecting
condition for Helmholtz: an undamped one-dimensional reflection retains unit
magnitude regardless of where the wall is placed.

CI generates `helmholtz_methods_1d.png`, `helmholtz_error_1d.png`, and
`benchmark.txt`. Generated images are not committed.

## Usage

```bash
fpm run --example helmholtz_open_boundary_comparison
```

## Source Code

```fortran
program helmholtz_open_boundary_comparison
    use fortfem_api, only: solve_scalar_helmholtz_pml_slab_1d
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve
    use fortplot, only: figure, legend, plot, savefig, title, xlabel, ylabel
    implicit none

    integer, parameter :: physical_nodes = 81
    integer, parameter :: pml_nodes = 161
    integer, parameter :: large_nodes = 297
    real(dp), parameter :: wave_number = 3.0_dp
    real(dp), parameter :: physical_end = 1.0_dp
    real(dp), parameter :: pml_end = 2.0_dp
    real(dp), parameter :: large_end = 3.7_dp
    real(dp), parameter :: sigma_max = 12.0_dp
    character(*), parameter :: output_directory = &
        "output/example/helmholtz_open_boundary_comparison"

    complex(dp), allocatable :: dtn(:), exact(:), large(:), pml(:)
    real(dp), allocatable :: large_grid(:), physical_grid(:), pml_grid(:)
    real(dp) :: dtn_error, dtn_seconds, end_time, large_error, large_seconds
    real(dp) :: pml_error, pml_seconds, start_time
    integer :: node, status, unit

    call execute_command_line("mkdir -p "//output_directory)
    call uniform_grid(physical_end, physical_nodes, physical_grid)
    call uniform_grid(pml_end, pml_nodes, pml_grid)
    call uniform_grid(large_end, large_nodes, large_grid)
    allocate(exact(physical_nodes))
    exact = exp(cmplx(0.0_dp, wave_number, dp)*physical_grid)

    call cpu_time(start_time)
    call solve_slab( &
        physical_grid, cmplx(0.0_dp, wave_number, dp), dtn, status)
    call cpu_time(end_time)
    if (status /= 0) error stop "slab DtN solve failed"
    dtn_seconds = end_time - start_time

    call cpu_time(start_time)
    call solve_scalar_helmholtz_pml_slab_1d( &
        pml_grid, physical_end, wave_number, sigma_max, 2, &
        cmplx(1.0_dp, 0.0_dp, dp), pml, status)
    call cpu_time(end_time)
    if (status /= 0) error stop "slab PML solve failed"
    pml_seconds = end_time - start_time

    call cpu_time(start_time)
    call solve_slab( &
        large_grid, cmplx(huge(1.0_dp), 0.0_dp, dp), large, status)
    call cpu_time(end_time)
    if (status /= 0) error stop "large-domain slab solve failed"
    large_seconds = end_time - start_time

    dtn_error = relative_error(dtn, exact)
    pml_error = relative_error(pml(:physical_nodes), exact)
    large_error = relative_error(large(:physical_nodes), exact)
    if (dtn_error >= 1.0e-3_dp) error stop "DtN accuracy regression"
    if (pml_error >= 2.0e-3_dp) error stop "PML accuracy regression"
    if (large_error <= 2.0e-1_dp) &
        error stop "large-domain reflection unexpectedly absent"

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(physical_grid, real(exact), label="analytical outgoing", &
        linestyle="-")
    call plot(physical_grid, real(dtn), label="exact DtN", linestyle="--")
    call plot(physical_grid, real(pml(:physical_nodes)), label="PML", &
        linestyle="-.")
    call plot(physical_grid, real(large(:physical_nodes)), &
        label="far Dirichlet boundary", linestyle=":")
    call xlabel("physical coordinate x")
    call ylabel("Re(u)")
    call title("Helmholtz open boundaries: DtN, PML, and larger domain")
    call legend()
    call savefig(output_directory//"/helmholtz_methods_1d.png")

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(physical_grid, abs(dtn - exact), label="DtN error")
    call plot(physical_grid, abs(pml(:physical_nodes) - exact), &
        label="PML error")
    call plot(physical_grid, abs(large(:physical_nodes) - exact), &
        label="far-boundary error")
    call xlabel("physical coordinate x")
    call ylabel("absolute error")
    call title("Open-boundary error in the physical slab")
    call legend()
    call savefig(output_directory//"/helmholtz_error_1d.png")

    open (newunit=unit, file=output_directory//"/benchmark.txt", &
        status="replace", action="write")
    write (unit, "(a,i0)") "physical nodes: ", physical_nodes
    write (unit, "(a,i0)") "PML nodes: ", pml_nodes
    write (unit, "(a,i0)") "larger-domain nodes: ", large_nodes
    write (unit, "(a,es14.6)") "DtN solve seconds: ", dtn_seconds
    write (unit, "(a,es14.6)") "PML solve seconds: ", pml_seconds
    write (unit, "(a,es14.6)") "larger-domain solve seconds: ", large_seconds
    write (unit, "(a,es14.6)") "DtN relative error: ", dtn_error
    write (unit, "(a,es14.6)") "PML relative error: ", pml_error
    write (unit, "(a,es14.6)") "larger-domain relative error: ", large_error
    close (unit)

contains

    subroutine uniform_grid(endpoint, count, grid)
        real(dp), intent(in) :: endpoint
        integer, intent(in) :: count
        real(dp), allocatable, intent(out) :: grid(:)

        integer :: index

        allocate(grid(count))
        do index = 1, count
            grid(index) = &
                endpoint*real(index - 1, dp)/real(count - 1, dp)
        end do
    end subroutine uniform_grid

    subroutine solve_slab(grid, right_dtn, solution, solve_status)
        real(dp), intent(in) :: grid(:)
        complex(dp), intent(in) :: right_dtn
        complex(dp), allocatable, intent(out) :: solution(:)
        integer, intent(out) :: solve_status

        complex(dp), allocatable :: matrix(:, :), rhs(:)
        complex(dp) :: local(2, 2)
        real(dp) :: length
        integer :: element, info

        allocate(matrix(size(grid), size(grid)), rhs(size(grid)))
        allocate(solution(size(grid)))
        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        rhs = cmplx(0.0_dp, 0.0_dp, dp)
        do element = 1, size(grid) - 1
            length = grid(element + 1) - grid(element)
            local = 1.0_dp/length*reshape( &
                [1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp], [2, 2]) - &
                wave_number**2*length/6.0_dp*reshape( &
                [2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
            matrix(element:element + 1, element:element + 1) = &
                matrix(element:element + 1, element:element + 1) + local
        end do
        if (abs(right_dtn) < 0.5_dp*huge(1.0_dp)) then
            matrix(size(grid), size(grid)) = &
                matrix(size(grid), size(grid)) - right_dtn
        else
            matrix(:, size(grid)) = cmplx(0.0_dp, 0.0_dp, dp)
            matrix(size(grid), :) = cmplx(0.0_dp, 0.0_dp, dp)
            matrix(size(grid), size(grid)) = cmplx(1.0_dp, 0.0_dp, dp)
        end if
        rhs = rhs - matrix(:, 1)
        matrix(:, 1) = cmplx(0.0_dp, 0.0_dp, dp)
        matrix(1, :) = cmplx(0.0_dp, 0.0_dp, dp)
        matrix(1, 1) = cmplx(1.0_dp, 0.0_dp, dp)
        rhs(1) = cmplx(1.0_dp, 0.0_dp, dp)
        call dense_solve(matrix, rhs, solution, info)
        solve_status = info
    end subroutine solve_slab

    pure real(dp) function relative_error(values, oracle) result(error)
        complex(dp), intent(in) :: values(:), oracle(:)

        error = sqrt(sum(abs(values - oracle)**2))/sqrt(sum(abs(oracle)**2))
    end function relative_error

end program helmholtz_open_boundary_comparison
```

## Generated Plots

*No plot artifact is produced by this example.*

---

[← Back to all examples](../index.html)
