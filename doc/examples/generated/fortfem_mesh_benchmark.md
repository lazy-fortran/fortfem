---
title: fortfem_mesh_benchmark Example
---

# fortfem_mesh_benchmark Example

# FortFEM Mesh Benchmark

This example benchmarks FortFEM's mesh generation performance against FreeFEM.

## Description

Compares mesh generation speed and quality between FortFEM and FreeFEM implementations for various mesh sizes and geometries.

## Usage

```bash
fpm run --example fortfem_mesh_benchmark
```

## What it does

- Generates meshes of different sizes using FortFEM
- Measures timing and mesh quality metrics
- Outputs benchmark results for comparison with FreeFEM
## Usage

```bash
fpm run --example fortfem_mesh_benchmark
```

## Source Code

```fortran
program fortfem_mesh_benchmark
    ! FortFEM mesh generation benchmark to compare with FreeFEM
    use fortfem_api
    use fortfem_kinds
    use fortplot, only: figure, legend, fortplot_plot => plot, savefig, &
        set_xscale, set_yscale, title, xlabel, ylabel
    implicit none

    character(*), parameter :: output_directory = &
        "output/example/fortfem_mesh_benchmark"
    integer, parameter :: n_sizes = 6
    integer :: mesh_sizes(n_sizes) = [5, 10, 15, 20, 25, 30]
    real(dp) :: fortfem_times(n_sizes)
    integer :: fortfem_vertices(n_sizes), fortfem_triangles(n_sizes)
    real(dp) :: mesh_sizes_real(n_sizes), time_plot(n_sizes)
    real(dp) :: vertices_plot(n_sizes), triangles_plot(n_sizes)

    type(mesh_t) :: mesh
    real(dp) :: start_time, end_time
    integer :: command_status, i, n

    call execute_command_line( &
        "mkdir -p "//output_directory, exitstat=command_status)
    if (command_status /= 0) error stop "cannot create example output directory"

    write(*,*) "=== FortFEM Mesh Generation Benchmark ==="
    write(*,*) ""

    ! Run benchmark for each mesh size
    do i = 1, n_sizes
        n = mesh_sizes(i)

        write(*,'(A,I0,A,I0,A)') "Testing ", n, "x", n, " mesh:"

        ! Time FortFEM mesh generation
        call cpu_time(start_time)
        mesh = unit_square_mesh(n)
        call cpu_time(end_time)

        fortfem_times(i) = end_time - start_time
        fortfem_vertices(i) = mesh%data%n_vertices
        fortfem_triangles(i) = mesh%data%n_triangles
        mesh_sizes_real(i) = real(n, dp)

        write(*,'(A,I0,A,I0,A,ES10.3,A)') "  FortFEM: ", fortfem_vertices(i), &
            " vertices, ", fortfem_triangles(i), " triangles, ", fortfem_times(i), "s"
    end do

    write(*,*) ""
    write(*,*) "=== Benchmark Results ==="
    write(*,*) "Size	FortFEM_vertices	FortFEM_triangles	FortFEM_time(s)"

    do i = 1, n_sizes
        write(*,'(I0,A,I0,A,I0,A,ES10.3)') mesh_sizes(i), "	", fortfem_vertices(i), &
            "		", fortfem_triangles(i), "		", fortfem_times(i)
    end do

    call plot(mesh, filename=output_directory//"/representative_mesh_2d.png", &
        title="Representative FortFEM benchmark mesh")

    ! Save results for comparison
    call save_benchmark_results()

    time_plot = max(fortfem_times, tiny(1.0_dp))
    vertices_plot = real(fortfem_vertices, dp)
    triangles_plot = real(fortfem_triangles, dp)
    call figure(figsize=[9.0_dp, 5.5_dp])
    call fortplot_plot(mesh_sizes_real, time_plot, &
        label="FortFEM mesh generation", marker="o")
    call set_xscale("log")
    call set_yscale("log")
    call xlabel("mesh resolution n")
    call ylabel("generation time [s]")
    call title("FortFEM mesh-generation scaling")
    call legend()
    call savefig(output_directory//"/generation_time_1d.png")

    call figure(figsize=[9.0_dp, 5.5_dp])
    call fortplot_plot(mesh_sizes_real, vertices_plot, &
        label="vertices", marker="o")
    call fortplot_plot(mesh_sizes_real, triangles_plot, &
        label="triangles", marker="s")
    call set_xscale("log")
    call set_yscale("log")
    call xlabel("mesh resolution n")
    call ylabel("entity count")
    call title("FortFEM mesh-size growth")
    call legend()
    call savefig(output_directory//"/entity_counts_1d.png")

    write(*,*) ""
    write(*,*) "FortFEM benchmark completed."
    write(*,*) "Results saved to fortfem_benchmark_results.dat"
    write(*,*) ""
    write(*,*) "To compare with FreeFEM:"
    write(*,*) "1. Run: FreeFem++ benchmark/mesh_comparison.edp"
    write(*,*) "2. Compare fortfem_benchmark_results.dat with freefem_benchmark_results.dat"

contains

    subroutine save_benchmark_results()
        integer :: unit_num, ios

        open(newunit=unit_num, file="fortfem_benchmark_results.dat", &
            status="replace", action="write", iostat=ios)

        if (ios /= 0) then
            write(*,*) "Warning: Could not save benchmark results"
            return
        end if

        write(unit_num, '(A)') "# mesh_size vertices triangles time_seconds"

        do i = 1, n_sizes
            write(unit_num, '(I0,1X,I0,1X,I0,1X,ES15.8)') mesh_sizes(i), &
                fortfem_vertices(i), fortfem_triangles(i), fortfem_times(i)
        end do

        close(unit_num)
    end subroutine

end program fortfem_mesh_benchmark
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/fortfem_mesh_benchmark/primary.png)

### entity_counts_1d.png

![entity_counts_1d.png](../../media/examples/fortfem_mesh_benchmark/entity_counts_1d.png)

### generation_time_1d.png

![generation_time_1d.png](../../media/examples/fortfem_mesh_benchmark/generation_time_1d.png)

### representative_mesh_2d.png

![representative_mesh_2d.png](../../media/examples/fortfem_mesh_benchmark/representative_mesh_2d.png)

---

[← Back to all examples](../index.html)
