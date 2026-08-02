---
title: biro_tree_cotree_benchmark Example
---

# biro_tree_cotree_benchmark Example

# Bíró tree--cotree curl--curl benchmark

This small manufactured case reproduces the structural part of the direct
tree--cotree gauge described by Bíró, Preis, and Richter, *On the use of the
magnetic vector potential in the nodal and edge finite element analysis of 3D
magnetostatic fields*, IEEE Transactions on Magnetics 32 (1996), 651--654
([DOI](https://doi.org/10.1109/20.497322)). It uses a cyclic graph-level
topology analogue, fixes the spanning-tree coefficients, forms a unit
topological curl--curl Gram block, solves the reduced cotree system directly,
and reconstructs the physical edge potential. The paper is a formulation
paper; its application-specific 3D geometry, metric, and source data are not
redistributed here.

The first output is `solution.png`: the graph and oriented edge-potential
arrows. `solution.csv` and `diagnostics.csv` are generated under
`output/example/biro_tree_cotree_benchmark`; no media are checked in. The
matrix and right-hand side are manufactured, so this is a reproducible
method/orientation oracle rather than a claim to ship the paper's source or
benchmark geometry.

Run with:

```text
fo run --example biro_tree_cotree_benchmark
```

## Usage

```bash
fpm run --example biro_tree_cotree_benchmark
```

## Source Code

```fortran
program biro_tree_cotree_benchmark
    !! Small, reproducible tree--cotree curl--curl fixture.
    !!
    !! This is a license-safe manufactured reproduction of the direct-gauge
    !! topology used by Biro, Preis, and Richter.  It intentionally exposes
    !! the physical edge potential and its loop-current proxy before any
    !! residual or convergence diagnostics.
    use fortfem_feec, only: apply_tree_cotree_prolongation, &
        apply_tree_cotree_restriction, build_tree_cotree_gauge, &
        reduce_tree_cotree_dense_system, tree_cotree_gauge_t
    use fortfem_kinds, only: dp
    use fortplot, only: add_scatter, figure, legend, plot, quiver, savefig, &
        title, xlabel, ylabel
    implicit none

    integer, parameter :: edge_count = 5
    integer, parameter :: incidence(4, edge_count) = reshape([ &
        -1, 1, 0, 0, 0, -1, 1, 0, 0, 0, -1, 1, 1, 0, 0, -1, &
        -1, 0, 1, 0], [4, edge_count])
    integer, parameter :: edge_nodes(2, edge_count) = reshape([ &
        1, 2, 2, 3, 3, 4, 4, 1, 1, 3], [2, edge_count])
    real(dp), parameter :: node_x(4) = [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
    real(dp), parameter :: node_y(4) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: mesh_color(3) = [0.65_dp, 0.65_dp, 0.65_dp]
    real(dp), parameter :: node_color(3) = [0.05_dp, 0.05_dp, 0.05_dp]
    real(dp), parameter :: arrow_color(3) = [0.12_dp, 0.35_dp, 0.75_dp]
    ! Add a discrete gradient to the physical cotree potential.  The curl
    ! block cannot see this component; the direct tree gauge removes it.
    real(dp), parameter :: target_solution(edge_count) = [ &
        -0.35_dp, 0.50_dp, -0.60_dp, 1.20_dp, -0.25_dp]
    real(dp), parameter :: physical_solution(edge_count) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.75_dp, -0.40_dp]
    ! Face-edge incidence for the two triangular faces of the square split.
    ! Its Gram matrix is the topological curl-curl block C^T C.
    real(dp), parameter :: curl_incidence(2, edge_count) = reshape([ &
        1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        0.0_dp, 1.0_dp, -1.0_dp, 1.0_dp], [2, edge_count])
    real(dp), parameter :: full_matrix(edge_count, edge_count) = reshape([ &
        1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, &
        1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
        -1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], &
        [edge_count, edge_count])
    character(*), parameter :: output_directory = &
        "output/example/biro_tree_cotree_benchmark"
    type(tree_cotree_gauge_t) :: gauge
    real(dp), allocatable :: expected_restricted(:), restricted(:), prolonged(:)
    real(dp), allocatable :: reduced_matrix(:, :), reduced_rhs(:)
    real(dp) :: face_circulation(2), full_rhs(edge_count), residual
    real(dp) :: midpoint_x(edge_count), midpoint_y(edge_count)
    real(dp) :: arrow_u(edge_count), arrow_v(edge_count)
    integer :: edge, first_node, second_node, status, unit

    call execute_command_line("mkdir -p "//output_directory)
    if (maxval(abs(full_matrix - matmul(transpose(curl_incidence), &
        curl_incidence))) > 1.0e-12_dp) error stop "curl-curl block mismatch"
    full_rhs = matmul(full_matrix, target_solution)
    call build_tree_cotree_gauge(incidence, gauge, status)
    if (status /= 0) error stop "Biro tree-cotree gauge construction failed"
    call apply_tree_cotree_restriction( &
        gauge, target_solution, expected_restricted, status)
    if (status /= 0) error stop "Biro tree-cotree restriction failed"
    call reduce_tree_cotree_dense_system( &
        gauge, full_matrix, full_rhs, reduced_matrix, reduced_rhs, status)
    if (status /= 0) error stop "Biro tree-cotree direct reduction failed"
    call solve_reduced_2x2(reduced_matrix, reduced_rhs, restricted, status)
    if (status /= 0) error stop "Biro reduced direct solve failed"
    if (maxval(abs(restricted - physical_solution(4:5))) > 1.0e-12_dp) &
        error stop "Biro reduced solve did not recover cotree solution"
    call apply_tree_cotree_prolongation(gauge, restricted, prolonged, status)
    if (status /= 0) error stop "Biro tree-cotree prolongation failed"
    face_circulation = matmul(curl_incidence, prolonged)
    residual = maxval(abs(matmul(full_matrix, prolonged) - full_rhs))
    if (residual > 1.0e-12_dp) error stop "Biro manufactured residual failed"

    open (newunit=unit, file=output_directory//"/solution.csv", &
        status="replace", action="write")
    write (unit, "(a)") "edge,x0,y0,x1,y1,A_edge"
    do edge = 1, edge_count
        first_node = edge_nodes(1, edge)
        second_node = edge_nodes(2, edge)
        midpoint_x(edge) = 0.5_dp*(node_x(first_node) + node_x(second_node))
        midpoint_y(edge) = 0.5_dp*(node_y(first_node) + node_y(second_node))
        arrow_u(edge) = 0.18_dp*prolonged(edge)* &
            (node_x(second_node) - node_x(first_node))
        arrow_v(edge) = 0.18_dp*prolonged(edge)* &
            (node_y(second_node) - node_y(first_node))
        write (unit, "(*(es24.16,:,','))") real(edge, dp), &
            node_x(first_node), node_y(first_node), node_x(second_node), &
            node_y(second_node), prolonged(edge)
    end do
    close (unit)
    open (newunit=unit, file=output_directory//"/diagnostics.csv", &
        status="replace", action="write")
    write (unit, "(a,es24.16)") "manufactured_residual,", residual
    write (unit, "(a,es24.16)") "face_circulation_1,", face_circulation(1)
    write (unit, "(a,es24.16)") "face_circulation_2,", face_circulation(2)
    write (unit, "(a,es24.16)") "direct_cotree_energy,", &
        0.5_dp*dot_product(restricted, matmul(reduced_matrix, restricted))
    close (unit)

    ! Physical solution first: solved edge potential and face circulation.
    call figure(figsize=[7.2_dp, 6.2_dp])
    do edge = 1, edge_count
        first_node = edge_nodes(1, edge)
        second_node = edge_nodes(2, edge)
        call plot([node_x(first_node), node_x(second_node)], &
            [node_y(first_node), node_y(second_node)], color=mesh_color, &
            linewidth=1.5_dp)
    end do
    call add_scatter(node_x, node_y, marker="o", color=node_color, &
        label="mesh vertices")
    call quiver(midpoint_x, midpoint_y, arrow_u, arrow_v, scale=1.0_dp, &
        scale_units="xy", angles="xy", color=arrow_color, width=0.004_dp, &
        headwidth=3.0_dp)
    call xlabel("x")
    call ylabel("y")
    call title("Bíró tree--cotree curl--curl solution")
    call legend()
    call savefig(output_directory//"/solution.png")

contains

    subroutine solve_reduced_2x2(matrix, rhs, solution, status)
        real(dp), intent(in) :: matrix(:, :), rhs(:)
        real(dp), allocatable, intent(out) :: solution(:)
        integer, intent(out) :: status
        real(dp) :: determinant

        status = 1
        allocate(solution(0))
        if (size(matrix, 1) /= 2 .or. size(matrix, 2) /= 2 .or. &
            size(rhs) /= 2) return
        determinant = matrix(1, 1)*matrix(2, 2) - &
            matrix(1, 2)*matrix(2, 1)
        if (abs(determinant) < 1.0e-14_dp) return
        deallocate(solution)
        allocate(solution(2))
        solution(1) = (matrix(2, 2)*rhs(1) - matrix(1, 2)*rhs(2))/determinant
        solution(2) = (matrix(1, 1)*rhs(2) - matrix(2, 1)*rhs(1))/determinant
        status = 0
    end subroutine solve_reduced_2x2
end program biro_tree_cotree_benchmark
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/biro_tree_cotree_benchmark/primary.png)

### solution.png

![solution.png](../../media/examples/biro_tree_cotree_benchmark/solution.png)

---

[← Back to all examples](../index.html)
