---
title: tetra_nedelec_p_convergence Example
---

# tetra_nedelec_p_convergence Example

Executable FortFEM tetra_nedelec_p_convergence.f90 example.

## Usage

```bash
fpm run --example tetra_nedelec_p_convergence
```

## Source Code

```fortran
program tetra_nedelec_p_convergence
    use fortfem_api, only: &
        evaluate_tetra_nedelec_first_kind, &
        initialize_tetra_nedelec_first_kind, &
        interpolate_reference_tetra_nedelec, tetra_duffy_quadrature, &
        tetra_nedelec_dof_count, tetra_nedelec_first_kind_t
    use fortfem_kinds, only: dp
    use fortplot, only: add_scatter, figure, plot, savefig, set_yscale, title, &
        xlabel, ylabel
    implicit none

    character(*), parameter :: output_directory = &
        "output/example/tetra_nedelec_p_convergence"
    type(tetra_nedelec_first_kind_t) :: basis
    real(dp) :: degrees(4), errors(4)
    integer :: command_status, order, status

    call execute_command_line( &
        "mkdir -p "//output_directory, exitstat=command_status)
    if (command_status /= 0) error stop "cannot create example output directory"

    print "(a)", "order  L2 error for grad(x^4)"
    do order = 1, 4
        degrees(order) = real(order, dp)
        call initialize_tetra_nedelec_first_kind(order, basis, status)
        if (status /= 0) error stop "basis initialization failed"
        call interpolation_error(basis, errors(order), status)
        if (status /= 0) error stop "interpolation failed"
        print "(i5,2x,es14.6)", order, errors(order)
    end do
    if (any(errors(2:4) >= errors(1:3))) then
        error stop "interpolation error did not decrease"
    end if
    if (errors(4) >= 2.0e-11_dp) then
        error stop "order four did not reproduce the cubic field"
    end if

    call render_field()

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(degrees, errors, label="Nedelec interpolation error", &
        marker="o")
    call set_yscale("log")
    call xlabel("Nedelec polynomial order")
    call ylabel("L2 vector error")
    call title("Tetrahedral Nedelec p-convergence")
    call savefig(output_directory//"/p_convergence_1d.png")

contains

    subroutine render_field()
        real(dp), parameter :: vertices(3, 8) = reshape([ &
            0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
            0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [3, 8])
        real(dp) :: field_magnitude(size(vertices, 2))
        integer :: vertex

        do vertex = 1, size(vertices, 2)
            field_magnitude(vertex) = &
                norm2(cubic_gradient_value(vertices(:, vertex)))
        end do
        call figure(figsize=[7.5_dp, 6.0_dp])
        call add_scatter( &
            vertices(1, :), vertices(2, :), vertices(3, :), &
            c=field_magnitude, cmap="viridis", marker="o", &
            markersize=9.0_dp, label="|grad(x^4)| samples")
        call title("Tetrahedral Nedelec vector field")
        call savefig(output_directory//"/nedelec_field_3d.png")
    end subroutine render_field

    subroutine interpolation_error(basis, error, status)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(out) :: error
        integer, intent(out) :: status

        real(dp), allocatable :: curls(:, :), dofs(:), values(:, :)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: difference(3), point(3)
        integer :: node

        call interpolate_reference_tetra_nedelec( &
            basis, cubic_gradient, dofs, status)
        if (status /= 0) return
        allocate( &
            values(3, tetra_nedelec_dof_count(basis)), &
            curls(3, tetra_nedelec_dof_count(basis)))
        call tetra_duffy_quadrature(12, x, y, z, weights, status)
        if (status /= 0) return
        error = 0.0_dp
        do node = 1, size(weights)
            point = [x(node), y(node), z(node)]
            call evaluate_tetra_nedelec_first_kind( &
                basis, point, values, curls, status)
            if (status /= 0) return
            difference = matmul(values, dofs) - cubic_gradient_value(point)
            error = error + weights(node) * dot_product(difference, difference)
        end do
        error = sqrt(error)
    end subroutine interpolation_error

    pure subroutine cubic_gradient(point, value)
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: value(3)

        value = cubic_gradient_value(point)
    end subroutine cubic_gradient

    pure function cubic_gradient_value(point) result(value)
        real(dp), intent(in) :: point(3)
        real(dp) :: value(3)

        value = [4.0_dp * point(1)**3, 0.0_dp, 0.0_dp]
    end function cubic_gradient_value

end program tetra_nedelec_p_convergence
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/tetra_nedelec_p_convergence/primary.png)

### nedelec_field_3d.png

![nedelec_field_3d.png](../../media/examples/tetra_nedelec_p_convergence/nedelec_field_3d.png)

### p_convergence_1d.png

![p_convergence_1d.png](../../media/examples/tetra_nedelec_p_convergence/p_convergence_1d.png)

---

[← Back to all examples](../index.html)
