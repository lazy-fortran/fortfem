---
title: circular_dtn_modes Example
---

# circular_dtn_modes Example

Boundary data on the unit circle:

## Usage

```bash
fpm run --example circular_dtn_modes
```

## Source Code

```fortran
program circular_dtn_modes
    use fortfem_api, only: apply_circular_helmholtz_dtn
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, title, xlabel, ylabel
    implicit none

    integer, parameter :: point_count = 32
    character(*), parameter :: output_directory = &
        "output/example/circular_dtn_modes"
    complex(dp) :: normal_derivative(point_count), trace(point_count)
    real(dp) :: angle, discarded_relative_norm, theta(point_count)
    integer :: command_status, point, status

    call execute_command_line( &
        "mkdir -p "//output_directory, exitstat=command_status)
    if (command_status /= 0) error stop "cannot create example output directory"

    ! Boundary data on the unit circle:
    ! u(theta) = exp(2 i theta) + 0.1 exp(9 i theta).
    do point = 1, point_count
        angle = 2.0_dp * acos(-1.0_dp) * &
            real(point - 1, dp) / real(point_count, dp)
        theta(point) = angle
        trace(point) = exp(cmplx(0.0_dp, 2.0_dp * angle, dp)) + &
            0.1_dp * exp(cmplx(0.0_dp, 9.0_dp * angle, dp))
    end do

    ! Retaining |mode| <= 4 removes the small ninth harmonic. The diagnostic
    ! is the relative discrete L2 norm of the discarded boundary trace.
    call apply_circular_helmholtz_dtn( &
        trace, 3.0_dp, 1.0_dp, 4, normal_derivative, &
        discarded_relative_norm, status)
    if (status /= 0) error stop "Circular DtN application failed"

    write(*, '(a, es12.4)') &
        "discarded relative trace norm: ", discarded_relative_norm
    write(*, '(a, 2es14.5)') &
        "normal derivative at theta=0: ", normal_derivative(1)

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(theta, real(trace), label="Re(trace)", linestyle="-")
    call plot(theta, aimag(trace), label="Im(trace)", linestyle="--")
    call xlabel("boundary angle theta [rad]")
    call ylabel("boundary trace")
    call title("Circular Helmholtz DtN input trace")
    call legend()
    call savefig(output_directory//"/circular_dtn_trace_1d.png")

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(theta, real(normal_derivative), label="Re(DtN trace)", &
        marker="o")
    call plot(theta, aimag(normal_derivative), label="Im(DtN trace)", &
        marker="s")
    call xlabel("boundary angle theta [rad]")
    call ylabel("normal derivative")
    call title("Circular Helmholtz DtN response")
    call savefig(output_directory//"/circular_dtn_response_1d.png")
end program circular_dtn_modes
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/circular_dtn_modes/primary.png)

### circular_dtn_response_1d.png

![circular_dtn_response_1d.png](../../media/examples/circular_dtn_modes/circular_dtn_response_1d.png)

### circular_dtn_trace_1d.png

![circular_dtn_trace_1d.png](../../media/examples/circular_dtn_modes/circular_dtn_trace_1d.png)

---

[← Back to all examples](../index.html)
