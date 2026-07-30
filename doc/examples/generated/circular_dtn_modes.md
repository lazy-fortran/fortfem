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
    implicit none

    integer, parameter :: point_count = 32
    complex(dp) :: normal_derivative(point_count), trace(point_count)
    real(dp) :: angle, discarded_relative_norm
    integer :: point, status

    ! Boundary data on the unit circle:
    ! u(theta) = exp(2 i theta) + 0.1 exp(9 i theta).
    do point = 1, point_count
        angle = 2.0_dp * acos(-1.0_dp) * &
            real(point - 1, dp) / real(point_count, dp)
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
end program circular_dtn_modes
```

## Generated Plots

*No plot artifact is produced by this example.*

---

[← Back to all examples](../index.html)
