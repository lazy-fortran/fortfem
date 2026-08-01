---
title: fci_parallel_diffusion Example
---

# fci_parallel_diffusion Example

# PARALLAX-aligned FCI parallel diffusion

This small fixture exercises FortFEM's field-coordinate-independent support
operator with identity interpolation maps. It applies

\[
L_\parallel u=-W_c^{-1}Q^T W_sK_\parallel Q u
\]

to a smooth cosine profile on an open field-line segment, then checks the
independent telescoping mass oracle and the negative weighted-energy identity.
The zero-gradient end behaviour keeps the boundary flux small in the plot.
The maps are trivial by design; field-line tracing and interpolation geometry
remain separate services, as in the PARALLAX architecture.

Outputs:

- `fci_parallel_profile_1d.png`
- `fci_parallel_dissipation_1d.png`
- `fci_parallel_profile.csv`
- `fci_parallel_benchmark.csv` with the measured matrix-free action time on
  the local runner.

The implementation follows the support-operator construction described by
[Stegmeir et al.](https://doi.org/10.1016/j.cpc.2015.09.016). FortFEM does not
copy PARALLAX source code or benchmark data.

## Usage

```bash
fpm run --example fci_parallel_diffusion
```

## Source Code

```fortran
program fci_parallel_diffusion
    ! Small PARALLAX-aligned support-operator fixture.
    ! Identity maps make the one-dimensional telescoping conservation oracle
    ! transparent while retaining the full FCI P K_parallel Q API.  The
    ! cosine profile has zero continuous gradient at both open-line ends, so
    ! the gallery plot does not obscure the support action with boundary flux.
    use fortfem_api, only: apply_fci_parallel_diffusion
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: segment_count = 64
    integer, parameter :: benchmark_repetitions = 1000
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(*), parameter :: output_directory = &
        "output/example/fci_parallel_diffusion"
    real(dp) :: forward_map(1, 1, segment_count)
    real(dp) :: backward_map(1, 1, segment_count)
    real(dp) :: line_lengths(1, segment_count)
    real(dp) :: parallel_coefficient(segment_count)
    real(dp) :: canonical_volumes(segment_count + 1)
    real(dp) :: staggered_volumes(segment_count)
    real(dp) :: coordinate(segment_count + 1), field(segment_count + 1)
    real(dp) :: diffusion_field(segment_count + 1)
    real(dp) :: total_mass_rate, dissipation
    integer :: point, unit, repetition, clock_rate, start_count, end_count
    real(dp) :: action_seconds
    type(fortsparse_status_t) :: status

    call execute_command_line("mkdir -p "//output_directory)
    forward_map = 1.0_dp
    backward_map = 1.0_dp
    line_lengths = 1.0_dp
    parallel_coefficient = 20.0_dp
    canonical_volumes = 1.0_dp
    staggered_volumes = 1.0_dp
    do point = 1, segment_count + 1
        coordinate(point) = real(point - 1, dp)/real(segment_count, dp)
        field(point) = cos(pi*coordinate(point))
    end do

    call apply_fci_parallel_diffusion( &
        forward_map, backward_map, line_lengths, parallel_coefficient, &
        canonical_volumes, staggered_volumes, field, diffusion_field, status)
    if (status%code /= 0) error stop "FCI parallel diffusion failed"
    total_mass_rate = sum(diffusion_field*canonical_volumes)
    dissipation = dot_product(field*canonical_volumes, diffusion_field)
    if (abs(total_mass_rate) > 2.0e-13_dp) &
        error stop "FCI support operator lost total mass"
    if (dissipation > 2.0e-13_dp) &
        error stop "FCI support operator gained energy"

    call system_clock(count_rate=clock_rate)
    if (clock_rate > 0) then
        call system_clock(start_count)
        do repetition = 1, benchmark_repetitions
            call apply_fci_parallel_diffusion( &
                forward_map, backward_map, line_lengths, parallel_coefficient, &
                canonical_volumes, staggered_volumes, field, diffusion_field, &
                status)
            if (status%code /= 0) error stop "FCI benchmark action failed"
        end do
        call system_clock(end_count)
        action_seconds = real(end_count - start_count, dp) / &
            real(clock_rate, dp) / real(benchmark_repetitions, dp)
    else
        action_seconds = -1.0_dp
    end if

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(coordinate, field, label="field u")
    call plot(coordinate, diffusion_field, label="P K_parallel Q u")
    call xlabel("normalized toroidal coordinate")
    call ylabel("field / diffusion action")
    call title("PARALLAX-aligned FCI parallel diffusion")
    call legend()
    call savefig(output_directory//"/fci_parallel_profile_1d.png")

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(coordinate, diffusion_field, label="parallel diffusion")
    call xlabel("normalized toroidal coordinate")
    call ylabel("P K_parallel Q u")
    call title("FCI support-operator dissipation profile")
    call legend()
    call savefig(output_directory//"/fci_parallel_dissipation_1d.png")

    open (newunit=unit, file=output_directory//"/fci_parallel_profile.csv", &
        status="replace", action="write")
    write (unit, "(a)") "coordinate,field,diffusion"
    do point = 1, segment_count + 1
        write (unit, "(3(es24.16,1x))") coordinate(point), field(point), &
            diffusion_field(point)
    end do
    write (unit, "(a,es24.16)") "total_mass_rate,", total_mass_rate
    write (unit, "(a,es24.16)") "dissipation,", dissipation
    write (unit, "(a,es24.16)") "action_seconds,", action_seconds
    write (unit, "(a,i0)") "action_repetitions,", benchmark_repetitions
    close (unit)

    open (newunit=unit, file=output_directory//"/fci_parallel_benchmark.csv", &
        status="replace", action="write")
    write (unit, "(a)") "segment_count,staggered_rows,repetitions,action_seconds"
    write (unit, "(i0,',',i0,',',i0,',',es24.16)") segment_count, &
        size(staggered_volumes), benchmark_repetitions, action_seconds
    close (unit)
end program fci_parallel_diffusion
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/fci_parallel_diffusion/primary.png)

### fci_parallel_dissipation_1d.png

![fci_parallel_dissipation_1d.png](../../media/examples/fci_parallel_diffusion/fci_parallel_dissipation_1d.png)

### fci_parallel_profile_1d.png

![fci_parallel_profile_1d.png](../../media/examples/fci_parallel_diffusion/fci_parallel_profile_1d.png)

---

[← Back to all examples](../index.html)
