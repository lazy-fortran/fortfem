---
title: iga_polar_feec Example
---

# iga_polar_feec Example

# Magnetic-axis isogeometric FEEC

This example constructs the Type-1 polar spline de Rham sequence

\[
H^1 \xrightarrow{\nabla} H(\mathrm{curl})
\xrightarrow{\mathrm{curl}} L^2
\]

on a periodic radial patch containing the magnetic axis. It assembles the
physical Piola-mapped Hodge operators and checks, for deterministic random
fields, that the discrete gradient and curl energies commute with their
corresponding Hodge maps.

The generated gallery artifacts show the polar control mesh, the energy
identity errors, and measured assembly time. No rendered image is committed.

## Provenance

- A. Toshniwal and T. J. R. Hughes, *Isogeometric discrete differential forms:
  Non-uniform degrees, Bézier extraction, polar splines and flows on surfaces*,
  CMAME 376 (2021), 113576,
  <https://doi.org/10.1016/j.cma.2020.113576>.
- F. Holderied, S. Possanner, and X. Wang, *MHD-kinetic hybrid code based on
  structure-preserving finite elements with particles-in-cell*, JCP 433
  (2021), 110143, <https://doi.org/10.1016/j.jcp.2021.110143>.

## Usage

```bash
fpm run --example iga_polar_feec
```

## Source Code

```fortran
program iga_polar_feec
    use fortfem_api, only: &
        assemble_bspline_polar_hcurl_operator_csc, &
        assemble_bspline_polar_h1_operator_csc, &
        assemble_bspline_polar_l2_mass_csc, &
        build_bspline_polar_feec_2d_operators
    use fortfem_kinds, only: dp
    use fortplot, only: figure, plot, savefig, set_yscale, title, xlabel, ylabel
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    integer, parameter :: azimuth_count = 16, radial_count = 5
    integer, parameter :: trial_count = 24
    real(dp), parameter :: radial_knots(7) = [ &
        0.0_dp, 0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp, 1.0_dp]
    character(*), parameter :: output_directory = &
        "output/example/iga_polar_feec"
    real(dp), allocatable :: control_points(:, :, :), curl(:, :)
    real(dp), allocatable :: edge_state(:), energy_error(:), face_state(:)
    real(dp), allocatable :: gradient(:, :), scalar_state(:), weights(:, :)
    real(dp) :: angle, assembly_seconds, curl_energy, h1_energy, start_time
    real(dp) :: hcurl_energy, l2_energy, trials(trial_count)
    type(csc_t) :: h1_stiffness, hcurl_curl, hcurl_mass, l2_mass
    type(fortsparse_status_t) :: sparse_status
    integer :: azimuth, radial, status, trial, unit

    call execute_command_line("mkdir -p "//output_directory)
    allocate( &
        control_points(2, azimuth_count, radial_count), &
        weights(azimuth_count, radial_count))
    weights = 1.0_dp
    do radial = 1, radial_count
        do azimuth = 1, azimuth_count
            angle = 2.0_dp*acos(-1.0_dp)* &
                real(azimuth, dp)/real(azimuth_count, dp)
            control_points(:, azimuth, radial) = &
                real(radial - 1, dp)/real(radial_count - 1, dp)* &
                [cos(angle), sin(angle)]
        end do
    end do
    call build_bspline_polar_feec_2d_operators( &
        azimuth_count, radial_count, gradient, curl, status)
    if (status /= 0) error stop "polar FEEC topology construction failed"

    call cpu_time(start_time)
    call assemble_bspline_polar_h1_operator_csc( &
        radial_knots, 1, azimuth_count, 1, control_points, weights, 3, &
        h1_stiffness, sparse_status, stiffness_coefficient=1.0_dp, &
        mass_coefficient=0.0_dp)
    if (sparse_status%code /= 0) error stop "polar H1 assembly failed"
    call assemble_bspline_polar_hcurl_operator_csc( &
        radial_knots, 1, azimuth_count, 1, control_points, weights, 3, &
        hcurl_mass, sparse_status, curl_coefficient=0.0_dp, &
        mass_coefficient=1.0_dp)
    if (sparse_status%code /= 0) error stop "polar Hcurl mass assembly failed"
    call assemble_bspline_polar_hcurl_operator_csc( &
        radial_knots, 1, azimuth_count, 1, control_points, weights, 3, &
        hcurl_curl, sparse_status, curl_coefficient=1.0_dp, &
        mass_coefficient=0.0_dp)
    if (sparse_status%code /= 0) error stop "polar curl-curl assembly failed"
    call assemble_bspline_polar_l2_mass_csc( &
        radial_knots, 1, azimuth_count, 1, control_points, weights, 3, &
        l2_mass, sparse_status)
    if (sparse_status%code /= 0) error stop "polar L2 assembly failed"
    call cpu_time(assembly_seconds)
    assembly_seconds = assembly_seconds - start_time

    allocate( &
        scalar_state(size(gradient, 2)), edge_state(size(gradient, 1)), &
        face_state(size(curl, 1)), energy_error(2*trial_count))
    call seed_random_numbers()
    do trial = 1, trial_count
        trials(trial) = real(trial, dp)
        call random_number(scalar_state)
        edge_state = matmul(gradient, scalar_state)
        h1_energy = dot_product( &
            scalar_state, csc_matvec(h1_stiffness, scalar_state))
        hcurl_energy = dot_product( &
            edge_state, csc_matvec(hcurl_mass, edge_state))
        energy_error(trial) = abs(h1_energy - hcurl_energy)
        call random_number(edge_state)
        face_state = matmul(curl, edge_state)
        curl_energy = dot_product( &
            edge_state, csc_matvec(hcurl_curl, edge_state))
        l2_energy = dot_product(face_state, csc_matvec(l2_mass, face_state))
        energy_error(trial_count + trial) = abs(curl_energy - l2_energy)
    end do
    if (maxval(energy_error) > 5.0e-11_dp) &
        error stop "polar FEEC commuting-energy regression"

    call render_mesh()
    call figure(figsize=[8.0_dp, 5.5_dp])
    call plot(trials, max(energy_error(:trial_count), epsilon(1.0_dp)), &
        label="grad/H(curl)")
    call plot(trials, max(energy_error(trial_count + 1:), epsilon(1.0_dp)), &
        label="curl/L2")
    call set_yscale("log")
    call xlabel("deterministic random trial")
    call ylabel("absolute energy identity error")
    call title("Magnetic-axis FEEC commuting Hodge identities")
    call savefig(output_directory//"/polar_feec_energy_1d.png")

    open (newunit=unit, file=output_directory//"/benchmark.txt", &
        status="replace", action="write")
    write (unit, "(a,i0)") "polar H1 degrees of freedom: ", size(gradient, 2)
    write (unit, "(a,i0)") "polar Hcurl degrees of freedom: ", size(gradient, 1)
    write (unit, "(a,i0)") "polar L2 degrees of freedom: ", size(curl, 1)
    write (unit, "(a,es14.6)") "four Hodge assemblies seconds: ", &
        assembly_seconds
    write (unit, "(a,es14.6)") "maximum commuting-energy error: ", &
        maxval(energy_error)
    close (unit)

contains

    subroutine render_mesh()
        real(dp) :: ring_x(azimuth_count + 1), ring_y(azimuth_count + 1)
        real(dp) :: spoke_x(radial_count), spoke_y(radial_count)

        call figure(figsize=[6.5_dp, 6.5_dp])
        do radial = 1, radial_count
            ring_x(:azimuth_count) = control_points(1, :, radial)
            ring_y(:azimuth_count) = control_points(2, :, radial)
            ring_x(azimuth_count + 1) = ring_x(1)
            ring_y(azimuth_count + 1) = ring_y(1)
            call plot(ring_x, ring_y)
        end do
        do azimuth = 1, azimuth_count
            spoke_x = control_points(1, azimuth, :)
            spoke_y = control_points(2, azimuth, :)
            call plot(spoke_x, spoke_y)
        end do
        call xlabel("R-like coordinate")
        call ylabel("Z-like coordinate")
        call title("Periodic polar spline patch with magnetic axis")
        call savefig(output_directory//"/polar_control_mesh_2d.png")
    end subroutine render_mesh

    subroutine seed_random_numbers()
        integer, allocatable :: seed(:)
        integer :: entry

        call random_seed(size=entry)
        allocate(seed(entry))
        seed = [(15485863 + 32452843*entry, entry = 1, size(seed))]
        call random_seed(put=seed)
    end subroutine seed_random_numbers

end program iga_polar_feec
```

## Generated Plots

*No plot artifact is produced by this example.*

---

[← Back to all examples](../index.html)
