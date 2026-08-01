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

The generated gallery artifacts lead with the computed scalar FEEC solution and
its physical gradient vectors. A separate plot evaluates the periodic spline
map at dense parameter samples, so it shows the curvilinear physical mesh
rather than only the control polygon. The energy identity errors and measured
assembly time are also reported. No rendered image is committed.

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
        build_bspline_polar_feec_2d_operators, &
        build_bspline_polar_h1_extraction, evaluate_bspline_basis, &
        evaluate_periodic_bspline_basis
    use fortfem_kinds, only: dp
    use fortplot, only: add_scatter, colorbar, figure, plot, quiver, savefig, &
        set_yscale, title, xlabel, ylabel
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    integer, parameter :: azimuth_count = 64, radial_count = 5
    integer, parameter :: trial_count = 24
    real(dp), parameter :: radial_knots(7) = [ &
        0.0_dp, 0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp, 1.0_dp]
    character(*), parameter :: output_directory = &
        "output/example/iga_polar_feec"
    integer, parameter :: plot_samples = 41, plot_line_samples = 121
    integer, parameter :: vector_stride = 4
    real(dp), allocatable :: control_points(:, :, :), curl(:, :)
    real(dp), allocatable :: edge_state(:), energy_error(:), face_state(:)
    real(dp), allocatable :: gradient(:, :), scalar_state(:), weights(:, :)
    real(dp), allocatable :: solution_state(:), tensor_solution(:)
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
        face_state(size(curl, 1)), energy_error(2*trial_count), &
        solution_state(size(gradient, 2)))
    call initialize_solution_state(solution_state, tensor_solution, status)
    if (status /= 0) error stop "polar FEEC solution construction failed"
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

    call render_solution()
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

    subroutine initialize_solution_state(polar_coefficients, tensor_coefficients, &
            local_status)
        real(dp), intent(out) :: polar_coefficients(:)
        real(dp), allocatable, intent(out) :: tensor_coefficients(:)
        integer, intent(out) :: local_status

        real(dp), allocatable :: extraction(:, :)
        real(dp) :: local_angle, local_radius
        integer :: local_azimuth, local_radial

        call build_bspline_polar_h1_extraction( &
            azimuth_count, radial_count, extraction, local_status)
        if (local_status /= 0) return
        if (size(polar_coefficients) /= size(extraction, 1)) then
            local_status = 1
            return
        end if
        polar_coefficients = 0.0_dp
        do local_radial = 1, radial_count - 2
            local_radius = real(local_radial, dp)/real(radial_count - 1, dp)
            do local_azimuth = 1, azimuth_count
                local_angle = 2.0_dp*acos(-1.0_dp)* &
                    (real(local_azimuth, dp) - 0.5_dp)/ &
                    real(azimuth_count, dp)
                polar_coefficients(3 + (local_radial - 1)*azimuth_count + &
                    local_azimuth) = local_radius**2*cos(2.0_dp*local_angle)
            end do
        end do
        tensor_coefficients = matmul(transpose(extraction), polar_coefficients)
        local_status = 0
    end subroutine initialize_solution_state

    subroutine render_solution()
        real(dp), allocatable :: x(:), y(:), values(:)
        real(dp), allocatable :: arrow_x(:), arrow_y(:), arrow_u(:), arrow_v(:)
        real(dp), allocatable :: radial_values(:), radial_derivatives(:)
        real(dp), allocatable :: angular_values(:), angular_derivatives(:)
        real(dp) :: coordinate_r, coordinate_a, jacobian(2, 2)
        real(dp) :: numerator, radial_derivative, angular_derivative
        real(dp) :: determinant, gradient_x, gradient_y, gradient_scale
        integer :: count, arrow_count, i, j, local_status, tensor_index

        allocate( &
            x(plot_samples**2), y(plot_samples**2), values(plot_samples**2), &
            arrow_x((plot_samples/vector_stride + 1)**2), &
            arrow_y((plot_samples/vector_stride + 1)**2), &
            arrow_u((plot_samples/vector_stride + 1)**2), &
            arrow_v((plot_samples/vector_stride + 1)**2))
        count = 0
        arrow_count = 0
        do j = 1, plot_samples
            coordinate_r = (real(j, dp) - 0.5_dp)/real(plot_samples, dp)
            do i = 1, plot_samples
                coordinate_a = (real(i, dp) - 0.5_dp)/real(plot_samples, dp)
                call evaluate_polar_point( &
                    coordinate_r, coordinate_a, x(count + 1), y(count + 1), &
                    jacobian, local_status)
                if (local_status /= 0) error stop "polar geometry evaluation failed"
                call evaluate_bspline_basis( &
                    radial_knots, 1, coordinate_r, radial_values, &
                    radial_derivatives, local_status)
                if (local_status /= 0) error stop "polar radial basis failed"
                call evaluate_periodic_bspline_basis( &
                    azimuth_count, 1, coordinate_a, angular_values, &
                    angular_derivatives, local_status)
                if (local_status /= 0) error stop "polar angular basis failed"
                numerator = 0.0_dp
                radial_derivative = 0.0_dp
                angular_derivative = 0.0_dp
                do tensor_index = 1, size(tensor_solution)
                    numerator = numerator + tensor_solution(tensor_index)* &
                        radial_values((tensor_index - 1)/azimuth_count + 1)* &
                        angular_values(modulo(tensor_index - 1, azimuth_count) + 1)
                    radial_derivative = radial_derivative + &
                        tensor_solution(tensor_index)* &
                        radial_derivatives((tensor_index - 1)/azimuth_count + 1)* &
                        angular_values(modulo(tensor_index - 1, azimuth_count) + 1)
                    angular_derivative = angular_derivative + &
                        tensor_solution(tensor_index)* &
                        radial_values((tensor_index - 1)/azimuth_count + 1)* &
                        angular_derivatives(modulo(tensor_index - 1, azimuth_count) + 1)
                end do
                count = count + 1
                values(count) = numerator
                determinant = jacobian(1, 1)*jacobian(2, 2) - &
                    jacobian(1, 2)*jacobian(2, 1)
                if (abs(determinant) > 128.0_dp*epsilon(1.0_dp)) then
                    gradient_x = (jacobian(2, 2)*radial_derivative - &
                        jacobian(2, 1)*angular_derivative)/determinant
                    gradient_y = (-jacobian(1, 2)*radial_derivative + &
                        jacobian(1, 1)*angular_derivative)/determinant
                else
                    gradient_x = 0.0_dp
                    gradient_y = 0.0_dp
                end if
                if (mod(i - 1, vector_stride) == 0 .and. &
                    j > plot_samples/8 .and. &
                    mod(j - 1, vector_stride) == 0) then
                    arrow_count = arrow_count + 1
                    arrow_x(arrow_count) = x(count)
                    arrow_y(arrow_count) = y(count)
                    arrow_u(arrow_count) = gradient_x
                    arrow_v(arrow_count) = gradient_y
                end if
            end do
        end do
        gradient_scale = maxval(sqrt(arrow_u(:arrow_count)**2 + &
            arrow_v(:arrow_count)**2))
        gradient_scale = max(gradient_scale, epsilon(1.0_dp))
        arrow_u(:arrow_count) = 0.12_dp*arrow_u(:arrow_count)/gradient_scale
        arrow_v(:arrow_count) = 0.12_dp*arrow_v(:arrow_count)/gradient_scale
        call figure(figsize=[8.0_dp, 7.0_dp])
        call add_scatter(x(:count), y(:count), c=values(:count), &
            cmap="coolwarm", marker=".", markersize=8.0_dp, &
            label="computed polar FEEC scalar solution")
        call quiver(arrow_x(:arrow_count), arrow_y(:arrow_count), &
            arrow_u(:arrow_count), arrow_v(:arrow_count), scale=1.0_dp, &
            scale_units="xy", angles="xy", color="black", width=0.003_dp, &
            headwidth=3.0_dp)
        call colorbar(label="scalar solution")
        call xlabel("physical x coordinate")
        call ylabel("physical y coordinate")
        call title("Computed polar FEEC solution and physical gradient")
        call savefig(output_directory//"/polar_feec_solution_2d.png")
    end subroutine render_solution

    subroutine render_mesh()
        real(dp) :: ring_x(plot_line_samples), ring_y(plot_line_samples)
        real(dp) :: spoke_x(plot_line_samples), spoke_y(plot_line_samples)
        real(dp) :: coordinate_r, coordinate_a
        real(dp) :: jacobian(2, 2)
        integer :: sample, local_status

        call figure(figsize=[6.5_dp, 6.5_dp])
        do radial = 2, radial_count
            coordinate_r = real(radial - 1, dp)/real(radial_count - 1, dp)
            do sample = 1, plot_line_samples
                coordinate_a = real(sample - 1, dp)/ &
                    real(plot_line_samples - 1, dp)
                call evaluate_polar_point( &
                    coordinate_r, coordinate_a, ring_x(sample), ring_y(sample), &
                    jacobian, local_status)
                if (local_status /= 0) error stop "polar ring evaluation failed"
            end do
            call plot(ring_x, ring_y)
        end do
        do azimuth = 1, azimuth_count
            coordinate_a = real(azimuth - 1, dp)/real(azimuth_count, dp)
            do sample = 1, plot_line_samples
                coordinate_r = real(sample - 1, dp)/ &
                    real(plot_line_samples - 1, dp)
                call evaluate_polar_point( &
                    coordinate_r, coordinate_a, spoke_x(sample), spoke_y(sample), &
                    jacobian, local_status)
                if (local_status /= 0) error stop "polar spoke evaluation failed"
            end do
            call plot(spoke_x, spoke_y)
        end do
        call xlabel("R-like coordinate")
        call ylabel("Z-like coordinate")
        call title("Curvilinear physical mesh from periodic polar splines")
        call savefig(output_directory//"/polar_curvilinear_mesh_2d.png")
    end subroutine render_mesh

    subroutine evaluate_polar_point(coordinate_r, coordinate_a, x, y, &
            jacobian, local_status)
        real(dp), intent(in) :: coordinate_r, coordinate_a
        real(dp), intent(out) :: x, y, jacobian(2, 2)
        integer, intent(out) :: local_status

        real(dp), allocatable :: radial_values(:), radial_derivatives(:)
        real(dp), allocatable :: angular_values(:), angular_derivatives(:)
        real(dp) :: denominator, denominator_derivative(2), factor
        real(dp) :: numerator(2), numerator_derivative(2, 2), point(2)
        integer :: local_azimuth, local_radial

        call evaluate_bspline_basis( &
            radial_knots, 1, coordinate_r, radial_values, radial_derivatives, &
            local_status)
        if (local_status /= 0) return
        call evaluate_periodic_bspline_basis( &
            azimuth_count, 1, coordinate_a, angular_values, &
            angular_derivatives, local_status)
        if (local_status /= 0) return
        denominator = 0.0_dp
        denominator_derivative = 0.0_dp
        numerator = 0.0_dp
        numerator_derivative = 0.0_dp
        do local_radial = 1, radial_count
            do local_azimuth = 1, azimuth_count
                factor = weights(local_azimuth, local_radial)* &
                    radial_values(local_radial)*angular_values(local_azimuth)
                denominator = denominator + factor
                numerator = numerator + factor* &
                    control_points(:, local_azimuth, local_radial)
                factor = weights(local_azimuth, local_radial)* &
                    radial_derivatives(local_radial)*angular_values(local_azimuth)
                denominator_derivative(1) = denominator_derivative(1) + factor
                numerator_derivative(:, 1) = numerator_derivative(:, 1) + &
                    factor*control_points(:, local_azimuth, local_radial)
                factor = weights(local_azimuth, local_radial)* &
                    radial_values(local_radial)*angular_derivatives(local_azimuth)
                denominator_derivative(2) = denominator_derivative(2) + factor
                numerator_derivative(:, 2) = numerator_derivative(:, 2) + &
                    factor*control_points(:, local_azimuth, local_radial)
            end do
        end do
        if (denominator <= 128.0_dp*epsilon(1.0_dp)) then
            local_status = 1
            return
        end if
        point = numerator/denominator
        do local_azimuth = 1, 2
            jacobian(:, local_azimuth) = ( &
                numerator_derivative(:, local_azimuth) - &
                point*denominator_derivative(local_azimuth))/denominator
        end do
        x = point(1)
        y = point(2)
        local_status = 0
    end subroutine evaluate_polar_point

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

### primary.png

![primary.png](../../media/examples/iga_polar_feec/primary.png)

### polar_control_mesh_2d.png

![polar_control_mesh_2d.png](../../media/examples/iga_polar_feec/polar_control_mesh_2d.png)

### polar_curvilinear_mesh_2d.png

![polar_curvilinear_mesh_2d.png](../../media/examples/iga_polar_feec/polar_curvilinear_mesh_2d.png)

### polar_feec_energy_1d.png

![polar_feec_energy_1d.png](../../media/examples/iga_polar_feec/polar_feec_energy_1d.png)

### polar_feec_solution_2d.png

![polar_feec_solution_2d.png](../../media/examples/iga_polar_feec/polar_feec_solution_2d.png)

---

[← Back to all examples](../index.html)
