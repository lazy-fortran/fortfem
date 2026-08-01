program mixed_elasticity_wave
    !! Structure-preserving mixed 1-D elastic bar modal example.
    !!
    !! The two sine modes are a manufactured compatible bar state.  The
    !! mixed midpoint step advances displacement-like coordinates and their
    !! velocity companions; stress is reconstructed from the strain map.
    use fortfem_api, only: &
        assemble_mixed_elasticity_residual, advance_mixed_wave_midpoint
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, set_yscale, title, &
        xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: mode_count = 2, step_count = 400, sample_count = 161
    real(dp), parameter :: time_step = 0.01_dp
    real(dp), parameter :: frequencies(mode_count) = [ &
        acos(-1.0_dp), 2.0_dp*acos(-1.0_dp)]
    real(dp), parameter :: initial_displacement(mode_count) = [0.7_dp, -0.25_dp]
    real(dp), parameter :: initial_velocity(mode_count) = [0.15_dp, 0.35_dp]
    real(dp), parameter :: blue(3) = [0.0_dp, 114.0_dp, 178.0_dp]/255.0_dp
    real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
    real(dp), parameter :: green(3) = [0.0_dp, 158.0_dp, 115.0_dp]/255.0_dp
    character(*), parameter :: output_directory = &
        "output/example/mixed_elasticity_wave"
    real(dp) :: mass_displacement(mode_count, mode_count)
    real(dp) :: mass_velocity(mode_count, mode_count)
    real(dp) :: coupling(mode_count, mode_count)
    real(dp) :: compliance(mode_count, mode_count)
    real(dp) :: strain_map(mode_count, mode_count)
    real(dp) :: divergence_map(mode_count, mode_count)
    real(dp) :: displacement(mode_count), velocity(mode_count)
    real(dp) :: displacement_reverse(mode_count), velocity_reverse(mode_count)
    real(dp) :: load(mode_count), stress_modes(mode_count)
    real(dp) :: constitutive_residual(mode_count), equilibrium_residual(mode_count)
    real(dp) :: x(sample_count), displacement_field(sample_count)
    real(dp) :: velocity_field(sample_count), stress_field(sample_count)
    real(dp) :: time(step_count + 1), energy_history(step_count + 1)
    real(dp) :: error_history(step_count + 1)
    real(dp) :: displacement_exact(mode_count), velocity_exact(mode_count)
    real(dp) :: maximum_error, maximum_energy_drift, maximum_residual
    real(dp) :: start_time, elapsed, energy_initial, energy
    integer :: mode, point, step, unit, residual_status
    type(fortsparse_status_t) :: status

    call execute_command_line("mkdir -p "//output_directory)
    mass_displacement = 0.0_dp
    mass_velocity = 0.0_dp
    coupling = 0.0_dp
    compliance = 0.0_dp
    strain_map = 0.0_dp
    divergence_map = 0.0_dp
    do mode = 1, mode_count
        mass_displacement(mode, mode) = 1.0_dp
        mass_velocity(mode, mode) = 1.0_dp
        coupling(mode, mode) = frequencies(mode)
        compliance(mode, mode) = 1.0_dp
        strain_map(mode, mode) = frequencies(mode)
        divergence_map(mode, mode) = 1.0_dp
    end do

    displacement = initial_displacement
    velocity = initial_velocity
    energy_initial = 0.5_dp*(dot_product(displacement, displacement) + &
        dot_product(velocity, velocity))
    time(1) = 0.0_dp
    energy_history(1) = 0.0_dp
    error_history(1) = 0.0_dp
    call cpu_time(start_time)
    do step = 1, step_count
        call advance_mixed_wave_midpoint( &
            mass_displacement, mass_velocity, coupling, time_step, &
            displacement, velocity, status)
        if (status%code /= 0) error stop "mixed elasticity midpoint failed"
        time(step + 1) = real(step, dp)*time_step
        energy = 0.5_dp*(dot_product(displacement, displacement) + &
            dot_product(velocity, velocity))
        energy_history(step + 1) = energy/energy_initial - 1.0_dp
        call exact_modal_state(time(step + 1), displacement_exact, velocity_exact)
        error_history(step + 1) = max(maxval(abs(displacement - displacement_exact)), &
            maxval(abs(velocity - velocity_exact)))

        stress_modes = matmul(strain_map, displacement)
        load = stress_modes
        call assemble_mixed_elasticity_residual( &
            compliance, strain_map, divergence_map, stress_modes, displacement, load, &
            constitutive_residual, equilibrium_residual, residual_status)
        if (residual_status /= 0) error stop "mixed elasticity residual failed"
        if (max(maxval(abs(constitutive_residual)), &
            maxval(abs(equilibrium_residual))) > 2.0e-13_dp) &
            error stop "mixed elasticity residual oracle failed"
    end do
    call cpu_time(elapsed)
    elapsed = elapsed - start_time

    maximum_error = maxval(error_history)
    maximum_energy_drift = maxval(abs(energy_history))
    displacement_reverse = displacement
    velocity_reverse = velocity
    do step = 1, step_count
        call advance_mixed_wave_midpoint( &
            mass_displacement, mass_velocity, coupling, -time_step, &
            displacement_reverse, velocity_reverse, status)
        if (status%code /= 0) error stop "mixed elasticity reverse step failed"
    end do
    if (maxval(abs(displacement_reverse - initial_displacement)) > 2.0e-12_dp .or. &
        maxval(abs(velocity_reverse - initial_velocity)) > 2.0e-12_dp) &
        error stop "mixed elasticity reversibility oracle failed"
    if (maximum_error > 5.0e-3_dp) &
        error stop "mixed elasticity analytical oracle failed"
    if (maximum_energy_drift > 2.0e-12_dp) error stop "mixed elasticity energy failed"
    maximum_residual = max(maxval(abs(constitutive_residual)), &
        maxval(abs(equilibrium_residual)))

    do point = 1, sample_count
        x(point) = real(point - 1, dp)/real(sample_count - 1, dp)
    end do
    call reconstruct_fields(x, displacement, velocity, displacement_field, &
        velocity_field, stress_field)
    call render_plots()
    open (newunit=unit, file=output_directory//"/benchmark.txt", &
        status="replace", action="write")
    write (unit, "(a)") "quantity,value"
    write (unit, "(a,es24.16)") "time_step,", time_step
    write (unit, "(a,i0)") "step_count,", step_count
    write (unit, "(a,es24.16)") "maximum_exact_state_error,", maximum_error
    write (unit, "(a,es24.16)") "maximum_relative_energy_drift,", maximum_energy_drift
    write (unit, "(a,es24.16)") "mixed_residual_norm,", maximum_residual
    write (unit, "(a,es24.16)") "midpoint_wall_time_seconds,", elapsed
    close (unit)

contains

    subroutine exact_modal_state(current_time, displacement_value, velocity_value)
        real(dp), intent(in) :: current_time
        real(dp), intent(out) :: displacement_value(:), velocity_value(:)
        integer :: local_mode

        do local_mode = 1, mode_count
            displacement_value(local_mode) = &
                initial_displacement(local_mode)* &
                cos(frequencies(local_mode)*current_time) - &
                initial_velocity(local_mode)*sin(frequencies(local_mode)*current_time)
            velocity_value(local_mode) = &
                initial_velocity(local_mode)* &
                cos(frequencies(local_mode)*current_time) + &
                initial_displacement(local_mode)* &
                sin(frequencies(local_mode)*current_time)
        end do
    end subroutine exact_modal_state

    subroutine reconstruct_fields(coordinates, displacement_value, velocity_value, &
            displacement_value_field, velocity_value_field, stress_value_field)
        real(dp), intent(in) :: coordinates(:), displacement_value(:), velocity_value(:)
        real(dp), intent(out) :: displacement_value_field(:), velocity_value_field(:)
        real(dp), intent(out) :: stress_value_field(:)
        integer :: local_point

        do local_point = 1, size(coordinates)
            displacement_value_field(local_point) = &
                displacement_value(1)*sin(frequencies(1)*coordinates(local_point)) + &
                displacement_value(2)*sin(frequencies(2)*coordinates(local_point))
            velocity_value_field(local_point) = &
                velocity_value(1)*sin(frequencies(1)*coordinates(local_point)) + &
                velocity_value(2)*sin(frequencies(2)*coordinates(local_point))
            stress_value_field(local_point) = &
                frequencies(1)*displacement_value(1)* &
                cos(frequencies(1)*coordinates(local_point)) + &
                frequencies(2)*displacement_value(2)* &
                cos(frequencies(2)*coordinates(local_point))
        end do
    end subroutine reconstruct_fields

    subroutine render_plots()
        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(x, displacement_field, label="displacement u(x)", color=blue)
        call plot(x, stress_field/max(1.0_dp, maxval(abs(stress_field))), &
            label="normalized stress sigma(x)", color=orange, linestyle="--")
        call xlabel("bar coordinate x")
        call ylabel("field amplitude")
        call title("Mixed elastic bar solution at final time")
        call legend()
        call savefig(output_directory//"/mixed_elasticity_solution_1d.png")

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(x, stress_field, label="stress sigma(x)", color=orange)
        call plot(x, velocity_field, label="velocity v(x)", color=blue, linestyle="--")
        call xlabel("bar coordinate x")
        call ylabel("physical amplitude")
        call title("Mixed elastic bar stress and velocity")
        call legend()
        call savefig(output_directory//"/mixed_elasticity_stress_1d.png")

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(time, energy_history, label="relative energy drift", color=green)
        call xlabel("time")
        call ylabel("E(t)/E(0) - 1")
        call title("Elastic midpoint energy invariant")
        call savefig(output_directory//"/mixed_elasticity_energy_1d.png")

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(time, max(error_history, tiny(1.0_dp)), &
            label="state error", color=orange)
        call set_yscale("log")
        call xlabel("time")
        call ylabel("absolute error")
        call title("Mixed elastic bar analytical comparison")
        call savefig(output_directory//"/mixed_elasticity_error_1d.png")
    end subroutine render_plots

end program mixed_elasticity_wave
