program test_mixed_wave_3d_structure_oracle
    !! Independent three-component manufactured wave structure oracle.
    !!
    !! The three diagonal blocks represent the x, y, and z components of a
    !! first-order port-Hamiltonian wave state.  The closed-form oscillator is
    !! deliberately independent of the time-step implementation.  The same
    !! fixture also checks the canonical two-form and proves that the separate
    !! dissipative Cayley block is energy-contracting but not symplectic.
    use check, only: check_condition, check_summary
    use fortfem_api, only: advance_dissipative_cayley, &
        advance_mixed_wave_midpoint, assemble_symplectic_map_defect
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: component_count = 3, state_count = 2*component_count
    integer, parameter :: step_count = 300
    real(dp), parameter :: time_step = 0.01_dp
    real(dp), parameter :: frequencies(component_count) = [ &
        0.75_dp, 1.25_dp, 1.80_dp]
    real(dp), parameter :: initial_q(component_count) = [ &
        0.80_dp, -0.40_dp, 0.20_dp]
    real(dp), parameter :: initial_v(component_count) = [ &
        -0.30_dp, 0.50_dp, 0.70_dp]
    real(dp), parameter :: tolerance = 3.0e-4_dp
    real(dp), parameter :: energy_tolerance = 3.0e-12_dp
    real(dp), parameter :: reversibility_tolerance = 3.0e-12_dp
    real(dp), parameter :: form_tolerance = 4.0e-13_dp
    real(dp), parameter :: damping_rate(component_count) = [ &
        0.20_dp, 0.35_dp, 0.50_dp]
    real(dp) :: mass_q(component_count, component_count)
    real(dp) :: mass_v(component_count, component_count)
    real(dp) :: coupling(component_count, component_count)
    real(dp) :: mass_state(state_count, state_count)
    real(dp) :: damping_state(state_count, state_count)
    real(dp) :: q(component_count), v(component_count)
    real(dp) :: q_reverse(component_count), v_reverse(component_count)
    real(dp) :: q_exact(component_count), v_exact(component_count)
    real(dp) :: state(state_count), state_next(state_count)
    real(dp) :: dissipative_state(state_count), dissipative_next(state_count)
    real(dp) :: map(state_count, state_count), dissipative_map(state_count, state_count)
    real(dp) :: symplectic_form(state_count, state_count)
    real(dp) :: defect(state_count, state_count)
    real(dp) :: dissipative_defect(state_count, state_count)
    real(dp) :: energy_initial, energy, maximum_error, maximum_energy_drift
    real(dp) :: dissipative_energy_initial, dissipative_energy_final, time
    integer :: component, step, column
    type(fortsparse_status_t) :: status

    mass_q = 0.0_dp
    mass_v = 0.0_dp
    coupling = 0.0_dp
    do component = 1, component_count
        mass_q(component, component) = 1.0_dp
        mass_v(component, component) = 1.0_dp
        coupling(component, component) = frequencies(component)
    end do

    q = initial_q
    v = initial_v
    energy_initial = wave_energy(q, v)
    maximum_error = 0.0_dp
    maximum_energy_drift = 0.0_dp
    do step = 1, step_count
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling, time_step, q, v, status)
        call check_condition(status%code == 0, &
            "3D manufactured midpoint accepts compatible wave blocks")
        time = real(step, dp)*time_step
        call exact_state(time, q_exact, v_exact)
        maximum_error = max(maximum_error, maxval(abs(q - q_exact)))
        maximum_error = max(maximum_error, maxval(abs(v - v_exact)))
        energy = wave_energy(q, v)
        maximum_energy_drift = max(maximum_energy_drift, &
            abs(energy - energy_initial))
    end do
    call check_condition(maximum_error < tolerance, &
        "3D midpoint follows the independent Cartesian oscillator solution")
    call check_condition(maximum_energy_drift < energy_tolerance, &
        "3D midpoint preserves the independent quadratic energy")

    q_reverse = q
    v_reverse = v
    do step = 1, step_count
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling, -time_step, q_reverse, v_reverse, status)
    end do
    call check_condition(maxval(abs(q_reverse - initial_q)) < &
        reversibility_tolerance .and. maxval(abs(v_reverse - initial_v)) < &
        reversibility_tolerance, &
        "3D midpoint is reversible under signed-step reversal")

    symplectic_form = 0.0_dp
    symplectic_form(:component_count, component_count + 1:state_count) = &
        eye(component_count)
    symplectic_form(component_count + 1:state_count, :component_count) = &
        -eye(component_count)
    map = 0.0_dp
    do column = 1, state_count
        state = 0.0_dp
        state(column) = 1.0_dp
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling, time_step, state(:component_count), &
            state(component_count + 1:state_count), status)
        map(:, column) = state
    end do
    call assemble_symplectic_map_defect( &
        map, symplectic_form, defect, status)
    call check_condition(status%code == 0 .and. maxval(abs(defect)) < form_tolerance, &
        "3D midpoint preserves the canonical symplectic two-form")

    mass_state = eye(state_count)
    damping_state = 0.0_dp
    do component = 1, component_count
        damping_state(component, component) = damping_rate(component)
        damping_state(component_count + component, component_count + component) = &
            damping_rate(component)
    end do
    dissipative_map = 0.0_dp
    do column = 1, state_count
        state = 0.0_dp
        state(column) = 1.0_dp
        call advance_dissipative_cayley( &
            mass_state, damping_state, time_step, state, state_next, status)
        dissipative_map(:, column) = state_next
    end do
    call assemble_symplectic_map_defect( &
        dissipative_map, symplectic_form, dissipative_defect, status)
    call check_condition(maxval(abs(dissipative_defect)) > 1.0e-6_dp, &
        "dissipative Cayley remains outside the symplectic contract")
    dissipative_state = [initial_q, initial_v]
    dissipative_energy_initial = 0.5_dp*dot_product( &
        dissipative_state, dissipative_state)
    call advance_dissipative_cayley( &
        mass_state, damping_state, time_step, dissipative_state, &
        dissipative_next, status)
    dissipative_energy_final = 0.5_dp*dot_product( &
        dissipative_next, dissipative_next)
    call check_condition(status%code == 0 .and. &
        dissipative_energy_final < dissipative_energy_initial, &
        "dissipative Cayley has a strictly decaying positive-time energy")

    call check_summary("3D manufactured mixed-wave structure oracle")

contains

    pure subroutine exact_state(current_time, q_value, v_value)
        real(dp), intent(in) :: current_time
        real(dp), intent(out) :: q_value(:), v_value(:)
        integer :: local_component

        do local_component = 1, component_count
            q_value(local_component) = initial_q(local_component)*cos( &
                frequencies(local_component)*current_time) - &
                initial_v(local_component)*sin( &
                frequencies(local_component)*current_time)
            v_value(local_component) = initial_v(local_component)*cos( &
                frequencies(local_component)*current_time) + &
                initial_q(local_component)*sin( &
                frequencies(local_component)*current_time)
        end do
    end subroutine exact_state

    pure function wave_energy(q_value, v_value) result(energy)
        real(dp), intent(in) :: q_value(:), v_value(:)
        real(dp) :: energy

        energy = 0.5_dp*(dot_product(q_value, q_value) + &
            dot_product(v_value, v_value))
    end function wave_energy

    pure function eye(dimension) result(identity)
        integer, intent(in) :: dimension
        real(dp) :: identity(dimension, dimension)
        integer :: diagonal

        identity = 0.0_dp
        do diagonal = 1, dimension
            identity(diagonal, diagonal) = 1.0_dp
        end do
    end function eye

end program test_mixed_wave_3d_structure_oracle
