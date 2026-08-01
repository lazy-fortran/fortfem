program test_mixed_wave_properties
    use, intrinsic :: iso_fortran_env, only: int32
    use check, only: check_condition, check_property, check_summary, &
        property_random_integer, property_random_unit, property_rng_initialize, &
        property_rng_t
    use fortfem_api, only: advance_mixed_wave_midpoint, advance_mixed_wave_strang
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    logical :: all_passed
    integer(int32) :: first_failed_seed, shrunk_seed

    call check_property( &
        "random mixed-wave midpoint and Strang structure preservation", &
        271828_int32, 16, wave_case, all_passed, first_failed_seed, shrunk_seed)
    call check_condition(all_passed .and. first_failed_seed == 0_int32 .and. &
        shrunk_seed == 0_int32, &
        "random mixed-wave property reports no failure seed")
    call check_summary("random mixed-wave structure properties")
    if (.not. all_passed) error stop 1

contains

    logical function wave_case(case_seed)
        integer(int32), intent(in) :: case_seed
        type(property_rng_t) :: rng
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: mass_q(:, :), mass_v(:, :)
        real(dp), allocatable :: coupling_a(:, :), coupling_b(:, :)
        real(dp), allocatable :: q(:), v(:), q_initial(:), v_initial(:)
        real(dp) :: time_step, energy_initial, energy_final
        integer :: dimension, row, column

        wave_case = .false.
        call property_rng_initialize(rng, case_seed)
        dimension = property_random_integer(rng, 1, 3)
        allocate(mass_q(dimension, dimension), mass_v(dimension, dimension), &
            coupling_a(dimension, dimension), coupling_b(dimension, dimension), &
            q(dimension), v(dimension), q_initial(dimension), v_initial(dimension))
        mass_q = 0.0_dp
        mass_v = 0.0_dp
        do row = 1, dimension
            mass_q(row, row) = 0.7_dp + 1.1_dp*property_random_unit(rng)
            mass_v(row, row) = 0.8_dp + 1.0_dp*property_random_unit(rng)
            q(row) = 2.0_dp*property_random_unit(rng) - 1.0_dp
            v(row) = 2.0_dp*property_random_unit(rng) - 1.0_dp
            do column = 1, dimension
                coupling_a(row, column) = 0.6_dp*( &
                    2.0_dp*property_random_unit(rng) - 1.0_dp)
                coupling_b(row, column) = 0.6_dp*( &
                    2.0_dp*property_random_unit(rng) - 1.0_dp)
            end do
        end do
        q_initial = q
        v_initial = v
        time_step = 0.03_dp + 0.2_dp*property_random_unit(rng)
        energy_initial = wave_energy(mass_q, mass_v, q, v)
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling_a, time_step, q, v, status)
        if (status%code /= 0) return
        energy_final = wave_energy(mass_q, mass_v, q, v)
        if (abs(energy_final - energy_initial) > 2.0e-11_dp) return
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling_a, -time_step, q, v, status)
        if (status%code /= 0) return
        if (maxval(abs(q - q_initial)) > 2.0e-11_dp .or. &
            maxval(abs(v - v_initial)) > 2.0e-11_dp) return

        q = q_initial
        v = v_initial
        energy_initial = wave_energy(mass_q, mass_v, q, v)
        call advance_mixed_wave_strang( &
            mass_q, mass_v, coupling_a, coupling_b, time_step, q, v, status)
        if (status%code /= 0) return
        energy_final = wave_energy(mass_q, mass_v, q, v)
        if (abs(energy_final - energy_initial) > 4.0e-11_dp) return
        call advance_mixed_wave_strang( &
            mass_q, mass_v, coupling_a, coupling_b, -time_step, q, v, status)
        if (status%code /= 0) return
        if (maxval(abs(q - q_initial)) > 4.0e-11_dp .or. &
            maxval(abs(v - v_initial)) > 4.0e-11_dp) return
        wave_case = .true.
    end function wave_case

    pure function wave_energy(mass_q, mass_v, q, v) result(energy)
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), q(:), v(:)
        real(dp) :: energy

        energy = 0.5_dp*(dot_product(q, matmul(mass_q, q)) + &
            dot_product(v, matmul(mass_v, v)))
    end function wave_energy

end program test_mixed_wave_properties
