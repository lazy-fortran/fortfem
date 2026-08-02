program test_mixed_wave_3d_properties
    use, intrinsic :: iso_fortran_env, only: int32
    use check, only: check_condition, check_property, check_summary, &
        property_random_unit, property_rng_initialize, property_rng_t
    use fortfem_api, only: &
        advance_mixed_wave_midpoint, advance_mixed_wave_symplectic_euler
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    logical :: all_passed
    integer(int32) :: first_failed_seed, shrunk_seed

    call check_property( &
        "random 3D mixed-wave matrix oracles preserve ideal structure", &
        1357911_int32, 16, wave_case, all_passed, first_failed_seed, shrunk_seed)
    call check_condition(all_passed .and. first_failed_seed == 0_int32 .and. &
        shrunk_seed == 0_int32, &
        "random 3D mixed-wave property reports no failure seed")
    call check_summary("random 3D mixed-wave matrix properties")
    if (.not. all_passed) error stop 1

contains

    logical function wave_case(case_seed)
        integer(int32), intent(in) :: case_seed
        integer, parameter :: component_count = 3
        integer, parameter :: state_count = 2*component_count
        type(property_rng_t) :: rng
        type(fortsparse_status_t) :: status
        real(dp) :: mass_q(component_count, component_count)
        real(dp) :: mass_v(component_count, component_count)
        real(dp) :: coupling(component_count, component_count)
        real(dp) :: state(state_count), state_reverse(state_count)
        real(dp) :: midpoint_map(state_count, state_count)
        real(dp) :: midpoint_reverse(state_count, state_count)
        real(dp) :: euler_map(state_count, state_count)
        real(dp) :: midpoint_oracle(state_count, state_count)
        real(dp) :: euler_oracle(state_count, state_count)
        real(dp) :: symplectic_form(state_count, state_count)
        real(dp) :: energy_metric(state_count, state_count)
        real(dp) :: skew_port(state_count, state_count)
        real(dp) :: energy_initial, energy_final, time_step
        real(dp) :: midpoint_defect, euler_defect, energy_defect
        real(dp) :: reversal_defect, matrix_defect
        integer :: row, column

        wave_case = .false.
        call property_rng_initialize(rng, case_seed)
        call random_spd_matrix(rng, mass_q)
        call random_spd_matrix(rng, mass_v)
        call random_port_coupling(rng, coupling)
        time_step = 0.02_dp + 0.12_dp*property_random_unit(rng)
        do row = 1, state_count
            state(row) = 2.0_dp*property_random_unit(rng) - 1.0_dp
        end do

        call build_midpoint_map( &
            mass_q, mass_v, coupling, time_step, midpoint_map, status)
        if (status%code /= 0) return
        call build_midpoint_map( &
            mass_q, mass_v, coupling, -time_step, midpoint_reverse, status)
        if (status%code /= 0) return
        call midpoint_matrix_oracle( &
            mass_q, mass_v, coupling, time_step, midpoint_oracle)
        call symplectic_geometry( &
            mass_q, mass_v, coupling, energy_metric, symplectic_form, skew_port)
        matrix_defect = maxval(abs(midpoint_map - midpoint_oracle))
        reversal_defect = maxval(abs(matmul(midpoint_reverse, midpoint_map) - &
            eye_matrix(state_count)))
        energy_defect = maxval(abs(matmul(transpose(midpoint_map), &
            matmul(energy_metric, midpoint_map)) - energy_metric))
        midpoint_defect = maxval(abs(matmul(transpose(midpoint_map), &
            matmul(symplectic_form, midpoint_map)) - symplectic_form))
        if (matrix_defect > 2.0e-10_dp .or. reversal_defect > 2.0e-10_dp .or. &
                energy_defect > 3.0e-10_dp .or. midpoint_defect > 5.0e-9_dp) return
        if (maxval(abs(skew_port + transpose(skew_port))) > 1.0e-14_dp) return

        state_reverse = matmul(midpoint_map, state)
        energy_initial = wave_energy(energy_metric, state)
        energy_final = wave_energy(energy_metric, state_reverse)
        if (abs(energy_final - energy_initial) > 3.0e-10_dp) return

        call build_euler_map( &
            mass_q, mass_v, coupling, time_step, euler_map, status)
        if (status%code /= 0) return
        call euler_matrix_oracle( &
            mass_q, mass_v, coupling, time_step, euler_oracle)
        matrix_defect = maxval(abs(euler_map - euler_oracle))
        euler_defect = maxval(abs(matmul(transpose(euler_map), &
            matmul(symplectic_form, euler_map)) - symplectic_form))
        if (matrix_defect > 2.0e-10_dp .or. euler_defect > 5.0e-9_dp) return
        wave_case = .true.
    end function wave_case

    subroutine random_spd_matrix(rng, matrix)
        type(property_rng_t), intent(inout) :: rng
        real(dp), intent(out) :: matrix(:, :)

        real(dp) :: lower(size(matrix, 1), size(matrix, 2))
        integer :: row, column, inner

        lower = 0.0_dp
        do row = 1, size(matrix, 1)
            do column = 1, row
                if (row == column) then
                    lower(row, column) = 0.9_dp + &
                        0.5_dp*property_random_unit(rng)
                else
                    lower(row, column) = 0.25_dp*( &
                        2.0_dp*property_random_unit(rng) - 1.0_dp)
                end if
            end do
        end do
        matrix = 0.0_dp
        do row = 1, size(matrix, 1)
            do column = 1, size(matrix, 2)
                do inner = 1, size(matrix, 2)
                    matrix(row, column) = matrix(row, column) + &
                        lower(row, inner)*lower(column, inner)
                end do
            end do
        end do
    end subroutine random_spd_matrix

    subroutine random_port_coupling(rng, coupling)
        type(property_rng_t), intent(inout) :: rng
        real(dp), intent(out) :: coupling(:, :)

        integer :: row, column

        do row = 1, size(coupling, 1)
            do column = 1, size(coupling, 2)
                coupling(row, column) = 0.12_dp*( &
                    2.0_dp*property_random_unit(rng) - 1.0_dp)
                if (row == column) coupling(row, column) = &
                    1.0_dp + 0.2_dp*property_random_unit(rng)
            end do
        end do
    end subroutine random_port_coupling

    subroutine build_midpoint_map( &
            mass_q, mass_v, coupling, time_step, map, status)
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), coupling(:, :)
        real(dp), intent(in) :: time_step
        real(dp), intent(out) :: map(:, :)
        type(fortsparse_status_t), intent(out) :: status

        real(dp) :: q(size(mass_q, 1)), v(size(mass_v, 1))
        integer :: column

        map = 0.0_dp
        do column = 1, size(map, 2)
            q = 0.0_dp
            v = 0.0_dp
            if (column <= size(q)) then
                q(column) = 1.0_dp
            else
                v(column - size(q)) = 1.0_dp
            end if
            call advance_mixed_wave_midpoint( &
                mass_q, mass_v, coupling, time_step, q, v, status)
            if (status%code /= 0) return
            map(:, column) = [q, v]
        end do
    end subroutine build_midpoint_map

    subroutine build_euler_map( &
            mass_q, mass_v, coupling, time_step, map, status)
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), coupling(:, :)
        real(dp), intent(in) :: time_step
        real(dp), intent(out) :: map(:, :)
        type(fortsparse_status_t), intent(out) :: status

        real(dp) :: q(size(mass_q, 1)), v(size(mass_v, 1))
        integer :: column

        map = 0.0_dp
        do column = 1, size(map, 2)
            q = 0.0_dp
            v = 0.0_dp
            if (column <= size(q)) then
                q(column) = 1.0_dp
            else
                v(column - size(q)) = 1.0_dp
            end if
            call advance_mixed_wave_symplectic_euler( &
                mass_q, mass_v, coupling, time_step, q, v, status)
            if (status%code /= 0) return
            map(:, column) = [q, v]
        end do
    end subroutine build_euler_map

    subroutine midpoint_matrix_oracle( &
            mass_q, mass_v, coupling, time_step, map)
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), coupling(:, :)
        real(dp), intent(in) :: time_step
        real(dp), intent(out) :: map(:, :)

        real(dp) :: system(size(map, 1), size(map, 2))
        real(dp) :: rhs(size(map, 1), size(map, 2))
        real(dp) :: energy_metric(size(map, 1), size(map, 2))
        real(dp) :: skew_port(size(map, 1), size(map, 2))
        integer :: n

        n = size(mass_q, 1)
        call block_matrices(mass_q, mass_v, coupling, energy_metric, skew_port)
        system = energy_metric + 0.5_dp*time_step*skew_port
        rhs = energy_metric - 0.5_dp*time_step*skew_port
        call solve_matrix(system, rhs, map)
    end subroutine midpoint_matrix_oracle

    subroutine euler_matrix_oracle( &
            mass_q, mass_v, coupling, time_step, map)
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), coupling(:, :)
        real(dp), intent(in) :: time_step
        real(dp), intent(out) :: map(:, :)

        real(dp) :: inverse_q(size(mass_q, 1), size(mass_q, 2))
        real(dp) :: inverse_v(size(mass_v, 1), size(mass_v, 2))
        real(dp) :: identity(size(map, 1), size(map, 2))
        integer :: n

        n = size(mass_q, 1)
        call solve_matrix(mass_q, eye_matrix(n), inverse_q)
        call solve_matrix(mass_v, eye_matrix(n), inverse_v)
        identity = eye_matrix(2*n)
        map = 0.0_dp
        map(:n, :n) = identity(:n, :n) - time_step**2* &
            matmul(inverse_q, matmul(transpose(coupling), &
            matmul(inverse_v, coupling)))
        map(:n, n + 1:2*n) = -time_step*matmul(inverse_q, transpose(coupling))
        map(n + 1:2*n, :n) = time_step*matmul(inverse_v, coupling)
        map(n + 1:2*n, n + 1:2*n) = identity(n + 1:2*n, n + 1:2*n)
    end subroutine euler_matrix_oracle

    subroutine symplectic_geometry( &
            mass_q, mass_v, coupling, energy_metric, symplectic_form, skew_port)
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), coupling(:, :)
        real(dp), intent(out) :: energy_metric(:, :), symplectic_form(:, :)
        real(dp), intent(out) :: skew_port(:, :)

        integer :: n
        real(dp) :: inverse_c(size(coupling, 1), size(coupling, 2))

        n = size(mass_q, 1)
        call solve_matrix(coupling, eye_matrix(n), inverse_c)
        call block_matrices(mass_q, mass_v, coupling, energy_metric, skew_port)
        symplectic_form = 0.0_dp
        symplectic_form(:n, n + 1:2*n) = matmul(mass_q, &
            matmul(inverse_c, mass_v))
        symplectic_form(n + 1:2*n, :n) = -transpose( &
            symplectic_form(:n, n + 1:2*n))
    end subroutine symplectic_geometry

    subroutine block_matrices( &
            mass_q, mass_v, coupling, energy_metric, skew_port)
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), coupling(:, :)
        real(dp), intent(out) :: energy_metric(:, :), skew_port(:, :)

        integer :: n

        n = size(mass_q, 1)
        energy_metric = 0.0_dp
        energy_metric(:n, :n) = mass_q
        energy_metric(n + 1:2*n, n + 1:2*n) = mass_v
        skew_port = 0.0_dp
        skew_port(:n, n + 1:2*n) = transpose(coupling)
        skew_port(n + 1:2*n, :n) = -coupling
    end subroutine block_matrices

    subroutine solve_matrix(matrix, rhs, solution)
        real(dp), intent(in) :: matrix(:, :), rhs(:, :)
        real(dp), intent(out) :: solution(:, :)

        real(dp) :: augmented(size(matrix, 1), size(matrix, 1) + size(rhs, 2))
        real(dp) :: pivot, factor, row_work(size(augmented, 2))
        integer :: dimension, rhs_count, row, column, pivot_row

        dimension = size(matrix, 1)
        rhs_count = size(rhs, 2)
        augmented(:, :dimension) = matrix
        augmented(:, dimension + 1:dimension + rhs_count) = rhs
        do column = 1, dimension
            pivot_row = column - 1 + maxloc( &
                abs(augmented(column:dimension, column)), dim=1)
            pivot = augmented(pivot_row, column)
            if (abs(pivot) < 1.0e-12_dp) error stop "matrix oracle singular"
            if (pivot_row /= column) then
                row_work = augmented(column, :)
                augmented(column, :) = augmented(pivot_row, :)
                augmented(pivot_row, :) = row_work
            end if
            augmented(column, :) = augmented(column, :)/pivot
            do row = 1, dimension
                if (row == column) cycle
                factor = augmented(row, column)
                augmented(row, :) = augmented(row, :) - factor*augmented(column, :)
            end do
        end do
        solution = augmented(:, dimension + 1:dimension + rhs_count)
    end subroutine solve_matrix

    pure function eye_matrix(dimension) result(identity)
        integer, intent(in) :: dimension
        real(dp) :: identity(dimension, dimension)
        integer :: diagonal

        identity = 0.0_dp
        do diagonal = 1, dimension
            identity(diagonal, diagonal) = 1.0_dp
        end do
    end function eye_matrix

    pure function wave_energy(metric, state) result(energy)
        real(dp), intent(in) :: metric(:, :), state(:)
        real(dp) :: energy

        energy = 0.5_dp*dot_product(state, matmul(metric, state))
    end function wave_energy

end program test_mixed_wave_3d_properties
